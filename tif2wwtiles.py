#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import os
import sys

from osgeo import gdal
from osgeo import osr

try:
    from PIL import Image
    import numpy
    import osgeo.gdal_array as gdalarray
except:
    # 'antialias' resampling is not available
    pass

__version__ = "$Id: tif2wwtiles.py 00001 2017-01-09 10:36:12Z goatbar $"

resampling_list = ('average','near','bilinear','cubic','cubicspline','lanczos','antialias')
#profile_list = ('mercator','geodetic','raster') #,'zoomify')
#webviewer_list = ('all','google','openlayers','leaflet','none')


MAXZOOMLEVEL = 21

class Global(object):
    def __init__(self, tileSize=512, size=36):
	"""初始化模具"""
    	self.tileSize = tileSize
	self.zerolevelsize = size
	
    def GetLztsd(self, n):
	return self.zerolevelsize / math.pow(2, n)
	
    def LatToRow(self, n, lat):
	return int( (lat + 90) / self.GetLztsd( n ))
		
    def LonToCol(self, n, lon):
	return int( (lon + 180) / self.GetLztsd( n ))
		
    def BLtoRowCol(self, n, degree, mid):
	return (degree + mid) / self.GetLztsd( n )
	
    def RowColtoBL(self, n, rowcol, mid):
	return rowcol * self.GetLztsd( n ) - mid
		
    def GetPixelsPerDegree(self, n):
        return self.tileSize / self.GetLztsd( n )
		
    def RowColtoBLH(self, n, row, col):
	b = self.RowColtoBL(n, row, 90)
	l = self.RowColtoBL(n, col, 180)
	return (b, l)
	
    def BLHtoRowCol(self, n, b, l):
	row = self.LatToRow(n, b)
	col = self.LonToCol(n, l)
	return (row, col)
	
class Tif2WWTiles(object):
    def process(self):
        """The main processing function, runs all the main steps of processing"""

        # Opening and preprocessing of the input file
        self.open_input()

        # Generation of the lowest tiles
        self.generate_base_tiles()

        # Generation of the overview tiles (higher in the pyramid)
        self.generate_overview_tiles()
		
	# Generation of xmldata
        self.generate_xmldata()

    def error(self, msg, details = "" ):
        """Print an error message and stop the processing"""

        if details:
            self.parser.error(msg + "\n\n" + details)
        else:
            self.parser.error(msg)

    def progressbar(self, complete = 0.0):
        """Print progressbar for float value 0..1"""

        gdal.TermProgress_nocb(complete)
            
    def __init__(self, arguments):
        """Constructor function - initialization"""
	self.input = None
        self.output = None
        
	# Tile format
        self.tilesize = 512
        self.tiledriver = 'GTiff'
        self.tileext = 'tif'
	self.querysize = 4 * self.tilesize
	self.tools = Global()

        # RUN THE ARGUMENT PARSER:
        self.optparse_init()
        self.options, self.args = self.parser.parse_args(args=arguments)

        if not self.args:
            self.error("No input file specified")
            
        # Test output directory, if it doesn't exist
        if (len(self.args) > 3):
            self.output = self.args[3]
            self.args = self.args[:3]

        self.input = self.args[0]

        # Default values for not given options

        if not self.output:
            # Directory with input filename without extension in actual directory
            self.output = os.path.splitext(os.path.basename( self.input ))[0]
        
	# User specified zoom levels
	if (len(self.args) != 3):
            self.error("please input parameters in this way: filename.tif maxlevel minlevel")
        else :
            import string
	    self.tmaxz = string.atoi(self.args[1],10)
	    self.tminz = string.atoi(self.args[2],10)
        
	
    def optparse_init(self):
        """Prepare the option parser for input (argv)"""

        from optparse import OptionParser
        usage = "Usage: %prog [options] input_file(s) maxlevel minlevel [output]"
        p = OptionParser(usage, version="%prog "+ __version__)
        
        p.add_option("-r", "--resampling", dest="resampling", type='choice', choices=resampling_list,
                        help="Resampling method (%s) - default 'average'" % ",".join(resampling_list))
        p.add_option('-z', '--zoom', dest="zoom",
                          help="Zoom levels to render (format:'2-5' or '10').")
        p.add_option('-e', '--resume', dest="resume", action="store_true",
                          help="Resume mode. Generate only missing files.")
        p.add_option("-v", "--verbose",
                          action="store_true", dest="verbose",
                          help="Print status messages to stdout")
        p.add_option("-q", "--quiet",
                          action="store_true", dest="quiet",
                          help="Disable messages and status to stdout")

        p.set_defaults(verbose=False, resampling='average', resume=False)

        self.parser = p
		
    def open_input(self):
	gdal.AllRegister()

	# Initialize necessary GDAL drivers
        self.out_drv = gdal.GetDriverByName( self.tiledriver )
        self.mem_drv = gdal.GetDriverByName( 'MEM' )

        if not self.out_drv:
            raise Exception("The '%s' driver was not found, is it available in this GDAL build?", self.tiledriver)
        if not self.mem_drv:
            raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?")
    
	# Open the input file
        self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
        
	# Get base image data & info
	self.width = self.in_ds.RasterXSize
	self.height = self.in_ds.RasterYSize
	self.bandnum = self.in_ds.RasterCount

        # Spatial Reference System of tiles
        out_srs = osr.SpatialReference()
        out_srs.ImportFromEPSG(3857)

        # Spatial Reference System of input
        self.wkt = self.in_ds.GetProjection()
        srs_in = osr.SpatialReference()
	srs_in.ImportFromWkt(self.wkt)

	# WGS84 Reference System
	srs4326 = osr.SpatialReference()
	srs4326.ImportFromEPSG(4326)

        # Are the reference systems the same? Reproject if necessary.
        self.out_ds = None
        if (srs_in.ExportToProj4() != out_srs.ExportToProj4()):

            # Generation of VRT dataset in tile projection, default 'nearest neighbour' warping
            self.out_ds = gdal.AutoCreateWarpedVRT( self.in_ds, self.wkt, out_srs.ExportToWkt() )

        if not self.out_ds:
            self.out_ds = self.in_ds

	# Get alpha band (either directly or from NODATA value)
        self.alphaband = self.in_ds.GetRasterBand(1).GetMaskBand()

        # Get reference systems transformer
	self.ctutmto84 = osr.CoordinateTransformation(srs_in, srs4326)
	self.ct84toutm = osr.CoordinateTransformation(srs4326, srs_in)
	self.ctmerto84 = osr.CoordinateTransformation(out_srs, srs4326)
	self.ct84tomer = osr.CoordinateTransformation(srs4326, out_srs)
		
	# Get GeoTransform
	self.transform = self.out_ds.GetGeoTransform()
		
	# Get input bounds in WGS84 Coordinate System
	lu = self.MERtoWGS84(self.transform[0], self.transform[3])
	rl = self.MERtoWGS84(self.transform[0] + self.transform[1] * self.out_ds.RasterXSize, self.transform[3] + self.transform[5] * self.out_ds.RasterYSize)
	self.omaxx = rl[0]
        self.ominx = lu[0]
        self.omaxy = lu[1]
        self.ominy = rl[1]
		
	# Generate table with min max tile coordinates for all zoomlevels
	self.tminmax = list(range(0,MAXZOOMLEVEL))
	for tz in range(0, MAXZOOMLEVEL):
	    tminy, tminx = self.tools.BLHtoRowCol(tz, self.ominy, self.ominx )
	    tmaxy, tmaxx = self.tools.BLHtoRowCol(tz, self.omaxy, self.omaxx )
	    # crop tiles extending world limits (+-180,+-90)
	    tminx, tminy = max(0, tminx), max(0, tminy)
	    tmaxx, tmaxy = min(10 * 2**tz-1, tmaxx), min(5 * 2**tz-1, tmaxy)
	    self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)
	
    def generate_base_tiles(self):

	print("Generating Base Tiles:")
	
	# Set the bounds
        tminx, tminy, tmaxx, tmaxy = self.tminmax[self.tmaxz]

	ds = self.out_ds
        tilebands = ds.RasterCount
        querysize = self.querysize
		
	#print tminx, tminy, tmaxx, tmaxy
        tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
	#print tcount
        ti = 0
		
	tz = self.tmaxz
        for ty in range(tmaxy, tminy-1, -1): #range(tminy, tmaxy+1):
            for tx in range(tminx, tmaxx+1):
			
                ti += 1
                tilefilename = os.path.join(self.output, str(tz), str(ty), "%s_%s.%s" % (str(ty), tx, self.tileext))
                #if self.options.verbose:
                #    print(ti,'/',tcount, tilefilename) #, "( TileMapService: z / x / y )"

                #if self.options.resume and os.path.exists(tilefilename):
                #    if self.options.verbose:
                #        print("Tile generation skipped because of --resume")
                #    else:
                #        self.progressbar( ti / float(tcount) )
                #    continue
				
		# Create directories for the tile
                if not os.path.exists(os.path.dirname(tilefilename)):
                    os.makedirs(os.path.dirname(tilefilename))
					
		# Get tile left-up corner LatLon
		tilelubl = self.tools.RowColtoBLH(tz, ty + 1, tx)
		tilerlbl = self.tools.RowColtoBLH(tz, ty , tx + 1)
				
		# Get read & write bounds
		rb, wb = self.geo_query(ds, tilelubl[1], tilelubl[0], tilerlbl[1], tilerlbl[0])
		# Tile bounds in raster coordinates for ReadRaster query
                rb, wb = self.geo_query(ds, tilelubl[1], tilelubl[0], tilerlbl[1], tilerlbl[0], querysize=querysize)

                rx, ry, rxsize, rysize = rb
                wx, wy, wxsize, wysize = wb
		
		# Tile dataset in memory
		dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)
                data = ds.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize)
                #alpha = self.alphaband.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize)

                if self.tilesize == querysize:
                    # Use the ReadRaster result directly in tiles ('nearest neighbour' query)
                    dstile.WriteRaster(wx, wy, wxsize, wysize, data)
                    #dstile.WriteRaster(wx, wy, wxsize, wysize, alpha, band_list=[tilebands])

                    # Note: For source drivers based on WaveLet compression (JPEG2000, ECW, MrSID)
                    # the ReadRaster function returns high-quality raster (not ugly nearest neighbour)
                    # TODO: Use directly 'near' for WaveLet files
                else:
                    # Big ReadRaster query in memory scaled to the tilesize - all but 'near' algo
                    dsquery = self.mem_drv.Create('', querysize, querysize, tilebands)
                    # TODO: fill the null value in case a tile without alpha is produced (now only png tiles are supported)
                    #for i in range(1, tilebands+1):
                    #   dsquery.GetRasterBand(1).Fill(tilenodata)
                    dsquery.WriteRaster(wx, wy, wxsize, wysize, data)
                    #dsquery.WriteRaster(wx, wy, wxsize, wysize, alpha, band_list=[tilebands])

                    self.scale_query_to_tile(dsquery, dstile, tilefilename)
                    del dsquery

                del data
                
                # Write a copy of tile to png/jpg
                self.out_drv.CreateCopy(tilefilename, dstile, strict=0)
                
                self.progressbar( ti / float(tcount) )
                
                del dstile

    def generate_overview_tiles(self):
        """Generation of the overview tiles (higher in the pyramid) based on existing tiles"""

        print("Generating Overview Tiles:")

        tilebands = self.in_ds.RasterCount

        # Usage of existing tiles: from 4 underlying tiles generate one as overview.

        tcount = 0
        for tz in range(self.tmaxz-1, self.tminz-1, -1):
            tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
            tcount += (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))

        ti = 0

        # querysize = tilesize * 2

        for tz in range(self.tmaxz-1, self.tminz-1, -1):
            tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
            for ty in range(tmaxy, tminy-1, -1): #range(tminy, tmaxy+1):
                for tx in range(tminx, tmaxx+1):

                    ti += 1
                    tilefilename = os.path.join( self.output, str(tz), str(ty), "%s_%s.%s" % (str(ty), tx, self.tileext) )

                    #if self.options.verbose:
                    #    print(ti,'/',tcount, tilefilename) #, "( TileMapService: z / x / y )"

                    #if self.options.resume and os.path.exists(tilefilename):
                    #    if self.options.verbose:
                    #        print("Tile generation skipped because of --resume")
                    #    else:
                    #        self.progressbar( ti / float(tcount) )
                    #    continue

                    # Create directories for the tile
                    if not os.path.exists(os.path.dirname(tilefilename)):
                        os.makedirs(os.path.dirname(tilefilename))

                    dsquery = self.mem_drv.Create('', 2*self.tilesize, 2*self.tilesize, tilebands)
                    # TODO: fill the null value
                    #for i in range(1, tilebands+1):
                    #   dsquery.GetRasterBand(1).Fill(tilenodata)
                    dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)

                    # TODO: Implement more clever walking on the tiles with cache functionality
                    # probably walk should start with reading of four tiles from top left corner
                    # Hilbert curve

                    children = []
                    # Read the tiles and write them to query window
                    for y in range(2*ty,2*ty+2):
                        for x in range(2*tx,2*tx+2):
                            minx, miny, maxx, maxy = self.tminmax[tz+1]
                            if x >= minx and x <= maxx and y >= miny and y <= maxy:
                                dsquerytile = gdal.Open( os.path.join( self.output, str(tz+1), str(y), "%s_%s.%s" % (str(y), x, self.tileext)), gdal.GA_ReadOnly)
                                if (ty==0 and y==1) or (ty!=0 and (y % (2*ty)) != 0):
                                    tileposy = 0
                                else:
                                    tileposy = self.tilesize
                                if tx:
                                    tileposx = x % (2*tx) * self.tilesize
                                elif tx==0 and x==1:
                                    tileposx = self.tilesize
                                else:
                                    tileposx = 0
                                dsquery.WriteRaster( tileposx, tileposy, self.tilesize, self.tilesize,
                                    dsquerytile.ReadRaster(0,0,self.tilesize,self.tilesize))
                                children.append( [x, y, tz+1] )

                    self.scale_query_to_tile(dsquery, dstile, tilefilename)
                   
                    # Write a copy of tile to png/jpg
                    self.out_drv.CreateCopy(tilefilename, dstile, strict=0)

                    #if not self.options.verbose and not self.options.quiet:
                    self.progressbar( ti / float(tcount) )
                    del dstile
				
    def geo_query(self, ds, ulx, uly, lrx, lry, querysize = 0):
        """For given dataset and query in cartographic coordinates
        returns parameters for ReadRaster() in raster coordinates and
        x/y shifts (for border tiles). If the querysize is not given, the
        extent is returned in the native resolution of dataset ds."""

	# Get read bounds
        geotran = ds.GetGeoTransform()
	ul = self.WGS84toMER(uly, ulx)
	lr = self.WGS84toMER(lry, lrx)
        rx= int((ul[0] - geotran[0]) / geotran[1] + 0.001)
        ry= int((ul[1] - geotran[3]) / geotran[5] + 0.001)
        rxsize= int((lr[0] - ul[0]) / geotran[1] + 0.5)
        rysize= int((lr[1] - ul[1]) / geotran[5] + 0.5)

        if not querysize:
            wxsize, wysize = rxsize, rysize
        else:
            wxsize, wysize = querysize, querysize

        # Coordinates should not go out of the bounds of the raster
        wx = 0
        if rx < 0:
            rxshift = abs(rx)
            wx = int( wxsize * (float(rxshift) / rxsize) )
            wxsize = wxsize - wx
            rxsize = rxsize - int( rxsize * (float(rxshift) / rxsize) )
            rx = 0
        if rx+rxsize > ds.RasterXSize:
            wxsize = int( wxsize * (float(ds.RasterXSize - rx) / rxsize) )
            rxsize = ds.RasterXSize - rx

        wy = 0
        if ry < 0:
            ryshift = abs(ry)
            wy = int( wysize * (float(ryshift) / rysize) )
            wysize = wysize - wy
            rysize = rysize - int( rysize * (float(ryshift) / rysize) )
            ry = 0
        if ry+rysize > ds.RasterYSize:
            wysize = int( wysize * (float(ds.RasterYSize - ry) / rysize) )
            rysize = ds.RasterYSize - ry

        return (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)

    def scale_query_to_tile(self, dsquery, dstile, tilefilename=''):
        """Scales down query dataset to the tile dataset"""

        querysize = dsquery.RasterXSize
        tilesize = dstile.RasterXSize
        tilebands = dstile.RasterCount

        # Function: gdal.RegenerateOverview()
        for i in range(1,tilebands+1):
            # Black border around NODATA
            #if i != 4:
            #   dsquery.GetRasterBand(i).SetNoDataValue(0)
            res = gdal.RegenerateOverview( dsquery.GetRasterBand(i),
                dstile.GetRasterBand(i), 'average' )
            if res != 0:
                self.error("RegenerateOverview() failed on %s, error %d" % (tilefilename, res))

    def generate_xmldata(self):
        """Generation world wind tiles xml data"""

        print("Generating wwt xml data:")

        f = open(os.path.join(self.output, 'szcut.xml'), 'w')
        f.write( self.xmldata())
        f.close()
        print("---- done.")

    def xmldata(self):
        """templet of world wind tiles xml data"""
        
        args = {}
        args['name'] = self.output
        args['format'] = self.tileext
        args['count'] = str(self.tmaxz - self.tminz + 1)
        args['topnum'] = str(self.tminz)
        args['swlat'] = self.ominy
        args['swlon'] = self.ominx
        args['nelat'] = self.omaxy
        args['nelon'] = self.omaxx
        args['tilesize'] = str(self.tilesize)
        args['lv0size'] = self.tools.GetLztsd(self.tminz)

        s = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Layer layerType="TiledImageLayer" version="1">
    <DisplayName>%(name)s</DisplayName>
    <DatasetName>%(name)s</DatasetName>
    <DataCacheName>%(name)s</DataCacheName>
    <Service serviceName="Offline">
            </Service>
    <FormatSuffix>.%(format)s</FormatSuffix>
    <NumLevels count="%(count)s" numEmpty="%(topnum)s"/>
    <Sector>
        <SouthWest>
          <LatLon latitude="%(swlat).8f" longitude="%(swlon).8f" units="degrees"/>
        </SouthWest>
        <NorthEast>
          <LatLon latitude="%(nelat).8f" longitude="%(nelon).8f" units="degrees"/>
        </NorthEast>
    </Sector>
    <TileOrigin>
        <LatLon latitude="-90.0" longitude="-180.0" units="degrees"/>
    </TileOrigin>
    <TileSize>
        <Dimension height="%(tilesize)s" width="%(tilesize)s"/>
    </TileSize>
    <LevelZeroTileDelta>
        <LatLon latitude="%(lv0size).8f" longitude="%(lv0size).8f" units="degrees"/>
    </LevelZeroTileDelta>
    <AvailableImageFormats>
        <ImageFormat>image/%(format)s</ImageFormat>
    </AvailableImageFormats>
</Layer>
""" % args

        return s
                
    def WGS84toUTM(self, b, l):
	if self.ct84toutm:
	    return self.ct84toutm.TransformPoint(l, b)
	else:
	    raise Exception("coordinate transform error!")

    def WGS84toMER(self, b, l):
	if self.ct84tomer:
	    return self.ct84tomer.TransformPoint(l, b)
	else:
	    raise Exception("coordinate transform error!")
	
    def UTMtoWGS84(self, x, y):
	if self.ctutmto84:
	    return self.ctutmto84.TransformPoint(x, y)[:2]
	else:
	    raise Exception("coordinate transform error!")

    def MERtoWGS84(self, x, y):
        if self.ctmerto84:
            return self.ctmerto84.TransformPoint(x, y)[:2]
        else:
	    raise Exception("coordinate transform error!")
		
def main():
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    #argv = ['gdal2tiles.py', '123.tif', '7', '0']
    if argv:
        tif2wwtiles = Tif2WWTiles( argv[1:] )
        tif2wwtiles.process()

if __name__=='__main__':
    main()
