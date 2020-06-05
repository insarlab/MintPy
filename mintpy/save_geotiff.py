#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, 2013                           #
############################################################

import gdal, ogr, os, osr
import numpy as np
import argparse
from mintpy.utils import readfile

EXAMPLE = """example:
  save_tif.py geo/geo_velocity.h5 
  save_tif.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120
  save_tif.py geo/geo_timeseries_ERA5_demErr.h5 -d 20200505_20200517
  save_tif.py geo/geo_ifgramStack.h5 -d 20101120_20110220
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate GeoTiff file from h5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dataset', dest='dset',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .tif')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    # Check 1: file in geo coord
    if 'X_FIRST' not in atr.keys():
        raise Exception('ERROR: Input file is not geocoded.')

    return inps

def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def main(iargs=None):
    inps = cmd_line_parse(iargs)

    if not inps.dset:
        array, attr = readfile.read(inps.file)
    array, attr = readfile.read(inps.file, datasetName=inps.dset)

    rasterOrigin = (float(attr['X_FIRST']),float(attr['Y_FIRST']))
    pixelWidth = float(attr['X_STEP'])
    pixelHeight = float(attr['Y_STEP'])

    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.file)[0]
    inps.outfile = os.path.abspath(inps.outfile)

    reversed_arr = array[::-1]
    array2raster(inps.outfile,rasterOrigin,pixelWidth,pixelHeight,reversed_arr) # convert array to raster

if __name__ == "__main__":
    main()
