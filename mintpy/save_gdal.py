#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Jun 2020                       #
############################################################


import os
import sys
from osgeo import gdal, osr
from mintpy.utils import readfile, utils0 as ut, plot as pp
from mintpy.utils.arg_utils import create_argument_parser


# link: https://gdal.org/drivers/raster/index.html
GDAL_DRIVER2EXT = {
    'GTiff' : '.tif',
    'ENVI'  : '',
    'GMT'   : '.grd',
    'GRIB'  : '.grb',
    'JPEG'  : '.jpg',
    'PNG'   : '.png',
}


EXAMPLE = """example:
  save_gdal.py geo/geo_velocity.h5 
  save_gdal.py geo/geo_timeseries_ERA5_demErr.h5 -d 20200505_20200517 --of ENVI
  save_gdal.py geo/geo_ifgramStack.h5 -d unwrapPhase-20101120_20110220 --of ISCE
  save_gdal.py geo/geo_ifgramStack.h5 -d coherence-20101120_20110220 --of ISCE
"""


def create_parser(subparsers=None):
    synopsis = 'Generate GDAL raster from MintPy h5 file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed by GDAL driver')
    parser.add_argument('--of', '--out-format', '--output-format', dest='out_format', default='GTiff',
                        help='file format as defined by GDAL driver name, e.g. GTiff, ENVI, default: %(default)s\n'
                             'GDAL driver names can be found at https://gdal.org/drivers/raster/index.html')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    # Check if file in geo coord
    if 'X_FIRST' not in atr.keys():
        raise Exception('ERROR: Input file is not geocoded.')

    return inps


##############################################################################
def array2raster(array, rasterName, rasterFormat, rasterOrigin, xStep, yStep, epsg=4326):

    # transform info
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]
    transform = (originX, xStep, 0, originY, 0, yStep)

    # write
    driver = gdal.GetDriverByName(rasterFormat)
    print('initiate GDAL driver: {}'.format(driver.LongName))

    print('create raster band')
    print('raster row / column number: {}, {}'.format(rows, cols))
    print('raster transform info: {}'.format(transform))
    outRaster = driver.Create(rasterName, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(transform)

    print('write data to raster band')
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)

    print('set projection as: EPSG {}'.format(epsg))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    print('finished writing to {}'.format(rasterName))

    return rasterName


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # check reference date
    print('read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)
    if attr['FILE_TYPE'] == 'timeseries' and inps.dset and '_' in inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date = None

    # read data
    print('read data     from file: {}'.format(inps.file))
    array, attr = readfile.read(inps.file, datasetName=inps.dset)
    if attr['FILE_TYPE'] == 'timeseries' and inps.ref_date:
        array -= readfile.read(inps.file, datasetName=inps.ref_date)[0]

    # output filename
    if not inps.outfile:
        fbase = pp.auto_figure_title(inps.file,
                                     datasetNames=inps.dset,
                                     inps_dict=vars(inps))
        inps.outfile = fbase + GDAL_DRIVER2EXT.get(inps.out_format, '')
    else:
        inps.outfile = os.path.abspath(inps.outfile)

    # coordinate info
    rasterOrigin = (float(attr['X_FIRST']),float(attr['Y_FIRST']))
    xStep = float(attr['X_STEP'])
    yStep = float(attr['Y_STEP'])
    kwargs = dict(xStep=xStep, yStep=yStep)

    epsg = attr.get('EPSG', None)
    if not epsg and 'UTM_ZONE' in attr.keys():
        epsg = ut.utm_zone2epsg_code(attr['UTM_ZONE'])
    if epsg:
        kwargs['epsg'] = int(epsg)

    # convert array to raster
    array2raster(array, inps.outfile, inps.out_format, rasterOrigin, **kwargs)

    return


##############################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
