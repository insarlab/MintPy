#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################
# Modified from _gmt.py, GIANT v1.0, Caltech.


import sys
import numpy as np
from scipy.io import netcdf
from mintpy.utils import readfile, plot as pp
from mintpy.utils.arg_utils import create_argument_parser


####################################################################################
EXAMPLE = """example:
  save_gmt.py  geo_velocity.h5
  save_gmt.py  geo_timeseries.h5  20071031
  save_gmt.py  geo_timeseries.h5
  save_gmt.py  geo_filt_100608-101024-sim_HDR_16rlks_c10.unw
  save_gmt.py  gsi10m.dem
"""


def create_parser(subparsers=None):
    synopsis = 'Export geocoded file to GMT grd file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    if 'Y_FIRST' not in atr.keys():
        raise Exception('ERROR: input file is not geocoded.')

    if not inps.dset and atr['FILE_TYPE'] in ['timeseries', 'ifgramStack']:
        raise Exception("No dataset input, it's required for {} file".format(atr['FILE_TYPE']))
    return inps


####################################################################################
def write_gmt_simple(lons, lats, z, fname, title='default', name='z', scale=1.0, offset=0, units='meters'):
    """Writes a simple GMT grd file with one array.
    This is based on the gdal2grd.py script found at:
        http://http://www.vso.cape.com/~nhv/files/python/gdal/gdal2grd.py

    Parameters: lons : 1D Array of lon values
                lats : 1D Array of lat values
                z : 2D slice to be saved
                fname : Output file name
    Kwargs:     title : Title for the grd file
                name : Name of the field in the grd file
                scale : Scale value in the grd file
                offset : Offset value in the grd file
    Returns:    fname
    """

    fid = netcdf.netcdf_file(fname, 'w')

    # Create a dimension variable
    fid.createDimension('side', 2)
    fid.createDimension('xysize', np.prod(z.shape))

    # Range variables
    fid.createVariable('x_range', 'd', ('side',))
    fid.variables['x_range'].units = 'degrees'

    fid.createVariable('y_range', 'd', ('side',))
    fid.variables['y_range'].units = 'degrees'

    fid.createVariable('z_range', 'd', ('side',))
    fid.variables['z_range'].units = units

    # Spacing
    fid.createVariable('spacing', 'd', ('side',))
    fid.createVariable('dimension', 'i4', ('side',))

    fid.createVariable('z', 'f', ('xysize',))
    fid.variables['z'].long_name = name
    fid.variables['z'].scale_factor = scale
    fid.variables['z'].add_offset = offset
    fid.variables['z'].node_offset = 0

    fid.title = title
    fid.source = 'MintPy'

    # Filling in the actual data
    fid.variables['x_range'][0] = lons[0]
    fid.variables['x_range'][1] = lons[-1]
    fid.variables['spacing'][0] = lons[1]-lons[0]

    fid.variables['y_range'][0] = lats[0]
    fid.variables['y_range'][1] = lats[-1]
    fid.variables['spacing'][1] = lats[1]-lats[0]

    # Range
    fid.variables['z_range'][0] = np.nanmin(z)
    fid.variables['z_range'][1] = np.nanmax(z)

    fid.variables['dimension'][:] = z.shape[::-1]
    fid.variables['z'][:] = np.flipud(z).flatten()
    fid.close()
    return fname
    ############################################################
    # Program is part of GIAnT v1.0                            #
    # Copyright 2012, by the California Institute of Technology#
    # Contact: earthdef@gps.caltech.edu                        #
    ############################################################


def get_geo_lat_lon(atr):
    X_FIRST = float(atr['X_FIRST'])
    Y_FIRST = float(atr['Y_FIRST'])
    X_STEP = float(atr['X_STEP'])
    Y_STEP = float(atr['Y_STEP'])
    W = int(atr['WIDTH'])
    L = int(atr['LENGTH'])
    Y_END = Y_FIRST + L*Y_STEP
    X_END = X_FIRST + W*X_STEP

    X = np.linspace(X_FIRST, X_END, W)
    Y = np.linspace(Y_FIRST, Y_END, L)
    #XI,YI = np.meshgrid(X,Y)

    return Y, X


def write_grd_file(data, atr, fname_out=None):
    """Write GMT .grd file for input data matrix, using giant._gmt module.
    Inputs:
        data - 2D np.array in int/float, data matrix to write
        atr  - dict, attributes of input data matrix
        fname_out - string, output file name
    Output:
        fname_out - string, output file name
    """
    # Get 1D array of lats and lons
    lats, lons = get_geo_lat_lon(atr)

    # writing
    print('writing >>> '+fname_out)
    write_gmt_simple(lons, np.flipud(lats), np.flipud(data), fname_out,
                     title='default', name=atr['FILE_TYPE'],
                     scale=1.0, offset=0, units=atr['UNIT'])
    return fname_out


####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # Read data
    data, atr = readfile.read(inps.file, datasetName=inps.dset) 

    # 2. Write GMT .grd file
    if not inps.outfile:
        outbase = pp.auto_figure_title(inps.file, datasetNames=inps.dset, inps_dict=vars(inps))
        inps.outfile = '{}.grd'.format(outbase)

    write_grd_file(data, atr, inps.outfile)

    print('Done.')
    return


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
