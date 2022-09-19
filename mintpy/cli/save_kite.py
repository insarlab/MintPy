############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils import arg_utils


#########################################################################################################
EXAMPLE = """example:
  ## displacement [event-type inversion]
  # option 1: use velocity file with step estimation from timeseries2velocity.py for co-seismic displacement
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  # option 2: use time-series / ifgramStack file with date1_date2 for the transient displacement:
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  save_kite.py geo/geo_ifgramStack.h5     -d unwrapPhase-20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  ## velocity [interseismic or tensile dislocation inversion]
  # https://pyrocko.org/beat/docs/current/examples/Rectangular_tensile.html
  save_kite.py geo/geo_velocity.h5 -d velocity -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  ## import to kite
  spool outfile_name    % /do quadtree,covariance/aps and then File>Save Scene and it is ready for GROND or BEAT
"""

KITE_URL = 'https://github.com/pyrocko/kite'


def create_parser(subparsers=None):
    synopsis = f'Generate KITE ({KITE_URL}) npz and yaml from MintPy HDF5 file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset', type=str, required=True,
                        help='dataset of interest to be converted.\n'+
                             'e.g.: velocity / stepYYYYMMDD for velocity HDF5 file,\n'+
                             '      date12 in YYYYMMDD_YYYYMMDD for time-series HDF5 file,\n'+
                             '      date12 in unwrapPhase-YYYYMMDD_YYYYMMDD for ifgramStack HDF5 file.')
    parser.add_argument('-g', '--geom', dest='geom_file', type=str,
                        help='geometry file for incidence /azimuth angle and height.')
    parser.add_argument('-m', '--mask', dest='mask_file', type=str,
                        help='mask file, or run mask.py to mask the input file beforehand.')
    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        help='output filename')
    parser = arg_utils.add_subset_argument(parser)
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps


#########################################################################################################
def main(iargs=None):
    import numpy as np
    from mintpy import subset
    from mintpy.utils import ptime, readfile, attribute
    from mintpy.save_kite import mintpy2kite

    inps = cmd_line_parse(iargs)

    print('\n-------------------READ INPUTS -------------------')
    print('Read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)

    #Extract subset if defined
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), attr)

    # output filename
    if not inps.outfile:
        inps.outfile = attr['PROJECT_NAME']

    # date1/2
    if attr['FILE_TYPE'] in ['timeseries', 'HDFEOS']:
        date1, date2 = inps.dset.split('_')
        inps.dset = date2

    elif attr['FILE_TYPE'] == 'ifgramStack':
        date1, date2 = inps.dset.split('-')[1].split('_')

    else:
        # velocity, unw
        date1, date2 = ptime.yyyymmdd(attr['DATE12'].replace('_','-').split('-'))
        if inps.dset.startswith('step'):
            date1 = inps.dset.split('step')[-1]
            date2 = date1
    print('First  InSAR date: {}'.format(date1))
    print('Second InSAR date: {}'.format(date2))

    # read data
    print('Read {} from file: {}'.format(inps.dset, inps.file))
    dis, attr = readfile.read(inps.file, datasetName=inps.dset, box=inps.pix_box)

    if attr['FILE_TYPE'] == 'timeseries':
        print('Read {} from file: {}'.format(date1, inps.file))
        dis -= readfile.read(inps.file, datasetName=date1, box=inps.pix_box)[0]

    # convert radians to meters
    if attr['UNIT'] == 'radian':
        dis *= (float(attr['WAVELENGTH']) / (-4*np.pi))

    # mask data
    if inps.mask_file is not None:
        mask = readfile.read(inps.mask_file, box=inps.pix_box)[0]
        print('Set data to NaN for pixels with zero value in file: {}'.format(inps.mask_file))
        dis[mask==0] = np.nan

    # read geometry incidence / azimuth angle
    print('\nread incidence / azimuth angle from file: {}'.format(inps.geom_file))
    inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle', box=inps.pix_box)[0]
    az_angle = readfile.read(inps.geom_file, datasetName='azimuthAngle', box=inps.pix_box)[0]
    print('Mean satellite incidence angle: {0:.2f}°'.format(np.nanmean(inc_angle)))
    print('Mean satellite heading   angle: {0:.2f}°\n'.format(90 - np.nanmean(az_angle)))

    # Update attributes
    if inps.subset_lat is not None or inps.subset_x is not None:
        attr = attribute.update_attribute4subset(attr, inps.pix_box)

    # create kite container
    mintpy2kite(dis, attr, date1, date2, inc_angle, az_angle, out_file=inps.outfile)


#########################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
