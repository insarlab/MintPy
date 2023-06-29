#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Marin Govorcin, Aug 2022      #
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
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_kite import save_kite

    # run
    save_kite(inps)


#########################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
