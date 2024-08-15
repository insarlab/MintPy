#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import glob
import os
import sys

from mintpy.utils import arg_utils

####################################################################################
EXAMPLE = """example:
  prep_fringe.py -u './PS_DS/unwrap/*.unw' -c ./PS_DS/tcorr_ds_ps.bin -g ./geometry -m '../reference/IW*.xml' -b ../baselines -o ./mintpy

  cd ~/data/SanAndreasSenDT42/fringe
  prep_fringe.py

  ## example commands after prep_fringe.py
  reference_point.py timeseries.h5 -y 500 -x 1150
  generate_mask.py temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5
  remove_ramp.py timeseries_ERA5.h5 -m maskTempCoh.h5 -s linear
  dem_error.py timeseries_ERA5_ramp.h5 -g inputs/geometryRadar.h5
  timeseries2velocity.py timeseries_ERA5_ramp_demErr.h5
  geocode.py velocity.h5 -l inputs/geometryRadar.h5
"""

def create_parser(subparsers=None):
    """Command Line Parser"""
    synopsis = "Prepare FRInGE products for MintPy"
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-u', '--unw-file', dest='unwFile', type=str, default='./PS_DS/unwrap/*.unw',
                        help='path pattern of unwrapped interferograms (default: %(default)s).')
    parser.add_argument('-c', '--coh-file', dest='cohFile', type=str, default='./PS_DS/tcorr_ds_ps.bin',
                        help='temporal coherence file (default: %(default)s).')
    parser.add_argument('--ps-mask', dest='psMaskFile', type=str, default='./ampDispersion/ps_pixels',
                        help='PS pixels file (default: %(default)s).')
    parser.add_argument('-g', '--geom-dir', dest='geomDir', type=str, default='./geometry',
                        help='FRInGE geometry directory (default: %(default)s).\n'
                             'This is used to grab 1) bounding box\n'
                             '                 AND 2) geometry source directory where the binary files are.')
    parser.add_argument('-w','--water-mask', dest='water_mask_file', type=str, default='./geometry/waterMask.rdr',
                        help='path to water mask file (default: %(default)s or None if not exist).')

    parser.add_argument('-m', '--meta-file', dest='metaFile', type=str, default='../reference/IW*.xml',
                        help='metadata file (default: %(default)s).\n'
                             'e.g.: ./reference/IW1.xml        for ISCE/topsStack OR\n'
                             '      ./referenceShelve/data.dat for ISCE/stripmapStack')
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default='../baselines',
                        help='baseline directory (default: %(default)s).')

    parser.add_argument('-o', '--out-dir', dest='outDir', type=str, default='./mintpy',
                        help='output directory (default: %(default)s).')

    parser.add_argument('-r','--range', dest='lks_x', type=int, default=1,
                        help='number of looks in range direction, for multilooking applied after fringe processing.\n'
                             'Only impacts metadata. (default: %(default)s).')
    parser.add_argument('-a','--azimuth', dest='lks_y', type=int, default=1,
                        help='number of looks in azimuth direction, for multilooking applied after fringe processing.\n'
                             'Only impacts metadata. (default: %(default)s).')

    parser.add_argument('--geom-only', action='store_true',
                        help='Only create the geometry file (useful for geocoding a watermask).')

    parser = arg_utils.add_subset_argument(parser, geo=False)

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: --meta-file option (if contains wildcard)
    inps.metaFile = sorted(glob.glob(inps.metaFile))[0]

    # check: --water-mask option (if default/input file exists)
    if not os.path.isfile(inps.water_mask_file):
        msg = f'WARNING: default/input water mask file ({inps.water_mask_file}) NOT found, '
        msg += 'continue without it (by setting water_mask_file to None).'
        print(msg)
        inps.water_mask_file = None

    return inps


####################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_fringe import load_fringe

    # run
    load_fringe(inps)


####################################################################################
if __name__=="__main__":
    main(sys.argv[1:])
