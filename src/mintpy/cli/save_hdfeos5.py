#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

################################################################
TEMPLATE = get_template_content('hdfeos5')

EXAMPLE = """example:
  save_hdfeos5.py geo/geo_timeseries_ERA5_ramp_demErr.h5
  save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5 --tc temporalCoherence.h5 --asc avgSpatialCoh.h5 -m maskTempCoh.h5 -g inputs/geometryGeo.h5
  save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5 --tc temporalCoherence.h5 --asc avgSpatialCoh.h5 -m maskTempCoh.h5 -g inputs/geometryRadar.h5
"""

NOTE = """
  https://earthdata.nasa.gov/esdis/eso/standards-and-references/hdf-eos5
  https://mintpy.readthedocs.io/en/latest/hdfeos5/
"""


def create_parser(subparsers=None):
    synopsis = 'Convert MintPy timeseries product into HDF-EOS5 format'
    epilog = TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file', default='timeseries.h5', help='Time-series file')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='Template file for 1) arguments/options and 2) missing metadata')

    parser.add_argument('--tc','--temp-coh', dest='tcoh_file',
                        help='Coherence/correlation file, i.e. temporalCoherence.h5')
    parser.add_argument('--asc','--avg-spatial-coh', dest='scoh_file',
                        help='Average spatial coherence file, i.e. avgSpatialCoh.h5')
    parser.add_argument('-m', '--mask', dest='mask_file', help='Mask file')
    parser.add_argument('-g', '--geometry', dest='geom_file', help='geometry file')
    parser.add_argument('--suffix', dest='suffix', help='suffix to be appended to file name (e.g. PS).')

    parser.add_argument('--update', action='store_true',
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')
    parser.add_argument('--subset', action='store_true',
                        help='Enable subset mode, a.k.a. put suffix _N31700_N32100_E130500_E131100')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check
    meta = readfile.read_attribute(inps.ts_file)

    # default: input file paths
    ts_dir = os.path.dirname(inps.ts_file)
    if os.path.basename(inps.ts_file).startswith('geo_'):
        tcoh_file = os.path.join(ts_dir, 'geo_temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'geo_avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'geo_maskTempCoh.h5')
        geom_file = os.path.join(ts_dir, 'geo_geometryRadar.h5')
    else:
        tcoh_file = os.path.join(ts_dir, 'temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'maskTempCoh.h5')
        geom_file = os.path.join(ts_dir, 'inputs/geometry')
        geom_file += 'Geo.h5' if 'Y_FIRST' in meta.keys() else 'Radar.h5'

    inps.tcoh_file = inps.tcoh_file if inps.tcoh_file else tcoh_file
    inps.scoh_file = inps.scoh_file if inps.scoh_file else scoh_file
    inps.mask_file = inps.mask_file if inps.mask_file else mask_file
    inps.geom_file = inps.geom_file if inps.geom_file else geom_file

    # check: existence of input files
    for fname in [inps.ts_file, inps.tcoh_file, inps.scoh_file, inps.mask_file, inps.geom_file]:
        if not os.path.isfile(fname):
            raise FileNotFoundError(fname)

    # check: --subset mode in conflict with input file in radar-coordinates
    if inps.subset and 'Y_FIRST' not in meta.keys():
        raise SystemExit('ERROR: --subset mode is NOT supported for time-series in radar-coordinates!')

    # check: coordinate
    if not 'Y_FIRST' in meta.keys():
        raise ValueError(f'Input file {inps.ts_file} is NOT geocoded!')

    return inps


################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_hdfeos5 import save_hdfeos5

    # run
    save_hdfeos5(inps)


################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
