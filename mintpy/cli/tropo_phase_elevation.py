#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import argparse
import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################################
REFERENCE = """reference:
  Doin, M. P., C. Lasserre, G. Peltzer, O. Cavalie, and C. Doubre (2009), Corrections of
    stratified tropospheric delays in SAR interferometry: Validation with global atmospheric
    models, J App. Geophy., 69(1), 35-50, doi:10.1016/j.jappgeo.2009.03.010.
"""

EXAMPLE = """example:
  tropo_phase_elevation.py  timeseries_demErr.h5      -g inputs/geometryRadar.h5  -m maskTempCoh.h5
  tropo_phase_elevation.py  geo_timeseries_demErr.h5  -g geo_geometryRadar.h5     -m geo_maskTempCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Correct Topo-correlated Stratified tropospheric delay'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file', help='time-series file to be corrected')
    parser.add_argument('-g', '--geometry', dest='geom_file', required=True,
                        help='DEM file used for correlation calculation.')
    parser.add_argument('-m', '--mask', dest='mask_file', required=True,
                        help='mask file for pixels used for correlation calculation')

    parser.add_argument('-t', '--threshold', type=float, default=0.,
                        help='correlation threshold to apply phase correction.\n'
                             'if not set, all dates will be corrected.')
    parser.add_argument('-l', '--looks', dest='num_multilook', type=int, default=8,
                        help='number of looks applied to data for empirical estimation (default: %(default)s).')

    parser.add_argument('--poly-order', '-p', dest='poly_order', type=int, default=1, choices=[1, 2, 3],
                        help='polynomial order of phase-height correlation (default: %(default)s).')
    parser.add_argument('-o', '--outfile', help='output corrected timeseries file name')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: -t / --threshold option (must be within [0,1])
    if inps.threshold and (not 0.0 <= inps.threshold <= 1.0):
        msg = f'correction threshold {inps.threshold} is NOT within [0.0, 1.0]'
        raise argparse.ArgumentTypeError(msg)

    return inps


############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.tropo_phase_elevation import run_tropo_phase_elevation

    # run
    run_tropo_phase_elevation(inps)


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
