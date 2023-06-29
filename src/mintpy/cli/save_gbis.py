#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

##############################################################################
REFERENCE = """references:
  Bagnardi, M., and A. Hooper (2018), Inversion of Surface Deformation Data for Rapid Estimates of
    Source Parameters and Uncertainties: A Bayesian Approach, Geochemistry, Geophysics, Geosystems,
    19, doi:10.1029/2018GC007585.
  Yunjun, Z., Amelung, F., & Aoki, Y. (2021), Imaging the hydrothermal system of Kirishima volcanic
    complex with L-band InSAR time series, Geophysical Research Letters, 48(11), e2021GL092879.
    doi:10.1029/2021GL092879
"""

EXAMPLE = """example:
  save_gbis.py velocity.h5 -g inputs/geometryGeo.h5 -o AlosDT73_20081012_20100302.mat
  save_gbis.py 20150223_20161031_msk.unw -g inputs/geometryGeo.h5 -o Alos2DT23_20150223_20161031.mat
  save_gbis.py 20150223_20161031_msk.unw -g inputs/geometryGeo.h5 --out-dir ../Model/data --ellipsoid2geoid

  #for non-velocity dataset, e.g. step or annualAmplitude
  save_gbis.py velocity.h5 step20080221    -o AlosDT73_step20080221.mat
  save_gbis.py velocity.h5 annualAmplitude -o AlosDT73_annualAmplitude.mat
"""

def create_parser(subparsers=None):
    synopsis = 'Convert MintPy product to GBIS .mat format.'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='deformation file.')
    parser.add_argument('dset', nargs='?',
                        help='date/date12 of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-g','--geometry', dest='geom_file', required=True, help='geometry file')
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file.')

    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2,
                        help='custom reference pixel in lat/lon')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                        help='do not display the figure')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('--out-dir', dest='outdir',
                        help='custom output directory, ONLY IF --output is not specified.')
    parser.add_argument('--ellipsoid2geoid', action='store_true',
                        help='Convert the height of ellipsoid to geoid using "geoidheight" module\n'+
                             'Download & install geoidheight as below:\n'+
                             'https://github.com/geodesymiami/2021_Kirishima')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.argv = iargs if iargs else sys.argv[1:]
    print('{} {}'.format(os.path.basename(__file__), ' '.join(inps.argv)))

    return inps


##############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_gbis import save_gbis

    # run
    save_gbis(inps)


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
