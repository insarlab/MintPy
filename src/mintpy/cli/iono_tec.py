#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

#####################################################################################
REFERENCE = """references:
  Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., & Aoki, Y. (2022).
    Range Geolocation Accuracy of C-/L-band SAR and its Implications for Operational
    Stack Coregistration. IEEE Trans. Geosci. Remote Sens., 60, doi:10.1109/TGRS.2022.3168509.
  Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format
    version 1.1. Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany.
"""

EXAMPLE = """example:
  iono_tec.py timeseriesRg.h5 -g inputs/geometryRadar.h5
  iono_tec.py timeseriesRg.h5 -g inputs/geometryRadar.h5 -s cod
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate ionospheric ramps using Global Iono Maps  from GNSS-based TEC products.'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('dis_file', help='displacement time-series HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g','--geomtry', dest='geom_file', type=str, required=True,
                        help='geometry file including incidence/azimuthAngle.')
    parser.add_argument('-s','--sol','--sol-code', dest='sol_code', default='jpl',
                        choices={'cod','esa','igs','jpl','upc','uqr'},
                        help='GIM solution center code (default: %(default)s).\n'
                             'https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html')
    parser.add_argument('--tec-dir', dest='tec_dir', default='${WEATHER_DIR}/IONEX',
                        help='directory of downloaded GNSS TEC data (default: %(default)s).')

    # output
    parser.add_argument('--update', dest='update_mode', action='store_true', help='Enable update mode.')
    parser.add_argument('--iono-file', dest='iono_file', help='calculated LOS iono ramp time-series file.')
    #parser.add_argument('-o', dest='cor_dis_file', help='Output file name for the corrected time-series.')

    # GIM extraction
    tec_cfg = parser.add_argument_group('GIM extraction', 'Parameters to extract TEC at point of interest from '
                                        'GIM (mainly for impact demonstration).')
    tec_cfg.add_argument('-i','--interp', dest='interp_method', default='linear3d',
                         choices={'nearest', 'linear2d', 'linear3d'},
                         help='Interpolation method to grab the GIM value at the point of interest'
                              ' (default: %(default)s).')
    tec_cfg.add_argument('--norotate', dest='rotate_tec_map', action='store_false',
                         help="Rotate TEC maps along the longitude direction to compensate the correlation\n"
                              "between the ionosphere and the Sun's position (Schaer et al. (1998).\n"
                              "For 'interp_method == linear3d' ONLY. (default: %(default)s).")
    tec_cfg.add_argument('--ratio', dest='sub_tec_ratio', type=str,
                         help='Ratio to calculate the sub-orbital TEC from the total TEC.\n'
                              'Set to "adaptive" for seasonally adaptive scaling.\n'
                              '     Based on equation (14) from Yunjun et al. (2022).\n'
                              'Set to "a value" within (0,1] for a fixed scaling\n'
                              'E.g. 0.75 for TerraSAR-X (Gisinger et al., 2021)\n'
                              '     0.90 for Sentinel-1 (Gisinger et al., 2021)\n'
                              '     0.69 for Sentinel-1 (Yunjun et al., 2022)\n')

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import utils0 as ut

    # check: --tec-dir option
    inps.tec_dir = os.path.expanduser(inps.tec_dir)
    inps.tec_dir = os.path.expandvars(inps.tec_dir)
    if not os.path.isdir(inps.tec_dir):
        print(f'WARNING: Input TEC dir "{inps.tec_dir}" does not exist!')
        inps.tec_dir = os.path.join(os.path.dirname(inps.dis_file), 'TEC')
        print(f'Use "{inps.tec_dir}" instead.')

    # check: --ratio option
    if inps.sub_tec_ratio is None:
        suffix = ''
    elif ut.is_number(inps.sub_tec_ratio):
        suffix = f'R{float(inps.sub_tec_ratio)*100:.0f}'
    elif inps.sub_tec_ratio.startswith('adap'):
        suffix = 'RA'
    else:
        raise ValueError('Input is neither a number nor startswith adap!')

    # default: input/output file paths
    inps.dis_file = os.path.abspath(inps.dis_file)
    inps.geom_file = os.path.abspath(inps.geom_file)

    if not inps.iono_file:
        geom_dir = os.path.dirname(inps.geom_file)
        inps.iono_file = os.path.join(geom_dir, f'TEC{inps.sol_code[0]}lr{suffix}.h5')

    #if not inps.cor_dis_file:
    #    dis_dir = os.path.dirname(inps.dis_file)
    #    fbase, fext = os.path.splitext(os.path.basename(inps.dis_file))
    #    suffix = os.path.splitext(os.path.basename(inps.iono_file))[0]
    #    inps.cor_dis_file = os.path.join(dis_dir, f'{fbase}_{suffix}{fext}')

    return inps


#####################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.iono_tec import run_iono_tec

    # run
    run_iono_tec(inps)


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
