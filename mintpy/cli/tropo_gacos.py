#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Sara Mirzaee, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################################
REFERENCE = """references:
  Yu, C., Li, Z., Penna, N. T., & Crippa, P. (2018). Generic atmospheric correction
    model for Interferometric Synthetic Aperture Radar observations. Journal of
    Geophysical Research: Solid Earth, 123(10), 9202-9222, doi:10.1029/2017JB015305
"""

DIR_DEMO = """--dir ./GACOS
  20060624.ztd
  20060624.ztd.rsc
  20061225.ztd
  20061225.ztd.rsc
  ...
  OR
  20060624.ztd.tif
  20061225.ztd.tif
  ...
"""

EXAMPLE = """example:
  tropo_gacos.py -f timeseries.h5         -g inputs/geometryRadar.h5  --dir ./GACOS
  tropo_gacos.py -f geo/geo_timeseries.h5 -g geo/geo_geometryRadar.h5 --dir ./GACOS
"""


def create_parser(subparsers=None):
    synopsis = 'Tropospheric correction using GACOS (http://www.gacos.net) delays'
    epilog = REFERENCE + '\n' + DIR_DEMO + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-f', '--file', dest='dis_file', required=True,
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g', '--geom', dest='geom_file', required=True,
                        help='geometry file.')
    parser.add_argument('--dir','--gacos-dir', dest='gacos_dir', default='./GACOS',
                        help='directory to downloaded GACOS delays data (default: %(default)s).')
    parser.add_argument('-o', dest='cor_dis_file',
                        help='Output file name for trospheric corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: --gacos-dir (use absolute path)
    inps.gacos_dir = os.path.abspath(inps.gacos_dir)
    print('Use GACOS products at directory:', inps.gacos_dir)

    # check: exsitence of input files
    for fname in [inps.dis_file, inps.geom_file]:
        if fname and not os.path.isfile(fname):
            raise FileNotFoundError(f'input file not exist: {fname}')

    # check: processors & coordinates of input files
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    coord1 = 'geo' if 'Y_FIRST' in atr1.keys() else 'radar'
    coord2 = 'geo' if 'Y_FIRST' in atr2.keys() else 'radar'
    proc = atr1.get('PROCESSOR', 'isce')

    # check: radar-coord product from gamma and roipac is not supported
    if coord1 == 'radar' and proc in ['gamma', 'roipac']:
        msg = f'Radar-coded file from {proc} is NOT supported!'
        msg += '\n    Try to geocode the time-series and geometry files and re-run with them instead.'
        raise ValueError(msg)

    # check: coordinate system must be consistent btw. displacement and geometry files
    if coord1 != coord2:
        n = max(len(os.path.basename(i)) for i in [inps.dis_file, inps.geom_file])
        msg = 'Input time-series and geometry file are NOT in the same coordinate!'
        msg += f'\n    file {os.path.basename(inps.dis_file):<{n}} coordinate: {coord1}'
        msg += f'\n    file {os.path.basename(inps.geom_file):<{n}} coordinate: {coord2}'
        raise ValueError(msg)

    # default: output tropo and corrected displacement file names
    inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), 'GACOS.h5')
    if not inps.cor_dis_file:
        inps.cor_dis_file = inps.dis_file.split('.')[0] + '_GACOS.h5'

    return inps


############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.tropo_gacos import run_tropo_gacos

    # run
    run_tropo_gacos(inps)


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
