#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: David Bekaert, March 2026                        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################################
REFERENCE = """references:
  Bekaert, D. et al., OPERA L4 Tropospheric Zenith Delay products
  derived from ECMWF HRES model (https://www.earthdata.nasa.gov/data/catalog/asf-opera-l4-tropo-zenith-v1-1).
"""

DIR_DEMO = """--dir ./OPERA
  Path to the directory for downloading and storing OPERA tropospheric delay
  products.  Files are automatically downloaded from ASF on demand.
"""

EXAMPLE = """example:
  tropo_opera.py -f timeseries.h5         -g inputs/geometryRadar.h5  --dir ./OPERA
  tropo_opera.py -f geo/geo_timeseries.h5 -g geo/geo_geometryRadar.h5 --dir ./OPERA
"""


def create_parser(subparsers=None):
    synopsis = 'Tropospheric correction using OPERA Zenith Tropospheric Delays from ECMWF HRES data (https://www.earthdata.nasa.gov/data/catalog/asf-opera-l4-tropo-zenith-v1-1)'
    epilog = REFERENCE + '\n' + DIR_DEMO + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)
    parser.add_argument('-f', '--file', dest='dis_file', required=True,
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g', '--geom', dest='geom_file', required=True,
                        help='geometry file.')
    parser.add_argument('--dir','--opera-dir', dest='opera_dir', default='./OPERA',
                        help='directory to store downloaded OPERA delays data (default: %(default)s).')
    parser.add_argument('-o', dest='cor_dis_file',
                        help='Output file name for tropospheric corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: --opera-dir (use absolute path)
    inps.opera_dir = os.path.abspath(inps.opera_dir)
    print('Use OPERA products at directory:', inps.opera_dir)

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
    inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), 'OPERA.h5')
    if not inps.cor_dis_file:
        inps.cor_dis_file = inps.dis_file.split('.')[0] + '_OPERA.h5'

    return inps


############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.tropo_opera import run_tropo_opera

    # run
    run_tropo_opera(inps)


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
