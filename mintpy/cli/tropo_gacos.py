############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
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
  tropo_gacos.py -f timeseries.h5 -g inputs/geometryRadar.h5 --dir ./GACOS
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
    parser.add_argument('--dir','--GACOS-dir', dest='GACOS_dir', default='./GACOS',
                        help='directory to downloaded GACOS delays data (default: %(default)s).')
    parser.add_argument('-o', dest='cor_dis_file',
                        help='Output file name for trospheric corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    from mintpy.utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.GACOS_dir = os.path.abspath(inps.GACOS_dir)
    print('Use GACOS products at directory:', inps.GACOS_dir)

    # check input files - existance
    for key in ['dis_file', 'geom_file']:
        fname = vars(inps)[key]
        if fname and not os.path.isfile(fname):
            raise FileNotFoundError('input file not exist: {}'.format(fname))

    # check input files - processors & coordinates
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    coord1 = 'geo' if 'Y_FIRST' in atr1.keys() else 'radar'
    coord2 = 'geo' if 'Y_FIRST' in atr2.keys() else 'radar'
    proc = atr1.get('PROCESSOR', 'isce')

    if coord1 == 'radar' and proc in ['gamma', 'roipac']:
        msg = 'Radar-coded file from {} is NOT supported!'.format(proc)
        msg += '\n    Try to geocode the time-series and geometry files and re-run with them instead.'
        raise ValueError(msg)

    if coord1 != coord2:
        n = max(len(os.path.basename(i)) for i in [inps.dis_file, inps.geom_file])
        msg = 'Input time-series and geometry file are NOT in the same coordinate!'
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.dis_file),  n=n, c=coord1)
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.geom_file), n=n, c=coord2)
        raise ValueError(msg)

    # default output filenames
    inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), 'GACOS.h5')
    if not inps.cor_dis_file:
        inps.cor_dis_file = inps.dis_file.split('.')[0] + '_GACOS.h5'

    return inps


############################################################################
def main(iargs=None):
    from mintpy.utils import readfile
    from mintpy.tropo_gacos import calculate_delay_timeseries, correct_timeseries, correct_single_ifgram

    inps = cmd_line_parse(iargs)

    # calculate tropo delay and savee to h5 file
    calculate_delay_timeseries(
        tropo_file=inps.tropo_file,
        dis_file=inps.dis_file,
        geom_file=inps.geom_file,
        GACOS_dir=inps.GACOS_dir)

    # correct tropo delay from dis time-series
    ftype = readfile.read_attribute(inps.dis_file)['FILE_TYPE']
    if ftype == 'timeseries':
        correct_timeseries(
            dis_file=inps.dis_file,
            tropo_file=inps.tropo_file,
            cor_dis_file=inps.cor_dis_file)

    elif ftype == '.unw':
        correct_single_ifgram(
            dis_file=inps.dis_file,
            tropo_file=inps.tropo_file,
            cor_dis_file=inps.cor_dis_file)
    else:
        print('input file {} is not timeseries nor .unw, correction is not supported yet.'.format(ftype))


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
