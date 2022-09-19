############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


####################################################################################
EXAMPLE = """example:
  save_gmt.py  geo_velocity.h5
  save_gmt.py  geo_timeseries.h5  20071031
  save_gmt.py  geo_timeseries.h5
  save_gmt.py  geo_filt_100608-101024-sim_HDR_16rlks_c10.unw
  save_gmt.py  gsi10m.dem
"""


def create_parser(subparsers=None):
    synopsis = 'Export geocoded file to GMT grd file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')
    return parser


def cmd_line_parse(iargs=None):
    from mintpy.utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    if 'Y_FIRST' not in atr.keys():
        raise Exception('ERROR: input file is not geocoded.')

    if not inps.dset and atr['FILE_TYPE'] in ['timeseries', 'ifgramStack']:
        raise Exception("No dataset input, it's required for {} file".format(atr['FILE_TYPE']))
    return inps


####################################################################################
def main(iargs=None):
    from mintpy.utils import readfile, plot as pp
    from mintpy.save_gmt import write_grd_file

    inps = cmd_line_parse(iargs)

    # Read data
    data, atr = readfile.read(inps.file, datasetName=inps.dset) 

    # 2. Write GMT .grd file
    if not inps.outfile:
        outbase = pp.auto_figure_title(inps.file, datasetNames=inps.dset, inps_dict=vars(inps))
        inps.outfile = '{}.grd'.format(outbase)

    write_grd_file(data, atr, inps.outfile)

    print('Done.')


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
