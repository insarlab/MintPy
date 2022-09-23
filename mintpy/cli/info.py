#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################
EXAMPLE = """example:
  info.py timeseries.h5
  info.py velocity.h5
  info.py ifgramStack.h5

  # Display dataset
  info.py timeseries.py --dset date
  info.py timeseries.py --dset bperp

  # Time / Date Info
  info.py ifgramStack.h5 --date                 #print date1_date2 info for all  interferograms
  info.py timeseries.h5  --num                  #print date list of timeseries with its number
  info.py ifgramStack.h5 --date --show kept     #print date1_date2 info for kept interferograms
  info.py ifgramStack.h5 --date --show dropped  #print date1_date2 info for dropped/excluded interferograms
  info.py LS-PARAMS.h5   --date > date_list.txt #print date list of timeseries and save it to txt file.
  info.py S1_IW12_128_0593_0597_20141213_20180619.h5 --date

  # save date1_date2 info of interferograms to a text file
  info.py ifgramStack.h5 --date --show kept > date12_list.txt

  # Slice / Dataset Info
  info.py timeseries.h5                              --slice
  info.py inputs/ifgramStack.h5                      --slice
  info.py S1_IW12_128_0593_0597_20141213_20180619.h5 --slice
  info.py LS-PARAMS.h5                               --slice
"""


def create_parser(subparsers=None):
    """Create command line parser."""
    synopsis = 'Display Metadata / Structure information of ANY File'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='File to check')
    parser.add_argument('--compact', action='store_true',
                        help='show compact info by displaying only the top 20 metadata')

    parser.add_argument('--dset', type=str, help='Show dataset')

    par_list = parser.add_argument_group('List','list date/slice info')
    par_list.add_argument('--date', dest='disp_date', action='store_true',
                          help='Show date/date12 info of input file')
    par_list.add_argument('--num', dest='disp_num', action='store_true',
                          help='Show date/date12 info with numbers')
    par_list.add_argument('--slice', dest='disp_slice', action='store_true',
                          help='Show slice list of the file')
    par_list.add_argument('--show','--show-ifgram', dest='disp_ifgram',
                          choices={'all','kept','dropped'}, default='all',
                          help='Show all / kept / dropped interferograms only. Default: all.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: input file existence
    if not os.path.isfile(inps.file):
        raise FileNotFoundError(inps.file)

    # default: --compact option and max number of metadata to show
    inps.max_meta_num = 10000
    if inps.compact:
        inps.max_meta_num = 20

    return inps


############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.info import print_info

    # run
    print_info(inps)


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
