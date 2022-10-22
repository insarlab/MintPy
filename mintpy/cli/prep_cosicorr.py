#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Forrest Williams, Aug 2022    #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

#########################################################################
EXAMPLE_META_FILE = """
offset1NS.tif  20160206 20161122
offset1EW.tif  20160206 20161122
offset1SNR.tif 20160206 20161122
offset2NS.tif  20160206 20170225
offset2EW.tif  20160206 20170225
offset2SNR.tif 20160206 20170225
...            ...   ...
"""

EXAMPLE = """example:
  prep_cosicorr.py offsets/*offset.tif -m metadata.txt
  prep_cosicorr.py snr/*snr.tif        -m metadata.txt
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare attributes file for COSI-Corr pixel offset product.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='cosicorr file(s)')
    parser.add_argument('-m', '--metadata', type=str, dest='meta_file',
                        help='metadata file with date info. E.g.:'+EXAMPLE_META_FILE)
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#########################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_cosicorr import prep_cosicorr

    # run
    prep_cosicorr(inps)


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
