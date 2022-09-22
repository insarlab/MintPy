#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

#########################################################################
EXAMPLE = """example:
  prep_gmtsar.py StHelensEnvDT156.txt
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare GMTSAR metadata files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('template_file', type=str, help='MintPy template file for GMTSAR products.')
    parser.add_argument('--mintpy-dir', dest='mintpy_dir', default='./',
                        help='MintPy directory (default: %(default)s).')
    parser.add_argument('--force', dest='update_mode', action='store_false',
                        help='Force to overwrite all .rsc metadata files.')
    return parser


def cmd_line_parse(iargs = None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: input paths (use full path and translate user symbol)
    inps.template_file = os.path.abspath(inps.template_file)
    inps.mintpy_dir = os.path.expanduser(inps.mintpy_dir)
    inps.mintpy_dir = os.path.abspath(inps.mintpy_dir)

    return inps


#########################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_gmtsar import prep_gmtsar

    # run
    prep_gmtsar(inps)


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
