#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Z. Yunjun, H. Fattahi, Antonio Valentino, 2013   #
############################################################


import os
import sys
import warnings

from mintpy.utils.arg_utils import create_argument_parser

# from mintpy.objects import RAMP_LIST
# copy as below to avoid importing the non-empty mintpy.objects.__init__.py
RAMP_LIST = [
    'linear',
    'linear_range',
    'linear_azimuth',
    'quadratic',
    'quadratic_range',
    'quadratic_azimuth',
]

# key configuration parameter name
config_keys = [
    'mintpy.deramp',
    'mintpy.deramp.maskFile',
]


###########################################################################################
EXAMPLE = """example:
  remove_ramp.py  timeseries.h5      -m maskTempCoh.h5
  remove_ramp.py  ifgramStack.h5     -m maskTempCoh.h5  -d unwrapPhase_bridging
  remove_ramp.py  090214_101120.unw  -m maskTempCoh.h5  -s quadratic
"""


def create_parser(subparsers=None):
    synopsis = 'Remove 2D ramp(s) from the input file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='File for ramp removal')
    parser.add_argument('-m', '--mask', dest='mask_file', default='maskTempCoh.h5',
                        help='mask for pixels used in ramp estimation\n'
                             'default - maskTempCoh.h5\n'
                             'no - use the whole area')
    parser.add_argument('-s', dest='surface_type', default='linear', choices=RAMP_LIST,
                        help='type of surface/ramp to remove, linear by default')
    parser.add_argument('-d','--dset', dest='dset',
                        help='dataset name to be derampped in ifgramStack file\n'
                             'e.g.: unwrapPhase\n'
                             '      unwrapPhase_bridging')

    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('--save-ramp-coeff', dest='save_ramp_coeff', action='store_true',
                        help='Save the estimated ramp coefficients into text file.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if:\n'
                             '1) output file already exists, readable '
                             'and newer than input file\n'
                             '2) all configuration parameters are the same.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: coupled options (--update requires --outfile)
    if inps.update_mode and not inps.outfile:
        inps.update_mode = False
        warnings.warn('--update is given but the required --outfile is NOT, ignore --update and continue.')

    return inps


###########################################################################################
def run_or_skip(inps, extra_meta):
    from mintpy.utils import readfile

    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print(f'1) output file {inps.outfile} NOT found.')
    else:
        print(f'1) output file {inps.outfile} already exists.')
        infiles = [inps.file]
        if inps.mask_file:
            infiles.append(inps.mask_file)
        ti = max(os.path.getmtime(i) for i in infiles)
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print(f'2) output file is NOT newer than input file: {infiles}.')
        else:
            print(f'2) output file is newer than input file: {infiles}.')

    # check configuration
    if flag == 'skip':
        atr = readfile.read_attribute(inps.outfile)
        if any(str(extra_meta[key]) != atr.get(key, 'None') for key in config_keys):
            flag = 'run'
            print(f'3) NOT all key configuration parameters are the same:{config_keys}')
        else:
            print(f'3) all key configuration parameters are the same:{config_keys}')

    # result
    print(f'run or skip: {flag}.')
    return flag


###########################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.utils import utils1 as ut

    # run or skip
    extra_meta = {'mintpy.deramp' : inps.surface_type,
                  'mintpy.deramp.maskFile' : inps.mask_file}
    if inps.update_mode and run_or_skip(inps, extra_meta) == 'skip':
        return

    # run
    ut.run_deramp(
        inps.file,
        ramp_type=inps.surface_type,
        mask_file=inps.mask_file,
        out_file=inps.outfile,
        datasetName=inps.dset,
        save_ramp_coeff=inps.save_ramp_coeff,
        extra_meta=extra_meta,
    )


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
