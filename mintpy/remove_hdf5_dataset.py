#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################


import os
import sys
import h5py
from mintpy.utils import writefile
from mintpy.utils.arg_utils import create_argument_parser


###########################################################################################
EXAMPLE = """Example:
  remove_hdf5_dataset.py  ifgramStack.h5  unwrapPhase_phaseClosure
  remove_hdf5_dataset.py  ifgramStack.h5  unwrapPhase_phaseClosure  unwrapPhase_bridging
  remove_hdf5_dataset.py  velocity.h5     velocityStd
"""


def create_parser(subparsers=None):
    synopsis = 'Remove an existing dataset from HDF5 file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='HDF5 file of interest')
    parser.add_argument('dset', type=str, nargs='+', help='dataset to be removed.')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if os.path.splitext(inps.file)[1] not in ['.h5', '.he5']:
        raise ValueError('input file is not HDF5: {}'.format(inps.file))
    return inps


###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    with h5py.File(inps.file, 'r') as f:
        dset_list = list(f.keys())
    if any(i not in dset_list for i in inps.dset):
        raise ValueError(('input dataset do not exists: {}'
                          '\navailable datasets:\n{}').format(inps.dset, dset_list))

    writefile.remove_hdf5_dataset(inps.file, inps.dset, print_msg=True)

    print('Done.')
    return


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
