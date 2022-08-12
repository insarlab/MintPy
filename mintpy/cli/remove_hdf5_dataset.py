############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################

import os
import sys
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
    from ..remove_hdf5_dataset import remove_hdf5_dataset
    inps = cmd_line_parse(iargs)
    remove_hdf5_dataset(inps.file, inps.dset)
    print('Done.')


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
