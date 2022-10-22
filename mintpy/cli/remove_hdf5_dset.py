#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Antonio Valentino, Aug 2018        #
############################################################

import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

###########################################################################################
EXAMPLE = """Example:
  remove_hdf5_dset.py  ifgramStack.h5  unwrapPhase_phaseClosure
  remove_hdf5_dset.py  ifgramStack.h5  unwrapPhase_phaseClosure  unwrapPhase_bridging
  remove_hdf5_dset.py  velocity.h5     velocityStd
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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: input file extension
    if os.path.splitext(inps.file)[1] not in ['.h5', '.he5']:
        raise ValueError(f'input file is NOT HDF5: {inps.file}')

    return inps


###########################################################################################
def run_remove_hdf5_dset(fname, ds_names):
    """Remove a dataset from the given HDF5 file.

    Parameters: fname    - str, path to the HDF5 data file
                ds_names - list(str), name of the HDF5 dataset to be removed.
    """
    import h5py

    from mintpy.utils import writefile

    # grab exiting dataset list
    with h5py.File(fname, 'r') as f:
        dset_list = list(f.keys())

    # check if given dataset exists
    if any(i not in dset_list for i in ds_names):
        msg = f'input dataset ({ds_names}) do not exist!'
        msg += f'\nAvailable datasets: {dset_list}'
        raise ValueError(msg)

    # update file
    writefile.remove_hdf5_dataset(fname, ds_names, print_msg=True)
    print('Done.')

    return


###########################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # run
    run_remove_hdf5_dset(inps.file, inps.dset)


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
