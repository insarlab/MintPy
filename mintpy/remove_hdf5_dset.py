############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################


import h5py
from mintpy.utils import writefile


def run_remove_hdf5_dset(fname, ds_names):
    """Remove a dataset from the given HDF5 file.

    Parameters: fname    - str, path to the HDF5 data file
                ds_names - list(str), name of the HDF5 dataset to be removed.
    """
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
