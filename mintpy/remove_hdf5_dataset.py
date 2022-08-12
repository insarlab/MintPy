############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################


import h5py
from mintpy.utils import writefile


###########################################################################################
def remove_hdf5_dataset(filename, dset):
    with h5py.File(filename, 'r') as f:
        dset_list = list(f.keys())
    if any(i not in dset_list for i in dset):
        raise ValueError(('input dataset do not exists: {}'
                          '\navailable datasets:\n{}').format(dset, dset_list))

    writefile.remove_hdf5_dataset(filename, dset, print_msg=True)
