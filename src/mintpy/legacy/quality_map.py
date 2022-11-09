#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import laplace

from mintpy.utils import ptime

##############################################################################
USAGE = """
Usage: quality_map.py  stack_file

Example: quality_map.py  inputs/ifgramStack.h5
"""


def usage():
    print(USAGE)
    return


##############################################################################
def main(argv):

    try:
        file = argv[0]
    except:
        usage()
        sys.exit(1)

    g_name = 'unwrapPhase'
    g_name_out = 'unwrapPhase_laplace'

    print(f'calculate Laplace filter of {g_name} based on approximate second derivatives.')

    f = h5py.File(file, 'a')
    ds = f[g_name]
    if g_name_out in f.keys():
        ds_out = f[g_name_out]
    else:
        ds_out = f.create_dataset(g_name_out, shape=ds.shape, dtype=np.float32, chunks=True, compression=None)
    print(f'write to dataset /{g_name_out}')

    num_ifgram = ds.shape[0]
    prog_bar = ptime.progressBar(maxValue=num_ifgram)
    for i in range(num_ifgram):
        unw = ds[i, :, :]
        ds_out[i, :, :] = laplace(unw)
        prog_bar.update(i+1, suffix=f'{i+1}/{num_ifgram}')
    prog_bar.close()
    f.close()
    print(f'finished writing to {file}')
    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
