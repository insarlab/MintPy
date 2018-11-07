#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import sys
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.ndimage.filters import laplace
from pysar.utils import ptime


##############################################################################
USAGE = """
Usage: quality_map.py  stack_file

Example: quality_map.py  INPUTS/ifgramStack.h5
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

    print('calculate Laplace filter of {} based on approximate second derivatives.'.format(g_name))

    f = h5py.File(file, 'a')
    ds = f[g_name]
    if g_name_out in f.keys():
        ds_out = f[g_name_out]
    else:
        ds_out = f.create_dataset(g_name_out, shape=ds.shape, dtype=np.float32, chunks=True, compression=None)
    print('write to dataset /{}'.format(g_name_out))

    num_ifgram = ds.shape[0]
    prog_bar = ptime.progressBar(maxValue=num_ifgram)
    for i in range(num_ifgram):
        unw = ds[i, :, :]
        ds_out[i, :, :] = laplace(unw)
        prog_bar.update(i+1, suffix='{}/{}'.format(i+1, num_ifgram))
    prog_bar.close()
    f.close()
    print('finished writing to {}'.format(file))
    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
