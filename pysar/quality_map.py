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
#from scipy.sparse.csgraph import laplacian
from scipy.ndimage.filters import laplace


##############################################################################
def usage():
    print("""
***************************************************************************
   Usage:
       quality_map.py  file.h5 

   Example:
       quality_map.py  Seeded_Loaded_data.h5

***************************************************************************
    """)
#def Sudo_Correlation():

#def Derivative_Variance():


##############################################################################
def main(argv):

    try:    file=argv[0]
    except: usage();sys.exit(1)
  
    h5file=h5py.File(file,'r')
    kh5=list(h5file.keys())
    ifgramList=list(h5file['interferograms'].keys())
  
    try:    OutName=argv[1]
    except: OutName='Laplacian.h5'
  
    h5laplace=h5py.File(OutName,'w')
    group=h5laplace.create_group('interferograms')
    print('Calculating the Discrete Laplacian Transform')   
    for ifgram in  ifgramList:
        print(ifgram)
        dset=h5file['interferograms'][ifgram].get(ifgram)
        unw=dset[0:dset.shape[0],0:dset.shape[1]]
        Lunw=laplace(unw)
        g=group.create_group(ifgram)
        g.create_dataset(ifgram,data=Lunw,compression='gzip')
        for key, value in h5file['interferograms'][ifgram].attrs.items():
            g.attrs[key] = value
  
    gm = h5laplace.create_group('mask')
    mask = h5file['mask'].get('mask')
    dset = gm.create_dataset('mask', data=mask)
  
    try:
        meanCoherence = h5file['meanCoherence'].get('meanCoherence')
        gc = h5laplace.create_group('meanCoherence')
        dset = gc.create_dataset('meanCoherence', data=meanCoherence)
    except:
        print('')   
  
    print('DONE!')

##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

