#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Modified from load_data.py written by Heresh Fattahi.
#

import os
import sys

import h5py
import numpy as np

import pysar._readfile as readfile


#################################  Usage  ####################################
def usage():
    print '''
******************************************************************************
  Calculate temporal average/mean of multi-temporal datasets.

  Usage:
      mean_temporal.py multi_temporal_file [output_filename]

  Example:
      mean_temporal.py Coherence.h5 average_spatial_coherence.h5

******************************************************************************
    '''


#############################  Main Function  ################################
def main(argv):

    try:
        File = argv[0]
        atr  = readfile.read_attributes(File)
    except: usage(); sys.exit(1)

    try:    outName = argv[1]
    except: outName = 'tempMean_'+File
    #except: outName = 'average_spatial_coherence.h5'

    print '\n*************** Temporal Average ******************'

    ##### Input File Info
    k = atr['FILE_TYPE']
    print 'Input file is '+k
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    h5file = h5py.File(File)
    epochList = h5file[k].keys()
    epochList = sorted(epochList)
    print 'number of epoch: '+str(len(epochList))

    ##### Calculation
    dMean = np.zeros((length,width))
    print 'calculating ...'
    for epoch in epochList:
        print epoch
        if k in ['interferogram','coherence','wrapped']:
            d = h5file[k][epoch].get(epoch)[:]
        elif k in ['timeseries']:
            d = h5file[k].get(epoch)[:]
        else: print k+' type is not supported currently.'; sys.exit(1)
        dMean += d
    dMean /= len(epochList)
    del d
    h5file.close()

    ##### Output
    print 'writing >>> '+outName
    h5mean = h5py.File(outName,'w')
    group  = h5mean.create_group('mask')
    dset = group.create_dataset(os.path.basename('mask'), data=dMean, compression='gzip')
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5mean.close()


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
