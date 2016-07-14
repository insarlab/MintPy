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


#################################  Usage  ####################################
def Usage():
    print '''
    ************************************************************************

       Calculate average spatial coherence from Coherence.h5

       Usage:
              mean_spatial_coherence.py coherence_file [output_filename]

       Example:
              mean_spatial_coherence.py Coherence.h5
              mean_spatial_coherence.py Coherence.h5 average_spatial_coherence.h5

    ************************************************************************
    '''


#############################  Main Function  ################################
def main(argv):

  try:
      file = argv[1]
      h5file = h5py.File(file)
      k = h5file.keys()
      if   'interferograms' in k: k[0] = 'interferograms'
      elif 'coherence'      in k: k[0] = 'coherence'
      print 'Input file is '+k[0]
  except: Usage(); sys.exit(1)

  try:    outName = argv[2]
  except: outName = 'average_spatial_coherence.h5'

  epochList = h5file[k[0]].keys(); print 'number of epoch: '+str(len(epochList))
  epochSet  = h5file[k[0]][epochList[0]].get(epochList[0])
  epochMean = np.zeros(epochSet.shape)
  print 'calculating ...'
  for epoch in epochList:
      print epoch
      epochSet  = h5file[k[0]][epoch].get(epoch)
      epochMean += epochSet[0:epochSet.shape[0],0:epochSet.shape[1]]
  epochMean /= len(epochList)

  print 'writing >>> '+outName
  h5mean = h5py.File(outName,'w')
  group  = h5mean.create_group('mask')
  dset = group.create_dataset(os.path.basename('mask'), data=epochMean, compression='gzip')
  for key,value in h5file[k[0]][epochList[0]].attrs.iteritems():   group.attrs[key] = value


  h5file.close()
  h5mean.close()


##############################################################################

if __name__ == '__main__':
  main(sys.argv[:])
