#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

def Usage():
  print '''
***************************************************************
***************************************************************
Generates the sum of two input files.
 
   Usage:
          add.py file1.h5 file2.h5

   output is file1_plus_file2.h5  (file1 + file2)
   
   example:
           
          add.py velocity_masked.h5 velocity_demCor_masked.h5
          add.py velocity_demCor_masked.h5 velocity_demCor_tropCor_masked.h5
          add.py timeseries_demCor.h5 timeseries_demCor_tropCor.h5
          add.py timeseries.h5 timeseries_demCor.h5
          add.py interferograms.h5 interferograms2.h5

***************************************************************
***************************************************************
'''

def main(argv):

  try:
    File1=sys.argv[1]
    File2=sys.argv[2]
  except:
    Usage();sys.exit(1)


  h5file1=h5py.File(File1,'r')
  k1=h5file1.keys()
  h5file2=h5py.File(File2,'r')
  k2=h5file2.keys()

  if k1[0]!=k2[0]:
    print 'Error'
    print 'Both input files should be the same type to calculate the difference'
    Usage();sys.exit(1)
  
  outName=File1.split('.')[0]+'_plus_'+File2.split('.')[0]+'.h5'

  if k1[0] in ('velocity','temporal_coherence','rmse','mask'):
     dset1 = h5file1[k1[0]].get(k1[0])
     data1=dset1[0:dset1.shape[0],0:dset1.shape[1]]
     dset2 = h5file2[k2[0]].get(k2[0])
     data2=dset2[0:dset2.shape[0],0:dset2.shape[1]]
     
     h5file = h5py.File(outName,'w')
     group=h5file.create_group(k1[0])
     dset = group.create_dataset(k1[0], data=data1+data2, compression='gzip')
     

     for key , value in h5file1[k1[0]].attrs.iteritems():
        group.attrs[key]=value
     h5file.close()
  
  elif 'timeseries' in k1:

     dateList1 = h5file1['timeseries'].keys() 
     dateList2 = h5file2['timeseries'].keys()

     h5timeseries = h5py.File(outName)
     group = h5timeseries.create_group('timeseries')
     for date in dateList1:
        dset1 = h5file1['timeseries'].get(date)
        data1 = dset1[0:dset1.shape[0],0:dset1.shape[1]]
        dset2 = h5file2['timeseries'].get(date)
        data2 = dset2[0:dset2.shape[0],0:dset2.shape[1]]
        dset = group.create_dataset(date, data=data1+data2, compression='gzip')

     for key,value in h5file1['timeseries'].attrs.iteritems():
        group.attrs[key] = value
  
     h5timeseries.close()
  
  elif 'interferograms' in k1:
     ifgramList = h5file1['interferograms'].keys()
     h5igrams = h5py.File(outName)   
     gg = h5igrams.create_group('interferograms')
     for igram in ifgramList:
        dset1=h5file1['interferograms'][igram].get(igram)
        data1 = dset1[0:dset1.shape[0],0:dset1.shape[1]]     
        dset2=h5file2['interferograms'][igram].get(igram)
        data2 = dset2[0:dset2.shape[0],0:dset2.shape[1]]
        group = gg.create_group(igram)
        dset = group.create_dataset(igram, data=data1+data2, compression='gzip')
        for key, value in h5file1['interferograms'][igram].attrs.iteritems():
           group.attrs[key] = value

     try:
       gm = h5igrams.create_group('mask')
       mask = h5file1['mask'].get('mask')
       dset = gm.create_dataset('mask', data=mask, compression='gzip')
     except:
       print 'mask not found' 

     h5igrams.close()

  h5file1.close()
  h5file2.close()

if __name__ == '__main__':

  main(sys.argv[1:])  

