#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import getopt
import time
import datetime

from  numpy import pi as pi
from numpy import shape,zeros,ones,hstack,dot,float32,reshape
import h5py

import pysar._pysar_utilities as ut


#####################################################################################
def reconstruct_igrams_from_timeseries(h5timeseries,h5igrams):

  dateList = h5timeseries['timeseries'].keys()

  dateIndex={}
  for ni in range(len(dateList)):
      dateIndex[dateList[ni]]=ni

  dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
  nrows,ncols=shape(dset)
  timeseries = zeros((len(h5timeseries['timeseries'].keys()),shape(dset)[0]*shape(dset)[1]),float32)
  for date in dateList:
      dset = h5timeseries['timeseries'].get(date)
      d = dset[0:dset.shape[0],0:dset.shape[1]]
      timeseries[dateIndex[date]][:]=d.flatten(0)
  del d

  lt,numpixels=shape(timeseries)
  A,B = ut.design_matrix(h5igrams)

  lam = float(h5timeseries['timeseries'].attrs['WAVELENGTH'])
  range2phase=-4*pi/lam
  timeseries=range2phase*timeseries

  p=-1*ones([A.shape[0],1])
  Ap=hstack((p,A))
  estData=dot(Ap,timeseries)

  return estData,nrows,ncols

#####################################################################################
def Usage():
  print '''
  ****************************************************************************************
  Reconstructs the interferograms from the time-series epochs

  usage:
       reconstruct_igrams.py simulatedIgrams.h5 timeseries.h5  [output_name]

       output reconstructed_simulatedIgrams.h5 by default.

  Example:
     
       reconstruct_igrams.py Seeded_LoadedData.h5 timeseries_ECMWF_demCor.h5
       reconstruct_igrams.py Seeded_LoadedData.h5 timeseries_ECMWF_demCor.h5  reconstructedIgrams.h5

  ****************************************************************************************
  '''


#####################################################################################
def main(argv):

  ##### Inputs
  try:
     igramFile = argv[0]
     tsFile    = argv[1]
  except:
     Usage() ; sys.exit(1)

  try:
     h5igrams     = h5py.File(igramFile,'r')
     h5timeseries = h5py.File(tsFile,'r')
  except:
     Usage() ; sys.exit(1)

  try:    outName = argv[2]
  except: outName = 'reconstructed_'+igramFile

  #####
  estData,nrows,ncols=reconstruct_igrams_from_timeseries(h5timeseries,h5igrams)
  
  h5estIgram=h5py.File(outName,'w')
  igramList=h5igrams['interferograms'].keys()
  gg = h5estIgram.create_group('interferograms')
  for i in range(len(igramList)):
      print igramList[i]
      data=reshape(estData[i,:],(nrows,ncols))
      group = gg.create_group(igramList[i])
      dset = group.create_dataset(igramList[i], data=data, compression='gzip')
      for key, value in h5igrams['interferograms'][igramList[i]].attrs.iteritems():
          group.attrs[key] = value     
  
  print '*****************************'    
  print 'writing to ' + outName
  print '*****************************'

  h5estIgram.close()   


#####################################################################################

if __name__ == '__main__':
  main(sys.argv[1:])


