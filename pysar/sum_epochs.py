#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
    
import sys
import os
import getopt
import time
import datetime
import h5py
import numpy as np

def Usage():
  print '''
  ********************************************
  Calculate the Sum of Time Series Displacement per epoch
      For each epoch, referencing it as the master date,
      get a new time series and calculate the temporal
      average displacement.

  Usage:
      sum_epochs.py timeseries_file [output_file]

      output filename is 'sum_'timeseries_file by default.

  Example:
      sum_epochs.py timeseries.h5
      sum_epochs.py timeseries.h5 ouput.h5

  ********************************************
  '''

def main(argv):
  try: timeSeriesFile=argv[0]
  except: Usage() ; sys.exit(1)

  try:    outname=argv[1]
  except: outname='sum_'+timeSeriesFile

  ########################################################
  print "\n*************** Calculating Sum of Epochs ****************"
  print "Loading time series: " + timeSeriesFile
  h5timeseries=h5py.File(timeSeriesFile) 
  dateList = h5timeseries['timeseries'].keys()
  dateList.sort()

  dateIndex={}
  for ni in range(len(dateList)):
      dateIndex[dateList[ni]]=ni

  dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
  nrows,ncols=np.shape(dset)
  D = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0]*np.shape(dset)[1]),np.float32)

  for date in dateList:
     print date
     dset = h5timeseries['timeseries'].get(date)
     d = dset[0:dset.shape[0],0:dset.shape[1]]
     D[dateIndex[date]][:]=d.flatten(0)

  lt=len(dateList)
  sumD=np.zeros(D.shape)
  print 'calculating epochs sum ...'
  for j in range(lt):
      print j
      sumD[j,:]=np.sum(np.abs(D-D[j,:]),0)/lt

  ##################################################
  #Normalize to 0 and 1 with high atmosphere equal to 0 and no atmosphere equal to 1
  sumD=-1*(sumD-np.max(sumD,0))
  sumD=sumD/np.max(sumD,0)
  ind=np.isnan(sumD)
  sumD[ind]=1

  ##################################################  
  print 'writing to >>> '+outname
  h5sum = h5py.File(outname,'w')
  group = h5sum.create_group('timeseries')
  for date in dateList:
    print date
    dset = group.create_dataset(date, data=np.reshape(sumD[dateIndex[date]][:],[nrows,ncols]), compression='gzip')

  for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group.attrs[key] = value

  print 'Done.'


##################################################  
if __name__ == '__main__':

  main(sys.argv[1:])
