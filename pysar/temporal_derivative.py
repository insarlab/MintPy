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
from numpy import sum,remainder,zeros,dot,reshape, float32, array, hstack, vstack, linalg, eye, ones
from scipy.stats import nanstd, nanmean
######################################
######################################
def Usage():
  print '''
    ***********************************************
    Usage:
    
      temporal_derivative.py  timeSeriesFile 

    Example:

      temporal_derivative.py  timeseries.h5 

   ***********************************************   
'''

######################################
def main(argv):
  try:
    timeSeriesFile = argv[0]
   # igramsFile = argv[1]
  except:
    Usage() ; sys.exit(1)


  ########################################################
  print "Loading time series: " + timeSeriesFile
  h5timeseries = h5py.File(timeSeriesFile)
  dateList = h5timeseries['timeseries'].keys()

  tbase=[]
  d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
  for ni in range(len(dateList)):
    print dateList[ni]
    d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
    diff = d2-d1
    tbase.append(diff.days)


 # dateList.remove(h5timeseries['timeseries'].attrs['ref_date'])
 # tbase.remove(0)
  tbase=array(tbase,float32)
   
  nrows=float(h5timeseries['timeseries'].attrs['FILE_LENGTH'])
  ncols=float(h5timeseries['timeseries'].attrs['WIDTH'])
  lt=len(dateList)
  npix=nrows*ncols
 # maskFile=argv[1]
 # h5mask=h5py.File(maskFile,'r')
 # maskSet=h5mask['mask'].get('mask')
 # Mask=maskSet[0:maskSet.shape[0],0:maskSet.shape[1]]
 # ndx = Mask !=0
 # npixels=sum(Mask)
 # n_all_pixels=Mask.size
 # nrows,ncols=Mask.shape 

  timeseries = zeros((lt,npix),float32)
  for i in range(lt):
    date=dateList[i]
    dset= h5timeseries['timeseries'].get(date)
    d=dset[0:dset.shape[0],0:dset.shape[1]]
    timeseries[i][:]=d.flatten(0)

  print "++++++++++++++++++++++++++++++++++++++"
  print "temporal fisrt derivative"
  timeseries_1st = zeros((lt-1,npix),float32)
  for i in range(lt-1):
    print i
    timeseries_1st[i][:]=timeseries[i+1][:]-timeseries[i][:]

  print "++++++++++++++++++++++++++++++++++++++"
  print "temporal second derivative"
  timeseries_2nd = zeros((lt-2,npix),float32)
  for i in range(lt-2):
    print i
    timeseries_2nd[i][:]=timeseries_1st[i+1][:]-timeseries_1st[i][:]
  print "++++++++++++++++++++++++++++++++++++++"
  print 'writing first_derivative.h5'
  h51=h5py.File('first_derivative.h5','w')
  h52=h5py.File('second_derivative.h5','w')
  group1=h51.create_group('timeseries')
  group2=h52.create_group('timeseries')

  for i in range(lt-1):
      date=dateList[i+1]
      print date
      dset = group1.create_dataset(date, data=reshape(timeseries_1st[i][:],[nrows,ncols]), compression='gzip')
  for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group1.attrs[key] = value
  print 'writing second_derivative.h5'
  for i in range(lt-2):
      date=dateList[i+2]
      print date
      dset = group2.create_dataset(date, data=reshape(timeseries_2nd[i][:],[nrows,ncols]), compression='gzip')
  for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group2.attrs[key] = value

  print "++++++++++++++++++++++++++++++++++++++"
if __name__ == '__main__':
  main(sys.argv[1:])
