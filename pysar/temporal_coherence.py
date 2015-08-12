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
import matplotlib.pyplot as plt


######################################
######################################
def date_list(h5file): 
  dateList = []
  tbase = []
  for ifgram in h5file['interferograms'].keys():
    dates = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
    if dates[0][0] == '9':
      dates[0] = '19'+dates[0]
    else:
      dates[0] = '20'+dates[0]
    if dates[1][0] == '9':
      dates[1] = '19'+dates[1]
    else:
      dates[1] = '20'+dates[1]
    if not dates[0] in dateList: dateList.append(dates[0])
    if not dates[1] in dateList: dateList.append(dates[1])
  dateList.sort()
  d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
  for ni in range(len(dateList)):
    d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
    diff = d2-d1
    tbase.append(diff.days)
  dateDict = {}
  for i in range(len(dateList)): dateDict[dateList[i]] = tbase[i]
  return tbase,dateList,dateDict

def design_matrix(h5file):
  '''Make the design matrix for the inversion.  '''
  tbase,dateList,dateDict = date_list(h5file)
  ifgramList = h5file['interferograms'].keys()
  numDates = len(dateDict)
  numIfgrams = len(ifgramList)
  A = np.zeros((numIfgrams,numDates))
  B = np.zeros(np.shape(A))
  daysList = []
  for day in tbase:
    daysList.append(day)
  tbase = np.array(tbase)
  t = np.zeros((numIfgrams,2))
  for ni in range(numIfgrams):
    date = h5file['interferograms'][ifgramList[ni]].attrs['DATE12'].split('-')
    if date[0][0] == '9':
      date[0] = '19'+date[0]
    else:
      date[0] = '20'+date[0]
    if date[1][0] == '9':
      date[1] = '19'+date[1]
    else:
      date[1] = '20'+date[1]
    ndxt1 = daysList.index(dateDict[date[0]])
    ndxt2 = daysList.index(dateDict[date[1]])
    A[ni,ndxt1] = -1
    A[ni,ndxt2] = 1
    B[ni,ndxt1:ndxt2] = tbase[ndxt1+1:ndxt2+1]-tbase[ndxt1:ndxt2]
    t[ni,:] = [dateDict[date[0]],dateDict[date[1]]]
  A = A[:,1:]
  B = B[:,:-1]
  return A,B


def Usage():
  
  print '''
  *************************************************************************
  *************************************************************************
  Generates a parameter called temporal coherence for every pixel.

    Usage:

        temporal_coherence.py inteferogramsfile timeseriesfile'

    Example:

        temporal_coherence.py Seeded_Loaded_interferograms.h5 timeseries.h5

  **************************************************************************
  **************************************************************************
'''
######################################
def main(argv):
  try:
    igramsFile = argv[0]
    timeSeriesFile=argv[1]
  except:
    Usage() ; sys.exit(1)

########################################################
  print "Loading time series: " + timeSeriesFile
  h5timeseries = h5py.File(timeSeriesFile)
  dateList = h5timeseries['timeseries'].keys()
  
  dateIndex={}
  for ni in range(len(dateList)):
    dateIndex[dateList[ni]]=ni 

  dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
  nrows,ncols=np.shape(dset)
  timeseries = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0]*np.shape(dset)[1]),np.float32)
  for date in dateList:
    dset = h5timeseries['timeseries'].get(date)
    d = dset[0:dset.shape[0],0:dset.shape[1]]
    timeseries[dateIndex[date]][:]=d.flatten(0)
  del d

#  h5timeseries.close()

  lt,numpixels=np.shape(timeseries)

######################################################
  print "Loading interferograms: " + igramsFile
  h5igrams = h5py.File(igramsFile)  
  ifgramList = h5igrams['interferograms'].keys()
  numIfgrams = len(ifgramList)

  data = np.zeros((numIfgrams,numpixels),np.float32)
  for ni in range(numIfgrams):
#    dset = h5igrams[ifgramList[ni]].get(h5igrams[ifgramList[ni]].keys()[0])
    dset=h5igrams['interferograms'][ifgramList[ni]].get(ifgramList[ni])
    d = dset[0:dset.shape[0],0:dset.shape[1]]
    data[ni] = d.flatten(0)
  del d
  
  A,B = design_matrix(h5igrams)
  h5igrams.close() 

###################################################
 # lam=0.056235646793737
  lam = float(h5timeseries['timeseries'].attrs['WAVELENGTH'])
  range2phase=-4*np.pi/lam
  timeseries=range2phase*timeseries

  p=-1*np.ones([A.shape[0],1])
  Ap=np.hstack((p,A))
  estData=np.dot(Ap,timeseries)
  dataDiff=data-estData
  
  #print list(data[:,0])
  #print '**********************'
  #print '**********************'
  #print list(estData[:,0])
  #print '**********************'
 # qq1=np.absolute(np.sum(np.exp(1j*(dataDiff)/data),0))/numIfgrams
 # print qq1
 # print '**********************'
  qq=np.absolute(np.sum(np.exp(1j*dataDiff),0))/numIfgrams
 # print qq
 # print '***********************'
#  C1=np.zeros([2,numIfgrams])
 # print np.shape(data)
 # print np.shape(estData)
 # C1[0][:]=data[:,0]
 # C1[1][:]=estData[:,0]
 # print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 # print ''
 # print 'Correlation of the velocity with the DEM:  '+ str(np.corrcoef(C1)[0][1])
 # print np.corrcoef(C1)

 # import matplotlib.pyplot as ply
 # plt.plot(data,'-^',ms=5, alpha=0.7, mfc='green')
 # plt.plot(estData,'--s',ms=5, alpha=0.7, mfc='red')
 # plt.show() 
 
  Temp_Coh=np.reshape(qq,[nrows,ncols])

  h5file = 'temporal_coherence.h5'
  h5TempCoh = h5py.File(h5file,'w')
  group=h5TempCoh.create_group('temporal_coherence')
  dset = group.create_dataset(os.path.basename('temporal_coherence'), data=Temp_Coh, compression='gzip')
  for key , value in h5timeseries['timeseries'].attrs.iteritems():
     group.attrs[key]=value 

#  group.attrs['date1'] = datevector[0]
#  group.attrs['date2'] = datevector[40]
  h5TempCoh.close()
  h5timeseries.close()
#  plt.imshow(Temp_Coh)
#  plt.colorbar()
#  plt.show() 

if __name__ == '__main__':
#  sys.argv.append('--clim=0.05')  
#  sys.argv.append('-i0,-1') 
#  sys.argv.append('-m0.3') 
#  sys.argv.append('-f/Users/sbaker/Desktop/data/TSSAR/timeseries_HawaiiT322EnvA1.h5')
#  sys.argv.append('-d/Users/sbaker/Desktop/data/TSSAR/hawaii_filled_0.00056.dem')

  main(sys.argv[1:])

  
