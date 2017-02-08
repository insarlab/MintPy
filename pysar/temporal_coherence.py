#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
#


import sys
import os
import getopt
import time
import datetime

import h5py
import numpy as np
import matplotlib.pyplot as plt

import pysar._readfile as readfile
import pysar._pysar_utilities as ut


######################################################################################################
def date_list(h5file): 
    dateList = []
    tbase    = []
    for ifgram in h5file['interferograms'].keys():
        dates = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
        if dates[0][0] == '9':  dates[0] = '19'+dates[0]
        else:                   dates[0] = '20'+dates[0]
        if dates[1][0] == '9':  dates[1] = '19'+dates[1]
        else:                   dates[1] = '20'+dates[1]
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


######################################################################################################
def design_matrix(h5file):
    '''Make the design matrix for the inversion.  '''
    tbase,dateList,dateDict = date_list(h5file)
    ifgramList = h5file['interferograms'].keys()
    numDates   = len(dateDict)
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
        if date[0][0] == '9':  date[0] = '19'+date[0]
        else:                  date[0] = '20'+date[0]
        if date[1][0] == '9':  date[1] = '19'+date[1]
        else:                  date[1] = '20'+date[1]
        ndxt1 = daysList.index(dateDict[date[0]])
        ndxt2 = daysList.index(dateDict[date[1]])
        A[ni,ndxt1] = -1
        A[ni,ndxt2] = 1
        B[ni,ndxt1:ndxt2] = tbase[ndxt1+1:ndxt2+1]-tbase[ndxt1:ndxt2]
        t[ni,:] = [dateDict[date[0]],dateDict[date[1]]]
    A = A[:,1:]
    B = B[:,:-1]
    return A,B


######################################################################################################
def usage():
    print '''
***************************************************************************************
  Generates a parameter called temporal coherence for every pixel.

  Usage:
      temporal_coherence.py inteferograms_file timeseries_file [output_name]

  Example:
      temporal_coherence.py Seeded_unwrapIfgram.h5 timeseries.h5
      temporal_coherence.py Seeded_unwrapIfgram.h5 timeseries.h5 temporal_coherence.h5

  Reference:
      Tizzani, P., P. Berardino, F. Casu, P. Euillades, M. Manzo, G. P. Ricciardi, G. Zeni,
      and R. Lanari (2007), Surface deformation of Long Valley Caldera and Mono Basin, 
      California, investigated with the SBAS-InSAR approach, Remote Sens. Environ., 108(3),
      277-289, doi:10.1016/j.rse.2006.11.015.

      Gourmelen, N., F. Amelung, and R. Lanari (2010), Interferometric synthetic aperture
      radar-GPS integration: Interseismic strain accumulation across the Hunter Mountain 
      fault in the eastern California shear zone, J. Geophys. Res., 115, B09408, 
      doi:10.1029/2009JB007064.

***************************************************************************************
    '''


######################################
def main(argv):
    try:
        igramsFile     = argv[0]
        timeSeriesFile = argv[1]
    except:
        usage() ; sys.exit(1)

    try:    tempCohFile = argv[2]
    except: tempCohFile = 'temporal_coherence.h5'

    ########################################################
    print '\n********** Temporal Coherence ****************'
    print "load time series: "+timeSeriesFile
    atr_ts = readfile.read_attribute(timeSeriesFile)
    h5timeseries = h5py.File(timeSeriesFile)
    dateList = h5timeseries['timeseries'].keys()
    numDates = len(dateList)

    print 'number of epoch: '+str(numDates)
    dateIndex={}
    for ni in range(numDates):
        dateIndex[dateList[ni]]=ni 

    dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
    nrows,ncols=np.shape(dset)
    timeseries = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0]*np.shape(dset)[1]),np.float32)

    for i in range(numDates):
        date = dateList[i]
        dset = h5timeseries['timeseries'].get(date)
        d = dset[0:dset.shape[0],0:dset.shape[1]]
        timeseries[dateIndex[date]][:]=d.flatten(0)
        ut.printProgress(i+1,numDates,'loading:',date)
    del d
    h5timeseries.close()

    lt,numpixels=np.shape(timeseries)
    range2phase = -4*np.pi/float(atr_ts['WAVELENGTH'])
    timeseries = range2phase*timeseries

    ######################################################
    print "load interferograms: " + igramsFile
    h5igrams   = h5py.File(igramsFile)  
    ifgramList = h5igrams['interferograms'].keys()
    numIfgrams = len(ifgramList)
    print 'number of epochs: '+str(numIfgrams)
    A,B = design_matrix(h5igrams)
    p   = -1*np.ones([A.shape[0],1])
    Ap  = np.hstack((p,A))

    print 'calculating temporal coherence ...'
    #data = np.zeros((numIfgrams,numpixels),np.float32)
    qq = np.zeros(numpixels)+0j
    for ni in range(numIfgrams):
        ## read interferogram
        igram = ifgramList[ni]
        data = h5igrams['interferograms'][igram].get(igram)[:].flatten(0)

        ## calculate difference between observed and estimated data
        ## interferogram by interferogram, less memory, Yunjun - 2016.06.10
        dataEst  = np.dot(Ap[ni,:], timeseries)
        dataDiff = data - dataEst
        qq += np.exp(1j*dataDiff)

        ## progress bar
        ut.printProgress(ni+1,numIfgrams,'calculating:',igram)
    del timeseries, data, dataEst, dataDiff
    h5igrams.close()

    #qq=np.absolute(np.sum(np.exp(1j*dataDiff),0))/numIfgrams
    qq = np.absolute(qq)/numIfgrams
    Temp_Coh=np.reshape(qq,[nrows,ncols])

    ##### write temporal coherence file ####################
    print 'writing >>> '+tempCohFile
    h5TempCoh = h5py.File(tempCohFile,'w')
    group=h5TempCoh.create_group('temporal_coherence')
    dset = group.create_dataset(os.path.basename('temporal_coherence'), data=Temp_Coh, compression='gzip')
    for key , value in atr_ts.iteritems():
        group.attrs[key]=value
    group.attrs['UNIT'] = '1'
    h5TempCoh.close()


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


