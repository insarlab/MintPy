#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#


import sys
import os

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


######################################################################################################
def temporal_coherence(timeseriesFile, ifgramFile):
    '''Calculate temporal coherence based on input timeseries file and interferograms file
    Inputs:
        timeseriesFile - string, path of time series file
        ifgramFile     - string, path of interferograms file
    Output:
        temp_coh - 2D np.array, temporal coherence in float32
    '''
    
    # Basic Info
    atr_ts = readfile.read_attribute(timeseriesFile)
    length = int(atr_ts['FILE_LENGTH'])
    width = int(atr_ts['WIDTH'])
    pixel_num = length * width

    # Read time series data
    h5timeseries = h5py.File(timeseriesFile, 'r')
    date_list = sorted(h5timeseries['timeseries'].keys())
    date_num = len(date_list)

    print "load time series: "+timeseriesFile
    print 'number of acquisitions: '+str(date_num)
    timeseries = np.zeros((date_num, pixel_num), np.float32)
    prog_bar = ptime.progress_bar(maxValue=date_num, prefix='loading: ')
    for i in range(date_num):
        date = date_list[i]
        d = h5timeseries['timeseries'].get(date)[:]
        timeseries[i][:] = d.flatten(0)
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    h5timeseries.close()

    # Convert displacement from meter to radian
    range2phase = -4*np.pi/float(atr_ts['WAVELENGTH'])
    timeseries *= range2phase

    # interferograms data
    print "interferograms file: " + ifgramFile
    atr_ifgram = readfile.read_attribute(ifgramFile)
    h5ifgram   = h5py.File(ifgramFile, 'r')
    ifgram_list = sorted(h5ifgram['interferograms'].keys())
    ifgram_list = ut.check_drop_ifgram(h5ifgram, atr_ifgram, ifgram_list)
    ifgram_num = len(ifgram_list)

    # Design matrix
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    A1,B = ut.design_matrix(ifgramFile, date12_list)
    A0 = -1*np.ones([ifgram_num,1])
    A = np.hstack((A0, A1))

    # Get reference pixel
    try:
        ref_x = int(atr_ts['ref_x'])
        ref_y = int(atr_ts['ref_y'])
        print 'find reference pixel in y/x: [%d, %d]'%(ref_y, ref_x)
    except ValueError:
        print 'No ref_x/y found! Can not calculate temporal coherence without it.'
    
    print 'calculating temporal coherence interferogram by interferogram ...'
    print 'number of interferograms: '+str(ifgram_num)
    temp_coh = np.zeros(pixel_num, dtype=np.float32)+0j
    prog_bar = ptime.progress_bar(maxValue=ifgram_num, prefix='calculating: ')
    for i in range(ifgram_num):
        ifgram = ifgram_list[i]
        # read interferogram
        data = h5ifgram['interferograms'][ifgram].get(ifgram)[:]
        data -= data[ref_y, ref_x]
        data = data.flatten(0)

        # calculate difference between observed and estimated data
        dataEst  = np.dot(A[i,:], timeseries)
        dataDiff = data - dataEst
        temp_coh += np.exp(1j*dataDiff)
        prog_bar.update(i+1, suffix=date12_list[i])
    prog_bar.close()
    del timeseries, data, dataEst, dataDiff
    h5ifgram.close()

    temp_coh = np.array((np.absolute(temp_coh)/ifgram_num).reshape((length,width)), dtype=np.float32)
    return temp_coh


######################################################################################################
USAGE='''usage: temporal_coherence.py [-h] interferograms_file timeseries_file [ output_file ]'''

DESCRIPTION='''Generates temporal coherence map.'''

REFERENCE='''reference:
  Tizzani, P., P. Berardino, F. Casu, P. Euillades, M. Manzo, G. P. Ricciardi, G. Zeni,
  and R. Lanari (2007), Surface deformation of Long Valley Caldera and Mono Basin, 
  California, investigated with the SBAS-InSAR approach, Remote Sens. Environ., 108(3),
  277-289, doi:10.1016/j.rse.2006.11.015.
'''

EXAMPLE='''example:
  temporal_coherence.py  unwrapIfgram.h5  timeseries.h5
  temporal_coherence.py  unwrapIfgram.h5  timeseries.h5  temporalCoherence.h5
'''

def usage():
    print USAGE+'\n\n'+DESCRIPTION+'\n\n'+REFERENCE+'\n'+EXAMPLE
    return


######################################################################################################
def main(argv):
    try:
        ifgramFile     = argv[0]
        timeseriesFile = argv[1]
    except:
        usage(); sys.exit()

    temp_coherence = temporal_coherence(timeseriesFile, ifgramFile)

    try:    tempCohFile = argv[2]
    except: tempCohFile = 'temporalCoherence.h5'
    print 'writing >>> '+tempCohFile
    
    atr = readfile.read_attribute(timeseriesFile)
    atr['FILE_TYPE'] = 'temporal_coherence'
    atr['UNIT'] = '1'
    writefile.write(temp_coherence, atr, tempCohFile)
    print 'Done.'
    return tempCohFile


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


