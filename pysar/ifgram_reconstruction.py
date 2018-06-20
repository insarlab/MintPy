#! /usr/bin/env python2
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import getopt
import time
import datetime

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._pysar_utilities as ut


#####################################################################################
def usage():
    print '''usage: ifgram_reconstruction.py  ifgram_file  timeseries_file  [output_name]

Reconstruct interferograms from time-series

arguments:
  ifgram_file     : original interferograms file
  timeseries_file : timeseries file, phase corrected timeseries file prefereably
  output_name     : file name of reconstructed interferograms file
                    default: add prefix 'reconstructed_' to input ifgram_file

example:
  ifgram_reconstruction.py  unwrapIfgram.h5  timeseries_ECMWF_demCor.h5
  ifgram_reconstruction.py  unwrapIfgram.h5  timeseries_ECMWF_demCor.h5  reconstructed_unwrapIfgram.h5
    '''
    return


#####################################################################################
def main(argv):

    ##### Inputs
    try:
        ifgram_file = argv[0]
        timeseries_file = argv[1]
    except:
        usage(); sys.exit(1)
  
    try:    outfile = argv[2]
    except: outfile = 'reconstructed_'+ifgram_file

    atr = readfile.read_attribute(timeseries_file)
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    ##### Read time-series file
    print 'loading timeseries ...'
    h5ts = h5py.File(timeseries_file, 'r')
    date_list = sorted(h5ts['timeseries'].keys())
    date_num = len(date_list)
    timeseries = np.zeros((date_num, length*width))

    print 'number of acquisitions: '+str(date_num)
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        d = h5ts['timeseries'].get(date)[:]
        timeseries[i,:] = d.flatten(0)
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    h5ts.close()
    del d

    range2phase = -4*np.pi/float(atr['WAVELENGTH'])
    timeseries = range2phase*timeseries

    #####  Estimate interferograms from timeseries
    print 'estimating interferograms from timeseries using design matrix from input interferograms'
    A,B = ut.design_matrix(ifgram_file)
    p = -1*np.ones([A.shape[0],1])
    Ap = np.hstack((p,A))
    estData = np.dot(Ap, timeseries)
    del timeseries

    ##### Write interferograms file
    print 'writing >>> '+outfile
    h5 = h5py.File(ifgram_file,'r')
    ifgram_list = sorted(h5['interferograms'].keys())
    ifgram_num = len(ifgram_list)
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    
    h5out = h5py.File(outfile,'w')
    group = h5out.create_group('interferograms')

    print 'number of interferograms: '+str(ifgram_num)
    prog_bar = ptime.progress_bar(maxValue=ifgram_num)
    for i in range(ifgram_num):
        ifgram = ifgram_list[i]
        data = np.reshape(estData[i,:],(length, width))

        gg = group.create_group(ifgram)
        dset = gg.create_dataset(ifgram, data=data, compression='gzip')
        for key, value in h5['interferograms'][ifgram].attrs.iteritems():
            gg.attrs[key] = value
        prog_bar.update(i+1, suffix=date12_list[i])
    prog_bar.close()
    h5.close()
    h5out.close()
    print 'Done.'
    return outfile


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


