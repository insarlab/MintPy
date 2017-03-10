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

######################################

def get_data(h5timeseries):

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
  
    lt,numpixels=np.shape(timeseries)

    ######################################################
    tbase=np.zeros((lt,1))
    d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
  
    for ni in range(len(dateList)):
        d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
        diff = d2-d1
        tbase[ni]=diff.days  
  
    return timeseries,tbase,dateList,lt,numpixels,nrows,ncols

def Usage():
    print '''
    Usage:
      
      filter_temporal.py -f timeseries_file -t time_window -o Name of output

      time_window : default is 0.3

    Example:

      filter_temporal.py -f timeseries_LODcor_demCor.h5

      filter_temporal.py -f timeseries_LODcor_demCor.h5 -t 0.5

      filter_temporal.py -f timeseries_LODcor_demCor.h5 -t 0.5 -o smoothed_timseries.h5

    '''

######################################
def main(argv):
    time_win=0.3
    try:  opts, args = getopt.getopt(argv,"h:f:t:o:")
    except getopt.GetoptError:  Usage() ; sys.exit(1)
  
    if opts==[]:    Usage() ; sys.exit(1)
    for opt,arg in opts:
        if opt in ("-h","--help"):  Usage();  sys.exit()
        elif opt == '-f':       file = arg
        elif opt == '-t':       time_win=float(arg)
        elif opt == '-o':       outname = arg      

    ########################################################
    print '-------------------------------' 
    print "Loading the time series: " + file
    h5File = h5py.File(file,'r')
    if 'timeseries' not in h5File.keys():
        print ''' ******************************
          ERROR!
          timeseries does not exist. This function can only be used for the time-sereis files.
          '''
        sys.exit(1)
  
    
    
    timeseries,t_days,dateList,lt,npix,nrows,ncols=get_data(h5File)
    t=t_days/365.
    timeseries_filt=np.zeros((lt,npix))
    print '-------------------------------'
    print 'smoothing the time-series using moving gaussian window '
    print 'time window : ' +str(time_win)
    print ''
    for i in range(lt):
        print str(i+1)+ ' of ' +str(lt)
        dt2=(t[i]-t)**2;
        weight=np.exp(-1*dt2/2/(time_win**2));
        weight=weight/(np.sum(weight));
        Wa=np.tile(weight,(1,npix));
        timeseries_filt[i,:]=np.sum(timeseries*Wa,0);
        timeseries_filt[i,:]=timeseries_filt[i,:]-timeseries_filt[0,:]
    
    #referencing back all epochs to the first epoch
  

    try:     outname
    except:  outname=file.replace('.h5','')+'_Smooth.h5'
    print '-------------------------------'
    print 'writing   >>>> '+outname
  
    h5filt = h5py.File(outname,'w')
    
    group = h5filt.create_group('timeseries')
    i=-1
    for date in dateList: 
        i=i+1   
        print date
        dset = group.create_dataset(date, data=np.reshape(timeseries_filt[i,:],[nrows,ncols]), compression='gzip') 
  
    #  group = h5timeseriesDEMcor.create_group('timeseries')
    for key,value in h5File['timeseries'].attrs.iteritems():
        group.attrs[key] = value
  
    h5File.close()
    h5filt.close()
 
  
if __name__ == '__main__':
    main(sys.argv[1:])

