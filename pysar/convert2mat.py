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
#import matplotlib.cm as cm
#from matplotlib import colors
import getopt
import h5py
import scipy.io as sio
import datetime
import time

def Usage():
    print '''
****************************************************************
****************************************************************
   This function converts the PySAR hdf5 file formats to the 
   matlab structure and saves to a .mat file.
   Current version only  timeseries, velocity, rmse and temporal coherence.

   usage: 
         convert2mat.py file.h5 

   example:
         convert2mat.py velocity.h5
         convert2mat.py timeseries.h5

****************************************************************
****************************************************************
    '''
def yyyymmdd2years(date):
    d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])
    yy = np.float(d.year) + np.float(d.month-1)/12 + np.float(d.day-1)/365
    return yy

def main(argv):
    try:
        File=argv[0]
    except:
        Usage();sys.exit(1)
   
    h5file=h5py.File(File,'r')
    k=h5file.keys()
    matFile=File.split('.')[0]+'.mat'
    print 'writing '+matFile
    if k[0] =='velocity' or k[0] =='rmse' or k[0] =='temporal_coherence':
        dset = h5file[k[0]].get(k[0])
        data = dset[0:dset.shape[0],0:dset.shape[1]]
    
        velocity = {}
        velocity['time_range']=''
        try:
            velocity['x_first']=float(h5file[k[0]].attrs['X_FIRST'])
            velocity['y_first']=float(h5file[k[0]].attrs['Y_FIRST'])
            velocity['x_step']=float(h5file[k[0]].attrs['X_STEP'])
            velocity['y_step']=float(h5file[k[0]].attrs['Y_STEP'])
            velocity['x_unit']=h5file[k[0]].attrs['X_UNIT']
            velocity['y_unit']=h5file[k[0]].attrs['Y_UNIT']
        except:
            velocity['x_first']=1
            velocity['y_first']=1
            velocity['x_step']=1
            velocity['y_step']=1
            velocity['x_unit']=''
            velocity['y_unit']=''
  
        try:  velocity['wavelength']=float(h5file[k[0]].attrs['WAVELENGTH'])
        except:  print 'WAVELENGTH was not found'
        try:  velocity['sat_height']=float(h5file[k[0]].attrs['HEIGHT'])
        except:  print 'HEIGHT was not found'
    
        try:  velocity['near_range']=float(h5file[k[0]].attrs['STARTING_RANGE'])
        except:  print 'STARTING_RANGE was not found'
    
        velocity['far_range']=''
    
        try:  velocity['near_LookAng']=float(h5file[k[0]].attrs['LOOK_REF1'])
        except:  print 'LOOK_REF1 was not found'
        try:  velocity['far_LookAng']=float(h5file[k[0]].attrs['LOOK_REF2'])
        except:  print 'LOOK_REF2 was not found'
       
        velocity['earth_radius']=''
        velocity['Unit']='m/yr'
        velocity['bperptop']=''
        velocity['bperpbot']=''
        velocity['sat']=''
        try:  velocity['width']=int(h5file[k[0]].attrs['WIDTH'])
        except:  print 'WIDTH was not found'
    
        try:  velocity['file_length']=int(h5file[k[0]].attrs['FILE_LENGTH'])
        except:  print 'FILE_LENGTH was not found'
        velocity['t']=''
        velocity['date']=''
        velocity['date_years']=''
        try:     velocity['sat']   = h5file[k[0]].attrs['satellite']
        except:  velocity['sat']   = ''

        ########################################################
        # matFile=File.split('.')[0]+'.mat'
        # print 'writing '+matFile

        if k[0] =='velocity':
            velocity['rate']=data
            sio.savemat(matFile, {'velocity': velocity})
            # sio.savemat('velocity.mat', {'velocity': velocity}) 
    
        elif k[0] =='rmse':
            velocity['data']=data
            sio.savemat(matFile, {'rmse': velocity})
            # sio.savemat('rmse.mat', {'rmse': velocity})
    
        elif k[0] =='temporal_coherence':
            velocity['data']=data
            sio.savemat(matFile, {'temporal_coherence': velocity})
            #  sio.savemat('temporal_coherence.mat', {'temporal_coherence': velocity})       

    elif 'timeseries' in k:
    
        epochList=h5file['timeseries'].keys()
        data_dict={}
        
        for epoch in epochList:
            print epoch
            d = h5file['timeseries'].get(epoch)
            ts={}
            ts['data'] = d[0:d.shape[0],0:d.shape[1]] 
            try:
                ts['x_first']=float(h5file[k[0]].attrs['X_FIRST'])
                ts['y_first']=float(h5file[k[0]].attrs['Y_FIRST'])
                ts['x_step']=float(h5file[k[0]].attrs['X_STEP'])
                ts['y_step']=float(h5file[k[0]].attrs['Y_STEP'])
                ts['x_unit']=h5file[k[0]].attrs['X_UNIT']
                ts['y_unit']=h5file[k[0]].attrs['Y_UNIT']
            except:
                ts['x_first']=1
                ts['y_first']=1
                ts['x_step']=1
                ts['y_step']=1
                ts['x_unit']=''
                ts['y_unit']=''        
    
            ts['wavelength']=float(h5file[k[0]].attrs['WAVELENGTH'])
            ts['sat_height']=float(h5file[k[0]].attrs['HEIGHT'])
            ts['near_range']=float(h5file[k[0]].attrs['STARTING_RANGE'])
            ts['far_range']=float(h5file[k[0]].attrs['STARTING_RANGE1'])
            ts['near_LookAng']=float(h5file[k[0]].attrs['LOOK_REF1'])
            ts['far_LookAng']=float(h5file[k[0]].attrs['LOOK_REF2'])
            ts['earth_radius']=float(h5file[k[0]].attrs['EARTH_RADIUS'])
            ts['Unit']='m'
            ts['bperptop']=float(h5file[k[0]].attrs['P_BASELINE_TOP_HDR'])
            ts['bperpbot']=float(h5file[k[0]].attrs['P_BASELINE_BOTTOM_HDR'])
            ts['sat']=h5file[k[0]].attrs['PLATFORM']
            ts['width']=int(h5file[k[0]].attrs['WIDTH'])
            ts['file_length']=int(h5file[k[0]].attrs['FILE_LENGTH'])
            ts['t']=np.round((yyyymmdd2years(epoch)-yyyymmdd2years(epochList[0]))*365)
            ts['date']=epoch
            ts['date_years']=yyyymmdd2years(epoch)
            
              
            data_dict['t'+str(epoch)]=ts  #
    
        data_dict['Number_of_epochs']=len(epochList)
        data_dict['epoch_dates']=epochList
        sio.savemat(matFile, {'timeseries': data_dict})

    h5file.close()

    
if __name__ == '__main__':
    main(sys.argv[1:])
