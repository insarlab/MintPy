#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import os
import sys
import time
import datetime

import h5py
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

import pysar._readfile as readfile
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


########################################################################################
def usage():
    print '''usage: save_mat.py  file  [outfile]

This function converts the PySAR hdf5 file formats to the matlab structure and saves to a .mat file.

example:
  save_mat.py  velocity.h5
  save_mat.py  timeseries.h5
    '''
    return


def yyyymmdd2years(date):
    d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])
    yy = np.float(d.year) + np.float(d.month-1)/12 + np.float(d.day-1)/365
    return yy


########################################################################################
def main(argv):
    try:
        File = argv[0]
    except:
        usage();sys.exit(1)

    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    print 'input is '+k+' file: '+File

    try:    matFile = argv[1]
    except: matFile = os.path.splitext(File)[0]+'.mat'
    print 'writing >>> '+matFile

    #####
    h5file = h5py.File(File,'r')
    if k in single_dataset_hdf5_file:
        data = h5file[k].get(k)[:]
    
        V = {}
        V['time_range']=''
        try:
            V['x_first'] = float(atr['X_FIRST'])
            V['y_first'] = float(atr['Y_FIRST'])
            V['x_step'] = float(atr['X_STEP'])
            V['y_step'] = float(atr['Y_STEP'])
            V['x_unit'] = atr['X_UNIT']
            V['y_unit'] = atr['Y_UNIT']
        except:
            V['x_first']=1
            V['y_first']=1
            V['x_step']=1
            V['y_step']=1
            V['x_unit']=''
            V['y_unit']=''
  
        try:  V['wavelength']=float(atr['WAVELENGTH'])
        except:  print 'WAVELENGTH was not found'
        try:  V['sat_height']=float(atr['HEIGHT'])
        except:  print 'HEIGHT was not found'
    
        try:  V['near_range']=float(atr['STARTING_RANGE'])
        except:  print 'STARTING_RANGE was not found'
    
        V['far_range']=''
    
        try:  V['near_LookAng']=float(atr['LOOK_REF1'])
        except:  print 'LOOK_REF1 was not found'
        try:  V['far_LookAng']=float(atr['LOOK_REF2'])
        except:  print 'LOOK_REF2 was not found'
       
        V['earth_radius']=''
        V['Unit']='m/yr'
        V['bperptop']=''
        V['bperpbot']=''
        V['sat']=''
        try:  V['width']=int(atr['WIDTH'])
        except:  print 'WIDTH was not found'
    
        try:  V['file_length']=int(atr['FILE_LENGTH'])
        except:  print 'FILE_LENGTH was not found'
        V['t']=''
        V['date']=''
        V['date_years']=''
        try:     V['sat']   = atr['satellite']
        except:  V['sat']   = ''

        ########################################################
        V['data']=data
        sio.savemat(matFile, {k: V})


    elif 'timeseries' in k:
        epochList = sorted(h5file['timeseries'].keys())
        data_dict={}
        for epoch in epochList:
            print epoch
            d = h5file['timeseries'].get(epoch)
            ts={}
            ts['data'] = d[0:d.shape[0],0:d.shape[1]] 
            try:
                ts['x_first']=float(atr['X_FIRST'])
                ts['y_first']=float(atr['Y_FIRST'])
                ts['x_step']=float(atr['X_STEP'])
                ts['y_step']=float(atr['Y_STEP'])
                ts['x_unit']=atr['X_UNIT']
                ts['y_unit']=atr['Y_UNIT']
            except:
                ts['x_first']=1
                ts['y_first']=1
                ts['x_step']=1
                ts['y_step']=1
                ts['x_unit']=''
                ts['y_unit']=''        
    
            ts['wavelength']=float(atr['WAVELENGTH'])
            ts['sat_height']=float(atr['HEIGHT'])
            ts['near_range']=float(atr['STARTING_RANGE'])
            ts['far_range']=float(atr['STARTING_RANGE1'])
            ts['near_LookAng']=float(atr['LOOK_REF1'])
            ts['far_LookAng']=float(atr['LOOK_REF2'])
            ts['earth_radius']=float(atr['EARTH_RADIUS'])
            ts['Unit']='m'
            ts['bperptop']=float(atr['P_BASELINE_TOP_HDR'])
            ts['bperpbot']=float(atr['P_BASELINE_BOTTOM_HDR'])
            ts['sat']=atr['PLATFORM']
            ts['width']=int(atr['WIDTH'])
            ts['file_length']=int(atr['FILE_LENGTH'])
            ts['t']=np.round((yyyymmdd2years(epoch)-yyyymmdd2years(epochList[0]))*365)
            ts['date']=epoch
            ts['date_years']=yyyymmdd2years(epoch)
              
            data_dict['t'+str(epoch)]=ts  #

        data_dict['Number_of_epochs']=len(epochList)
        data_dict['epoch_dates']=epochList
        sio.savemat(matFile, {k: data_dict})

    h5file.close()
    return


########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

