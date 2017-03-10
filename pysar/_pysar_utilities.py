#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

# timeseries_inversion and Remove_plaane are modified 
# from a software originally written by Scott Baker with 
# the following licence:

###############################################################################
#  Copyright (c) 2011, Scott Baker 
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
############################################################################### 
# Yunjun, Oct 2015: Add glob2radar() and radar2glob() (modified from radar2geo.py
#                       written by Heresh)
# Yunjun, Dec 2015: Use k[0] instead of 'interferograms' in some functions for 
#                       better support of interferograms, coherence and wrapped
# Yunjun, Jan 2016: Add yyyymmdd() and yymmdd()
# Yunjun, Jun 2016: Add print_progress() written by Greenstick from Stack Overflow
#                   Removed remove_plane functions since a better version in _remove_plane
#                   Add inner function ts_inverse() to faster time series inversion
#                   Add P_BASELINE_TIMESERIES attribute to timeseries file.
# Yunjun, Jul 2016: add get_file_list() to support multiple files input
# Yunjun, Aug 2016: add spatial_average()
# Yunjun, Jan 2017: add temporal_average(), nonzero_mask()


import os
import sys
import re
import time
import datetime
import glob
import warnings

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime
import pysar._network as pnet
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file

    


############################################################
def incidence_angle(atr, dimension=2):
    '''Calculate 2D matrix of incidence angle from ROI_PAC attributes, very accurate.
    Input:
        dictionary - ROI_PAC attributes including the following items:
                     STARTING_RANGE
                     RANGE_PIXEL_SIZE
                     EARTH_RADIUS
                     HEIGHT
                     FILE_LENGTH
                     WIDTH
        dimension - int,
                    2 for 2d matrix
                    1 for 1d array
                    0 for one center value
    Output: 2D np.array - incidence angle in degree for each pixel
    '''
    ## Read Attributes
    near_range = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    r  = float(atr['EARTH_RADIUS'])
    H  = float(atr['HEIGHT'])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    
    ## Calculation
    far_range = near_range+dR*width
    incidence_n = (np.pi-np.arccos((r**2+near_range**2-(r+H)**2)/(2*r*near_range)))*180.0/np.pi
    incidence_f = (np.pi-np.arccos((r**2+ far_range**2-(r+H)**2)/(2*r*far_range)))*180.0/np.pi
    incidence_c = (incidence_f+incidence_n)/2
    if dimension == 0:
        return incidence_c
    
    print 'near    incidence angle : '+ str(incidence_n)
    print 'far     incidence angle : '+ str(incidence_f)
    print 'average incidence angle : '+ str(incidence_c)
    
    incidence_x = np.linspace(incidence_n, incidence_f, num=width, endpoint='FALSE')
    if dimension == 1:
        return incidence_x
    else:
        incidence_xy = np.tile(incidence_x,(length,1))
        return incidence_xy


def which(program):
    '''Test if executable exists'''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def get_file_stack(File, maskFile=None):
    '''Get stack file of input File and return the stack 2D matrix
    Input:   File/maskFile - string
    Output:  stack - 2D np.array matrix
    '''
    stack = None
    atr = readfile.read_attribute(File)
    stackFile = os.path.splitext(File)[0]+'_stack.h5'
    
    # Read stack from existed file
    if os.path.isfile(stackFile):
        atrStack = readfile.read_attribute(stackFile)
        if atrStack['WIDTH'] == atr['WIDTH'] and atrStack['FILE_LENGTH'] == atr['FILE_LENGTH']:
            print 'reading stack from existed file: '+stackFile
            stack = readfile.read(stackFile)[0]
    # Calculate stack
    else:
        print 'calculating stack of input file ...'
        stack = stacking(File)
        # Write stack file is input file is multi-dataset (large file size usually)
        if atr['FILE_TYPE'] in multi_group_hdf5_file+multi_dataset_hdf5_file:
            atrStack = atr.copy()
            atrStack['FILE_TYPE'] = 'mask'
            print 'writing stack file >>> '+stackFile
            writefile.write(stack, atrStack, stackFile)

    # set masked out area into NaN
    if maskFile:
        mask = readfile.read(maskFile)[0]
        stack[mask==0] = np.nan
    
    return stack


def nonzero_mask(File, outFile='Mask.h5'):
    '''Generate mask file for non-zero value of input multi-group hdf5 file'''
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    
    mask = np.ones([length, width])
    
    h5 = h5py.File(File,'r')
    igramList = sorted(h5[k].keys())
    for i in range(len(igramList)):
        igram = igramList[i]
        data = h5[k][igram].get(igram)[:]
        
        mask[data==0] = 0
        print_progress(i+1, len(igramList))

    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, atr, outFile)
    
    return outFile


######################################################################################################
def spatial_average(File, mask=None, box=None, saveList=False):
    '''Calculate  Spatial Average.
        Only non-nan pixel is considered.
    Input:
        File : string, path of input file
        mask : 2D np.array, mask file 
        box  : 4-tuple defining the left, upper, right, and lower pixel coordinate
        saveList: bool, save (list of) mean value into text file
    Output:
        meanList : list for float, average value in space for each epoch of input file
    Example:
        meanList = spatial_average('coherence.h5')
        meanList = spatial_average('coherence.h5', mask, saveList=True)
        refList = spatial_average('unwrapIfgram.h5', box=(100,200,101,201))
    '''

    # Baic File Info
    atr  = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    if not box:
        box = (0,0,width,length)
    if not mask is None:
        mask = mask[box[1]:box[3],box[0]:box[2]]

    # Calculate mean coherence list
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum  = len(epochList)

        meanList   = []
        for i in range(epochNum):
            epoch = epochList[i]
            if k in multi_group_hdf5_file:
                dset = h5file[k][epoch].get(epoch)
            elif k in multi_dataset_hdf5_file:
                dset = h5file[k].get(epoch)
            else:  print 'Unrecognized group type: '+k
            
            data = dset[box[1]:box[3],box[0]:box[2]]
            if not mask is None:
                data[mask==0] = np.nan
            ## supress warning 
            ## url - http://stackoverflow.com/questions/29688168/mean-nanmean-and-warning-mean-of-empty-slice
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                meanList.append(np.nanmean(data))
            print_progress(i+1, epochNum, suffix=epoch)
        del data
        h5file.close()
    else:
        data,atr = readfile.read(File, box)
        if not mask is None:
            data[mask==0] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            meanList = [np.nanmean(data)]

    # Write mean coherence list into text file
    if saveList:
        txtFile = os.path.splitext(File)[0]+'_spatialAverage.list'
        print 'write average coherence in space into text file: '+txtFile
        fl = open(txtFile, 'w')
        # 1st column of file
        if k in multi_group_hdf5_file:
            str1_list = pnet.get_date12_list(File)
        elif k in multi_dataset_hdf5_file:
            str1_list = epochList
        else:
            str1_list = [os.path.basename(File)]
        for i in range(len(str1_list)):
            line = str1_list[i]+'    '+str(meanList[i])+'\n'
            fl.write(line)
        fl.close()

    if len(meanList) == 1:
        meanList = meanList[0]
    return meanList


def temporal_average(File, outFile=None):
    '''Calculate temporal average.'''
    # Input File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    h5file = h5py.File(File)
    epochList = sorted(h5file[k].keys())
    epochNum = len(epochList)

    # Calculation
    dMean = np.zeros((length,width))
    for i in range(epochNum):
        epoch = epochList[i]
        if k in multi_group_hdf5_file:
            d = h5file[k][epoch].get(epoch)[:]
        elif k in ['timeseries']:
            d = h5file[k].get(epoch)[:]
        else: print k+' type is not supported currently.'; sys.exit(1)
        dMean += d
        print_progress(i+1, epochNum, suffix=epoch)
    dMean /= float(len(epochList))
    del d
    h5file.close()

    # Output
    if not outFile:
        outFile = os.path.splitext(File)[0]+'_tempAverage.h5'
    print 'writing >>> '+outFile
    h5mean = h5py.File(outFile, 'w')
    group  = h5mean.create_group('mask')
    dset = group.create_dataset(os.path.basename('mask'), data=dMean, compression='gzip')
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5mean.close()

    return outFile


######################################################################################################
def get_file_list(fileList):
    '''Get all existed files matching the input list of file pattern
    Example:
    fileList = get_file_list(['*velocity*.h5','timeseries*.h5'])
    '''
  
    fileListOut = []
    for i in range(len(fileList)):
        file0 = fileList[i]
        fileList0 = glob.glob(file0)
        fileListOut += list(set(fileList0) - set(fileListOut))
    return fileListOut


######################################################################################################
def print_progress(iteration, total, prefix='calculating:', suffix='complete', decimals=1, barLength=50, elapsed_time=None):
    """Print iterations progress - Greenstick from Stack Overflow
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int) 
        barLength   - Optional  : character length of bar (Int) 
        elapsed_time- Optional  : elapsed time in seconds (Int/Float)
    
    Reference: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    if elapsed_time:
        sys.stdout.write('%s [%s] %s%s    %s    %s secs\r' % (prefix, bar, percents, '%', suffix, int(elapsed_time)))
    else:
        sys.stdout.write('%s [%s] %s%s    %s\r' % (prefix, bar, percents, '%', suffix))
    sys.stdout.flush()
    if iteration == total:
        print("\n")

    '''
    Sample Useage:
    for i in range(len(dateList)):
        print_progress(i+1,len(dateList))
    '''
    return


def range_resolution(atr):
    '''Get range resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord'''
    if 'X_FIRST' in atr.keys():
        print 'Input file is in geo coord, no range resolution info.'
        return
    inc_angle = incidence_angle(atr, 0)
    rg_step = float(atr['RANGE_PIXEL_SIZE'])/np.sin(inc_angle/180.0*np.pi)
    return rg_step

def azimuth_resolution(atr):
    '''Get azimuth resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord'''
    if 'X_FIRST' in atr.keys():
        print 'Input file is in geo coord, no azimuth resolution info.'
        return
    Re = float(atr['EARTH_RADIUS'])
    Height = float(atr['HEIGHT'])
    az_step = float(atr['AZIMUTH_PIXEL_SIZE']) *Re/(Re+Height)
    return az_step


#########################################################################
def glob2radar(lat, lon, transFile='geomap*.trans', atr_rdr=dict()):
    '''Convert geo coordinates into radar coordinates.
    Inputs:
        lat/lon    - np.array, float, latitude/longitude
        transFile - string, trans/look up file
        atr_rdr    - dict, attributes of file in radar coord, optional but recommended.
    Output:
        az/rg     - np.array, float, range/azimuth pixel number
        az/rg_res - float, residul/uncertainty of coordinate conversion
    '''

    try:    transFile = glob.glob(transFile)[0]
    except: transFile = None

    ########## Precise conversion using geomap.trans file, if it exists.
    if transFile:
        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3;    # in meter
        print 'reading file: '+transFile
        trans_rg, trans_atr = readfile.read(transFile, (), 'range')
        trans_az = readfile.read(transFile, (), 'azimuth')[0]
        lat_first = float(trans_atr['Y_FIRST'])
        lon_first = float(trans_atr['X_FIRST'])
        lat_center = lat_first + float(trans_atr['Y_STEP'])*float(trans_atr['FILE_LENGTH'])/2
        lat_step_deg = float(trans_atr['Y_STEP'])
        lon_step_deg = float(trans_atr['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0*earth_radius
        lon_step = lon_step_deg*np.pi/180.0*earth_radius*np.sin(lat_center*np.pi/180)
        
        # Get range/azimuth ground resolution/step in meter
        if atr_rdr:
            az_step = azimuth_resolution(atr_rdr)
            rg_step = range_resolution(atr_rdr)
            try:    az0 = int(atr_rdr['subset_y0'])
            except: az0 = 0
            try:    rg0 = int(atr_rdr['subset_x0'])
            except: rg0 = 0
            
            x_factor = np.ceil(abs(lon_step)/rg_step).astype(int)
            y_factor = np.ceil(abs(lat_step)/az_step).astype(int)
        else:
            x_factor = 10
            y_factor = 10
            az0 = 0
            rg0 = 0
        
        width  = int(trans_atr['WIDTH'])
        row = np.rint((lat - lat_first)/lat_step_deg).astype(int)
        col = np.rint((lon - lon_first)/lon_step_deg).astype(int)
        rg = np.rint(trans_rg[row, col]).astype(int) - rg0
        az = np.rint(trans_az[row, col]).astype(int) - az0
        rg_resid = x_factor
        az_resid = y_factor

    ########## Simple conversion using 2D linear transformation, with 4 corners' lalo info
    elif atr_rdr:
        print 'finding approximate radar coordinate with 2D linear transformation estimation.'
        print '    using four corner lat/lon info from '+rdrFile+' file.'
        # This method only works for whole frame/track, thus input file cannot be subsetted before.
        if 'subset_x0' in atr_rdr.keys():
            print 'WARNING: Cannot use subset file as input! No coordinate converted.'
            return None

        LAT_REF1=float(atr_rdr['LAT_REF1'])
        LAT_REF2=float(atr_rdr['LAT_REF2'])
        LAT_REF3=float(atr_rdr['LAT_REF3'])
        LAT_REF4=float(atr_rdr['LAT_REF4'])
        LON_REF1=float(atr_rdr['LON_REF1'])
        LON_REF2=float(atr_rdr['LON_REF2'])
        LON_REF3=float(atr_rdr['LON_REF3'])
        LON_REF4=float(atr_rdr['LON_REF4'])
        W =      float(atr_rdr['WIDTH'])
        L =      float(atr_rdr['FILE_LENGTH'])

        LAT_REF = np.array([LAT_REF1,LAT_REF2,LAT_REF3,LAT_REF4]).reshape(4,1)
        LON_REF = np.array([LON_REF1,LON_REF2,LON_REF3,LON_REF4]).reshape(4,1)
        X = np.array([1,W,1,W]).reshape(4,1)
        Y = np.array([1,1,L,L]).reshape(4,1)
    
        ### estimate 2D tranformation from Lease Square
        A = np.hstack([LAT_REF,LON_REF,np.ones((4,1))])
        B = np.hstack([X,Y])
        affine_par = np.linalg.lstsq(A,B)[0]
        res = B - np.dot(A,affine_par)
        res_mean = np.mean(np.abs(res),0)
        rg_resid = (res_mean[0]+0.5).astype(int)
        az_resid = (res_mean[1]+0.5).astype(int)
        print 'Residul - rg: '+str(rg_resid)+', az: '+str(az_resid)
    
        ### calculate radar coordinate of inputs
        N = len(lat)
        A = np.hstack([lat.reshape(N,1), lon.reshape(N,1), np.ones((N,1))])
        rg = np.rint(np.dot(A, affine_par[:,0])).astype(int)
        az = np.rint(np.dot(A, affine_par[:,1])).astype(int)
  
    return az, rg, az_resid, rg_resid


def radar2glob(az, rg, transFile='geomap*.trans', atr_rdr=dict()):
    '''Convert radar coordinates into geo coordinates
    Inputs:
        rg/az      - np.array, int, range/azimuth pixel number
        transFile - string, trans/look up file
        atr_rdr    - dict, attributes of file in radar coord, optional but recommended.
    Output:
        lon/lat    - np.array, float, longitude/latitude of input point (rg,az)
        latlon_res - float, residul/uncertainty of coordinate conversion
    '''
    try:    transFile = glob.glob(transFile)[0]
    except: transFile = None
    
    lon = np.zeros(rg.shape)
    lat = np.zeros(rg.shape)    

    ##### Use geomap*.trans file for precious (pixel-level) coord conversion
    ## by searching pixels in trans file with value falling buffer lat/lon value
    if transFile:
        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3;    # in meter
        print 'reading file: '+transFile
        trans_rg, trans_atr = readfile.read(transFile, (), 'range')
        trans_az = readfile.read(transFile, (), 'azimuth')[0]
        lat_first = float(trans_atr['Y_FIRST'])
        lon_first = float(trans_atr['X_FIRST'])
        lat_center = lat_first + float(trans_atr['Y_STEP'])*float(trans_atr['FILE_LENGTH'])/2
        lat_step_deg = float(trans_atr['Y_STEP'])
        lon_step_deg = float(trans_atr['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0*earth_radius
        lon_step = lon_step_deg*np.pi/180.0*earth_radius*np.sin(lat_center*np.pi/180)

        # Get range/azimuth ground resolution/step
        if atr_rdr:
            az_step = azimuth_resolution(atr_rdr)
            rg_step = range_resolution(atr_rdr)
        
            x_factor = np.ceil(abs(lon_step)/rg_step)
            y_factor = np.ceil(abs(lat_step)/az_step)
        else:
            x_factor = 10
            y_factor = 10
        
        for i in range(len(rg)):
            mask_rg = np.multiply(trans_rg>=rg[i]-x_factor, trans_rg<=rg[i]+x_factor)
            mask_az = np.multiply(trans_az>=az[i]-y_factor, trans_az<=az[i]+y_factor)
            idx = np.where(np.multiply(mask_rg, mask_az))
            trans_row, trans_col = np.mean(idx,1)
            
            lat[i] = trans_row*lat_step_deg + lat_first
            lon[i] = trans_col*lon_step_deg + lon_first
        
        lat_resid = y_factor*lat_step_deg
        lon_resid = x_factor*lon_step_deg

    ##### Use corner lat/lon for rough (ten-pixels or more level) coord conversion
    ## by estimating a simple 2D linear transformation
    elif atr_rdr:
        # This method only works for whole frame/track, thus input file cannot be subsetted before.
        if 'subset_x0' in atr_rdr.keys():
            print 'WARNING: Cannot use subset file as input! No coordinate converted.'
            return None
        
        LAT_REF1=float(atr_rdr['LAT_REF1'])
        LAT_REF2=float(atr_rdr['LAT_REF2'])
        LAT_REF3=float(atr_rdr['LAT_REF3'])
        LAT_REF4=float(atr_rdr['LAT_REF4'])
        LON_REF1=float(atr_rdr['LON_REF1'])
        LON_REF2=float(atr_rdr['LON_REF2'])
        LON_REF3=float(atr_rdr['LON_REF3'])
        LON_REF4=float(atr_rdr['LON_REF4'])
        W =      float(atr_rdr['WIDTH'])
        L =      float(atr_rdr['FILE_LENGTH'])
        
        LAT_REF = np.array([LAT_REF1,LAT_REF2,LAT_REF3,LAT_REF4]).reshape(4,1)
        LON_REF = np.array([LON_REF1,LON_REF2,LON_REF3,LON_REF4]).reshape(4,1)
        X = np.array([1,W,1,W]).reshape(4,1)
        Y = np.array([1,1,L,L]).reshape(4,1)
        
        ### estimate 2D tranformation from Lease Square
        A = np.hstack([X,Y,np.ones((4,1))])
        B = np.hstack([LAT_REF,LON_REF])
        affine_par = np.linalg.lstsq(A,B)[0]
        res = B - np.dot(A,affine_par)
        res_mean = np.mean(np.abs(res),0)
        lat_resid = res_mean[0]
        lon_resid = res_mean[1]
        print 'Residul - lat: '+str(lat_resid)+', lon: '+str(lon_resid)
        
        ### calculate geo coordinate of inputs
        N = len(rg)
        A = np.hstack([rg.reshape(N,1), az.reshape(N,1), np.ones((N,1))])
        lat = np.dot(A, affine_par[:,0])
        lon = np.dot(A, affine_par[:,1])
        
    else:
        print 'No geomap*.trans or radar coord file found!'
        return None
        
    return lat, lon, lat_resid, lon_resid


#########################################################################
def check_variable_name(path):
    s=path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        p0=os.getenv(s[1:])
        path=path.replace(path.split("/")[0],p0)
    return path


#########################################################################
def hillshade(data,scale):
    '''from scott baker, ptisk library '''
    azdeg  = 315.0
    altdeg = 45.0
    az  = azdeg  * np.pi/180.0
    alt = altdeg * np.pi/180.0
    dx, dy = np.gradient(data/scale)
    slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dx, dy)
    data = np.sin(alt)*np.sin(slope) + np.cos(alt)*np.cos(slope)*np.cos(-az - aspect - 0.5*np.pi)
    return data


#################################################################
def date_list(h5file):
    dateList = []
    tbase = []
    k = h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    ifgramList = sorted(h5file[k[0]].keys())
    for ifgram in  ifgramList:
        dates = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
        dates1= h5file[k[0]][ifgram].attrs['DATE12'].split('-')
        if dates[0][0] == '9':      dates[0] = '19'+dates[0]
        else:                       dates[0] = '20'+dates[0]
        if dates[1][0] == '9':      dates[1] = '19'+dates[1]
        else:                       dates[1] = '20'+dates[1]
        if not dates[0] in dateList: dateList.append(dates[0])
        if not dates[1] in dateList: dateList.append(dates[1])
      
    dateList.sort()
    dateList1=[]
    for ni in range(len(dateList)):
        dateList1.append(dateList[ni][2:])
  
    d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
    for ni in range(len(dateList)):
        d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
        diff = d2-d1
        tbase.append(diff.days)
    dateDict = {}
    for i in range(len(dateList)): dateDict[dateList[i]] = tbase[i]
    return tbase,dateList,dateDict,dateList1


#####################################
def YYYYMMDD2years(d):
    dy = datetime.datetime(*time.strptime(d,"%Y%m%d")[0:5])
    dyy=np.float(dy.year) + np.float(dy.month-1)/12 + np.float(dy.day-1)/365
    return dyy

######################################
def design_matrix(h5file):
    '''Make the design matrix for the inversion.  '''
    tbase,dateList,dateDict,dateList1 = date_list(h5file)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    ifgramList = h5file[k[0]].keys()
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
        date = h5file[k[0]][ifgramList[ni]].attrs['DATE12'].split('-')
        if date[0][0] == '9':      date[0] = '19'+date[0]
        else:                      date[0] = '20'+date[0]
        if date[1][0] == '9':      date[1] = '19'+date[1]
        else:                      date[1] = '20'+date[1]
        ndxt1 = daysList.index(dateDict[date[0]])
        ndxt2 = daysList.index(dateDict[date[1]])
        A[ni,ndxt1] = -1
        A[ni,ndxt2] = 1
        B[ni,ndxt1:ndxt2] = tbase[ndxt1+1:ndxt2+1]-tbase[ndxt1:ndxt2]
        t[ni,:] = [dateDict[date[0]],dateDict[date[1]]]
    A = A[:,1:]
    B = B[:,:-1]
    return A,B

######################################
def timeseries_inversion(igramsFile, timeseriesFile):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 
    
    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms 
      h5timeseries: hdf5 file with the output from the inversion
    '''
    total = time.time()
  
    global B1, dt, numDates
    h5flat = h5py.File(igramsFile,'r')
    A,B = design_matrix(h5flat)
    tbase,dateList,dateDict,dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    B1 = np.linalg.pinv(B)
    B1 = np.array(B1,np.float32)
    numDates = len(dateList)
  
    ##### Basic Info
    ifgramList = h5flat['interferograms'].keys()
    numIfgrams = len(ifgramList)
    atr = readfile.read_attribute(igramsFile)
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    numPixels = length * width
    print 'number of interferograms: '+str(numIfgrams)
    print 'number of pixels        : '+str(numPixels)
    #numPixels_step = int(numPixels/10)

    ##### Inversion Function
    def ts_inverse_point(dataPoint):
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):
            B1tmp = np.dot(B1,np.diag(fin_ndx))
            tmp_rate = np.dot(B1tmp,dataPoint)
            zero = np.array([0.],np.float32)
            defo = np.concatenate((zero,np.cumsum([tmp_rate*dt])))
        else: defo = np.zeros((modelDimension+1,1),np.float32)
        return defo
  
    def ts_inverse(dataLine):
        numPoint = dataLine.shape[1]
        tmp_rate = np.dot(B1,dataLine)
        defo1 = tmp_rate * np.tile(dt.reshape((numDates-1,1)),(1,numPoint))
        defo0 = np.array([0.]*numPoint,np.float32)
        defo  = np.vstack((defo0, np.cumsum(defo1,axis=0)))
        return defo
  
    ##### Read Interferograms
    print 'Reading interferograms ...'
    data = np.zeros((numIfgrams,numPixels),np.float32)
    for j in range(numIfgrams):
        ifgram = ifgramList[j]
        print_progress(j+1, numIfgrams, prefix='loading: ', suffix=ifgram)
        d = h5flat['interferograms'][ifgramList[j]].get(ifgramList[j])[:]
        data[j] = d.flatten(1)
    h5flat.close()

    ##### Inversion
    print 'Inversing time series ...'
    dataPoint = np.zeros((numIfgrams,1),np.float32)
    dataLine  = np.zeros((numIfgrams,width),np.float32)
    modelDimension  = np.shape(B)[1]
    tempDeformation = np.zeros((modelDimension+1,numPixels),np.float32)
    for i in range(length):
        dataLine = data[:,i*width:(i+1)*width]
        defoLine = ts_inverse(dataLine)
        tempDeformation[:,i*width:(i+1)*width] = defoLine
  
        #for j in range(width):
        #    dataPoint = data[:,j]
        #    try: tempDeformation[:,i*length+j] = point_inverse(dataPoint)
        #    except: pass
  
        print_progress(i+1,length,prefix='calculating:')
    del data
  
    ##### Time Series Data Preparation
    print 'converting phase to range'
    timeseries = np.zeros((modelDimension+1,length,width),np.float32)
    phase2range = -1*float(atr['WAVELENGTH'])/(4.*np.pi)
    for ni in range(modelDimension+1):
        timeseries[ni] = tempDeformation[ni].reshape(width,length).T
        timeseries[ni] = timeseries[ni]*phase2range
    del tempDeformation
  
    ##### Output Time Series File
    print 'writing >>> '+timeseriesFile
    print 'number of dates: '+str(numDates)
    h5timeseries = h5py.File(timeseriesFile,'w')
    group = h5timeseries.create_group('timeseries')
    #dateIndex = ptime.date_index(dateList)
    for i in range(numDates):
        date = dateList[i]
        print_progress(i+1, numDates, suffix=date)
        dset = group.create_dataset(date, data=timeseries[i], compression='gzip')

    ## Attributes
    print 'calculating perpendicular baseline timeseries'
    Bperp = Baseline_timeseries(igramsFile)
    Bperp = str(Bperp.tolist()).translate(None,'[],')
    atr['P_BASELINE_TIMESERIES'] = Bperp
    atr['ref_date'] = dateList[0]
    for key,value in atr.iteritems():   group.attrs[key] = value
    h5timeseries.close()
  
    print 'Done.\nTime series inversion took ' + str(time.time()-total) +' secs'
    
###################################################
def timeseries_inversion_FGLS(h5flat,h5timeseries):
    '''Implementation of the SBAS algorithm.
    
    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms 
      h5timeseries: hdf5 file with the output from the inversion
    ##################################################'''
  
    total = time.time()
    A,B = design_matrix(h5flat)
    tbase,dateList,dateDict,dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    B1 = np.linalg.pinv(B)
    B1 = np.array(B1,np.float32)
    ifgramList = h5flat['interferograms'].keys()
    numIfgrams = len(ifgramList)
    #dset = h5flat[ifgramList[0]].get(h5flat[ifgramList[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset=h5flat['interferograms'][ifgramList[0]].get(ifgramList[0])
    data = dset[0:dset.shape[0],0:dset.shape[1]] 
    numPixels = np.shape(data)[0]*np.shape(data)[1]
    print 'Reading in the interferograms'
    #print numIfgrams,numPixels
    print 'number of interferograms: '+str(numIfgrams)
    print 'number of pixels: '+str(numPixels)
    numPixels_step = int(numPixels/10)
  
    data = np.zeros((numIfgrams,numPixels),np.float32)
    for ni in range(numIfgrams):
        dset=h5flat['interferograms'][ifgramList[ni]].get(ifgramList[ni])
        #dset = h5flat[ifgramList[ni]].get(h5flat[ifgramList[ni]].keys()[0])
        d = dset[0:dset.shape[0],0:dset.shape[1]]
        #print np.shape(d)

    del d
    dataPoint = np.zeros((numIfgrams,1),np.float32)
    modelDimension = np.shape(B)[1]
    tempDeformation = np.zeros((modelDimension+1,numPixels),np.float32)
    for ni in range(numPixels):
        dataPoint = data[:,ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):
            B1tmp = np.dot(B1,np.diag(fin_ndx))
            tmpe_ratea = np.dot(B1tmp,dataPoint)
            zero = np.array([0.],np.float32)
            defo = np.concatenate((zero,np.cumsum([tmpe_ratea*dt])))
            tempDeformation[:,ni] = defo
        #if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,numPixels)
        if not np.remainder(ni,numPixels_step):
            print 'Processing point: %8d of %8d, %3d' % (ni,numPixels,(10*ni/numPixels_step))+'%'
    del data
    timeseries = np.zeros((modelDimension+1,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgramList[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(modelDimension+1):
        timeseries[ni] = tempDeformation[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del tempDeformation
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgramList[0]].attrs.iteritems():
        timeseriesDict[key] = value 
  
    dateIndex={}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]]=ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key,value in timeseriesDict.iteritems():
            group.attrs[key] = value
    
    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]], compression='gzip')
    print 'Time series inversion took ' + str(time.time()-total) +' secs'



def timeseries_inversion_L1(h5flat,h5timeseries):
    try:
        from l1 import l1
        from cvxopt import normal,matrix
    except:
        print '-----------------------------------------------------------------------'
        print 'cvxopt should be installed to be able to use the L1 norm minimization.'
        print '-----------------------------------------------------------------------'
        sys.exit(1)
        #modified from sbas.py written by scott baker, 2012 
  
    
    total = time.time()
    A,B = design_matrix(h5flat)
    tbase,dateList,dateDict,dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    BL1 = matrix(B)
    B1 = np.linalg.pinv(B)
    B1 = np.array(B1,np.float32)
    ifgramList = h5flat['interferograms'].keys()
    numIfgrams = len(ifgramList)
    #dset = h5flat[ifgramList[0]].get(h5flat[ifgramList[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset=h5flat['interferograms'][ifgramList[0]].get(ifgramList[0]) 
    data = dset[0:dset.shape[0],0:dset.shape[1]] 
    numPixels = np.shape(data)[0]*np.shape(data)[1]
    print 'Reading in the interferograms'
    print numIfgrams,numPixels
  
    #data = np.zeros((numIfgrams,numPixels),np.float32)
    data = np.zeros((numIfgrams,numPixels))
    for ni in range(numIfgrams):
        dset=h5flat['interferograms'][ifgramList[ni]].get(ifgramList[ni])
        #dset = h5flat[ifgramList[ni]].get(h5flat[ifgramList[ni]].keys()[0])
        d = dset[0:dset.shape[0],0:dset.shape[1]]
        #print np.shape(d)
    
        data[ni] = d.flatten(1)
    del d
    dataPoint = np.zeros((numIfgrams,1),np.float32)
    modelDimension = np.shape(B)[1]
    tempDeformation = np.zeros((modelDimension+1,numPixels),np.float32)
    print data.shape
    DataL1=matrix(data)
    L1ORL2=np.ones((numPixels,1))
    for ni in range(numPixels):
        print ni
        dataPoint = data[:,ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):
          
            B1tmp = np.dot(B1,np.diag(fin_ndx))
            #tmpe_ratea = np.dot(B1tmp,dataPoint)
            try:
                tmpe_ratea=np.array(l1(BL1,DataL1[:,ni]))
                zero = np.array([0.],np.float32)
                defo = np.concatenate((zero,np.cumsum([tmpe_ratea[:,0]*dt])))
            except:
                tmpe_ratea = np.dot(B1tmp,dataPoint)
                L1ORL2[ni]=0      
                zero = np.array([0.],np.float32)
                defo = np.concatenate((zero,np.cumsum([tmpe_ratea*dt])))
    
            tempDeformation[:,ni] = defo
        if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,numPixels)
    del data
    timeseries = np.zeros((modelDimension+1,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgramList[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(modelDimension+1):
        timeseries[ni] = tempDeformation[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del tempDeformation
    L1ORL2=np.reshape(L1ORL2,(np.shape(dset)[1],np.shape(dset)[0])).T
    
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgramList[0]].attrs.iteritems():
            timeseriesDict[key] = value
  
    dateIndex={}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]]=ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key,value in timeseriesDict.iteritems():
            group.attrs[key] = value
  
    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]], compression='gzip')
    print 'Time series inversion took ' + str(time.time()-total) +' secs'
    L1orL2h5=h5py.File('L1orL2.h5','w')
    gr=L1orL2h5.create_group('mask') 
    dset=gr.create_dataset('mask',data=L1ORL2,compression='gzip')
    L1orL2h5.close()

def Baseline_timeseries(igramsFile):
    h5file = h5py.File(igramsFile,'r')
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    igramList = h5file[k[0]].keys()
    Bp_igram=[]
    for igram in igramList:
        Bp_igram.append((float(h5file[k[0]][igram].attrs['P_BASELINE_BOTTOM_HDR'])+\
                         float(h5file[k[0]][igram].attrs['P_BASELINE_TOP_HDR']))/2)
    
    A,B=design_matrix(h5file)
    dateList       = ptime.igram_date_list(igramsFile)
    tbase,dateDict = ptime.date_list2tbase(dateList)
    dt = np.diff(tbase)
  
    Bp_rate=np.dot(np.linalg.pinv(B),Bp_igram)
    zero = np.array([0.],np.float32)
    Bperp = np.concatenate((zero,np.cumsum([Bp_rate*dt])))
    h5file.close()
    
    return Bperp


def dBh_dBv_timeseries(igramsFile):
    h5file = h5py.File(igramsFile)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    igramList = h5file[k[0]].keys()
    dBh_igram=[]
    dBv_igram=[]
    for igram in igramList:
        dBh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_RATE_HDR']))
        dBv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_RATE_HDR']))
    
  
    A,B=design_matrix(h5file)
    tbase,dateList,dateDict,dateList1 = date_list(h5file)
    dt = np.diff(tbase)
  
    Bh_rate=np.dot(np.linalg.pinv(B),dBh_igram)
    zero = np.array([0.],np.float32)
    dBh = np.concatenate((zero,np.cumsum([Bh_rate*dt])))
    
    Bv_rate=np.dot(np.linalg.pinv(B),dBv_igram)
    zero = np.array([0.],np.float32)
    dBv = np.concatenate((zero,np.cumsum([Bv_rate*dt])))
  
    h5file.close()
  
    return dBh,dBv

def Bh_Bv_timeseries(igramsFile):
    h5file = h5py.File(igramsFile)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    igramList = h5file[k[0]].keys()
    Bh_igram=[]
    Bv_igram=[]
    for igram in igramList:
        Bh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_TOP_HDR']))
        Bv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_TOP_HDR']))
  
  
    A,B=design_matrix(h5file)
    tbase,dateList,dateDict,dateList1 = date_list(h5file)
    dt = np.diff(tbase)
  
    Bh_rate=np.dot(np.linalg.pinv(B),Bh_igram)
    zero = np.array([0.],np.float32)
    Bh = np.concatenate((zero,np.cumsum([Bh_rate*dt])))
  
    Bv_rate=np.dot(np.linalg.pinv(B),Bv_igram)
    zero = np.array([0.],np.float32)
    Bv = np.concatenate((zero,np.cumsum([Bv_rate*dt])))
  
    h5file.close()
  
    return Bh,Bv

def stacking(File):
    '''Stack multi-temporal dataset into one
       equivalent to temporal sum
    '''

    ## File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])

    ## Calculation
    stack  = np.zeros([length,width])
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Input File Info
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum  = len(epochList)
        for i in range(epochNum):
            epoch = epochList[i]
            if k == 'timeseries':  data = h5file[k].get(epoch)[:]
            else:                  data = h5file[k][epoch].get(epoch)[:]
            stack += data
            print_progress(i+1,epochNum)
        h5file.close()

    else:
        try: stack, atrStack = readfile.read(File)
        except: print 'Cannot read file: '+File; sys.exit(1)

    return stack


def yymmdd2YYYYMMDD(date):
    if date[0] == '9':      date = '19'+date
    else:                   date = '20'+date
    return date
 

def yyyymmdd(dates):
    datesOut = []
    for date in dates:
        if len(date) == 6:
            if date[0] == '9':  date = '19'+date
            else:               date = '20'+date
        datesOut.append(date)
    return datesOut


def yymmdd(dates):
    datesOut = []
    for date in dates:
        if len(date) == 8:  date = date[2:8]
        datesOut.append(date)
    return datesOut


def make_triangle(dates12,igram1,igram2,igram3):
    dates=[]
    dates.append(igram1.split('-')[0])
    dates.append(igram1.split('-')[1])
    dates.append(igram2.split('-')[1])
    datesyy=[]
    for d in dates:
        datesyy.append(yymmdd2YYYYMMDD(d))
  
    datesyy.sort()
    Igramtriangle=[]
    Igramtriangle.append(datesyy[0][2:]+'-'+datesyy[1][2:])
    Igramtriangle.append(datesyy[0][2:]+'-'+datesyy[2][2:])
    Igramtriangle.append(datesyy[1][2:]+'-'+datesyy[2][2:])
  
    IgramtriangleIndexes=[dates12.index(Igramtriangle[0]),dates12.index(Igramtriangle[1]),dates12.index(Igramtriangle[2])]
    return Igramtriangle,IgramtriangleIndexes

def get_triangles(h5file):
    k=h5file.keys()
    igramList=h5file[k[0]].keys()
   
    dates12=[]
    for igram in igramList:
        dates12.append(h5file[k[0]][igram].attrs['DATE12'])
    Triangles=[]
    Triangles_indexes=[]
    for igram1 in dates12:
        igram1_date1=igram1.split('-')[0]
        igram1_date2=igram1.split('-')[1]
  
        igram2=[]
        igram2_date2=[]
        for d in dates12:
            if igram1_date2==d.split('-')[0]:
                igram2.append(d)
                igram2_date2.append(d.split('-')[1])
    
        igram3=[]
        igram3_date2=[]
        for d in dates12:
            if igram1_date1==d.split('-')[0] and d != igram1:
                igram3.append(d)
                igram3_date2.append(d.split('-')[1])
  
        for date in  igram2_date2:
            if date in igram3_date2:
                Igramtriangle,IgramtriangleIndexes = make_triangle(dates12,igram1,igram2[igram2_date2.index(date)],\
                                                                                  igram3[igram3_date2.index(date)])
                if not Igramtriangle in Triangles:
                    Triangles.append(Igramtriangle)
                    Triangles_indexes.append(IgramtriangleIndexes)
  
    numTriangles = np.shape(Triangles_indexes)[0]
    curls=np.zeros([numTriangles,3],dtype=np.int)
    for i in range(numTriangles):
        curls[i][:]=Triangles_indexes[i]
 
    numIgrams=len(igramList)
    C=np.zeros([numTriangles,numIgrams])
    for ni in range(numTriangles):
        C[ni][curls[ni][0]]=1
        C[ni][curls[ni][1]]=-1
        C[ni][curls[ni][2]]=1
 
    return curls,Triangles,C


def generate_curls(curlfile,h5file,Triangles,curls):
    ifgramList = h5file['interferograms'].keys()
    h5curlfile=h5py.File(curlfile,'w')
    gg = h5curlfile.create_group('interferograms')
    lcurls=np.shape(curls)[0]
    for i in range(lcurls):
        d1=h5file['interferograms'][ifgramList[curls[i,0]]].get(ifgramList[curls[i,0]])
        d2=h5file['interferograms'][ifgramList[curls[i,1]]].get(ifgramList[curls[i,1]])
        d3=h5file['interferograms'][ifgramList[curls[i,2]]].get(ifgramList[curls[i,2]])
        data1=d1[0:d1.shape[0],0:d1.shape[1]]
        data2=d2[0:d2.shape[0],0:d2.shape[1]]
        data3=d3[0:d3.shape[0],0:d3.shape[1]]
 
        print i
        group = gg.create_group(Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2])
        dset = group.create_dataset(Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2],\
                                    data=data1+data3-data2, compression='gzip')
        for key, value in h5file['interferograms'][ifgramList[curls[i,0]]].attrs.iteritems():
            group.attrs[key] = value
 
    h5curlfile.close()


