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
# Yunjun, Jun 2016: Removed remove_plane functions since a better version in _remove_plane
#                   Add inner function ts_inverse() to faster time series inversion
#                   Add P_BASELINE_TIMESERIES attribute to timeseries file.
# Yunjun, Jul 2016: add get_file_list() to support multiple files input
# Yunjun, Aug 2016: add spatial_average()
# Yunjun, Jan 2017: add temporal_average(), nonzero_mask()


import os
import sys
import time
import datetime
import glob
import warnings

import h5py
import numpy as np
import multiprocessing

import pysar
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime
import pysar._network as pnet
import pysar._remove_surface as rm
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file



def circle_index(atr,circle_par):
    ## Return Index of Elements within a Circle
    ## Inputs:
    ##     atr        : (dictionary) attributes containging width, length info
    ##     circle_par : (string) center_x,center_y,radius

    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    cir_par = circle_par.split(',')
    #import pdb; pdb.set_trace()
    try:
        c_y    = int(cir_par[0])
        c_x    = int(cir_par[1])
        radius = int(float(cir_par[2]))
    except:
        try:
            c_lat  = float(cir_par[0])
            c_lon  = float(cir_par[1])
            radius = int(float(cir_par[2]))
            c_y = subset.coord_geo2radar(c_lat,atr,'lat')
            c_x = subset.coord_geo2radar(c_lon,atr,'lon')
        except:
            print '\nERROR: Unrecognized circle index format: '+circle_par
            print 'Supported format:'
            print '--circle 200,300,20            for radar coord input'
            print '--circle 31.0214,130.5699,20   for geo   coord input\n'
            return 0

    y,x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx


def update_template_file(template_file, template_dict):
    '''Update option value in template_file with value from input template_dict'''
    tmp_file = template_file+'.tmp'
    f_orig = open(template_file, 'r')
    f_tmp = open(tmp_file, 'w')
    for line in f_orig:
        line = line.strip()
        c = [i.strip() for i in line.split('=', 1)]
        if not line.startswith('%') and not line.startswith('#') and len(c) > 1:
            key = c[0]
            value = str.replace(c[1],'\n','').split("#")[0].strip()
            if key in template_dict.keys() and template_dict[key] != value:
                line = line.replace(value, template_dict[key], 1)
                print '    '+key+': '+value+' --> '+template_dict[key]
        f_tmp.write(line+'\n')
    f_orig.close()
    f_tmp.close()

    # Overwrite exsting original template file
    mvCmd = 'mv '+tmp_file+' '+template_file
    os.system(mvCmd)

    return template_file


def get_residual_std(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    '''Calculate deramped standard deviation in space for each epoch of input timeseries file.
    Inputs:
        timeseries_resid_file - string, timeseries HDF5 file, e.g. timeseries_ECMWF_demErrInvResid.h5
        mask_file - string, mask file, e.g. maskTempCoh.h5
        ramp_type - string, ramp type, e.g. plane, quadratic, no for do not remove ramp
    outputs:
        std_list  - list of float, standard deviation of deramped input timeseries file
        date_list - list of string in YYYYMMDD format, corresponding dates
    Example:
        import pysar._pysar_utilities as ut
        std_list, date_list = ut.get_residual_std('timeseries_ECMWF_demErrInvResid.h5', 'maskTempCoh.h5')
    '''
    # Intermediate files name
    if ramp_type == 'no':
        print 'No ramp removal'
        deramp_file = timeseries_resid_file
    else:
        deramp_file = os.path.splitext(timeseries_resid_file)[0]+'_'+ramp_type+'.h5'
    std_file = os.path.splitext(deramp_file)[0]+'_std.txt'

    # Get residual std text file
    if update_file(std_file, [deramp_file,mask_file], check_readable=False):
        if update_file(deramp_file, timeseries_resid_file):
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print 'removing a '+ramp_type+' ramp from file: '+timeseries_resid_file
                deramp_file = rm.remove_surface(timeseries_resid_file, ramp_type, mask_file, deramp_file)
        print 'Calculating residual standard deviation for each epoch from file: '+deramp_file
        std_file = timeseries_std(deramp_file, mask_file, std_file)

    # Read residual std text file
    print 'read timeseries RSD from file: '+std_file
    std_fileContent = np.loadtxt(std_file, dtype=str)
    std_list = std_fileContent[:,1].astype(np.float).tolist()
    date_list = list(std_fileContent[:,0]) 
    
    return std_list, date_list


def timeseries_std(inFile, maskFile='maskTempCoh.h5', outFile=None):
    '''Calculate the standard deviation for each epoch of input timeseries file
    and output result to a text file.
    '''
    try:
        mask = readfile.read(maskFile)[0]
        print 'read mask from file: '+maskFile
    except:
        maskFile = None
        print 'no mask input, use all pixels'

    if not outFile:
        outFile = os.path.splitext(inFile)[0]+'_std.txt'

    atr = readfile.read_attribute(inFile)
    k = atr['FILE_TYPE']
    if not k in ['timeseries']:
        raise Exception('Only timeseries file is supported, input file is: '+k)

    h5 = h5py.File(inFile, 'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)

    f = open(outFile, 'w')
    f.write('# Residual Standard Deviation in space for each epoch of timeseries\n')
    f.write('# Timeseries file: '+inFile+'\n')
    f.write('# Mask file: '+maskFile+'\n')
    f.write('# Date      STD(m)\n')
    for i in range(date_num):
        date = date_list[i]
        data = h5[k].get(date)[:]
        if maskFile:
            data[mask==0] = np.nan
        std = np.nanstd(data)
        msg = '%s    %.4f' % (date, std)
        f.write(msg+'\n')
        print msg
    h5.close()
    f.close()
    print 'write to '+outFile

    return outFile


def get_residual_rms(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    '''Calculate deramped Root Mean Square in space for each epoch of input timeseries file.
    Inputs:
        timeseries_resid_file - string, timeseries HDF5 file, e.g. timeseries_ECMWF_demErrInvResid.h5
        mask_file - string, mask file, e.g. maskTempCoh.h5
        ramp_type - string, ramp type, e.g. plane, quadratic, no for do not remove ramp
    outputs:
        rms_list  - list of float, Root Mean Square of deramped input timeseries file
        date_list - list of string in YYYYMMDD format, corresponding dates
    Example:
        import pysar._pysar_utilities as ut
        rms_list, date_list = ut.get_residual_rms('timeseries_ECMWF_demErrInvResid.h5', 'maskTempCoh.h5')
    '''
    # Intermediate files name
    if ramp_type == 'no':
        print 'No ramp removal'
        deramp_file = timeseries_resid_file
    else:
        deramp_file = os.path.splitext(timeseries_resid_file)[0]+'_'+ramp_type+'.h5'
    rms_file = os.path.splitext(deramp_file)[0]+'_rms.txt'

    # Get residual RMS text file
    if update_file(rms_file, [deramp_file,mask_file], check_readable=False):
        if update_file(deramp_file, timeseries_resid_file):
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print 'removing a '+ramp_type+' ramp from file: '+timeseries_resid_file
                deramp_file = rm.remove_surface(timeseries_resid_file, ramp_type, mask_file, deramp_file)
        print 'Calculating residual RMS for each epoch from file: '+deramp_file
        rms_file = timeseries_rms(deramp_file, mask_file, rms_file)

    # Read residual RMS text file
    print 'read timeseries RSD from file: '+rms_file
    rms_fileContent = np.loadtxt(rms_file, dtype=str)
    rms_list = rms_fileContent[:,1].astype(np.float).tolist()
    date_list = list(rms_fileContent[:,0]) 
    
    return rms_list, date_list


def timeseries_rms(inFile, maskFile='maskTempCoh.h5', outFile=None, dimension=2):
    '''Calculate the Root Mean Square for each epoch of input timeseries file
    and output result to a text file.
    '''
    try:
        mask = readfile.read(maskFile)[0]
        print 'read mask from file: '+maskFile
    except:
        maskFile = None
        print 'no mask input, use all pixels'

    if not outFile:
        outFile = os.path.splitext(inFile)[0]+'_rms.txt'

    atr = readfile.read_attribute(inFile)
    k = atr['FILE_TYPE']
    if not k in ['timeseries']:
        raise Exception('Only timeseries file is supported, input file is: '+k)

    h5 = h5py.File(inFile, 'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)

    f = open(outFile, 'w')
    f.write('# Root Mean Square in space for each epoch of timeseries\n')
    f.write('# Timeseries file: '+inFile+'\n')
    f.write('# Mask file: '+maskFile+'\n')
    if dimension == 2:
        f.write('# Date      RMS(m)\n')
        for i in range(date_num):
            date = date_list[i]
            data = h5[k].get(date)[:]
            if maskFile:
                data[mask==0] = np.nan
            rms = np.sqrt(np.nanmean(np.square(data)))
            msg = '%s    %.4f' % (date, rms)
            f.write(msg+'\n')
            print msg
        h5.close()
        f.close()
        print 'write to '+outFile
        return outFile

    elif dimension == 3:
        length = int(atr['FILE_LENGTH'])
        width = int(atr['WIDTH'])
        ts_data = np.zeros((date_num, length*width))
        for i in range(date_num):
            data = h5[k].get(date_list[i])[:]
            if maskFile:
                data[mask==0] = np.nan
            ts_data[i,:] = data.flatten()
        rms = np.sqrt(np.nanmean(np.square(ts_data)))
        return rms


def timeseries_coherence(inFile, maskFile='maskTempCoh.h5', outFile=None):
    '''Calculate spatial average coherence for each epoch of input time series file
    Inputs:
        inFile   - string, timeseries HDF5 file
        maskFile - string, mask file 
        outFile  - string, output text file 
    Example:
        txtFile = timeseries_coherence('timeseries_ECMWF_demErrInvResid_quadratic.h5')
    '''
    try:
        mask = readfile.read(maskFile)[0]
        print 'read mask from file: '+maskFile
    except:
        maskFile = None
        print 'no mask input, use all pixels'

    if not outFile:
        outFile = os.path.splitext(inFile)[0]+'_coh.txt'

    atr = readfile.read_attribute(inFile)
    k = atr['FILE_TYPE']
    if not k in ['timeseries']:
        raise Exception('Only timeseries file is supported, input file is: '+k)
    range2phase = -4*np.pi/float(atr['WAVELENGTH'])

    h5 = h5py.File(inFile, 'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)

    f = open(outFile, 'w')
    f.write('# Date      spatial_average_coherence\n')
    for i in range(date_num):
        date = date_list[i]
        data = h5[k].get(date)[:]
        data = np.exp(1j*range2phase*data)
        if maskFile:
            data[mask==0] = np.nan
        coh = np.absolute(np.nanmean(data))
        msg = '%s    %.4f' % (date, coh)
        f.write(msg+'\n')
        print msg
    h5.close()
    f.close()
    print 'write to '+outFile

    return outFile


def normalize_timeseries(ts_mat, nanValue=0):
    '''Normalize timeseries of 2D matrix in time domain'''
    ts_mat -= np.min(ts_mat, 0)
    ts_mat *= 1/np.max(ts_mat, 0)
    ts_mat[np.isnan(ts_mat)] = 0
    return ts_mat

def normalize_timeseries_old(ts_mat, nanValue=0):
    ts_mat -= np.max(ts_mat, 0)
    ts_mat *= -1
    ts_mat /= np.max(ts_mat, 0)
    ts_mat[np.isnan(ts_mat)] = 1
    return ts_mat


############################################################
def update_file(outFile, inFile=None, overwrite=False, check_readable=True):
    '''Check whether to update outFile or not.
    return True if any of the following meets:
        1. if overwrite option set to True
        2. outFile is empty, e.g. None, []
        3. outFile is not existed
        4. outFile is not readable by readfile.read_attribute() when check_readable=True
        5. outFile is older than inFile, if inFile is not None
    Otherwise, return False.
    
    If inFile=None and outFile exists and readable, return False
    
    Inputs:
        inFile - string or list of string, input file(s)
    '''
    if overwrite:
        return True

    if not outFile or not os.path.isfile(outFile):
        return True

    if check_readable:
        try:
            atr = readfile.read_attribute(outFile)
        except:
            print outFile+' exists, but can not read, remove it.'
            rmCmd = 'rm '+outFile;  print rmCmd;  os.system(rmCmd)
            return True

    if inFile:
        # Convert string to list
        if isinstance(inFile, basestring):
            inFile = [inFile]

        # Check existance of each item
        fileList = list(inFile)
        for File in fileList:
            if not os.path.isfile(File):
                inFile.remove(File)

        # Check modification time
        if inFile:
            if any(os.path.getmtime(outFile) < os.path.getmtime(File) for File in inFile):
                return True
            else:
                print outFile+' exists and is newer than '+str(inFile)+', skip updating.'
                return False

    return False

def update_attribute_or_not(atr_new, atr_orig, update=False):
    '''Compare new attributes with exsiting ones'''
    for key in atr_new.keys():
        value = str(atr_new[key])
        if (key in atr_orig.keys() and value == str(atr_orig[key]) or\
            key not in atr_orig.keys() and value == 'None'):
            next
        else:
            update = True
    return update


def add_attribute(File, atr_new=dict()):
    '''Add/update input attribute into File
    Inputs:
        File - string, path/name of file
        atr_new - dict, attributes to be added/updated
                  if value is None, delete the item from input File attributes
    Output:
        File - string, path/name of updated file
    '''
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']

    # Compare new attributes with exsiting ones
    update = False
    if k in multi_dataset_hdf5_file+single_dataset_hdf5_file:
        update = update_attribute_or_not(atr_new, atr, update)
    elif k in multi_group_hdf5_file:
        h5 = h5py.File(File, 'r')
        epochList = h5[k].keys()
        for epoch in epochList:
            atr = h5[k][epoch].attrs
            update = update_attribute_or_not(atr_new, atr, update)
        h5.close()
    else:
        raise Exception('Un-recognized file type: '+k)

    if not update:
        print 'All updated (removed) attributes already exists (do not exists) and have the same value, skip update.'
        return File

    # Update attributes
    h5 = h5py.File(File,'r+')
    if k in multi_dataset_hdf5_file+single_dataset_hdf5_file:
        for key, value in atr_new.iteritems():
            # delete the item is new value is None
            if value == 'None':
                try: h5[k].attrs.pop(key)
                except: pass
            else:
                h5[k].attrs[key] = value
    elif k in multi_group_hdf5_file:
        epochList = h5[k].keys()
        for epoch in epochList:
            for key, value in atr_new.iteritems():
                if value == 'None':
                    try: h5[k][epoch].attrs.pop(key)
                    except: pass
                else:
                    h5[k][epoch].attrs[key] = value
    h5.close()
    return File


def check_parallel(file_num=1):
    '''Check parallel option based on pysar setting, file num and installed module'''
    enable_parallel = True

    # Disable parallel option for one input file
    if file_num <= 1:
        enable_parallel = False
        print 'parallel processing is diabled for one input file'
        return 1, enable_parallel, None, None

    # Check required python module
    try:
        from joblib import Parallel, delayed
    except:
        print 'Can not import joblib'
        print 'parallel is disabled.'
        enable_parallel = False
        return 1, enable_parallel, None, None

    # Find proper number of cores for parallel processing
    num_cores = min(multiprocessing.cpu_count(), file_num, pysar.parallel_num)
    if num_cores <= 1:
        enable_parallel = False
        print 'parallel processing is disabled because min of the following two numbers <= 1:'
        print 'available cpu number of the computer: '+str(multiprocessing.cpu_count())
        print 'pysar.__init__.py: parallel_num: '+str(pysar.parallel_num)
    else:
        print 'parallel processing using %d cores ...'%(num_cores)
    
    try:
        return num_cores, enable_parallel, Parallel, delayed
    except:
        return num_cores, enable_parallel, None, None


def perp_baseline_timeseries(atr, dimension=1):
    '''Calculate perpendicular baseline for each acquisition within timeseries
    Inputs:
        atr - dict, including the following PySAR attribute
              FILE_LENGTH
              P_BASELINE_TIMESERIES
              P_BASELINE_TOP_TIMESERIES (optional)
              P_BASELINE_BOTTOM_TIMESERIES (optional)
        dimension - int, choices = [0, 1]
                    0 for constant P_BASELINE in azimuth direction
                    1 for linear P_BASELINE in azimuth direction, for radar coord only
    Output:
        pbase - np.array, with shape = [date_num, 1] or [date_num, length]
    '''
    # return constant value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        print 'input file is geocoded, return constant P_BASELINE in azimuth direction within one interferogram'

    length = int(atr['FILE_LENGTH'])
    pbase_list = [float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]
    date_num = len(pbase_list)
    pbase = np.array(pbase_list).reshape(date_num, 1)

    if dimension > 0:
        try:
            pbase_top = np.array([float(i) for i in atr['P_BASELINE_TOP_TIMESERIES'].split()]).reshape(date_num, 1)
            pbase_bottom = np.array([float(i) for i in atr['P_BASELINE_BOTTOM_TIMESERIES'].split()]).reshape(date_num, 1)
            pbase = np.zeros((date_num, length))
            for i in range(date_num):
                pbase[i,:] = np.linspace(pbase_top[i], pbase_bottom[i], num=length, endpoint='FALSE')
        except:
            dimension = 0
            print 'Can not find P_BASELINE_TOP/BOTTOM_TIMESERIES in input attribute'
            print 'return constant P_BASELINE in azimuth direction for each acquisition instead'
    return pbase


def range_distance(atr, dimension=2):
    '''Calculate range distance from input attribute dict
    Inputs:
        atr - dict, including the following ROI_PAC attributes:
              STARTING_RANGE
              RANGE_PIXEL_SIZE
              FILE_LENGTH
              WIDTH
        dimension - int, choices = [0,1,2]
                    2 for 2d matrix, vary in range direction, constant in az direction, for radar coord only
                    1 for 1d matrix, in range direction, for radar coord file
                    0 for center value
    Output: np.array (0, 1 or 2 D) - range distance between antenna and ground target in meters
    '''
    # return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        print 'input file is geocoded, return center range distance for the whole area'
    
    near_range = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    
    far_range = near_range + dR*width
    center_range = (far_range + near_range)/2.0
    print 'center range : %.2f m' % (center_range)
    if dimension == 0:
        return np.array(center_range)
    
    print 'near   range : %.2f m' % (near_range)
    print 'far    range : %.2f m' % (far_range)
    range_x = np.linspace(near_range, far_range, num=width, endpoint='FALSE')
    if dimension == 1:
        return range_x
    else:
        range_xy = np.tile(range_x, (length, 1))
        return range_xy


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
    # Return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        print 'input file is geocoded, return center incident angle only'
    
    ## Read Attributes
    near_range = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    r  = float(atr['EARTH_RADIUS'])
    H  = float(atr['HEIGHT'])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    
    ## Calculation
    far_range = near_range+dR*width
    incidence_n = (np.pi - np.arccos((r**2+near_range**2-(r+H)**2)/(2*r*near_range))) * 180.0/np.pi
    incidence_f = (np.pi - np.arccos((r**2+ far_range**2-(r+H)**2)/(2*r*far_range))) * 180.0/np.pi
    incidence_c = (incidence_f+incidence_n)/2.0
    print 'center incidence angle : %.4f degree' % (incidence_c)
    if dimension == 0:
        return np.array(incidence_c)
    
    print 'near   incidence angle : %.4f degree' % (incidence_n)
    print 'far    incidence angle : %.4f degree' % (incidence_f)
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
        print 'read mask from file: '+maskFile
        mask = readfile.read(maskFile)[0]
        stack[mask==0] = np.nan
    
    return stack


def check_drop_ifgram(h5, atr, ifgram_list, print_message=True):
    '''Update ifgram_list based on 'drop_ifgram' attribute
    Inputs:
        h5          - HDF5 file object
        atr         - dict, file attribute
        ifgram_list - list of string, all group name existed in file
    Outputs:
        ifgram_list_out  - list of string, group name with drop_ifgram = 'yes'
        ifgram_list_drop - list of string, group name with drop_ifgram = 'no'
    '''
    # Return all interferogram list if 'drop_ifgram' do not exist
    if 'drop_ifgram' not in atr.keys():
        return ifgram_list

    ifgram_list_out = list(ifgram_list)
    k = atr['FILE_TYPE']
    for ifgram in ifgram_list:
        if h5[k][ifgram].attrs['drop_ifgram'] == 'yes':
            ifgram_list_out.remove(ifgram)
    
    if len(ifgram_list) > len(ifgram_list_out) and print_message:
        print "remove interferograms with 'drop_ifgram'='yes'"
    return ifgram_list_out


def nonzero_mask(File, outFile='mask.h5'):
    '''Generate mask file for non-zero value of input multi-group hdf5 file'''
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    
    mask = np.ones([length, width])
    
    h5 = h5py.File(File,'r')
    igramList = sorted(h5[k].keys())
    igramList = check_drop_ifgram(h5, atr, igramList)
    date12_list = ptime.list_ifgram2date12(igramList)
    prog_bar = ptime.progress_bar(maxValue=len(igramList), prefix='loading: ')
    for i in range(len(igramList)):
        igram = igramList[i]
        data = h5[k][igram].get(igram)[:]
        mask[data==0] = 0
        prog_bar.update(i+1, suffix=date12_list[i])
    prog_bar.close()

    atr['FILE_TYPE'] = 'mask'
    print 'writing >>> '+outFile
    writefile.write(mask, atr, outFile)
    
    return outFile


######################################################################################################

def get_spatial_average(File, maskFile=None, box=None, saveList=True):
    '''Get spatial average info from input File.
    Inputs:
        File     - string, path of HDF5 file or txt file
        maskFile - string, path of mask file, e.g. maskTempCoh.h5
        box      - 4-tuple defining the left, upper, right, and lower pixel coordinate
        saveList - bool, save (list of) mean value into text file
    outputs:
        mean_list - list of float, spatial average value of file
        date_list - list of string for date info
    Example:
        import pysar._pysar_utilities as ut
        mean_list, date_list = ut.get_spatial_average('coherence.h5', 'maskTempCoh.h5')
    '''
    suffix='_spatialAverage.txt'
    if File.endswith(suffix):
        print 'Input file is spatial average txt already, read it directly'
        txtFile = File

    else:
        txtFile = os.path.splitext(File)[0]+suffix
        if os.path.isfile(txtFile):
            print 'spatial average file exists: '+txtFile
            print 'read it directly, or delete it and re-run the script to re-calculate the list'
        else:
            print 'calculating spatial average from file: '+File
            if maskFile:
                mask = readfile.read(maskFile)[0]
                print 'read mask from file: '+maskFile
            else:
                mask = None
            mean_list = spatial_average(File, mask, box, saveList)

    # Read txt file
    txtContent = np.loadtxt(txtFile, dtype=str)
    mean_list = [float(i) for i in txtContent[:,1]]
    date_list = [i for i in txtContent[:,0]]
    return mean_list, date_list


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
        prog_bar = ptime.progress_bar(maxValue=epochNum, prefix='calculating: ')
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
            prog_bar.update(i+1)
        prog_bar.close()
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
        txtFile = os.path.splitext(File)[0]+'_spatialAverage.txt'
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
    epochList = check_drop_ifgram(h5file, atr, epochList)
    epochNum = len(epochList)

    # Calculation
    dMean = np.zeros((length,width))
    prog_bar = ptime.progress_bar(maxValue=epochNum, prefix='calculating: ')
    for i in range(epochNum):
        epoch = epochList[i]
        if k in multi_group_hdf5_file:
            d = h5file[k][epoch].get(epoch)[:]
        elif k in ['timeseries']:
            d = h5file[k].get(epoch)[:]
        else: print k+' type is not supported currently.'; sys.exit(1)
        dMean += d
        prog_bar.update(i+1)
    prog_bar.close()
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
def get_file_list(fileList, abspath=False):
    '''Get all existed files matching the input list of file pattern
    Example:
        fileList = get_file_list(['*velocity*.h5','timeseries*.h5'])
        fileList = get_file_list('timeseries*.h5')
    '''
    if not fileList:
        return []

    if isinstance(fileList, basestring):
        fileList = [fileList]

    fileListOut = []
    for i in range(len(fileList)):
        file0 = fileList[i]
        fileList0 = glob.glob(file0)
        fileListOut += sorted(list(set(fileList0) - set(fileListOut)))
    if abspath:
        fileListOut = [os.path.abspath(i) for i in fileListOut]
    return fileListOut


def mode (thelist):
    '''Find Mode (most common) item in the list'''
    if not thelist:
        return None
    if len(thelist) == 1:
        return thelist[0]

    counts = {}
    for item in thelist:
        counts[item] = counts.get(item, 0) + 1
    maxcount = 0
    maxitem  = None
    for k, v in counts.items():
        if v > maxcount:
            maxitem  = k
            maxcount = v

    if maxcount == 1:
        print "All values only appear once"
        return None
    elif counts.values().count(maxcount) > 1:
        print "List has multiple modes"
        return None
    else:
        return maxitem


######################################################################################################
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

    ##### Use geomap*.trans file for precious (pixel-level) coord conversion
    def get_trans_row_col4radar(az, rg, trans_az, trans_rg, x_factor=10, y_factor=10):
        mask_rg = np.multiply(trans_rg>=max(rg-x_factor,0.5), trans_rg<=rg+x_factor)
        mask_az = np.multiply(trans_az>=max(az-y_factor,0.5), trans_az<=az+y_factor)
        idx = np.where(np.multiply(mask_rg, mask_az))
        trans_row, trans_col = np.mean(idx,1)
        return trans_row, trans_col


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

            x_factor = 2*np.ceil(abs(lon_step)/rg_step)
            y_factor = 2*np.ceil(abs(lat_step)/az_step)
        else:
            x_factor = 10
            y_factor = 10

        trans_row = np.zeros(rg.shape)
        trans_col = np.zeros(rg.shape)
        if rg.size == 1:
            trans_row, trans_col = get_trans_row_col4radar(az, rg, trans_az, trans_rg, x_factor, y_factor)
        else:
            for i in range(rg.size):
                trans_row[i], trans_col[i] = get_trans_row_col4radar(az[i], rg[i], trans_az, trans_rg, x_factor, y_factor)

        lat = trans_row*lat_step_deg + lat_first
        lon = trans_col*lon_step_deg + lon_first
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
    ifgram_list = sorted(h5file[k[0]].keys())
    for ifgram in  ifgram_list:
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


######################################
def design_matrix(ifgramFile=None, date12_list=[]):
    '''Make the design matrix for the inversion based on date12_list.
    Input:
        ifgramFile  - string, name/path of interferograms file
        date12_list - list of string, date12 used in calculation in YYMMDD-YYMMDD format
                      use all date12 from ifgramFile if input is empty
    Outputs:
        A - 2D np.array in size (igram_num, date_num-1)
            representing date combination for each interferogram
        B - 2D np.array in size (igram_num, date_num-1)
            representing temporal baseline timeseries between master and slave date for each interferogram
    '''
    # Check Inputs
    if not date12_list:
        if ifgramFile:
            date12_list = pnet.get_date12_list(ifgramFile)
        else:
            raise ValueError

    # date12_list to date6_list
    m_dates = [i.split('-')[0] for i in date12_list]
    s_dates = [i.split('-')[1] for i in date12_list]
    date6_list = sorted(list(set(m_dates + s_dates)))
    tbase = np.array(ptime.date_list2tbase(date6_list)[0])
    date_num = len(date6_list)
    igram_num = len(date12_list)

    A = np.zeros((igram_num, date_num))
    B = np.zeros(np.shape(A))
    #t = np.zeros((igram_num, 2))
    for i in range(igram_num):
        m_idx, s_idx = [date6_list.index(j) for j in date12_list[i].split('-')]
        A[i, m_idx] = -1
        A[i, s_idx] = 1
        B[i, m_idx:s_idx] = tbase[m_idx+1:s_idx+1] - tbase[m_idx:s_idx]
        #t[i,:] = [tbase[m_idx], tbase[s_idx]]
    # Remove the 1st date assuming it's zero
    A = A[:,1:]
    B = B[:,:-1]

    return A,B


######################################
def timeseries_inversion(ifgramFile, timeseriesFile):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 
    
    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms 
      h5timeseries: hdf5 file with the output from the inversion
    '''
    total = time.time()

    # Basic Info
    atr = readfile.read_attribute(ifgramFile)
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    pixel_num = length * width

    h5ifgram = h5py.File(ifgramFile,'r')
    ifgram_list = sorted(h5ifgram['interferograms'].keys())
    ifgram_list = check_drop_ifgram(h5ifgram, atr, ifgram_list)
    ifgram_num = len(ifgram_list)

    # Convert ifgram_list to date12/8_list
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    m_dates = [i.split('-')[0] for i in date12_list]
    s_dates = [i.split('-')[1] for i in date12_list]
    date8_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    date_num = len(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    dt = np.diff(tbase_list).reshape((date_num-1,1))

    print 'number of interferograms : '+str(ifgram_num)
    print 'number of pixels in space: '+str(pixel_num)
    print 'number of acquisitions   : '+str(date_num)

    # Design matrix
    A,B = design_matrix(ifgramFile, date12_list)
    B_inv = np.array(np.linalg.pinv(B), np.float32)

    # Reference pixel in space
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print 'reference pixel in y/x: [%d, %d]'%(ref_y, ref_x)
    except:
        print 'ERROR: No ref_x/y found! Can not inverse interferograms without reference in space.'
        print 'run seed_data.py '+ifgramFile+' --mark-attribute for a quick referencing.'
        sys.exit(1)

    ##### Inversion Function
    def ts_inverse(dataLine, B_inv, dt, date_num):
        numPoint = dataLine.shape[1]
        tmp_rate = np.dot(B_inv, dataLine)
        defo1 = tmp_rate * np.tile(dt,(1,numPoint))
        defo0 = np.array([0.]*numPoint,np.float32)
        defo  = np.vstack((defo0, np.cumsum(defo1,axis=0)))
        return defo

    ##### Read Interferograms
    print 'Reading interferograms ...'
    data = np.zeros((ifgram_num,pixel_num), np.float32)
    prog_bar = ptime.progress_bar(maxValue=ifgram_num, prefix='loading: ')
    for j in range(ifgram_num):
        ifgram = ifgram_list[j]
        group = h5ifgram['interferograms'][ifgram]
        d = group.get(ifgram)[:]
        d -= d[ref_y, ref_x]
        data[j] = d.flatten(1)
        prog_bar.update(j+1, suffix=date12_list[j])
    h5ifgram.close()
    prog_bar.close()

    ##### Inversion
    print 'Inversing time series ...'
    dataPoint = np.zeros((ifgram_num,1),np.float32)
    dataLine  = np.zeros((ifgram_num,width),np.float32)
    tempDeformation = np.zeros((date_num,pixel_num),np.float32)

    prog_bar = ptime.progress_bar(maxValue=length, prefix='calculating: ')
    for i in range(length):
        dataLine = data[:,i*width:(i+1)*width]
        defoLine = ts_inverse(dataLine, B_inv, dt, date_num)
        tempDeformation[:,i*width:(i+1)*width] = defoLine
        prog_bar.update(i+1, every=length/100)
    prog_bar.close()
    del data

    ##### Time Series Data Preparation
    print 'converting phase to range'
    timeseries = np.zeros((date_num,length,width),np.float32)
    phase2range = -1*float(atr['WAVELENGTH'])/(4.*np.pi)
    for i in range(date_num):
        timeseries[i] = tempDeformation[i].reshape(width,length).T
        timeseries[i] *= phase2range
    del tempDeformation
  
    ##### Output Time Series File
    print 'writing >>> '+timeseriesFile
    print 'number of dates: '+str(date_num)
    h5timeseries = h5py.File(timeseriesFile,'w')
    group = h5timeseries.create_group('timeseries')
    prog_bar = ptime.progress_bar(maxValue=date_num, prefix='writing: ')
    for i in range(date_num):
        date = date8_list[i]
        dset = group.create_dataset(date, data=timeseries[i], compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    ## Attributes
    print 'calculating perpendicular baseline timeseries'
    pbase, pbase_top, pbase_bottom = perp_baseline_ifgram2timeseries(ifgramFile, ifgram_list)
    # convert np.array into string with each item separated by white space
    pbase = str(pbase.tolist()).translate(None,'[],')
    pbase_top = str(pbase_top.tolist()).translate(None,'[],')
    pbase_bottom = str(pbase_bottom.tolist()).translate(None,'[],')
    atr['P_BASELINE_TIMESERIES'] = pbase
    atr['P_BASELINE_TOP_TIMESERIES'] = pbase_top
    atr['P_BASELINE_BOTTOM_TIMESERIES'] = pbase_bottom
    atr['ref_date'] = date8_list[0]
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5timeseries.close()
    print 'Time series inversion took ' + str(time.time()-total) +' secs\nDone.'
    return timeseriesFile

    
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
    ifgram_list = h5flat['interferograms'].keys()
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset=h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0])
    data = dset[0:dset.shape[0],0:dset.shape[1]] 
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print 'Reading in the interferograms'
    #print ifgram_num,pixel_num
    print 'number of interferograms: '+str(ifgram_num)
    print 'number of pixels: '+str(pixel_num)
    pixel_num_step = int(pixel_num/10)
  
    data = np.zeros((ifgram_num,pixel_num),np.float32)
    for ni in range(ifgram_num):
        dset=h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0],0:dset.shape[1]]
        #print np.shape(d)

    del d
    dataPoint = np.zeros((ifgram_num,1),np.float32)
    modelDimension = np.shape(B)[1]
    tempDeformation = np.zeros((date_num,pixel_num),np.float32)
    for ni in range(pixel_num):
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
        #if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,pixel_num)
        if not np.remainder(ni,pixel_num_step):
            print 'Processing point: %8d of %8d, %3d' % (ni,pixel_num,(10*ni/pixel_num_step))+'%'
    del data
    timeseries = np.zeros((date_num,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = tempDeformation[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del tempDeformation
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.iteritems():
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
    ifgram_list = h5flat['interferograms'].keys()
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset=h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0]) 
    data = dset[0:dset.shape[0],0:dset.shape[1]] 
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print 'Reading in the interferograms'
    print ifgram_num,pixel_num
  
    #data = np.zeros((ifgram_num,pixel_num),np.float32)
    data = np.zeros((ifgram_num,pixel_num))
    for ni in range(ifgram_num):
        dset=h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0],0:dset.shape[1]]
        #print np.shape(d)
    
        data[ni] = d.flatten(1)
    del d
    dataPoint = np.zeros((ifgram_num,1),np.float32)
    modelDimension = np.shape(B)[1]
    tempDeformation = np.zeros((date_num,pixel_num),np.float32)
    print data.shape
    DataL1=matrix(data)
    L1ORL2=np.ones((pixel_num,1))
    for ni in range(pixel_num):
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
        if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,pixel_num)
    del data
    timeseries = np.zeros((date_num,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = tempDeformation[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del tempDeformation
    L1ORL2=np.reshape(L1ORL2,(np.shape(dset)[1],np.shape(dset)[0])).T
    
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.iteritems():
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


def perp_baseline_ifgram2timeseries(ifgramFile, ifgram_list=[]):
    '''Calculate perpendicular baseline timeseries from input interferograms file
    Input:
        ifgramFile - string, file name/path of interferograms file
        ifgram_list - list of string, group name that is used for calculation
                      use all if it's empty
    Outputs:
        pbase        - 1D np.array, P_BASELINE_TIMESERIES
        pbase_top    - 1D np.array, P_BASELINE_TOP_TIMESERIES
        pbase_bottom - 1D np.array, P_BASELINE_BOTTOM_TIMESERIES
    '''
    k = readfile.read_attribute(ifgramFile)['FILE_TYPE']
    h5file = h5py.File(ifgramFile,'r')

    if not ifgram_list:
        ifgram_list = sorted(h5file[k].keys())

    # P_BASELINE of all interferograms
    pbase_ifgram = []
    pbase_top_ifgram = []
    pbase_bottom_ifgram = []
    for ifgram in ifgram_list:
        pbase_top    = float(h5file[k][ifgram].attrs['P_BASELINE_TOP_HDR'])
        pbase_bottom = float(h5file[k][ifgram].attrs['P_BASELINE_BOTTOM_HDR'])
        pbase_ifgram.append((pbase_bottom+pbase_top)/2.0)
        pbase_top_ifgram.append(pbase_top)
        pbase_bottom_ifgram.append(pbase_bottom)
    h5file.close()

    # Temporal baseline velocity
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    m_dates = [i.split('-')[0] for i in date12_list]
    s_dates = [i.split('-')[1] for i in date12_list]
    date8_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    tbase_v = np.diff(tbase_list)

    A,B = design_matrix(ifgramFile, date12_list)
    B_inv = np.linalg.pinv(B)

    pbase_rate        = np.dot(B_inv, pbase_ifgram)
    pbase_top_rate    = np.dot(B_inv, pbase_top_ifgram)
    pbase_bottom_rate = np.dot(B_inv, pbase_bottom_ifgram)

    zero = np.array([0.],np.float32)
    pbase        = np.concatenate((zero,np.cumsum([pbase_rate*tbase_v])))
    pbase_top    = np.concatenate((zero,np.cumsum([pbase_top_rate*tbase_v])))
    pbase_bottom = np.concatenate((zero,np.cumsum([pbase_bottom_rate*tbase_v])))

    return pbase, pbase_top, pbase_bottom


def dBh_dBv_timeseries(ifgramFile):
    h5file = h5py.File(ifgramFile)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    igramList = h5file[k[0]].keys()
    dBh_igram=[]
    dBv_igram=[]
    for igram in igramList:
        dBh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_RATE_HDR']))
        dBv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_RATE_HDR']))

    A,B=design_matrix(ifgramFile)
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

def Bh_Bv_timeseries(ifgramFile):
    h5file = h5py.File(ifgramFile)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    igramList = h5file[k[0]].keys()
    Bh_igram=[]
    Bv_igram=[]
    for igram in igramList:
        Bh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_TOP_HDR']))
        Bv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_TOP_HDR']))
  
    A,B=design_matrix(ifgramFile)
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
        prog_bar = ptime.progress_bar(maxValue=epochNum, prefix='calculating: ')
        for i in range(epochNum):
            epoch = epochList[i]
            if k == 'timeseries':  data = h5file[k].get(epoch)[:]
            else:                  data = h5file[k][epoch].get(epoch)[:]
            stack += data
            prog_bar.update(i+1)
        prog_bar.close()
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
    ifgram_list = h5file['interferograms'].keys()
    h5curlfile=h5py.File(curlfile,'w')
    gg = h5curlfile.create_group('interferograms')
    lcurls=np.shape(curls)[0]
    for i in range(lcurls):
        d1=h5file['interferograms'][ifgram_list[curls[i,0]]].get(ifgram_list[curls[i,0]])
        d2=h5file['interferograms'][ifgram_list[curls[i,1]]].get(ifgram_list[curls[i,1]])
        d3=h5file['interferograms'][ifgram_list[curls[i,2]]].get(ifgram_list[curls[i,2]])
        data1=d1[0:d1.shape[0],0:d1.shape[1]]
        data2=d2[0:d2.shape[0],0:d2.shape[1]]
        data3=d3[0:d3.shape[0],0:d3.shape[1]]
 
        print i
        group = gg.create_group(Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2])
        dset = group.create_dataset(Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2],\
                                    data=data1+data3-data2, compression='gzip')
        for key, value in h5file['interferograms'][ifgram_list[curls[i,0]]].attrs.iteritems():
            group.attrs[key] = value
 
    h5curlfile.close()


