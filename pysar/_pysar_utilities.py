#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
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


###############################################################################
def touch(fname_list, times=None):
    '''python equivalent function to Unix utily - touch
    It sets the modification and access times of files to the current time of day.
    If the file doesn't exist, it is created with default permissions.
    Inputs/Output:
        fname_list - string / list of string
    '''
    if not fname_list:
        return None

    if isinstance(fname_list, basestring):
        fname_list = [fname_list]

    fname_list = filter(lambda x: x!=None, fname_list)
    for fname in fname_list:
        with open(fname, 'a'):
            os.utime(fname, times)
            print 'touch '+fname

    if len(fname_list) == 1:
        fname_list = fname_list[0]
    return fname_list


def get_lookup_file(filePattern=None, abspath=False):
    '''Find lookup table file with/without input file pattern'''
    ##Files exists
    if not filePattern:
        filePattern = ['geometryRadar.h5',\
                       'geometryGeo_tight.h5', 'geometryGeo.h5',\
                       'geomap*lks_tight.trans', 'geomap*lks.trans',\
                       'sim*_tight.UTM_TO_RDC', 'sim*.UTM_TO_RDC']
    existFiles = []
    try:
        existFiles = get_file_list(filePattern)
    except:
        print 'ERROR: No geometry / lookup table file found!'
        print 'It should be like:'
        print filePattern
        return None

    ##Files with lookup table info
    lookupFile = None
    for fname in existFiles:
        atr = readfile.read_attribute(fname)
        if 'Y_FIRST' in atr.keys():
            epoch2check = 'rangeCoord'
        else:
            epoch2check = 'latitude'
        try:
            dset = readfile.read(fname, epoch=epoch2check, print_msg=False)[0]
            lookupFile = fname
            break
        except:
            pass

    if not lookupFile:
        print 'No lookup table info range/lat found in files.'
        return None

    if abspath:
        lookupFile = os.path.abspath(lookupFile)
    return lookupFile


def get_dem_file(coordType='radar', filePattern=None, abspath=False):
    '''Find DEM file with/without input file pattern'''
    ##Files exists
    if not filePattern:
        if coordType == 'radar':
            filePattern = ['geometryRadar.h5',\
                           'demRadar.h5',\
                           'radar*.hgt']
        else:
            filePattern = ['geometryGeo_tight.h5', 'geometryGeo.h5',\
                           'demRadar.h5', 'demGeo.h5',\
                           '*.dem', '*.dem.wgs84']
    existFiles = []
    try:
        existFiles = get_file_list(filePattern)
    except:
        print 'ERROR: No DEM file found!'
        print 'It should be like:'
        print filePattern
        return None

    ##Files with lookup table info
    demFiles = []
    for fname in existFiles:
        atr = readfile.read_attribute(fname)
        try:
            dset = readfile.read(fname, epoch='height')[0]
            demFiles.append(fname)
        except: pass

    if not demFiles:
        print 'No height info found in files.'
        return None

    demFile = demFiles[0]
    if abspath:
        demFile = os.path.abspath(demFile)
    return demFile


def check_loaded_dataset(work_dir='./', inps=None, print_msg=True):
    '''Check the result of loading data for the following two rules:
        1. file existance
        2. file attribute readability

    If inps is valid/not_empty: return updated inps;
    Otherwise, return True/False if all recommended file are loaded and readably or not

    Inputs:
        work_dir : string, PySAR working directory
        inps     : Namespace, optional, variable for pysarApp.py. Not needed for check loading result.
    Outputs:
        load_complete  : bool, complete loading or not
        ifgram_file    : string, file name/path of unwrapped interferograms
        coherence_file : string, file name/path of spatial coherence
        dem_file_radar : string, file name/path of DEM file in radara coord (for interferograms in radar coord)
        dem_file_geo   : string, file name/path of DEM file in geo coord
        lookup_file    : string, file name/path of lookup table file (for interferograms in radar coord)
    Example:
        from pysar.pysarApp import check_loaded_dataset
        True = check_loaded_dataset($SCRATCHDIR+'/SinabungT495F50AlosA/PYSAR') #if True, PROCESS, SLC folder could be removed.
        inps = check_loaded_dataset(inps.work_dir, inps)
    '''
    ##### Find file name/path of all loaded files
    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    if inps:
        inps.ifgram_file    = None
        inps.coherence_file = None
        inps.dem_radar_file = None
        inps.dem_geo_file   = None
        inps.lookup_file    = None

    # Required files - 1. unwrapped interferograms
    file_list = [work_dir+'/Modified_unwrapIfgram.h5',\
                 work_dir+'/unwrapIfgram.h5',\
                 work_dir+'/Modified_LoadedData.h5',\
                 work_dir+'/LoadedData.h5']
    ifgram_file = is_file_exist(file_list, abspath=True)

    if not ifgram_file:
        if inps:
            return inps
        else:
            return False
    else:
        atr = readfile.read_attribute(ifgram_file)

    if print_msg:
        print 'Loaded dataset are processed by %s InSAR software' % atr['INSAR_PROCESSOR']

    if 'X_FIRST' in atr.keys():
        geocoded = True
        if print_msg:
            print 'Loaded dataset are in geo coordinates'
    else:
        geocoded = False
        if print_msg:
            print 'Loaded dataset are in radar coordinates'

    if print_msg:
        print 'Unwrapped interferograms: '+ifgram_file

    # Recommended files (None if not found)
    # 2. Spatial coherence for each interferogram
    file_list = [work_dir+'/Modified_coherence.h5',\
                 work_dir+'/coherence.h5',\
                 work_dir+'/Modified_Coherence.h5',\
                 work_dir+'/Coherence.h5']
    coherence_file = is_file_exist(file_list, abspath=True)
    if print_msg:
        if coherence_file:
            print 'Spatial       coherences: '+coherence_file
        else:
            print 'WARNING: No coherences file found. Cannot use coherence-based network modification without it.'
            print "It's supposed to be like: "+str(file_list)

    # 3. DEM in radar coord
    dem_radar_file = get_dem_file(coordType='radar', abspath=True)
    if print_msg:
        if dem_radar_file:
            print 'DEM in radar coordinates: '+dem_radar_file
        elif not geocoded:
            print 'WARNING: No DEM file in radar coord found.'
            print "It's supposed to be like: "+str(file_list)

    # 4. DEM in geo coord
    dem_geo_file = get_dem_file(coordType='geo', abspath=True)
    if print_msg:
        if dem_geo_file:
            print 'DEM in geo   coordinates: '+dem_geo_file
        else:
            print 'WARNING: No DEM file in geo coord found.'
            print "It's supposed to be like: "+str(file_list)

    # 5. Lookup table file for geocoding
    lookup_file = get_lookup_file(inps.lookup_file, abspath=True)
    if print_msg:
        if lookup_file:
            print 'Lookup table        file: '+lookup_file
        elif not geocoded:
            print 'No lookup file found! Can not geocode without it!'
            print "It's supposed to be like: "+str(file_list)

    ##### Update namespace inps if inputed
    load_complete = True
    if None in [ifgram_file, coherence_file, dem_geo_file]:
        load_complete = False
    if not geocoded and None in [dem_radar_file, lookup_file]:
        load_complete = False
    if load_complete and print_msg:
        print '-----------------------------------------------------------------------------------'
        print 'All data needed found/loaded/copied. Processed 2-pass InSAR data can be removed.'
    print '-----------------------------------------------------------------------------------'

    if inps:
        inps.ifgram_file    = ifgram_file
        inps.coherence_file = coherence_file
        inps.dem_radar_file = dem_radar_file
        inps.dem_geo_file   = dem_geo_file
        inps.lookup_file    = lookup_file
        return inps

    ##### Check 
    else:
        return load_complete


def is_file_exist(file_list, abspath=True):
    '''Check if any file in the file list 1) exists and 2) readable
    Inputs:
        file_list : list of string, file name with/without wildcards
        abspath   : bool, return absolute file name/path or not
    Output:
        file_path : string, found file name/path; None if not.
    '''
    try:
        file_path = get_file_list(file_list, abspath=abspath)[0]
        atr_temp = readfile.read_attribute(file_path)
    except:
        file_path = None
    return file_path


def four_corners(atr):
    '''Return 4 corners lat/lon'''
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    west  = float(atr['X_FIRST'])
    north = float(atr['Y_FIRST'])
    south = north + lat_step*length
    east  = west  + lon_step*width

    return west, east, south, north


def circle_index(atr,circle_par):
    '''Return Index of Elements within a Circle centered at input pixel
    Inputs: atr : dictionary
                containging the following attributes:
                WIDT
                FILE_LENGTH
            circle_par : string in the format of 'y,x,radius'
                i.e. '200,300,20'          for radar coord
                     '31.0214,130.5699,20' for geo   coord
    Output: idx : 2D np.array in bool type
                mask matrix for those pixel falling into the circle defined by circle_par
    Examples: idx_mat = ut.circle_index(atr, '200,300,20')
              idx_mat = ut.circle_index(atr, '31.0214,130.5699,20')
    '''

    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    if type(circle_par) == tuple:
        cir_par = circle_par
    elif type(circle_par) == list:
        cir_par = circle_par
    else:
        cir_par = circle_par.replace(',',' ').split()

    try:
        c_y    = int(cir_par[0])
        c_x    = int(cir_par[1])
        radius = int(float(cir_par[2]))
    except:
        try:
            c_lat  = float(cir_par[0])
            c_lon  = float(cir_par[1])
            radius = int(float(cir_par[2]))
            c_y = np.rint((c_lat-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
            c_x = np.rint((c_lon-float(atr['X_FIRST']))/float(atr['X_STEP']))
        except:
            print '\nERROR: Unrecognized circle index format: '+circle_par
            print 'Supported format:'
            print '--circle 200,300,20            for radar coord input'
            print '--circle 31.0214,130.5699,20   for geo   coord input\n'
            return 0

    y,x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx


def update_template_file(template_file, extra_dict):
    '''Update option value in template_file with value from input extra_dict'''
    extra_key_list = extra_dict.keys()

    ## Compare and skip updating template_file is no new option value found.
    update = False
    orig_dict = readfile.read_template(template_file)
    for key, value in orig_dict.iteritems():
        if key in extra_key_list and extra_dict[key] != value:
            update = True
    if not update:
        print 'No new option value found, skip updating '+template_file
        return template_file

    ## Update template_file with new value from extra_dict
    tmp_file = template_file+'.tmp'
    f_orig = open(template_file, 'r')
    f_tmp = open(tmp_file, 'w')
    for line in f_orig:
        line = line.strip()
        c = [i.strip() for i in line.split('=', 1)]
        if not line.startswith('%') and not line.startswith('#') and len(c) > 1:
            key = c[0]
            value = str.replace(c[1],'\n','').split("#")[0].strip()
            if key in extra_key_list and extra_dict[key] != value:
                line = line.replace(value, extra_dict[key], 1)
                print '    '+key+': '+value+' --> '+extra_dict[key]
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
        mask = readfile.read(maskFile, epoch='mask')[0]
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
    rms_file = os.path.dirname(os.path.abspath(deramp_file))+'/rms_'+os.path.splitext(deramp_file)[0]+'.txt'

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
    print 'read timeseries residual RMS from file: '+rms_file
    rms_fileContent = np.loadtxt(rms_file, dtype=str)
    rms_list = rms_fileContent[:,1].astype(np.float).tolist()
    date_list = list(rms_fileContent[:,0]) 
    
    return rms_list, date_list


def timeseries_rms(inFile, maskFile='maskTempCoh.h5', outFile=None, dimension=2):
    '''Calculate the Root Mean Square for each epoch of input timeseries file
    and output result to a text file.
    '''
    try:
        mask = readfile.read(maskFile, epoch='mask')[0]
        print 'read mask from file: '+maskFile
    except:
        maskFile = None
        print 'no mask input, use all pixels'

    if not outFile:
        outFile = os.path.dirname(os.path.abspath(inFile))+'/rms_'+os.path.splitext(inFile)[0]+'.txt'

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
        mask = readfile.read(maskFile, epoch='mask')[0]
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
    '''Check whether to update outFile/outDir or not.
    return True if any of the following meets:
        1. if overwrite option set to True
        2. outFile is empty, e.g. None, []
        3. outFile is not existed
        4. outFile is not readable by readfile.read_attribute() when check_readable=True
        5. outFile is older than inFile, if inFile is not None
    Otherwise, return False.
    
    If inFile=None and outFile exists and readable, return False
    
    Inputs:
        inFile - string or list of string, input file(s)/directories
    Output:
        True/False - bool, whether to update output file or not
    Example:
        if ut.update_file('timeseries_ECMWF_demErr.h5', 'timeseries_ECMWF.h5'):
        if ut.update_file('exclude_date.txt', ['timeseries_ECMWF_demErrInvResid.h5','maskTempCoh.h5','pysar_template.txt'],\
                          check_readable=False):
    '''
    if overwrite:
        return True

    if not outFile or (not os.path.isfile(outFile) and not os.path.isdir(outFile)):
        return True

    if check_readable:
        try:
            atr = readfile.read_attribute(outFile)
        except:
            print outFile+' exists, but can not read, remove it.'
            rmCmd = 'rm '+outFile;  print rmCmd;  os.system(rmCmd)
            return True

    if inFile:
        inFile = get_file_list(inFile)

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


def check_parallel(file_num=1, print_msg=True):
    '''Check parallel option based on pysar setting, file num and installed module
    Examples:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(1000)
    '''
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
    elif print_msg:
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
    if dimension > 0 and 'Y_FIRST' in atr.keys():
        dimension = 0
        print 'file is in geo coordinates, return constant P_BASELINE for one interferogram'
    if dimension > 0 and any(i not in atr.keys() for i in ['P_BASELINE_TOP_TIMESERIES',\
                                                           'P_BASELINE_BOTTOM_TIMESERIES']):
        dimension = 0
        print 'No P_BASELINE_TOP/BOTTOM_TIMESERIES attributes found, return constant P_BASELINE for one interferogram'

    pbase_center = np.array([float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]).reshape(-1,1)
    if dimension == 0:
        pbase = pbase_center
    elif dimension == 1:
        pbase_top    = np.array([float(i) for i in atr['P_BASELINE_TOP_TIMESERIES'].split()]).reshape(-1, 1)
        pbase_bottom = np.array([float(i) for i in atr['P_BASELINE_BOTTOM_TIMESERIES'].split()]).reshape(-1, 1)
        length = int(atr['FILE_LENGTH'])
        date_num = pbase_center.shape[0]
        pbase = np.zeros((date_num, length))
        for i in range(date_num):
            pbase[i,:] = np.linspace(pbase_top[i], pbase_bottom[i], num=length, endpoint='FALSE')
    else:
        raise ValueError('Input pbase dimension: %s, only support 0 and 1.' % (str(dimension)))

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


def incidence_angle(atr, dimension=2, print_msg=True):
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
        if print_msg:
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
    if print_msg:
        print 'center incidence angle : %.4f degree' % (incidence_c)
    if dimension == 0:
        return np.array(incidence_c)

    if print_msg:
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


def check_drop_ifgram(h5, print_msg=True):
    '''Update ifgram_list based on 'drop_ifgram' attribute
    Input:
        h5          - HDF5 file object
    Output:
        dsListOut  - list of string, group name with drop_ifgram = 'yes'
    Example:
        h5 = h5py.File('unwrapIfgram.h5','r')
        ifgram_list = ut.check_drop_ifgram(h5)
    '''
    # Return all interferogram list if 'drop_ifgram' do not exist
    k = h5.keys()[0]
    dsList = sorted(h5[k].keys())
    atr = h5[k][dsList[0]].attrs
    if 'drop_ifgram' not in atr.keys():
        return dsList

    dsListOut = list(dsList)
    for ds in dsList:
        if h5[k][ds].attrs['drop_ifgram'] == 'yes':
            dsListOut.remove(ds)

    if len(dsList) > len(dsListOut) and print_msg:
        print "remove interferograms with 'drop_ifgram'='yes'"
    return dsListOut


def nonzero_mask(File, outFile='mask.h5'):
    '''Generate mask file for non-zero value of input multi-group hdf5 file'''
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    
    mask = np.ones([length, width])
    
    h5 = h5py.File(File,'r')
    igramList = sorted(h5[k].keys())
    igramList = check_drop_ifgram(h5)
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
def spatial_average(File, maskFile=None, box=None, saveList=False, checkAoi=True):
    '''Read/Calculate Spatial Average of input file.

    If input file is text file, read it directly;
    If input file is data matrix file:
        If corresponding text file exists with the same mask file/AOI info, read it directly;
        Otherwise, calculate it from data file.

        Only non-nan pixel is considered.
    Input:
        File     : string, path of input file
        maskFile : string, path of mask file, e.g. maskTempCoh.h5
        box      : 4-tuple defining the left, upper, right, and lower pixel coordinate
        saveList : bool, save (list of) mean value into text file
    Output:
        mean_list : list for float, average value in space for each epoch of input file
        date_list : list of string for date info
                    date12_list, e.g. 101120-110220, for interferograms/coherence
                    date8_list, e.g. 20101120, for timeseries
                    file name, e.g. velocity.h5, for all the other file types
    Example:
        mean_list = spatial_average('coherence.h5')[0]
        ref_list  = spatial_average('unwrapIfgram.h5', box=(100,200,101,201))[0]
        mean_list, date12_list = spatial_average('coherence.h5', 'maskTempCoh.h5', saveList=True)
        
        stack = ut.get_file_stack('unwrapIfgram.h5', 'mask.h5')
        mask = ~np.isnan(stack)
        ref_list = ut.spatial_average('unwrapIfgram.h5', mask, (100,200,101,201))
    '''
    suffix='_spatialAverage.txt'
    if File.endswith(suffix):
        print 'Input file is spatial average txt already, read it directly'
        txtFile = File
        txtContent = np.loadtxt(txtFile, dtype=str)
        mean_list = [float(i) for i in txtContent[:,1]]
        date_list = [i for i in txtContent[:,0]]
        return mean_list, date_list


    # Baic File Info
    atr  = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    if not box:
        box = (0,0,width,length)

    # Convert input mask argument (maskFile) to mask file name (maskFile) and matrix (mask)
    if maskFile is None:
        maskFile = None
        mask = None
        print 'no mask input, use all pixels available'
    elif type(maskFile) is str:
        print 'mask from file: '+maskFile
        mask = readfile.read(maskFile, epoch='mask')[0]
        mask = mask[box[1]:box[3],box[0]:box[2]]
    elif type(maskFile) is np.ndarray:
        mask = maskFile
        mask = mask[box[1]:box[3],box[0]:box[2]]
        maskFile = 'np.ndarray matrix'
        print 'mask from input matrix'
    else:
        print 'Unsupported mask input format: '+str(type(maskFile))
        return None, None

    # Read existing txt file only if 1) data file is older AND 2) same AOI
    read_txt = False
    txtFile = os.path.splitext(File)[0]+suffix
    file_line = '# Data file: %s\n' % os.path.basename(str(File))
    mask_line = '# Mask file: %s\n' % os.path.basename(str(maskFile))
    aoi_line = '# AOI box: %s\n' % str(box)

    try:
        # Read AOI line from existing txt file
        fl = open(txtFile,'r')
        lines = fl.readlines()
        fl.close()
        if checkAoi:
            try:    aoi_line_orig = [i for i in lines if '# AOI box:' in i][0]
            except: aoi_line_orig = ''
        else:
            aoi_line_orig = aoi_line
        try:    mask_line_orig = [i for i in lines if '# Mask file:' in i][0]
        except: mask_line_orig = ''
        if (aoi_line_orig == aoi_line \
            and mask_line_orig == mask_line \
            and not update_file(txtFile, [File, maskFile], check_readable=False)):
            read_txt = True
    except: pass

    if read_txt:
        print txtFile+' already exists, read it directly'
        txtContent = np.loadtxt(txtFile, dtype=str)
        mean_list = [float(i) for i in txtContent[:,1]]
        date_list = [i for i in txtContent[:,0]]
        return mean_list, date_list


    # Calculate mean coherence list
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        print 'calculating spatial average of file: '+os.path.basename(File)
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum  = len(epochList)

        mean_list   = []
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
                mean_list.append(np.nanmean(data))
            prog_bar.update(i+1)
        prog_bar.close()
        del data
    else:
        data,atr = readfile.read(File, box)
        if not mask is None:
            data[mask==0] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_list = [np.nanmean(data)]

    # Get date/pair list
    if k in multi_group_hdf5_file:
        date_list = pnet.get_date12_list(File)
        # Temp/Perp Baseline
        m_dates = [date12.split('-')[0] for date12 in date_list]
        s_dates = [date12.split('-')[1] for date12 in date_list]
        date6_list = ptime.yymmdd(sorted(list(set(m_dates + s_dates))))
        tbase_ts_list = ptime.date_list2tbase(date6_list)[0]
        tbase_list = []
        pbase_list = []
        for i in range(epochNum):
            ifgram = epochList[i]
            pbase_top    = float(h5file[k][ifgram].attrs['P_BASELINE_TOP_HDR'])
            pbase_bottom = float(h5file[k][ifgram].attrs['P_BASELINE_BOTTOM_HDR'])
            pbase = (pbase_bottom+pbase_top)/2.0
            pbase_list.append(pbase)

            m_idx = date6_list.index(m_dates[i])
            s_idx = date6_list.index(s_dates[i])
            tbase = tbase_ts_list[s_idx] - tbase_ts_list[m_idx]
            tbase_list.append(tbase)

    elif k in multi_dataset_hdf5_file:
        date_list = epochList
    else:
        date_list = [os.path.basename(File)]

    try: h5file.close()
    except: pass


    # Write mean coherence list into text file
    if saveList:
        print 'write average coherence in space into text file: '+txtFile
        fl = open(txtFile, 'w')
        # Write comments
        fl.write(file_line)
        fl.write(mask_line)
        fl.write(aoi_line)

        # Write data list
        line_num = len(date_list)
        if k in multi_group_hdf5_file:
            fl.write('#   DATE12        Mean      Btemp/days   Bperp/m\n')
            for i in range(line_num):
                line = '%s    %.4f    %8.0f    %8.1f\n' % (date_list[i], mean_list[i], tbase_list[i], pbase_list[i])
                fl.write(line)
        else:
            fl.write('#   DATE12        Mean\n')
            for i in range(line_num):
                line = '%s    %.4f\n' % (date_list[i], mean_list[i])
                fl.write(line)
        fl.close()

    if len(mean_list) == 1:
        mean_list = mean_list[0]
        date_list = date_list[0]
    return mean_list, date_list


def temporal_average(File, outFile=None):
    '''Calculate temporal average.'''
    # Input File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    h5file = h5py.File(File)
    epochList = sorted(h5file[k].keys())
    epochList = check_drop_ifgram(h5file)
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
def get_file_list(fileList, abspath=False, coord=None):
    '''Get all existed files matching the input list of file pattern
    Inputs:
        fileList - string or list of string, input file/directory pattern
        abspath  - bool, return absolute path or not
        coord    - string, return files with specific coordinate type: geo or radar
                   if none, skip the checking and return all files
    Output:
        fileListOut - list of string, existed file path/name, [] if not existed
    Example:
        fileList = get_file_list(['*velocity*.h5','timeseries*.h5'])
        fileList = get_file_list('timeseries*.h5')
    '''
    if not fileList:
        return []

    if isinstance(fileList, basestring):
        fileList = [fileList]

    # Get rid of None element
    fileList = filter(lambda x: x!=None, fileList)
    fileListOut = []
    for i in range(len(fileList)):
        file0 = fileList[i]
        fileList0 = glob.glob(file0)
        fileListOut += sorted(list(set(fileList0) - set(fileListOut)))

    if abspath:
        fileListOut = [os.path.abspath(i) for i in fileListOut]

    if coord is not None:
        fileListOutBk = list(fileListOut)
        for fname in fileListOutBk:
            atr = readfile.read_attribute(fname)
            if coord in ['geo']:
                if 'Y_FIRST' not in atr.keys():
                    fileListOut.remove(fname)
            elif coord in ['radar','rdr','rdc']:
                if 'Y_FIRST' in atr.keys():
                    fileListOut.remove(fname)
            else:
                raise ValueError('Input coord type: '+str(coord)+'\n. Only support geo, radar, rdr, rdc inputs.')

    return fileListOut


##################################################################
def check_file_size(fname_list, mode_width=None, mode_length=None):
    '''Check file size in the list of files, and drop those not in the same size with majority.'''
    # If input file list is empty
    if not fname_list:
        return fname_list, None, None

    # Read Width/Length list
    width_list = []
    length_list = []
    for fname in fname_list:
        atr = readfile.read_attribute(fname)
        width_list.append(atr['WIDTH'])
        length_list.append(atr['FILE_LENGTH'])

    # Mode of Width and Length
    if not mode_width:
        mode_width = mode(width_list)
    if not mode_length:
        mode_length = mode(length_list)
    
    # Update Input List
    fname_list_out = list(fname_list)
    if width_list.count(mode_width)!=len(width_list) or length_list.count(mode_length)!=len(length_list):
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print 'WARNING: Some files may have the wrong dimensions!'
        print 'All files should have the same size.'
        print 'The width and length of the majority of files are: %s, %s' % (mode_width, mode_length)
        print 'But the following files have different dimensions and thus will not be loaded:'
        for i in range(len(fname_list)):
            if width_list[i] != mode_width or length_list[i] != mode_length:
                print '%s    width: %s  length: %s' % (fname_list[i], width_list[i], length_list[i])
                fname_list_out.remove(fname_list[i])
        print '\nNumber of files left: '+str(len(fname_list_out))
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    return fname_list_out, mode_width, mode_length


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
def range_ground_resolution(atr, print_msg=False):
    '''Get range resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord'''
    if 'X_FIRST' in atr.keys():
        print 'Input file is in geo coord, no range resolution info.'
        return
    inc_angle = incidence_angle(atr, 0, print_msg)
    rg_step = float(atr['RANGE_PIXEL_SIZE'])/np.sin(inc_angle/180.0*np.pi)
    return rg_step

def azimuth_ground_resolution(atr):
    '''Get azimuth resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord'''
    if 'X_FIRST' in atr.keys():
        print 'Input file is in geo coord, no azimuth resolution info.'
        return
    try:    processor = atr['INSAR_PROCESSOR']
    except: processor = atr['PROCESSOR']
    if processor in ['roipac','isce']:
        Re = float(atr['EARTH_RADIUS'])
        Height = float(atr['HEIGHT'])
        az_step = float(atr['AZIMUTH_PIXEL_SIZE']) *Re/(Re+Height)
    elif processor == 'gamma':
        try: atr = readfile.attribute_gamma2roipac(atr)
        except: pass
        az_step = float(atr['AZIMUTH_PIXEL_SIZE'])
    return az_step


#########################################################################
##### Use geomap*.trans file for precious (pixel-level) coord conversion
def get_lookup_row_col(y, x, lut_y, lut_x, y_factor=10, x_factor=10, geoCoord=False):
    '''Get row/col number in y/x value matrix from input y/x
    Use overlap mean value between y and x buffer;
    To support point outside of value pool/matrix, could use np.polyfit to fit a line
    for y and x value buffer and return the intersection point row/col
    '''
    ymin = y - y_factor
    xmin = x - x_factor
    if not geoCoord:
        ymin = max(ymin, 0.5)
        xmin = max(xmin, 0.5)
    mask_y = np.multiply(lut_y >= ymin, lut_y <= (y+y_factor))
    mask_x = np.multiply(lut_x >= xmin, lut_x <= (x+x_factor))
    row, col = np.nanmean(np.where(np.multiply(mask_y, mask_x)), axis=1)
    return row, col

def glob2radar(lat, lon, lookupFile=None, atr_rdr=dict(), print_msg=True):
    '''Convert geo coordinates into radar coordinates.
    Inputs:
        lat/lon    - np.array, float, latitude/longitude
        lookupFile - string, trans/look up file
        atr_rdr    - dict, attributes of file in radar coord, optional but recommended.
    Output:
        az/rg     - np.array, float, range/azimuth pixel number
        az/rg_res - float, residul/uncertainty of coordinate conversion
    '''
    lookupFile = get_lookup_file(lookupFile)
    if not lookupFile:
        print('WARNING: No lookup table found! Can not convert coordinates without it.')
        return None
    atr_lut = readfile.read_attribute(lookupFile)
    if print_msg:
        print 'reading file: '+lookupFile

    #####For lookup table in geo-coord, read value directly
    if 'Y_FIRST' in atr_lut.keys():
        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3
        lut_x = readfile.read(lookupFile, epoch='rangeCoord')[0]
        lut_y = readfile.read(lookupFile, epoch='azimuthCoord')[0]
        lat0 = float(atr_lut['Y_FIRST'])
        lon0 = float(atr_lut['X_FIRST'])
        lat_center = lat0 + float(atr_lut['Y_STEP'])*float(atr_lut['FILE_LENGTH'])/2
        lat_step_deg = float(atr_lut['Y_STEP'])
        lon_step_deg = float(atr_lut['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0*earth_radius
        lon_step = lon_step_deg*np.pi/180.0*earth_radius*np.cos(lat_center*np.pi/180)

        # Get range/azimuth ground resolution/step in meter
        x_factor = 2
        y_factor = 2
        az0 = 0
        rg0 = 0
        if 'Y_FIRST' not in atr_rdr.keys():
            az_step = azimuth_ground_resolution(atr_rdr)
            rg_step = range_ground_resolution(atr_rdr, print_msg)
            x_factor = np.ceil(abs(lon_step)/rg_step).astype(int)
            y_factor = np.ceil(abs(lat_step)/az_step).astype(int)
            if 'subset_y0' in atr_rdr.keys():
                az0 = int(atr_rdr['subset_y0'])
            if 'subset_x0' in atr_rdr.keys():
                rg0 = int(atr_rdr['subset_x0'])

        width  = int(atr_lut['WIDTH'])
        row = np.rint((lat - lat0)/lat_step_deg).astype(int)
        col = np.rint((lon - lon0)/lon_step_deg).astype(int)
        rg = np.rint(lut_x[row, col]).astype(int) - rg0
        az = np.rint(lut_y[row, col]).astype(int) - az0


    #####For lookup table in radar-coord, search the buffer and use center pixel
    else:
        lut_x = readfile.read(lookupFile, epoch='lon')[0]
        lut_y = readfile.read(lookupFile, epoch='lat')[0]
        az = np.zeros(lat.shape)
        rg = np.zeros(lat.shape)
        x_factor = 10
        y_factor = 10
        try:    earth_radius = float(atr_rdr['EARTH_RADIUS'])
        except: earth_radius = 6371.0e3
        az_step = azimuth_ground_resolution(atr_rdr)
        rg_step = range_ground_resolution(atr_rdr)
        lat0 = np.nanmax(lat)
        lat1 = np.nanmin(lat)
        az_step_deg = 180. / np.pi * az_step / earth_radius
        rg_step_deg = 180. / np.pi * rg_step / (earth_radius*np.cos((lat0+lat1)/2*np.pi/180.))

        if lat.size == 1:
            az, rg = get_lookup_row_col(lat, lon, lut_y, lut_x,\
                                        y_factor*az_step_deg, x_factor*rg_step_deg, geoCoord=True)
        else:
            for i in range(rg.size):
                az[i], rg[i] = get_lookup_row_col(lat[i], lon[i], lut_y, lut_x,\
                                                  y_factor*az_step_deg, x_factor*rg_step_deg, geoCoord=True)
        az = np.rint(az).astype(int)
        rg = np.rint(rg).astype(int)
    rg_resid = x_factor
    az_resid = y_factor
    return az, rg, az_resid, rg_resid


def radar2glob(az, rg, lookupFile=None, atr_rdr=dict(), print_msg=True):
    '''Convert radar coordinates into geo coordinates
    Inputs:
        rg/az      - np.array, int, range/azimuth pixel number
        lookupFile - string, trans/look up file
        atr_rdr    - dict, attributes of file in radar coord, optional but recommended.
    Output:
        lon/lat    - np.array, float, longitude/latitude of input point (rg,az); nan if not found.
        latlon_res - float, residul/uncertainty of coordinate conversion
    '''
    lookupFile = get_lookup_file(lookupFile)
    if not lookupFile:
        print('WARNING: No lookup table found! Can not convert coordinates without it.')
        return None
    atr_lut = readfile.read_attribute(lookupFile)
    if print_msg:
        print 'reading file: '+lookupFile

    #####For lookup table in geo-coord, search the buffer and use center pixel
    if 'Y_FIRST' in atr_lut.keys():
        if 'subset_x0' in atr_rdr.keys():
            rg += int(atr_rdr['subset_x0'])
            az += int(atr_rdr['subset_y0'])        

        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3;    # in meter
        lut_x = readfile.read(lookupFile, epoch='rangeCoord')[0]
        lut_y = readfile.read(lookupFile, epoch='azimuthCoord')[0]
        lat0 = float(atr_lut['Y_FIRST'])
        lon0 = float(atr_lut['X_FIRST'])
        lat_center = lat0 + float(atr_lut['Y_STEP'])*float(atr_lut['FILE_LENGTH'])/2
        lat_step_deg = float(atr_lut['Y_STEP'])
        lon_step_deg = float(atr_lut['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0*earth_radius
        lon_step = lon_step_deg*np.pi/180.0*earth_radius*np.cos(lat_center*np.pi/180)

        # Get range/azimuth ground resolution/step
        x_factor = 10
        y_factor = 10
        if 'Y_FIRST' not in atr_rdr.keys():
            az_step = azimuth_ground_resolution(atr_rdr)
            rg_step = range_ground_resolution(atr_rdr, print_msg)
            x_factor = 2*np.ceil(abs(lon_step)/rg_step)
            y_factor = 2*np.ceil(abs(lat_step)/az_step)

        lut_row = np.zeros(rg.shape)
        lut_col = np.zeros(rg.shape)
        if rg.size == 1:
            lut_row, lut_col = get_lookup_row_col(az, rg, lut_y, lut_x, y_factor, x_factor)
        else:
            for i in range(rg.size):
                lut_row[i], lut_col[i] = get_lookup_row_col(az[i], rg[i], lut_y, lut_x, y_factor, x_factor)
        lat = lut_row*lat_step_deg + lat0
        lon = lut_col*lon_step_deg + lon0
        lat_resid = abs(y_factor*lat_step_deg)
        lon_resid = abs(x_factor*lon_step_deg)

    #####For lookup table in radar-coord, read the value directly.
    else:
        lut_x = readfile.read(lookupFile, epoch='lon')[0]
        lut_y = readfile.read(lookupFile, epoch='lat')[0]
        lat = lut_y[az, rg]
        lon = lut_x[az, rg]

        x_factor = 2
        y_factor = 2
        try:    earth_radius = float(atr_rdr['EARTH_RADIUS'])
        except: earth_radius = 6371.0e3
        az_step = azimuth_ground_resolution(atr_rdr)
        rg_step = range_ground_resolution(atr_rdr)
        lat0 = np.nanmax(lat)
        lat1 = np.nanmin(lat)
        az_step_deg = 180. / np.pi * az_step / earth_radius
        rg_step_deg = 180. / np.pi * rg_step / (earth_radius*np.cos((lat0+lat1)/2*np.pi/180.))
        lat_resid = abs(y_factor * az_step_deg)
        lon_resid = abs(x_factor * rg_step_deg)

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
def design_matrix(ifgramFile=None, date12_list=[], zero_first=True):
    '''Make the design matrix for the inversion based on date12_list.
    Reference:
        Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002).
        A new algorithm for surface deformation monitoring based on small
        baseline differential SAR interferograms. IEEE TGRS, 40(11), 2375-2383.

    Input:
        ifgramFile  - string, name/path of interferograms file
        date12_list - list of string, date12 used in calculation in YYMMDD-YYMMDD format
                      use all date12 from ifgramFile if input is empty
    Outputs:
        A - 2D np.array in size of (ifgram_num, date_num-1)
            representing date combination for each interferogram (-1 for master, 1 for slave, 0 for others)
        B - 2D np.array in size of (ifgram_num, date_num-1)
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
    ifgram_num = len(date12_list)

    A = np.zeros((ifgram_num, date_num))
    B = np.zeros(np.shape(A))
    #t = np.zeros((ifgram_num, 2))
    for i in range(ifgram_num):
        m_idx, s_idx = [date6_list.index(j) for j in date12_list[i].split('-')]
        A[i, m_idx] = -1
        A[i, s_idx] = 1
        B[i, m_idx:s_idx] = tbase[m_idx+1:s_idx+1] - tbase[m_idx:s_idx]
        #t[i,:] = [tbase[m_idx], tbase[s_idx]]

    # Remove the 1st date assuming it's zero
    if zero_first:
        A = A[:,1:]
        B = B[:,:-1]

    return A,B


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
    ts_data = np.zeros((date_num,pixel_num),np.float32)
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
            ts_data[:,ni] = defo
        #if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,pixel_num)
        if not np.remainder(ni,pixel_num_step):
            print 'Processing point: %8d of %8d, %3d' % (ni,pixel_num,(10*ni/pixel_num_step))+'%'
    del data
    timeseries = np.zeros((date_num,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
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
    ts_data = np.zeros((date_num,pixel_num),np.float32)
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
    
            ts_data[:,ni] = defo
        if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,pixel_num)
    del data
    timeseries = np.zeros((date_num,np.shape(dset)[0],np.shape(dset)[1]),np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1],np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
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
    if stack is None:
        print 'calculating stack of input file ...'
        stack = stacking(File)

    # set masked out area into NaN
    if maskFile:
        print 'read mask from file: '+maskFile
        mask = readfile.read(maskFile, epoch='mask')[0]
        stack[mask==0] = np.nan

    return stack


def stacking(File):
    '''Stack multi-temporal dataset into one equivalent to temporal sum
    For interferograms, the averaged velocity is calculated.
    '''

    ## File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])
    if k in ['interferograms']:
        phase2range = -1 * float(atr['WAVELENGTH']) / (4.0 * np.pi)
        atr['FILE_TYPE'] = 'velocity'
        atr['UNIT'] = 'm/yr'
    else:
        atr['FILE_TYPE'] = 'mask'

    ## Calculation
    stack = np.zeros([length,width])
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Input File Info
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum  = len(epochList)
        prog_bar = ptime.progress_bar(maxValue=epochNum, prefix='calculating: ')
        for i in range(epochNum):
            epoch = epochList[i]
            if k == 'timeseries':
                data = h5file[k].get(epoch)[:]
            else:
                data = h5file[k][epoch].get(epoch)[:]
                if k in ['interferograms']:
                    m_date, s_date = h5file[k][epoch].attrs['DATE12'].split('-')
                    t1 = datetime.datetime(*time.strptime(m_date, "%y%m%d")[0:5])
                    t2 = datetime.datetime(*time.strptime(s_date, "%y%m%d")[0:5])
                    dt = float((t2-t1).days)/365.25
                    data *= phase2range / dt
            stack += data
            prog_bar.update(i+1)
        stack *= 1.0/float(epochNum)
        prog_bar.close()
        h5file.close()

        # Write stack file is input file is multi-dataset (large file size usually)
        stackFile = os.path.splitext(File)[0]+'_stack.h5'
        print 'writing stack file >>> '+stackFile
        writefile.write(stack, atr, stackFile)

    else:
        try:
            stack, atrStack = readfile.read(File)
        except:
            print 'Cannot read file: '+File; sys.exit(1)
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


def generate_curls(curlfile, h5file, Triangles, curls):
    ifgram_list = h5file['interferograms'].keys()
    h5curlfile = h5py.File(curlfile,'w')
    gg = h5curlfile.create_group('interferograms')

    curl_num = np.shape(curls)[0]
    prog_bar = ptime.progress_bar(maxValue=curl_num)
    for i in range(curl_num):
        ifgram1 = ifgram_list[curls[i,0]]
        ifgram2 = ifgram_list[curls[i,1]]
        ifgram3 = ifgram_list[curls[i,2]]
        d1 = h5file['interferograms'][ifgram1].get(ifgram1)[:]
        d2 = h5file['interferograms'][ifgram2].get(ifgram2)[:]
        d3 = h5file['interferograms'][ifgram3].get(ifgram3)[:]

        triangle_date = Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2]
        group = gg.create_group(triangle_date)
        dset = group.create_dataset(triangle_date, data=d1+d3-d2, compression='gzip')
        for key, value in h5file['interferograms'][ifgram1].attrs.iteritems():
            group.attrs[key] = value
        prog_bar.update(i+1)

    h5curlfile.close()
    prog_bar.close()
    return curlfile

