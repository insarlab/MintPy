############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#This program is modified from the software originally written by Scott Baker with 
#the following licence:
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
# Yunjun, Sep 2015: Add read_par_file()
#                   Add read_gamma_float() and read_gamma_scomplex()
# Yunjun, Oct 2015: Add box option for read_float32()
# Heresh, Nov 2015: Add ISCE xml reader
# Yunjun, Jan 2016: Add read()
# Yunjun, May 2016: Add read_attribute() and 'PROCESSOR','FILE_TYPE','UNIT' attributes


import os
import sys
import re

import h5py
import numpy as np
import xml.etree.ElementTree as ET
from PIL import Image
import json


#########################################################################
'''Three types of HDF5 files in PySAR
multi_group   : multiple groups with one      dataset and one attribute dict per group (Ngroup-1dset-1atr)
multi_dataset : one      group  with multiple dataset and one attribute dict per group (1group-Ndset-1atr)
single_dataset: one      group  with one      dataset and one attribute dict per gropu (1group-1dset-1atr)

Recommend usage:
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file
'''
multi_group_hdf5_file=['interferograms','coherence','wrapped','snaphu_connect_component']
multi_dataset_hdf5_file=['timeseries']
single_dataset_hdf5_file=['dem','mask','rmse','temporal_coherence', 'velocity']


#########################################################################
def read(File, box=(), epoch=None):
    '''Read one dataset and its attributes from input file.
    
    Read one dataset, i.e. interferogram, coherence, velocity, dem ...
    return 0 if failed.

    Inputs:
        File  : str, path of file to read
                PySAR   file: interferograms, timeseries, velocity, etc.
                ROI_PAC file: .unw .cor .hgt .dem .trans
                Gamma   file: .mli .slc
                Image   file: .jpeg .jpg .png .ras .bmp
        box   : 4-tuple of int, area to read, defined in (x0, y0, x1, y1) in pixel coordinate
        epoch : string, epoch to read, for multi-dataset files
                for .trans file:
                '' - return both dataset
                rg, range   - for geomap_*.trans file
                az, azimuth - for geomap_*.trans file

    Outputs:
        data : 2-D matrix in numpy.array format, return None if failed
        atr  : dictionary, attributes of data, return None if failed

    Examples:
        data, atr = read('velocity.h5')
        data, atr = read('100120-110214.unw', (100,1100, 500, 2500))
        data, atr = read('timeseries.h5', (), '20101120')
        data, atr = read('timeseries.h5', (100,1100, 500, 2500), '20101120')
        az,   atr = read('geomap*.trans', (), 'azimuth')
        rg,az,atr = read('geomap*.trans')
    '''

    # Basic Info
    ext = os.path.splitext(File)[1].lower()
    atr = read_attribute(File, epoch)
    processor = atr['PROCESSOR']

    ## Update attributes if subset
    #if box:
    #    width = float(atr['WIDTH'])
    #    length = float(atr['FILE_LENGTH'])
    #    if (box[2]-box[0])*(box[3]-box[1]) < width*length:
    #        atr = subset_attribute(atr, box)

    ##### HDF5
    if ext in ['.h5','.he5']:
        h5file = h5py.File(File,'r')
        k = atr['FILE_TYPE']

        # Read Dataset
        if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
            # Check input epoch exists or not
            epoch_list = sorted(h5file[k].keys())
            try:    epoch2read = [i for i in epoch_list if epoch in i][0]
            except: epoch2read = None
            if not epoch2read:
                print 'ERROR: no input epoch found!'
                print 'input epoch: '+str(epoch)
                print 'available epoches: '+str(epoch_list)
                sys.exit(1)

            elif k in multi_dataset_hdf5_file:
                dset = h5file[k].get(epoch2read)
            else:
                dset = h5file[k][epoch2read].get(epoch2read)

        elif k in single_dataset_hdf5_file:
            dset = h5file[k].get(k)
        else:
            print 'ERROR: Unrecognized h5 file type: '+k
            sys.exit(1)

        # Crop
        if box:
            data = dset[box[1]:box[3],box[0]:box[2]]
        else:
            data = dset[:,:]

        h5file.close()
        return data, atr

    ##### Image
    elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        atr = read_roipac_rsc(File+'.rsc')
        data  = Image.open(File)
        if box:  data = data.crop(box)
        return data, atr

    ##### ISCE
    elif processor == 'isce':
        if   ext in ['.flat']:
            amp, data, atr = read_complex_float32(File)
        elif ext in ['.cor']:
            data, atr = read_real_float32(File)
        elif ext in ['.slc']:
            data, pha, atr = read_complex_float32(File)
            #ind = np.nonzero(data)
            #data[ind] = np.log10(data[ind])     # dB
            #atr['UNIT'] = 'dB'
        else:
            print 'Un-supported '+processfor+' file format: '+ext
            sys.exit(1)
  
        if box:  data = data[box[1]:box[3],box[0]:box[2]]
        return data, atr

    ##### ROI_PAC
    elif processor in ['roipac']:
        if ext in ['.unw','.cor','.hgt', '.msk']:
            if box:
                amp,pha,atr = read_float32(File,box)
            else:
                amp,pha,atr = read_float32(File)
            return pha, atr

        elif ext in ['.dem']:
            dem,atr = read_real_int16(File)
            if box:  dem = dem[box[1]:box[3],box[0]:box[2]]
            return dem, atr
  
        elif ext in ['.int']:
            amp, pha, atr = read_complex_float32(File)
            if box:
                pha = pha[box[1]:box[3],box[0]:box[2]]
            return pha, atr
        elif ext in ['.amp']:
            masterAmplitude, slaveAmplitude, atr = read_complex_float32(File, real_imag=True)
            if box:
                masterAmplitude = masterAmplitude[box[1]:box[3],box[0]:box[2]]
                slaveAmplitude = slaveAmplitude[box[1]:box[3],box[0]:box[2]]
            return masterAmplitude, slaveAmplitude, atr
        elif ext in ['.flg', '.byt']:
            flag, atr = read_flag(File)
            return flag, atr
  
        elif ext == '.trans':
            if box:
                rg,az,atr = read_float32(File,box)
            else:
                rg,az,atr = read_float32(File)

            if not epoch:
                #print 'read range and azimuth from '+File
                return rg, az, atr
            elif epoch in ['rg','range']:
                #print 'read range from '+File
                return rg, atr
            elif epoch in ['az','azimuth']:
                #print 'read azimuth from '+File
                return az, atr
            else:
                print 'Un-recognized epoch input: '+epoch
                sys.exit(1)

    ##### Gamma
    elif processor == 'gamma':
        if ext in ['.unw','.cor','.hgt_sim']:
            data, atr = read_real_float32(File, byteorder='ieee-be')
            if box: data = data[box[1]:box[3],box[0]:box[2]]
            return data, atr

        elif ext == '.mli':
            data,atr = read_real_float32(File)
            if box: data = data[box[1]:box[3],box[0]:box[2]]
            return data, atr

        elif ext == '.slc':
            if box:
                amplitude, phase,atr = read_complex_int16(File, box)
            else:
                amplitude, phase, atr = read_complex_int16(File)
            del phase
            return amplitude, atr

        else:
            print 'Un-supported '+process+' for file format: '+ext
            sys.exit(1)
    else:
        print 'Unrecognized file format: '+ext
        sys.exit(1)


#########################################################################
def read_attribute(File, epoch=None):
    '''Read attributes of input file into a dictionary
    Input  : string, file name and epoch (optional)
    Output : dictionary, attributes dictionary
    '''
    ext = os.path.splitext(File)[1].lower()
    if not os.path.isfile(File):
        print 'Input file not existed: '+File
        print 'Current directory: '+os.getcwd()
        sys.exit(1)

    ##### PySAR
    if ext in ['.h5','.he5']:
        h5 = h5py.File(File,'r')
        k = h5.keys()
        if   'interferograms' in k: k[0] = 'interferograms'
        elif 'coherence'      in k: k[0] = 'coherence'
        elif 'timeseries'     in k: k[0] = 'timeseries'

        if k[0] in multi_group_hdf5_file:
            # Check input epoch exists or not
            epoch_list = sorted(h5[k[0]].keys())
            try:    epoch2read = [i for i in epoch_list if epoch in i][0]
            except: epoch2read = None
            if not epoch2read:
                print 'ERROR: no input epoch found!'
                print 'input epoch: '+str(epoch)
                print 'available epoches: '+str(epoch_list)
                sys.exit(1)

            if epoch2read:
                attrs  = h5[k[0]][epoch2read].attrs
            else:
                attrs  = h5[k[0]][h5[k[0]].keys()[0]].attrs

        elif k[0] in multi_dataset_hdf5_file+single_dataset_hdf5_file:
            attrs  = h5[k[0]].attrs
        else:
            sys.exit('Unrecognized h5 file key: '+k[0])

        atr = dict()
        for key, value in attrs.iteritems():  atr[key] = str(value)
        atr['PROCESSOR'] = 'pysar'
        atr['FILE_TYPE'] = str(k[0])

        if k[0] == 'timeseries':
            try:    atr['ref_date']
            except: atr['ref_date'] = sorted(h5[k[0]].keys())[0]

        h5.close()

    else:
        # attribute file list
        try:
            potentialRscFileList = [File+'.rsc', File.split('_snap_connect.byt')[0]+'.unw.rsc']
            rscFile = [rscFile for rscFile in potentialRscFileList if os.path.isfile(rscFile)][0]
        except: pass

        ##### GAMMA
        if os.path.isfile(File+'.par'):
            atr = read_gamma_par(File+'.par')
            atr['FILE_TYPE'] = ext
            #if 'FILE_TYPE' not in atr.keys():
            #    atr['FILE_TYPE'] = ext
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'gamma'

        ##### ROI_PAC
        elif rscFile:
            atr = read_roipac_rsc(rscFile)
            atr['FILE_TYPE'] = ext
            #if 'FILE_TYPE' not in atr.keys():
            #    atr['FILE_TYPE'] = ext
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'roipac'

        ##### ISCE
        elif os.path.isfile(File+'.xml'):
            atr = read_isce_xml(File+'.xml')
            atr['PROCESSOR'] = 'isce'
            atr['FILE_TYPE'] = ext
    
        else: print 'Unrecognized file extension: '+ext; sys.exit(1)

    # Unit - str
    #if 'UNIT' not in atr.keys():
    if atr['FILE_TYPE'] in ['interferograms','wrapped','.unw','.int','.flat']:
        atr['UNIT'] = 'radian'
    elif atr['FILE_TYPE'] in ['timeseries','dem','.dem','.hgt']:
        atr['UNIT'] = 'm'
    elif atr['FILE_TYPE'] in ['velocity']:
        atr['UNIT'] = 'm/yr'
    else:
        atr['UNIT'] = '1'

    atr['FILE_PATH'] = os.path.abspath(File)

    return atr


#########################################################################
def check_variable_name(path):
    s=path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        p0=os.getenv(s[1:])
        path=path.replace(path.split("/")[0],p0)
    return path

def is_plot_attribute(attribute):
    tokens = attribute.split(".")
    if tokens is None:
        return False

    return tokens[0] == "plot" and len(tokens) > 1

def read_template(File, delimiter='='):
    '''Reads the template file into a python dictionary structure.
    Input : string, full path to the template file
    Output: dictionary, pysar template content
    Example:
        tmpl = read_template(KyushuT424F610_640AlosA.template)
        tmpl = read_template(R1_54014_ST5_L0_F898.000.pi, ':')
    '''
    template_dict = {}
    plotAttributeDict = {}
    insidePlotObject = False
    plotAttributes = []
    # the below logic for plotattributes object can be made much more simple
    # if we assume that any plot attribute coming after a > belongs to the
    # same object. Must Ask Falk and Yunjung if we can assume this to eliminate
    # all these conditionals
    for line in open(File):
        line = line.strip()
        c = [i.strip() for i in line.split(delimiter, 1)]  #split on the 1st occurrence of delimiter
        if len(c) < 2 or line.startswith(('%','#')):
            if line.startswith(">"):
                plotAttributeDict = {}
                insidePlotObject = True
            # otherwise, if previously inside attributes object, we are now outsid    e
            # unless the line is a comment
            elif insidePlotObject and not line.startswith('%') and not line.startswith('#'):
                # just came from being inside plot object, but now we are outside
                insidePlotObject = False
                plotAttributes.append(plotAttributeDict)
            next #ignore commented lines or those without variables
        else:
            atrName  = c[0]
            atrValue = str.replace(c[1],'\n','').split("#")[0].strip()
            atrValue = check_variable_name(atrValue)

            if insidePlotObject:
                if is_plot_attribute(atrName):
                    plotAttributeDict[atrName] = atrValue
                else:
                    # just came from being inside plot object, but now we are outside
                    insidePlotObject = False
                    plotAttributes.append(plotAttributeDict)
                    template_dict[atrName] = atrValue

            elif atrValue != '':
                template_dict[atrName] = atrValue

    # what if no \n at end of file? write out last plot attributes dict
    if insidePlotObject:
        plotAttributes.append(plotAttributeDict)

    if len(plotAttributes) > 0:
        template_dict["plotAttributes"] = json.dumps(plotAttributes)

    return template_dict


def read_roipac_rsc(File):
    '''Read ROI_PAC .rsc file into a python dictionary structure.'''
    rsc_dict = dict(np.loadtxt(File, dtype=str, usecols=(0,1)))
    return rsc_dict


def load_roipac_attribute(fname):
    '''Read/extract attributes for PySAR from ROI_PAC product
    Parameters: fname : str
                    ROIPAC interferogram filename or path, i.e. /PopoSLT143TsxD/filt_130118-130129_4rlks.unw
    Returns:    atr : dict
                    Attributes dictionary
    '''
    atr = {}

    ## Get info: dir, date12, num of loooks
    file_dir = os.path.dirname(fname)
    file_basename = os.path.basename(fname)
    date12 = str(re.findall('\d{6}[-_]\d{6}', file_basename)[0]).replace('_','-')
    m_date, s_date = date12.split('-')

    data_rsc_file = fname+'.rsc'
    baseline_rsc_file = file_dir+'/'+m_date+'_'+s_date+'_baseline.rsc'

    data_rsc_dict     = read_roipac_rsc(data_rsc_file)
    baseline_rsc_dict = read_roipac_rsc(baseline_rsc_file)

    atr.update(data_rsc_dict)
    atr.update(baseline_rsc_dict)

    return atr


def read_gamma_par(fname, delimiter=':', skiprows=3, convert2roipac=True):
    '''Read GAMMA .par/.off file into a python dictionary structure.
    Parameters: fname : file, str, or path. 
                    File path of .par, .off file.
                delimiter : str, optional
                    String used to separate values.
                skiprows : int, optional
                    Skip the first skiprows lines.
    Returns:    par_dict : dict
                    Attributes dictionary
    '''
    par_dict = {}

    # Read txt file
    f = open(fname,'r')
    lines = f.readlines()[skiprows:]
    for line in lines:
        line = line.strip()
        c = [i.strip() for i in line.split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%','#')):
            next
        else:
            key = c[0]
            value = str.replace(c[1],'\n','').split("#")[0].split()[0].strip()
            par_dict[key] = value
    f.close()

    return par_dict


def read_isce_xml(File):
    '''Read ISCE .xml file input a python dictionary structure.'''
    #from lxml import etree as ET
    tree = ET.parse(File)
    root = tree.getroot()
    xmldict={}
    for child in root.findall('property'):
        attrib = child.attrib['name']
        value  = child.find('value').text
        xmldict[attrib]=value

    xmldict['WIDTH']       = xmldict['width']
    xmldict['FILE_LENGTH'] = xmldict['length']
    #xmldict['WAVELENGTH'] = '0.05546576'
    #Date1=os.path.dirname(File).split('/')[-1].split('_')[0][2:]
    #Date2=os.path.dirname(File).split('/')[-1].split('_')[1][2:]
    #xmldict['DATE12'] = Date1 + '-' + Date2
    #xmldict['DATE1'] = Date1
    #xmldict['DATE2'] = Date2
    return xmldict


def merge_attribute(atr1,atr2):
    atr = dict()
    for key, value in atr1.iteritems():  atr[key] = str(value)
    for key, value in atr2.iteritems():  atr[key] = str(value)

    return atr


#########################################################################
def read_float32(File, box=None):
    '''Reads roi_pac data (RMG format, interleaved line by line)
    should rename it to read_rmg_float32()
    
    ROI_PAC file: .unw, .cor, .hgt, .trans, .msk
    
    RMG format (named after JPL radar pionner Richard M. Goldstein): made
    up of real*4 numbers in two arrays side-by-side. The two arrays often
    show the magnitude of the radar image and the phase, although not always
    (sometimes the phase is the correlation). The length and width of each 
    array are given as lines in the metadata (.rsc) file. Thus the total
    width width of the binary file is (2*width) and length is (length), data
    are stored as:
    magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    ......
    
       box  : 4-tuple defining the left, upper, right, and lower pixel coordinate.
    Example:
       a,p,r = read_float32('100102-100403.unw')
       a,p,r = read_float32('100102-100403.unw',(100,1200,500,1500))
    '''

    atr = read_attribute(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0,0,width,length]

    data = np.fromfile(File,np.float32,box[3]*2*width).reshape(box[3],2*width)
    amplitude = data[box[1]:box[3],box[0]:box[2]]
    phase     = data[box[1]:box[3],width+box[0]:width+box[2]]

    #oddindices = np.where(np.arange(length*2)&1)[0]
    #data = np.fromfile(File,np.float32,length*2*width).reshape(length*2,width)
    #amplitude = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
    #phase     = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

    return amplitude, phase, atr


def read_complex_float32(File, real_imag=False):
    '''Read complex float 32 data matrix, i.e. roi_pac int or slc data.
    old name: read_complex64()
    
    ROI_PAC file: .slc, .int, .amp
    
    Data is sotred as:
    real, imaginary, real, imaginary, ...
    real, imaginary, real, imaginary, ...
    ...
    
    Usage:
        File : input file name
        real_imag : flag for output format, 
                    0 for amplitude and phase [by default], 
                    non-0 : for real and imagery
    
    Example:
        amp, phase, atr = read_complex_float32('geo_070603-070721_0048_00018.int')
        data, atr       = read_complex_float32('150707.slc', 1)
    '''

    atr = read_attribute(File)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    ##data = np.fromfile(File,np.complex64,length*2*width).reshape(length,width)
    data = np.fromfile(File,np.complex64,length*width).reshape(length,width)

    if not real_imag:
        amplitude = np.array([np.hypot(  data.real,data.imag)]).reshape(length,width)
        phase     = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
        return amplitude, phase, atr
    else:
        return data, atr


def read_real_float32(fname, byteorder=None):
    '''Read real float 32 data matrix, i.e. GAMMA .mli file
    Parameters: fname     : str, path, filename to be read
                byteorder : str, optional, order of reading byte in the file
    Returns: data : 2D np.array, data matrix 
             atr  : dict, attribute dictionary
    Usage: data, atr = read_real_float32('20070603.mli')
           data, atr = read_real_float32('diff_filt_130118-130129_4rlks.unw')
    '''
    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))

    if byteorder in ['big-endian','b','ieee-be']:
        data = np.fromfile(fname, dtype='>f4').reshape(length, width)
    else:
        data = np.fromfile(fname, dtype=np.float32).reshape(length, width)
    return data, atr


def read_complex_int16(File, box=None, real_imag=False):
    '''Read complex int 16 data matrix, i.e. GAMMA SCOMPLEX file (.slc)
    
    Gamma file: .slc
    
    Inputs:
       file: complex data matrix (cpx_int16)
       box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
    Example:
       data,rsc = read_complex_int16('100102.slc')
       data,rsc = read_complex_int16('100102.slc',(100,1200,500,1500))
    '''

    atr = read_attribute(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0,0,width,length]

    data = np.fromfile(File,np.int16,box[3]*2*width).reshape(box[3],2*width)
    data = data[box[1]:box[3],2*box[0]:2*box[2]].flatten()
    odd_idx = np.arange(1, len(data), 2)
    real = data[odd_idx-1].reshape(box[3]-box[1],box[2]-box[0])
    imag = data[odd_idx].reshape(box[3]-box[1],box[2]-box[0])

    if real_imag:
        return real, imag, atr
    else:
        amplitude = np.array([np.hypot(imag,real)]).reshape(length,width)
        phase = np.array([np.arctan2(imag,real)]).reshape(length,width)
        return amplitude, phase, atr

    #data = np.fromfile(File,np.int16,length*2*width).reshape(length*2,width)
    #oddindices = np.where(np.arange(length*2)&1)[0]
    #real = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
    #imag = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

    #amplitude = np.array([np.hypot(  real,imag)]).reshape(length,width)
    #phase     = np.array([np.arctan2(imag,real)]).reshape(length,width)
    #return amplitude, phase, parContents


def read_dem(File):
    '''Read real int 16 data matrix, i.e. ROI_PAC .dem file.
    Input:  roi_pac format dem file
    Usage:  dem, atr = read_real_int16('gsi10m_30m.dem')
    '''
    atr = read_attribute(File)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH'])) 
    dem = np.fromfile(File, dtype=np.int16).reshape(length, width)
    return dem, atr


def read_real_int16(File):
    '''Same as read_dem() above'''
    atr = read_attribute(File)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH'])) 
    dem = np.fromfile(File, dtype=np.int16).reshape(length, width)
    return dem, atr


def read_flag(File):
    '''Read binary file with flags, 1-byte values with flags set in bits
    For ROI_PAC .flg, *_snap_connect.byt file.    
    '''
    # Read attribute
    if File.endswith('_snap_connect.byt'):
        rscFile = File.split('_snap_connect.byt')[0]+'.unw.rsc'
    else:
        rscFile = File+'.rsc'
    atr = read_attribute(rscFile.split('.rsc')[0])
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH'])) 
    
    flag = np.fromfile(File, dtype=bool).reshape(length, width)
    
    return flag, atr


def read_GPS_USGS(File):  
    yyyymmdd= np.loadtxt(File,dtype=str,usecols = (0,1))[:,0]
    YYYYMMDD=[]
    for y in yyyymmdd:
        YYYYMMDD.append(y)
    data=np.loadtxt(File,usecols = (1,2,3,4))
    dates=data[:,0]
    north=np.array(data[:,1])
    east=np.array(data[:,2])
    up=np.array(data[:,3])
 
    return east,north,up,dates,YYYYMMDD


#########################################################################
def read_multiple(File,box=''):  # Not ready yet
    '''Read multi-temporal 2D datasets into a 3-D data stack
    Inputs:
        File  : input file, interferograms,coherence, timeseries, ...
        box   : 4-tuple defining the left, upper, right, and lower pixel coordinate [optional]
    Examples:
        stack = stacking('timeseries.h5',(100,1200,500,1500))
    '''

    ##### File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    length = int(float(atr['FILE_LENGTH']))
    width  = int(float(atr['WIDTH']))

    ##### Bounding Box
    if box == '':  box = [0,0,width,length]

    epochList = h5file[k].keys()
    epochNum  = len(epochList)
    if epochNum == 0:   print "There is no data in the file";  sys.exit(1)
 
    data = np.zeros([length,width])
    for igram in igramList:
        print igram
        
        dset = h5file[k][igram].get(igram)
        ##### Crop
        try:    data = dset[box[1]:box[3],box[0]:box[2]]
        except: data = dset[:,:]
        unw=dset[0:dset.shape[0],0:dset.shape[1]]
        stack=stack+unw
    return stack




