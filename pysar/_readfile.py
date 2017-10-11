#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
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
multi_dataset_hdf5_file=['timeseries','geometry']
single_dataset_hdf5_file=['dem','mask','rmse','temporal_coherence', 'velocity']


#########################################################################
def read(File, box=None, epoch=None):
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
    fbase = os.path.splitext(os.path.basename(File))[0]
    atr = read_attribute(File, epoch)
    k = atr['FILE_TYPE']
    processor = atr['INSAR_PROCESSOR']
    length = int(float(atr['FILE_LENGTH']))
    width = int(float(atr['WIDTH']))
    if not box:
        box = (0, 0, width, length)

    ##### HDF5
    if ext in ['.h5','.he5']:
        h5file = h5py.File(File,'r')

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

        data = dset[box[1]:box[3],box[0]:box[2]]

        h5file.close()
        return data, atr

    ##### Image
    elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        atr = read_roipac_rsc(File+'.rsc')
        data  = Image.open(File).crop(box)
        return data, atr

    ##### ISCE
    elif processor in ['isce']:
        if k in ['.unw','unw']:
            amp, pha, atr = read_float32(File, box=box)
            return pha, atr

        elif k in ['.cor','cor']:
            data, atr = read_real_float32(File, box=box)
            return data, atr

        elif k in ['.int','int']:
            amp, pha, atr = read_complex_float32(File, box=box, cpx=False)
            return pha, atr

        elif k in ['.slc']:
            amp, pha, atr = read_complex_float32(File, box=box, cpx=False)
            return amp, atr

        elif k in ['.flat','cpx']:
            amp, pha, atr = read_complex_float32(File, box=box, cpx=False)
            return pha, atr

        elif fbase.startswith('los'):
            incAngle, azAngle, atr = read_float32(File, box=box)
            if not epoch:
                return incAngle, azAngle, atr
            elif epoch.startswith('inc'):
                return incAngle, atr
            elif epoch.startswith(('az','head')):
                return azAngle, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

        elif atr['DATA_TYPE'].lower() in ['float64','double']:
            data, atr = read_real_float64(File, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['cfloat32']:
            data, atr = read_complex_float32(File, box=box, cpx=True)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['float32','float']:
            data, atr = read_real_float32(File, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['int16','short']:
            data, atr = read_real_int16(File, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['bool','byte','flag']:
            data, atr = read_bool(File, box=box)
            return data, atr
        else:
            sys.exit('Un-supported '+processor+' file format: '+os.path.basename(File))

    ##### ROI_PAC
    elif processor in ['roipac']:
        if ext in ['.unw','.cor','.hgt', '.msk']:
            amp, pha, atr = read_float32(File, box=box)
            return pha, atr

        elif ext in ['.dem']:
            dem, atr = read_real_int16(File, box=box)
            return dem, atr

        elif ext in ['.int']:
            amp, pha, atr = read_complex_float32(File, box=box, cpx=False)
            return pha, atr

        elif ext in ['.amp']:
            amp, atr = read_complex_float32(File, box=box, cpx=True)
            m_amp = amp.real
            s_amp = amp.imag
            return m_amp, s_amp, atr

        elif ext in ['.flg', '.byt']:
            flag, atr = read_bool(File, box=box)
            return flag, atr

        elif ext in ['.trans']:
            rg, az, atr = read_float32(File, box=box)
            if not epoch:
                return rg, az, atr
            elif epoch in ['rg','range']:
                return rg, atr
            elif epoch in ['az','azimuth']:
                return az, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

    ##### Gamma
    elif processor == 'gamma':
        if ext in ['.unw','.cor','.hgt_sim','.dem']:
            data, atr = read_real_float32(File, box=box, byte_order='ieee-be')
            return data, atr

        elif ext in ['.UTM_TO_RDC', '.utm_to_rdc']:
            data, atr = read_complex_float32(File, box=box, byte_order='ieee-be', cpx=True)
            if not epoch:
                return data.real, data.imag, atr
            elif epoch in ['rg','range']:
                return data.real, atr
            elif epoch in ['az','azimuth']:
                return data.imag, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

        elif ext in ['.mli']:
            data, atr = read_real_float32(File, box=box)
            return data, atr

        elif ext in ['.amp','.ramp']:
            data, atr = read_real_float32(File, box=box, byte_order='ieee-be')
            return data, atr

        elif ext == '.slc':
            amp, pha, atr = read_complex_int16(File, box=box, cpx=False)
            return amp, atr

        else:
            sys.exit('Un-supported '+processor+' for file format: '+ext)
    else:
        sys.exit('Unrecognized file format: '+ext)


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
            if epoch:
                # Check input epoch exists or not
                epoch_list = sorted(h5[k[0]].keys())
                try:    epoch = [i for i in epoch_list if epoch in i][0]
                except: epoch = None

            if not epoch:
                epoch = h5[k[0]].keys()[0]
            attrs = h5[k[0]][epoch].attrs

        elif k[0] in multi_dataset_hdf5_file+single_dataset_hdf5_file:
            attrs  = h5[k[0]].attrs
        elif k[0] in ['HDFEOS']:
            attrs = h5.attrs
        else:
            sys.exit('Unrecognized h5 file key: '+k[0])

        atr = dict()
        for key, value in attrs.iteritems():
            atr[key] = str(value)
        atr['FILE_TYPE'] = str(k[0])
        atr['PROCESSOR'] = 'pysar'

        if k[0] == 'timeseries':
            try:    atr['ref_date']
            except: atr['ref_date'] = sorted(h5[k[0]].keys())[0]

        h5.close()

    else:
        # attribute file list
        try:
            rscFileList = [File+'.rsc', File.split('_snap_connect.byt')[0]+'.unw.rsc']
            rscFile = [i for i in rscFileList if os.path.isfile(i)][0]
        except:
            rscFile = None

        ##### RSC File
        if rscFile:
            atr = read_roipac_rsc(rscFile)
            atr['FILE_TYPE'] = ext
            #if 'FILE_TYPE' not in atr.keys():
            #    atr['FILE_TYPE'] = ext
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'roipac'
            if 'INSAR_PROCESSOR' not in atr.keys():
                atr['INSAR_PROCESSOR'] = 'roipac'

        ##### PAR File
        elif os.path.isfile(File+'.par'):
            atr = read_gamma_par(File+'.par')
            atr['FILE_TYPE'] = ext
            #if 'FILE_TYPE' not in atr.keys():
            #    atr['FILE_TYPE'] = ext
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'gamma'
            if 'INSAR_PROCESSOR' not in atr.keys():
                atr['INSAR_PROCESSOR'] = 'gamma'

        ##### XML File
        elif os.path.isfile(File+'.xml'):
            atr = read_isce_xml(File+'.xml')
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = ext
            atr['PROCESSOR'] = 'isce'
            if 'INSAR_PROCESSOR' not in atr.keys():
                atr['INSAR_PROCESSOR'] = 'isce'

        else:
            sys.exit('Unrecognized file extension: '+ext)

    # Unit - str
    #if 'UNIT' not in atr.keys():
    if atr['FILE_TYPE'] in ['interferograms','wrapped','.unw','.int','.flat','unw']:
        atr['UNIT'] = 'radian'
    elif atr['FILE_TYPE'] in ['timeseries','dem','.dem','.hgt']:
        atr['UNIT'] = 'm'
    elif atr['FILE_TYPE'] in ['velocity']:
        atr['UNIT'] = 'm/yr'
    else:
        atr['UNIT'] = '1'

    atr['FILE_PATH'] = os.path.abspath(File)
    if 'INSAR_PROCESSOR' not in atr.keys():
        if atr['PROCESSOR'] == 'pysar':
            atr['INSAR_PROCESSOR'] = 'roipac'
        else:
            atr['INSAR_PROCESSOR'] = atr['PROCESSOR']
    return atr


#########################################################################
def check_variable_name(path):
    s=path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        try:
            p0 = os.getenv(s[1:])
            path = path.replace(path.split("/")[0], p0)
        except:
            print 'WARNING: Un-recognized environmental variable: '+s
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

    par_dict = attribute_gamma2roipac(par_dict)

    return par_dict


def read_isce_xml(File):
    '''Read ISCE .xml file input a python dictionary structure.'''
    tree = ET.parse(File)
    root = tree.getroot()
    xml_dict={}

    for child in root.findall('property'):
        key = child.attrib['name']
        value  = child.find('value').text
        xml_dict[key] = value

    ## Read lat/lon info for geocoded file
    try:
        comp1 = root.find("./component[@name='coordinate1']")
        xml_dict['X_STEP']  = comp1.find("./property[@name='delta']/value").text
        xml_dict['X_FIRST'] = comp1.find("./property[@name='startingvalue']/value").text
        xml_dict['X_LAST']  = comp1.find("./property[@name='endingvalue']/value").text
    except: pass

    try:
        comp2 = root.find("./component[@name='coordinate2']")
        xml_dict['Y_STEP']  = comp2.find("./property[@name='delta']/value").text
        xml_dict['Y_FIRST'] = comp2.find("./property[@name='startingvalue']/value").text
        xml_dict['Y_LAST']  = comp2.find("./property[@name='endingvalue']/value").text
    except: pass

    xml_dict = attribute_isce2roipac(xml_dict)
    return xml_dict


def attribute_gamma2roipac(par_dict_in):
    '''Convert Gamma par attribute into ROI_PAC format'''
    par_dict = dict()
    for key, value in par_dict_in.iteritems():
        par_dict[key] = value
    key_list = par_dict.keys()

    # Length - number of rows
    key = 'azimuth_lines'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    key = 'interferogram_azimuth_lines'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    key = 'nlines'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    key = 'az_samp_1'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    # Width - number of columns
    key = 'range_samples'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    key = 'interferogram_width'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    key = 'width'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    key = 'range_samp_1'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    # WAVELENGTH
    speed_of_light = 299792458.0   # meter/second
    key = 'radar_frequency'
    if key in key_list:
        par_dict['WAVELENGTH'] = str(speed_of_light/float(par_dict[key]))

    # HEIGHT & EARTH_RADIUS
    key = 'earth_radius_below_sensor'
    if key in key_list:
        par_dict['EARTH_RADIUS'] = par_dict[key]

        key2 = 'sar_to_earth_center'
        if key2 in key_list:
            par_dict['HEIGHT'] = str(float(par_dict[key2]) - float(par_dict[key]))

    # UTC TIME
    key = 'center_time'
    if key in key_list:
        par_dict['CENTER_LINE_UTC'] = par_dict[key]

    # STARTING_RANGE
    key = 'near_range_slc'
    if key in key_list:
        par_dict['STARTING_RANGE'] = par_dict[key]

    # RANGE_PIXEL_SIZE
    key = 'range_pixel_spacing'
    if key in key_list:
        par_dict['RANGE_PIXEL_SIZE'] = par_dict[key]

    key = 'interferogram_range_pixel_spacing'
    if key in key_list:
        par_dict['RANGE_PIXEL_SIZE'] = par_dict[key]

    key = 'range_pixel_spacing_1'
    if key in key_list:
        par_dict['RANGE_PIXEL_SIZE'] = par_dict[key]

    # PLATFORM
    key = 'sensor'
    if key in key_list:
        par_dict['PLATFORM'] = par_dict[key]

    # ORBIT_DIRECTION
    key = 'heading'
    if key in key_list:
        value = float(par_dict[key])
        if 270 < value < 360 or -90 < value < 90:
            par_dict['ORBIT_DIRECTION'] = 'ascending'
        else:
            par_dict['ORBIT_DIRECTION'] = 'descending'

        par_dict['HEADING'] = str(value)


    ##### attributes in geo coordinates
    key = 'corner_lat'
    if key in key_list:
        par_dict['Y_FIRST'] = par_dict[key]

    key = 'corner_lon'
    if key in key_list:
        par_dict['X_FIRST'] = par_dict[key]

    key = 'post_lat'
    if key in key_list:
        par_dict['Y_STEP'] = par_dict[key]

    key = 'post_lon'
    if key in key_list:
        par_dict['X_STEP'] = par_dict[key]


    ##### Optional attributes for PySAR from ROI_PAC
    # ANTENNA_SIDE
    key = 'azimuth_angle'
    if key in key_list:
        value = float(par_dict[key])
        if 0 < value < 180:
            par_dict['ANTENNA_SIDE'] = '-1'
        else:
            par_dict['ANTENNA_SIDE'] = '1'

    # AZIMUTH_PIXEL_SIZE
    key = 'azimuth_pixel_spacing'
    if key in key_list:
        par_dict['AZIMUTH_PIXEL_SIZE'] = par_dict[key]

    key = 'interferogram_azimuth_pixel_spacing'
    if key in key_list:
        par_dict['AZIMUTH_PIXEL_SIZE'] = par_dict[key]

    key = 'az_pixel_spacing_1'
    if key in key_list:
        par_dict['AZIMUTH_PIXEL_SIZE'] = par_dict[key]

    # RLOOKS
    key = 'interferogram_range_looks'
    if key in key_list:
        par_dict['RLOOKS'] = par_dict[key]

    # ALOOKS
    key = 'interferogram_azimuth_looks'
    if key in key_list:
        par_dict['ALOOKS'] = par_dict[key]

    # PRF
    key = 'prf'
    if key in key_list:
        par_dict['PRF'] = par_dict['prf']

    return par_dict


def attribute_isce2roipac(metaDict, dates=[], baselineDict={}):
    '''Convert ISCE xml attribute into ROI_PAC format'''

    rscDict={}
    for key in metaDict.keys():
        rscDict[key] = str(metaDict[key]).strip().split()[0]
    keyList = rscDict.keys()

    rscDict['WIDTH'] = rscDict['width']
    rscDict['FILE_LENGTH'] = rscDict['length']

    if 'dataType' in keyList:
        rscDict['DATA_TYPE'] = rscDict['dataType']
    if 'data_type' in keyList:
        rscDict['DATA_TYPE'] = rscDict['data_type']

    if 'radarWavelength' in keyList:
        rscDict['WAVELENGTH'] = rscDict['radarWavelength']

    rscDict['PROCESSOR'] = 'isce'
    rscDict['INSAR_PROCESSOR'] = 'isce'
    rscDict['PLATFORM'] = 'Sentinel1'

    if 'rangePixelSize' in keyList:
        rscDict['RANGE_PIXEL_SIZE'] = rscDict['rangePixelSize']

    if 'azimuthPixelSize' in keyList:
        rscDict['AZIMUTH_PIXEL_SIZE'] = rscDict['azimuthPixelSize']

    if 'earthRadius' in keyList:
        rscDict['EARTH_RADIUS'] = rscDict['earthRadius']

    if 'altitude' in keyList:
        rscDict['HEIGHT'] = rscDict['altitude']

    if 'startingRange' in keyList:
        rscDict['STARTING_RANGE']  = rscDict['startingRange']
        rscDict['STARTING_RANGE1'] = rscDict['startingRange']

    rscDict['ANTENNA_SIDE'] = '-1'
    if 'passDirection' in keyList:
        rscDict['ORBIT_DIRECTION'] = rscDict['passDirection']

    if dates:
        rscDict['DATE12'] = str(dates[0][2:]+'-'+dates[1][2:])
        #rscDict['DATE'] = str(dates[0])

    if dates and baselineDict:
        bperp = baselineDict['bperp'][dates[1]] - baselineDict['bperp'][dates[0]]
        bpar  = baselineDict['bpar'][dates[1]]  - baselineDict['bpar'][dates[0]]
        rscDict['P_BASELINE_TOP_HDR']    = str(bperp)
        rscDict['P_BASELINE_BOTTOM_HDR'] = str(bperp)
        rscDict['H_BASELINE_TOP_HDR']    = str(bpar)
        rscDict['H_BASELINE_BOTTOM_HDR'] = str(bpar)

    return rscDict


#########################################################################
def read_float32(File, box=None, byte_order='l'):
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

    data_type = 'f4'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f4'

    data = np.fromfile(File, dtype=data_type, count=box[3]*2*width).reshape(box[3],2*width)
    amplitude = data[box[1]:box[3], box[0]:box[2]]
    phase = data[box[1]:box[3], width+box[0]:width+box[2]]

    return amplitude, phase, atr


def read_real_float64(fname, box=None, byte_order='l'):
    '''Read real float64/double data matrix, i.e. isce lat/lon.rdr
    '''
    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f8'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f8'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[1]:box[2]]
    return data, atr

def read_complex_float32(fname, box=None, byte_order='l', cpx=False):
    '''Read complex float 32 data matrix, i.e. roi_pac int or slc data.
    old name: read_complex64()

    ROI_PAC file: .slc, .int, .amp
    
    Data is sotred as:
    real, imaginary, real, imaginary, ...
    real, imaginary, real, imaginary, ...
    ...

    Inputs:
        fname      : str, input file name
        box        : 4-tuple defining (left, upper, right, lower) pixel coordinate.
        byte_order : str, optional, order of reading byte in the file
        cpx        : flag for output format, 
                    0 for amplitude and phase [by default], 
                    non-0 : for real and imagery
    Output:
        data : 2D np.array in complex float32 
    Example:
        amp, phase, atr = read_complex_float32('geo_070603-070721_0048_00018.int')
        data, atr       = read_complex_float32('150707.slc', 1)
    '''

    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'c8'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>c8'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[1]:box[2]]

    if cpx:
        return data, atr
    else:
        amplitude = np.hypot(data.imag, data.real)
        phase = np.arctan2(data.imag, data.real)
        return amplitude, phase, atr


def read_real_float32(fname, box=None, byte_order='l'):
    '''Read real float 32 data matrix, i.e. GAMMA .mli file
    Parameters: fname     : str, path, filename to be read
                byte_order : str, optional, order of reading byte in the file
    Returns: data : 2D np.array, data matrix 
             atr  : dict, attribute dictionary
    Usage: data, atr = read_real_float32('20070603.mli')
           data, atr = read_real_float32('diff_filt_130118-130129_4rlks.unw')
    '''
    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f4'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f4'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[1]:box[2]]
    return data, atr


def read_complex_int16(File, box=None, byte_order='l', cpx=False):
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

    data_type = 'i2'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>i2'

    data = np.fromfile(File, dtype=data_type, count=box[3]*2*width).reshape(box[3],2*width)
    data = data[box[1]:box[3],2*box[0]:2*box[2]].flatten()

    odd_idx = np.arange(1, len(data), 2)
    real = data[odd_idx-1].reshape(box[3]-box[1],box[2]-box[0])
    imag = data[odd_idx].reshape(box[3]-box[1],box[2]-box[0])

    if cpx:
        return real, imag, atr
    else:
        amplitude = np.hypot(imag,real)
        phase = np.arctan2(imag,real)
        return amplitude, phase, atr


def read_real_int16(File, box=None, byte_order='l'):
    atr = read_attribute(File)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    if not box:
        box = [0,0,width,length]

    data_type = 'i2'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>i2'

    data = np.fromfile(File, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[1]:box[2]]
    return data, atr


def read_bool(File, box=None):
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
    if not box:
        box = [0,0,width,length]
    
    data = np.fromfile(File, dtype=bool, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[1]:box[2]]
    return data, atr


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




