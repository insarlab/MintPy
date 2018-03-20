###############################################################################
# Program is part of PySAR v2.0 
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun
# Author:  Heresh Fattahi, Zhang Yunjun
###############################################################################
#  This program is modified from the software originally written by Scott Baker
#  with the following licence:
#
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
#
# Recommend usage:
#   import pysar.utils.readfile as readfile
#


import os
import sys
import re
from datetime import datetime as dt

import h5py
import numpy as np
from lxml import objectify
#from PIL import Image
import json


standardMetadatKeys={'width':'WIDTH','Width':'WIDTH','samples':'WIDTH',
                     'length':'LENGTH','FILE_LENGTH':'LENGTH','lines':'WIDTH',
                     'wavelength':'WAVELENGTH','Wavelength':'WAVELENGTH','radarWavelength':'WAVELENGTH',
                     'prf':'PRF',
                     'post_lat':'Y_STEP',
                     'post_lon':'X_STEP',
                     'range_looks':'RLOOKS',
                     'azimuth_looks':'ALOOKS',
                     'dataType':'DATA_TYPE','data_type':'DATA_TYPE',
                     'rangePixelSize':'RANGE_PIXEL_SIZE',
                     'range_pixel_spacing':'RANGE_PIXEL_SIZE','rg_pixel_spacing':'RANGE_PIXEL_SIZE',
                     'azimuthPixelSize':'AZIMUTH_PIXEL_SIZE',
                     'azimuth_pixel_spacing':'AZIMUTH_PIXEL_SIZE','az_pixel_spacing':'AZIMUTH_PIXEL_SIZE',
                     'earthRadius':'EARTH_RADIUS','earth_radius_below_sensor':'EARTH_RADIUS',
                     'altitude':'HEIGHT',
                     'startingRange':'STARTING_RANGE',
                     'center_time':'CENTER_LINE_UTC'
                     }

GDAL2NUMPY_DATATYPE = {

1 : np.uint8,
2 : np.uint16,
3 : np.int16,
4 : np.uint32,
5 : np.int32,
6 : np.float32,
7 : np.float64,
10: np.complex64,
11: np.complex128,

}


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
geometry_dataset=['rangeCoord','azimuthCoord','latitude','longitude','height',\
                  'incidenceAngle','headingAngle','slantRangeDistance','waterMask','shadowMask']


#########################################################################
def read(fname, box=None, epoch=None, print_msg=True):
    '''Read one dataset and its attributes from input file.

    Read one dataset, i.e. interferogram, coherence, velocity, dem ...
    return 0 if failed.

    Inputs:
        fname : str, path of file to read
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
        data, atr = read('100120-110214.unw', box=(100,1100, 500, 2500))
        data, atr = read('timeseries.h5', epoch='20101120')
        data, atr = read('timeseries.h5', box=(100,1100, 500, 2500), epoch='20101120')
        az,   atr = read('geomap*.trans', epoch='azimuth')
        rg,az,atr = read('geomap*.trans')
    '''

    # Basic Info
    ext = os.path.splitext(fname)[1].lower()
    fbase = os.path.splitext(os.path.basename(fname))[0]
    atr = read_attribute(fname, epoch)
    k = atr['FILE_TYPE']
    processor = atr['INSAR_PROCESSOR']
    length = int(float(atr['LENGTH']))
    width = int(float(atr['WIDTH']))
    if not box:
        box = (0, 0, width, length)

    ##### HDF5
    if ext in ['.h5','.he5']:
        h5file = h5py.File(fname,'r')

        # Read Dataset
        if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
            # Check input epoch exists or not
            epoch_list = sorted(h5file[k].keys())
            try:    epoch2read = [i for i in epoch_list if epoch.lower() in i.lower()][0]
            except: epoch2read = None
            if not epoch2read:
                if print_msg:
                    print('ERROR: no input epoch found!')
                    print('input epoch: '+str(epoch))
                    print('available epoches: '+str(epoch_list))
                sys.exit(1)

            elif k in multi_dataset_hdf5_file:
                dset = h5file[k].get(epoch2read)
            else:
                dset = h5file[k][epoch2read].get(epoch2read)

        elif k in ['GIANT_TS']:
            dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in h5file['dates'][:].tolist()]
            dateIndx = dateList.index(epoch)
            if 'rawts' in list(h5file.keys()):
                dset = h5file['rawts'][dateIndx,:,:]
            elif 'recons' in list(h5file.keys()):
                dset = h5file['recons'][dateIndx,:,:]
        else:
            dset = h5file[k].get(k)
        #else:
        #    print 'ERROR: Unrecognized h5 file type: '+k
        #    sys.exit(1)

        data = dset[box[1]:box[3],box[0]:box[2]]

        h5file.close()
        return data, atr

    ###### Image
    #elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
    #    atr = read_roipac_rsc(fname+'.rsc')
    #    data  = Image.open(fname).crop(box)
    #    return data, atr

    ##### ISCE
    elif processor in ['isce']:
        if k in ['.unw','unw']:
            amp, pha, atr = read_float32(fname, box=box)
            return pha, atr

        elif k in ['.cor','cor']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif k in ['.int','int']:
            amp, pha, atr = read_complex_float32(fname, box=box, cpx=False)
            return pha, atr

        elif k in ['.slc']:
            amp, pha, atr = read_complex_float32(fname, box=box, cpx=False)
            return amp, atr

        elif k in ['.flat','cpx']:
            amp, pha, atr = read_complex_float32(fname, box=box, cpx=False)
            return pha, atr

        elif fbase.startswith('los'):
            incAngle, azAngle, atr = read_float32(fname, box=box)
            if not epoch:
                return incAngle, azAngle, atr
            elif epoch.startswith('inc'):
                return incAngle, atr
            elif epoch.startswith(('az','head')):
                return azAngle, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

        elif atr['DATA_TYPE'].lower() in ['float64','double']:
            data, atr = read_real_float64(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['cfloat32']:
            data, atr = read_complex_float32(fname, box=box, cpx=True)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['float32','float']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['int16','short']:
            data, atr = read_real_int16(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['bool','byte','flag']:
            data, atr = read_bool(fname, box=box)
            return data, atr
        else:
            sys.exit('Un-supported '+processor+' file format: '+os.path.basename(fname))

    ##### ROI_PAC
    elif processor in ['roipac']:
        if ext in ['.unw','.cor','.hgt', '.msk']:
            amp, pha, atr = read_float32(fname, box=box)
            return pha, atr

        elif ext in ['.dem','.wgs84']:
            dem, atr = read_real_int16(fname, box=box)
            return dem, atr

        elif ext in ['.int']:
            amp, pha, atr = read_complex_float32(fname, box=box, cpx=False)
            return pha, atr

        elif ext in ['.amp']:
            amp, atr = read_complex_float32(fname, box=box, cpx=True)
            m_amp = amp.real
            s_amp = amp.imag
            return m_amp, s_amp, atr

        elif ext in ['.flt']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif ext in ['.flg', '.byt']:
            flag, atr = read_bool(fname, box=box)
            return flag, atr

        elif ext in ['.trans']:
            rg, az, atr = read_float32(fname, box=box)
            if not epoch:
                return rg, az, atr
            elif epoch.startswith(('rg','range')):
                return rg, atr
            elif epoch.startswith(('az','azimuth')):
                return az, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

    ##### Gamma
    elif processor == 'gamma':
        if ext in ['.unw','.cor','.hgt_sim','.dem']:
            data, atr = read_real_float32(fname, box=box, byte_order='ieee-be')
            return data, atr

        elif ext in ['.UTM_TO_RDC', '.utm_to_rdc']:
            data, atr = read_complex_float32(fname, box=box, byte_order='ieee-be', cpx=True)
            if not epoch:
                return data.real, data.imag, atr
            elif epoch.startswith(('rg','range')):
                return data.real, atr
            elif epoch.startswith(('az','azimuth')):
                return data.imag, atr
            else:
                sys.exit('Un-recognized epoch input: '+epoch)

        elif ext in ['.int']:
            amp, pha, atr = read_complex_float32(fname, box=box, byte_order='ieee-be', cpx=False)
            return pha, atr

        elif ext in ['.mli']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif ext in ['.amp','.ramp']:
            data, atr = read_real_float32(fname, box=box, byte_order='ieee-be')
            return data, atr

        elif ext == '.slc':
            amp, pha, atr = read_complex_int16(fname, box=box, cpx=False)
            return amp, atr

        else:
            sys.exit('Un-supported '+processor+' for file format: '+ext)
    else:
        sys.exit('Unrecognized file format: '+ext)


#########################################################################
def read_attribute(fname, epoch=None):
    '''Read attributes of input file into a dictionary
    Input  : string, file name and epoch (optional)
    Output : dictionary, attributes dictionary
    '''
    ext = os.path.splitext(fname)[1].lower()
    if not os.path.isfile(fname):
        print('Input file not existed: '+fname)
        print('Current directory: '+os.getcwd())
        sys.exit(1)

    ##### PySAR
    if ext in ['.h5','.he5']:
        h5 = h5py.File(fname,'r')
        if   'interferograms' in list(h5.keys()): k = 'interferograms'
        elif 'coherence'      in list(h5.keys()): k = 'coherence'
        elif 'timeseries'     in list(h5.keys()): k = 'timeseries'
        else: k = list(h5.keys())[0]

        attrs = None
        if k in multi_group_hdf5_file:
            if epoch:
                # Check input epoch exists or not
                epoch_list = sorted(h5[k].keys())
                try:    epoch = [i for i in epoch_list if epoch in i][0]
                except: epoch = None

            if not epoch:
                epoch = list(h5[k].keys())[0]
            attrs = h5[k][epoch].attrs

        #elif k in multi_dataset_hdf5_file+single_dataset_hdf5_file:
        else:
            key = 'WIDTH'
            if key in h5.attrs.keys():
                attrs = h5.attrs
            else:
                for groupK in h5.keys():
                    if key in h5[groupK].attrs.keys():
                        attrs = h5[groupK].attrs
                        break
            if fname.endswith('PARAMS.h5'):
                #dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in h5['dates'][:].tolist()]
                attrs = dict()
                attrs['LENGTH'] = h5['cmask'].shape[0]
                attrs['WIDTH'] = h5['cmask'].shape[1]
                #attrs['ORBIT_DIRECTION'] = 'descending'
                #attrs['ref_y'] = '134'
                #attrs['ref_x'] = '637'
                #attrs['REF_DATE'] = '20141225'
                k = 'GIANT_TS'
            if attrs is None:
                raise ValueError('No attribute '+key+' found in 1/2 group level!')

        atr = dict()
        for key, value in attrs.items():
            try:     atr[key] = value.decode('utf-8')
            except:  atr[key] = value
        atr['FILE_TYPE'] = str(k)
        #atr['PROCESSOR'] = 'pysar'

        if k == 'timeseries':
            try:    atr['REF_DATE']
            except: atr['REF_DATE'] = sorted(h5[k].keys())[0]

        h5.close()

    else:
        # attribute file list
        try:
            rscFileList = [fname+'.rsc', fname.split('_snap_connect.byt')[0]+'.unw.rsc']
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
        elif os.path.isfile(fname+'.par'):
            atr = read_gamma_par(fname+'.par')
            atr['FILE_TYPE'] = ext
            #if 'FILE_TYPE' not in atr.keys():
            #    atr['FILE_TYPE'] = ext
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'gamma'
            if 'INSAR_PROCESSOR' not in atr.keys():
                atr['INSAR_PROCESSOR'] = 'gamma'

        ##### XML File
        elif os.path.isfile(fname+'.xml'):
            atr = read_isce_xml(fname+'.xml')
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = ext
            atr['PROCESSOR'] = 'isce'
            if 'INSAR_PROCESSOR' not in atr.keys():
                atr['INSAR_PROCESSOR'] = 'isce'

        elif os.path.isfile(fname+'.hdr'):
            atr = read_template(fname+'.hdr')
            atr = attribute_envi2roipac(atr)
            atr['FILE_TYPE'] = atr['file type']
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
    elif atr['FILE_TYPE'] in ['GIANT_TS']:
        atr['UNIT'] = 'mm'
    else:
        if 'UNIT' not in atr.keys():
            atr['UNIT'] = '1'

    atr['FILE_PATH'] = os.path.abspath(fname)
    if 'INSAR_PROCESSOR' not in atr.keys():
        if atr['PROCESSOR'] == 'pysar':
            atr['INSAR_PROCESSOR'] = 'roipac'
        else:
            atr['INSAR_PROCESSOR'] = atr['PROCESSOR']

    if atr['PROCESSOR'] == 'isce' and ext == '.wgs84':
        atr['FILE_TYPE'] = 'dem'

    atr = standardize_metadata(atr, standardMetadataKeys)
    return atr


#########################################################################
def check_variable_name(path):
    s=path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        try:
            p0 = os.getenv(s[1:])
            path = path.replace(path.split("/")[0], p0)
        except:
            print('WARNING: Un-recognized environmental variable: '+s)
    return path


def is_plot_attribute(attribute):
    tokens = attribute.split(".")
    if tokens is None:
        return False
    return tokens[0] == "plot" and len(tokens) > 1


def read_template(fname, delimiter='='):
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
    for line in open(fname):
        line = line.strip()
        c = [i.strip() for i in line.split(delimiter, 1)]  #split on the 1st occurrence of delimiter
        if len(c) < 2 or line.startswith(('%','#')):
            if line.startswith(">"):
                plotAttributeDict = {}
                insidePlotObject = True
            # otherwise, if previously inside attributes object, we are now outside
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


def read_roipac_rsc(fname):
    '''Read ROI_PAC style RSC file.
    Parameters: fname : str.
                    File path of .rsc file.
    Returns:    rscDict : dict
                    Dictionary of keys and values in RSC file.
    Examples:
        from pysar.utils import readfile
        atr = readfile.read_roipac_rsc('filt_101120_110220_c10.unw.rsc')
    '''
    #rsc_dict = dict(np.loadtxt(fname, dtype=bytes, usecols=(0,1)).astype(str))
    #return rsc_dict
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    rscDict = {}
    for line in lines:
        key, value = line.strip().split()[0:2]
        rscDict[key] = value
    return rscDict


def read_gamma_par(fname, delimiter=':', skiprows=3, convert2roipac=True):
    '''Read GAMMA .par/.off file into a python dictionary structure.
    Parameters: fname : str. 
                    File path of .par, .off file.
                delimiter : str, optional
                    String used to separate values.
                skiprows : int, optional
                    Skip the first skiprows lines.
    Returns:    parDict : dict
                    Attributes dictionary
    '''
    parDict = {}

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
            parDict[key] = value
    f.close()

    parDict = attribute_gamma2roipac(parDict)

    return parDict


def read_isce_xml(fname):
    '''Read ISCE .xml file input a python dictionary structure.'''
    xmlDict={}
    fObj = objectify.parse(fname)
    root = fObj.getroot()

    for child in root.findall('property'):
        xmlDict[child.attrib['name']] = str(child.value)

    ## Read lat/lon info for geocoded file
    try:
        comp1 = root.find("./component[@name='coordinate1']")
        x_step = comp1.find("./property[@name='delta']/value").text
        if x_step not in  ['1','-1']:
            xmlDict['X_STEP']  = x_step
            xmlDict['X_FIRST'] = comp1.find("./property[@name='startingvalue']/value").text
            xmlDict['X_LAST']  = comp1.find("./property[@name='endingvalue']/value").text
    except: pass

    try:
        comp2 = root.find("./component[@name='coordinate2']")
        y_step = comp2.find("./property[@name='delta']/value").text
        if y_step not in ['1','-1']:
            xmlDict['Y_STEP']  = y_step
            xmlDict['Y_FIRST'] = comp2.find("./property[@name='startingvalue']/value").text
            xmlDict['Y_LAST']  = comp2.find("./property[@name='endingvalue']/value").text
    except: pass

    xmlDict = attribute_isce2roipac(xmlDict)
    return xmlDict


def standardize_metadata(metaDict, standardMetadatKeys):
    metaDict_standard = {}
    for k in metaDict.keys():
        if k in standardMetadataKeys.keys():
            metaDict_standard[standardMetadatKeys[k]] = metaDict[k]
        else:
            metaDict_standard[k] = metaDict[k]
    return metaDict_standard


def attribute_gamma2roipac(par_dict_in):
    '''Convert Gamma par attribute into ROI_PAC format'''
    par_dict = dict()
    for key, value in iter(par_dict_in.items()):
        par_dict[key] = value

    # Length - number of rows
    for key in par_dict.keys():
        if any(key.startswith(i) for i in ['azimuth_lines','nlines','az_samp','interferogram_azimuth_lines']):
            par_dict['LENGTH'] = par_dict[key]

    # Width - number of columns
    for key in par_dict.keys():
        if any(key.startswith(i) for i in ['width','range_samp','interferogram_width']):
            par_dict['WIDTH'] = par_dict[key]

    # WAVELENGTH
    speed_of_light = 299792458.0   # meter/second
    key = 'radar_frequency'
    if key in par_dict.keys():
        par_dict['WAVELENGTH'] = str(speed_of_light/float(par_dict[key]))

    # HEIGHT & EARTH_RADIUS
    key = 'earth_radius_below_sensor'
    if key in par_dict.keys():
        par_dict['EARTH_RADIUS'] = par_dict[key]

        key2 = 'sar_to_earth_center'
        if key2 in par_dict.keys():
            par_dict['HEIGHT'] = str(float(par_dict[key2]) - float(par_dict[key]))

    # STARTING_RANGE
    key = 'near_range_slc'
    if key in par_dict.keys():
        par_dict['STARTING_RANGE'] = par_dict[key]

    # PLATFORM
    key = 'sensor'
    if key in par_dict.keys():
        par_dict['PLATFORM'] = par_dict[key]

    # ORBIT_DIRECTION
    key = 'heading'
    if key in par_dict.keys():
        value = float(par_dict[key])
        if 270 < value < 360 or -90 < value < 90:
            par_dict['ORBIT_DIRECTION'] = 'ascending'
        else:
            par_dict['ORBIT_DIRECTION'] = 'descending'

        par_dict['HEADING'] = str(value)


    ##### attributes in geo coordinates
    key = 'corner_lat'
    if key in par_dict.keys():
        par_dict['Y_FIRST'] = par_dict[key]

    key = 'corner_lon'
    if key in par_dict.keys():
        par_dict['X_FIRST'] = par_dict[key]

    ##### Optional attributes for PySAR from ROI_PAC
    # ANTENNA_SIDE
    key = 'azimuth_angle'
    if key in par_dict.keys():
        value = float(par_dict[key])
        if 0 < value < 180:
            par_dict['ANTENNA_SIDE'] = '-1'
        else:
            par_dict['ANTENNA_SIDE'] = '1'

    return par_dict


def attribute_isce2roipac(metaDict, dates=[], baselineDict={}):
    '''Convert ISCE xml attribute into ROI_PAC format'''

    rscDict={}
    for key in metaDict.keys():
        rscDict[key] = str(metaDict[key]).strip().split()[0]

    rscDict['WIDTH'] = rscDict['width']
    rscDict['LENGTH'] = rscDict['length']

    rscDict['PROCESSOR'] = 'isce'
    rscDict['INSAR_PROCESSOR'] = 'isce'
    rscDict['PLATFORM'] = 'Sentinel1'

    rscDict['ANTENNA_SIDE'] = '-1'
    if 'passDirection' in rscDict.keys():
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

def attribute_envi2roipac(metaDict):
    '''Convert ISCE xml attribute into ROI_PAC format'''

    rscDict={}
    for key in metaDict.keys():
        rscDict[key] = str(metaDict[key]).strip().split()[0]

    enviDataType = rscDict['data type']
    if enviDataType == '4':
        rscDict['DATA_TYPE'] = 'float32'
    return rscDict


#########################################################################
def read_float32(fname, box=None, byte_order='l'):
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

    atr = read_attribute(fname)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0,0,width,length]

    data_type = 'f4'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f4'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*2*width).reshape(box[3],2*width)
    amplitude = data[box[1]:box[3], box[0]:box[2]]
    phase = data[box[1]:box[3], width+box[0]:width+box[2]]

    return amplitude, phase, atr


def read_real_float64(fname, box=None, byte_order='l'):
    '''Read real float64/double data matrix, i.e. isce lat/lon.rdr
    '''
    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f8'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f8'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[0]:box[2]]
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
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'c8'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>c8'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[0]:box[2]]

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
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f4'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>f4'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[0]:box[2]]
    return data, atr


def read_complex_int16(fname, box=None, byte_order='l', cpx=False):
    '''Read complex int 16 data matrix, i.e. GAMMA SCOMPLEX file (.slc)
    
    Gamma file: .slc
    
    Inputs:
       file: complex data matrix (cpx_int16)
       box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
    Example:
       data,rsc = read_complex_int16('100102.slc')
       data,rsc = read_complex_int16('100102.slc',(100,1200,500,1500))
    '''

    atr = read_attribute(fname)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0,0,width,length]

    data_type = 'i2'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>i2'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*2*width).reshape(box[3],2*width)
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


def read_real_int16(fname, box=None, byte_order='l'):
    atr = read_attribute(fname)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0,0,width,length]

    data_type = 'i2'
    if byte_order in ['b','big','big-endian','ieee-be']:
        data_type = '>i2'

    data = np.fromfile(fname, dtype=data_type, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[0]:box[2]]
    return data, atr


def read_bool(fname, box=None):
    '''Read binary file with flags, 1-byte values with flags set in bits
    For ROI_PAC .flg, *_snap_connect.byt file.    
    '''
    # Read attribute
    if fname.endswith('_snap_connect.byt'):
        rscFile = fname.split('_snap_connect.byt')[0]+'.unw.rsc'
    else:
        rscFile = fname+'.rsc'
    atr = read_attribute(rscFile.split('.rsc')[0])
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0,0,width,length]
    
    data = np.fromfile(fname, dtype=bool, count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3], box[0]:box[2]]
    return data, atr


def read_GPS_USGS(fname):  
    yyyymmdd= np.loadtxt(fname,dtype=bytes,usecols = (0,1)).astype(str)[:,0]
    YYYYMMDD=[]
    for y in yyyymmdd:
        YYYYMMDD.append(y)
    data=np.loadtxt(fname,usecols = (1,2,3,4))
    dates=data[:,0]
    north=np.array(data[:,1])
    east=np.array(data[:,2])
    up=np.array(data[:,3])
 
    return east,north,up,dates,YYYYMMDD


#########################################################################


