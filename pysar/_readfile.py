
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
#
# Yunjun, Sep 2015: Add read_par_file()
#                   Add read_gamma_float() and read_gamma_scomplex()
# Yunjun, Oct 2015: Add box option for read_float32()
# Heresh, Nov 2015: Add ISCE xml reader
# Yunjun, Jan 2016: Add read()
# Yunjun, May 2016: Add read_attribute() and 'PROCESSOR','FILE_TYPE','UNIT' attributes


import os
import sys

import h5py
import numpy as np
import xml.etree.ElementTree as ET
from PIL import Image


#########################################################################
'''Three types of HDF5 files in PySAR
multi_group   : multiple groups with one      dataset and one attribute dict per group (Ngroup-1dset-1atr)
multi_dataset : one      group  with multiple dataset and one attribute dict per group (1group-Ndset-1atr)
single_dataset: one      group  with one      dataset and one attribute dict per gropu (1group-1dset-1atr)
'''
multi_group_hdf5_file=['interferograms','coherence','wrapped','snaphu_connect_component']
multi_dataset_hdf5_file=['timeseries']
single_dataset_hdf5_file=['dem','mask','rmse','temporal_coherence', 'velocity']

#########################################################################
#######################  Read Attributes  ######################
def read_attribute(File, epoch=''):
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
        h5f = h5py.File(File,'r')
        k = h5f.keys()
        if   'interferograms' in k: k[0] = 'interferograms'
        elif 'coherence'      in k: k[0] = 'coherence'
        elif 'timeseries'     in k: k[0] = 'timeseries'
        if   k[0] in multi_group_hdf5_file:
            if epoch:
                attrs  = h5f[k[0]][epoch].attrs
            else:
                attrs  = h5f[k[0]][h5f[k[0]].keys()[0]].attrs
        elif k[0] in multi_dataset_hdf5_file+single_dataset_hdf5_file:
            attrs  = h5f[k[0]].attrs
        else: print 'Unrecognized h5 file key: '+k[0]

        atr = dict()
        for key, value in attrs.iteritems():  atr[key] = str(value)
        atr['PROCESSOR'] = 'pysar'
        atr['FILE_TYPE'] = str(k[0])

        if k[0] == 'timeseries':
            try: atr['ref_date']
            except: atr['ref_date'] = sorted(h5f[k[0]].keys())[0]

        h5f.close()
    
    else:
        # attribute file list
        potentialRscFileList = [File+'.rsc', File.split('_snap_connect.byt')[0]+'.unw.rsc']
        rscFile = [rscFile for rscFile in potentialRscFileList if os.path.isfile(rscFile)][0]

        ##### GAMMA
        if os.path.isfile(File+'.par'):
            atr = read_gamma_par(File+'.par')
            atr['PROCESSOR'] = 'gamma'
            atr['FILE_TYPE'] = ext
        
        ##### ROI_PAC
        elif rscFile:
            atr = read_roipac_rsc(rscFile)
            atr['PROCESSOR'] = 'roipac'
            atr['FILE_TYPE'] = ext
    
        ##### ISCE
        elif os.path.isfile(File+'.xml'):
            atr = read_isce_xml(File+'.xml')
            atr['PROCESSOR'] = 'isce'
            atr['FILE_TYPE'] = ext
    
        else: print 'Unrecognized file extension: '+ext; sys.exit(1)

    # Unit - str
    if atr['FILE_TYPE'] in ['interferograms','wrapped','.unw','.int','.flat']:
        atr['UNIT'] = 'radian'
    elif atr['FILE_TYPE'] in ['timeseries','dem','.dem','.hgt']:
        atr['UNIT'] = 'm'
    elif atr['FILE_TYPE'] in ['velocity']:
        atr['UNIT'] = 'm/yr'
    else:
        atr['UNIT'] = '1'

    return atr


#########################################################################
def check_variable_name(path):
    s=path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        p0=os.getenv(s[1:])
        path=path.replace(path.split("/")[0],p0)
    return path

#########################################################################
def read_template(File):
    '''Reads the template file into a python dictionary structure.
    Input : string, full path to the template file
    Output: dictionary, pysar template content
    '''
    template_dict = {}
    for line in open(File):
        c = line.split("=")
        if len(c) < 2 or line.startswith('%') or line.startswith('#'):
            next #ignore commented lines or those without variables
        else:
            atrName  = c[0].strip()
            atrValue = str.replace(c[1],'\n','').split("#")[0].strip()
            atrValue = check_variable_name(atrValue)
            template_dict[atrName] = atrValue
    return template_dict

#########################################################################
def read_roipac_rsc(File):
    '''Read ROI_PAC .rsc file into a python dictionary structure.'''
    rsc_dict = dict(np.loadtxt(File,dtype=str))
    return rsc_dict

#########################################################################
def read_gamma_par(File):
    '''Read GAMMA .par file into a python dictionary structure.'''
    par_dict = {}
    
    f = open(File)
    next(f); next(f);
    for line in f:
        l = line.split()
        if len(l) < 2 or line.startswith('%') or line.startswith('#'):
            next #ignore commented lines or those without variables
        else:
            par_dict[l[0].strip()] = str.replace(l[1],'\n','').split("#")[0].strip()
    
    par_dict['WIDTH']       = par_dict['range_samples:']
    par_dict['FILE_LENGTH'] = par_dict['azimuth_lines:']
    
    return par_dict

#########################################################################
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

#########################################################################
def merge_attribute(atr1,atr2):
    atr = dict()
    for key, value in atr1.iteritems():  atr[key] = str(value)
    for key, value in atr2.iteritems():  atr[key] = str(value)

    return atr


#########################################################################
##############################  Read Data  ##############################
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

#########################################################################
##def read_complex64(File, real_imag=0):
def read_complex_float32(File, real_imag=False):
    '''Read complex float 32 data matrix, i.e. roi_pac int or slc data.
    should rename it to read_complex_float32()
    
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

#########################################################################
def read_real_float32(File):
    '''Read real float 32 data matrix, i.e. GAMMA .mli file
    Usage:  data, atr = read_real_float32('20070603.mli')
    '''
    atr = read_attribute(File)
    width = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))
    data = np.fromfile(File,dtype=np.float32).reshape(length,width)
    return data, atr

#########################################################################
#def read_complex_int16(*args):
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

#########################################################################
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


#########################################################################
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


#########################################################################
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
#def read(*args):
def read(File, box=(), epoch=''):
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
        if k in multi_group_hdf5_file:
            epochList = h5file[k].keys()
            epochList = sorted(epochList)

            if not epoch in epochList:
                print 'input epoch is not included in file: '+File
                print 'input epoch: '+epoch
                print 'epoch in file '+File
                print epochList

            if k in multi_dataset_hdf5_file:
                dset = h5file[k].get(epoch)
            else:
                dset = h5file[k][epoch].get(epoch)

        elif k in single_dataset_hdf5_file:
            dset = h5file[k].get(k)
        else: print 'Unrecognized h5 file type: '+k

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
        else: print 'Un-supported '+processfor+' file format: '+ext
  
        if box:  data = data[box[1]:box[3],box[0]:box[2]]
        return data, atr

    ##### ROI_PAC
    elif processor in ['roipac','gamma']:
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
                print 'read range and azimuth from '+File
                return rg, az, atr
            elif epoch in ['rg','range']:
                print 'read range from '+File
                return rg, atr
            elif epoch in ['az','azimuth']:
                print 'read azimuth from '+File
                return az, atr

        ##### Gamma
        #elif processor == 'gamma':
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

        else: print 'Un-supported '+processfor+' file format: '+ext
    else: print 'Unrecognized file format: '+ext; return 0


#########################################################################
## Not ready yet
def read_multiple(File,box=''):
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
    print 'number of epochs: '+str(epochNum)
 
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




