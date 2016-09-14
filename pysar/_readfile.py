
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
# Yunjun, May 2016: Add read_attributes() and 'PROCESSOR','FILE_TYPE','UNIT' attributes


import os
import sys

import numpy as np
import h5py


#########################################################################
#######################  Read Attributes  ######################
def read_attributes(File):
    ## Read attributes of input file into a dictionary
    ## Input  : file name
    ## Output : atr  - attributes dictionary

    ext = os.path.splitext(File)[1].lower()
    if not os.path.isfile(File):
        print 'Input file not existed: '+File
        print 'Current directory: '+os.getcwd()
        sys.exit(1)

    ##### PySAR
    if ext == '.h5':
        h5f = h5py.File(File,'r')
        k = h5f.keys()
        if   'interferograms' in k: k[0] = 'interferograms'
        elif 'coherence'      in k: k[0] = 'coherence'
        elif 'timeseries'     in k: k[0] = 'timeseries'
        if   k[0] in ('interferograms','coherence','wrapped'):
            attrs  = h5f[k[0]][h5f[k[0]].keys()[0]].attrs
        elif k[0] in ('dem','velocity','mask','temporal_coherence','rmse','timeseries'):
            attrs  = h5f[k[0]].attrs
        else: print 'Unrecognized h5 file key: '+k[0]

        atr = dict()
        for key, value in attrs.iteritems():  atr[key] = str(value)
        atr['PROCESSOR'] = 'pysar'
        atr['FILE_TYPE'] = str(k[0])
        #try: atr['UNIT']
        #except:
        if   k[0] in ['interferograms','wrapped']:                atr['UNIT'] = 'radian'
        elif k[0] in ['coherence','temporal_coherence','rmse']:   atr['UNIT'] = '1'
        elif k[0] in ['dem','timeseries']:                        atr['UNIT'] = 'm'
        elif k[0] in ['velocity']:                                atr['UNIT'] = 'm/yr'
        else:                                                     atr['UNIT'] = '1'

        h5f.close()

    ##### ROI_PAC
    elif os.path.isfile(File + '.rsc'):
        atr = read_roipac_rsc(File + '.rsc')
        atr['PROCESSOR'] = 'roipac'
        atr['FILE_TYPE'] = ext
        if   ext in ['.unw','.int']:
            atr['UNIT'] = 'radian'
        elif ext in ['.cor','.msk','.trans','.slc','.mli']:
            atr['UNIT'] = '1'
        elif ext in ['.dem','.hgt']:
            atr['UNIT'] = 'm'

    ##### ISCE
    elif os.path.isfile(File + '.xml'):
        atr = read_isce_xml(File + '.xml')
        atr['PROCESSOR'] = 'isce'
        atr['FILE_TYPE'] = ext
        if   ext in ['.slc','.cor']:
            atr['UNIT'] = '1'
        elif ext in ['.flat']:
            atr['UNIT'] = 'radian'

    ##### GAMMA
    elif os.path.isfile(File + '.par'):
        atr = read_gamma_par(File + '.par')
        atr['PROCESSOR'] = 'gamma'
        atr['FILE_TYPE'] = ext
        if  ext in ['.mli','.slc']:
            atr['UNIT'] = '1'
    # obselete
    elif os.path.isfile(os.path.splitext(File)[0] + '.par'):
        atr = read_gamma_par(os.path.splitext(File)[0] + '.par')
        atr['PROCESSOR'] = 'gamma'
        atr['FILE_TYPE'] = ext

    else: print 'Unrecognized file extension: '+ext; sys.exit(1)

    ##### Common Attributes
    atr['FILE_EXTENSION'] = ext

    return atr

#########################################################################
def read_template(File):
    ## Reads the template file into a python dictionary structure.
    ## Input: full path to the template file
    ##  

    template_dict = {}
    for line in open(File):
        c = line.split("=")
        if len(c) < 2 or line.startswith('%') or line.startswith('#'): 
            next #ignore commented lines or those without variables
        else: 
            template_dict[c[0].strip()] = str.replace(c[1],'\n','').split("#")[0].strip()

    return template_dict

#########################################################################
def read_roipac_rsc(File):
    ## Read ROI_PAC .rsc file into a python dictionary structure.
    ##

    rsc_dict = dict(np.loadtxt(File,dtype=str))
    return rsc_dict

#########################################################################
def read_gamma_par(File):
    ## Read GAMMA .par file into a python dictionary structure.
    ##

    par_dict = dict(np.loadtxt(File,dtype=str,skiprows=2))
    par_dict['WIDTH']       = par_dict['range_samples:']
    par_dict['FILE_LENGTH'] = par_dict['azimuth_lines:']
    return par_dict

#########################################################################
def read_isce_xml(File):
    ## Read ISCE .xml file input a python dictionary structure.
    ##

    #from lxml import etree as ET
    import xml.etree.ElementTree as ET
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
#def read_float32(File):
def read_float32(*args):
    ## Reads roi_pac data (RMG format, interleaved line by line)
    ## should rename it to read_rmg_float32()
    ##
    ## RMG format (named after JPL radar pionner Richard M. Goldstein): made
    ## up of real*4 numbers in two arrays side-by-side. The two arrays often
    ## show the magnitude of the radar image and the phase, although not always
    ## (sometimes the phase is the correlation). The length and width of each 
    ## array are given as lines in the metadata (.rsc) file. Thus the total
    ## width width of the binary file is (2*width) and length is (length), data
    ## are stored as:
    ## magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    ## magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    ## ......
    ##
    ##    file: .unw, .cor, .hgt, .trans
    ##    box:  4-tuple defining the left, upper, right, and lower pixel coordinate.
    ## Example:
    ##    a,p,r = read_float32('100102-100403.unw')
    ##    a,p,r = read_float32('100102-100403.unw',(100,1200,500,1500))

    File = args[0]
    atr = read_attributes(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))

    if   len(args)==1:     box = [0,0,width,length]
    elif len(args)==2:     box = args[1]
    else: print 'Error: only support 1/2 inputs.'; return 0

    data = np.fromfile(File,np.float32,box[3]*2*width).reshape(box[3],2*width)
    amplitude = data[box[1]:box[3],box[0]:box[2]]
    phase     = data[box[1]:box[3],width+box[0]:width+box[2]]

    #oddindices = np.where(np.arange(length*2)&1)[0]
    #data = np.fromfile(File,np.float32,length*2*width).reshape(length*2,width)
    #amplitude = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
    #phase     = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

    return amplitude, phase, atr

#########################################################################
def read_complex64(File):
    ## Read complex float 32 data matrix, i.e. roi_pac int or slc data.
    ## should rename it to read_complex_float32()
    ##
    ## Data is sotred as:
    ## real, imaginary, real, imaginary, ...
    ## real, imaginary, real, imaginary, ...
    ## ...
    ##
    ## Usage:
    ##     amp, phase, atr = read_complex64('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.int')
    ## 

    atr = read_attributes(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))

    data = np.fromfile(File,np.complex64,length*2*width).reshape(length,width)
    amplitude = np.array([np.hypot(  data.real,data.imag)]).reshape(length,width)
    phase     = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
    return amplitude, phase, atr

#########################################################################
def read_real_float32(File):
    ## Read real float 32 data matrix, i.e. GAMMA .mli file
    ##
    ## Usage:
    ##     data, atr = read_real_float32('20070603.mli')

    atr = read_attributes(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))

    data = np.fromfile(File,np.float32,length*width).reshape(length,width)
    return data, atr

#########################################################################
#def read_gamma_scomplex(File,box):
def read_complex_int16(*args):
    ## Read complex int 16 data matrix, i.e. GAMMA SCOMPLEX file (.slc)
    ##
    ## Inputs:
    ##    file: complex data matrix (cpx_int16)
    ##    box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
    ## Example:
    ##    data,rsc = read_complex_int16('100102.slc')
    ##    data,rsc = read_complex_int16('100102.slc',(100,1200,500,1500))

    File = args[0]
    atr = read_attributes(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH']))

    if   len(args)==1:     box = [0,0,width,length]
    elif len(args)==2:     box = args[1]
    else: print 'Error: only support 1/2 inputs.'; return 0
    nlines = box[3]-box[1]
    WIDTH  = box[2]-box[0]

    data = np.fromfile(File,np.int16,box[3]*2*width).reshape(box[3],2*width)
    data = data[box[1]:box[3],2*box[0]:2*box[2]].reshape(2*nlines*WIDTH,1)
    id1 = range(0,2*nlines*WIDTH,2)
    id2 = range(1,2*nlines*WIDTH,2)
    real = data[id1].reshape(nlines,WIDTH)
    imag = data[id2].reshape(nlines,WIDTH)

    data_cpx = real + imag*1j
    return data_cpx, atr

    #data = np.fromfile(File,np.int16,length*2*width).reshape(length*2,width)
    #oddindices = np.where(np.arange(length*2)&1)[0]
    #real = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
    #imag = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

    #amplitude = np.array([np.hypot(  real,imag)]).reshape(length,width)
    #phase     = np.array([np.arctan2(imag,real)]).reshape(length,width)
    #return amplitude, phase, parContents

#########################################################################
def read_dem(File):
    ## Read real int 16 data matrix, i.e. ROI_PAC .dem file.
    ## should rename it to read_real_int16()
    ##
    ## Input:
    ##     roi_pac format dem file
    ## Usage:
    ##     dem, atr = read_dem('gsi10m_30m.dem')

    atr = read_attributes(File)
    width  = int(float(atr['WIDTH']))
    length = int(float(atr['FILE_LENGTH'])) 
    dem = np.fromfile(File,dtype=np.int16).reshape(length,width)
    return dem, atr

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
def read(*args):
    ## Read one dataset, i.e. interferogram, coherence, velocity, dem ...
    ##     return 0 if failed.
  
    ## Things to think: how to output multiple dataset, like all timeseries, rg and az coord of
    ##     .trans file, multiple dataset and attrs of interferograms h5 file.
  
    ## Usage:
    ##     [ data,atr] = read(file, [box|epoch_date|epoch_num])
    ##     [rg,az,atr] = read(file, [box|epoch_date|epoch_num])   ## .trans file
    ##
    ##   Inputs:
    ##     file : file path, supported format:
    ##            PySAR   file: interferograms, timeseries, velocity, etc.
    ##            ROI_PAC file: .unw .cor .hgt .dem .trans
    ##            Gamma   file: .mli .slc
    ##            Image   file: .jpeg .jpg .png .ras .bmp
    ##
    ##     box  : 4-tuple defining the left, upper, right, and lower pixel coordinate [optional]
    ##     epoch_date : date (pair) [optional]
    ##     epoch_num  : epoch number [starting from 0, optional]
  
    ##   Outputs:
    ##     data : 2D data matrix
    ##     rg/az: 2D matrix of range and azimuth pixel location of SAR (radar coord) for each DEM pixel
    ##     atr  : attribute object
  
    ## Examples:
    ##     data,atr = read('velocity.h5')
    ##     data,atr = read('temporal_coherence.h5')
    ##     data,atr = read('timeseries.h5','20100120')
    ##     data,atr = read('timeseries.h5',3)
    ##     data,atr = read('LoadedData.h5','100120-110214')
    ##     data,atr = read('LoadedData.h5',3)
    ##
    ##     data,atr = read('100120-110214.unw')
    ##     data,atr = read('strm1.dem')
    ##     data,atr = read('strm1.dem.jpg')
    ##     data,atr = read('100120.mli')
    ##     rg,az,atre = read('geomap*.trans')
    ##
    ##     data,atr = read('velocity.h5',             (100,1200,500,1500))
    ##     data,atr = read('100120-110214.unw',       (100,1200,500,1500))
    ##     data,atr = read('timeseries.h5','20100120',(100,1200,500,1500))


    ########## Check Inputs ##########
    File = args[0]
    for i in range(1,len(args)):
        if   isinstance(args[i], tuple):       box        = args[i]
        elif isinstance(args[i], basestring):  epoch_date = args[i]
        elif isinstance(args[i], int):         epoch_num  = args[i]
  
    ############### Read ###############
    ext = os.path.splitext(File)[1].lower()
    atr = read_attributes(File)
    processor = atr['PROCESSOR']
  
    ##### PySAR HDF5
    if processor == 'pysar':
        h5file = h5py.File(File,'r')
        k = atr['FILE_TYPE']

        ##### Read Data
        if k == 'timeseries':
            epochList = h5file[k].keys()
            try:     epoch_num
            except:
                try: epoch_num = epochList.index(epoch_date)
                except: print 'Unrecognized / No epoch input!';  return 0;
            dset = h5file[k].get(epochList[epoch_num])

        elif k in ['interferograms','coherence','wrapped']:
            epochList = h5file[k].keys()
            try:    epoch_num
            except:
                try:
                    epoch_date
                    for i in range(len(epochList)):
                        if epoch_date in epochList[i]:   epoch_num = i
                except: print 'Unrecognized / No epoch input!';  return 0;
            dset = h5file[k][epochList[epoch_num]].get(epochList[epoch_num])

        elif k in ['dem','mask','rmse','temporal_coherence', 'velocity']:
            dset = h5file[k].get(k)
        else: print 'Unrecognized h5 file type: '+k

        ##### Crop
        try:    data = dset[box[1]:box[3],box[0]:box[2]]
        except: data = dset[:,:]

        h5file.close()
        return data, atr

    ##### ISCE
    elif processor == 'isce':
        if   ext in ['.flat']:
            amp, data, atr = read_complex64(File)
        elif ext in ['.cor']:
            data, atr = read_real_float32(File)
        elif ext in ['.slc']:
            data, pha, atr = read_complex64(File)
            ind = np.nonzero(data)
            data[ind] = np.log10(data[ind])     # dB
            atr['UNIT'] = 'dB'
  
        try: data = data[box[1]:box[3],box[0]:box[2]]
        except: pass
        return data, atr
  
    ##### ROI_PAC
    elif processor == 'roipac':
        if ext in ['.unw','.cor','.hgt']:
            try:    amp,pha,atr = read_float32(File,box)
            except: amp,pha,atr = read_float32(File)
            return pha, atr
  
        elif ext == '.dem':
            dem,atr = read_dem(File)
            try: dem = dem[box[1]:box[3],box[0]:box[2]]
            except: pass
            return dem, atr
  
        elif ext == '.int':
            amp, pha, atr = read_complex64(File)
            try: pha = pha[box[1]:box[3],box[0]:box[2]]
            except: pass
            return pha, atr
            #ind = np.nonzero(amp)
            #amp[ind] = np.log10(amp[ind])
            #atr['UNIT'] = 'dB'
            #return amp, atr
  
        elif ext == '.trans':
            try:    rg,az,atr = read_float32(file,box)
            except: rg,az,atr = read_float32(file)
            return rg, az, atr
  
    ##### Gamma
    elif processor == 'gamma':
        if ext == '.mli':
            data,atr = read_real_float32(file)
            try: data = data[box[1]:box[3],box[0]:box[2]]
            except: pass
            data = np.nanlog10(data)     # dB
            atr['UNIT'] = 'dB'
            return data, atr
  
        elif ext == '.slc':
            try:    data,atr = read_complex_int16(File,box)
            except: data,atr = read_complex_int16(File)
            data = np.nanlog10(data)     # dB
            atr['UNIT'] = 'dB'
            return data, atr

    ##### Image
    elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        atr = read_roipac_rsc(File+'.rsc')
        import Image
        data  = Image.open(File)
        try: data = data.crop(box)
        except: pass
        return data, atr
  
    else: print 'Unrecognized file format: '+ext; return 0


#########################################################################
## Not ready yet
def read_multiple(File,box=''):
    ## Read multi-temporal 2D datasets into a 3-D data stack
    ## Inputs:
    ##     File  : input file, interferograms,coherence, timeseries, ...
    ##     box   : 4-tuple defining the left, upper, right, and lower pixel coordinate [optional]
    ## Examples:
    ##     stack = stacking('timeseries.h5',(100,1200,500,1500))

    ##### File Info
    atr = readfile.read_attributes(File)
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




