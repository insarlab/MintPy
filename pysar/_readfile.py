
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
def read_attributes(file):
  ## Read attributes of input file into a dictionary
  ## Input  : file name
  ## Output : atr  - attributes dictionary

  ext = os.path.splitext(file)[1].lower()
  if not os.path.isfile(file):
      print 'Input file not existed: '+file
      print 'Current directory: '+os.getcwd()
      sys.exit(1)

  ##### PySAR
  if ext == '.h5':
      h5f = h5py.File(file,'r')
      k = h5f.keys()
      if 'interferograms' in k: k[0] = 'interferograms'
      elif 'coherence'    in k: k[0] = 'coherence'
      elif 'timeseries'   in k: k[0] = 'timeseries'
      if   k[0] in ('interferograms','coherence','wrapped'):
          attrs  = h5f[k[0]][h5f[k[0]].keys()[0]].attrs
      elif k[0] in ('dem','velocity','mask','temporal_coherence','rmse','timeseries'):
          attrs  = h5f[k[0]].attrs
      else: print 'Unrecognized h5 file key: '+k[0]

      atr = dict()
      for key, value in attrs.iteritems():  atr[key] = str(value)
      atr['PROCESSOR'] = 'pysar'
      atr['FILE_TYPE']      = str(k[0])
      #try: atr['UNIT']
      #except:
      if   k[0] in ['interferograms','wrapped']:                atr['UNIT'] = 'radian'
      elif k[0] in ['coherence','temporal_coherence','rmse']:   atr['UNIT'] = '1'
      elif k[0] in ['dem','timeseries']:                        atr['UNIT'] = 'm'
      elif k[0] in ['velocity']:                                atr['UNIT'] = 'm/yr'
      else:                                                     atr['UNIT'] = '1'

      h5f.close()

 ##### ROI_PAC
  elif os.path.isfile(file + '.rsc'):
      atr = read_roipac_rsc(file + '.rsc')
      atr['PROCESSOR'] = 'roipac'
      atr['FILE_TYPE'] = ext
      if   ext in ['.unw','.int']:
          atr['UNIT'] = 'radian'
      elif ext in ['.cor','.msk','.trans','.slc','.mli']:
          atr['UNIT'] = '1'
      elif ext in ['.dem','.hgt']:
          atr['UNIT'] = 'm'

  ##### ISCE
  elif os.path.isfile(file + '.xml'):
      atr = read_isce_xml(file + '.xml')
      atr['PROCESSOR'] = 'isce'
      atr['FILE_TYPE'] = ext
      if   ext in ['.slc','.cor']:
          atr['UNIT'] = '1'
      elif ext in ['.flat']:
          atr['UNIT'] = 'radian'

  ##### GAMMA
  elif os.path.isfile(file + '.par'):
      atr = read_gamma_par(file + '.par')
      atr['PROCESSOR'] = 'gamma'
      atr['FILE_TYPE'] = ext
      if  ext in ['.mli','.slc']:
          atr['UNIT'] = '1'
  # obselete
  elif os.path.isfile(os.path.splitext(file)[0] + '.par'):
      atr = read_gamma_par(os.path.splitext(file)[0] + '.par')
      atr['PROCESSOR'] = 'gamma'
      atr['FILE_TYPE'] = ext

  else: print 'Unrecognized file extension: '+ext; sys.exit(1)

  ##### Common Attributes
  atr['FILE_EXTENSION'] = ext

  return atr

#########################################################################
def read_template(file):
  '''Reads the template file into a python dictionary structure.
     Input: full path to the template file
  '''  

  template_dict = {}
  for line in open(file):
    c = line.split("=")
    if len(c) < 2 or line.startswith('%') or line.startswith('#'): 
      next #ignore commented lines or those without variables
    else: 
      template_dict[c[0].strip()] = str.replace(c[1],'\n','').split("#")[0].strip()

  return template_dict

#########################################################################
def read_roipac_rsc(file):
  '''Read the .rsc file into a python dictionary structure.
  '''

  rsc_dict = dict(np.loadtxt(file,dtype=str))
  return rsc_dict

#########################################################################
def read_gamma_par(file):
  '''Read the .par file into a python dictionary structure.
  '''

  par_dict = dict(np.loadtxt(file,dtype=str,skiprows=2))
  par_dict['WIDTH']       = par_dict['range_samples:']
  par_dict['FILE_LENGTH'] = par_dict['azimuth_lines:']
  return par_dict

#########################################################################
def read_isce_xml(file):
  #from lxml import etree as ET
  import xml.etree.ElementTree as ET
  tree = ET.parse(file)
  root = tree.getroot()
  xmldict={}
  for child in root.findall('property'):
    attrib = child.attrib['name']
    value = child.find('value').text
    xmldict[attrib]=value

  xmldict['WIDTH']       = xmldict['width']
  xmldict['FILE_LENGTH'] = xmldict['length']
  #xmldict['WAVELENGTH'] = '0.05546576'
  #Date1=os.path.dirname(file).split('/')[-1].split('_')[0][2:]
  #Date2=os.path.dirname(file).split('/')[-1].split('_')[1][2:]
  #xmldict['DATE12'] = Date1 + '-' + Date2
  #xmldict['DATE1'] = Date1
  #xmldict['DATE2'] = Date2
  return xmldict

#########################################################################
#def read_float32(file):
def read_float32(*args):
  '''Reads roi_pac data (RMG format, interleaved line by line)

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

        file: .unw, .cor, .hgt, .trans
        box:  4-tuple defining the left, upper, right, and lower pixel coordinate.
     Example:
        a,p,r = read_float32('100102-100403.unw')
        a,p,r = read_float32('100102-100403.unw',(100,1200,500,1500))
  '''

  file = args[0]
  atr = read_attributes(file)
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])

  if   len(args)==1:     box = [0,0,width,length]
  elif len(args)==2:     box = args[1]
  else: print 'Error: only support 1/2 inputs.'; return

  data = np.fromfile(file,np.float32,box[3]*2*width).reshape(box[3],2*width)
  amplitude = data[box[1]:box[3],box[0]:box[2]]
  phase     = data[box[1]:box[3],width+box[0]:width+box[2]]

  #oddindices = np.where(np.arange(length*2)&1)[0]
  #data = np.fromfile(file,np.float32,length*2*width).reshape(length*2,width)
  #amplitude = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  #phase     = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

  return amplitude, phase, atr

#########################################################################
def read_complex64(file):
  '''Reads roi_pac int or slc data.
     Requires the file path and returns amplitude and phase
     Usage:
       amp, phase, rscDictionary = readInt('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.int')
  '''

  atr = read_attributes(file)
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])

  data = np.fromfile(file,np.complex64,length*2*width).reshape(length,width)
  amplitude = np.array([np.hypot(  data.real,data.imag)]).reshape(length,width)
  phase     = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
  return amplitude, phase, atr

#########################################################################
def read_real_float32(file):
  '''Read GAMMA FLOAT file (.mli)
  '''

  atr = read_attributes(file)
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])

  data = np.fromfile(file,np.float32,length*width).reshape(length,width)
  return data, atr

#########################################################################
#def read_gamma_scomplex(file,box):
def read_complex_int16(*args):
  '''Read GAMMA SCOMPLEX file (.slc)
        file: complex data matrix (cpx_int16)
        box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
     Example:
        data,rsc = read_gamma_scomplex('100102.slc')
        data,rsc = read_gamma_scomplex('100102.slc',(100,1200,500,1500))
  '''

  file = args[0]
  atr = read_attributes(file)
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])

  if   len(args)==1:     box = [0,0,width,length]
  elif len(args)==2:     box = args[1]
  else: print 'Error: only support 1/2 inputs.'; return
  nlines = box[3]-box[1]
  WIDTH  = box[2]-box[0]

  data = np.fromfile(file,np.int16,box[3]*2*width).reshape(box[3],2*width)
  data = data[box[1]:box[3],2*box[0]:2*box[2]].reshape(2*nlines*WIDTH,1)
  id1 = range(0,2*nlines*WIDTH,2)
  id2 = range(1,2*nlines*WIDTH,2)
  real = data[id1].reshape(nlines,WIDTH)
  imag = data[id2].reshape(nlines,WIDTH)

  data_cpx = real + imag*1j
  return data_cpx, atr

  #data = np.fromfile(file,np.int16,length*2*width).reshape(length*2,width)
  #oddindices = np.where(np.arange(length*2)&1)[0]
  #real = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  #imag = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

  #amplitude = np.array([np.hypot(  real,imag)]).reshape(length,width)
  #phase     = np.array([np.arctan2(imag,real)]).reshape(length,width)
  #return amplitude, phase, parContents

#########################################################################
def read_dem(file):
  '''Read a roipac dem file.
     Input:
       roi_pac format dem file
  '''  

  atr = read_attributes(file)
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])  
  dem = np.fromfile(file,dtype=np.int16).reshape(length,width)
  return dem, atr

#########################################################################
def read_GPS_USGS(file):  
   yyyymmdd= np.loadtxt(file,dtype=str,usecols = (0,1))[:,0]
   YYYYMMDD=[]
   for y in yyyymmdd:
      YYYYMMDD.append(y)
   data=np.loadtxt(file,usecols = (1,2,3,4))
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
  file = args[0]
  for i in range(1,len(args)):
      if   isinstance(args[i], tuple):       box        = args[i]
      elif isinstance(args[i], basestring):  epoch_date = args[i]
      elif isinstance(args[i], int):         epoch_num  = args[i]

  ############### Read ###############
  ext = os.path.splitext(file)[1].lower()
  atr = read_attributes(file)
  processor = atr['PROCESSOR']
  type      = atr['FILE_TYPE']

  ##### PySAR HDF5
  if processor == 'pysar':
      h5file = h5py.File(file,'r')
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
      if   type in ['.flat']:
          amp, data, atr = read_complex64(file)
      elif type in ['.cor']:
          data, atr = read_real_float32(file)
      elif type in ['.slc']:
          data, pha, atr = read_complex64(file)
          ind = np.nonzero(data)
          data[ind] = np.log10(data[ind])     # dB
          atr['UNIT'] = 'dB'

      try: data = data[box[1]:box[3],box[0]:box[2]]
      except: pass
      return data, atr

  ##### ROI_PAC
  elif processor == 'roipac':
      if ext in ['.unw','.cor','.hgt']:
          try:    amp,pha,atr = read_float32(file,box)
          except: amp,pha,atr = read_float32(file)
          return pha, atr

      elif ext == '.dem':
          dem,atr = read_dem(file)
          try: dem = dem[box[1]:box[3],box[0]:box[2]]
          except: pass
          return dem, atr

      elif ext == '.int':
          amp, pha, atr = read_complex64(file)
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
          try:    data,atr = read_complex_int16(file,box)
          except: data,atr = read_complex_int16(file)
          data = np.nanlog10(data)     # dB
          atr['UNIT'] = 'dB'
          return data, atr

  ##### Image
  elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
      atr = read_roipac_rsc(file+'.rsc')
      import Image
      data  = Image.open(file)
      try: data = data.crop(box)
      except: pass
      return data, atr

  else: print 'Unrecognized file format: '+ext; return 0


