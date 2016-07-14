
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

import numpy as np

def read_template(file):
  '''Reads the template file into a python dictionary structure.
  
  Input:
    full path to the template file
  '''  
  template_dict = {}
  for line in open(file):
    c = line.split("=")
    if len(c) < 2 or line.startswith('%') or line.startswith('#'): 
      next #ignore commented lines or those without variables
    else: 
     # template_dict[c[0].strip()] = str.replace(c[1],'\n','').strip()
      template_dict[c[0].strip()] = str.replace(c[1],'\n','').split("#")[0].strip()
  return template_dict
  
def read_rsc_file(file):
  '''Read the .rsc file into a python dictionary structure.
  
  '''
  rsc_dict = dict(np.loadtxt(file,dtype=str))
  return rsc_dict

def read_par_file(file):
  '''Read the .par file into a python dictionary structure.

  '''
  par_dict = dict(np.loadtxt(file,dtype=str,skiprows=2))
  return par_dict

#def read_float32(file):
def read_float32(*args):
  '''Reads roi_pac data (RMG format, interleaved line by line)
        file: .unw, .cor, .hgt, .trans
        box:  4-tuple defining the left, upper, right, and lower pixel coordinate.
     Example:
        a,p,r = read_float32('100102-100403.unw')
        a,p,r = read_float32('100102-100403.unw',(100,1200,500,1500))
  '''

  file = args[0]
  rscContents = read_rsc_file(file + '.rsc')
  width  = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
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

  return amplitude, phase, rscContents

def read_complex64(file):
  '''Reads roi_pac int or slc data.
  
  Requires the file path and returns amplitude and phase
  Usage:
    amp, phase, rscDictionary = readInt('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.int')
  '''
  rscContents = read_rsc_file(file + '.rsc')
  width  = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  data = np.fromfile(file,np.complex64,length*2*width).reshape(length,width)
  amplitude = np.array([np.hypot(  data.real,data.imag)]).reshape(length,width)
  phase     = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
  return amplitude, phase, rscContents

def read_gamma_float(file):
  '''Read GAMMA FLOAT file (.mli)

  '''
  try:
     parContents = read_rsc_file(file + '.rsc')
     width  = int(parContents['WIDTH'])
     length = int(parContents['FILE_LENGTH'])
  except:
     parContents = read_par_file(file + '.par')
     width  = int(parContents['range_samples:'])
     length = int(parContents['azimuth_lines:'])
  data = np.fromfile(file,np.float32,length*width).reshape(length,width)
  return data, parContents

#def read_gamma_scomplex(file,box):
def read_gamma_scomplex(*args):
  '''Read GAMMA SCOMPLEX file (.slc)
        file: complex data matrix (cpx_int16)
        box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
     Example:
        data,rsc = read_gamma_scomplex('100102.slc')
        data,rsc = read_gamma_scomplex('100102.slc',(100,1200,500,1500))
  '''

  file = args[0]
  try:
     parContents = read_rsc_file(file + '.rsc')
     width  = int(parContents['WIDTH'])
     length = int(parContents['FILE_LENGTH'])
  except:
     parContents = read_par_file(file + '.par')
     width  = int(parContents['range_samples:'])
     length = int(parContents['azimuth_lines:'])

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
  return data_cpx, parContents

  #data = np.fromfile(file,np.int16,length*2*width).reshape(length*2,width)
  #oddindices = np.where(np.arange(length*2)&1)[0]
  #real = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  #imag = np.array([data.take(oddindices,  axis=0)]).reshape(length,width)

  #amplitude = np.array([np.hypot(  real,imag)]).reshape(length,width)
  #phase     = np.array([np.arctan2(imag,real)]).reshape(length,width)
  #return amplitude, phase, parContents

def read_dem(file):
  '''Read a roipac dem file.
  
  Input:
    roi_pac format dem file
  '''  
  rscContents = read_rsc_file(file + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])  
  dem = np.fromfile(file,dtype=np.int16).reshape(length,width)
  return dem, rscContents

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

def read_isce_xml(file):
  from lxml import etree as ET
  tree = ET.parse(file)
  root = tree.getroot()
  xmldict={}
  for child in root.findall('property'):
    attrib = child.attrib['name']
    value = child.find('value').text
    xmldict[attrib]=value
  xmldict['WIDTH']=xmldict['width']
  xmldict['FILE_LENGTH']=xmldict['length']
  xmldict['WAVELENGTH'] = '0.05546576'
  Date1=os.path.dirname(file).split('/')[-1].split('_')[0][2:]
  Date2=os.path.dirname(file).split('/')[-1].split('_')[1][2:]
  xmldict['DATE12'] = Date1 + '-' + Date2
  xmldict['DATE1'] = Date1
  xmldict['DATE2'] = Date2
  return xmldict
