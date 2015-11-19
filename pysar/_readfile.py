
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

def read_float32(file):
  '''Reads roi_pac unw, cor, or hgt data.
  
  Requires the file path and returns amplitude and phase
  Usage:
    amplitude, phase, rscDictionary = readUnw('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.unw')
  '''
  rscContents = read_rsc_file(file + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  oddindices = np.where(np.arange(length*2)&1)[0]
  data = np.fromfile(file,np.float32,length*2*width)
#  print np.shape(data)
  data=data.reshape(length*2,width)
  amplitude = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  phase = np.array([data.take(oddindices,axis=0)]).reshape(length,width)
  return amplitude, phase, rscContents

def read_complex64(file):
  '''Reads roi_pac int or slc data.
  
  Requires the file path and returns amplitude and phase
  Usage:
    amp, phase, rscDictionary = readInt('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.int')
  '''
  rscContents = read_rsc_file(file + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  data = np.fromfile(file,np.complex64,length*2*width).reshape(length,width)
  amplitude = np.array([np.hypot(data.real,data.imag)]).reshape(length,width)
  phase = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
  return amplitude, phase, rscContents

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
   print "isce"
