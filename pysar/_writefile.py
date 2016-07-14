############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Sep 2015: Add write_gamma_float() and write_gamma_scomplex()
# Yunjun, Oct 2015: Add support for write_float32(amp, phase, outname)

import numpy as np

#def write_float32(data,outname):
def write_float32(*args):
   # To write an array to a binary file with float32 precision
   # Format of the binary file is same as roi_pac unw, cor, or hgt data.
   # Exmaple:
   #         write_float32(phase, outname)
   #         write_float32(amp, phase, outname)

   if len(args)==2:
      amp     = args[0]
      pha     = args[0]
      outname = args[1]
   elif len(args)==3:
      amp     = args[0]
      pha     = args[1]
      outname = args[2]
   else:
      print 'Error while getting args: support 2/3 args only.'
      return

   nlines = pha.shape[0]
   WIDTH  = pha.shape[1]
   F=np.zeros([2*nlines*WIDTH,1],np.float32)

   for line in range(nlines):
       F[(2*WIDTH)*(line) :       (2*WIDTH)*(line)+WIDTH]=np.reshape(amp[line][:],[WIDTH,1])
       F[(2*WIDTH)*(line)+WIDTH : (2*WIDTH)*(line+1)]    =np.reshape(pha[line][:],[WIDTH,1])

   F.tofile(outname)

def write_complex64(data,outname):   
  # writes roi_pac .int data

   nlines=data.shape[0]
   WIDTH=data.shape[1]
   R=np.cos(data)
   Im=np.sin(data)
#   F=np.zeros([2*nlines*WIDTH,1],np.complex64) 
   F=np.zeros([2*nlines*WIDTH,1],np.float32)  
   id1=range(0,2*nlines*WIDTH,2)
   id2=range(1,2*nlines*WIDTH,2)
   F[id1]=np.reshape(R,(nlines*WIDTH,1))
   F[id2]=np.reshape(Im,(nlines*WIDTH,1))
   F.tofile(outname)

def write_dem(data,outname):
   data=np.array(data,dtype=np.int16)
   data.tofile(outname)

def write_gamma_float(data,outname):
   data=np.array(data,dtype=np.float32)
   data.tofile(outname)

def write_gamma_scomplex(data,outname):
   # write gamma scomplex data, i.e. .slc file.
   nlines = data.shape[0]
   WIDTH  = data.shape[1]
   id1 = range(0,2*nlines*WIDTH,2)
   id2 = range(1,2*nlines*WIDTH,2)
   F=np.zeros([2*nlines*WIDTH,1],np.int16)
   F[id1]=np.reshape(data.real,(nlines*WIDTH,1))
   F[id2]=np.reshape(data.imag,(nlines*WIDTH,1))
   F.tofile(outname)




