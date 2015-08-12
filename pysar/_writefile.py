############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import numpy as np

def write_float32(data,outname):
   # To write an array to a binary file with float32 precision
   # Format of the binary file is same as roi_pac unw, cor, or hgt data.
  
   F=np.zeros([2*data.shape[0]*data.shape[1],1],np.float32)
   nlines=data.shape[0]
   WIDTH=data.shape[1]

   for line in range(nlines):

       F[(line)*(2*WIDTH):(line)*(2*WIDTH)+WIDTH]=np.reshape(data[line][:],[WIDTH,1])
       F[(line)*(2*WIDTH)+WIDTH : (line+1)*(WIDTH*2)]=np.reshape(data[line][:],[WIDTH,1])
   
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




