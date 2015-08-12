#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
import sys 
import os
import getopt
import h5py
import numpy as np
import matplotlib.pyplot as plt
import _readfile as readfile

def Usage():
   print '''

************************************************************************
************************************************************************

   Calculates the correlation of the dem with the InSAR velocity field.
  
   Usage:
     
     correlation_with_dem.py dem velocity

   Example:

     correlation_with_dem.py radar_8rlks.hgt  velocity.h5

***********************************************************************
***********************************************************************
'''

try:
  demFile=sys.argv[1]
  File=sys.argv[2]
except:
  Usage()
  sys.exit(1)


if os.path.basename(demFile).split('.')[1]=='hgt':
       amp,dem,demRsc = readfile.read_float32(demFile)

elif os.path.basename(demFile).split('.')[1]=='dem':
       dem,demRsc = readfile.read_dem(demFile)

#amp,dem,demRsc = readfile.read_float32(demFile)
h5data = h5py.File(File)
dset = h5data['velocity'].get('velocity')
data = dset[0:dset.shape[0],0:dset.shape[1]]

try:
  suby=sys.argv[3].split(':')
  subx=sys.argv[4].split(':')
  data = data[int(suby[0]):int(suby[1]),int(subx[0]):int(subx[1])]
  dem = dem[int(suby[0]):int(suby[1]),int(subx[0]):int(subx[1])]
except:
  print 'no subset'

dem=dem.flatten(1)
data=data.flatten(1)
ndx = ~np.isnan(data)
C1=np.zeros([2,len(dem[ndx])])
C1[0][:]=dem[ndx]
C1[1][:]=data[ndx]
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print ''
print 'Correlation of the velocity with the DEM:  '+ str(np.corrcoef(C1)[0][1])
print ''
print'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print 'DEM info:'
print ''
print 'Maximum height difference (m) : ' + str(np.max(dem[ndx])-np.min(dem[ndx]))
print 'Average height (m) :'+str(np.mean(dem[ndx]))
print 'Height Std: '+str(np.std(dem[ndx]))

