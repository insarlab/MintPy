#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
from numpy import pi
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib import colors
#import getopt
import h5py
#import scipy.io as sio
import _writefile as writefile

def Usage():
  print '''
****************************************************************
****************************************************************
   To converts the PySAR hdf5 file formats to the roipac unw format 
   usage: 
         save_unw.py file.h5 
        
   example:
         save_unw.py velocity.h5
         save_unw.py timeseries.h5 20050601
         save_unw.py LoadedData_ChamanT256EnvA6.h5 filt_091225-100723-sim_HDR_8rlks_c10.unw


    If the input is velocity, the ouput will be a one year interferogram.
    for the timeseries if date is not specified, by default the last date is used
         save_unw.py timeseries.h5     
****************************************************************
****************************************************************
  '''


def main(argv):
  try:
    File=argv[0]
  except:
    Usage();sys.exit(1)
 
  h5file=h5py.File(File,'r')

  k=h5file.keys()
  if 'velocity' in k:
    
    dset = h5file['velocity'].get('velocity')
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    print "converting velocity to a 1 year interferogram."
    wvl=float(h5file[k[0]].attrs['WAVELENGTH'])
    data=(-4*pi/wvl)*data

    outname=File.split('.')[0]+'.unw'
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]].attrs.iteritems():
      f.write(key+'    '+str(value)+'\n')
    f.close()

  if k[0] =='rmse' or k[0] =='temporal_coherence' or k[0]=='mask':
    dset = h5file[k[0]].get(k[0])
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    outname=File.split('.')[0]+'.unw'
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]].attrs.iteritems():
      f.write(key+'    '+str(value)+'\n')
    f.close()
  #  outname=File.split('.')[0]+'.unw'
  #  writefile.write_float32(data,outname) 
  #  f = open(outname+'.rsc','w')
  #  for key , value in h5file[k[0]].attrs.iteritems():
  #    f.write(key+'    '+str(value)+'\n')
  #  f.close()

  elif 'timeseries' in k:
    dateList=h5file['timeseries'].keys() 
    try:
      d=sys.argv[2]
    except:
      print 'No input date specified >>> continue with the last date'
      dateList=h5file['timeseries'].keys()
      d=dateList[-1]
    print 'reading '+d + ' ... '
    dset=h5file['timeseries'].get(d)    
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    wvl=float(h5file[k[0]].attrs['WAVELENGTH'])
    data=(-4*pi/wvl)*data

    outname=File.split('.')[0]+'.unw'
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]].attrs.iteritems():
      # update 'DATE12', Yunjun, Aug 3rd 2015
      if key=='DATE12':
        try:      master_d=h5file[k[0]].attrs['ref_date']
        except:   master_d=h5file[k[0]].attrs['DATE']
        if len(master_d)==8:  master_d=master_d[2:8]
        if len(d)==8:         d=d[2:8]
        f.write(key+'    '+master_d+'-'+d+'\n')
      else:
        f.write(key+'    '+str(value)+'\n')
    f.close()
  

  elif 'interferograms' in k:
    igramList=h5file['interferograms'].keys()
    try:
      d=sys.argv[2]
    except:
      print 'No input date specified >>> continue with the last date'
      d=igramList[-1]
    print 'reading '+d + ' ... '
    dset=h5file['interferograms'][d].get(d)
    data = dset[0:dset.shape[0],0:dset.shape[1]]
#    wvl=float(h5file['interferograms'][d].attrs['WAVELENGTH'])
#    data=(-4*pi/wvl)*data
    
   # outname=File.split('.')[0]+'.unw'
    outname=d
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file['interferograms'][d].attrs.iteritems():
      f.write(key+'    '+str(value)+'\n')
    f.close()    

  h5file.close()
if __name__ == '__main__':

  main(sys.argv[1:])
