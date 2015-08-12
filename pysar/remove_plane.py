#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import numpy as np
import h5py
import sys
import _remove_surface as rm

######################################
def Usage():
  print '''
********************************************************
********************************************************

    usage:

          remove_plane.py file  method Maskfile

     file: interferograms or time-series saved in HDF5 file format.

     method: quadratic, plane, quardatic_range, quadratic_azimiuth, plane_range, plane_azimuth

     Maskfile: a mask file with 0 values for those pixels which are not considered in
               plane estimation.


    example:
          
          remove_plane.py  timeseries.h5 plane Mask.h5
          remove_plane.py  LoadedData_SanAndreasT356EnvD.h5 plane Mask.h5
          remove_plane.py  LoadedData_SanAndreasT356EnvD.h5 quadratic_range Mask.h5
          remove_plane.py  timeseries.h5 quadratic_azimuth Mask.h5

********************************************************
********************************************************
  '''

######################################
def main(argv):
  try:
    igramsFile = argv[0]
    surfType=argv[1]
    
  except:
    Usage() ; sys.exit(1)

  h5file = h5py.File(igramsFile)
  try:
    maskFile=argv[2]
    h5Mask = h5py.File(maskFile,'r')
    kMask=h5Mask.keys()
    dset1 = h5Mask[kMask[0]].get(kMask[0])
    Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
    Masking='yes'
  except:
    Masking='no'
  
  if Masking=='no':
    try:
      dset1 = h5file['mask'].get('mask')
      Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
      Masking=='yes'
    except:
      Masking=='no'

  h5flat = h5py.File(igramsFile.split('.')[0]+'_'+surfType+'.h5','w')
  k=h5file.keys()
  if 'interferograms' in k:  

     if Masking=='no':
       igramList=h5file['interferograms'].keys()
       W=int(h5file['interferograms'][igramList[0]].attrs['WIDTH'])
       L=int(h5file['interferograms'][igramList[0]].attrs['FILE_LENGTH'])
       Mask=np.ones((L,W))

     rm.remove_surface_igrams(surfType,h5file,h5flat,Mask)
     group=h5flat.create_group('mask')
     dset = group.create_dataset('mask', data=Mask, compression='gzip')
  elif 'timeseries' in k:
     if Masking=='no':
       W=int(h5file['timeseries'].attrs['WIDTH'])
       L=int(h5file['timeseries'].attrs['FILE_LENGTH'])
       Mask=np.ones((L,W))
     rm.remove_surface_timeseries(surfType,h5file,h5flat,Mask)
     group=h5flat.create_group('mask')
     dset = group.create_dataset('mask', data=Mask, compression='gzip')
  elif 'velocity' in k:
     if Masking=='no':
       W=int(h5file['velocity'].attrs['WIDTH'])
       L=int(h5file['velocity'].attrs['FILE_LENGTH'])
       Mask=np.ones((L,W))
     rm.remove_surface_velocity(surfType,h5file,h5flat,Mask)
  else:
     print 'input should be interferograms, timeseries or velocity!'


  h5file.close()
  h5flat.close()
  
######################################
if __name__ == '__main__':

  main(sys.argv[1:])


