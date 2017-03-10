#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

#########################################################################################
#                                                                                       #
# The empiriocal model in this program to correct the Local Oscilator Frequency Decay   #
# of Envisat ASAR instrument was suggested by Petar Marinkovic and Yngvar Larsen, 2013. #
#                                                                                       #
#########################################################################################

import os
import sys
import h5py
#import getopt
import numpy as np
import time
import datetime

def Usage():
  print '''
  *****************************************************************
  Applying an empirical model to correct the Local Oscilator Drift 
  of Envisat ASAR instrument. The empiriocal model was suggested 
  by Petar Marinkovic and Yngvar Larsen, 2013.

  Usage:

      lod.py 'time-series in radar coordinate' outname

  Example:
  
      lod.py timeseries.h5
 
      lod.py timeseries.h5 timeseries_LODcor.h5


  *****************************************************************
  '''

def tempbase(dateList):
   tbase=[]
   d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
   for ni in range(len(dateList)):
     d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
     diff = d2-d1
    # print diff
     tbase.append(diff.days)
   return tbase

def correct_igram_LOD(h5file,igram,unw):
   width=float(h5file['interferograms'][igram].attrs['WIDTH'])
   length=float(h5file['interferograms'][igram].attrs['FILE_LENGTH'])
   range_resolution=float(h5file['interferograms'][igram].attrs['RANGE_PIXEL_SIZE'])
   wvl=float(h5file['interferograms'][igram].attrs['WAVELENGTH'])
   r=np.linspace(0,width-1,width)
   R=range_resolution*r*(3.87e-7)
   Ramp=np.tile(R,[length,1])
   yref=int(h5file['interferograms'][igram].attrs['ref_y'])
   xref=int(h5file['interferograms'][igram].attrs['ref_x'])

   date1,date2=h5file['interferograms'][igram].attrs['DATE12'].split('-')
   date1=YYYYMMDD2years(date1)
   date2=YYYYMMDD2years(date2)
   dt=date1-date2;
   Ramp=Ramp*dt*(4*np.pi/wvl)
   Ramp=Ramp-Ramp[yref][xref]    
   unw=unw-Ramp
   return unw

def YYYYMMDD2years(d):
  if d[0] == '9':
      d = '19'+d
  else:
      d = '20'+d

  dy = datetime.datetime(*time.strptime(d,"%Y%m%d")[0:5])
  dyy=np.float(dy.year) + np.float(dy.month-1)/12 + np.float(dy.day-1)/365
  return dyy


def main(argv):

  try:
    File = argv[0]
  except:
    Usage();sys.exit(1)

  try:
    outName=argv[2]
  except:
    outName=File.split('.')[0]+'_LODcor.h5'

  h5file=h5py.File(File,'r')
  k=h5file.keys()
  
  if k[0]=='velocity':
    dset = h5file[k[0]].get(k[0])
    data=dset[0:dset.shape[0],0:dset.shape[1]]
    width=dset.shape[1]
    length=dset.shape[0]
    range_resolution=float(h5file[k[0]].attrs['RANGE_PIXEL_SIZE'])
    r=np.linspace(1,width,width)
    R=range_resolution*r*(3.87e-7)
    Ramp=np.tile(R,[length,1])
    yref=int(h5file[k[0]].attrs['ref_y'])
    xref=int(h5file[k[0]].attrs['ref_x'])
    Ramp=Ramp-Ramp[yref][xref]
    data=data-Ramp
    
    try:
      outName=argv[1]
    except:
      outName=File.split('.')[0]+'_LODcor.h5'

    h5velocity = h5py.File(outName,'w')
    group=h5velocity.create_group('velocity')
    dset = group.create_dataset(os.path.basename('velocity'), data=data, compression='gzip')
    for key , value in h5file[k[0]].attrs.iteritems():
        group.attrs[key]=value
    h5velocity.close()
    h5file.close()
    
  elif 'timeseries' in k:
    width=float(h5file['timeseries'].attrs['WIDTH'])
    length=float(h5file['timeseries'].attrs['FILE_LENGTH'])
    range_resolution=float(h5file['timeseries'].attrs['RANGE_PIXEL_SIZE'])
    r=np.linspace(0,width-1,width)
    R=range_resolution*r*(3.87e-7)
    Ramp=np.tile(R,[length,1])
    yref=int(h5file['timeseries'].attrs['ref_y'])
    xref=int(h5file['timeseries'].attrs['ref_x'])
    dateList = h5file['timeseries'].keys() 
    tbase=tempbase(dateList)
    t1=datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5]) 

    h5timeseries = h5py.File(outName,'w')
    group = h5timeseries.create_group('timeseries')
    for date in dateList:
       print date
       dset = h5file['timeseries'].get(date)
       data = dset[0:dset.shape[0],0:dset.shape[1]]
       
       t=datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])
       dt=(t-t1)
       dt=float(dt.days)/365.0
       Rampt=Ramp*dt
       Rampt=Rampt-Rampt[yref][xref]
       dset = group.create_dataset(date, data=data-Rampt, compression='gzip')
    
    for key,value in h5file['timeseries'].attrs.iteritems():
      group.attrs[key] = value

    try:
        dset1 = h5file['mask'].get('mask')
        group=h5timeseries.create_group('mask')
        dset = group.create_dataset('mask', data=dset1, compression='gzip')
    except: pass
    h5timeseries.close()
    h5file.close()    
  
  elif 'interferograms' in k:
    print 'LOD correction for interferograms'
    h5out = h5py.File(outName,'w')
    gg=h5out.create_group('interferograms')
    igramList=h5file['interferograms'].keys()
    for igram in igramList:
      print igram
      unwSet = h5file['interferograms'][igram].get(igram)
      unw = unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]
      platform=h5file['interferograms'][igram].attrs['PLATFORM']
      if platform in ['ENVISAT','Env','env','Envisat','envisat']:
         unw=correct_igram_LOD(h5file,igram,unw)         
      else:
         print 'No LOD correction for '+platform
      group = gg.create_group(igram)
      dset = group.create_dataset(igram, data=unw, compression='gzip')
      for key, value in h5file['interferograms'][igram].attrs.iteritems():
          group.attrs[key] = value      

    try:
        mask = h5file['mask'].get('mask')
        gm = h5out.create_group('mask')
        dset = gm.create_dataset('mask', data=mask, compression='gzip')    
    except: pass

    try:
       meanCoherence = h5file['meanCoherence'].get('meanCoherence')
       gc = h5out.create_group('meanCoherence')
       dset = gc.create_dataset('meanCoherence', data=meanCoherence, compression='gzip')
    except:
       print 'The Loaded interferograms does not contain the average coherence'

if __name__ == '__main__':
     
  main(sys.argv[1:])



