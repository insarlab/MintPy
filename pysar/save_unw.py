#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: update DATE12 for timeseries option
# Yunjun, Oct 2015: add coherence/wrapped option
#                   add two dates option for timeseries


import sys
import os

def Usage():
  print '''
****************************************************************
****************************************************************
   To converts the PySAR hdf5 file formats to the roipac unw format 

   Usage: 
         save_unw.py file.h5 

   Example:
         save_unw.py velocity.h5
         save_unw.py timeseries.h5 20050601
         save_unw.py timeseries.h5 20040728 20050601
         save_unw.py LoadedData.h5 filt_091225-100723-sim_HDR_8rlks_c10.unw
         save_unw.py temporal_coherence.h5


      for velocity:   the ouput will be a one year interferogram.
      for timeseries: if date is not specified, the last date will be used
                      if two dates are specified, the earlier date will be
                         used as the reference date.
****************************************************************
****************************************************************
  '''


def main(argv):
  try:    File=argv[0]
  except: Usage();sys.exit(1)

  from numpy import pi
  import h5py
  import pysar._writefile as writefile
 
  h5file=h5py.File(File,'r')
  k=h5file.keys()
  if 'interferograms' in k: k[0] = 'interferograms'
  elif 'coherence'    in k: k[0] = 'coherence'
  elif 'timeseries'   in k: k[0] = 'timeseries'

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

  elif k[0] in ['rmse','temporal_coherence','mask']:
    dset = h5file[k[0]].get(k[0])
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    if k[0] == 'temporal_coherence': outname=File.split('.')[0]+'.cor'
    else:                            outname=File.split('.')[0]+'.unw'
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]].attrs.iteritems():
      f.write(key+'    '+str(value)+'\n')
    f.close()

  elif 'timeseries' in k:
    dateList=h5file['timeseries'].keys() 
    ## Input
    if   len(sys.argv)==2:
      print 'No input date specified >>> continue with the last date'
      dateList=h5file['timeseries'].keys()
      d=dateList[-1]
    elif len(sys.argv)==3:
      d=sys.argv[2]
    elif len(sys.argv)==4:
      ds=sys.argv[2:4]; ds.sort()
      d_ref = ds[0]
      d     = ds[1]
    else: Usage(); sys.exit(1)

    ## Data Operation
    print 'reading '+d+' ... '
    dset = h5file['timeseries'].get(d)    
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    try:
      dset_ref = h5file['timeseries'].get(d_ref)
      print 'reading '+d_ref+' ... '
      data_ref = dset_ref[0:dset_ref.shape[0],0:dset_ref.shape[1]]
      data = data - data_ref
    except: pass
    wvl=float(h5file[k[0]].attrs['WAVELENGTH'])
    data=(-4*pi/wvl)*data

    ## outName
    try:      master_d = d_ref
    except:
      try:    master_d = h5file[k[0]].attrs['ref_date']
      except: master_d = h5file[k[0]].attrs['DATE']
    if len(master_d)==8:  master_d=master_d[2:8]
    if len(d)==8:         d=d[2:8]
    outname = master_d+'_'+d+'.unw'

    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]].attrs.iteritems():
      if key=='DATE12':
        f.write(key+'    '+master_d+'-'+d+'\n')
      else:
        f.write(key+'    '+str(value)+'\n')
    f.close()


  elif k[0] in ['interferograms','coherence','wrapped']:
    igramList=h5file[k[0]].keys()
    try:     d=sys.argv[2]
    except:  d=igramList[-1];   print 'No input date specified >>> continue with the last date'
    print 'reading '+d + ' ... '
    dset=h5file[k[0]][d].get(d)
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    outname=d
    print 'writing >>> '+ outname
    writefile.write_float32(data,outname)
    f = open(outname+'.rsc','w')
    for key , value in h5file[k[0]][d].attrs.iteritems():
      f.write(key+'    '+str(value)+'\n')
    f.close()    

  h5file.close()


##########################################################################

if __name__ == '__main__':

  main(sys.argv[1:])
