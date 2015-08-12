#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import getopt
import h5py
import _readfile  as readfile 
import _writefile as writefile


def Usage():

   print '''

   ***************************************************************
   Usage:

    multi_looking.py  file  azimuth_looks  range_looks

    file: PySAR h5 files [interferogram, time-series, velocity], 
          roi_pac files [dem, unw, hgt, ]
          or image files [jpeg, jpg, png]

   Example:
       
     multi_looking.py  velocity.h5 15 15
     multi_looking.py  Loaded_Igrams.h5 15 20
     multi_looking.py  timeseries.h5 15 15
     multi_looking.py  SanAndreas.dem 10 10
     multi_looking.py  SanAndreas.dem.jpeg 10 10

   ***************************************************************


'''

def multilook(ifg,lksy,lksx):

    
    rows,cols=ifg.shape
    rows_lowres=np.floor(rows/lksy)
    cols_lowres=np.floor(cols/lksx)


   # %ifg_lowres=NaN(rows_lowres,cols_lowres)
    
    thr = np.floor(lksx*lksy/2)
    
    
    ifg_Clowres=np.zeros((rows,cols_lowres))
    ifg_lowres=np.zeros((rows_lowres,cols_lowres))
    
    
    for c in range(int(cols_lowres)):
        ifg_Clowres[:,c]=np.sum(ifg[:,(c)*lksx:(c+1)*lksx],1)
        
    
    
    for r in range(int(rows_lowres)):
        
        ifg_lowres[r,:]=np.sum(ifg_Clowres[(r)*lksy:(r+1)*lksy,:],0)
    
    
    
    ifg_lowres=ifg_lowres/(lksy*lksx) 
    return ifg_lowres

def main(argv):

  try:  
    file=argv[0]
    alks=float(argv[1])
    rlks=float(argv[2])
  except:
    Usage();sys.exit(1)



  ext = os.path.splitext(file)[1]


  outName=file.split('.')[0]+'_a'+str(int(alks))+'lks_r'+str(int(rlks))+'lks'+ext
  if ext == '.int' or ext == '.slc':
    a,p,r = readfile.read_complex64(file)
    plks=multilook(p,alks,rlks)
    alks=multilook(a,alks,rlks)


    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    try:
       r['Y_STEP']=str(float(r['Y_STEP'])*alks)
       r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()   

  elif ext == '.unw' or ext == '.cor' or ext == '.hgt':
    a,p,r = readfile.read_float32(file)
    plks=multilook(p,alks,rlks)
    alks=multilook(a,alks,rlks)
    
    writefile.write_float32(plks,outName)

    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    
    try:
       r['Y_STEP']=str(float(r['Y_STEP'])*alks)
       r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext == ('.dem'):
    d,r = readfile.read_dem(file)
    dlks=multilook(d,alks,rlks)

    print 'writing '+outName
    writefile.write_dem(dlks,outName)
    
    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)

    try:
      r['Y_STEP']=str(float(r['Y_STEP'])*alks)
      r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
      Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext in ['.jpeg','jpg','png']:

    import Image
    im = Image.open(file)

    width = im.size[0] / int(rlks)
    height = im.size[1] / int(alks)

    imlks = im.resize((width, height), Image.NEAREST)
    print 'writing ' + outName
    imlks.save(outName)

    try:
      r=readfile.read_rsc_file(file+'.rsc')
    except:
      sys.exit(1)

    r['FILE_LENGTH']=str(height)
    r['WIDTH']=str(width)
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    try:
      r['Y_STEP']=str(float(r['Y_STEP'])*alks)
      r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
      Geo=0
    
    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()


  elif ext == ('.h5'):

    h5file=h5py.File(file,'r')
   # outName=file.split('.')[0]+'_a'+str(int(alks))+'lks_r'+str(int(rlks))+'lks.h5'
    h5file_lks=h5py.File(outName,'w')
  
    if 'interferograms' in h5file.keys():
      print 'Multilooking the interferograms'
      gg = h5file_lks.create_group('interferograms')
      igramList=h5file['interferograms'].keys()
      for igram in igramList:
        print igram
        unw = h5file['interferograms'][igram].get(igram)
        unwlks=multilook(unw,alks,rlks)
        group = gg.create_group(igram)
        dset = group.create_dataset(igram, data=unwlks, compression='gzip')
        for key, value in h5file['interferograms'][igram].attrs.iteritems():
            group.attrs[key] = value
        group.attrs['WIDTH']=unwlks.shape[1]
        group.attrs['FILE_LENGTH']=unwlks.shape[0]
        try:
           group.attrs['Y_STEP']=alks*float(group.attrs['Y_STEP'])
           group.attrs['X_STEP']=rlks*float(group.attrs['X_STEP'])
        except:
           group.attrs['AZIMUTH_PIXEL_SIZE']=alks*float(group.attrs['AZIMUTH_PIXEL_SIZE'])
           group.attrs['RANGE_PIXEL_SIZE']=rlks*float(group.attrs['RANGE_PIXEL_SIZE'])         

      dset1=h5file['mask'].get('mask')
      mask=dset1[0:dset1.shape[0],0:dset1.shape[1]]
      masklks=multilook(mask,alks,rlks)
      group=h5file_lks.create_group('mask')
      dset = group.create_dataset('mask', data=masklks, compression='gzip')

    elif 'timeseries' in h5file.keys():
      print 'Multilooking the time-series'
      group = h5file_lks.create_group('timeseries')
      dateList=h5file['timeseries'].keys()
      for d in dateList:
        print d
        unw = h5file['timeseries'].get(d)
        unwlks=multilook(unw,alks,rlks)
        dset = group.create_dataset(d, data=unwlks, compression='gzip')      

      for key,value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value
      group.attrs['WIDTH']=unwlks.shape[1]
      group.attrs['FILE_LENGTH']=unwlks.shape[0]

      try:
           group.attrs['Y_STEP']=alks*float(group.attrs['Y_STEP'])
           group.attrs['X_STEP']=rlks*float(group.attrs['X_STEP'])
      except:
           group.attrs['AZIMUTH_PIXEL_SIZE']=alks*float(group.attrs['AZIMUTH_PIXEL_SIZE'])
           group.attrs['RANGE_PIXEL_SIZE']=rlks*float(group.attrs['RANGE_PIXEL_SIZE'])


      try:
        dset1 = h5file['mask'].get('mask')
        Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
        Masklks=multilook(Mask,alks,rlks)
        group=h5file_lks.create_group('mask')
        dset = group.create_dataset('mask', data=Masklks, compression='gzip')
      except:
        print 'Multilooked file does not include the maske'


    elif 'temporal_coherence' in h5file.keys() or 'velocity' in h5file.keys() or 'mask' in h5file.keys():
      k=h5file.keys()
      print 'multi looking the '+ k[0]
   
      group=h5file_lks.create_group(k[0])    
      dset1 = h5file[k[0]].get(k[0])
      Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
      Masklks=multilook(Mask,alks,rlks)
      dset = group.create_dataset(k[0], data=Masklks, compression='gzip')
      for key , value in h5file[k[0]].attrs.iteritems():
         group.attrs[key]=value

      try:
           group.attrs['Y_STEP']=alks*float(group.attrs['Y_STEP'])
           group.attrs['X_STEP']=rlks*float(group.attrs['X_STEP'])
      except:
           group.attrs['AZIMUTH_PIXEL_SIZE']=alks*float(group.attrs['AZIMUTH_PIXEL_SIZE'])
           group.attrs['RANGE_PIXEL_SIZE']=rlks*float(group.attrs['RANGE_PIXEL_SIZE'])
    group.attrs['WIDTH']=Masklks.shape[1]
    group.attrs['FILE_LENGTH']=Masklks.shape[0]
    h5file.close()
    h5file_lks.close()

if __name__ == '__main__':

  main(sys.argv[1:])



