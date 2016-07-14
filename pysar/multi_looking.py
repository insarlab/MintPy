#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Oct 2015: Merge timeseries/velocity into one
#                   Merge all non-hdf5 into one

import sys
import os
import getopt
import numpy as np

##################################################################################################

def Usage():
   print '''
   ***************************************************************

   Usage:

     multi_looking.py  file  azimuth_looks  [range_looks]

     file: PySAR h5 files [interferogram, time-series, velocity] 
           roi_pac  files [dem, unw, cor, hgt]
           Image    files [jpeg, jpg, png]

   Example:
       
     multi_looking.py  velocity.h5         15 15
     multi_looking.py  velocity.h5         15
     multi_looking.py  LoadedData.h5       15 20
     multi_looking.py  timeseries.h5       15 15
     multi_looking.py  filt*.unw           10 10
     multi_looking.py  SanAndreas.dem      10 10
     multi_looking.py  SanAndreas.dem.jpeg 10 10

   ***************************************************************
'''

##################################################################################################

def multilook(ifg,lksy,lksx):
    rows,cols=ifg.shape
    lksx = int(lksx)
    lksy = int(lksy)
    rows_lowres=int(np.floor(rows/lksy))
    cols_lowres=int(np.floor(cols/lksx))
    #thr = np.floor(lksx*lksy/2)
    ifg_Clowres=np.zeros((rows,cols_lowres))
    ifg_lowres =np.zeros((rows_lowres,cols_lowres))
    
    for c in range(int(cols_lowres)):  ifg_Clowres[:,c]=np.sum(ifg[:,(c)*lksx:(c+1)*lksx],1)
    for r in range(int(rows_lowres)):  ifg_lowres[r,:] =np.sum(ifg_Clowres[(r)*lksy:(r+1)*lksy,:],0)
    #for c in range(lksx):   ifg_Clowres = ifg_Clowres +         ifg[:,range(c,cols_lowres*lksx,lksx)]
    #for r in range(lksy):   ifg_lowres  = ifg_lowres  + ifg_Clowres[  range(r,rows_lowres*lksy,lksy),:]
    ifg_lowres=ifg_lowres/(lksy*lksx) 

    return ifg_lowres


##################################################################################################

def main(argv):

  try:  
    file=argv[0]
    alks=float(argv[1])
    rlks=float(argv[2])
  except:
    try:
       file=argv[0]
       alks=float(argv[1])
       rlks=float(argv[1])
    except:  Usage();sys.exit(1)

  ext = os.path.splitext(file)[1]
  outName=file.split('.')[0]+'_a'+str(int(alks))+'lks_r'+str(int(rlks))+'lks'+ext

  ################################################################################

  if ext == '.h5':
    import h5py
    h5file_mli=h5py.File(outName,'w')
    h5file=h5py.File(file,'r')
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    elif 'timeseries'   in k: k[0] = 'timeseries'

    if k[0] in ['interferograms','coherence','wrapped']:
      print 'Multilooking the interferograms'
      gg = h5file_mli.create_group(k[0])
      igramList=h5file[k[0]].keys()
      for igram in igramList:
        print igram
        unw = h5file[k[0]][igram].get(igram)
        unwlks=multilook(unw,alks,rlks)
        group = gg.create_group(igram)
        dset = group.create_dataset(igram, data=unwlks, compression='gzip')
        for key, value in h5file[k[0]][igram].attrs.iteritems():
           group.attrs[key] = value
        group.attrs['WIDTH']       = unwlks.shape[1]
        group.attrs['FILE_LENGTH'] = unwlks.shape[0]
        try:
           group.attrs['Y_STEP']=alks*float(group.attrs['Y_STEP'])
           group.attrs['X_STEP']=rlks*float(group.attrs['X_STEP'])
        except:
           group.attrs['AZIMUTH_PIXEL_SIZE'] = alks*float(group.attrs['AZIMUTH_PIXEL_SIZE'])
           group.attrs['RANGE_PIXEL_SIZE']   = rlks*float(group.attrs['RANGE_PIXEL_SIZE'])

      try:
        dset1=h5file['mask'].get('mask')
        mask=dset1[0:dset1.shape[0],0:dset1.shape[1]]
        masklks=multilook(mask,alks,rlks)
        group=h5file_mli.create_group('mask')
        dset = group.create_dataset('mask', data=masklks, compression='gzip')
      except: print 'No mask group found.'

    elif k[0] in ['timeseries','temporal_coherence', 'velocity', 'mask', 'rmse']:
      print 'Multilooking '+k[0]
      group = h5file_mli.create_group(k[0])

      if k[0] == 'timeseries':
        dateList=h5file[k[0]].keys()
        for d in dateList:
          print d
          unw = h5file[k[0]].get(d)
          unwlks=multilook(unw,alks,rlks)
          dset = group.create_dataset(d, data=unwlks, compression='gzip')
      elif k[0] in ['temporal_coherence', 'velocity', 'mask', 'rmse']:
        dset1 = h5file[k[0]].get(k[0])
        unw = dset1[0:dset1.shape[0],0:dset1.shape[1]]
        unwlks=multilook(unw,alks,rlks)
        dset = group.create_dataset(k[0], data=unwlks, compression='gzip')

      try:
        dset1 = h5file['mask'].get('mask')
        Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
        Masklks=multilook(Mask,alks,rlks)
        group=h5file_mli.create_group('mask')
        dset = group.create_dataset('mask', data=Masklks, compression='gzip')
      except:  print 'No mask group found.'

      ## Update attributes
      for key,value in h5file[k[0]].attrs.iteritems():      group.attrs[key] = value
      group.attrs['WIDTH']       = unwlks.shape[1]
      group.attrs['FILE_LENGTH'] = unwlks.shape[0]

      try:
        group.attrs['Y_STEP']=alks*float(group.attrs['Y_STEP'])
        group.attrs['X_STEP']=rlks*float(group.attrs['X_STEP'])
      except:
        group.attrs['AZIMUTH_PIXEL_SIZE'] = alks*float(group.attrs['AZIMUTH_PIXEL_SIZE'])
        group.attrs['RANGE_PIXEL_SIZE']   = rlks*float(group.attrs['RANGE_PIXEL_SIZE'])


    h5file.close()
    h5file_mli.close()

  ################################################################################

  elif ext in ['.unw','.cor','.hgt','.dem','.trans'] or\
       ext in ['.jpeg','.jpg','.png','.ras','.bmp'] or\
       ext in ['.mli','.slc']:
    import pysar._readfile  as readfile 
    import pysar._writefile as writefile
    r = readfile.read_rsc_file(file + '.rsc')

    if ext == '.int' or ext == '.slc':
       a,p,r = readfile.read_complex64(file)
       pmli=multilook(p,alks,rlks)
       amli=multilook(a,alks,rlks)
       r['FILE_LENGTH']=str(pmli.shape[0])
       r['WIDTH']      =str(pmli.shape[1])
    elif ext == '.unw' or ext == '.cor' or ext == '.hgt':
       a,p,r = readfile.read_float32(file)
       pmli=multilook(p,alks,rlks)
       amli=multilook(a,alks,rlks)
       print 'writing >>>'+outName 
       writefile.write_float32(pmli,outName)
       r['FILE_LENGTH']=str(pmli.shape[0])
       r['WIDTH']      =str(pmli.shape[1])
    elif ext == ('.dem'):
       d,r = readfile.read_dem(file)
       dmli=multilook(d,alks,rlks)
       print 'writing >>>'+outName
       writefile.write_dem(dmli,outName)
       r['FILE_LENGTH']=str(dmli.shape[0])
       r['WIDTH']      =str(dmli.shape[1])
    elif ext in ['.jpeg','.jpg','.png','.bmp']:
       import Image
       im = Image.open(file)
       width  = im.size[0] / int(rlks)
       height = im.size[1] / int(alks)
       imlks = im.resize((width, height), Image.NEAREST)
       print 'writing >>>'+outName
       imlks.save(outName)
       r['FILE_LENGTH']=str(height)
       r['WIDTH']      =str(width)

    ## Update attributes
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    try:
       r['Y_STEP']=str(float(r['Y_STEP'])*alks)
       r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:  
       r['AZIMUTH_PIXEL_SIZE'] = alks*float(r['AZIMUTH_PIXEL_SIZE'])
       r['RANGE_PIXEL_SIZE']   = rlks*float(r['RANGE_PIXEL_SIZE'])

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()



###################################################################################################

if __name__ == '__main__':

  main(sys.argv[1:])



