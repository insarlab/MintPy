#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Oct 2015: add support for ROI_PAC product

import sys
import os
import getopt

def Usage():
  print '''
  ******************************************
  ******************************************
  Masking a PySAR product using a mask file.

  Usage:
        masking.py -f file -m MaskFile
        masking.py -f file -m MaskFile -t threshold
        masking.py -f file -m MaskFile -o outFile

  file: an hdf5/roi_pac file format to be masked.
  MaskFile: the mask file 
  threshold: the threshold value used for masking. if not 
             specified then only pixels with mask value 
             equal to zero is masked out.

  if user wants to mask out an area with coordinates y1<y<y2 and x1<x<x2:
          masking.py -f file -m MaskFile -t threshold -y 'y1:y2' -x 'x1:x2'
   
  Example: 
    masking velocity:
          masking.py -f velocity.h5 -m temporal_coherence.h5 -t 0.7
          masking.py -f velocity.h5 -m Mask.h5
          masking.py -f velocity.h5 -m temporal_coherence.h5 -t 0.9 -y '200:300' -x '300:400'

    masking interferograms (unwrapped/wrapped): 
          masking.py -f LoadedData.h5 -m temporal_coherence.h5 -t 0.7
          masking.py -f LoadedData.h5 -m Coherence.h5 -t 0.7
          masking.py -f Wrapped.h5    -m temporal_coherence.h5 -t 0.7
          masking.py -f Wrapped.h5    -m Coherence.h5 -t 0.7
          masking.py -f geo_100102_101120.unw -m geo_100102_101120.cor -t 0.7
          masking.py -f geo_100102_101120.unw -m geo_temporal_coherence.h5 -t 0.7

    masking time-series:
          masking.py -f timeseries.h5 -m temporal_coherence.h5 -t 0.7    
  ******************************************
  ******************************************
  '''

######################################
def main(argv):

  #lineWidth=4
  #fontSize=32
  #markerColor='orange'
  #markerSize=20
  #if len(sys.argv)>2:

  try:
    opts, args = getopt.getopt(argv,"h:f:m:t:x:y:o:")
  except getopt.GetoptError:
    Usage() ; sys.exit(1)

  for opt,arg in opts:
    if opt in ("-h","--help"):
      Usage()
      sys.exit()
    elif opt == '-f':        File = arg
    elif opt == '-m':        maskFile = arg
    elif opt == '-t':        thr  = float(arg)
    elif opt == '-y':        ysub = [int(i) for i in arg.split(':')];     ysub.sort()
    elif opt == '-x':        xsub = [int(i) for i in arg.split(':')];     xsub.sort()
    elif opt == '-o':        MaskedFile = arg

  try:
    File
    maskFile
  except:
     Usage() ; sys.exit(1)

######################################

  ext = os.path.splitext(File)[1]
  import h5py
  import numpy as np

  ################################  PySAR HDF5  #################################

  if ext == '.h5':
    h5file=h5py.File(File,'r')
    h5mask=h5py.File(maskFile,'r')
    kf=h5file.keys()

    if 'coherence' not in h5mask.keys():
       Mset=h5mask[h5mask.keys()[0]].get(h5mask.keys()[0])
       M=Mset[0:Mset.shape[0],0:Mset.shape[1]]

    try:     MaskedFile
    except:  MaskedFile = File.split('.')[0]+'_masked.h5'
    h5file2 = h5py.File(MaskedFile,'w')

    if len(kf)==1 and kf[0] in ('velocity','temporal_coherence','rmse','mask'):
      Vset=h5file[h5file.keys()[0]].get(h5file.keys()[0])
      V=Vset[0:Vset.shape[0],0:Vset.shape[1]]
 
      try:
        xsub
        ysub
        M[ysub[0]:ysub[1],xsub[0]:xsub[1]]=0
      except:   print 'No subset'

      try:     V[M<thr] = np.nan
      except:  V[M==0]  = np.nan
    
      group=h5file2.create_group(kf[0])
      dset = group.create_dataset(os.path.basename(kf[0]), data=V, compression='gzip')

      for key, value in h5file[kf[0]].attrs.iteritems():
            group.attrs[key] = value


    elif 'timeseries' in h5file.keys():
      print 'Masking the time-series'
      group = h5file2.create_group('timeseries')
      dateList=h5file['timeseries'].keys()
      for d in dateList:
        print d
        unwset = h5file['timeseries'].get(d)
        unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

        try:     unw[M<thr] = np.nan
        except:  unw[M==0]  = np.nan

        dset = group.create_dataset(d, data=unw, compression='gzip')

      for key,value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value


    elif kf[0] in ('interferograms','wrapped') and 'coherence' in h5mask.keys():
      print 'Masking each '+kf[0]+' using its coherence file'
      igramList=h5file[kf[0]].keys()
      cohList=h5mask['coherence'].keys()
      gg = h5file2.create_group(kf[0])
      for igram in igramList:
         print igram
         date12 = h5file[kf[0]][igram].attrs['DATE12']
         for cohFile in cohList:
             if h5mask['coherence'][cohFile].attrs['DATE12']==date12:
                 igramCoh=cohFile
         print igramCoh

         unwset = h5file[kf[0]][igram].get(igram)
         unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

         cohset=h5mask['coherence'][igramCoh].get(igramCoh)      
         coh=cohset[0:cohset.shape[0],0:cohset.shape[1]]
         unw[coh<thr]=np.nan
         group = gg.create_group(igram)
         dset = group.create_dataset(igram, data=unw, compression='gzip')
         for key, value in h5file[kf[0]][igram].attrs.iteritems():
            group.attrs[key] = value
      try:
      #if kf[0] == 'interferograms':
         mask = h5file['mask'].get('mask')
         gm = h5file2.create_group('mask')
         dset = gm.create_dataset('mask', data=mask, compression='gzip')
      except: print 'no mask group found.'


    elif kf[0] in ('interferograms','wrapped') and 'coherence' not in h5mask.keys():
      print 'Masking the '+kf[0]
      igramList=h5file[kf[0]].keys()     
      gg = h5file2.create_group(kf[0])
      for igram in igramList:
         print igram
         unwset = h5file[kf[0]][igram].get(igram)
         unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]
         try:     unw[M<thr] = np.nan
         except:  unw[M==0]  = np.nan

         group = gg.create_group(igram)
         dset = group.create_dataset(igram, data=unw, compression='gzip')
         for key, value in h5file[kf[0]][igram].attrs.iteritems():
            group.attrs[key] = value

      try:
      #if kf[0] == 'interferograms':
         mask = h5file['mask'].get('mask')
         gm = h5file2.create_group('mask')
         dset = gm.create_dataset('mask', data=mask, compression='gzip') 
      except: print 'no mask group found.'


    h5file.close()
    h5mask.close()
    h5file2.close()    


  ################################  ROI_PAC  ######################################

  elif ext in ['.unw','.int']:
    import pysar._readfile as readfile
    import pysar._writefile as writefile
    Mext = os.path.splitext(maskFile)[1]
    if Mext == '.cor':
       a,M,atr = readfile.read_float32(maskFile)
    elif Mext == '.h5':
       h5mask=h5py.File(maskFile,'r')
       Mset=h5mask[h5mask.keys()[0]].get(h5mask.keys()[0])
       M=Mset[0:Mset.shape[0],0:Mset.shape[1]]
    else: sys.exit(['Un-recoganized mask file extension: '+Mext])

    try:     MaskedFile
    except:  MaskedFile = File.split('.')[0]+'_masked'+ext

    a,unw,atr = readfile.read_float32(File)
    try:     unw[M<thr] = np.nan
    except:  unw[M==0]  = np.nan

    print 'writting >>> '+MaskedFile
    writefile.write_float32(unw,MaskedFile)
    f = open(MaskedFile+'.rsc','w')
    for k in atr.keys():  f.write(k+'    '+atr[k]+'\n')
    f.close()

    try: h5mask.close()
    except: pass

############################################################

if __name__ == '__main__':
  main(sys.argv[1:])

