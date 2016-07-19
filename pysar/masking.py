#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Oct 2015: add support for ROI_PAC product
# Yunjun, Jul 2016: add mask_data(), mask_file()
#                   add parallel processing using joblib


import os
import sys
import getopt
import glob

import h5py
import numpy as np

import pysar._pysar_utilities as ut
import pysar._readfile as readfile


############################################################
def mask_data(data,mask):
  ## Masking a 2D matrxi data with mask
  try:
      xsub
      ysub
      mask[ysub[0]:ysub[1],xsub[0]:xsub[1]]=0
  except:   pass

  try:     data[mask<thr] = np.nan
  except:  data[mask==0]  = np.nan

  return data


############################################################
def mask_with_multi_masks(in_file,mask_file):

  h5file=h5py.File(in_file,'r')
  h5mask=h5py.File(mask_file,'r')
  kf=h5file.keys()

  ext = os.path.splitext(in_file)[1]
  out_file = in_file.split('.')[0]+'_masked'+ext
  h5out = h5py.File(out_file,'w')

  if kf[0] in ('interferograms','wrapped','coherence') and 'coherence' in h5mask.keys():
      print 'file type: '+kf[0]
      print 'Masking each '+kf[0]+' using its coherence file'
      igramList = h5file[kf[0]].keys();         igramList = sorted(igramList)
      cohList   = h5mask['coherence'].keys();   cohList   = sorted(cohList)
      gg = h5out.create_group(kf[0])
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

          unw = mask_data(unw,coh)

          group = gg.create_group(igram)
          dset = group.create_dataset(igram, data=unw, compression='gzip')
          for key, value in h5file[kf[0]][igram].attrs.iteritems():   group.attrs[key] = value
          del unw, coh
      try:
          mask = h5file['mask'].get('mask')
          gm = h5out.create_group('mask')
          dset = gm.create_dataset('mask', data=mask, compression='gzip')
      except: print 'no mask group found.'

  else: print 'Only support multiple dataset file with multiple masks.'; sys.exit(1)


############################################################
def mask_file(in_file,M):
  ## Mask input file with mask matrix M

  atr = readfile.read_attributes(in_file)
  k = atr['FILE_TYPE']
  print 'file type: '+k

  ext      = os.path.splitext(in_file)[1]
  out_file = os.path.basename(in_file).split('.')[0]+'_masked'+ext

  if k in ['timeseries','interferograms','wrapped','coherence']:
      h5file = h5py.File(in_file,'r')
      epochList = h5file[k].keys()
      epochList = sorted(epochList)
      print 'number of epochs: '+str(len(epochList))

      h5out = h5py.File(out_file,'w')
      print 'writing >>> '+out_file

  ##### Multiple Dataset File
  if k == 'timeseries':
      group = h5out.create_group(k)
      for d in epochList:
          print d
          unwset = h5file[k].get(d)
          unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

          unw = mask_data(unw,M)

          dset = group.create_dataset(d, data=unw, compression='gzip')
      for key,value in atr.iteritems():   group.attrs[key] = value

  elif k in ['interferograms','wrapped','coherence']:
      gg = h5out.create_group(k)
      for igram in epochList:
          print igram
          unwset = h5file[kf[0]][igram].get(igram)
          unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

          unw = mask_data(unw,M)

          group = gg.create_group(igram)
          dset = group.create_dataset(igram, data=unw, compression='gzip')
          for key, value in h5file[k][igram].attrs.iteritems():
              group.attrs[key] = value
      try:
          mask = h5file['mask'].get('mask')
          gm = h5out.create_group('mask')
          dset = gm.create_dataset('mask', data=mask, compression='gzip')
      except: print 'no mask group found.'

  ##### Single Dataset File
  else:
      import pysar._writefile as writefile
      unw,atr = readfile.read(in_file)
      unw     = mask_data(unw,M)
      writefile.write(unw,atr,out_file)

  try:
      h5file.close()
      h5out.close()
  except: pass


############################################################
def Usage():
  print '''
**************************************************************************
  Masking File(s) using a mask

  Usage:
      masking.py file MaskFile
      masking.py -f file -m MaskFile [ -t threshold ]

      -f : file (list) that needed to be masked
      -m : mask file 
      -o : output file name

      -t : threshold value used for masking. if not 
           specified then only pixels with mask value 
           equal to zero is masked out.
      -x : masking subset in range   / column / longtitude direction
      -y : masking subset in azimuth / row    / latitude   direction

      --parallel : enable parallel computing

  Example: 
      masking.py velocity Mask.h5

      masking.py -f geo_100102_101120.unw -m Mask.h5
      masking.py -f timeseries.h5         -m temporal_coherence.h5 -t 0.7
      masking.py -f LoadedData.h5         -m 100102_101120.cor     -t 0.9 -y '200:300' -x '300:400'

      masking.py -f 'timeseries*.h5'              -m Mask.h5
      masking.py -f 'timeseries*.h5'              -m Mask.h5 --parallel
      masking.py -f 'timeseries*.h5,*velocity.h5' -m Mask.h5 --parallel

**************************************************************************
  '''

############################################################
def main(argv):

  global xsub, ysub, thr
  #lineWidth=4
  #fontSize=32
  #markerColor='orange'
  #markerSize=20
  #if len(sys.argv)>2:
  parallel = 'no'

  ######################################
  try:    opts, args = getopt.getopt(argv,"h:f:m:t:x:y:o:",['parallel'])
  except getopt.GetoptError:    Usage() ; sys.exit(1)

  if len(sys.argv) > 3:
      for opt,arg in opts:
          if opt in ("-h","--help"):     Usage();  sys.exit()
          elif opt == '-f':        File     = arg.split(',')
          elif opt == '-m':        maskFile = arg
          elif opt == '-t':        thr  = float(arg)
          elif opt == '-y':        ysub = [int(i) for i in arg.split(':')];     ysub.sort()
          elif opt == '-x':        xsub = [int(i) for i in arg.split(':')];     xsub.sort()
          elif opt == '-o':        outFile = arg
          elif opt == '--parallel':   parallel = 'yes'

  elif len(sys.argv)==3:
      if os.path.isfile(argv[0]) and os.path.isfile(argv[1]):
          File     = argv[0].split(',')
          maskFile = argv[1]
      else:  print 'Input file does not existed: '+argv[0]+' / '+argv[1];  sys.exit(1)
  else:   Usage();  sys.exit(1)

  try:
      File
      maskFile
  except:    Usage() ; sys.exit(1)

  ##### Check Input File List
  print '\n****************** Masking *********************'
  fileList = ut.get_file_list(File)
  print 'number of file to mask: '+str(len(fileList))
  print fileList

  ######################################
  atr_mask = readfile.read_attributes(maskFile)
  k_mask = atr_mask['FILE_TYPE']

  ##### Multiple Mask
  if k_mask == 'coherence':
      parallel = 'no'
      mask_with_multi_masks(fileList[0],maskFile)

  ##### Single Mask
  else:
      M,Matr = readfile.read(maskFile)
      print 'mask file: '+maskFile

      if parallel == 'no':
          for file in fileList:
              print '-------------------------------------------'
              print 'masking  : '+file
              mask_file(file,M)

      else:
          print '-----------------------'
          print 'parallel masking ...'
          print '-----------------------'
          from joblib import Parallel, delayed
          import multiprocessing
          num_cores = multiprocessing.cpu_count()
          Parallel(n_jobs=num_cores)(delayed(mask_file)(file,M) for file in fileList)


############################################################
if __name__ == '__main__':
  main(sys.argv[1:])


