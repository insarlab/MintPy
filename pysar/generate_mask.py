#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jan 2016: support ROI_PAC files
# Yunjun, Jun 2016: use readfile.read()
#                   Add nonzero method, equivalent to Mask.h5


import sys
import os
import getopt

import numpy as np
import h5py

import pysar._readfile as readfile


def Usage():
  print '''
  ******************************************
  ******************************************
  Generating a mask file with the same size of the input file

  Usage:

     generate_mask.py -f file [ -m min -M max -y ymin:ymax -x xmin:xmax -o output_file ]

     -f: file used to generate mask, supported file format:
             PySAR HDF5 files: temporal_coherence, velocity, dem, mask, rmse
             ROI_PAC    files: .unw .cor .hgt .dem
     -d: epoch date   for timeseries/interferograms file
     -e: epoch number for timeseries/interferograms file (start from 1)
     -y: bounding box in y direction
     -x: boudning box in x direction
     -m: minimum value
     -M: maximum value
     -o: output file, in HDF5 format [default name is mask.h5]

     --nonzero : mask of all nonzero pixels, equivalent to Mask.h5 from LoadedData.h5


  Example:
     generate_mask.py -f temporal_coherence.h5 -m 0.3 -M 0.8 -o mask_temp_coh.h5
     generate_mask.py -f temporal_coherence.h5 -m 0.3 -M 0.8 -y 100:700 -x 200:800
     generate_mask.py -f temporal_coherence.h5 -y 100:700 -x 200:800
     generate_mask.py -f Mask.h5               -y 100:700 -x 200:800
     generate_mask.py -f 081018_090118.unw     -m 2    -M 4
     generate_mask.py -f 081018_090118.cor     -m 0.5
     generate_mask.py -f srtm1.dem             -m 1000 -M 1500
     generate_mask.py -f Seed_LoadedData_mask.h5 -e 50 -m 4

     generate_mask.py -f LoadedData.h5 --nonzero

  ******************************************
  ******************************************
  '''

def main(argv):

  outName = 'mask.h5'
  method  = 'threshold'

  ##### Check Inputs
  if len(sys.argv)>2:
      try:   opts, args = getopt.getopt(argv,'h:f:m:M:x:y:o:d:e:',['nonzero'])
      except getopt.GetoptError:      Usage() ; sys.exit(1)

      for opt,arg in opts:
          if opt in ("-h","--help"):   Usage();   sys.exit()
          elif opt == '-f':         File = arg
          elif opt == '-m':         minV = float(arg)
          elif opt == '-M':         maxV = float(arg)
          elif opt == '-y':         ysub = [int(i) for i in arg.split(':')];        ysub.sort()
          elif opt == '-x':         xsub = [int(i) for i in arg.split(':')];        xsub.sort()
          elif opt == '-o':         outName    = arg
          elif opt == '-d':         epoch_date = arg
          elif opt == '-e':         epoch_num  = int(arg) - 1
          elif opt == '--nonzero':  method     = 'nonzero'

  elif len(sys.argv)==2:
      if   argv[0] in ['-h','--help']:    Usage(); sys.exit(1)
      elif os.path.isfile(argv[0]):       File = argv[0]
      else:    print 'Input file does not existed: '+argv[0];  sys.exit(1)
  else:                                   Usage(); sys.exit(1)


  ##### Input File Info
  atr = readfile.read_attributes(File)
  print '\n****************** Generate Mask *******************'
  print 'Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']+': '+File
  mask = np.ones([int(atr['FILE_LENGTH']),int(atr['WIDTH'])])
  print 'Create initial mask with the same size as the input file and all = 1'


  ##### Non-zero Mask #######
  if method == 'nonzero':
      k = atr['FILE_TYPE']
      MaskZero = np.ones([int(atr['FILE_LENGTH']),int(atr['WIDTH'])])

      ext = os.path.splitext(File)[1].lower()
      if ext == '.h5' and k in ['interferograms','coherence','wrapped','timeseries']:
          h5file = h5py.File(File,'r')
          epochList = h5file[k].keys()

          for epoch in epochList:
              print epoch
              if k in ['interferograms','coherence','wrapped']:
                  data = h5file[k][epoch].get(epoch)[:]
              elif k in ['timeseries']:
                  data = h5file[k].get(epoch)
              MaskZero *= data
              MaskZero[np.isnan(data)] = 0
          h5file.close()

      else:
          data,atr = readfile.read(File)
          MaskZero *= data
          MaskZero[np.isnan(data)] = 0

      mask = np.ones([int(atr['FILE_LENGTH']),int(atr['WIDTH'])])
      mask[MaskZero==0] = 0


  ##### Threshold ##########
  else:
      ##### Read and Initiate Mask
      try:        V, atr = readfile.read(File,epoch_date)
      except:
          try:    V, atr = readfile.read(File,epoch_num)
          except: V, atr = readfile.read(File)

      ##### Calculating Mask
      ## threshold
      try:
          mask[V<minV]=0
          print 'all value < '+str(minV)+' = 0'
      except:  print 'No min threshold'
      try:
          mask[V>maxV]=0
          print 'all value > '+str(maxV)+' = 0'
      except:  print 'No max threshold'  
      ## nan value
      mask[np.isnan(V)]=0

  ## subset
  try:
      mask[0:ysub[0],:]=0
      mask[ysub[1]:mask.shape[0],:]=0
      print 'all y in [0,'+str(ysub[0])+'] and ['+str(ysub[1])+',end] = 0'
  except:  print 'No subset in y direction'
  try:
      mask[:,0:xsub[0]]=0
      mask[:,xsub[1]:mask.shape[1]]=0
      print 'all x in [0,'+str(xsub[0])+'] and ['+str(xsub[1])+',end] = 0'
  except:  print 'No subset in x direction'
 

  ##### Writing mask file
  h5mask = h5py.File(outName,'w');  print 'writing >>> '+outName
  group = h5mask.create_group('mask')
  dset = group.create_dataset('mask', data=mask, compression='gzip')
  for key, value in atr.iteritems():    group.attrs[key] = value
  h5mask.close()  


############################################################
if __name__ == '__main__':
  main(sys.argv[1:])

