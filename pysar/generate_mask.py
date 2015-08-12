#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib import colors
import getopt
import h5py

def Usage():
  print '''
  ******************************************
  ******************************************
  Generating a mask file with the same size of the input file

  Usage:

     generate_mask.py -f file.h5 -m min -M max -y 'ymin:ymax' -x 'xmin:xmax' -o output name

     -y 'ymin:ymax' -x 'xmin:xmax'  generates a mask file bounded by this box.

  Example:
     generate_mask.py -f temporal_coherence.h5 -m 0.3 -M 0.8 -y -o maskfile.h5
     generate_mask.py -f temporal_coherence.h5 -m 0.3 -M 0.8 -y '100:700' -x '200:800'
     generate_mask.py -f temporal_coherence.h5  -y '100:700' -x '200:800'
     generate_mask.py -f velocity.h5 -y '100:700' -x '200:800'
     generate_mask.py -f Mask.h5 -y '100:700' -x '200:800'
   

  ******************************************
  ******************************************
  '''

def main(argv):

  try:
      opts, args = getopt.getopt(argv,"h:f:m:M:x:y:o:")

  except getopt.GetoptError:
      Usage() ; sys.exit(1)

  for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        File = arg
      elif opt == '-m':
        minV = float(arg)
      elif opt == '-M':
        maxV=float(arg)
      elif opt=='-y':
        ysub=[int(i) for i in arg.split(':')]
        ysub.sort()
      elif opt=='-x':
        xsub = [int(i) for i in arg.split(':')]
        xsub.sort()
      elif opt=='-o':
        outName=arg


  try:
     
     h5file=h5py.File(File,'r') 
      
  except:
     Usage() ; sys.exit(1)

  try:  
    h5mask=h5py.File(outName,'w')
  except:
    h5mask=h5py.File('mask.h5','w')

  kf=h5file.keys()

  Vset=h5file[h5file.keys()[0]].get(h5file.keys()[0])
  V=Vset[0:Vset.shape[0],0:Vset.shape[1]]
  
  if 'mask' in h5file.keys():
         mask=V
  else:
         mask=np.ones([Vset.shape[0],Vset.shape[1]])

  try:
      mask[0:ysub[0],:]=0
      mask[ysub[1]:mask.shape[0],:]=0
      mask[:,0:xsub[0]]=0
      mask[:,xsub[1]:mask.shape[1]]=0
     
  except:
      print 'No subset'

  try:
      mask[V<minV]=0
      mask[V>maxV]=0
  except:
      print 'No threshold'  

  group=h5mask.create_group('mask')
  dset = group.create_dataset('mask', data=mask, compression='gzip')

  for key, value in h5file[kf[0]].attrs.iteritems():
          group.attrs[key] = value
 
  h5file.close()
  h5mask.close()  

############################################################

if __name__ == '__main__':
  main(sys.argv[1:])

