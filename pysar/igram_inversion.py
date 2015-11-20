#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import getopt
#import re
#import time
#import datetime
import h5py
import _pysar_utilities as ut

######################################
def Usage():
  print '********************************************************************************\n'
  print 'Inversion of interferograms using L1 or L2 norm minimization (Default is L2)'
  print 'Usage:'
  print '          igram_inversion.py interferogramsFile.h5  method'
  print 'Example:'
  print '          igram_inversion.py interferograms.h5'
  print '          igram_inversion.py -f interferograms.h5 -l L1'
  print '          igram_inversion.py -f interferograms.h5 -m average_spatial_coherence.h5'
  print '\n********************************************************************************'


######################################
def main(argv):

  inversion_method = 'L2'

  try:
    opts, args = getopt.getopt(argv,"h:f:m:l:")
  except getopt.GetoptError:
    Usage() ; sys.exit(1)

  if len(sys.argv)>2:
    for opt,arg in opts:
      if opt in ("-h","--help"):  Usage();   sys.exit()
      elif opt == '-f':           igramsFile        = arg
      elif opt == '-m':           maskFile          = arg
      elif opt == '-l':           inversion_method  = arg

  elif len(sys.argv)==2:
    if os.path.isfile(argv[0]):   igramsFile = argv[0]
    else:  Usage(); sys.exit(1)
  elif len(sys.argv)<2:           Usage(); sys.exit(1)

  #try:     igramsFile = argv[0]
  #except:  Usage() ; sys.exit(1)
  #try:     inversion_method = argv[1]
  #except:  inversion_method = 'L2'

  #try:
  h5file = h5py.File(igramsFile,'r')
  if not 'interferograms' in h5file.keys():
      print '**********************************************************************'
      print 'ERROR:'
      print '       '+igramsFile+ '  was not found or the file is not readable!'
      print '**********************************************************************'
      Usage();sys.exit(1)

  numIfgrams = len(h5file['interferograms'].keys())
  if numIfgrams == 0.:
      print "load interferograms first by running: load_data.py TEMPLATEFILE"
      sys.exit(1)
  h5timeseries = h5py.File('timeseries.h5','w')
  if   inversion_method in ['L2','l2']: 
      ut.timeseries_inversion(h5file,h5timeseries)
  elif inversion_method in ['L1','l1']:
      ut.timeseries_inversion_L1(h5file,h5timeseries)
  else:
      print '******************************'
      print 'WARNING:'
      print '  Inversion method not recognized!'
      print '  L2 norm minimization is used.'
      print '******************************'
      ut.timeseries_inversion(h5file,h5timeseries)

  ## generate 'mask' for timeseries.h5
  try:
      maskFile
      h5filem = h5py.File(maskFile,'r')
      dset1 = h5filem['mask'].get('mask')
      Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
      group=h5timeseries.create_group('mask')
      dset = group.create_dataset('mask', data=Mask, compression='gzip')
      h5filem.close()
  except: 
      dset1 = h5file['mask'].get('mask')
      Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]] 
      group=h5timeseries.create_group('mask')
      dset = group.create_dataset('mask', data=Mask, compression='gzip')

  h5file.close()
  h5timeseries.close()
  #h5flat.close()

######################################
if __name__ == '__main__':

  main(sys.argv[1:])

