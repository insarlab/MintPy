#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jun 2016: Add template input option
#                   Add multiple files support
#


import os
import sys
import getopt
import glob

import numpy as np
import h5py

import pysar._remove_surface as rm
import pysar._readfile as readfile


######################################
def Usage():
  print '''
********************************************************
********************************************************

    Remove phase ramp

    Usage:
        remove_plane.py file method [Maskfile]

        -f : input file (list) that need to remove ramp
        -s : quadratic, plane, quardatic_range, quadratic_azimiuth, plane_range, plane_azimuth
        -m : (optional) a mask file with 0 values for those pixels which are not considered in
               plane estimation.
        -t : template file
        -o : output name

    example:
        remove_plane.py  timeseries.h5 plane
        remove_plane.py  timeseries.h5 plane             Mask.h5
        remove_plane.py  LoadedData.h5 quadratic_range   Mask.h5

        remove_plane.py  -f timeseries.h5 -t KyushuT424F640AlosA.template

        remove_plane.py  -f 'geo_100102_*.unw'  -s plane -m Mask_tempCoh.h5

********************************************************
********************************************************
  '''

######################################
def main(argv):

  ########################## Check Inputs ################################################
  ## Default value
  Masking  = 'no'
  surfType = 'plane'

  if len(sys.argv) > 4:
      try: opts, args = getopt.getopt(argv,'h:f:m:o:s:t:')
      except getopt.GetoptError:  print 'Error while getting args!\n';  Usage(); sys.exit(1)

      for opt,arg in opts:
          if   opt in ['-h','--help']:    Usage(); sys.exit()
          elif opt in '-f':    File     = arg
          elif opt in '-m':    maskFile = arg
          elif opt in '-o':    outName  = arg
          elif opt in '-s':    surfType = arg.lower()
          elif opt in '-t':    templateFile = arg

  elif len(sys.argv) in [3,4]:
      File          = argv[0]
      surfType      = argv[1].lower()
      try: maskFile = argv[2]
      except: pass
  else: Usage(); sys.exit(1)

  ##### Tempate File
  try:
      templateFile
      templateContents = readfile.read_template(templateFile)
  except: pass

  try: surfType
  except:
      try: surfType = templateContents['pysar.orbitError.method']
      except: print 'No ramp type found!'; sys.exit(1)

  ##### Read Mask File 
  ## Priority:
  ## Input mask file > pysar.mask.file > existed Modified_Mask.h5 > existed Mask.h5
  try:      maskFile
  except:
      try:  maskFile = templateContents['pysar.mask.file']
      except: pass
          #if   os.path.isfile('Modified_Mask.h5'):  maskFile = 'Modified_Mask.h5'
          #elif os.path.isfile('Mask.h5'):           maskFile = 'Mask.h5'
          #else: print 'No mask found!'; sys.exit(1)
  print '\n*************** Phase Ramp Removal ***********************'
  try:
      Mask,Matr = readfile.read(maskFile)
      print 'mask: '+maskFile
      Masking = 'yes'
  except:
      print 'No mask'
      Masking = 'no'
      pass
      #sys.exit(1)

  if Masking=='no':
      atr = readfile.read_attributes(File)
      length = int(atr['FILE_LENGTH'])
      width  = int(atr['WIDTH'])
      Mask=np.ones((length,width))


  ############################## Removing Phase Ramp #######################################
  fileList = glob.glob(File)
  fileList = sorted(fileList)
  print 'number of file to de-ramp: '+str(len(fileList))
  print fileList
  for file in fileList:
      print '------------------------------------------'
      print 'input files : '+file
      rm.remove_surface(file,Mask,surfType)

###########################################################################################
if __name__ == '__main__':
  main(sys.argv[1:])


