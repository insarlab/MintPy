#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################

import os
import sys
import glob
import time

import numpy as np
import h5py

import _readfile as readfile
import matplotlib.pyplot as plt

def main(argv):
  try:
    templateFile = argv[1]
    print templateFile
  except:
    print '''
    *******************************************

       Copy DEM (radar and geo coordinate) and 
           geomap*.trans file from PROCESS folder
           to TSSAR folder
       
       Usage: copy_dem_trans.py TEMPLATEFILE  

    *******************************************         
    '''
    sys.exit(1)

  templateContents = readfile.read_template(templateFile)
  projectName = os.path.basename(templateFile.partition('.')[0])
  masterDate = templateContents['masterDate']
  Rlooks = templateContents['Rlooks_unw']

  try: 
     processDir = os.getenv('PROCESSDIR') +'/'+ projectName
     tssarDir = os.getenv('TSSARDIR') +'/'+projectName
  except:
     processDir = os.getenv('SCRATCHDIR') + '/' + projectName + '/PROCESS'                           # FA 7/2015: adopted for new directory structure
     tssarDir   = os.getenv('SCRATCHDIR') + '/' + projectName + '/TSSAR'                             # FA 7/2015: adopted for new directory structure
     #tssarDir  = scratchDirectory + '/' + projectName + "/TSSAR"

  print "QQ: AFTER: " + processDir
  print "QQ: AFTER: " + tssarDir
  sys.exit(1)


  """copy geomap_*rlks.trans"""
  geomapPath = processDir+'/GEO/geo_'+masterDate+'*'
  geomapList = glob.glob(geomapPath)
  numGeomap = len(geomapList)
  if numGeomap>0:
    geomapDir = geomapList[0]
    geomapFile = geomapDir+'/geomap_'+Rlooks+'rlks.trans'
    geomapRscFile = geomapFile+'.rsc'
#    print geomapFile
#    print geomapRscFile
    geomapCmd = 'cp '+geomapFile+' '+geomapRscFile+' '+tssarDir
    print geomapCmd
    os.system(geomapCmd)
  else:
    print 'WARNING: no geomap file of master interferograms found!'
    print 'Check the PROCESS/Project/GEO folder.'
    sys.exit(1)

  """copy radar_*rlks.hgt"""
  doneIgramPath = processDir+'/DONE/IFGRAM_*'
  doneIgramList = glob.glob(doneIgramPath)
  numDoneIgram = len(doneIgramList)
  if numDoneIgram>0:
    doneIgramDir = doneIgramList[0]
    radarHgt = doneIgramDir+'/radar_'+Rlooks+'rlks.hgt'
    radarHgtRsc = radarHgt + '.rsc'
#    print radarHgt
#    print radarHgtRsc
    radarCmd = 'cp '+radarHgt+' '+radarHgtRsc+' '+tssarDir
    print radarCmd
    os.system(radarCmd)
  else:
    print 'WARNING: no radar_*rlks.hgt file found!'
    print 'Check the PROCESS/Project/DONE folder.'
    sys.exit(1)


  """copy DEM file from $DEMDIR"""
  demPath=templateContents['DEM']
  demRscPath=demPath+'.rsc'
  demCmd = 'cp '+demPath+' '+demRscPath+' '+tssarDir
  print demCmd
  os.system(demCmd)


if __name__ == '__main__':
  main(sys.argv[:])
