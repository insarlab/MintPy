#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Reference:
# Fattahi, H. and F. Amelung, (2013), DEM-error correction in 
# InSAR time-series analysis, IEEE TGRS, vol. no.99,
# doi: 10.1109/TGRS.2012.2227761.
#
# Yunjun, Jun 2016: Add phase velocity approach from the paper.
#                   Use different range and look angle for each column
#


import sys
import os
import getopt
import time
import datetime

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._pysar_utilities as ut
#import pysar._writefile as writefile
import pysar._readfile as readfile


######################################
def Usage():
  print '''
  **********************************************************

  DEM Error Correction
      Reference: Fattahi and Amelung, 2013, IEEE TGRS.

  Usage:

      dem_error.py timeSeriesFile interferogramsFile [ outputFile ]

      -f : timeseries file
      -F : interferograms file for baseline time series calculation
      -o : output filename [add _demCor by default]

      --phase-velocity: use phase velocity history approach instead of
                            phase history approach.
      --no-timeseries-update: calculate DEM error file only and do not
                              update the timeseries file.


  Example:
      dem_error.py timeseries_ECMWF.h5
      dem_error.py -f timeseries_ECMWF.h5 --phase-velocity
      dem_error.py -f timeseries_NARR.h5  -o timeseries_demCor.h5
      dem_error.py -f timeseries_MERRA.h5 -F Seeded_LoadedData.h5

  **********************************************************
  '''

######################################
def main(argv):

  ## Default value
  phase_velocity    = 'no'        # 'no' means use 'phase history'
  update_timeseries = 'yes'

  if len(sys.argv)>2:
      try:   opts, args = getopt.getopt(argv,'h:f:F:o:v:',['phase-velocity','no-timeseries-update'])
      except getopt.GetoptError:
          print 'Error in reading input options!';  Usage() ; sys.exit(1)

      for opt,arg in opts:
          if opt in ['-h','--help']:    Usage() ; sys.exit()
          elif opt == '-f':    timeSeriesFile = arg
          elif opt == '-F':    igramsFile     = arg
          elif opt == '-o':    outname        = arg
          elif opt == '--phase-velocity'      :  phase_velocity = 'yes'
          elif opt == '--no-timeseries-update':  update_timeseries = 'no'

  elif len(sys.argv)==2:
      if argv[0] in ['-h','--help']:  Usage(); sys.exit(1)
      else:  timeSeriesFile = argv[0]
  else:  Usage(); sys.exit(1)

  try:    outname
  except: outname = timeSeriesFile.replace('.h5','')+'_demCor.h5'

  ##### Read Time Series
  print '\n*************** Topographic Error Correction ****************'
  print "Loading time series: " + timeSeriesFile
  atr = readfile.read_attributes(timeSeriesFile)
  h5timeseries = h5py.File(timeSeriesFile)
  dateList = h5timeseries['timeseries'].keys()
  dateList.sort()
  lt = len(dateList)
  print 'number of epochs: '+str(lt)

  dateIndex={}
  for ni in range(len(dateList)):   dateIndex[dateList[ni]]=ni

  nrows = int(atr['FILE_LENGTH'])
  ncols = int(atr['WIDTH'])
  timeseries = np.zeros((len(dateList),nrows*ncols),np.float32)
  for date in dateList:
      print date
      d = h5timeseries['timeseries'].get(date)[:]
      timeseries[dateIndex[date]][:]=d.flatten('F')
  del d
  h5timeseries.close()   
  print '**************************************'

  ##### Temporal Baseline
  print 'read temporal baseline'
  tbase,date_dict = ptime.date_list2tbase(dateList)
  tbase = np.array(tbase).reshape(lt,1)

  ##### Perpendicular Baseline
  try:
      Bp = [float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]
      Bp = np.array(Bp).reshape(lt,1)
  except:
      print 'Cannot find P_BASELINE_TIMESERIES from timeseries file.'
      print 'Trying to calculate it from interferograms file'
      try:
          Bp = ut.Baseline_timeseries(igramsFile)
          Bp = np.array(Bp).reshape(lt,1)
      except:
          print 'Error in calculating baseline time series!'
          sys.exit(1)
  Bp_v = (Bp[1:lt] - Bp[0:lt-1]) / (tbase[1:lt] - tbase[0:lt-1])


  ##### Cubic Temporal Deformation Model
  ## Formula (10) in (Fattahi and Amelung, 2013, TGRS)
  if phase_velocity == 'yes':
      print 'using phase velocity history'
      M1 = np.ones((lt-1,1))
      M2 = (tbase[1:lt]+tbase[0:lt-1])/2
      M3 = (tbase[1:lt]**2 + tbase[1:lt]*tbase[0:lt-1] + tbase[0:lt-1]**2)/6
      M  = np.hstack((M1,M2,M3))
  else:
      print 'using phase history'
      M  = np.hstack((.5*tbase**2,tbase,np.ones((lt,1))))

  ## Testing
  #teta = (tetaN+tetaF)/2
  #r = (rN+rF)/2
  #teta=19.658799999999999*np.pi/180
  #r=846848.2
  #Bperp=1000*np.random.random((lt,1))

  ##### Range and Look Angle
  near_range = float(atr['STARTING_RANGE1'])
  dR         = float(atr['RANGE_PIXEL_SIZE'])
  r          = float(atr['EARTH_RADIUS'])
  H          = float(atr['HEIGHT'])
  far_range  = near_range + dR*(ncols-1)
  incidence_n = np.pi-np.arccos((r**2+near_range**2-(r+H)**2)/(2*r*near_range))
  incidence_f = np.pi-np.arccos((r**2+ far_range**2-(r+H)**2)/(2*r*far_range))

  various_range = 'yes'
  if various_range == 'yes':
      range_x      = np.linspace(near_range, far_range,  num=ncols,endpoint='FALSE')
      look_angle_x = np.linspace(incidence_n,incidence_f,num=ncols,endpoint='FALSE')
  else:
      print 'using center range and look angle to represent the whole area'
      center_range = (near_range + far_range)/2
      center_look_angle = np.pi-np.arccos((r**2+center_range**2-(r+H)**2)/(2*r*center_range))
      range_x      = np.tile(center_range,     ncols)
      look_angle_x = np.tile(center_look_angle,ncols)
  #C1_v = Bp_v / (center_range * np.sin(center_look_angle))
  #C1   = Bp   / (center_range * np.sin(center_look_angle))
  #timeseries_v = (timeseries[1:lt,:] - timeseries[0:lt-1,:]) / (tbase[1:lt] - tbase[0:lt-1])

  ##### Inversion column by column
  print 'inversing using L2-norm minimization (unweighted least squares)...'
  dz = np.zeros([1,nrows*ncols])

  for i in range(ncols):
      ## Design Matrix Inversion
      C1_v = Bp_v / (range_x[i] * np.sin(look_angle_x[i]))
      C1   = Bp   / (range_x[i] * np.sin(look_angle_x[i]))
      if phase_velocity == 'yes':  C = np.hstack((M,C1_v))
      else:                        C = np.hstack((M,C1))

      #print '    rank of the design matrix : '+str(np.linalg.matrix_rank(C))
      #if np.linalg.matrix_rank(C) == 4:  print '    design matrix has full rank'
      Cinv = np.linalg.pinv(C)

      ## (Phase) Velocity History
      ts_x  = timeseries[:,i*nrows:(i+1)*nrows]
      ts_xv = (ts_x[1:lt,:] - ts_x[0:lt-1,:]) / (tbase[1:lt] - tbase[0:lt-1])

      ## DEM error
      if phase_velocity == 'yes':    par  = np.dot(Cinv,ts_xv)
      else:                          par  = np.dot(Cinv,ts_x)
      dz_x = par[3].reshape((1,nrows))

      ## Update DEM error matrix and timeseries matrix
      dz[0][i*nrows:(i+1)*nrows]         = dz_x
      timeseries[:,i*nrows:(i+1)*nrows] -= np.dot(C1,dz_x)

      ut.printProgress(i+1,ncols)

  #dz[0][:] = par[3][:]
  dz = np.reshape(dz,[nrows,ncols],order='F')


  ########## Output - DEM error #######################
  #print '**************************************'
  #print 'writing DEM_error.hgt'
  #writefile.write_float32(dz,'DEM_error.hgt')
  #f = open('DEM_error.hgt.rsc','w')
  #f.write('FILE_LENGTH       '+str(int(nrows))+'\n')
  #f.write('WIDTH             '+str(int(ncols))+'\n')  
  #print '**************************************'

  h5fileDEM = 'DEM_error.h5'
  print 'writing >>> '+h5fileDEM
  h5rmse = h5py.File(h5fileDEM,'w')
  group=h5rmse.create_group('dem')
  dset = group.create_dataset(os.path.basename('dem'), data=dz, compression='gzip')
  for key , value in atr.iteritems():
      group.attrs[key]=value
  group.attrs['UNIT']='m'
  print '**************************************'

  ########### Output - Corrected Time Series ##########
  if update_timeseries == 'yes':
      print 'writing >>> '+outname
      print 'number of dates: '+str(len(dateList))

      h5timeseriesDEMcor = h5py.File(outname,'w')
      group = h5timeseriesDEMcor.create_group('timeseries')
      for i in range(len(dateList)):
          print dateList[i]
          d = np.reshape(timeseries[i][:],[nrows,ncols],order='F')
          dset = group.create_dataset(dateList[i], data=d, compression='gzip')
      #for date in dateList:
      #    print date
      #    if not date in h5timeseriesDEMcor['timeseries']:
      #        d = np.reshape(timeseries[dateIndex[date]][:],[nrows,ncols],order='F')
      #        dset = group.create_dataset(date, data=d, compression='gzip') 
      for key,value in atr.iteritems():  group.attrs[key] = value
      h5timeseriesDEMcor.close()


################################################################################
if __name__ == '__main__':
  main(sys.argv[1:])  



