#! /usr/bin/env python

############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import sys
import os
import getopt
import time
import datetime
import h5py
import numpy as np
import matplotlib.pyplot as plt
import _pysar_utilities as ut
import _writefile as writefile

######################################
def Usage():
  print '''
    Usage:

      dem_error.py timeSeriesFile interferogramsFile

    Example:

      dem_error.py timeseries.h5 Loaded_igrams.h5

'''

######################################
def main(argv):
  try:
    timeSeriesFile = argv[0]
    igramsFile = argv[1]
  except:
    Usage() ; sys.exit(1)

#  projectName = os.path.basename(templateFile.partition('.')[0])
#  tssarDirectory = os.getenv('TSSARDIR')
#  projectDir=tssarDirectory+'/'+projectName
#  os.chdir(projectDir)
#  timeSeriesFile='timeseries_'+projectName+'.h5'

  ########################################################
  print "Loading time series: " + timeSeriesFile
  h5timeseries = h5py.File(timeSeriesFile)
  dateList = h5timeseries['timeseries'].keys()
  
  dateIndex={}
  for ni in range(len(dateList)):
    dateIndex[dateList[ni]]=ni

  dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
  nrows,ncols=np.shape(dset)
  timeseries = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0]*np.shape(dset)[1]),np.float32)
  for date in dateList:
    dset = h5timeseries['timeseries'].get(date)
    d = dset[0:dset.shape[0],0:dset.shape[1]]
    timeseries[dateIndex[date]][:]=d.flatten(0)
  del d

#  h5timeseries.close()

  lt,numpixels=np.shape(timeseries)

######################################################
  tbase=np.zeros((lt,1))
  d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])

  for ni in range(len(dateList)):
    d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
    diff = d2-d1
    tbase[ni]=diff.days

#####################################################
  Bperp=np.zeros([lt,1])
  #Bh,Bv = ut.Bh_Bv_timeseries(igramsFile)
  #Theta_Near = float(h5timeseries['timeseries'].attrs['LOOK_REF1'])*np.pi/180.
  #Theta_Far = float(h5timeseries['timeseries'].attrs['LOOK_REF2'])*np.pi/180.
  #Theta=(Theta_Near+Theta_Far)/2
  #Bp=Bh*np.sin(Theta)-Bv*np.cos(Theta) 
  Bp = ut.Baseline_timeseries(igramsFile)
  for i in range(lt):
     Bperp[i][0]=Bp[i]
####################################################
#  bl_list=np.loadtxt(baselineFile,usecols = (0,1),dtype=str)

  
  
#  Bperp=np.zeros([lt,1])
#  print np.shape(bl_list)
#  for i in range(np.shape(bl_list)[0]):
#      Bperp[i]=float(bl_list[i][1])

#  Bperp=Bperp-Bperp[0][0]  
#  print Bperp
#  print Heresh
####################################################
#  teta = (tetaN+tetaF)/2
#  r = (rN+rF)/2
#  teta=19.658799999999999*np.pi/180
#  r=846848.2
#  Bperp=1000*np.random.random((lt,1))
  h5file = h5py.File(igramsFile)
  igramList = h5file['interferograms'].keys()
  try:
     r=float(h5file['interferograms'][igramList[0]].attrs['STARTING_RANGE1'])
     theta=(float(h5file['interferograms'][igramList[0]].attrs['LOOK_REF1'])+float(h5file['interferograms'][igramList[0]].attrs['LOOK_REF2']))/2
     theta=theta*np.pi/180
     C1=Bperp/r/np.sin(theta)
  except:
     C1=Bperp

#  tbase=np.array(tbase)
#################################################  
#  print np.shape(C1)

  M=np.hstack((.5*tbase**2,tbase,np.ones((lt,1))))
  C=np.hstack((M,C1))
  print 'rank of the design matrix : '+str(np.linalg.matrix_rank(C))
  if np.linalg.matrix_rank(C) ==4:
      print 'Design matrix has full rank'
  Cinv=np.linalg.pinv(C)
  par=np.dot(Cinv,timeseries)
  dz=np.zeros([1,numpixels])
  dz[0][:]=par[3][:]
#  print np.shape(dz)
  timeseries=timeseries-np.dot(C1,dz)
  dz=np.reshape(dz,[nrows,ncols])

  print '**************************************'
  print 'writing DEM_error.hgt'
  writefile.write_float32(dz,'DEM_error.hgt')
  f = open('DEM_error.hgt.rsc','w')
  f.write('FILE_LENGTH       '+str(int(nrows))+'\n')
  f.write('WIDTH             '+str(int(ncols))+'\n')  
  print '**************************************'
###################################################
  print 'writing DEM_error.h5'
  h5fileDEM = 'DEM_error.h5'
  h5rmse = h5py.File(h5fileDEM,'w')
  group=h5rmse.create_group('dem')
  dset = group.create_dataset(os.path.basename('dem'), data=dz, compression='gzip')
  for key , value in h5timeseries['timeseries'].attrs.iteritems():
     group.attrs[key]=value
  print '**************************************'
##################################################

  try:
     outname=argv[2]
  except:
     outname=timeSeriesFile.replace('.h5','')+'_demCor.h5'
  print 'writing '+outname
  h5timeseriesDEMcor = h5py.File(outname,'w')
  
  group = h5timeseriesDEMcor.create_group('timeseries')
  for date in dateList:
    print date
    if not date in h5timeseriesDEMcor['timeseries']:
      dset = group.create_dataset(date, data=np.reshape(timeseries[dateIndex[date]][:],[nrows,ncols]), compression='gzip') 

#  group = h5timeseriesDEMcor.create_group('timeseries')
  for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group.attrs[key] = value

  dset1 = h5timeseries['mask'].get('mask')
  Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]] 
  group=h5timeseriesDEMcor.create_group('mask')
  dset = group.create_dataset('mask', data=Mask, compression='gzip')
  
  h5timeseriesDEMcor.close()
  h5timeseries.close()   
if __name__ == '__main__':

  main(sys.argv[1:])  



