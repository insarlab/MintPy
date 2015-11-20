#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Add -M option, Yunjun, Aug 2015
#

import sys
import os
import getopt
import time
import datetime
#import h5py
import numpy as np
#import matplotlib.pyplot as plt

def yyyymmdd2years(date):
  d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])   
  yy = np.float(d.year) + np.float(d.month-1)/12 + np.float(d.day-1)/365
  return yy
######################################
def Usage():

  print '''
****************************************************************
****************************************************************
    Estimating displacement velocity for each pixel.
    It also generates the standadrd deviation of the velocity and the RMSE.

    Usage:

         timeseries2velocity.py -f timeseries.h5 -o OutputName.h5 -m minimum_date -M maximum_date -d exclude_dates
    
    Example:

         timeseries2velocity.py timeSeriesFile.h5
         timeseries2velocity.py -f timeseries.h5 -m 20080201
         timeseries2velocity.py -f timeseries.h5 -m 20080201 -M 20100508
         timeseries2velocity.py -f timeseries.h5 -d '20040502 20060708 20090103'
         timeseries2velocity.py -f timeseries.h5 -o velocity_demCor.h5  
         timeseries2velocity.py -f timeseries.h5 -m 20080201 -M 20100508 -d '20090703'

****************************************************************
****************************************************************
   '''


######################################
def main(argv):

  if len(sys.argv)>2:

    try:
      opts, args = getopt.getopt(argv,"f:d:m:M:h:o:")
       
    except getopt.GetoptError:
      Usage() ; sys.exit(1)
       
    for opt,arg in opts:
      if opt == '-f':
        timeSeriesFile = arg
      elif opt == '-d':
        datesNot2include = arg.split()
      elif opt == '-m':
        minDate = arg
      elif opt == '-M':
        maxDate = arg
      elif opt == '-o':
        outName = arg
      
  elif len(sys.argv)==2:
    if argv[0]=='-h':
       Usage(); sys.exit(1)
    elif os.path.isfile(argv[0]):
       timeSeriesFile = argv[0]
    else:
       Usage(); sys.exit(1)

  else:
    Usage(); sys.exit(1)    

##############################################################
  print "Loading time series file: " + timeSeriesFile
  import h5py
  h5timeseries = h5py.File(timeSeriesFile)
  dateList1 = h5timeseries['timeseries'].keys()

##############################################################
  print 'All dates existed:'
  print dateList1
  print '*******************'

  try:
    datesNot2include
    print 'exclude dates: '+str(datesNot2include)
  except:
    datesNot2include=[]

  try:
    minDate
    minDateyy=yyyymmdd2years(minDate)
    print 'minimum date: '+minDate
    for date in dateList1:
       yy=yyyymmdd2years(date)
       if yy < minDateyy:
           print '  remove date: '+date
           datesNot2include.append(date)
  except:
    print ''

 # maxDate='20100521'
  try:
    maxDate
    maxDateyy=yyyymmdd2years(maxDate) 
    print 'maximum date: '+maxDate
    for date in dateList1:
       yy=yyyymmdd2years(date)
       if yy > maxDateyy:
           print '  remove date: '+date
           datesNot2include.append(date)
  except:
    print ''

  try:
    # datesNot2include = '20100903 20100730 20100625 20100521 20100416'
     dateList=[]
     for date in dateList1:
        if date not in datesNot2include:
           dateList.append(date)
  except:
     dateList=dateList1
     print 'using all dates to calculate the vlocity'
  print '--------------------------------------------'
  print 'dates used to estimate the velocity:'
  print dateList
  print '--------------------------------------------'
##############################################################
  dateIndex={}
  for ni in range(len(dateList)):
    dateIndex[dateList[ni]]=ni
  tbase=[]
  d1 = datetime.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
  
  for ni in range(len(dateList)):
    d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
    diff = d2-d1
    tbase.append(diff.days)

  dates=[]
  for ni in range(len(dateList)):
    d = datetime.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
    dates.append(d)
#  print 'Index and dates from ' + timeSeriesFile
#  for ni in range(len(dates)):
#    print ni,dates[ni]

###########################################
  print 'Calculating Velocity'

  datevector=[]
  for i in range(len(dates)):
    datevector.append(np.float(dates[i].year) + np.float(dates[i].month-1)/12 + np.float(dates[i].day-1)/365)

  B=np.ones([len(datevector),2])
  B[:,0]=datevector
  #B1 = np.linalg.pinv(B)
  B1 = np.dot(np.linalg.inv(np.dot(B.T,B)),B.T)
  B1 = np.array(B1,np.float32)
#########################################
    

  dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[0])
#  timeseries = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0],np.shape(dset)[1]),np.float32)
#  for date in h5timeseries['timeseries'].keys():
#    timeseries[dateIndex[date]] = h5timeseries['timeseries'].get(date)

  timeseries = np.zeros((len(dateList),np.shape(dset)[0],np.shape(dset)[1]),np.float32)
  for date in dateList:
    timeseries[dateIndex[date]] = h5timeseries['timeseries'].get(date)


  lt,rows,cols=np.shape(timeseries)
  numpixels=rows*cols
  
  Data=np.zeros([lt,numpixels])
  for i in range(lt):
     Data[i,:]=np.reshape(timeseries[i],[1,numpixels])

  x=np.dot(B1,Data)
  velocity=np.reshape(x[0,:],[rows,cols])
#  import matplotlib.pyplot as plt
#  plt.imshow(velocity,vmin=-0.02, vmax=.02)
#  plt.colorbar()
#  plt.show()
#####################################################
  print 'Calculating rmse'
  Data_linear=np.dot(B,x)
  rmse=np.reshape(np.sqrt((np.sum((Data_linear-Data)**2,0))/lt),[rows,cols])
 # se=np.reshape((np.sum(np.abs(Data_linear-Data),0)/lt),[rows,cols])
 # rmse=np.reshape((np.sum((Data_linear-Data)**2,0))/lt,[rows,cols])
######################################################
  print 'Calculating the standard deviation of the estimated velocities'
  residual=Data_linear-Data
  s1=np.sqrt(np.sum(residual**2,0)/(lt-2))
  s2=np.sqrt(np.sum((datevector-np.mean(datevector))**2))
  se=np.reshape(s1/s2,[rows,cols])
######################################################
   
 # SSt=np.sum((Data-np.mean(Data,0))**2,0)
 # SSres=np.sum(residual**2,0)
 # SS_REG=SSt-SSres
 # Rsquared=np.reshape(SS_REG/SSt,[rows,cols])
######################################################  
  # covariance of the velocities
  
######################################################
#  h5file = projectDir+'/velocity_'+projectName+'.h5'
 # print 'saving results to hdf5 file'


  try:
    outName
    outName_rmse='rmse_'+outName
    outName_se='std_'+outName
    outName_Rsquared='R2_'+outName
  except:
    outName='velocity.h5'
    outName_rmse='rmse_velocity.h5'
    outName_se='std_velocity.h5'
    outName_Rsquared='R2_velocity.h5'


 # try:
   # h5file = argv[1]
  #  print 'writing velocity to '+argv[1]
 # except:
   # h5file = 'velocity.h5'
   # print 'writing to velocity.h5'
  print '--------------------------------------'
  print 'writing to '+outName
  h5velocity = h5py.File(outName,'w')
  group=h5velocity.create_group('velocity')
  dset = group.create_dataset('velocity', data=velocity, compression='gzip')
  group.attrs['date1'] = datevector[0]
  group.attrs['date2'] = datevector[lt-1]
  
  for key , value in h5timeseries['timeseries'].attrs.iteritems():
     group.attrs[key]=value
  h5velocity.close()  
#  h5timeseries.close()

  print '--------------------------------------'
  print 'writing to '+outName_rmse
  h5file = outName_rmse
  h5rmse = h5py.File(h5file,'w')
  group=h5rmse.create_group('rmse')
  dset = group.create_dataset(os.path.basename('rmse'), data=rmse, compression='gzip')
  group.attrs['date1'] = datevector[0]
  group.attrs['date2'] = datevector[lt-1]


  for key , value in h5timeseries['timeseries'].attrs.iteritems():
     group.attrs[key]=value  

  print '--------------------------------------'
  print 'writing to '+outName_se
  h5se = h5py.File(outName_se,'w')
  group=h5se.create_group('rmse')
  dset = group.create_dataset('rmse', data=se, compression='gzip')
  group.attrs['date1'] = datevector[0]
  group.attrs['date2'] = datevector[lt-1]

  for key , value in h5timeseries['timeseries'].attrs.iteritems():
     group.attrs[key]=value

  print '--------------------------------------'
 # print 'writing to '+outName_Rsquared
 # h5rsquared = h5py.File(outName_Rsquared,'w')
 # group=h5rsquared.create_group('rmse')
 # dset = group.create_dataset('rmse', data=Rsquared, compression='gzip')
 # group.attrs['date1'] = datevector[0]
 # group.attrs['date2'] = datevector[lt-1]

#  for key , value in h5timeseries['timeseries'].attrs.iteritems():
#     group.attrs[key]=value

#  h5rsquared.close()
  h5se.close()
  h5rmse.close()
  h5timeseries.close()
if __name__ == '__main__':
    
  main(sys.argv[1:])

