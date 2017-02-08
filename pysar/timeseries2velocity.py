#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: Add -m/M/d option
# Yunjun, Jun 2016: Add -t option
# Yunjun, Aug 2015: Support drop_date txt file input


import sys
import os
import getopt
import time
import datetime

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._datetime as ptime


############################################################################
def yyyymmdd2years(date):
    d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])   
    yy = np.float(d.year) + np.float(d.month-1)/12 + np.float(d.day-1)/365
    return yy

############################################################################
def usage():
    print '''
****************************************************************
  Estimating displacement velocity for each pixel.
      It also generates the standadrd deviation of the velocity and the RMSE.

  Usage:
      timeseries2velocity.py -f timeseries.h5 [-o OutputName.h5 -m minimum_date -M maximum_date -E exclude_dates]

      -f: time series h5 file
      -E: excluded dates for velocity estimation, or drop_date.txt
          e.g.: -E 20040502,20060708,20090103
                -E drop_date.txt
                drop_date.txt: 20040502
                               20060708
                               20090103
      -m: earliest date  for velocity estimation
      -M: latest   date  for velocity estimation
      -o: output file name
      -t: template file
          pysar.drop.date  = 20070107,20080712,20090530,20090830,20101203

  Example:
      timeseries2velocity.py timeSeriesFile.h5
      timeseries2velocity.py -f timeseries_ECMWF_demCor_plane.h5 -t $TE/KyushuT73F2980_2990AlosD.template
      timeseries2velocity.py -f timeseries.h5 -m 20080201
      timeseries2velocity.py -f timeseries.h5 -m 20080201 -M 20100508
      timeseries2velocity.py -f timeseries.h5 -m 20080201 -M 20100508 -E 20090703
      timeseries2velocity.py -f timeseries.h5 -E 20040502,20060708,20090103
      timeseries2velocity.py -f timeseries_ECMWF_demCor.h5 -E drop_date.txt

****************************************************************
    '''

############################################################################
def main(argv):

    if len(sys.argv)>2:
        try:   opts, args = getopt.getopt(argv,"f:E:m:M:h:o:t:")
        except getopt.GetoptError:   usage() ; sys.exit(1)
    
        for opt,arg in opts:
            if   opt == '-f':    timeSeriesFile   = arg
            elif opt == '-E':    datesNot2include = arg.replace(' ','').split(',')
            elif opt == '-m':    minDate          = arg
            elif opt == '-M':    maxDate          = arg
            elif opt == '-o':    outName          = arg
            elif opt == '-t':    templateFile     = arg
  
    elif len(sys.argv)==2:
        if   argv[0]=='-h':  usage(); sys.exit(1)
        elif os.path.isfile(argv[0]):   timeSeriesFile = argv[0]
        else:  usage(); sys.exit(1)
    else:  usage(); sys.exit(1)    
  
    ##### Read excluded date list Input
    try:  datesNot2include
    except:
        try:
            templateFile
            templateContents = readfile.read_template(templateFile)
            datesNot2include = templateContents['pysar.drop.date'].replace(' ','').split(',')
        except: pass

    ##############################################################
    print '\n********** Inversion: Time Series to Velocity ***********'
    atr = readfile.read_attribute(timeSeriesFile)
    k = atr['FILE_TYPE']
    print 'input file: '+k
    if not k == 'timeseries':  print 'Input file is not timeseries!'; sys.exit(1)
    print "Loading time series file: " + timeSeriesFile
    h5timeseries = h5py.File(timeSeriesFile)
    dateList1 = sorted(h5timeseries[k].keys())
  
    ##############################################################
    print '--------------------------------------------'
    print 'Dates from input file: '+str(len(dateList1))
    print dateList1
  
    try:
        datesNot2include
        if os.path.isfile(datesNot2include[0]):
            try:  datesNot2include = ptime.read_date_list(datesNot2include[0])
            except:  print 'Can not read date list file: '+datesNot2include[0]
        print '--------------------------------------------'
        print 'Date excluded: '+str(len(datesNot2include))
        print datesNot2include
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
    except: pass
  
    try:
        maxDate
        maxDateyy=yyyymmdd2years(maxDate) 
        print 'maximum date: '+maxDate
        for date in dateList1:
            yy=yyyymmdd2years(date)
            if yy > maxDateyy:
                print '  remove date: '+date
                datesNot2include.append(date)
    except: pass

    try:
        dateList=[]
        for date in dateList1:
            if date not in datesNot2include:
                dateList.append(date)
    except:  pass

    print '--------------------------------------------'
    if len(dateList) == len(dateList1):
        print 'using all dates to calculate the velocity'
    else:
        print 'Dates used to estimate the velocity: '+str(len(dateList))
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
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    lt     = len(dateList)
    timeseries = np.zeros((lt,length,width),np.float32)
    for date in dateList:
        timeseries[dateIndex[date]] = h5timeseries[k].get(date)
  
    numpixels=length*width
    
    Data=np.zeros([lt,numpixels])
    for i in range(lt):
        Data[i,:]=np.reshape(timeseries[i],[1,numpixels])
  
    x=np.dot(B1,Data)
    velocity=np.reshape(x[0,:],[length,width])

    #####################################################
    print 'Calculating rmse'
    Data_linear=np.dot(B,x)
    rmse=np.reshape(np.sqrt((np.sum((Data_linear-Data)**2,0))/lt),[length,width])
    # se=np.reshape((np.sum(np.abs(Data_linear-Data),0)/lt),[length,width])
    # rmse=np.reshape((np.sum((Data_linear-Data)**2,0))/lt,[length,width])
    ######################################################
    print 'Calculating the standard deviation of the estimated velocities'
    residual=Data_linear-Data
    s1=np.sqrt(np.sum(residual**2,0)/(lt-2))
    s2=np.sqrt(np.sum((datevector-np.mean(datevector))**2))
    se=np.reshape(s1/s2,[length,width])
    ######################################################
     
    # SSt=np.sum((Data-np.mean(Data,0))**2,0)
    # SSres=np.sum(residual**2,0)
    # SS_REG=SSt-SSres
    # Rsquared=np.reshape(SS_REG/SSt,[length,width])
    ######################################################  
    # covariance of the velocities
  
    ##### Output File Name
    try:    outName
    except:
        if not datesNot2include == []: outName = 'velocity_ex.h5'
        else:                          outName = 'velocity.h5'
    outName_rmse='rmse_'+outName
    outName_se='std_'+outName
    outName_Rsquared='R2_'+outName
  
    #####################################
    print '--------------------------------------'
    print 'writing to '+outName
    h5velocity = h5py.File(outName,'w')
    group=h5velocity.create_group('velocity')
    dset = group.create_dataset('velocity', data=velocity, compression='gzip')
    group.attrs['date1'] = datevector[0]
    group.attrs['date2'] = datevector[lt-1]
    
    for key , value in atr.iteritems():
        group.attrs[key]=value
    h5velocity.close()  
  
    #####################################
    print 'writing to '+outName_rmse
    h5file = outName_rmse
    h5rmse = h5py.File(h5file,'w')
    group=h5rmse.create_group('rmse')
    dset = group.create_dataset(os.path.basename('rmse'), data=rmse, compression='gzip')
    group.attrs['date1'] = datevector[0]
    group.attrs['date2'] = datevector[lt-1]
  
    for key , value in atr.iteritems():
        group.attrs[key]=value  

    #####################################
    print 'writing to '+outName_se
    h5se = h5py.File(outName_se,'w')
    group=h5se.create_group('rmse')
    dset = group.create_dataset('rmse', data=se, compression='gzip')
    group.attrs['date1'] = datevector[0]
    group.attrs['date2'] = datevector[lt-1]
  
    for key , value in atr.iteritems():
        group.attrs[key]=value
  
    # print 'writing to '+outName_Rsquared
    # h5rsquared = h5py.File(outName_Rsquared,'w')
    # group=h5rsquared.create_group('rmse')
    # dset = group.create_dataset('rmse', data=Rsquared, compression='gzip')
    # group.attrs['date1'] = datevector[0]
    # group.attrs['date2'] = datevector[lt-1]
  
  
    # h5rsquared.close()
    h5se.close()
    h5rmse.close()
    h5timeseries.close()
    print 'Done.'

############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

