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
import _readfile as readfile
######################################
def Usage():
  print '''
***************************************************************
***************************************************************

  Tropospheric correction with height-correlation approach

  Usage:

    tropospheric_correction.py -f timeSeriesFile -d demfile -p polynomial order -m maskFile -M maskThreshold -t corelation threshold between dem and data

  Example:
      
    tropospheric_correction.py -f timeseries_demCor.h5 -d radar_8rlks.hgt -p 1 -m temporal_coherence.h5 -M 0.9 -t 0.5
     tropospheric_correction.py -f timeseries_demCor.h5 -d radar_8rlks.hgt -p 2 -m mask.h5 -t 0.5
***************************************************************
***************************************************************
'''

######################################
def main(argv):
    maskThr=0.7
    try:
      opts, args = getopt.getopt(argv,"f:d:p:m:M:t:")

    except getopt.GetoptError:
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if opt == '-f':
        timeSeriesFile = arg
      elif opt == '-d':
        demFile = arg
      elif opt == '-p':
        p = int(arg)
      elif opt == '-m':
        maskFile=arg
        print maskFile
      elif opt == '-M':
        maskThr=float(arg)
      elif opt == '-t':
        corThr=float(arg)
    

    try:
      timeSeriesFile
      demFile
    except:
      Usage() ; sys.exit(1)
    
    try:
       p
    except:
       p=1

    
    try:
       maskFile
    except:
       maskFile='Mask.h5'

    print 'Mask file: ' + maskFile 
###################################################
    h5Mask=h5py.File(maskFile)
    kMask=h5Mask.keys()
    dset = h5Mask[kMask[0]].get(kMask[0])
    Mask = dset[0:dset.shape[0],0:dset.shape[1]]
   # Mask[600:dset.shape[0],:]=0
    Mask=Mask.flatten(1)

    print maskFile
    print maskThr

    if kMask[0]=='mask':
       ndx = Mask !=0
    elif kMask[0]=='temporal_coherence':
       ndx = Mask >maskThr
    else:
       print 'Mask file not recognized!'
       Usage();sys.exit(1)    

    h5Mask.close()
###################################################
    h5timeseries = h5py.File(timeSeriesFile)
    yref=h5timeseries['timeseries'].attrs['ref_y']
    xref=h5timeseries['timeseries'].attrs['ref_x']
###################################################
    if os.path.basename(demFile).split('.')[1]=='hgt':
       amp,dem,demRsc = readfile.read_float32(demFile)

    elif os.path.basename(demFile).split('.')[1]=='dem':
       dem,demRsc = readfile.read_dem(demFile) 

    dem=dem-dem[yref][xref]
   # try:
    print 'considering the look angle of each resolution cell...'
    near_LA=float(h5timeseries['timeseries'].attrs['LOOK_REF1'])
    far_LA=float(h5timeseries['timeseries'].attrs['LOOK_REF2'])
    Length,Width=np.shape(dem)
    LA=np.linspace(near_LA,far_LA,Width)
    LA=np.tile(LA,[Length,1])
    dem=dem/np.cos(LA*np.pi/180.0)       
       
#    except:
 #      print 'Look angle is not considered' 
    dem=dem.flatten(1)
    print np.shape(dem)
###################################################
    if p==1:
        A=np.vstack((dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem,np.ones(len(dem)))).T
    elif p==2: 
        A=np.vstack((dem[ndx]**2,dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem**2,dem,np.ones(len(dem)))).T  
    elif p==3:
        A = np.vstack((dem[ndx]**3,dem[ndx]**2,dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem**3,dem**2,dem,np.ones(len(dem)))).T
    print np.shape(A)

    Ainv=np.linalg.pinv(A)
###################################################
    print 'Estimating the tropospheric effect using the differences of the subsequent epochs and DEM'
    
    dateList = h5timeseries['timeseries'].keys()
    nrows,ncols=np.shape(h5timeseries['timeseries'].get(dateList[0]))
    PAR_EPOCH_DICT_2={} 
    par_diff_Dict={}
    Correlation_Dict={}
    Correlation_Dict[dateList[0]]=0
    Correlation_diff_Dict={}

    for i in range(len(dateList)-1):
       dset1 = h5timeseries['timeseries'].get(dateList[i])
       dset2 = h5timeseries['timeseries'].get(dateList[i+1])
       data1 = dset1[0:dset1.shape[0],0:dset1.shape[1]]
       data2 = dset2[0:dset2.shape[0],0:dset2.shape[1]]
       d = dset2[0:dset2.shape[0],0:dset2.shape[1]] - dset1[0:dset1.shape[0],0:dset1.shape[1]]
        
       del dset1
       del dset2
       d=d.flatten(1)
       data1=data1.flatten(1)
       data2=data2.flatten(1)
       ##############################
       print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       print 'correlation of dem with :'
       
       print '******************************'
       print dateList[i]
       C1=np.zeros([2,len(dem[ndx])])
       C1[0][:]=dem[ndx]
       C1[1][:]=data1[ndx]
       print np.corrcoef(C1)[0][1]
      # Correlation_Dict[dateList[i]]=np.corrcoef(C1)[0][1]
      # print '******************************'
       print dateList[i+1]
       C2=np.zeros([2,len(dem[ndx])])
       C2[0][:]=dem[ndx]
       C2[1][:]=data2[ndx]
       print np.corrcoef(C2)[0][1]
       Correlation_Dict[dateList[i+1]]=np.corrcoef(C2)[0][1]
      # print '******************************'
       print dateList[i]+'-'+dateList[i+1]
       C=np.zeros([2,len(dem[ndx])])
       C[0][:]=dem[ndx]
       C[1][:]=d[ndx]
       print np.corrcoef(C)[0][1]
       print '******************************'
       Correlation_diff_Dict[dateList[i]+'-'+dateList[i+1]]=np.corrcoef(C)[0][1]
       ##############################
      # try:
        # if np.corrcoef(C)[0][1] >= corThr:
        #    par=np.dot(Ainv,d[ndx])
       #  else:
         #   par=[0,0]
      # except:
       par=np.dot(Ainv,d[ndx])
    #   del d
       par_diff_Dict[dateList[i]+'-'+dateList[i+1]]=par  
   #    print par

       try:
         if np.abs(np.corrcoef(C2)[0][1]) >= corThr:        
            PAR2=np.dot(Ainv,data2[ndx])
         else:
            #PAR2=[0,0]
            PAR2=list(np.zeros(p+1))
       except:
            PAR2=np.dot(Ainv,data2[ndx])
 
       PAR_EPOCH_DICT_2[dateList[i+1]]=PAR2
  #     print PAR2
###################################################
    print'****************************************'
    print 'Correlation of DEM with each time-series epoch:'
    average_phase_height_cor=0
    for date in dateList:
       print date + ' : '+str(Correlation_Dict[date])
      # print date + ' : ' + str(Correlation_Dict[date])
       average_phase_height_cor=average_phase_height_cor+np.abs(Correlation_Dict[date])
    print'****************************************'
    print'****************************************'
    print ''
    print 'Average Correlation of DEM with time-series epochs: ' + str(average_phase_height_cor/(len(dateList)-1))
    print ''
    print '****************************************'
    print'****************************************'

   # print 'Correlation of DEM with epoch differences'
   # for key , value in Correlation_diff_Dict.iteritems():
    #   print key+' : '+str(value)
   # print'****************************************'
#    print 'Estimated parameters for each time-series epoch:'
   # for date in dateList:
   #    print date + ' : ' + str(PAR_EPOCH_DICT_2[date])
   # print'****************************************'


###################################################

    #fig=plt.figure(1)
    #ax = fig.add_subplot(3,1,1)
    #ax.plot(dem[ndx],data1[ndx],'o',ms=1)
    #ax = fig.add_subplot(3,1,2)
    #ax.plot(dem[ndx],data2[ndx],'o',ms=1)
    #ax = fig.add_subplot(3,1,3)
    #ax.plot(dem[ndx],d[ndx],'o',ms=1)
    #plt.show()
###################################################
   # print par_diff_Dict
    par_epoch_Dict={}
    par_epoch_Dict[dateList[1]]=par_diff_Dict[dateList[0]+'-'+dateList[1]]

    for i in range(2,len(dateList)):
        par_epoch_Dict[dateList[i]]=par_epoch_Dict[dateList[i-1]]+par_diff_Dict[dateList[i-1]+'-'+dateList[i]]

    print '###############################'
    #for key , value in PAR_EPOCH_DICT_2.iteritems():
    #   print key+' : '+str(value) +' | '+str(par_epoch_Dict[key])
       
  #  print PAR_EPOCH_DICT_2
  #  print par_epoch_Dict
    print '###############################'

    yref=h5timeseries['timeseries'].attrs['ref_y']
    xref=h5timeseries['timeseries'].attrs['ref_x']
    print 'removing the tropospheric delay from each epoch and writing '+timeSeriesFile.split('.')[0]+'_tropCor.h5:'
    h5tropCor = h5py.File(timeSeriesFile.split('.')[0]+'_tropCor.h5','w')
    group = h5tropCor.create_group('timeseries')
    dset = group.create_dataset(dateList[0], data=h5timeseries['timeseries'].get(dateList[0]), compression='gzip')
    for date in dateList:
       if not date in h5tropCor['timeseries']:
         dset = h5timeseries['timeseries'].get(date) 
         data = dset[0:dset.shape[0],0:dset.shape[1]]
        # par=par_epoch_Dict[date]
         par=PAR_EPOCH_DICT_2[date]
        # print np.shape(B)
        # print np.shape(par)

         tropo_effect=np.reshape(np.dot(B,par),[dset.shape[1],dset.shape[0]]).T
         tropo_effect=tropo_effect-tropo_effect[yref][xref]
         dset = group.create_dataset(date, data=data-tropo_effect, compression='gzip')

    for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group.attrs[key] = value
   
 
    dset1 = h5timeseries['mask'].get('mask')
    group=h5tropCor.create_group('mask')
    dset = group.create_dataset('mask', data=dset1, compression='gzip')
   
    h5tropCor.close()
    h5timeseries.close()
    
if __name__ == '__main__':

  main(sys.argv[1:])


