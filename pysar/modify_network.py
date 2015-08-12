#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Add: Update Coherence file, Yunjun, Jul 2015
#


import _pysar_utilities as ut
import sys
import os
import getopt
import time
import datetime

from numpy import *
from scipy.io import loadmat
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import h5py

import random
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import _readfile as readfile


def nearest_neighbor(x,y, tbase, pbase):
  """ find nearest neighbour """
  dist = sqrt((tbase -x)**2+(pbase -y)**2) 
  indx=dist==min(dist)
  return indx
  
######################################
def Usage():
  print '''
  ******************************************
  ******************************************
  Modify the network of interferograms and/or coherence.

  usage:
       modify_network.py -f interferogramsFile -s fontsize -w linewidth

    -f : interferograms file stored in hdf5 file format
    -C : coherence file stored in hdf5 file format
    -s : the font size used for x and y labels (default is 12)
    -w : line width used for plotting (default is 2)
    -m : marker size (default is 16)
    -c : marker face color (default is orange) 
    -t : temporal threshold
    -b : baseline threshold
    -d : date (all interferograms that includes date is removed)
    -l : list of interferograms to remove  
    -N : interferogram numbers to remove
    -n : (yes or no)network display to manually choose igrams from the network to remove(default is yes)

    
  Example:

       modify_network.py LoadedData_SanAndreasT356EnvD.h5
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -b 400 -n no
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -b 500 -t 600 -d 080307
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -d '080307 091023'           
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -l 'filt_080307-091023-sim_HDR_8rlks_c10.unw filt_080307-090814-sim_HDR_8rlks_c10.unw' 
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -N '0 4 20 76 89 100'
       modify_network.py -f LoadedData_SanAndreasT356EnvD.h5 -C Coherence_SanAndreasT356EnvD.h5 -N '0 4 20 76 89 100'

  ******************************************
  ******************************************  
  '''


def yymmdd2yyyymmdd(date):
  if date[0][0] == '9':
      date[0] = '19'+date[0]
  else:
      date[0] = '20'+date[0]
  return date
######################################
def main(argv):

  lineWidth=2
  fontSize=12
  markerColor='orange'
  markerSize=16
  networkDisplay='yes'

  if len(sys.argv)>2:

    try:
      opts, args = getopt.getopt(argv,"h:f:C:s:w:m:c:t:b:d:l:n:N:")
      
    except getopt.GetoptError:
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        igramsFile = arg
      elif opt == '-C':
        corFile = arg
      elif opt == '-s':
        fontSize = int(arg)
      elif opt == '-w':
        lineWidth=int(arg)
      elif opt == '-m':
        markerSize=int(arg)
      elif opt == '-c':
        markerColor=arg
      elif opt == '-t':
        temp_thr=float(arg)
      elif opt == '-b':
        base_thr=float(arg)
      elif opt == '-d':
        dates2Rmv=arg
      elif opt == '-l':
        ifgrams_to_rmv=arg
      elif opt == '-n':
        networkDisplay=arg
      elif opt == '-N':
        ifgrams_Number_to_rmv=arg.split()


    try:
      igramsFile
    except:
       Usage() ; sys.exit(1)

  elif len(sys.argv)==2:
     igramsFile = argv[0]
  else:
     Usage() ; sys.exit(1)

############################################
  h5file = h5py.File(igramsFile)
#  import pdb;  pdb.set_trace()
  if h5file.keys()[0] != 'interferograms':
      print 'Input file should be interferograms'
      Usage() ; sys.exit(1)
  ifgramList=h5file['interferograms'].keys()


  try:
    ifgrams_to_rmv
  except:     
    ifgrams_to_rmv=[]
############################################

  try:
    for i in ifgrams_Number_to_rmv:
        print i
        print ifgramList[int(i)]
        ifgrams_to_rmv.append(ifgramList[int(i)])
  except:
    print ''
    
  
  try:
    base_thr
    print 'interferograms with the spatial baseline longer than '+ str(base_thr)+' m is removed'
    for ifgram in  ifgramList:
       Baseline = (float(h5file['interferograms'][ifgram].attrs['P_BASELINE_BOTTOM_HDR'])+float(h5file['interferograms'][ifgram].attrs['P_BASELINE_TOP_HDR']))/2
       if abs(Baseline) > base_thr:
         if not ifgram in ifgrams_to_rmv:
            ifgrams_to_rmv.append(ifgram)
      
  except:
    print 'No Spatial Baseline threshold applied'

#  print ifgrams_to_rmv
  ###########################################################
  #Check if interferograms made of specific dates should be removed
#  print dates2Rmv
  try:
    dates2Rmv
    print 'interferograms with any of following dates will be removed: '+ dates2Rmv
    for ifgram in  ifgramList:
      date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
      
      if (date1 in dates2Rmv) or (date2 in dates2Rmv):
         if not ifgram in ifgrams_to_rmv:
            ifgrams_to_rmv.append(ifgram)
  except:
    print 'No specific dates selected to remove'

###########################################################
  #Check the temporal threshold
  tbase,dateList,dateDict,dateList1=ut.date_list(h5file)
  try:
    temp_thr
    print 'Applying the temporal baseline threshold with threshold of '+str(temp_thr)+' days'
    for ifgram in  ifgramList:
       date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')      
       ind1 = dateList1.index(date1)
       ind2 = dateList1.index(date2)
       dt=tbase[ind2]-tbase[ind1]
       if dt>temp_thr:
          if not ifgram in ifgrams_to_rmv:
            ifgrams_to_rmv.append(ifgram)
  except:
    print 'No Temporal Baseline threshold applied'
############################################################
############################################################
  if networkDisplay=='yes':
  
    tbase,dateList,dateDict,dateList1=ut.date_list(h5file)
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

    datevector=[]
    for i in range(len(dates)):
      datevector.append(np.float(dates[i].year) + np.float(dates[i].month-1)/12 + np.float(dates[i].day-1)/365)
    datevector2=[round(i,2) for i in datevector]

##################################################  
    Bp = ut.Baseline_timeseries(igramsFile)
#############################################################
 
    ifgramList = h5file['interferograms'].keys()
    igram_pairs=np.zeros([len(ifgramList),2],np.int)
    i=0
    for ifgram in  ifgramList:
      date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
      igram_pairs[i][0]=dateList1.index(date1)
      igram_pairs[i][1]=dateList1.index(date2)
      i=i+1

#  h5file.close()
############################################################
    fig1 = plt.figure(1)
    ax1=fig1.add_subplot(111)

    ax1.cla()
 # ax1.plot(dates,Bp, 'o',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    print tbase
    ax1.plot(tbase,Bp, 'o',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    for ni in range(len(ifgramList)):
    #  ax1.plot(array([dates[igram_pairs[ni][0]],dates[igram_pairs[ni][1]]]),array([Bp[igram_pairs[ni][0]],Bp[igram_pairs[ni][1]]]),'k',lw=4)
      ax1.plot(array([tbase[igram_pairs[ni][0]],tbase[igram_pairs[ni][1]]]),array([Bp[igram_pairs[ni][0]],Bp[igram_pairs[ni][1]]]),'k',lw=4) 
 # ax1.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
    ax1.set_ylabel('Bperp [m]',fontsize=fontSize)
    ax1.set_xlabel('Time [years]',fontsize=fontSize)
    ts=datevector[0]+0.2
    te=datevector[-1]+0.2
    ys=int(ts)
    ye=int(te)
    ms=int((ts-ys)*12)
    me=int((te-ye)*12)
    if ms>12:
       ys =ys+1
       ms=1
    if me>12:
       ye =ye+1
       me=1
    
    if ms<1:
       ys =ys-1
       ms=12
    if me<1:
       ye =ye-1
       me=12

    dss=datetime.datetime(ys,ms,1,0,0)
    dee=datetime.datetime(ye,me,1,0,0)
#  ax1.set_xlim(dss,dee)
    ax1.set_ylim(min(Bp)-0.4*abs(min(Bp)),max(Bp)+0.4*max(Bp))

    xticklabels = getp(gca(), 'xticklabels')
    yticklabels = getp(gca(), 'yticklabels')
    setp(yticklabels, 'color', 'k', fontsize=fontSize)
    setp(xticklabels, 'color', 'k', fontsize=fontSize)
 # fig1.autofmt_xdate()
#  plt.show() 
#  ax1.plot(array([tbase[igram_pairs[ni][0]],tbase[igram_pairs[ni][1]]]),array([Bp[igram_pairs[ni][0]],Bp[igram_pairs[ni][1]]]),'r',lw=10)
##########################################  
    x=[]
    y=[]
    Master_index_torremove=[]
    Slave_index_torremove=[]
    a_tbase=array(tbase)
    a_Bp=array(Bp)
#  print igram_pairs
    def onclick(event):
      if event.button==1:
        print 'click'
        xClick = event.xdata
        yClick = event.ydata
        idx=nearest_neighbor(xClick,yClick, a_tbase, a_Bp)       
        xr = a_tbase[idx]
        yr = a_Bp[idx]
        ix=tbase.index(xr)+1
    #  iy=Bp.index(yr)
        print ix
     # print iy
        x.append(xr)
        y.append(yr)
        if mod(len(x),2)==0:
           Master_index_torremove.append(tbase.index(xr))
           ax1.plot([x[len(x)-1],x[len(x)-2]],[y[len(x)-1],y[len(x)-2]],'r',lw=4)
        else:
           Slave_index_torremove.append(tbase.index(xr))
      plt.show()
    cid = fig1.canvas.mpl_connect('button_press_event', onclick)


    plt.show()
    print Master_index_torremove
    print Slave_index_torremove

    if len(Master_index_torremove) == len(Slave_index_torremove):
       R=np.vstack((Master_index_torremove,Slave_index_torremove))
    else:
       R=np.vstack((Master_index_torremove[:-1],Slave_index_torremove))

    R=np.vstack((Master_index_torremove,Slave_index_torremove)) 
    R.sort(0)
    print R
    print dateList1
    numIgrams_rmv=np.shape(R)[1]
    for ifgram in  ifgramList:
       date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
       for i in range(numIgrams_rmv):
           if dateList1[R[0][i]]==date1 and dateList1[R[1][i]]==date2:
               ifgrams_to_rmv.append(ifgram)

  else:
    print 'No network display.'
############################################################
############################################################

#  import pdb;  pdb.set_trace()

  print 'The list of interferograms to remove:' 
  print ifgrams_to_rmv
  igramsFile_modified='Modified_'+igramsFile
  h5filem = h5py.File(igramsFile_modified,'w')
  gg = h5filem.create_group('interferograms')
  ifgram=ifgramList[0]
  unw = h5file['interferograms'][ifgram].get(ifgram)
  MaskZero=np.ones([unw.shape[0],unw.shape[1]])

  print 'writing the modified interferogram file ...'
  for ifgram in  ifgramList:
     if not ifgram in ifgrams_to_rmv:
        print ifgram
        unwSet = h5file['interferograms'][ifgram].get(ifgram)
        unw = unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]        
        MaskZero=unw*MaskZero
        group = gg.create_group(ifgram)
        dset = group.create_dataset(ifgram, data=unw, compression='gzip')
        for key, value in h5file['interferograms'][ifgram].attrs.iteritems():
           group.attrs[key] = value

  Mask=np.ones([unwSet.shape[0],unwSet.shape[1]])
  Mask[MaskZero==0]=0

  # updating Coherence file
  # convert ifgrams_to_rmv to cor_to_rmv
  date12_to_rmv=[]
  for igram in ifgrams_to_rmv:
     date12_to_rmv.append(igram.split('-sim')[0].split('filt_')[-1])

  try:
     corFile
     h5fileCor=h5py.File(corFile)
     corList=h5fileCor['coherence'].keys()

     corFile_modified='Modified_'+corFile
     h5fileCorm=h5py.File(corFile_modified,'w')
     gc = h5fileCorm.create_group('coherence')
     print 'writing the modified coherence file ...'
     for cor in corList:
        date12=cor.split('-sim')[0].split('filt_')[-1]
        if not date12 in date12_to_rmv:
           print cor
           unwSet = h5fileCor['coherence'][cor].get(cor)
           unw = unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]
           group = gc.create_group(cor)
           dset = group.create_dataset(cor, data=unw, compression='gzip')
           for key, value in h5fileCor['coherence'][cor].attrs.iteritems():
              group.attrs[key] = value  
  except:
     print 'No coherence file to be updated.'

########################################################################

  print 'writing Modified_Mask.h5'
  
  h5mask = h5py.File('Modified_Mask.h5','w')
  group=h5mask.create_group('mask')
  dset = group.create_dataset(os.path.basename('mask'), data=Mask, compression='gzip')
  h5mask.close()      

  gm = h5filem.create_group('mask')
  dset = gm.create_dataset('mask', data=Mask, compression='gzip')

  h5file.close()
  h5filem.close()

  
############################################################

if __name__ == '__main__':
  main(sys.argv[1:])



