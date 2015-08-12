#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

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


  
######################################
def Usage():
  print '''
  ******************************************
  ******************************************
  Ploting the network of interferograms and 
  the baseline history of SAR acquisitions.

  usage:
       plot_network.py -f interferogramsFile -s fontsize -w linewidth

    -f : interferograms file stored in hdf5 file format
    -s : the font size used for x and y labels (default is 12)
    -w : line width used for plotting (default is 2)
    -m : marker size (default is 16)
    -c : marker face color (default is orange)
    -t : temporal threshold
    -b : baseline threshold
    -d : date (all interferograms with master or slave using the specified date is removed) 
    
  ******************************************
  ******************************************  
  '''



  
######################################
def main(argv):

  lineWidth=2
  fontSize=12
  markerColor='orange'
  markerSize=16
  saveFig='no'
  if len(sys.argv)>2:

    try:
      opts, args = getopt.getopt(argv,"h:f:s:w:m:c:S:")
      
#      for i in range(opts):
 #        if '-f' in opts[0][i]:
  #          fileCheck
   #      Usage() ; sys.exit(1)

    except getopt.GetoptError:
#      print 'No input option by user'
#      print 'runing with default options'
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        igramsFile = arg
      elif opt == '-s':
        fontSize = int(arg)
      elif opt == '-w':
        lineWidth=int(arg)
      elif opt == '-m':
        markerSize=int(arg)
      elif opt == '-c':
        markerColor=arg
      elif opt == '-S':
        saveFig=arg
    try:
      igramsFile
    except:
       Usage() ; sys.exit(1)

  elif len(sys.argv)==2:
     igramsFile = argv[0]
  else:
     Usage() ; sys.exit(1)

#  igramsFile=argv[0]

#############################################
#  Bp = ut.Baseline_timeseries(igramsFile)
 
  h5file = h5py.File(igramsFile)
  if h5file.keys()[0] != 'interferograms':
      print 'Inout file should be interferograms'
      Usage() ; sys.exit(1)
   
 
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
  
  fig2 = plt.figure(2)
  ax2=fig2.add_subplot(111) 
  
  ax2.cla()
  ax2.plot(dates,Bp, '-ko',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
  
  ax2.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
  ax2.set_ylabel('Bperp [m]',fontsize=fontSize)
  ax2.set_xlabel('Time [years]',fontsize=fontSize)
  ts=datevector[0]-0.2
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
  ax2.set_xlim(dss,dee)
  ax2.set_ylim(min(Bp)-0.4*abs(min(Bp)),max(Bp)+0.4*max(Bp))

  xticklabels = getp(gca(), 'xticklabels')
  yticklabels = getp(gca(), 'yticklabels')
  setp(yticklabels, 'color', 'k', fontsize=fontSize)
  setp(xticklabels, 'color', 'k', fontsize=fontSize)
  
  fig2.autofmt_xdate()
  if saveFig=='yes':
     figName='baselineHistory.png'
     plt.savefig(figName)
#############################################################
 

  ifgramList = h5file['interferograms'].keys()
  print 'Number of interferograms: '+str(len(ifgramList))
  igram_pairs=np.zeros([len(ifgramList),2],np.int)
  i=0
  for ifgram in  ifgramList:
     date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
     igram_pairs[i][0]=dateList1.index(date1)
     igram_pairs[i][1]=dateList1.index(date2)
     i=i+1

##########################################################################
# For simulated interferograms only
# To plot the interferograms with unwrapping errors with a different color
  N_unw_err=0
  try:
    for ifgram in  ifgramList:
      if h5file['interferograms'][ifgram].attrs['unwrap_error']=='yes':
         N_unw_err=N_unw_err+1

  except:
    print ""
      
  if N_unw_err>0:
     igram_pairs_ue=np.zeros([N_unw_err,2],np.int)
     i=0
     for ifgram in  ifgramList:
       if h5file['interferograms'][ifgram].attrs['unwrap_error']=='yes':
         date1,date2 = h5file['interferograms'][ifgram].attrs['DATE12'].split('-')
         igram_pairs_ue[i][0]=dateList1.index(date1)
         igram_pairs_ue[i][1]=dateList1.index(date2)
         i=i+1


  h5file.close()
##########################################################################
  fig1 = plt.figure(1)
  ax1=fig1.add_subplot(111)

  ax1.cla()
  ax1.plot(dates,Bp, 'o',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)

  for ni in range(len(ifgramList)):
    ax1.plot(array([dates[igram_pairs[ni][0]],dates[igram_pairs[ni][1]]]),array([Bp[igram_pairs[ni][0]],Bp[igram_pairs[ni][1]]]),'k',lw=lineWidth)
  
  if N_unw_err>0:
     for ni in range(N_unw_err):
        ax1.plot(array([dates[igram_pairs_ue[ni][0]],dates[igram_pairs_ue[ni][1]]]),array([Bp[igram_pairs_ue[ni][0]],Bp[igram_pairs_ue[ni][1]]]),'r',lw=lineWidth)
  
  ax1.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
  ax1.set_ylabel('Bperp [m]',fontsize=fontSize)
  ax1.set_xlabel('Time [years]',fontsize=fontSize)
  ts=datevector[0]-0.2
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
  ax1.set_xlim(dss,dee)
  ax1.set_ylim(min(Bp)-0.4*abs(min(Bp)),max(Bp)+0.4*max(Bp))

  xticklabels = getp(gca(), 'xticklabels')
  yticklabels = getp(gca(), 'yticklabels')
  setp(yticklabels, 'color', 'k', fontsize=fontSize)
  setp(xticklabels, 'color', 'k', fontsize=fontSize)
  fig1.autofmt_xdate()
 
  if saveFig=='yes':
     figName='igramsNetwork.png'
     plt.savefig(figName)

  plt.show() 

############################################################

if __name__ == '__main__':
  main(sys.argv[1:])



