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
from numpy import *
from scipy.io import loadmat
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import h5py
from mpl_toolkits.axes_grid.inset_locator import inset_axes


def Usage():
    print '''
***************************************************************
***************************************************************    
Time-series viewer

Usage:

         tsviewer.py -f timeseriesFile.h5 -v velocityFile.h5 -l lower bound -h higher bound -s fontsize -m Marker Size -c marker color -w linewidth -u unit


        -f : file of the timeseries
        -v : velocity file (if not specified then the last time-series epoch is displayed)
        -l : lower bound of the displacement (default is minimum of the displacemen)
        -h : higher bound of the displacemet (default is max of the displacemen)
        -s : size of font used x and y labels (default = 22)
        -m : marker size (default = 16)
        -c : color of the markers (default = red). some options are: orange, black, yellow, blue, green...
        -w : width of lines to connect the points (default = 2 ). set to 0 (-l 0) if you don't want any line connecting the points
        -u : unit of the displacement (default = cm). Other optons are: mm and m
        -e : event dates 
        -a : lower bound of the colorscale to display the velocity to display
        -b : higher bound of the colorscale to display the velocity to display
        -F : another  timeseries file (can be used to compare 2 time-series)
        -S : save to matlab [default: no]
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Example:
         tsviewer.py timeseries.h5
         tsviewer.py -f timeseries.h5
         tsviewer.py -f timeseries_demCor.h5 -v velocity_masked.h5 -u m -c blue
         tsviewer.py -f timeseries.h5  -v velocity.h5 -s 24 -m 12 -c orange -l -10 -h 10 -w 4 -u mm 
         tsviewer.py -f timeseries.h5  -F timeseries_tropCor.h5 
         tsviewer.py -f timeseries.h5  -v velocity.h5 -a -0.01 -b 0.01
         tsviewer.py -f timeseries.h5 -S yes
***************************************************************
***************************************************************
    '''


def main(argv):

  markerSize=16
  markerSize2=16
  markerColor='g'
  markerColor2='red'
  lineWidth=2
  fontSize=22
  unit='cm'
  Save_timeseries='no'

  if len(sys.argv)>2:

    try:
      opts, args = getopt.getopt(argv,"f:F:v:a:b:s:m:c:w:u:l:h:S:")

    except getopt.GetoptError:
      Usage() ; sys.exit(1)
  
    for opt,arg in opts:
      if opt == '-f':
        timeSeriesFile = arg
      elif opt == '-F':
        timeSeriesFile_2 = arg
      elif opt == '-v':
        velocityFile = arg
      elif opt == '-a':
        vmin = float(arg)
      elif opt == '-b':
        vmax = float(arg)
      elif opt == '-s':
        fontSize = int(arg)
      elif opt == '-m':
        markerSize=int(arg)
        markerSize2=int(arg)
      elif opt == '-S':
        Save_timeseries=arg
      elif opt == '-c':
        markerColor=arg
      elif opt == '-w':
        lineWidth=int(arg)
      elif opt == '-u':
        unit=arg
      elif opt == '-l':
        lbound=float(arg)
      elif opt == '-h':
        hbound=float(arg)


  elif len(sys.argv)==2:
    if argv[0]=='-h':
       Usage(); sys.exit(1)
    elif os.path.isfile(argv[0]):
       timeSeriesFile = argv[0]
       h5timeseries = h5py.File(timeSeriesFile)
       if not 'timeseries' in h5timeseries.keys():
          print 'ERROR'
          Usage(); sys.exit(1)
    else:
       
       Usage(); sys.exit(1)
       
  elif len(sys.argv)<2:
    Usage(); sys.exit(1)



  if unit in ('m','M'):
     unitFac=1
  elif unit in ('cm','Cm','CM'):
     unitFac=100
  elif unit in ('mm','Mm','MM','mM'):
     unitFac=1000
  else:
     print 'Warning:'
     print 'wrong unit input!'
     print 'cm is considered to display the displacement'
############################################

  
  if not os.path.isfile(timeSeriesFile):
       Usage();sys.exit(1)

  h5timeseries = h5py.File(timeSeriesFile)
  if not 'timeseries' in h5timeseries.keys():
          Usage(); sys.exit(1)

 
  dateList = h5timeseries['timeseries'].keys()

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
###########################################  
#  eventDates=['20041223','20051003']
 # try:
 #   eventDates
 #   events=[]
  #  for ni in range(len(eventDates)):
  #    d = datetime.datetime(*time.strptime(eventDates[ni],"%Y%m%d")[0:5])
  #    events.append(d)
 # except:
  #  print ''
  #print events
###########################################
  try:
     velocityFile
     h5file=h5py.File(velocityFile,'r')
     k=h5file.keys()
     dset= h5file[k[0]].get(k[0])
     print 'The file to display is: ' + k[0]
  except:
     dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[-1])
 # timeseries = np.zeros((len(h5timeseries['timeseries'].keys()),np.shape(dset)[0],np.shape(dset)[1]),np.float32)
 # for date in h5timeseries['timeseries'].keys():
 #   timeseries[dateIndex[date]] = h5timeseries['timeseries'].get(date)
  
###########################################
  
  fig = plt.figure()
  ax=fig.add_subplot(111)
  try:
    vmin
    vmax
    img=ax.imshow(dset,vmin=vmin,vmax=vmax)
  except:
    img=ax.imshow(dset)
  fig2 = plt.figure(2)
  ax2=fig2.add_subplot(111) 
 # print dates
 # print dateList

  try:
     timeSeriesFile_2
     h5timeseries_2=h5py.File(timeSeriesFile_2)
  except:
     print""   
##########################################  
  def onclick(event):
    if event.button==1:
      print 'click'
      xClick = int(event.xdata)
      yClick = int(event.ydata)
      Dis=[]
      for date in h5timeseries['timeseries'].keys():
             Dis.append( h5timeseries['timeseries'].get(date)[yClick][xClick])    
      ax2.cla()
      
      
      try:
         Dis2=[]
         for date in dateList:
             Dis2.append( h5timeseries_2['timeseries'].get(date)[yClick][xClick])
         dis2=array(Dis2)
         dis2=dis2*unitFac
         ax2.plot(dates,dis2, '^',ms=markerSize2, alpha=0.7, mfc=markerColor2)
      except:
         Dis2=[]
      
    #  ax2.plot(dates,dis, '-ko',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
      dis=array(Dis)
      if Save_timeseries in ['yes','y','YES','Yes']:
         import scipy.io as sio
         Delay={}
         Delay['displacement']=dis
         Delay['unit']='m'
         Delay['time']=datevector
         sio.savemat('displacement.mat', {'displacement': Delay})
      dis=dis*unitFac
      ax2.plot(dates,dis, '-ko',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
     # print dis
     # print dates
      print dset[yClick][xClick]

      ax2.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
      if unitFac==100:
        ax2.set_ylabel('Displacement [cm]',fontsize=fontSize)
      elif unitFac==1000:
        ax2.set_ylabel('Displacement [mm]',fontsize=fontSize)
      else:
        ax2.set_ylabel('Displacement [m]',fontsize=fontSize)

      ax2.set_xlabel('Time [years]',fontsize=fontSize)
      ds=datevector[0]-0.2
      de=datevector[-1]+0.2
      ys=int(ds)
      ye=int(de)
      ms=int((ds-ys)*12)+1
      me=int((de-ye)*12)+1
      
      
      
      dss=datetime.datetime(ys,ms,1,0,0)
      dee=datetime.datetime(ye,me,1,0,0)
      ax2.set_xlim(dss,dee)
      
      try:
        lbound
        hbound
        ax2.set_ylim(lbound,hbound)
          
      except: 
        ax2.set_ylim(min(dis)-0.4*abs(min(dis)),max(dis)+0.4*max(dis))

      for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontSize)
      for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontSize)
                # specify integer or one of preset strings, e.g.
                #tick.label.set_fontsize('x-small')
               # tick.label.set_rotation('vertical')


      fig2.autofmt_xdate()
 
      plt.show()
      
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  plt.show()

if __name__ == '__main__':

  main(sys.argv[1:])
