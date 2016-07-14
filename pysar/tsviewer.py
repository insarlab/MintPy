#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Add minDate/maxDate/excludeDateList/contour option, Yunjun, Aug 2015
# Add radius option and errorbar in the plot, Yunjun, Aug 2015
# Add x/y/dispVelFig/dispTsFig option, Yunjun, Sep 2015

import sys
import os
import getopt
import time
import datetime
from numpy import *
#from scipy.io import loadmat
#import matplotlib
#import matplotlib.pyplot as plt
from pylab import *
import h5py
#from mpl_toolkits.axes_grid.inset_locator import inset_axes

def yyyymmdd2years(date):
  d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])
  yy = float(d.year) + float(d.month-1)/12 + float(d.day-1)/365
  return yy

def Usage():
    print '''
***************************************************************
***************************************************************    
  Time-series Viewer

  Usage:
         tsviewer.py -f timeseriesFile.h5 -v velocityFile.h5 -l lower bound -h higher bound -s fontsize -m Marker Size -c marker color -w linewidth -u unit

        -f : file of the timeseries
        -v : velocity file (if not specified then the last time-series epoch is displayed)
        -l : lower bound of the displacement [default is min of the displacemen]
        -h : higher bound of the displacemet [default is max of the displacemen]
        -e : event dates (not enabled yet) 
        -a : lower bound of the colorscale to display the velocity to display
        -b : higher bound of the colorscale to display the velocity to display
        -F : another  timeseries file (can be used to compare 2 time-series)

        -s : size of font used x and y labels [default = 22]
        -m : marker size [default = 16]
        -c : color of the markers [default = green]. some options are: orange, black, yellow, blue, green...
        -w : width of lines to connect the points [default = 2]. set to 0 (-l 0) if you don't want any line connecting the points
        -u : unit of the displacement [default = cm]. Other optons are: mm and m

        -r : radius of selecting square, display mean value from (+/-radius,+/-radius). [default is 0 - one point]
        -x : x coordiante (range) of selection
        -y : y coordinate (range) of selection

        -D : dem file
        -C : display Contour: no, yes, or only ('no' for DEM only, 'yes' for DEM basemap and Contour, 'only' for contour only). [Default is 'only']
        -V : contour step for display [default is 200 m]

        -t : minimum date for display
        -T : maximum date for display
        -d : exclude dates list for display

        -S : save data to matlab and plot to pdf [default is no]
        -P : display time series plot [default is yes. Set -P no -S yes to save to pdf file without display]
        -p : display the map of velocity or last timeseries interferogram [default is yes, except when x/y are setted]

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Example:
         tsviewer.py timeseries.h5
         tsviewer.py -f timeseries.h5
         tsviewer.py -f timeseries.h5 -v velocity_masked.h5 -u m -c blue
         tsviewer.py -f timeseries.h5 -v velocity.h5 -s 24 -m 12 -c orange -l -10 -h 10 -w 4 -u mm 
         tsviewer.py -f timeseries.h5 -F timeseries_tropCor.h5 
         tsviewer.py -f timeseries.h5 -v velocity.h5 -a -0.01 -b 0.01
         tsviewer.py -f timeseries.h5 -S yes
         tsviewer.py -f timeseries.h5 -v velocity.h5 -a -0.05 -b 0.05 -D Andreas.dem -C yes -l -5 -h 15
         tsviewer.py -f timeseries.h5 -v velocity.h5 -D Andreas.dem
         tsviewer.py -f timeseries.h5 -v velocity.h5 -t 20100102 -T 20101120 -d '20100520 20100705'
         tsviewer.py -f timeseries.h5 -v velocity.h5 -r 10
         tsviewer.py -f timeseries.h5 -x 300 -y 500 -r 10
         tsviewer.py -f timeseries.h5 -x 300:330 -y 500:530
         tsviewer.py -f timeseries.h5 -x 300:330 -y 500:530 -S yes -P no

***************************************************************
***************************************************************
    '''


def main(argv):

  #default settings
  markerSize=16
  markerSize2=16
  markerColor='g'
  markerColor2='red'
  lineWidth=2
  fontSize=16
  unit='cm'
  Save_timeseries='no'
  dispTsFig='yes'
  dispVelFig='yes'
  dispContour='only'
  contour_step=200
  smoothContour='no'
  radius=0;
  edgeWidth=1.5
  fig_dpi=300

  if len(sys.argv)>2:
    try:
      opts, args = getopt.getopt(argv,"f:F:v:a:b:s:m:c:w:u:l:h:S:D:C:V:t:T:d:r:x:y:P:p:")
    except getopt.GetoptError:
      Usage() ; sys.exit(1)
 
    for opt,arg in opts:
      if   opt == '-f':     timeSeriesFile = arg
      elif opt == '-F':     timeSeriesFile_2 = arg
      elif opt == '-v':     velocityFile = arg
      elif opt == '-a':     vmin = float(arg)
      elif opt == '-b':     vmax = float(arg)
      elif opt == '-s':     fontSize = int(arg)
      elif opt == '-m':     markerSize=int(arg);       markerSize2=int(arg)
      elif opt == '-S':     Save_timeseries=arg
      elif opt == '-c':     markerColor=arg
      elif opt == '-w':     lineWidth=int(arg)
      elif opt == '-u':     unit=arg
      elif opt == '-l':     lbound=float(arg)
      elif opt == '-h':     hbound=float(arg)
      elif opt == '-D':     demFile=arg
      elif opt == '-C':     dispContour=arg
      elif opt == '-V':     contour_step=float(arg)
      elif opt == '-t':     minDate=arg
      elif opt == '-T':     maxDate=arg
      elif opt == '-d':     datesNot2show = arg.split()
      elif opt == '-r':     radius=abs(int(arg))
      elif opt == '-x':     xsub = [int(i) for i in arg.split(':')];   xsub.sort();   dispVelFig='no'
      elif opt == '-y':     ysub = [int(i) for i in arg.split(':')];   ysub.sort();   dispVelFig='no'
      elif opt == '-P':     dispTsFig=arg
      elif opt == '-p':     dispVelFig=arg


  elif len(sys.argv)==2:
    if argv[0]=='-h':
       Usage(); sys.exit(1)
    elif os.path.isfile(argv[0]):
       timeSeriesFile = argv[0]
       h5timeseries = h5py.File(timeSeriesFile)
       if not 'timeseries' in h5timeseries.keys():
          print 'ERROR'
          Usage(); sys.exit(1)
    else:  Usage(); sys.exit(1)
  elif len(sys.argv)<2:
    Usage(); sys.exit(1)

  if   unit in ('m','M'):              unitFac=1
  elif unit in ('cm','Cm','CM'):       unitFac=100
  elif unit in ('mm','Mm','MM','mM'):  unitFac=1000
  else:
     print 'Warning:'
     print 'wrong unit input!'
     print 'cm is considered to display the displacement'

##############################################################
# Read time series file info

  if not os.path.isfile(timeSeriesFile):
     Usage();sys.exit(1)

  h5timeseries = h5py.File(timeSeriesFile)
  if not 'timeseries' in h5timeseries.keys():
     Usage(); sys.exit(1)
 
  dateList1 = h5timeseries['timeseries'].keys()

##############################################################
# Dates to show time series plot

  import matplotlib.dates as mdates
  years    = mdates.YearLocator()   # every year
  months   = mdates.MonthLocator()  # every month
  yearsFmt = mdates.DateFormatter('%Y')

  print '*******************'
  print 'All dates existed:'
  print dateList1
  print '*******************'

  try:
     datesNot2show
     print 'dates not to show: '+str(datesNot2show)
  except:  datesNot2show=[]

  try:
    minDate
    minDateyy=yyyymmdd2years(minDate)
    print 'minimum date: '+minDate
    for date in dateList1:
       yy=yyyymmdd2years(date)
       if yy < minDateyy:
           datesNot2show.append(date)
  except:  pass
  try:
    maxDate
    maxDateyy=yyyymmdd2years(maxDate)
    print 'maximum date: '+maxDate
    for date in dateList1:
       yy=yyyymmdd2years(date)
       if yy > maxDateyy:
           datesNot2show.append(date)
  except:  pass

  try:
     dateList=[]
     for date in dateList1:
        if date not in datesNot2show:
           dateList.append(date)
     print '--------------------------------------------'
     print 'dates used to show time series displacements:'
     print dateList
     print '--------------------------------------------'
  except:
     dateList=dateList1
     print 'using all dates to show time series displacement'

###################################################################
# Date info

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
# Plot Fig 1 - Velocity / last epoch of time series / DEM

  import matplotlib.pyplot as plt
  if dispVelFig in ('yes','Yes','y','Y','YES'):
     fig = plt.figure()
     ax=fig.add_subplot(111)

     try:
        velocityFile
        h5file=h5py.File(velocityFile,'r')
        k=h5file.keys()
        dset= h5file[k[0]].get(k[0])
        print 'display: ' + k[0]
     except:
        dset = h5timeseries['timeseries'].get(h5timeseries['timeseries'].keys()[-1])
        print 'display: last epoch of timeseries'

     #DEM/contour option
     try:
        demFile
        import _readfile as readfile
        if   os.path.basename(demFile).split('.')[1]=='hgt':  amp,dem,demRsc = readfile.read_float32(demFile)
        elif os.path.basename(demFile).split('.')[1]=='dem':  dem,demRsc = readfile.read_dem(demFile)

        if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):
           print 'show DEM as basemap'
           cmap_dem=plt.get_cmap('gray')
           import _pysar_utilities as ut
           plt.imshow(ut.hillshade(dem,50.0),cmap=cmap_dem)
        if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):
           print 'show contour'
           if smoothContour in ('yes','Yes','y','Y','YES'):
              import scipy.ndimage as ndimage
              dem=ndimage.gaussian_filter(dem,sigma=10.0,order=0)
           contour_sequence=np.arange(-6000,9000,contour_step)
           plt.contour(dem,contour_sequence,origin='lower',colors='black',alpha=0.5)
     except: print 'No DEM file' 

     try:     img=ax.imshow(dset,vmin=vmin,vmax=vmax)
     except:  img=ax.imshow(dset)

     import matplotlib.patches as patches      # need for draw rectangle of points selected on VelFig

########################################## 
# Plot Fig 2 - Time series plot
  import scipy.stats as stats
  fig2 = plt.figure(2)
  ax2=fig2.add_subplot(111) 

  try:
     timeSeriesFile_2
     h5timeseries_2=h5py.File(timeSeriesFile_2)
     print 'plot 2nd time series'
  except:  pass   

  ########### Plot Time Series with x/y ##########
  try:
     xsub
     ysub
     try:     xmin=xsub[0];         xmax=xsub[1]+1;         print 'x='+str(xsub[0])+':'+str(xsub[1])
     except:  xmin=xsub[0]-radius;  xmax=xsub[0]+radius+1;  print 'x='+str(xsub[0])+'+/-'+str(radius)
     try:     ymin=ysub[0];         ymax=ysub[1]+1;         print 'y='+str(ysub[0])+':'+str(ysub[1])
     except:  ymin=ysub[0]-radius;  ymax=ysub[0]+radius+1;  print 'y='+str(ysub[0])+'+/-'+str(radius)
     try:
        fig
        rectSelect=patches.Rectangle((xmin,ymin),radius*2+1,radius*2+1,fill=False,lw=edgeWidth)
        ax.add_patch(rectSelect)
     except: pass

     Dis=[]
     for date in dateList:  Dis.append(h5timeseries['timeseries'].get(date)[ymin:ymax,xmin:xmax])
     Dis0=array(Dis)
     dis=Dis0*unitFac
     dis=reshape(dis,(len(dateList),-1))
     dis_mean=stats.nanmean(dis,1)
     if (xmax-xmin)*(ymax-ymin)==1:  dis_std=[0]*len(dateList)
     else:                           dis_std=stats.nanstd(dis,1)
     (_, caps, _)=ax2.errorbar(dates,dis_mean,yerr=dis_std,fmt='-ko',\
                               ms=markerSize, lw=lineWidth, alpha=1, mfc=markerColor,\
                               elinewidth=edgeWidth,ecolor='black',capsize=markerSize*0.5)
     for cap in caps:  cap.set_markeredgewidth(edgeWidth)
     print dis_mean

     # x axis format
     ax2.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
     if unitFac==100:     ax2.set_ylabel('Displacement [cm]',fontsize=fontSize)
     elif unitFac==1000:  ax2.set_ylabel('Displacement [mm]',fontsize=fontSize)
     else:                ax2.set_ylabel('Displacement [m]' ,fontsize=fontSize)
     ax2.set_xlabel('Time [years]',fontsize=fontSize)
     ax2.set_title('x='+str(xmin)+':'+str(xmax-1)+', y='+str(ymin)+':'+str(ymax-1))
     ax2.xaxis.set_major_locator(years)
     ax2.xaxis.set_major_formatter(yearsFmt)
     ax2.xaxis.set_minor_locator(months)
     datemin = datetime.date(int(datevector[0]),1,1)
     datemax = datetime.date(int(datevector[-1])+1,1,1)
     ax2.set_xlim(datemin, datemax)

     # y axis format
     try:
        lbound
        hbound
        ax2.set_ylim(lbound,hbound)
     except:
        ax2.set_ylim(nanmin(dis_mean-dis_std)-0.4*abs(nanmin(dis_mean)),\
                     nanmax(dis_mean+dis_std)+0.4*abs(nanmax(dis_mean)))

     for tick in ax2.xaxis.get_major_ticks():  tick.label.set_fontsize(fontSize)
     for tick in ax2.yaxis.get_major_ticks():  tick.label.set_fontsize(fontSize)
     #fig2.autofmt_xdate()     #adjust x overlap by rorating, may enble again

     if Save_timeseries in ('yes','Yes','Y','y','YES'):
        import scipy.io as sio
        Delay={}
        Delay['displacement']=Dis0
        Delay['unit']='m'
        Delay['time']=datevector
        tsNameBase='ts_x'+str(xmin)+'_'+str(xmax-1)+'y'+str(ymin)+'_'+str(ymax-1)
        sio.savemat(tsNameBase+'.mat', {'displacement': Delay})
        print 'saved data to '+tsNameBase+'.mat'
        plt.savefig(tsNameBase+'.pdf',dpi=fig_dpi)
        print 'saved plot to '+tsNameBase+'.pdf'

  except:  print 'No x/y input' ; pass

  ########### Plot Time Series with Click ##########
  def onclick(event):
    if event.button==1:
      xClick = int(event.xdata)
      yClick = int(event.ydata)
      print 'x='+str(xClick)+'+/-'+str(radius)+', y='+str(yClick)+'+/-'+str(radius)
      xmin=xClick-radius;  xmax=xClick+radius+1;
      ymin=yClick-radius;  ymax=yClick+radius+1;
      try:
         fig
         rectSelect=patches.Rectangle((xmin,ymin),radius*2+1,radius*2+1,fill=False,lw=edgeWidth)
         ax.add_patch(rectSelect)
      except: pass

      ax2.cla()

      #plot 1st time series
      Dis=[]
      for date in dateList:  Dis.append(h5timeseries['timeseries'].get(date)[ymin:ymax,xmin:xmax])
      Dis0=array(Dis)
      dis=Dis0*unitFac
      dis=reshape(dis,(len(dateList),-1))
      dis_mean=stats.nanmean(dis,1)
      if (xmax-xmin)*(ymax-ymin)==1:  dis_std=[0]*len(dateList)
      else:                           dis_std=stats.nanstd(dis,1)
      (_, caps, _)=ax2.errorbar(dates,dis_mean,yerr=dis_std,fmt='-ko',\
                                ms=markerSize, lw=lineWidth, alpha=1, mfc=markerColor,\
                                elinewidth=edgeWidth,ecolor='black',capsize=markerSize*0.5)
      for cap in caps:  cap.set_markeredgewidth(edgeWidth)
      print dis_mean

      #plot 2nd time series
      try:
         timeSeriesFile_2
         Dis2=[]
         for date in dateList:  Dis2.append(h5timeseries_2['timeseries'].get(date)[ymin:ymax,xmin:xmax])
         dis2=array(Dis2)
         dis2=dis2*unitFac
         dis2=reshape(dis2,(len(dateList),-1))
         dis2_mean=stats.nanmean(dis2,1)
         if (xmax-xmin)*(ymax-ymin)==1:  dis2_std=[0]*len(dateList)
         else:                           dis2_std=stats.nanstd(dis2,1)
         (_, caps, _)=ax2.errorbar(dates,dis2_mean,yerr=dis2_std,fmt='^',\
                                   ms=markerSize2, lw=lineWidth, alpha=1, mfc=markerColor2,\
                                   elinewidth=edgeWidth,ecolor='black',capsize=markerSize*0.5)
         for cap in caps:  cap.set_markeredgewidth(edgeWidth)
      except:  Dis2=[]

      #axis formating
      ax2.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
      if unitFac==100:     ax2.set_ylabel('Displacement [cm]',fontsize=fontSize)
      elif unitFac==1000:  ax2.set_ylabel('Displacement [mm]',fontsize=fontSize)
      else:                ax2.set_ylabel('Displacement [m]' ,fontsize=fontSize)
      ax2.set_xlabel('Time [years]',fontsize=fontSize)
      ax2.set_title('x='+str(xClick)+'+/-'+str(radius)+', y='+str(yClick)+'+/-'+str(radius))
      #ds=datevector[0]-0.2
      #de=datevector[-1]+0.2
      #ys=int(ds)
      #ye=int(de)
      #ms=int((ds-ys)*12)+1
      #me=int((de-ye)*12)+1
      #dss=datetime.datetime(ys,ms,1,0,0)
      #dee=datetime.datetime(ye,me,1,0,0)
      #ax2.set_xlim(dss,dee)
      ax2.xaxis.set_major_locator(years)
      ax2.xaxis.set_major_formatter(yearsFmt)
      ax2.xaxis.set_minor_locator(months)
      datemin = datetime.date(int(datevector[0]),1,1)
      datemax = datetime.date(int(datevector[-1])+1,1,1)
      ax2.set_xlim(datemin, datemax)

      try:
        lbound
        hbound
        ax2.set_ylim(lbound,hbound)
      except:
        ax2.set_ylim(nanmin(dis_mean-dis_std)-0.4*abs(nanmin(dis_mean)),\
                     nanmax(dis_mean+dis_std)+0.4*abs(nanmax(dis_mean)))

      for tick in ax2.xaxis.get_major_ticks():  tick.label.set_fontsize(fontSize)
      for tick in ax2.yaxis.get_major_ticks():  tick.label.set_fontsize(fontSize)
      #fig2.autofmt_xdate()     #adjust x overlap by rorating, may enble again

      if Save_timeseries in ('yes','Yes','Y','y','YES'):
         import scipy.io as sio
         Delay={}
         Delay['displacement']=Dis0
         Delay['unit']='m'
         Delay['time']=datevector
         tsNameBase='ts_x'+str(xmin)+'_'+str(xmax-1)+'y'+str(ymin)+'_'+str(ymax-1)
         sio.savemat(tsNameBase+'.mat', {'displacement': Delay})
         print 'saved data to '+tsNameBase+'.mat'
         plt2.savefig(tsNameBase+'.pdf',dpi=fig_dpi)
         print 'saved plot to '+tsNameBase+'.pdf'

      if dispTsFig in ('yes','Yes','Y','y','YES'):  plt.show()
  try:
     cid = fig.canvas.mpl_connect('button_press_event', onclick)       # Click function is available when VelFig is shown
  except: pass

  if dispTsFig in ('yes','Yes','Y','y','YES'):  plt.show()

if __name__ == '__main__':

  main(sys.argv[1:])
