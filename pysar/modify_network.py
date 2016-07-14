#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: Update Coherence file
# Yunjun, Oct 2015: Add 'T' option for template file input
#                   and pysar.dropIfgIndex in template content


import sys
import os
import getopt
import time
import datetime

import matplotlib
from numpy import *
from pylab import *
import h5py

import pysar._pysar_utilities as ut
import pysar._datetime as ptime
import pysar._readfile as readfile


###########################  Sub Function  #############################
######################################
def nearest_neighbor(x,y, tbase, pbase):
  """ find nearest neighbour """
  dist = sqrt((tbase -x)**2+(pbase -y)**2)
  indx=dist==min(dist)
  return indx

##############  Usage  ###############
def Usage():
  print '''
  ******************************************
  ******************************************
  Modify the network of interferograms and/or coherence.

  usage:
       modify_network.py -f interferogramsFile -s fontsize -w linewidth -T templateFile

    -f : interferograms file stored in hdf5 file format
    -C : coherence file stored in hdf5 file format

    Remove Pairs:
    -t : temporal threshold
    -b : baseline threshold
    -d : date (all interferograms that includes date is removed)
    -l : list of interferograms to remove  
    -N : interferogram numbers to remove (1 as the first)
    -n : (yes or no)network display to manually choose igrams from the network to remove(default is yes)
    -T : template file with pysar.dropIfgIndex setted (recommend)
         Example for template option:
             pysar.drop.ifgIndex   = 7:9,15,25,26,31,35,39,48,53,62,66,67,72,73,77,85,86,90,102,107,111
             pysar.drop.date       = 20080102          # not implemented yet

    -r : list file with interferograms to remove

    -s : the font size used for x and y labels (default is 12)
    -w : line width used for plotting (default is 2)
    -m : marker size (default is 16)
    -c : marker face color (default is orange) 
    
  Example:

       modify_network.py LoadedData.h5
       modify_network.py -f Seeded_LoadedData.h5 -T $TE/KirishimaT246EnvD2.template
       modify_network.py -f LoadedData.h5 -b 400 -n no
       modify_network.py -f LoadedData.h5 -b 500 -t 600 -d 080307
       modify_network.py -f LoadedData.h5 -d '080307 091023'           
       modify_network.py -f LoadedData.h5 -l 'filt_080307-091023-sim_HDR_8rlks_c10.unw filt_080307-090814-sim_HDR_8rlks_c10.unw' 
       modify_network.py -f LoadedData.h5 -N '1 4 20 76 89 100'
       modify_network.py -f LoadedData.h5 -C Coherence.h5 -N '1 4 20 76 89 100'

  ******************************************
  ******************************************  
  '''


#########################  Main Function  ##############################
def main(argv):

  lineWidth   = 2
  fontSize    = 12
  markerColor = 'orange'
  markerSize  = 16
  networkDisplay = 'no'

  if len(sys.argv)>2:

    try:
      opts, args = getopt.getopt(argv,"h:f:C:s:w:m:c:t:b:d:l:n:N:T:l:")
    except getopt.GetoptError:
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage();  sys.exit()
      elif opt == '-f':        file           = arg
      elif opt == '-C':        corFile        = arg
      elif opt == '-s':        fontSize       = int(arg)
      elif opt == '-w':        lineWidth      = int(arg)
      elif opt == '-m':        markerSize     = int(arg)
      elif opt == '-c':        markerColor    = arg
      elif opt == '-t':        temp_thr       = float(arg)
      elif opt == '-b':        base_thr       = float(arg)
      elif opt == '-d':        dates2Rmv      = arg
      elif opt == '-l':        ifgrams_to_rmv = arg
      elif opt == '-n':        networkDisplay = arg
      elif opt == '-N':        ifgrams_Number_to_rmv = arg.split()
      elif opt == '-T':        templateFile   = arg
      elif opt == '-l':        list_fileList  = arg.split(',')


    try:  file
    except:  Usage() ; sys.exit(1)

  elif len(sys.argv)==2:
    file = argv[0]
    networkDisplay = 'yes'
  else:   Usage() ; sys.exit(1)

  ## display network for modification, if no other limit setted
  try:
    temp_thr
    base_trh
    dates2Rmv
    ifgrams_to_rmv
    ifgrams_Number_to_rmv
    networkDisplay = 'yes'
  except: pass

###########################################################
  h5file = h5py.File(file)
  k=h5file.keys()
  if 'interferograms' in k: k[0] = 'interferograms'
  elif 'coherence'    in k: k[0] = 'coherence'
  print 'Input file is '+k[0]
  #if h5file.keys()[0] != 'interferograms':
  #    print 'Input file should be interferograms'; sys.exit(1)
  ifgramList=h5file[k[0]].keys()

  try:     ifgrams_to_rmv
  except:  ifgrams_to_rmv=[]

###########################################################
  ##### L - list file of interferograms
  try:
    ifgrams_to_rmv = list_file.read()
  except: pass 

  #####  T - templateFile, pysar.dropIfgIndex
  try:
    templateFile
    template = readfile.read_template(templateFile)
    drop_ifg_index = template['pysar.drop.ifgIndex'].split(',')
    print 'drop interferogram index:'
    print drop_ifg_index
    try:    ifgrams_Number_to_rmv
    except: ifgrams_Number_to_rmv = []
    for index in drop_ifg_index:
       index_temp = [int(i) for i in index.split(':')];    index_temp.sort()
       if   len(index_temp)==2:
           for j in range(index_temp[0],index_temp[1]+1):  ifgrams_Number_to_rmv.append(str(j))
       elif len(index_temp)==1:                            ifgrams_Number_to_rmv.append(index)
       else: print 'Unrecoganized input: '+index
  except: pass

  #####  N - interferogram number list
  try:
    for i in ifgrams_Number_to_rmv:
       print i+'    '+ifgramList[int(i)-1]
       ifgrams_to_rmv.append(ifgramList[int(i)-1])
  except: pass

  #####  b - perpendicular baseline limit
  try:
    base_thr
    print 'interferograms with the spatial baseline longer than '+ str(base_thr)+' m is removed'
    for ifgram in  ifgramList:
       Baseline = (float(h5file[k[0]][ifgram].attrs['P_BASELINE_BOTTOM_HDR'])+\
                   float(h5file[k[0]][ifgram].attrs['P_BASELINE_TOP_HDR']))/2
       if abs(Baseline) > base_thr:
         if not ifgram in ifgrams_to_rmv:   ifgrams_to_rmv.append(ifgram)
  except:    print 'No Spatial Baseline threshold applied'

  ##### d - dates to remove
  try:
    dates2Rmv
    print 'interferograms with any of following dates will be removed: '+ dates2Rmv
    for ifgram in  ifgramList:
      date1,date2 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
      if (date1 in dates2Rmv) or (date2 in dates2Rmv):
         if not ifgram in ifgrams_to_rmv:   ifgrams_to_rmv.append(ifgram)
  except:   print 'No specific dates selected to remove'

  ##### t - temporal baseline limit
  tbase,dateList,dateDict,dateList6=ut.date_list(h5file)
  try:
    temp_thr
    print 'Applying the temporal baseline threshold with threshold of '+str(temp_thr)+' days'
    for ifgram in  ifgramList:
       date1,date2 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')      
       ind1 = dateList6.index(date1)
       ind2 = dateList6.index(date2)
       dt=tbase[ind2]-tbase[ind1]
       if dt>temp_thr:
          if not ifgram in ifgrams_to_rmv:
            ifgrams_to_rmv.append(ifgram)
  except:
    print 'No Temporal Baseline threshold applied'

############################################################
############################################################
  if networkDisplay=='yes':

    dates,datevector = ptime.date_list2vector(dateList)

    ##################################################  
    Bp = ut.Baseline_timeseries(file)
    #############################################################
 
    ifgramList = h5file[k[0]].keys()
    igram_pairs=np.zeros([len(ifgramList),2],np.int)
    i=0
    for ifgram in  ifgramList:
      date1,date2 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
      igram_pairs[i][0]=dateList6.index(date1)
      igram_pairs[i][1]=dateList6.index(date2)
      i=i+1

    ############################################################
    import matplotlib.pyplot as plt
    fig1 = plt.figure(1)
    ax1=fig1.add_subplot(111)

    ax1.cla()
    # ax1.plot(dates,Bp, 'o',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    print tbase
    ax1.plot(tbase,Bp, 'o',ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    for ni in range(len(ifgramList)):
      ax1.plot(array([tbase[igram_pairs[ni][0]],tbase[igram_pairs[ni][1]]]),\
               array([Bp[igram_pairs[ni][0]],Bp[igram_pairs[ni][1]]]),'k',lw=4) 
    # ax1.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
    ax1.set_ylabel('Bperp [m]',fontsize=fontSize)
    ax1.set_xlabel('Time [years]',fontsize=fontSize)
    ts=datevector[0]+0.2
    te=datevector[-1]+0.2
    ys=int(ts)
    ye=int(te)
    ms=int((ts-ys)*12)
    me=int((te-ye)*12)
    if ms>12:       ys =ys+1;       ms=1
    if me>12:       ye =ye+1;       me=1
    if ms<1:        ys =ys-1;       ms=12
    if me<1:        ye =ye-1;       me=12

    dss=datetime.datetime(ys,ms,1,0,0)
    dee=datetime.datetime(ye,me,1,0,0)
    ax1.set_ylim(min(Bp)-0.4*abs(min(Bp)),max(Bp)+0.4*max(Bp))

    xticklabels = getp(gca(), 'xticklabels')
    yticklabels = getp(gca(), 'yticklabels')
    setp(yticklabels, 'color', 'k', fontsize=fontSize)
    setp(xticklabels, 'color', 'k', fontsize=fontSize)

    ##########################################  
    x=[]
    y=[]
    Master_index_torremove=[]
    Slave_index_torremove=[]
    a_tbase=array(tbase)
    a_Bp=array(Bp)
    def onclick(event):
      if event.button==1:
        print 'click'
        xClick = event.xdata
        yClick = event.ydata
        idx=nearest_neighbor(xClick,yClick, a_tbase, a_Bp)       
        xr = a_tbase[idx]
        yr = a_Bp[idx]
        ix=tbase.index(xr)+1
        print ix
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
    print dateList6
    numIgrams_rmv=np.shape(R)[1]
    for ifgram in  ifgramList:
       date1,date2 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
       for i in range(numIgrams_rmv):
           if dateList6[R[0][i]]==date1 and dateList6[R[1][i]]==date2:
               ifgrams_to_rmv.append(ifgram)

  else:
    print 'No network display.'
############################################################
############################################################

  print 'Number of interferograms to remove: '+str(len(ifgrams_to_rmv))
  print 'List   of interferograms to remove:' 
  print ifgrams_to_rmv
  file_modified='Modified_'+file
  h5filem = h5py.File(file_modified,'w')
  gg = h5filem.create_group(k[0])
  ifgram=ifgramList[0]
  unw = h5file[k[0]][ifgram].get(ifgram)
  MaskZero=np.ones([unw.shape[0],unw.shape[1]])

  print 'writing >>> modified interferogram file ...'
  print 'Number of interferograms: '+str(len(ifgramList)-len(ifgrams_to_rmv))
  for ifgram in  ifgramList:
     if not ifgram in ifgrams_to_rmv:
        print ifgram
        unwSet = h5file[k[0]][ifgram].get(ifgram)
        unw = unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]        
        MaskZero=unw*MaskZero
        group = gg.create_group(ifgram)
        dset = group.create_dataset(ifgram, data=unw, compression='gzip')
        for key, value in h5file[k[0]][ifgram].attrs.iteritems():
           group.attrs[key] = value

  Mask=np.ones([unwSet.shape[0],unwSet.shape[1]])
  Mask[MaskZero==0]=0
  atrMask = h5file[k[0]][ifgram].attrs

  gm = h5filem.create_group('mask')
  dset = gm.create_dataset('mask', data=Mask, compression='gzip')
  for key, value in atrMask.iteritems():
      gm.attrs[key] = value

  h5file.close()
  h5filem.close()


  ####################### Coherence ########################

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
     print 'writing >>> modified coherence file ...'
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
     h5fileCor.close()
     h5fileCorm.close()
  except:
     print 'No coherence file to be updated.'

  ############################################################

  print 'writing >>> Modified_Mask.h5'
  
  h5mask = h5py.File('Modified_Mask.h5','w')
  group=h5mask.create_group('mask')
  dset = group.create_dataset(os.path.basename('mask'), data=Mask, compression='gzip')
  for key, value in atrMask.iteritems():
      group.attrs[key] = value
  h5mask.close()      


########################################################################
if __name__ == '__main__':
  main(sys.argv[1:])



