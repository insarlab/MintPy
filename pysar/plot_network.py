#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Dec 2015: Add support for coherence/wrapped, update display
# Yunjun, Jun 2016: Add plot_network(), plot_bperp_hist(),
#                   axis_adjust_date_length(), igram_pairs()


import sys
import os
import getopt
import time
import datetime

import h5py
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from pylab import *
from numpy import *

import pysar._pysar_utilities as ut
import pysar._datetime as ptime
import pysar._network  as pnet
import pysar._readfile as readfile

  
######################################
def usage():
    print '''
******************************************************************************
  Ploting the network of interferograms and 
      the baseline history of SAR acquisitions.

  Usage:
      plot_network.py -f interferogramsFile -s fontsize -w linewidth
      plot_network.py -b baselineFile -l pairsFile

      -f : interferograms file stored in hdf5 file format, 
               supported files: interferograms, coherence and wrapped
      -b : baseline file
      -l : pairs info list file
      -s : the font size used for x and y labels (default is 12)
      -w : line width used for plotting (default is 2)
      -m : marker size (default is 16)
      -c : marker face color (default is orange)
      -t : temporal threshold
      -d : date (all interferograms with master or slave using the specified date is removed) 
  
      Save and Output
      --save      : save                    figure
      --nodisplay : save and do not display figure
      --list      : output pairs list file named Pairs.list

  Example:
      plot_network.py LoadedData.h5
      plot_network.py -f Coherence.h5           --nodisplay
      plot_network.py -f Modified_LoadedData.h5 --nodisplay --list
      plot_network.py -f LoadedData.h5 -l pairs.list --save
      plot_network.py -b bl_list.txt   -l Pairs.list

******************************************************************************  
    '''


##########################  Main Function  ##############################
def main(argv):

    ## Global Variables
    global fontSize, lineWidth, markerColor, markerSize

    ## Default Values
    fontSize    = 12
    lineWidth   = 2
    markerColor = 'orange'
    markerSize  = 16
    saveFig  = 'no'
    dispFig  = 'yes'
    saveList = 'no'

    if len(sys.argv)>2:
        try:  opts, args = getopt.getopt(argv,"h:b:f:s:w:l:m:c:o:",['save','nodisplay','list'])
        except getopt.GetoptError:   usage() ; sys.exit(1)

        for opt,arg in opts:
            if opt in ("-h","--help"):  usage();  sys.exit()
            elif opt == '-b':        baselineFile= arg
            elif opt == '-f':        igramsFile  = arg
            elif opt == '-l':        listFile    = arg
            elif opt == '-s':        fontSize    = int(arg)
            elif opt == '-w':        lineWidth   = int(arg)
            elif opt == '-m':        markerSize  = int(arg)
            elif opt == '-c':        markerColor = arg
            #elif opt == '-o':        figName2  = arg;   saveFig = 'yes'
            elif opt == '--save'     : saveFig   = 'yes'
            elif opt == '--nodisplay': dispFig   = 'no';  saveFig = 'yes'
            elif opt == '--list'     : saveList  = 'yes'

        try:  igramsFile
        except:
            try:  baselineFile
            except:  usage() ; sys.exit(1)

    elif len(sys.argv)==2:   igramsFile = argv[0]
    else:                    usage() ; sys.exit(1)

    ##### Output figure name
    figName1 = 'BperpHist.pdf'
    figName2 = 'Network.pdf'
    try:
        igramsFile
        if 'Modified' in igramsFile:
            figName1 = 'BperpHist_Modified.pdf'
            figName2 = 'Network_Modified.pdf'
    except: pass

    ############# Read Time and Bperp ################ 
    print '\n******************** Plot Network **********************'
    try:
        igramsFile
        atr = readfile.read_attribute(igramsFile)
        k = atr['FILE_TYPE']
        if k not in  ['interferograms','coherence','wrapped']:
            print 'Only interferograms / coherence / wrapped are supported.';  sys.exit(1)

        print 'reading date and perpendicular baseline from '+k
        dateList  = ptime.date_list(igramsFile)
        dateList6 = ptime.yymmdd(dateList)
        print 'number of acquisitions: '+str(len(dateList))
        Bp = ut.Baseline_timeseries(igramsFile)
    except:
        baselineFile
        print 'reading date and perpendicular baseline from '+baselineFile
        dateList, Bp = pnet.read_baseline_file(baselineFile)[0:2]
        dateList6 = ptime.yymmdd(dateList)

    ############# Read Pairs Info ####################
    print 'reading pairs info'
    try:
        listFile
        pairs = pnet.read_pairs_list(listFile,dateList)
    except:
        pairs = pnet.read_igram_pairs(igramsFile)
    print 'number of pairs       : '+str(len(pairs))

    if saveList == 'yes':
        pnet.write_pairs_list(pairs,dateList,'Pairs.list')
        print 'save pairs info to Pairs.list'

    ############# Read Unwrapping Error Info #######################
    ## For simulated interferograms only
    ## To plot the interferograms with unwrapping errors with a different color
    #N_unw_err=0
    #try:
    #  for ifgram in  ifgramList:
    #    if h5file[k[0]][ifgram].attrs['unwrap_error']=='yes':
    #       N_unw_err=N_unw_err+1
    #except: pass
    #    
    #if N_unw_err>0:
    #   pairs_ue=np.zeros([N_unw_err,2],np.int)
    #   i=0
    #   for ifgram in  ifgramList:
    #     if h5file[k[0]][ifgram].attrs['unwrap_error']=='yes':
    #       date1,date2 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
    #       pairs_ue[i][0]=dateList6.index(date1)
    #       pairs_ue[i][1]=dateList6.index(date2)
    #       i=i+1
    #

    ############### Fig 1 - Interferogram Network ##################
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    fig1 = pnet.plot_network(ax1, pairs, dateList, Bp)

    if saveFig=='yes':
        plt.savefig(figName2,bbox_inches='tight')
        print 'save figure to '+figName2

    ############## Fig 2 - Baseline History ###################  
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    fig2 = pnet.plot_perp_baseline_hist(ax2,dateList,Bp)

    if saveFig=='yes':
        plt.savefig(figName1,bbox_inches='tight')
        print 'save figure to '+figName1

    if dispFig == 'yes':  plt.show() 

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



