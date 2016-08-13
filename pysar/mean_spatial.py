#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# 

import os
import sys
import getopt

import h5py
import numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt

import pysar._readfile as readfile
import pysar._datetime as ptime
import pysar.subset as subset


##############################################################################
def circle_index(atr,circle_par):
    ## Return Index of Elements within a Circle
    ## Inputs:
    ##     atr        : (dictionary) attributes containging width, length info
    ##     circle_par : (string) center_x,center_y,radius

    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    c_x,c_y,radius = [int(i) for i in circle_par.split(',')]
    y,x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx


#################################  Usage  ####################################
def Usage():
    print '''
******************************************************************************

  Calculate Spatial average/mean of multi-temporal 2D datasets.

  Usage:
      mean_spatial.py    multi_temporal_file    mask_file
      mean_spatial.py -f multi_temporal_file -m mask_file --circle x,y,rad -x sub -y sub

  Example:
      mean_spatial.py sum_ts.h5 mask.h5

******************************************************************************
    '''


#############################  Main Function  ################################
def main(argv):

    ##### Default
    fontSize    = 12
    lineWidth   = 2
    markerColor = 'crimson'
    markerSize  = 16

    ##### Check Inputs
    if len(sys.argv)>3:
        try:
            opts, args = getopt.getopt(argv,'h:f:m:o:x:y:',['help','circle='])
        except getopt.GetoptError:
            print 'Error in reading input options!';  Usage() ; sys.exit(1)

        for opt,arg in opts:
            if opt in ("-h","--help"):    Usage() ; sys.exit()
            elif opt == '-f':  File      = arg
            elif opt == '-m':  maskFile  = arg
            elif opt == '-x':  xsub = [int(i) for i in arg.split(':')];  xsub.sort()
            elif opt == '-y':  ysub = [int(i) for i in arg.split(':')];  ysub.sort()
            elif opt == '--circle' :  cir_par = [i for i in arg.split(';')]
            #elif opt == '-o':  outName   = arg
            
    else:
        File = argv[1]
        try:  maskFile = argv[2]
        except: pass

    try:  atr  = readfile.read_attributes(File)
    except: Usage(); sys.exit(1)
    ext      = os.path.splitext(File)[1].lower()
    FileBase = os.path.basename(File).split(ext)[0]
    outNameBase = 'spatialMean_'+FileBase
    print '\n*************** Spatial Average ******************'

    ##### Input File Info
    k = atr['FILE_TYPE']
    print 'Input file is '+k
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])

    h5file = h5py.File(File)
    epochList = h5file[k].keys();
    epochList = sorted(epochList)
    epochNum  = len(epochList)
    print 'number of epoch: '+str(epochNum)
    dates,datevector = ptime.date_list2vector(epochList)

    ##### Mask Info
    try:
        Mask,Matr = readfile.read(maskFile)
        print 'mask file: '+maskFile
        Masking = 'yes'
    except:
        print 'No mask. Use the whole area for ramp estimation.'
        Masking = 'no'
        Mask=np.ones((length,width))

    ## Bounding Subset
    try:
        xsub, ysub
        ysub,xsub = subset.check_subset_range(ysub,xsub,atr)
        Mask[0:ysub[0],:]      = 0
        Mask[ysub[1]:length,:] = 0
        Mask[:,0:xsub[0]]      = 0
        Mask[:,xsub[1]:width]  = 0
    except: print 'No subset input.'

    ## Circle Inputs
    try:
        cir_par
        for i in range(len(cir_par)):
            cir_idx = circle_index(atr,cir_par[i])
            Mask[cir_idx] = 0
            print 'Circle '+str(i)+': '+cir_par[i]
    except: print 'No circle of interest input.'
    fig = plt.figure()
    plt.imshow(Mask,vmin=0,vmax=1,cmap='gray')
    plt.savefig(outNameBase+'_mask.png',bbox_inches='tight')
    print 'save mask to '+outNameBase+'_mask.png'

    ##### Calculation
    meanList = np.zeros(epochNum)
    print 'calculating ...'
    for i in range(epochNum):
        epoch = epochList[i]
        d      = h5file[k].get(epoch)[:]

        d[Mask==0]  = np.nan
        dMean       = np.nanmean(d)
        print epoch+' : '+str(dMean)
        meanList[i] = dMean
    del d
    h5file.close()

    ##### Max Value
    top3 = sorted(zip(meanList,epochList), reverse=True)[:3]
    print '------------ Top 3 Dates ------------------'
    print top3
    
    meanT = 0.7
    idxT  = meanList < meanT
    print '------------ Below Threshold - '+str(meanT)+' --------'
    print np.array(epochList)[idxT]
    print meanList[idxT]
    print '-------------------------------------------'

    ##### Display
    fig = plt.figure(2,figsize=(12,6))
    ax  = fig.add_subplot(111)
    ax.plot(dates, meanList, '-ko', ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    ax.plot([dates[0],dates[-1]],[meanT,meanT], '--b', lw=lineWidth)
    ax = ptime.adjust_xaxis_date(ax,datevector)
    ax.set_ylim(0,1)
    ax.set_title('Normalized Sum Epochs', fontsize=fontSize)
    ax.set_xlabel('Time [years]',         fontsize=fontSize)
    ax.set_ylabel('Normalized Sum Epochs',fontsize=fontSize)
    plt.savefig(outNameBase+'.png',bbox_inches='tight')
    print 'save figure to '+outNameBase+'.png'

    ##### Output
    epochList6 = ptime.yymmdd(epochList)
    fl = open(outNameBase+'.txt','w')
    for i in range(epochNum):
        str_line = epochList6[i]+'    '+str(meanList[i])+'\n'
        fl.write(str_line)
    fl.close()
    print 'write data to '+outNameBase+'.txt\n'


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
