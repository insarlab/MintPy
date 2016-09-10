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
    cir_par = circle_par.split(',')
    #import pdb; pdb.set_trace()
    try:
        c_y    = int(cir_par[0])
        c_x    = int(cir_par[1])
        radius = int(float(cir_par[2]))
    except:
        try:
            c_lat  = float(cir_par[0])
            c_lon  = float(cir_par[1])
            radius = int(float(cir_par[2]))
            c_y = subset.coord_geo2radar(c_lat,atr,'lat')
            c_x = subset.coord_geo2radar(c_lon,atr,'lon')
        except:
            print '\nERROR: Unrecognized circle index format: '+circle_par
            print 'Supported format:'
            print '--circle 200,300,20            for radar coord input'
            print '--circle 31.0214,130.5699,20   for geo   coord input\n'
            return 0

    y,x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx



#################################  Usage  ####################################
def Usage():
    print '''
******************************************************************************

  Calculate Spatial average/mean of multi-temporal 2D datasets.
      write the optimal reference date to reference_date.txt, and
      write turbulent date to drop_date.txt

  Usage:
      mean_spatial.py    multi_temporal_file    mask_file
      mean_spatial.py -f multi_temporal_file -m mask_file [--circle x,y,rad -x sub -y sub]

      -f : multi-temporal file

      ## Mark area included
      -m : mask file
      -x : subset in x/column/azimuth/longitude direction
      -y : subset in y/row   /range  /latitude  direction

      ## Mark area excluded
      --circle : circle shape subset in y/x or lat/lon, [center_y/lat,center_x/lon,radius;center_y/lat2,...]
                 i.e.: --circle 361,567,50;610,536,25          [int   coord input for radar coordinate]
                       --circle 31.0212,131.0569,30            [float coord input for geo   coordinate]

  Example:
      mean_spatial.py sum_Seeded_ts.h5 mask.h5
      mean_spatial.py -f sum_Seeded_ts.h5 -m Mask.h5 -x 200:800 -y 230:700 --circle 361,567,50;610,536,25
      mean_spatial.py -f sum_Seeded_ts.h5 -m Mask.h5 --circle 31.0212,131.0569,30

******************************************************************************
    '''


#############################  Main Function  ################################
def main(argv):

    ##### Default
    fontSize    = 12
    lineWidth   = 2
    markerColor = 'crimson'
    markerSize  = 16

    disp_fig  = 'no'
    save_fig  = 'yes'
    save_list = 'yes'

    ref_file  = 'reference_date.txt'
    drop_file = 'drop_date.txt'

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
            elif opt == '--circle'   :  cir_par   = [i for i in arg.split(';')]
            #elif opt == '-o':  outName   = arg
            
    else:
        try:  File = argv[0]
        except: Usage(); sys.exit(1)
        try:  maskFile = argv[1]
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
        Mask_orig,Matr = readfile.read(maskFile)
        print 'mask file: '+maskFile
        Masking = 'yes'
    except:
        print 'No mask. Use the whole area for ramp estimation.'
        Masking = 'no'
        Mask_orig=np.ones((length,width))
    Mask = np.zeros((length,width))
    Mask[:] = Mask_orig[:]

    ## Bounding Subset
    try:
        xsub
        ysub
        ysub,xsub = subset.check_subset_range(ysub,xsub,atr)
        Mask[ysub[0]:ysub[1],xsub[0]:xsub[1]] = Mask_orig[ysub[0]:ysub[1],xsub[0]:xsub[1]]*2
        #Mask[0:ysub[0],:]      = 0
        #Mask[ysub[1]:length,:] = 0
        #Mask[:,0:xsub[0]]      = 0
        #Mask[:,xsub[1]:width]  = 0
    except:
        Mask = Mask_orig*2
        print 'No subset input.'

    ## Circle Inputs
    try:
        cir_par
        for i in range(len(cir_par)):
            cir_idx = circle_index(atr,cir_par[i])
            Mask[cir_idx] = Mask_orig[cir_idx]
            print 'Circle '+str(i)+': '+cir_par[i]
    except: print 'No circle of interest input.'
    
    ## Mask output
    idx = Mask == 2
    idxNum = float(sum(sum(idx)))
    
    fig = plt.figure()
    plt.imshow(Mask,cmap='gray')
    plt.savefig(outNameBase+'_mask.png',bbox_inches='tight')
    print 'save mask to '+outNameBase+'_mask.png'
    #fig.clf()

    ##### Calculation
    meanList   = np.zeros(epochNum)
    pixPercent = np.zeros(epochNum)
    pixT = 0.7
    print 'calculating ...'
    print '  Date       Mean   Percentage'
    for i in range(epochNum):
        epoch = epochList[i]
        d      = h5file[k].get(epoch)[:]
        #d[Mask==0]  = np.nan
        
        meanList[i]   = np.nanmean(d[idx])
        pixPercent[i] = np.sum(d[idx] >= pixT)/idxNum
        
        print epoch+' :   %.2f    %.1f%%'%(meanList[i],pixPercent[i]*100)
    del d
    h5file.close()

    ##### Reference date - Max Value
    top3 = sorted(zip(meanList,epochList), reverse=True)[:3]
    print '------------ Top 3 Mean ------------------'
    print top3
    ## Write to txt file
    fref = open(ref_file,'w')
    fref.write(str(top3[0][1])+'\n')
    fref.close()
    print 'write optimal reference date to '+ref_file
    idxMean = meanList == np.nanmax(meanList)

    ##### Drop dates - mean threshold
    #meanT = 0.7
    #idxMean  = meanList < meanT
    #print '------------ Mean Value < '+str(meanT)+' --------'
    #print np.array(epochList)[idxMean]
    #print meanList[idxMean]

    ##### Drop dates - good pixel percentage
    pixNumT = 0.7
    print '------------ Good Pixel Percentage < %.0f%% -------'%(pixNumT*100)
    idxPix = pixPercent < pixNumT
    dropEpochList = np.array(epochList)[idxPix]
    print dropEpochList
    print pixPercent[idxPix]
    ## Write to txt file
    fdrop = open(drop_file,'w')
    for i in range(len(dropEpochList)):
        fdrop.write(str(dropEpochList[i])+'\n')
    fdrop.close()
    print 'write drop dates to '+drop_file
    print '-------------------------------------------'

    ##### Display
    fig = plt.figure(figsize=(12,12))
    ax  = fig.add_subplot(211)
    ax.plot(dates, meanList, '-ko', ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    #ax.plot([dates[0],dates[-1]],[meanT,meanT], '--b', lw=lineWidth)
    #sc = ax.scatter(dates, np.tile(0.5,epochNum), c=meanList, s=22**2, alpha=0.3, vmin=0.0, vmax=1.0)
    #ax.scatter(np.array(dates)[idxMean], 0.5, c=meanList[idxMean], s=22**2, alpha=1.0, vmin=0.0, vmax=1.0)
    ax = ptime.adjust_xaxis_date(ax,datevector)
    ax.set_ylim(0,1)
    ax.set_title('Spatial Average Value', fontsize=fontSize)
    ax.set_xlabel('Time [years]',         fontsize=fontSize)
    #cbar = plt.colorbar(sc)
    #cbar.set_label('Spatial Mean of Normalized Sum Epochs')

    ax  = fig.add_subplot(212)
    ax.plot(dates, pixPercent, '-ko', ms=markerSize, lw=lineWidth, alpha=0.7, mfc=markerColor)
    ax.plot([dates[0],dates[-1]],[pixNumT,pixNumT], '--b', lw=lineWidth)
    ax = ptime.adjust_xaxis_date(ax,datevector)
    ax.set_ylim(0,1)
    ax.set_title('Percenrage of Pixels with Value > '+str(pixNumT), fontsize=fontSize)
    ax.set_xlabel('Time [years]',         fontsize=fontSize)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:3.0f}%'.format(i*100) for i in vals])

    if save_fig == 'yes':
        plt.savefig(outNameBase+'.png',bbox_inches='tight')
        print 'save figure to '+outNameBase+'.png'

    if disp_fig == 'yes':
        plt.show()

    ##### Output
    if save_list == 'yes':
        epochList6 = ptime.yymmdd(epochList)
        fl = open(outNameBase+'.txt','w')
        for i in range(epochNum):
            str_line = epochList6[i]+'    %.2f    %.2f\n'%(meanList[i],pixPercent[i])
            fl.write(str_line)
        fl.close()
        print 'write data to '+outNameBase+'.txt\n'


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
