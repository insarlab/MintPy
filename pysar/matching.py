#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import sys
import os
import getopt

import numpy as np
import matplotlib.pyplot as plt
import h5py


############################################################################################################
def corners(h5V1):
    k=h5V1.keys()
    WIDTH  = int(h5V1[k[0]].attrs['WIDTH'])
    LENGTH = int(h5V1[k[0]].attrs['FILE_LENGTH'])
    West   = float(h5V1[k[0]].attrs['X_FIRST'])
    North  = float(h5V1[k[0]].attrs['Y_FIRST'])
    lon_step = float(h5V1[k[0]].attrs['X_STEP'])
    lat_step = float(h5V1[k[0]].attrs['Y_STEP'])
    South  = North+lat_step*(LENGTH-1)
    East   = West +lon_step*(WIDTH -1)
 
    lon=np.arange(West, West +WIDTH *lon_step,lon_step)
    lat=np.arange(North,North+LENGTH*lat_step,lat_step)

    return West,East,North,South,lon,lat,WIDTH,LENGTH


#################################
def nearest_neighbor(x,y, tbase, pbase):
    """ find nearest neighbour """
    dist = np.sqrt((tbase -x)**2+(pbase -y)**2)
    indx=dist==min(dist)
    return indx

##################################
def nearest(x, X):
    """ find nearest neighbour """
    dist = np.sqrt((X -x)**2)
    indx=np.where(dist==min(dist))
  
    return indx[0]


############################################################################################################
def Usage():
    print'''
    ***************************************************************
    To match two geocoded velocity maps with a common area.
    Function automatically finds the common area and calculates 
    the average offset between the two velocity.
    
    Usage:
       matching.py Velocity1.h5  Velocity2.h5
    
    ***************************************************************  
    '''

############################################################################################################
def main(argv):

    ## default value
    fig_size   = [16.0,16.0]
    manual_matching = 'no'

    try: 
        V1file=sys.argv[1]
        V2file=sys.argv[2]
    except:  Usage();sys.exit(1)

    try:     manual_matching = sys.argv[3]
    except:  manual_matching = 'no' 

    ####################################################
    h5V1=h5py.File(V1file,'r')
    h5V2=h5py.File(V2file,'r')

    k=h5V1.keys()
    V1set = h5V1[k[0]].get(k[0])
    V2set = h5V2[k[0]].get(k[0])

    V1=V1set[0:V1set.shape[0],0:V1set.shape[1]]
    V2=V2set[0:V2set.shape[0],0:V2set.shape[1]]

    ####################################################
    West1,East1,North1,South1,lon1,lat1,WIDTH1,LENGTH1 = corners(h5V1)
    West2,East2,North2,South2,lon2,lat2,WIDTH2,LENGTH2 = corners(h5V2)
    print 'finding the corners of the whole area'
    West  = min(West1, West2)
    East  = max(East1, East2)
    North = max(North1,North2)
    South = min(South1,South2)

    # width and length of the whole area
    lon_step = float(h5V1[k[0]].attrs['X_STEP'])
    lat_step = float(h5V1[k[0]].attrs['Y_STEP'])

    WIDTH  = int(round((East  - West )/lon_step + 1.0))
    LENGTH = int(round((South - North)/lat_step + 1.0))

    ####################################################
    LON = np.arange(West, West +WIDTH *lon_step,lon_step) 
    LAT = np.arange(North,North+LENGTH*lat_step,lat_step)

    indx1 = nearest(West1,  LON)[0]
    indy1 = nearest(North1, LAT)[0]
    indx2 = nearest(West2,  LON)[0]
    indy2 = nearest(North2, LAT)[0]

    ####################################################
    VV1=np.zeros([LENGTH,WIDTH]);   VV1[:,:]=np.nan
    VV2=np.zeros([LENGTH,WIDTH]);   VV2[:,:]=np.nan

    #import pdb; pdb.set_trace()
    VV=np.zeros([LENGTH,WIDTH])
    VV[:,:]=np.nan
    VV1[indy1:indy1+LENGTH1,indx1:indx1+WIDTH1] = V1   
    VV2[indy2:indy2+LENGTH2,indx2:indx2+WIDTH2] = V2

    M=VV2-VV1
    offset=np.nansum(M) / np.sum(np.isfinite(M))  

    if np.isnan(offset):
        print '**************************************************'
        print 'WARNING:'
        print ''
        print 'No common area found between two velocity maps'
        print 'At least one common pixel is required.'
        print 'No matching applied. Continue without matching ...'
        print '**************************************************'
 
    fig = plt.figure(figsize=fig_size)
    if manual_matching=='yes':
        ax=fig.add_subplot(111)
        ax.imshow(V1)
        xc=[] 
        yc=[] 
        print 'please click on start and end point of the desired profile'
        def onclick(event):
            if event.button==1:
                print 'click'
                xc.append(int(event.xdata))
                yc.append(int(event.ydata))
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()    
        x0=xc[0];x1=xc[1]
        y0=yc[0];y1=yc[1]
  
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.imshow(V2)
        xc=[]
        yc=[]
        print 'please click on start and end point of the desired profile'
        def onclick(event):
            if event.button==1:
                print 'click'
                xc.append(int(event.xdata))
                yc.append(int(event.ydata))
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        x00=xc[0];x11=xc[1]
        y00=yc[0];y11=yc[1]
  
        offset = np.nansum(V2[y00:y11,x00:x11]) / np.sum(np.isfinite(V2[y00:y11,x00:x11]))\
                 - np.nansum(V1[y0:y1,x0:x1]) / np.sum(np.isfinite(V1[y0:y1,x0:x1]))
  
      #  offset=V2[y00:y11,x00:x11]-V2[y0:y1,x0:x1]
      
    if np.isnan(offset):
        print '**************************************************'
        print 'WARNING:'
        print ''
        print 'No common area found between two velocity maps'
        print 'At least one common pixel is required.'
        print 'No matching applied. Continue without matching ...'
        print '**************************************************'   
    else:
        print 'Average offset between two velocity in the common area is: ' + str(offset)
        V2=V2-offset
 
    indv2=np.isfinite(V2)
    VV[indy1:indy1+LENGTH1,indx1:indx1+WIDTH1] = V1
    VV[indy2:indy2+LENGTH2,indx2:indx2+WIDTH2][indv2] = V2[indv2]
    
    ####################################
    outName='Matched.h5'
    print 'writing '+outName
    h5velocity = h5py.File(outName,'w')
    group=h5velocity.create_group('velocity')
    dset = group.create_dataset(os.path.basename('velocity'), data=VV, compression='gzip')
    
    for key , value in h5V1[k[0]].attrs.iteritems():
        group.attrs[key]=value
 
    group.attrs['WIDTH']       = WIDTH
    group.attrs['FILE_LENGTH'] = LENGTH
    group.attrs['X_FIRST']     = West
    group.attrs['Y_FIRST']     = North
    group.attrs['X_STEP']      = lon_step
    group.attrs['Y_STEP']      = lat_step
  
    h5velocity.close()
    h5V1.close()
    h5V2.close()
 
    ####################################
    

    fig=plt.subplot(2,2,1)
    plt.imshow(VV1)
    plt.colorbar()
 
    fig=plt.subplot(2,2,2)
    plt.imshow(VV2)
    plt.colorbar()
 
    fig=plt.subplot(2,2,3)
    plt.imshow(VV)
    plt.colorbar()
    
    fig=plt.subplot(2,2,4)
    plt.imshow(M)
    plt.colorbar()
 
    plt.savefig('match.png',bbox_inches='tight')
 
    plt.show()


############################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



