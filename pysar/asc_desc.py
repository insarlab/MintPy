#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib import colors
import getopt
import h5py


def usage():
    print'''
***************************************************************
  Projecting Asc and Desc LOS velocities to Horizontal and Vertical 
  components. Horizontal component is parallel to the azimuth angle.


  Usage: asc_des.py  V1.h5 V2.h5 azimuth incidence1 incidence2

  Example: 
      asc_des.py seeded_T134_masked.h5 seeded_T256.h5 16  23 38
      asc_des.py seeded_T134_masked.h5 seeded_T256.h5 16

***************************************************************  
    '''


def corners(h5V1):
    k=h5V1.keys()
    WIDTH=int(h5V1[k[0]].attrs['WIDTH'])
    LENGTH=int(h5V1[k[0]].attrs['FILE_LENGTH'])
    West=float(h5V1[k[0]].attrs['X_FIRST'])
    North=float(h5V1[k[0]].attrs['Y_FIRST'])
    lon_step=float(h5V1[k[0]].attrs['X_STEP'])
    lat_step=float(h5V1[k[0]].attrs['Y_STEP'])
    South=North+lat_step*(LENGTH-1)
    East=West+lon_step*(WIDTH-1)
 
    lon=np.arange(West,West+WIDTH*lon_step,lon_step)
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

#################################

def find_row_column(Lon,Lat,h5file):
    ################################################
    # finding row and column numbers of the GPS point
  
    lat,lon,lat_step,lon_step = get_lat_lon(h5file)
    idx= nearest(Lon, lon, lon_step)
    idy= nearest(Lat, lat, lat_step)
    if idx !=[] and idy != []:
        IDX=np.where(idx==True)[0][0]
        IDY=np.where(idy==True)[0][0]
    else:
        IDX=np.nan
        IDY=np.nan
    return IDY, IDX

###############################################

def get_lat_lon(h5file):

    k=h5file.keys()
 
    if 'interferograms' in k:
 
        ifgramList = h5file['interferograms'].keys()
        Width=float(h5file['interferograms'][ifgramList[0]].attrs['WIDTH'])
        Length= float(h5file['interferograms'][ifgramList[0]].attrs['FILE_LENGTH'])
        ullon=float(h5file['interferograms'][ifgramList[0]].attrs['X_FIRST'])
        ullat=float(h5file['interferograms'][ifgramList[0]].attrs['Y_FIRST'])
        lon_step=float(h5file['interferograms'][ifgramList[0]].attrs['X_STEP'])
        lat_step=float(h5file['interferograms'][ifgramList[0]].attrs['Y_STEP'])
        lon_unit=h5file['interferograms'][ifgramList[0]].attrs['Y_UNIT']
        lat_unit=h5file['interferograms'][ifgramList[0]].attrs['X_UNIT']
 
    elif 'timeseries' in k:
        Width=float(h5file['timeseries'].attrs['WIDTH'])
        Length= float(h5file['timeseries'].attrs['FILE_LENGTH'])
        ullon=float(h5file['timeseries'].attrs['X_FIRST'])
        ullat=float(h5file['timeseries'].attrs['Y_FIRST'])
        lon_step=float(h5file['timeseries'].attrs['X_STEP'])
        lat_step=float(h5file['timeseries'].attrs['Y_STEP'])
        lon_unit=h5file['timeseries'].attrs['Y_UNIT']
        lat_unit=h5file['timeseries'].attrs['X_UNIT']
 
    elif 'velocity' in k:
        Width=float(h5file['velocity'].attrs['WIDTH'])
        Length= float(h5file['velocity'].attrs['FILE_LENGTH'])
        ullon=float(h5file['velocity'].attrs['X_FIRST'])
        ullat=float(h5file['velocity'].attrs['Y_FIRST'])
        lon_step=float(h5file['velocity'].attrs['X_STEP'])
        lat_step=float(h5file['velocity'].attrs['Y_STEP'])
        lon_unit=h5file['velocity'].attrs['Y_UNIT']
        lat_unit=h5file['velocity'].attrs['X_UNIT']
 
    lllat=ullat+Length*lat_step
    urlon=ullon+Width*lon_step
    lat=np.arange(ullat,lllat,lat_step)
    lon=np.arange(ullon,urlon,lon_step)
    return lat,lon,lat_step,lon_step





def main(argv):

    try: 
        V1file=sys.argv[1]
        V2file=sys.argv[2]
    except:
        usage();sys.exit(1)
 
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
    ####################################################
    print 'finding the corners of the whole area'
    West = max(West1,West2)
    East = min(East1,East2)
    North = min(North1,North2)
    South = max(South1,South2)
    # width and length of the whole area
    lon_step=float(h5V1[k[0]].attrs['X_STEP'])
    lat_step=float(h5V1[k[0]].attrs['Y_STEP'])
 
    WIDTH  = int(round((East  - West )/lon_step + 1.0))
    LENGTH = int(round((South - North)/lat_step + 1.0))

    ####################################################
    indx11=nearest(West, lon1)
    indy11=nearest(North, lat1)
    indx12=nearest(East, lon1)
    indy12=nearest(South, lat1)
 
    indx21=nearest(West, lon2)
    indy21=nearest(North, lat2)
    indx22=nearest(East, lon2)
    indy22=nearest(South, lat2)


    ####################################################

    VV1=V1[indy11:indy12+1,indx11:indx12+1]
    VV2=V2[indy21:indy22+1,indx21:indx22+1]
 
    
    data=np.zeros((2,WIDTH*LENGTH))
    data[0,:]=VV1.flatten(0)
    data[1,:]=VV2.flatten(0)
 
    # theta1=np.pi*41./180.
    # theta2=np.pi*23./180.
 
    heading1 = float(h5V1['velocity'].attrs['HEADING'])
    if heading1 < 0:
        heading1=heading1+360.
    heading1=heading1*np.pi/180.
 
    heading2 = float(h5V2['velocity'].attrs['HEADING'])
    if heading2 < 0:
        heading2=heading2+360.
    heading2=heading2*np.pi/180.

    #   heading1=np.pi*346.8/180.
    #   heading2=np.pi*193.4/180.
    azimuth=float(sys.argv[3])
    print 'azimuth = '+str(azimuth)
    azimuth=np.pi*azimuth/180.
 
    try:     theta1=float(sys.argv[4])
    except:  theta1=float(h5V1['velocity'].attrs['LOOK_REF2'])
    print 'Look angle 1: '+str(theta1)
 
    try:     theta2=float(sys.argv[5])
    except:  theta2=float(h5V2['velocity'].attrs['LOOK_REF2'])
    print 'Look angle 2: '+str(theta2)  
    theta1=theta1*np.pi/180.
    theta2=theta2*np.pi/180.
    # azimuth=np.pi*29.5/180.

    A=np.zeros((2,2));
    A[0,0]=np.cos(theta1);A[0,1]=np.sin(theta1)*np.sin(heading1-azimuth);
    A[1,0]=np.cos(theta2);A[1,1]=np.sin(theta2)*np.sin(heading2-azimuth);
 
    A1=np.linalg.pinv(A)
    Luh=np.dot(A1,data)
    
    ######
    outName='Up.h5'
    print 'writing '+outName
    h5velocity = h5py.File(outName,'w')
    group=h5velocity.create_group('velocity')
    dset = group.create_dataset('velocity', data=np.reshape(Luh[0,:],(LENGTH,WIDTH)), compression='gzip')
    
    for key , value in h5V1[k[0]].attrs.iteritems():
        group.attrs[key]=value
 
    group.attrs['WIDTH']=WIDTH
    group.attrs['FILE_LENGTH']=LENGTH
    group.attrs['X_FIRST']=West
    group.attrs['Y_FIRST']=North
    group.attrs['X_STEP']=lon_step
    group.attrs['Y_STEP']=lat_step
  
    h5velocity.close()

    ########
    outName='Hz.h5'
    print 'writing '+outName
    h5velocity = h5py.File(outName,'w')
    group=h5velocity.create_group('velocity')
    dset = group.create_dataset('velocity', data=np.reshape(Luh[1,:],(LENGTH,WIDTH)), compression='gzip')
 
    for key , value in h5V1[k[0]].attrs.iteritems():
        group.attrs[key]=value

    group.attrs['WIDTH']=WIDTH
    group.attrs['FILE_LENGTH']=LENGTH
    group.attrs['X_FIRST']=West
    group.attrs['Y_FIRST']=North
    group.attrs['X_STEP']=lon_step
    group.attrs['Y_STEP']=lat_step
 
    h5velocity.close()
 
    h5V1.close()
    h5V2.close()


####################################
if __name__ == '__main__':
    main(sys.argv[1:])



