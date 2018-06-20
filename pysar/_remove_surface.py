#! /usr/bin/env python2
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################
# remove_plane are modified from a software originally
# written by Scott Baker with the following licence:
###############################################################################
#  Copyright (c) 2011, Scott Baker 
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
############################################################################### 
# Recommend usage:
#     import pysar._remove_surface as rm


import os
import time

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile


##################################################################
def remove_data_surface(data, mask, surf_type='plane'):
    '''Remove surface from input data matrix based on pixel marked by mask'''
    mask[np.isnan(data)] = 0
    mask = mask.flatten(1) 
    z = data.flatten(1)
    ndx= mask !=0
    x = range(0,np.shape(data)[1])
    y = range(0,np.shape(data)[0])
    x1,y1 = np.meshgrid(x,y)
    points = np.vstack((y1.flatten(1),x1.flatten(1))).T
    if surf_type=='quadratic':
        G = np.array([points[:,0]**2,points[:,1]**2,points[:,0],points[:,1],points[:,0]*points[:,1],\
                     np.ones(np.shape(points)[0])],np.float32).T
    elif surf_type=='plane':
        G = np.array([points[:,0],points[:,1],\
                     np.ones(np.shape(points)[0])],np.float32).T
    elif surf_type == 'quadratic_range':
        G = np.array([points[:,1]**2,points[:,1],\
                     np.ones(np.shape(points)[0])],np.float32).T
    elif surf_type == 'quadratic_azimuth':
        G = np.array([points[:,0]**2,points[:,0],\
                     np.ones(np.shape(points)[0])],np.float32).T
    elif surf_type=='plane_range':
        G = np.array([points[:,1],\
                     np.ones(np.shape(points)[0])],np.float32).T
    elif surf_type=='plane_azimuth':
        G = np.array([points[:,0],\
                     np.ones(np.shape(points)[0])],np.float32).T

    z = z[ndx]
    G = G[ndx]
    G1=np.linalg.pinv(G)
    plane = np.dot(G1,z)

    if   surf_type == 'quadratic':
        zplane = plane[0]*y1**2 + plane[1]*x1**2 + plane[2]*y1 + plane[3]*x1 + plane[4]*y1*x1 + plane[5]
    elif surf_type =='plane':
        zplane = plane[0]*y1 + plane[1]*x1 + plane[2]
    elif surf_type == 'quadratic_range':
        zplane = plane[0]*x1**2  + plane[1]*x1 + plane[2]
    elif surf_type == 'quadratic_azimuth':
        zplane = plane[0]*y1**2  + plane[1]*y1 + plane[2]
    elif surf_type == 'plane_range':
        zplane = plane[0]*x1 + plane[1]
    elif surf_type == 'plane_azimuth':
        zplane = plane[0]*y1 + plane[1]

    '''
    ## Some notes from _pysar_utilities.py remove_surface_velocity()
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print 'Plane parameters:'
    if surf_type == 'plane_range':
       print 'range gradient = ' + str(1000*plane[0][0]) + ' mm/yr/pixel'
       width= float(h5file['velocity'].attrs['WIDTH'])
       MaxRamp=width*1000*plane[0][0]
       print 'Maximum ramp in range direction = ' + str(MaxRamp) + ' mm/yr'
       h5flat['velocity'].attrs['Range_Gradient'] = str(1000*plane[0][0]) + '   mm/yr/pixel'
       h5flat['velocity'].attrs['Range_Ramp'] = str(MaxRamp) + '   mm/yr'
    elif surf_type == 'plane_azimuth':
       print 'azimuth gradient = ' + str(1000*plane[0][0]) + ' mm/yr/pixel'
       length= float(h5file['velocity'].attrs['FILE_LENGTH'])
       MaxRamp=length*1000*plane[0][0]
       h5flat['velocity'].attrs['Azimuth_Gradient'] = str(1000*plane[0][0]) + '   mm/yr/pixel'
       h5flat['velocity'].attrs['Azimuth_Ramp'] = str(MaxRamp) +'   mm/yr'
       print 'Maximum ramp in azimuth direction = '+ str(MaxRamp) + ' mm/yr'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    '''

    data_n = data - zplane
    data_n[data == 0.] = 0.
    data_n = np.array(data_n, data.dtype)
    zplane = np.array(zplane, data.dtype)

    return data_n, zplane


##################################################################
def remove_data_multiple_surface(data, mask, surf_type, ysub):
    ## ysub = [0,2400,2000,6800]
    dataOut = np.zeros(data.shape,data.dtype)
    dataOut[:] = np.nan

    surfaceNum = len(ysub)/2
    ## 1st Mask
    print 'removing 1st surface ...'
    i = 0
    mask_i = np.zeros(data.shape,data.dtype)
    mask_i[ysub[2*i]:ysub[2*i+1],:] = mask[ysub[2*i]:ysub[2*i+1],:]

    dataOut_i,ramp_i = remove_data_surface(data,mask_i,surf_type)
    dataOut[ysub[2*i]:ysub[2*i+1],:] = dataOut_i[ysub[2*i]:ysub[2*i+1],:]

    ## 2 - last Masks
    for i in range(1,surfaceNum):
        print 'removing '+str(i+1)+'th surface ...'
        mask_i = np.zeros(data.shape,data.dtype)
        mask_i[ysub[2*i]:ysub[2*i+1],:] = mask[ysub[2*i]:ysub[2*i+1],:]

        dataOut_i,ramp_i = remove_data_surface(data,mask_i,surf_type)

        if ysub[2*i] < ysub[2*i-1]:
            dataOut[ysub[2*i]:ysub[2*i-1],:]  += dataOut_i[ysub[2*i]:ysub[2*i-1],:]
            dataOut[ysub[2*i]:ysub[2*i-1],:]  /= 2
            dataOut[ysub[2*i-1]:ysub[2*i+1],:] = dataOut_i[ysub[2*i-1]:ysub[2*i+1],:]
        else:
            dataOut[ysub[2*i]:ysub[2*i+1],:]   = dataOut_i[ysub[2*i]:ysub[2*i+1],:]

    return dataOut


##################################################################
def remove_surface(File, surf_type, maskFile=None, outFile=None, ysub=None):
    start = time.time()
    atr = readfile.read_attribute(File)
    
    # Output File Name
    if not outFile:
        outFile = os.path.splitext(File)[0]+'_'+surf_type+os.path.splitext(File)[1]
    
    if maskFile:
        Mask = readfile.read(maskFile, epoch='mask')[0]
        print 'read mask file: '+maskFile
    else:
        Mask = np.ones((int(atr['FILE_LENGTH']), int(atr['WIDTH'])))
        print 'use mask of the whole area'
    
    ##### Input File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    print 'Input file is '+k
    print 'remove ramp type: '+surf_type
    
    ## Multiple Datasets File
    if k in ['interferograms','coherence','wrapped','timeseries']:
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epoch_num = len(epochList)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

        h5flat = h5py.File(outFile,'w')
        group  = h5flat.create_group(k)
        print 'writing >>> '+outFile

    if k in ['timeseries']:
        print 'number of acquisitions: '+str(len(epochList))
        for i in range(epoch_num):
            epoch = epochList[i]
            data = h5file[k].get(epoch)[:]
            
            if not ysub:
                data_n,ramp = remove_data_surface(data, Mask, surf_type) 
            else:
                data_n = remove_data_multiple_surface(data, Mask, surf_type, ysub)
  
            dset = group.create_dataset(epoch, data=data_n, compression='gzip')
            prog_bar.update(i+1, suffix=epoch)
        for key,value in h5file[k].attrs.iteritems():
            group.attrs[key] = value
  
    elif k in ['interferograms','wrapped','coherence']:
        print 'number of interferograms: '+str(len(epochList))
        date12_list = ptime.list_ifgram2date12(epochList)

        if k == 'interferograms':
            mask_bk = np.zeros(Mask.shape)
            mask_bk = Mask
            print 'do not consider zero value pixel for interferograms'

        for i in range(epoch_num):
            epoch = epochList[i]
            data = h5file[k][epoch].get(epoch)[:]
            if k == 'interferograms':
                Mask = mask_bk
                Mask[data == 0.] = 0

            if not ysub:
                data_n,ramp = remove_data_surface(data,Mask,surf_type)
            else:
                data_n = remove_data_multiple_surface(data, Mask, surf_type, ysub)

            gg   = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data_n, compression='gzip')
            for key,value in h5file[k][epoch].attrs.iteritems():
                gg.attrs[key] = value
            prog_bar.update(i+1, suffix=date12_list[i])

    ## Single Dataset File
    else:
        data,atr = readfile.read(File)
        print 'Removing '+surf_type+' from '+k
  
        if not ysub:
            data_n,ramp = remove_data_surface(data, Mask, surf_type)
        else:
            data_n = remove_data_multiple_surface(data, Mask, surf_type, ysub)

        print 'writing >>> '+outFile
        writefile.write(data_n,atr,outFile)
  
    try:
        h5file.close()
        h5flat.close()
        prog_bar.close()
    except: pass
  
    print 'Remove '+surf_type+' took ' + str(time.time()-start) +' secs'
    return outFile

