###############################################################################
# Program is part of PySAR      
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun
# Author:  Heresh Fattahi, Zhang Yunjun 
###############################################################################
#  deramp are modified from a software originally
#  written by Scott Baker with the following licence:
#
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
#     import pysar.utils.deramp as deramp
#


import os
import time
import h5py
import numpy as np
from pysar.utils import datetime as ptime, readfile, writefile
from pysar.objects import timeseries, ifgramStack


##################################################################
def remove_data_surface(data, mask, surf_type='plane'):
    '''Remove surface from input data matrix based on pixel marked by mask'''
    mask[np.isnan(data)] = 0
    mask = mask.flatten(1) 
    z = data.flatten(1)
    ndx= mask !=0
    x = list(range(0,np.shape(data)[1]))
    y = list(range(0,np.shape(data)[0]))
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
       length= float(h5file['velocity'].attrs['LENGTH'])
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
    ## 1st mask
    print('removing 1st surface ...')
    i = 0
    mask_i = np.zeros(data.shape,data.dtype)
    mask_i[ysub[2*i]:ysub[2*i+1],:] = mask[ysub[2*i]:ysub[2*i+1],:]

    dataOut_i,ramp_i = remove_data_surface(data,mask_i,surf_type)
    dataOut[ysub[2*i]:ysub[2*i+1],:] = dataOut_i[ysub[2*i]:ysub[2*i+1],:]

    ## 2 - last masks
    for i in range(1,surfaceNum):
        print('removing '+str(i+1)+'th surface ...')
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
    k = atr['FILE_TYPE']
    print('Input file is '+k)
    print('remove ramp type: '+surf_type)

    if not outFile:
        outFile = os.path.splitext(File)[0]+'_'+surf_type+os.path.splitext(File)[1]

    if maskFile:
        mask = readfile.read(maskFile, datasetName='mask')[0]
        print('read mask file: '+maskFile)
    else:
        mask = np.ones((int(atr['LENGTH']), int(atr['WIDTH'])))
        print('use mask of the whole area')

    if k == 'timeseries':
        obj = timeseries(File)
        data = obj.read()
        numDate = data.shape[0]
        print('estimating ramp for each acquisition ...')
        prog_bar = ptime.progressBar(maxValue=numDate)
        for i in range(numDate):
            if not ysub:
                data[i,:,:] = remove_data_surface(np.squeeze(data[i,:,:]), mask, surf_type)[0]
            else:
                data[i,:,:] = remove_data_multiple_surface(np.squeeze(data), mask, surf_type, ysub)
            prog_bar.update(i+1, suffix=obj.dateList[i])
        prog_bar.close()
        objOut = timeseries(outFile)
        objOut.write2hdf5(data=data, refFile=File)
  
    elif k in ['interferograms','wrapped','coherence']:
        print('number of interferograms: '+str(len(epochList)))
        date12_list = ptime.list_ifgram2date12(epochList)

        if k == 'interferograms':
            mask_bk = np.zeros(mask.shape)
            mask_bk = mask
            print('do not consider zero value pixel for interferograms')

        for i in range(epoch_num):
            epoch = epochList[i]
            data = h5file[k][epoch].get(epoch)[:]
            if k == 'interferograms':
                mask = mask_bk
                mask[data == 0.] = 0

            if not ysub:
                data_n,ramp = remove_data_surface(data,mask,surf_type)
            else:
                data_n = remove_data_multiple_surface(data, mask, surf_type, ysub)

            gg   = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data_n)
            for key,value in h5file[k][epoch].attrs.items():
                gg.attrs[key] = value
            prog_bar.update(i+1, suffix=date12_list[i])

    ## Single Dataset File
    else:
        data,atr = readfile.read(File)
        print('Removing '+surf_type+' from '+k)
        data_n,ramp = remove_data_surface(data, mask, surf_type)
        print('writing >>> '+outFile)
        writefile.write(data_n,atr,outFile)
  
    print('Remove '+surf_type+' took ' + str(time.time()-start) +' secs')
    return outFile

