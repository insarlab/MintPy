#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
import numpy as np
import time 

def remove_surface(data,Mask,surf_type):
      Mask[np.isnan(data)]=0
      Mask=Mask.flatten(1) 
      z = data.flatten(1)
      ndx= Mask !=0
      x = range(0,np.shape(data)[1])
      y = range(0,np.shape(data)[0])
      x1,y1 = np.meshgrid(x,y)
      points = np.vstack((y1.flatten(1),x1.flatten(1))).T
      if surf_type=='quadratic':
           G = np.array([points[:,0]**2,points[:,1]**2,points[:,0],points[:,1],points[:,0]*points[:,1],np.ones(np.shape(points)[0])],np.float32).T
      elif surf_type=='plane':
           G = np.array([points[:,0],points[:,1],np.ones(np.shape(points)[0])],np.float32).T
      elif surf_type == 'quadratic_range':
           G = np.array([points[:,1]**2,points[:,1],np.ones(np.shape(points)[0])],np.float32).T
      elif surf_type == 'quadratic_azimuth':
           G = np.array([points[:,0]**2,points[:,0],np.ones(np.shape(points)[0])],np.float32).T
      elif surf_type=='plane_range':
           G = np.array([points[:,1],np.ones(np.shape(points)[0])],np.float32).T
      elif surf_type=='plane_azimuth':
           G = np.array([points[:,0],np.ones(np.shape(points)[0])],np.float32).T

      z = z[ndx]
      G = G[ndx]
      G1=np.linalg.pinv(G)
      plane = np.dot(G1,z)

      if surf_type=='quadratic':
           zplane=plane[0]*y1**2 + plane[1]*x1**2 + plane[2]*y1 + plane[3]*x1 + plane[4]*y1*x1 + plane[5]
      elif surf_type=='plane':
           zplane= plane[0]*y1 + plane[1]*x1 + plane[2]
      elif surf_type == 'quadratic_range':
           zplane= plane[0]*x1**2  + plane[1]*x1 + plane[2]
      elif surf_type == 'quadratic_azimuth':
           zplane= plane[0]*y1**2  + plane[1]*y1 + plane[2]
      elif surf_type == 'plane_range':
           zplane=  plane[0]*x1 + plane[1]
      elif surf_type == 'plane_azimuth':
           zplane= plane[0]*y1 + plane[1]

      data_n = data - zplane
      data_n[data == 0.] = 0.
      data_n = np.array(data_n,np.float32) 
      return data_n

##################################################################
def remove_surface_timeseries(surf_type,h5file,h5flat,Mask):
#  Mask=Mask.flatten(1)
  start = time.time()
  ifgramList = h5file['timeseries'].keys()
  group = h5flat.create_group('timeseries')
  for key,value in h5file['timeseries'].attrs.iteritems():
         group.attrs[key] = value

  for ifgram in ifgramList:
        print "Removing " + surf_type  +" from " + ifgram
        dset1 = h5file['timeseries'].get(ifgram)
        data = dset1[0:dset1.shape[0],0:dset1.shape[1]]
        data_n = remove_surface(data,Mask,surf_type) 
        data_n = np.array(data_n,np.float32)
        # referencing to ref point
        dset = group.create_dataset(ifgram, data=data_n, compression='gzip')

  for key,value in h5file['timeseries'].attrs.iteritems():
         group.attrs[key] = value
  print 'Remove Plane took ' + str(time.time()-start) +' secs'


def remove_surface_velocity(surf_type,h5file,h5flat,Mask):
#  Mask=Mask.flatten(1)
  start = time.time()
  dset1 = h5file['velocity'].get('velocity')
  data = dset1[0:dset1.shape[0],0:dset1.shape[1]]
  print "Removing " + surf_type  +" from velocity"
  data_n = remove_surface(data,Mask,surf_type) 
  data_n = np.array(data_n,np.float32)
  group = h5flat.create_group('velocity')

  dset = group.create_dataset('velocity', data=data_n, compression='gzip')
  for key,value in h5file['velocity'].attrs.iteritems():
         group.attrs[key] = value

  print 'Remove Plane took ' + str(time.time()-start) +' secs'


def remove_surface_igrams(surf_type,h5file,h5flat,Mask):
  start = time.time()
  ifgramList = h5file['interferograms'].keys()
  gg = h5flat.create_group('interferograms')
#  Mask=Mask.flatten(1)  
  for ifgram in ifgramList:
      group = gg.create_group(ifgram)
      for key,value in h5file['interferograms'][ifgram].attrs.iteritems():
           group.attrs[key] = value
      print "Removing " + surf_type  +" from " + ifgram
      dset1 = h5file['interferograms'][ifgram].get(ifgram)
      data = dset1[0:dset1.shape[0],0:dset1.shape[1]]

      data_n = remove_surface(data,Mask,surf_type)
      data_n = np.array(data_n,np.float32)
      dset = group.create_dataset(ifgram, data=data_n, compression='gzip')


  print 'Remove Plane took ' + str(time.time()-start) +' secs'

##################################################################


