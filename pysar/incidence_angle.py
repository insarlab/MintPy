#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import getopt
import h5py


def Usage():
  print '''
***************************************************************
***************************************************************

  reads required attributes from the h5 file and generates incidence angles for each pixel
  Output is incidence_angle.h5

  Usage:
       
       incidence_angle.py -f  h5 file

  Example:

       incidence_angle.py -f velocity.h5
       incidence_angle.py -f timeseries.h5
       incidence_angle.py -f temporal_coherence.h5
       
***************************************************************
***************************************************************
'''

def main(argv):
  try:
    opts, args = getopt.getopt(argv,"f:h")

  except getopt.GetoptError:
    Usage() ; sys.exit(1)


  if opts==[]:
    Usage() ; sys.exit(1)
  for opt,arg in opts:
    if opt in ("-h","--help"):
       Usage()
       sys.exit()
    elif opt == '-f':
       File = arg

    h5file=h5py.File(File,'r')
    kk=h5file.keys()
    if 'timeseries' in kk:
      k=['timeseries']
    near_range=float(h5file[k[0]].attrs['STARTING_RANGE1'])
    dR=float(h5file[k[0]].attrs['RANGE_PIXEL_SIZE'])
    r=float(h5file[k[0]].attrs['EARTH_RADIUS'])
    H=float(h5file[k[0]].attrs['HEIGHT'])
    length=float(h5file[k[0]].attrs['FILE_LENGTH'])
    width = float(h5file[k[0]].attrs['WIDTH'])
    
    far_range=near_range+dR*width
    incidence_n=np.pi-np.arccos((r**2+near_range**2-(r+H)**2)/(2*r*near_range))
    incidence_f=np.pi-np.arccos((r**2+far_range**2-(r+H)**2)/(2*r*far_range))

    look_step=(incidence_f-incidence_n)/width
    lookx=np.arange(incidence_n,incidence_f,look_step)
    
    look_angle = np.tile(lookx,(length,1))
    look_angle =  look_angle*(180./np.pi)

    print '*********************************************'
    print ''
    print 'near incidence angle : '+ str(incidence_n*180./np.pi)
    print 'far incidence angle : '+ str(incidence_f*180./np.pi)
    print 'average incidence angle : '+ str(((incidence_f+incidence_n)/2)*180./np.pi)
    print ''
    print '*********************************************'
    print ''
    print 'writing incidence_angle.h5 ...'
    print ''
    print '*********************************************'

    h5file2 = h5py.File('incidence_angle.h5','w')
    group=h5file2.create_group('mask')
    dset = group.create_dataset('mask', data=look_angle, compression='gzip')

    for key, value in h5file[k[0]].attrs.iteritems():
          group.attrs[key] = value

    h5file.close()
    h5file2.close()

if __name__ == '__main__':

  main(sys.argv[1:])





