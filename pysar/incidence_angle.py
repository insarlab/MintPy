#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jun 2016: use read() and write() for IO
#                   add incidence_angle()


import sys
import os

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile


############################################################
def incidence_angle(atr):
    ## Read Attributes
    near_range = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    r  = float(atr['EARTH_RADIUS'])
    H  = float(atr['HEIGHT'])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    
    ## Calculation
    far_range = near_range+dR*width
    incidence_n = np.pi-np.arccos((r**2+near_range**2-(r+H)**2)/(2*r*near_range))
    incidence_f = np.pi-np.arccos((r**2+ far_range**2-(r+H)**2)/(2*r*far_range))
    
    print 'near    incidence angle : '+ str(incidence_n*180./np.pi)
    print 'far     incidence angle : '+ str(incidence_f*180./np.pi)
    print 'average incidence angle : '+ str(((incidence_f+incidence_n)/2)*180./np.pi)
    print 'writing incidence_angle.h5 ...'
    
    angle_x = np.linspace(incidence_n, incidence_f, num=width, endpoint='FALSE')
    angle_xy = np.tile(angle_x,(length,1))
    angle_xy *= 180.0/np.pi
  
    return angle_xy

############################################################
def usage():
    print '''
***************************************************************
  Generates incidence angles (in Radar Coordinate) for each pixel,
     with required attributes read from the h5 file

  Usage:
      incidence_angle.py  hdf5_file  [output_file]

  Example:
      incidence_angle.py velocity.h5
      incidence_angle.py timeseries.h5
      incidence_angle.py temporal_coherence.h5
       
***************************************************************
    '''

############################################################
def main(argv):
    try:
        File = argv[0]
        atr = readfile.read_attribute(File)
    except:
        usage();  sys.exit(1)
    
    try:    outFile = argv[1]
    except: outFile = 'incidence_angle.h5'
    
    print '\n*************** Generate Incidence Angle *****************'
    ##### Calculate look angle
    angle = incidence_angle(atr)
    
    ##### Output
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 'degree'
    writefile.write(angle, atr, outFile)

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])





