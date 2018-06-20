#! /usr/bin/env python2
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import os
import sys 

import h5py
import numpy as np
import matplotlib.pyplot as plt

import pysar._readfile as readfile


################################################################################
def usage():
    print '''usage: correlation_with_dem.py  dem_file  file  [y0:y1 x0:x1]

Calculates the correlation of the DEM and file

example:
  correlation_with_dem.py   radar_8rlks.hgt  velocity_masked.h5
    '''
    return


################################################################################
def main(argv):
    try:
        demFile = argv[0]
        File    = argv[1]
    except:
        usage(); sys.exit(1)

    dem, demRsc = readfile.read(demFile)
    data, atr   = readfile.read(File)
    print 'Input file is '+atr['FILE_TYPE']

    # Subset
    try:
        y0,y1 = [int(i) for i in argv[2].split(':')]
        x0,x1 = [int(i) for i in argv[3].split(':')]
        data = data[y0:y1, x0:x1]
        dem = dem[y0:y1, x0:x1]
    except: pass

    # Calculation
    dem = dem.flatten(1)
    data = data.flatten(1)
    ndx = ~np.isnan(data)
    C1 = np.zeros([2,len(dem[ndx])])
    C1[0][:] = dem[ndx]
    C1[1][:] = data[ndx]

    # Display
    print '-------------------------------------------'
    print 'Correlation with the DEM:  %.2f' % np.corrcoef(C1)[0][1]
    print '-------------------------------------------'
    print 'DEM info:'
    print '    Max height difference: %.2f m' % (np.max(dem[ndx])-np.min(dem[ndx]))
    print '    Average        height: %.2f m' % np.mean(dem[ndx])
    print '    Height            Std: %.2f m' % np.std(dem[ndx])
    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:]) 
