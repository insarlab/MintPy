#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os, sys
import h5py
import numpy as np
from pysar.utils import readfile, writefile, datetime as ptime, utils as ut


def usage():
    print("""
usage:  perp_baseline.py  timeseries_file   date  [outfile]

Generates perpendicular baseline (in Radar Coordinate) for each pixel
  with required attributes read from the h5 file

input arguments:
  timeseries_file : string, input file name/path
  epoch   : string, date in YYMMDD/YYYYMMDD format
  outfile : string, output file name/path for 2D incidence angle 
            calculated from file in radar coord

example:
  perp_baseline.py  timeseries.h5  20101020
    """)
    return

def main(argv):
    try:
        File = argv[0]
        atr = readfile.read_attribute(File)
        epoch = argv[1]
    except:
        usage();  sys.exit(1)

    try:    outFile = argv[2]
    except: outFile = None
    
    # Calculate look angle
    pbase = ut.perp_baseline_timeseries(atr, dimension=1)

    if pbase.shape[1] == 1:
        print(pbase)
        return pbase
    
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    
    h5 = h5py.File(File, 'r')
    epochList = sorted(h5[k].keys())
    epoch = ptime.yyyymmdd(epoch)
    epoch_idx = epochList.index(epoch)
    
    pbase_y = pbase[epoch_idx,:].reshape(length,1)
    pbase_xy = np.tile(pbase_y, (1, width))
    
    if not outFile:
        outFile = 'perpBaseline_'+epoch+'.h5'

    print('writing >>> '+outFile)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 'm'
    writefile.write(pbase_xy, out_file=outFile, metadata=atr)
    return outFile

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])





