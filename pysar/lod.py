#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

#########################################################################################
#                                                                                       #
# The empiriocal model in this program to correct the Local Oscilator Frequency Decay   #
# of Envisat ASAR instrument was suggested by Petar Marinkovic and Yngvar Larsen, 2013. #
#                                                                                       #
#########################################################################################
# Yunjun, Jan 2017: using pysar._readfile/_writefile/_datetime


import os
import sys
import time
import datetime

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


def correct_lod_file(File, outFile=None):
    # Check Sensor Type
    print 'input file: '+File
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    platform = atr['PLATFORM']
    print 'platform: '+platform
    if not platform.lower() in ['env','envisat']:
        print 'No need to correct LOD for '+platform
        sys.exit(1)

    # Output Filename
    if not outFile:
        ext = os.path.splitext(File)[1]
        outFile = os.path.splitext(File)[0]+'_LODcor'+ext
    
    # Get LOD phase ramp from empirical model
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    range_resolution = float(atr['RANGE_PIXEL_SIZE'])

    r = np.linspace(0, width-1, width)
    R = range_resolution*r*(3.87e-7)
    Ramp = np.tile(R,[length,1])
    
    yref=int(atr['ref_y'])
    xref=int(atr['ref_x'])
    Ramp -= Ramp[yref][xref]

    # Correct LOD Ramp for Input File
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5 = h5py.File(File,'r')
        epochList = sorted(h5[k].keys())
        
        h5out = h5py.File(outFile,'w')
        group = h5out.create_group(k)

        if k in ['interferograms','wrapped']:
            print 'number of interferograms: '+str(len(epochList))
            wvl = float(atr['WAVELENGTH'])
            Ramp *= -4*np.pi/wvl
            for epoch in epochList:
                print epoch
                data = h5[k][epoch].get(epoch)[:]
                atr = h5[k][epoch].attrs
                
                dates = ptime.yyyymmdd(atr['DATE12'].split('-'))
                dates = ptime.yyyymmdd2years(dates)
                dt = date[1] - date[0]
                data -= Ramp*dt
                 
                gg = group.create_group(epoch)
                dset = gg.create_dataset(epoch, data=data, compression='gzip')
                for key, value in atr.iteritems():
                    gg.attrs[key] = value

        elif k == 'timeseries':
            print 'number of acquisitions: '+str(len(epochList))
            tbase = [float(dy)/365.25 for dy in ptime.date_list2tbase(epochList)[0]]
            for i in range(len(epochList)):
                epoch = epochList[i]
                print epoch
                data = h5[k].get(epoch)[:]
                
                data -= Ramp*tbase[i]
                
                dset = group.create_dataset(epoch, data=data, compression='gzip')
            for key, value in atr.iteritems():
                group.attrs[key] = value
        else:
            print 'No need to correct for LOD for '+k+' file'
            sys.exit(1)

        h5.close()
        h5out.close()

    else:
        data, atr = readfile.read(File)
        data -= Ramp
        writefile.write(data, atr, outFile)

    return outFile


def usage():
    print '''
*****************************************************************
  Applying an empirical model to correct the Local Oscilator Drift 
  of Envisat ASAR instrument. The empiriocal model was suggested 
  by Petar Marinkovic and Yngvar Larsen, 2013.

  Usage:
      lod.py file_radarCoord [out_name]

  Example:
      lod.py timeseries.h5
      lod.py timeseries.h5 timeseries_LODcor.h5
      lod.py Seeded_unwrapIfgram.h5

*****************************************************************
    '''


#########################################################################################
def main(argv):

    # Check Inputs
    try:     File = argv[0]
    except:  usage();  sys.exit(1)
    try:     outName = argv[2]
    except:  outName = os.path.splitext(File)[0]+'_LODcor'+os.path.splitext(File)[1]

    #print '\n***************** Correct Local Oscilator Drift *******************'    
    outFile = correct_lod_file(File, outName)
    
    print 'Done.'

#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



