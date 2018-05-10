#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import sys
import os
import argparse

from numpy import pi
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


##############################################################################
def usage():
    print '''usage: save_roipac.py  file  [date_info]

Convert PySAR hdf5 file to ROI_PAC format 

argument:
  file : file to be converted.
         for velocity  : the ouput will be a one year interferogram.
         for timeseries: if date is not specified, the last date will be used
                         if two dates are specified, the earlier date will be
                             used as the reference date.

example:
    '''
    return

############################################################
EXAMPLE='''example:
  save_roipac.py  velocity.h5
  save_roipac.py  timeseries.h5    20050601
  save_roipac.py  timeseries.h5    050601    --ref-date 040728
  save_roipac.py  unwrapIfgram.h5  filt_091225-100723-sim_HDR_8rlks_c10.unw
  save_roipac.py  unwrapIfgram.h5  091225-100723
  save_roipac.py  temporal_coherence.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Convert PySAR HDF5 file to ROI_PAC format.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='HDF5 file to be converted.\n'+\
                        'for velocity  : the ouput will be a one year interferogram.\n'+\
                        'for timeseries: if date is not specified, the last date will be used.')
    parser.add_argument('epoch', nargs='?', help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o','--output', dest='outfile', help='output file name.')
    parser.add_argument('-r','--ref-date', dest='ref_date', help='Reference date for timeseries file')

    inps = parser.parse_args()
    return inps


##############################################################################
def main(argv):
    inps = cmdLineParse()
  
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    atr['PROCESSOR'] = 'roipac'
    atr['INSAR_PROCESSOR'] = 'roipac'

    h5file = h5py.File(inps.file,'r')
  
    if k == 'velocity':
        dset = h5file['velocity'].get('velocity')
        data = dset[0:dset.shape[0],0:dset.shape[1]]
        print "converting velocity to a 1 year interferogram."
        wvl=float(h5file[k].attrs['WAVELENGTH'])
        data=(-4*pi/wvl)*data

        inps.outfile=inps.file.split('.')[0]+'.unw'
        print 'writing >>> '+inps.outfile
        writefile.write(data,atr,inps.outfile)

    elif k in multi_dataset_hdf5_file:
        dateList = sorted(h5file[k].keys())
        try:
            inps.epoch = [date for date in dateList if inps.epoch in date][0]
        except:
            print 'No input date specified >>> continue with the last date'
            inps.epoch = dateList[-1]
        if k in ['timeseries']:
            inps.epoch = ptime.yyyymmdd(inps.epoch)

        ## Data
        print 'reading %s and %s ...' % (inps.ref_date, inps.epoch)
        data = h5file[k].get(inps.epoch)[:]
        if inps.ref_date:
            inps.ref_date = ptime.yyyymmdd(inps.ref_date)
            data -= h5file[k].get(inps.ref_date)[:]

        ## Attributes
        if k in ['timeseries']:
            wvl = float(atr['WAVELENGTH'])
            data *= -4*pi/wvl
            atr['FILE_TYPE']             = '.unw'
            atr['P_BASELINE_TIMESERIES'] = '0.0'
            atr['UNIT']                  = 'radian'
        if inps.ref_date:
            atr['DATE']              = inps.ref_date[2:8]
            atr['DATE12']            = '%s-%s' % (inps.ref_date[2:8],inps.epoch[2:8])

        ## Writing
        if not inps.outfile:
            if k in ['timeseries']:
                inps.outfile = '%s_%s.unw' % (inps.ref_date[2:8],inps.epoch[2:8])
            else:
                inps.outfile = '%s.cor' % (inps.epoch)
        print 'writing >>> '+inps.outfile
        writefile.write(data,atr,inps.outfile)

    elif k in ['interferograms','coherence','wrapped']:
        ## Check input
        igramList = sorted(h5file[k].keys())
        try:
            inps.epoch = [igram for igram in igramList if inps.epoch in igram][0]
        except:
            print 'No input interferogram specified >>> continue with the last one'
            inps.epoch = igramList[-1]

        ## Read and Write
        print 'reading '+inps.epoch+' ... '
        atr = dict(h5file[k][inps.epoch].attrs)
        data = h5file[k][inps.epoch].get(inps.epoch)[:]
        if k == 'interferograms':
            try:
                ref_y = int(atr['ref_y'])
                ref_x = int(atr['ref_x'])
                data -= data[ref_y,ref_x]
                print 'consider the reference pixel in y/x: %d/%d' % (ref_y, ref_x)
            except:
                print 'No ref_y/x info found in attributes.'
        atr['PROCESSOR'] = 'roipac'
        atr['INSAR_PROCESSOR'] = 'roipac'

        inps.outfile = inps.epoch
        print 'writing >>> '+ inps.outfile
        writefile.write(data, atr, inps.outfile)  

    else:
        data = h5file[k].get(k)[:]
        if not inps.outfile:
            if k in ['temporal_coherence']:
                inps.outfile=inps.file.split('.')[0]+'.cor'
            elif k in ['dem','.hgt','.dem']:
                atr['FILE_TYPE'] = '.dem'
                inps.outfile=os.path.splitext(inps.file)[0]+'.dem'
            else:
                inps.outfile=inps.file.split('.')[0]+'.unw'
        print 'writing >>> '+ inps.outfile
        writefile.write(data,atr,inps.outfile)


    h5file.close()
    return


##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
