#! /usr/bin/env python

############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
#import re
import h5py
from numpy import pi,round
#import getopt

def usage():
    print '''
************************************************************************

  Usage: rewrap.py  interferograms  output

  Example:
         rewrap.py  interferograms.h5

************************************************************************
    '''

def rewrap(unw):
    rewrapped = unw - round(unw/(2*pi)) * 2*pi
    return rewrapped

def main(argv):

    try:     file=argv[0]
    except:  usage();sys.exit(1)
 
    h5file=h5py.File(file)
    try:     OutName=argv[1]
    except:  OutName='rewrapped_'+file
    h5file_rewarap=h5py.File(OutName,'w')
    gg = h5file_rewarap.create_group('interferograms')
    ifgramList = h5file['interferograms'].keys()
    for ifgram in ifgramList:
        print ifgram
        unwset=h5file['interferograms'][ifgram].get(ifgram)
        unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]
        rewrapped=rewrap(unw)
        group = gg.create_group(ifgram)
        dset = group.create_dataset(ifgram, data=rewrapped, compression='gzip')
        for key, value in h5file['interferograms'][ifgram].attrs.iteritems():
            group.attrs[key] = value
 
    try:
        gm = h5file_rewarap.create_group('mask')
        mask = h5file['mask'].get('mask')
        dset = gm.create_dataset('mask', data=mask, compression='gzip')
    except:
        print 'mask not found'


if __name__ == '__main__':
    main(sys.argv[1:])


