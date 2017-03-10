#! /usr/bin/env python
#! /usr/bin/env python
###############################################################################
# 
# Project: PySAR 
# Purpose: Python Module for InSAR Time-series Analysis
# Author: Heresh Fattahi
# Created: Nov 2013
#
###############################################################################

import h5py
import sys


def Usage():
    print '''
  ######################################

  remove_dates.py timeseries_file  dates_to_remove

  remove_dates.py timeseries.h5 '20050708,20060304,20101101'

  ######################################
    '''

def main(argv):

    try:
        tsFile=sys.argv[1]
        dates2rmv=sys.argv[2]
    except:
        Usage();sys.exit(1)
  
    h5file=h5py.File(tsFile,'r')
    k=h5file.keys()
    if not 'timeseries' in k:
        sys.exit(1)
  
    dateList = h5file['timeseries'].keys()
    
    h5modified=h5py.File('modified_'+tsFile,'w')
    group=h5modified.create_group('timeseries')
    for d in dateList:
        if not d in dates2rmv:
            dataSet=h5file['timeseries'].get(d)
            dset = group.create_dataset(d, data=dataSet, compression='gzip')
        else:
            print 'removing '+ d
  
    for key,value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value
  
    try:
        dataSet=h5file['mask'].get('mask')
        group=h5modified.create_group('mask')
        dset = group.create_dataset('mask', data=dataSet, compression='gzip')
    except:
        print 'mask not found!'



if __name__ == '__main__':
    main(sys.argv[1:])


