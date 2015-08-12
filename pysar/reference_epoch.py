#! /usr/bin/env python

############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import h5py


def Usage():

  print '''
  ***************************************************************

  Referencing all time-series epochs to the specified date.

  reference_epoch.py   time-series  ref_Date Outname

  Example:

    reference_epoch.py timeseries.h5 20050107 timeseries_ref.h5

    reference_epoch.py timeseries.h5 20050107

  ***************************************************************

'''

def main(argv):
  try:
    timeSeriesFile = argv[0]
    refDate = argv[1]
  except:
    Usage() ; sys.exit(1)

  try:
    outName=argv[2]
  except:
    outName='Ref_changed_'+timeSeriesFile


  h5t=h5py.File(timeSeriesFile)
  dateList=h5t['timeseries'].keys()
  if not refDate in dateList:
     print '''**********************
     Error:  Reference date was not found.
             Choose a date available in the time-series.
             Exit without any action.
     '''
     sys.exit(1)
  
  refDataSet=h5t['timeseries'].get(refDate)
  refData=refDataSet[0:refDataSet.shape[0],0:refDataSet.shape[1]]
  print 'referencing all epochs to ' + refDate

  h5t2=h5py.File(outName,'w')
  group = h5t2.create_group('timeseries')
  for d in dateList:
     print d
     ds=h5t['timeseries'].get(d)
     data=ds[0:ds.shape[0],0:ds.shape[1]]
     dset = group.create_dataset(d, data=data-refData, compression='gzip')

  for key,value in h5t['timeseries'].attrs.iteritems():
      group.attrs[key] = value
  
  group.attrs['ref_date']=refDate
  try:
      h5t['mask'].get('mask')
      dset1 = h5t['mask'].get('mask')
      Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
      group=h5t2.create_group('mask')
      dset = group.create_dataset('mask', data=Mask, compression='gzip')
  except:
      print 'no mask in the file.'


if __name__ == '__main__':

  main(sys.argv[1:])


