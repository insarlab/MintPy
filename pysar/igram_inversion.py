#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import sys
import os
import getopt

import h5py

import pysar._pysar_utilities as ut
import pysar._readfile as readfile


######################################
def usage():
    print '''
********************************************************************************
  Inversion of interferograms using L1 or L2 norm minimization (Default is L2)
  
  Usage:
      igram_inversion.py interferograms_file
      igram_inversion.py -f interferograms_file [ -l method -o timeseries_file]
  
      -f: stacked interferograms file
      -l: inverse method, L2 (default) or L1
      -o: output timeseries file name
  
  Example:
      igram_inversion.py Seeded_unwrapIfgram.h5
      igram_inversion.py -f Seeded_unwrapIfgram.h5 -l L1

********************************************************************************
    '''

######################################
def main(argv):

    inversion_method = 'l2'
    #maskFile = 'Mask.h5'
  
    if len(sys.argv)>2:
        try:   opts, args = getopt.getopt(argv,"h:f:l:o:")
        except getopt.GetoptError:
            usage() ; sys.exit(1)
  
        for opt,arg in opts:
            if opt in ("-h","--help"):  usage();   sys.exit()
            elif opt == '-f':           igramsFile        = arg
            elif opt == '-l':           inversion_method  = arg.lower()
            elif opt == '-o':           timeseriesFile    = arg
  
    elif len(sys.argv)==2:
        if os.path.isfile(argv[0]):     igramsFile = argv[0]
        else:  usage(); sys.exit(1)
    else:  usage(); sys.exit(1)
    
    #try:     igramsFile = argv[0]
    #except:  usage() ; sys.exit(1)
    #try:     inversion_method = argv[1]
    #except:  inversion_method = 'L2'
  
    try:    timeseriesFile
    except: timeseriesFile = 'timeseries.h5'
  
    #h5file = h5py.File(igramsFile,'r')
    atr = readfile.read_attribute(igramsFile)
    if not atr['FILE_TYPE'] == 'interferograms':
        print '**********************************************************************'
        print 'ERROR:'
        print '       '+igramsFile+ '  was not found or the file is not readable!'
        print '**********************************************************************'
        usage();sys.exit(1)

    #numIfgrams = len(h5file['interferograms'].keys())
    #if numIfgrams == 0.:
    #    print "load interferograms first by running: load_data.py TEMPLATEFILE"
    #    sys.exit(1)
    #h5timeseries = h5py.File(timeseriesFile,'w')
  
    #print '\n************** Inverse Time Series ****************'
    if not inversion_method == 'l1':
        print 'Inverse time series using L2 norm minimization'
        ut.timeseries_inversion(igramsFile,timeseriesFile)
    else:
        print 'Inverse time series using L1 norm minimization'
        ut.timeseries_inversion_L1(igramsFile,timeseriesFile)
  
    ## generate 'mask' for timeseries.h5
    #print 'Generate mask group in timeseries file.'
    #try:
    #    maskFile
    #    h5filem = h5py.File(maskFile,'r'); kM = h5filem.keys()
    #    dset1 = h5filem[kM[0]].get(kM[0])
    #    Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
    #    group=h5timeseries.create_group('mask')
    #    dset = group.create_dataset('mask', data=Mask, compression='gzip')
    #    h5filem.close()
    #    print 'mask: '+maskFile
    #except:
    #    try:
    #        dset1 = h5file['mask'].get('mask')
    #        Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]] 
    #        group=h5timeseries.create_group('mask')
    #        dset = group.create_dataset('mask', data=Mask, compression='gzip')
    #        print 'mask: mask group in '+igramsFile
    #    except:
    #        print 'No mask found, cannot inverse timeseries without mask!'
    #        sys.exit(1)
  
    #h5file.close()
    #h5timeseries.close()
    #h5flat.close()

######################################
if __name__ == '__main__':
    main(sys.argv[1:])

