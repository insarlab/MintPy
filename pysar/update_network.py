#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Yunjun, Jun 2016: Support multiple Pairs list file
#


import os
import sys
import getopt

import h5py

import pysar._network as pnet


######################################
def Usage():
    print '''
  *********************************************************

  Update the network of coherence/interferograms based on 
    the reference network of interferograms/coherence.

  usage:
       updateNetwork.py -f inputFile -r referenceFile

    -f : input file stored in hdf5 file format
    -r : reference file stored in hdf5 file format
         or pairs list file (list)
    -o : output file name

  Example:
       update_network.py -f Coherence.h5 -r Modified_LoadedData.h5
       update_network.py -f LoadedData.h5 -r Pairs.list
       update_network.py -f LoadedData.h5 -r Pairs1.list,Pairs2.list
       update_network.py -f LoadedData.h5 -r Paris_pre1995.list -o LoadedData_pre1995.h5

  *********************************************************
    '''

######################################
def main(argv):

    if len(sys.argv)>4:
        try:        opts, args = getopt.getopt(argv,"h:f:o:r:")
        except getopt.GetoptError:       Usage();  sys.exit(1)
   
        for opt,arg in opts:
            if opt in ("-h","--help"):  Usage(); sys.exit()
            elif opt == '-f':           inFile  = arg
            elif opt == '-r':           refFile = arg
            elif opt == '-o':           outFile = arg
           
        try:
            inFile
            refFile
        except:  Usage(); sys.exit(1)
    else:       Usage(); sys.exit(1)
  
    try:     outFile
    except:  outFile = 'Modified_'+inFile
  
    ######################################
    h5fileIn = h5py.File(inFile)
    inK0=h5fileIn.keys()[0]
    if not inK0  in ('interferograms','coherence'):
        print 'Input file           should be interferograms/coherence';
        Usage(); sys.exit(1)
    print '\n************* Update Network ***************************'
    print 'Match/Update epoch info in '+inFile+' to '+refFile
  
    ######## Date List Info ##############
    ##### Input File
    inList  = h5fileIn[inK0].keys();     inList=sorted(inList);
    inDate12List  = [i.split('-sim')[0].split('filt_')[-1] for i in inList]
    print 'number of pairs in input     file: '+str(len(inList))
  
    ##### Reference File
    refDate12List=[]
    refExt = os.path.splitext(refFile)[1].lower()

    if refExt == '.h5':
        h5fileRef = h5py.File(refFile)
        refK0=h5fileRef.keys()[0]
        refList = h5fileRef[refK0].keys();  refList=sorted(refList);
        refDate12List = [i.split('-sim')[0].split('filt_')[-1] for i in refList]
        h5fileRef.close()
    else:
        for ref_file in refFile.split(','):
            print 'reading paris info in '+ref_file
            fl = open(ref_file)
            ref_list = fl.read().splitlines()
            new_list = list(set(ref_list) - set(refDate12List))
            print 'number of pairs: '+str(len(ref_list))
            print 'number of unique pairs to add: '+str(len(new_list))
            refDate12List += new_list
            fl.close()
        refDate12List = sorted(refDate12List)
  
    print 'number of pairs in reference file: '+str(len(refDate12List))
  
    rmDate12List   = list( set(inDate12List)  - set(refDate12List) )
    lackDate12List = list( set(refDate12List) - set(inDate12List) )
    if not lackDate12List==[]:
        print 'Warning: the following DATE12 are missing in input file, may be a problem for later processing:'
        print lackDate12List
    else: print 'Input file has all dates of reference file.'
    print '-------------------------------------------------'
    print 'DATE12 to be removed: '+str(len(rmDate12List))
    print rmDate12List
    print '-------------------------------------------------\n'
  
    ######### Write Updated File #########
    h5fileOut = h5py.File(outFile,'w')
    gg = h5fileOut.create_group(inK0)
  
    print 'writing the updated '+inK0+' file into '+outFile+' ...'
    print 'number of epochs: '+str(len(refDate12List))
    i = 1
    for inEpoch in inList:
        inDate12 = inEpoch.split('-sim')[0].split('filt_')[-1]
        if inDate12 in refDate12List:
            print inEpoch+'   '+str(i)
            unw = h5fileIn[inK0][inEpoch].get(inEpoch)[:]
    
            group = gg.create_group(inEpoch)
            dset = group.create_dataset(inEpoch, data=unw, compression='gzip')
            for key, value in h5fileIn[inK0][inEpoch].attrs.iteritems():
                group.attrs[key] = value
            i += 1
   
    h5fileIn.close()
    h5fileOut.close()


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



