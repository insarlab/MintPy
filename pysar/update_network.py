#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################


import sys
import os
import h5py
import getopt

######################################
def Usage():
  print '''
*********************************************************
*********************************************************
  Update the network of coherence/interferograms based on 
  the reference network of interferograms/coherence.

  usage:
       updateNetwork.py -f inFile -r refFile

    -f : coherence/interferograms file stored in hdf5 file format
    -r : reference coherence/interferograms file stored in hdf5 file format


  Example:
       update_network.py -f Coherence.h5 -r Modified_LoadedData.h5
  '''

######################################

def main(argv):

  if len(sys.argv)>4:

     try:
        opts, args = getopt.getopt(argv,"h:f:r:")
     except getopt.GetoptError:
        Usage();  sys.exit(1)

     for opt,arg in opts:
        if opt in ("-h","--help"):
           Usage(); sys.exit()
        elif opt == '-f':
           inFile = arg
        elif opt == '-r':
           refFile = arg
        
     try:
        inFile
        refFile
     except:
        Usage(); sys.exit(1)

  else:
     Usage(); sys.exit(1)

######################################
#  import pdb;  pdb.set_trace()

  h5fileIn = h5py.File(inFile)
  h5fileRef = h5py.File(refFile)
  inK0=h5fileIn.keys()[0]
  refK0=h5fileRef.keys()[0]
  if not inK0 in ('interferograms','coherence'):
     print 'Input file should be interferograms/coherence'
     Usage(); sys.exit(1)
  if not refK0 in ('interferograms','coherence'):
     print 'Input reference file should be interferograms/coherence'
     Usage(); sys.exit(1)

  try:
     refDate12List
  except:
     refDate12List=[]
######################################

  print ''
  print '****************************************'
  print 'DATE12 included in reference file:'
  refList=h5fileRef[refK0].keys()
  for refData in refList:
     refDate12=refData.split('-sim')[0].split('filt_')[-1]
     print refDate12
     refDate12List.append(refDate12)
  h5fileRef.close()
  lackDate12List=refDate12List
  #import pdb;  pdb.set_trace()

  outFile = 'Modified_'+inFile
  h5fileOut = h5py.File(outFile,'w')
  gg = h5fileOut.create_group(inK0)

  inList=h5fileIn[inK0].keys()
  rmDate12List=[]
  print 'writing the updated '+inK0+' file into '+outFile+' ...'
  for inData in inList:
     inDate12=inData.split('-sim')[0].split('filt_')[-1]
     if inDate12 in refDate12List:
        print inData
        lackDate12List.remove(inDate12)
        unwSet = h5fileIn[inK0][inData].get(inData)
        unw = unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]
        group = gg.create_group(inData)
        dset = group.create_dataset(inData, data=unw, compression='gzip')
        for key, value in h5fileIn[inK0][inData].attrs.iteritems():
           group.attrs[key] = value
     else:
        rmDate12List.append(inDate12)
  
  print 'Removed DATE12:'
  print rmDate12List
  print ''
  if not lackDate12List==[]:
     print 'Warning: the following DATE12 are missing in input file, may be a problem for later processing:'
     print lackDate12List
  else:
     print 'Update successfully!'
     print ''

  h5fileIn.close()
  h5fileOut.close()



############################################################

if __name__ == '__main__':
  main(sys.argv[1:])



