#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Mar 2016: add diff()
#
#

import sys
import os
import getopt

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile


#####################################################################################

def diff(data1,data2):
    data = data1 - data2;
    #data[np.isnan(data2)] = data1[np.isnan(data2)];
  
    return data

def usage():
    print '''
***************************************************************
  Generates the difference of two input files.

  Usage:
      diff.py file1 file2 [ output_file ]

      output is file1_diff_file2.h5 by default

  Example:
      diff.py velocity_masked.h5 velocity_demCor_masked.h5    demCor.h5
      diff.py timeseries.h5      timeseries_demCor.h5         demCor.h5
      diff.py LoadedData.h5      reconstruct_LoadedData.h5
      diff.py -f velocity.h5,velocity_2.h5  -o velocity_diff.h5

***************************************************************
    '''


#####################################################################################

def main(argv):

    ####################### Inputs Check ########################
    try:    opts, args = getopt.getopt(argv,"h:f:o:",['help'])
    except getopt.GetoptError:    usage() ; sys.exit(1)
  
    if len(sys.argv) > 4:
        for opt,arg in opts:
            if opt in ("-h","--help"):  usage();  sys.exit()
            elif opt == '-f':   fileList = arg.split(',')
            elif opt == '-o':   outName  = arg
  
    elif len(sys.argv) <= 4 and len(sys.argv) >= 3:
        fileList = [sys.argv[1],sys.argv[2]]
        try: outName = sys.argv[3]
        except: pass
    else: usage();  sys.exit(1)
  
    print '\n**************** Diff *******************'
    print 'Input files: '
    print fileList
  
    ext = os.path.splitext(fileList[0])[1].lower()
    try:     outName
    except:  outName = fileList[0].split('.')[0]+'_diff_'+fileList[1].split('.')[0]+ext
  
    ##### Read File Info / Attributes
    atr  = readfile.read_attribute(fileList[0])
    print 'Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']
    k = atr['FILE_TYPE']
  
    ##### File Type Check
    if k in ['timeseries','interferograms','coherence','wrapped']:
        for i in range(1,len(fileList)):
            File = fileList[i]
            r = readfile.read_attribute(File)
            if not r['FILE_TYPE'] == k:
                print 'Input file type is not the same: '+r['FILE_TYPE']
                sys.exit(1)
  
        h5out = h5py.File(outName,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+outName
  
        h5in  = h5py.File(fileList[0])
        epochList = h5in[k].keys()
        print 'number of epochs: '+str(len(epochList))

    ########################### Diff file by file ########################
    if k in ['timeseries']:
        for epoch in epochList:
            print epoch
            data = h5in[k].get(epoch)[:]
            for i in range(1,len(fileList)):
                #print fileList[i]
                h5file = h5py.File(fileList[i],'r')
                d = h5file[k].get(epoch)[:]
  
                data = diff(data,d)
  
            dset = group.create_dataset(epoch, data=data, compression='gzip')
        for key,value in atr.iteritems():   group.attrs[key] = value
  
        h5out.close()
        h5in.close()
  
    elif k in ['interferograms','coherence','wrapped']:
        for epoch in epochList:
            print epoch
            data = h5in[k][epoch].get(epoch)[:]
            for i in range(1,len(fileList)):
                #print fileList[i]
                h5file = h5py.File(fileList[i],'r')
                d = h5file[k][epoch].get(epoch)[:]
  
                data = diff(data,d)
  
            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in h5in[k][epoch].attrs.iteritems():
                gg.attrs[key] = value
  
        h5out.close()
        h5in.close()
  
    ## All the other file types
    else:
        data,atr = readfile.read(fileList[0])
        for i in range(1,len(fileList)):
            #print fileList[i]
            d,r = readfile.read(fileList[i])
            data = diff(data,d)
        writefile.write(data,atr,outName)


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

