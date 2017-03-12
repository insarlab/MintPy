#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jan 2016: add out_name option, support ROI_PAC product
#                   support coherence/wrapped
#                   nan + value = value for ROI_PAC product
# Yunjun, Jun 2016: support multiple input files


import sys
import os
import getopt

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile


def add(data1,data2):
    data = data1 + data2;
    data[np.isnan(data1)] = data2[np.isnan(data1)];
    data[np.isnan(data2)] = data1[np.isnan(data2)];
  
    return data


def usage():
    print '''
***************************************************************
  Generates the sum of two input files.

  Usage:
      add.py file1 file2 [out_name]
      add.py -f file1,file2,... -o out_name

      -f : files to be added, need to be the same format, supported file format:
      -o : output filename, optional [default is file1_plus_file2.h5]

  example:
           
      add.py velocity_masked.h5    velocity_demCor_masked.h5
      add.py timeseries.h5         timeseries_demCor.h5
      add.py unwrapIfgram.h5       unwrapIfgram2.h5
      add.py 081008_100220.unw     100220_110417.unw      081008_110417.unw

      add.py -f mask_1.h5,mask_2.h5,mask_3.h5       -o mask_all.h5

***************************************************************
    '''

################################################################################

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
  
    print '\n****************** Add **********************'
    print 'Input files: '
    print fileList
  
    ext = os.path.splitext(fileList[0])[1].lower()
    try:     outName
    except:  outName = File1.split('.')[0]+'_plus_'+File2.split('.')[0]+ext
  
  
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
  
        h5in  = h5py.File(fileList[0])
        epochList = sorted(h5in[k].keys())

    ########################### Add file by file ########################
    if k in ['timeseries']:
        for epoch in epochList:
            print epoch
            data = np.zeros((int(atr['FILE_LENGTH']),int(atr['WIDTH'])))
            for File in fileList:
                print File
                h5file = h5py.File(File,'r')
                d = h5file[k].get(epoch)[:]
  
                data = add(data,d)
  
            dset = group.create_dataset(epoch, data=data, compression='gzip')
        for key,value in atr.iteritems():   group.attrs[key] = value
  
        h5out.close()
        h5in.close()
  
    elif k in ['timeseries','interferograms','coherence','wrapped']:
        for epoch in epochList:
            print epoch
            data = np.zeros((int(atr['FILE_LENGTH']),int(atr['WIDTH'])))
            for File in fileList:
                print File
                h5file = h5py.File(File,'r')
                d = h5file[k][epoch].get(epoch)[:]
  
                data = add(data,d)
  
            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in h5in[k][epoch].attrs.iteritems():
                gg.attrs[key] = value
  
        h5out.close()
        h5in.close()
  
    ## All the other file types
    else:
        data = np.zeros((int(atr['FILE_LENGTH']),int(atr['WIDTH'])))
        for File in fileList:
            print 'loading '+File
            d,r = readfile.read(File)
            data = add(data,d)
        writefile.write(data,atr,outName)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

