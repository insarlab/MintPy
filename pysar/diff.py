#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Mar 2016: add diff_data()
# Yunjun, Apr 2017: add diff_file()


import sys
import os

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime


#####################################################################################
def diff_data(data1,data2):
    '''data1 - data2'''
    data = data1 - data2
    #data[np.isnan(data2)] = data1[np.isnan(data2)];
    return data


def diff_file(file1, file2, outName=None):
    '''Subtraction/difference of two input files'''
    if not outName:
        outName = os.path.splitext(file1)[0]+'_diff_'+os.path.splitext(os.path.basename(file2))[0]+\
                  os.path.splitext(file1)[1]
    
    print file1+' - '+file2
    # Read basic info
    atr  = readfile.read_attribute(file1)
    print 'Input first file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']
    k = atr['FILE_TYPE']
  
    # Multi-dataset/group file
    if k in ['timeseries','interferograms','coherence','wrapped']:
        # Check input files type for multi_dataset/group files
        atr2 = readfile.read_attribute(file2)
        k2 = atr2['FILE_TYPE']
  
        h5out = h5py.File(outName,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+outName
  
        h5_1  = h5py.File(file1)
        h5_2  = h5py.File(file2)
        epochList = sorted(h5_1[k].keys())
        epochList2 = sorted(h5_2[k2].keys())
        if not all(i in epochList2 for i in epochList):
            print file2+' does not contain all group of '+file1
            sys.exit(1)
        #if not len(epochList) == len(epochList2):
        #    print 'ERROR: input files have different number of datasets: '
        #    print 'number of datasets in '+file1+' : '+str(len(epochList))
        #    print 'number of datasets in '+file2+' : '+str(len(epochList2))
        #    sys.exit(1)
        epoch_num = len(epochList)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

    if k in ['timeseries']:
        print 'number of acquisitions: '+str(len(epochList))
        for i in range(epoch_num):
            date = epochList[i]
            data1 = h5_1[k].get(date)[:]
            data2 = h5_2[k2].get(date)[:]
            if atr2['ref_date'] != atr['ref_date']:
                data2 -= h5_2[k2].get(atr['ref_date'])[:]
            data = diff_data(data1, data2)
            dset = group.create_dataset(date, data=data, compression='gzip')
            prog_bar.update(i+1, suffix=date)
        for key,value in atr.iteritems():
            group.attrs[key] = value

        prog_bar.close()
        h5out.close()
        h5_1.close()
        h5_2.close()

    elif k in ['interferograms','coherence','wrapped']:
        print 'number of interferograms: '+str(len(epochList))
        date12_list = ptime.list_ifgram2date12(epochList)
        for i in range(epoch_num):
            epoch1 = epochList[i]
            epoch2 = epochList2[i]
            data1 = h5_1[k][epoch1].get(epoch1)[:]
            data2 = h5_2[k2][epoch2].get(epoch2)[:]
            data = diff_data(data1, data2)  
            gg = group.create_group(epoch1)
            dset = gg.create_dataset(epoch1, data=data, compression='gzip')
            for key, value in h5_1[k][epoch1].attrs.iteritems():
                gg.attrs[key] = value
            prog_bar.update(i+1, suffix=date12_list[i])

        prog_bar.close()
        h5out.close()
        h5_1.close()
        h5_2.close()
  
    # Sing dataset file
    else:
        data1, atr1 = readfile.read(file1)
        data2, atr2 = readfile.read(file2)
        data = diff_data(data1, data2)
        print 'writing >>> '+outName
        writefile.write(data, atr1, outName)

    return outName


#####################################################################################
def usage():
    print '''
usage:  diff.py  file1  file2  [ outfile ]

Generates the difference of two input files.

positional arguments:
  file1/2               file 1 and 2 used for differencing

optional argument:
  outfile               output file name, default is file1_diff_file2.h5

example:
  diff.py velocity_masked.h5  velocity_demCor_masked.h5    demCor.h5
  diff.py timeseries.h5       timeseries_demCor.h5         demCor.h5
  diff.py unwrapIfgram.h5     reconstruct_unwrapIfgram.h5
    '''
    return


#####################################################################################
def main(argv):
    if 3 <= len(sys.argv) <= 4:
        file1 = argv[0]
        file2 = argv[1]
        try:    outName = argv[2]
        except: outName = None
    else:
        usage(); sys.exit(1)
  
    outName = diff_file(file1, file2, outName)
    return outName


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

