#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jan 2016: add out_name option, support ROI_PAC product
#                   support coherence/wrapped
#                   nan + value = value for ROI_PAC product
# Yunjun, Jun 2016: support multiple input files


import os
import sys
import argparse

import h5py
import numpy as np

import _datetime as ptime
import _readfile as readfile
import _writefile as writefile
from _readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


################################################################################
def add_matrix(data1,data2):
    '''Sum of 2 input matrix'''
    data = data1 + data2
    data[np.isnan(data1)] = data2[np.isnan(data1)]
    data[np.isnan(data2)] = data1[np.isnan(data2)]
    return data


def add_files(fname_list, fname_out=None):
    '''Generate sum of all input files
    Inputs:
        fname_list - list of string, path/name of input files to be added
        fname_out  - string, optional, path/name of output file
    Output:
        fname_out  - string, path/name of output file
    Example:
        'mask_all.h5' = add_file(['mask_1.h5','mask_2.h5','mask_3.h5'], 'mask_all.h5')
    '''
    # Default output file name
    ext = os.path.splitext(fname_list[0])[1]
    if not fname_out:
        fname_out = os.path.splitext(fname_list[0])[0]
        for i in range(1,len(fname_list)):
            fname_out += '_plus_'+os.path.splitext(os.path.basename(fname_list[i]))[0]
        fname_out += ext

    # Basic Info
    atr  = readfile.read_attribute(fname_list[0])
    k = atr['FILE_TYPE']
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])
    print(('First input file is '+atr['PROCESSOR']+' '+k))

    ## Multi-dataset/group file
    if k in multi_group_hdf5_file + multi_dataset_hdf5_file:
        # File Type Check
        for i in range(1,len(fname_list)):
            ki = readfile.read_attribute(fname_list[i])['FILE_TYPE']
            if (k in multi_dataset_hdf5_file and ki in multi_dataset_hdf5_file or
                k in multi_group_hdf5_file   and ki in multi_group_hdf5_file):
                pass
            else:
                print(('Input files structure are not the same: '+k+' v.s. '+ki))
                sys.exit(1)

        print(('writing >>> '+fname_out))
        h5out = h5py.File(fname_out, 'w')
        group = h5out.create_group(k)

        h5  = h5py.File(fname_list[0], 'r')
        epoch_list = sorted(h5[k].keys())
        epoch_num = len(epoch_list)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

    if k in multi_dataset_hdf5_file:
        print(('number of acquisitions: %d' % epoch_num))
        for i in range(epoch_num):
            epoch = epoch_list[i]
            data = np.zeros((length,width))
            for fname in fname_list:
                h5file = h5py.File(fname, 'r')
                d = h5file[k].get(epoch)[:]
                data = add_matrix(data, d)

            dset = group.create_dataset(epoch, data=data, compression='gzip')
            prog_bar.update(i+1, suffix=epoch)

        for key,value in list(atr.items()):
            group.attrs[key] = value
        h5out.close()
        h5.close()
        prog_bar.close()
  
    elif k in multi_group_hdf5_file:
        print(('number of interferograms: %d' % epoch_num))
        date12_list = ptime.list_ifgram2date12(epoch_list)
        for i in range(epoch_num):
            epoch = epoch_list[i]
            data = np.zeros((length,width))
            for fname in fname_list:
                h5file = h5py.File(fname,'r')
                temp_k = list(h5file.keys())[0]
                temp_epoch_list = sorted(h5file[temp_k].keys())
                d = h5file[temp_k][temp_epoch_list[i]].get(temp_epoch_list[i])[:]
                data = add_matrix(data,d)

            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in list(h5[k][epoch].attrs.items()):
                gg.attrs[key] = value
            prog_bar.update(i+1, suffix=date12_list[i])
        h5out.close()
        h5.close()
        prog_bar.close()

    ## Single dataset files
    else:
        data = np.zeros((length,width))
        for fname in fname_list:
            print(('loading '+fname))
            d,r = readfile.read(fname)
            data = add_matrix(data,d)

        print(('writing >>> '+fname_out))
        writefile.write(data,atr,fname_out)

    return fname_out


################################################################################
EXAMPLE='''example:
  add.py  mask_1.h5 mask_2.h5 mask_3.h5           -o mask_all.h5
  add.py  081008_100220.unw    100220_110417.unw  -o 081008_110417.unw
  add.py  timeseries_ECMWF.h5  ECMWF.h5           -o  timeseries.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Generate sum of multiple input files.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='files (2 or more) to be added')
    parser.add_argument('-o','--output', dest='outfile', help='output file name')

    inps = parser.parse_args()
    if len(inps.file) < 2:
        parser.print_usage()
        sys.exit(os.path.basename(sys.argv[0])+': error: number of input files is < 2')
    return inps


################################################################################
def main(argv):
    inps = cmdLineParse()
    print('Input files to be added: ')
    print((inps.file))

    inps.outfile = add_files(inps.file, inps.outfile)
    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

