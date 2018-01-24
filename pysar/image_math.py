#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2015, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import argparse

import h5py

import _datetime as ptime
import _readfile as readfile
import _writefile as writefile
from _readfile import multi_group_hdf5_file, multi_dataset_hdf5_file


#######################################################################################
def data_operation(data, operator, operand):
    '''Mathmatic operation of 2D matrix'''
    if   operator == '+':  data += operand
    elif operator == '-':  data -= operand
    elif operator == '*':  data *= operand
    elif operator == '/':  data *= 1.0/operand
    elif operator == '^':  data = data**operand
    return data


def file_operation(fname, operator, operand, fname_out=None):
    '''Mathmathic operation of file'''

    # Basic Info
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    print('input is '+k+' file: '+fname)
    print('operation: file %s %f' % (operator, operand))

    # default output filename
    if not fname_out:
        if   operator in ['+','plus',  'add',      'addition']:        suffix = 'plus'
        elif operator in ['-','minus', 'substract','substraction']:    suffix = 'minus'
        elif operator in ['*','times', 'multiply', 'multiplication']:  suffix = 'multiply'
        elif operator in ['/','obelus','divide',   'division']:        suffix = 'divide'
        elif operator in ['^','pow','power']:                          suffix = 'pow'
        ext = os.path.splitext(fname)[1]
        fname_out = os.path.splitext(fname)[0]+'_'+suffix+str(operand)+ext

    ##### Multiple Dataset HDF5 File
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5 = h5py.File(fname,'r')
        epoch_list = sorted(h5[k].keys())
        epoch_num = len(epoch_list)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

        h5out = h5py.File(fname_out,'w')
        group = h5out.create_group(k)
        print('writing >>> '+fname_out)

        if k == 'timeseries':
            print('number of acquisitions: '+str(epoch_num))
            for i in range(epoch_num):
                date = epoch_list[i]
                data = h5[k].get(date)[:]

                data_out = data_operation(data, operator, operand)

                dset = group.create_dataset(date, data=data_out, compression='gzip')
                prog_bar.update(i+1, suffix=date)
            for key,value in atr.items():
                group.attrs[key] = value

        elif k in ['interferograms','wrapped','coherence']:
            print('number of interferograms: '+str(epoch_num))
            date12_list = ptime.list_ifgram2date12(epoch_list)
            for i in range(epoch_num):
                ifgram = epoch_list[i]
                data = h5[k][ifgram].get(ifgram)[:]

                data_out = data_operation(data, operator, operand)

                gg = group.create_group(ifgram)
                dset = gg.create_dataset(ifgram, data=data_out, compression='gzip')
                for key, value in h5[k][ifgram].attrs.items():
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])

        h5.close()
        h5out.close()
        prog_bar.close()

    ##### Duo datasets non-HDF5 File
    elif k in ['.trans']:
        rg, az, atr = readfile.read(fname)
        rg_out = data_operation(rg, operator, operand)
        az_out = data_operation(az, operator, operand)
        print('writing >>> '+fname_out)
        writefile.write(rg_out, az_out, atr, fname_out)

    ##### Single Dataset File
    else:
        data, atr = readfile.read(fname)
        data_out = data_operation(data, operator, operand)
        print('writing >>> '+fname_out)
        writefile.write(data_out, atr, fname_out)

    return fname_out


#######################################################################################
EXAMPLE='''example:
  image_math.py  velocity.h5            '+'  0.5
  image_math.py  geo_080212_101120.cor  '-'  0.2
  image_math.py  timeseries.h5          '*'  1.5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Basic Mathmatic Operation of file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='input file')
    parser.add_argument('-o','--output', dest='outfile', help='output file name.')
    parser.add_argument('operator', choices=['+','-','*','/','^'], help='mathmatical operator')
    parser.add_argument('operand', metavar='VALUE', type=float, help='value to be operated with input file')

    inps = parser.parse_args()
    return inps


#######################################################################################
def main(argv):
    inps = cmdLineParse()

    inps.outfile = file_operation(inps.file, inps.operator, inps.operand, inps.outfile)  
    print('Done.')
    return inps.outfile


#######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

