#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################


import argparse
import os
import sys

import numpy as np

from mintpy.utils import arg_utils, readfile, utils as ut, writefile

DATA_TYPE_STR2OBJ = {
    "bool"       : np.bool_,
    "byte"       : np.byte,
    "int16"      : np.int16,
    "float32"    : np.float32,
    "float64"    : np.float64,
    "complex64"  : np.complex64,
    "complex128" : np.complex128,
}


############################################################
EXAMPLE = """example:
  load2hdf5.py waterMask.rdr --dtype bool --dname mask -o waterMask.h5
  load2hdf5.py SanAndreas.dem --dtype float32 --dname height -o dem.h5
  load2hdf5.py './20*/*.slc.full' --dtype complex64 --dname timeseries -o slcStack.h5
  load2hdf5.py './20*/*.slc.full' --dtype complex64 --dname timeseries -o slcStack.h5 --sub-x 4400 6000 --sub-y 1000 1600
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Load binary data file(s) into an HDF5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    # input
    parser.add_argument('file', help='file to be loaded.')
    parser.add_argument('--dtype','--data-type', dest='data_type', choices=DATA_TYPE_STR2OBJ.keys(),
                        help='output data type')
    parser = arg_utils.add_subset_argument(parser, geo=False)

    # output
    parser.add_argument('--dname','--dset-name', dest='dset_name', help='output dataset name(s)')
    parser.add_argument('--meta', dest='metadata', nargs='*', help='add custom metadata')
    parser.add_argument('--comp','--compression', dest='compression', choices={None, 'lzf', 'gzip'},
                        default=None, help='compression while writing to HDF5 file (default: %(default)s).')

    parser.add_argument('-o', '--output', dest='outfile', required=True, help='output HDF5 file name')
    parser.add_argument('--force', dest='force', action='store_true', help='enforce output data overwrite.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # input file list
    inps.file = ut.get_file_list(inps.file)
    if len(inps.file) < 1:
        raise FileNotFoundError('Input file(s) NOT found!')

    # --dtype
    inps.data_type = DATA_TYPE_STR2OBJ[inps.data_type]

    # --output
    if os.path.isfile(inps.outfile) and not inps.force:
        raise FileExistsError('Output file already exists! To overwrite, re-run with "--force" option.')

    return inps


def read_inps2box(inps, atr):
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    box = [0, 0, width, length]

    if inps.subset_x:
        box[0] = inps.subset_x[0]
        box[2] = inps.subset_x[1]

    if inps.subset_y:
        box[1] = inps.subset_y[0]
        box[3] = inps.subset_y[1]

    print(f'input file length / width: {length} / {width}')
    print(f'read bounding box: {box}')

    return box


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    atr = readfile.read_attribute(inps.file[0])
    box = read_inps2box(inps, atr)
    num_file = len(inps.file)

    # initiate output
    dsDict = {}
    if num_file == 1:
        dsDict[inps.dset_name] = np.zeros((box[3]-box[1], box[2]-box[0]), dtype=inps.data_type)
    else:
        dsDict[inps.dset_name] = np.zeros((num_file, box[3]-box[1], box[2]-box[0]), dtype=inps.data_type)

        # add "date" dataset for timeseries
        if inps.dset_name and inps.dset_name == 'timeseries':
            date_list = [os.path.basename(os.path.dirname(x)) for x in inps.file]
            dsDict['date'] = np.array(date_list, dtype=np.bytes_)

    # metadata
    if inps.metadata:
        print('add/update the following metadata:')
        for meta_str in inps.metadata:
            key, value = meta_str.split('=')
            atr[key] = value
            print(f'{key} : {value}')

    # read
    for i, fname in enumerate(inps.file):
        print(f'reading file {i+1} / {num_file}: {fname}')
        ds_names = readfile.get_dataset_list(fname)
        ds_name = inps.dset_name if inps.dset_name in ds_names else None
        data = readfile.read(fname, datasetName=ds_name, box=box)[0]
        if num_file == 1:
            dsDict[inps.dset_name][:] = data
        else:
            dsDict[inps.dset_name][i, :, :] = data

    # write
    atr['LENGTH'] = box[3] - box[1]
    atr['WIDTH'] = box[2] - box[0]
    writefile.write(dsDict, out_file=inps.outfile, metadata=atr, compression=inps.compression)

    return inps.outfile


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
