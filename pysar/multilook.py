#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse
import warnings
import numpy as np
from pysar.utils import readfile, writefile, utils as ut


##################################################################################################
EXAMPLE = """example:
  multilook.py  velocity.h5  -r 15 -a 15
  multilook.py  srtm30m.dem  -r 10 -a 10  -o srtm30m_300m.dem

  To interpolate input file into larger size file:
  multilook.py  bperp.rdr  -10 -2 -o bperp_full.rdr
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Multilook.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to multilook')
    parser.add_argument('-r','--range','-x', dest='lks_x', type=int,
                        help='number of multilooking in azimuth/y direction')
    parser.add_argument('-a','--azimuth','-y', dest='lks_y', type=int,
                        help='number of multilooking in range  /x direction')
    parser.add_argument('-o', '--outfile',
                        help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false', default=True,
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file)

    if len(inps.file) > 1 and inps.outfile:
        inps.outfile = None
        print('more than one file is input, disable custom output filename.')
    return inps


######################################## Sub Functions ############################################
def multilook_matrix(matrix, lks_y, lks_x):
    """Obsolete multilooking functions"""
    rows, cols = matrix.shape
    lks_x = int(lks_x)
    lks_y = int(lks_y)
    if lks_x == 1 and lks_y == 1:
        return matrix

    rows_mli = int(np.floor(rows/lks_y))
    cols_mli = int(np.floor(cols/lks_x))
    #thr = np.floor(lks_x*lks_y/2)
    matrix_Cmli = np.zeros((rows,    cols_mli))
    matrix_mli = np.zeros((rows_mli, cols_mli))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for c in range(cols_mli):
            matrix_Cmli[:, c] = np.nanmean(matrix[:, (c)*lks_x:(c+1)*lks_x], 1)
        for r in range(rows_mli):
            matrix_mli[r, :] = np.nanmean(matrix_Cmli[(r)*lks_y:(r+1)*lks_y, :], 0)
    del matrix, matrix_Cmli
    return matrix_mli


def multilook_data(data, lksY, lksX):
    """Modified from Praveen on StackOverflow:
    https://stackoverflow.com/questions/34689519/how-to-coarser-the-2-d-array-data-resolution
    Parameters: data : 2D / 3D np.array
                lksY : int, number of multilook in y/azimuth direction
                lksX : int, number of multilook in x/range direction
    Returns:    coarseData : 2D / 3D np.array after multilooking in last two dimension
    """
    shape = np.array(data.shape, dtype=float)
    if len(shape) == 2:
        newShape = np.floor(shape / (lksY, lksX)).astype(int) * (lksY, lksX)
        cropData = data[:newShape[0], :newShape[1]]
        temp = cropData.reshape((newShape[0] // lksY, lksY,
                                 newShape[1] // lksX, lksX))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            coarseData = np.nanmean(temp, axis=(1, 3))
    elif len(shape) == 3:
        newShape = np.floor(shape / (1, lksY, lksX)).astype(int) * (1, lksY, lksX)
        cropData = data[:newShape[0], :newShape[1], :newShape[2]]
        temp = cropData.reshape((newShape[0],
                                 newShape[1] // lksY, lksY,
                                 newShape[2] // lksX, lksX))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            coarseData = np.nanmean(temp, axis=(2, 4))
    return coarseData


def multilook_attribute(atr_dict, lks_y, lks_x, print_msg=True):
    atr = dict()
    for key, value in iter(atr_dict.items()):
        atr[key] = str(value)

    length_mli = int(np.floor(int(atr['LENGTH']) / lks_y))
    width_mli = int(np.floor(int(atr['WIDTH']) / lks_x))

    # Update attributes
    atr['LENGTH'] = str(length_mli)
    atr['WIDTH'] = str(width_mli)
    atr['XMIN'] = '0'
    atr['YMIN'] = '0'
    atr['XMAX'] = str(width_mli-1)
    atr['YMAX'] = str(length_mli-1)
    if print_msg:
        print('update LENGTH, WIDTH, YMIN, YMAX, XMIN, XMAX')

    if 'Y_STEP' in atr.keys():
        atr['Y_STEP'] = str(lks_y * float(atr['Y_STEP']))
        atr['X_STEP'] = str(lks_x * float(atr['X_STEP']))
        if print_msg:
            print('update Y/X_STEP')

    if 'RANGE_PIXEL_SIZE' in atr.keys():
        atr['AZIMUTH_PIXEL_SIZE'] = str(lks_y * float(atr['AZIMUTH_PIXEL_SIZE']))
        atr['RANGE_PIXEL_SIZE'] = str(lks_x * float(atr['RANGE_PIXEL_SIZE']))
        if print_msg:
            print('update AZIMUTH/RANGE_PIXEL_SIZE')

    if 'Y_FIRST' not in atr.keys() and 'ALOOKS' in atr.keys():
        atr['RLOOKS'] = str(int(atr['RLOOKS']) * lks_x)
        atr['ALOOKS'] = str(int(atr['ALOOKS']) * lks_y)
        if print_msg:
            print('update R/ALOOKS')

    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str(int(int(atr['REF_Y']) / lks_y))
        atr['REF_X'] = str(int(int(atr['REF_X']) / lks_x))
        if print_msg:
            print('update REF_Y/X')

    if 'SUBSET_XMIN' in atr.keys():
        atr['SUBSET_YMIN'] = str(int(int(atr['SUBSET_YMIN'])/lks_y))
        atr['SUBSET_YMAX'] = str(int(int(atr['SUBSET_YMAX'])/lks_y))
        atr['SUBSET_XMIN'] = str(int(int(atr['SUBSET_XMIN'])/lks_x))
        atr['SUBSET_XMAX'] = str(int(int(atr['SUBSET_XMAX'])/lks_x))
        if print_msg:
            print('update SUBSET_XMIN/XMAX/YMIN/YMAX')
    return atr


def multilook_file(infile, lks_y, lks_x, outfile=None):
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    # input file info
    atr = readfile.read_attribute(infile)
    k = atr['FILE_TYPE']
    print('multilooking {} {} file: {}'.format(atr['PROCESSOR'], k, infile))
    print('number of looks in y / azimuth direction: %d' % lks_y)
    print('number of looks in x / range   direction: %d' % lks_x)

    # output file name
    if not outfile:
        if os.getcwd() == os.path.dirname(os.path.abspath(infile)):
            ext = os.path.splitext(infile)[1]
            outfile = os.path.splitext(infile)[0]+'_'+str(lks_y)+'alks_'+str(lks_x)+'rlks'+ext
        else:
            outfile = os.path.basename(infile)
    #print('writing >>> '+outfile)

    # read source data and multilooking
    dsNames = readfile.get_dataset_list(infile)
    maxDigit = max([len(i) for i in dsNames])
    dsDict = dict()
    for dsName in dsNames:
        print('multilooking {d:<{w}} from {f} ...'.format(
            d=dsName, w=maxDigit, f=os.path.basename(infile)))
        data = readfile.read(infile, datasetName=dsName, print_msg=False)[0]
        data = multilook_data(data, lks_y, lks_x)
        dsDict[dsName] = data
    atr = multilook_attribute(atr, lks_y, lks_x)
    writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)
    return outfile


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    for infile in inps.file:
        multilook_file(infile, lks_y=inps.lks_y, lks_x=inps.lks_x, outfile=inps.outfile)

    print('Done.')
    return


###################################################################################################
if __name__ == '__main__':
    main()
