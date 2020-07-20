#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import sys
import argparse
import warnings
import numpy as np

# suppress numpy.RuntimeWarning message
import logging
np_logger = logging.getLogger('numpy')
np_logger.setLevel(logging.WARNING)

from mintpy.utils import readfile, writefile, utils as ut


##################################################################################################
EXAMPLE = """example:
  multilook.py  velocity.h5  -r 15 -a 15
  multilook.py  srtm30m.dem  -r 10 -a 10  -o srtm30m_300m.dem

  # Ignore / skip marginal pixels
  multilook.py ../../geom_reference/hgt.rdr.full -r 300 -a 100 --margin 58 58 58 58 -o hgt.rdr
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Multilook.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to multilook')
    parser.add_argument('-r','--range','-x', dest='lks_x', type=int,
                        help='number of multilooking in range  /x direction')
    parser.add_argument('-a','--azimuth','-y', dest='lks_y', type=int,
                        help='number of multilooking in azimuth/y direction')
    parser.add_argument('-o', '--outfile',
                        help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--margin', dest='margin', type=int, nargs=4, metavar=('TOP','BOTTOM','LEFT','RIGHT'),
                        default=[0,0,0,0], help='number of pixels on the margin to skip, default: 0 0 0 0.')
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


def multilook_data(data, lks_y, lks_x):
    """Modified from Praveen on StackOverflow:

    link: https://stackoverflow.com/questions/34689519/how-to-coarser-the-2-d-array-data-resolution

    Parameters: data        : 2D / 3D np.array in real or complex
                lks_y       : int, number of multilook in y/azimuth direction
                lks_x       : int, number of multilook in x/range direction
    Returns:    coarse_data : 2D / 3D np.array after multilooking in last two dimension
    """

    shape = np.array(data.shape, dtype=float)
    if len(shape) == 2:
        # crop data to the exact multiple of the multilook number
        new_shape = np.floor(shape / (lks_y, lks_x)).astype(int) * (lks_y, lks_x)
        crop_data = data[:new_shape[0],
                         :new_shape[1]]

        # reshape to more dimensions and collapse the extra dimensions with mean
        temp = crop_data.reshape((new_shape[0] // lks_y, lks_y,
                                  new_shape[1] // lks_x, lks_x))
        coarse_data = np.nanmean(temp, axis=(1, 3))

    elif len(shape) == 3:
        # crop data to the exact multiple of the multilook number
        new_shape = np.floor(shape / (1, lks_y, lks_x)).astype(int) * (1, lks_y, lks_x)
        crop_data = data[:new_shape[0],
                         :new_shape[1],
                         :new_shape[2]]

        # reshape to more dimensions and collapse the extra dimensions with mean
        temp = crop_data.reshape((new_shape[0],
                                  new_shape[1] // lks_y, lks_y,
                                  new_shape[2] // lks_x, lks_x))
        coarse_data = np.nanmean(temp, axis=(2, 4))

    return coarse_data


def multilook_attribute(atr_dict, lks_y, lks_x, box=None, print_msg=True):
    # make a copy of original meta dict
    atr = dict()
    for key, value in iter(atr_dict.items()):
        atr[key] = str(value)

    if box is None:
        box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    length, width = box[3] - box[1], box[2] - box[0]

    length_mli = length // lks_y
    width_mli = width // lks_x
    print('output data in size: {}, {}'.format(length_mli, width_mli))

    # Update attributes
    atr['LENGTH'] = str(length_mli)
    atr['WIDTH'] = str(width_mli)
    atr['XMIN'] = str(box[0])
    atr['YMIN'] = str(box[1])
    atr['XMAX'] = str(width_mli - 1 + box[0])
    atr['YMAX'] = str(length_mli - 1 + box[1])
    atr['RLOOKS'] = str(int(atr.get('RLOOKS', '1')) * lks_x)
    atr['ALOOKS'] = str(int(atr.get('ALOOKS', '1')) * lks_y)
    if print_msg:
        print('update LENGTH, WIDTH, Y/XMIN/MAX, A/RLOOKS')

    if 'Y_STEP' in atr.keys():
        atr['Y_STEP'] = str(lks_y * float(atr['Y_STEP']))
        atr['X_STEP'] = str(lks_x * float(atr['X_STEP']))
        if print_msg:
            print('update Y/X_STEP')

    if 'AZIMUTH_PIXEL_SIZE' in atr.keys():
        atr['AZIMUTH_PIXEL_SIZE'] = str(lks_y * float(atr['AZIMUTH_PIXEL_SIZE']))
        if print_msg:
            print('update AZIMUTH_PIXEL_SIZE')

    if 'RANGE_PIXEL_SIZE' in atr.keys():
        atr['RANGE_PIXEL_SIZE'] = str(lks_x * float(atr['RANGE_PIXEL_SIZE']))
        if print_msg:
            print('update RANGE_PIXEL_SIZE')

    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str( (int(atr['REF_Y']) - box[1]) // lks_y )
        atr['REF_X'] = str( (int(atr['REF_X']) - box[0]) // lks_x )
        if print_msg:
            print('update REF_Y/X')

    if 'SUBSET_XMIN' in atr.keys():
        atr['SUBSET_YMIN'] = str( (int(atr['SUBSET_YMIN']) - box[1]) // lks_y )
        atr['SUBSET_YMAX'] = str( (int(atr['SUBSET_YMAX']) - box[1]) // lks_y )
        atr['SUBSET_XMIN'] = str( (int(atr['SUBSET_XMIN']) - box[0]) // lks_x )
        atr['SUBSET_XMAX'] = str( (int(atr['SUBSET_XMAX']) - box[0]) // lks_x )
        if print_msg:
            print('update SUBSET_XMIN/XMAX/YMIN/YMAX')
    return atr


def multilook_file(infile, lks_y, lks_x, outfile=None, margin=[0,0,0,0]):
    """ Multilook input file
    Parameters: infile - str, path of input file to be multilooked.
                lks_y  - int, number of looks in y / row direction.
                lks_x  - int, number of looks in x / column direction.
                margin - list of 4 int, number of pixels to be skipped during multilooking.
                         useful for offset product, where the marginal pixels are ignored during
                         cross correlation matching.
                outfile - str, path of output file
    Returns:    outfile - str, path of output file
    """
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    # input file info
    atr = readfile.read_attribute(infile)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    k = atr['FILE_TYPE']
    print('multilooking {} {} file: {}'.format(atr['PROCESSOR'], k, infile))
    print('number of looks in y / azimuth direction: %d' % lks_y)
    print('number of looks in x / range   direction: %d' % lks_x)

    # margin --> box
    if margin is not [0,0,0,0]:    # top, bottom, left, right
        box = (margin[2], margin[0], width - margin[3], length - margin[1])
        print('number of pixels to skip in top/bottom/left/right boundaries: {}'.format(margin))
    else:
        box = (0, 0, width, length)

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
        data = readfile.read(infile, datasetName=dsName, box=box, print_msg=False)[0]

        # keep timeseries data as 3D matrix when there is only one acquisition
        # because readfile.read() will squeeze it to 2D
        if atr['FILE_TYPE'] == 'timeseries' and len(data.shape) == 2:
            data = np.reshape(data, (1, data.shape[0], data.shape[1]))

        data = multilook_data(data, lks_y, lks_x)
        dsDict[dsName] = data

    # update metadata
    atr = multilook_attribute(atr, lks_y, lks_x, box=box)

    # for binary file with 2 bands, always use BIL scheme
    if (len(dsDict.keys()) == 2
            and os.path.splitext(infile)[1] not in ['.h5','.he5']
            and atr.get('scheme', 'BIL').upper() != 'BIL'):
        print('the input binary file has 2 bands with band interleave as: {}'.format(atr['scheme']))
        print('for the output binary file, change the band interleave to BIL as default.')
        atr['scheme'] = 'BIL'

    writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)
    return outfile


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    for infile in inps.file:
        multilook_file(infile,
                       lks_y=inps.lks_y,
                       lks_x=inps.lks_x,
                       outfile=inps.outfile,
                       margin=inps.margin)

    print('Done.')
    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
