#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2019, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import argparse
import numpy as np
from pysar.utils import readfile, writefile


############################################################
EXAMPLE = """example:
  mask.py  velocity.h5     -m Mask.h5
  mask.py  timeseries.h5   -m temporalCoherence.h5  -t 0.7
  mask.py  ifgramStack.h5  -m 100102_101120.cor     -t 0.9  -y  200 300  -x 300 400

  mask.py  filt_20060924_20090214.int -m waterMask.h5 -o filt_20060924_20090214_msk.int
  mask.py  filt_20060924_20090214.cor -m waterMask.h5 -o filt_20060924_20090214_msk.cor
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Mask file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='File to be masked')
    parser.add_argument('-m', '--mask', dest='mask_file', required=True,
                        help='mask for pixels used in ramp estimation')
    parser.add_argument('-o', '--outfile', help='Output file name.')

    # modify input mask
    parser.add_argument('-t', dest='threshold', type=float,
                        help='threshold value used for masking.\n' +
                        'if not specified, only pixels with mask value equal to zero is masked out.')
    parser.add_argument('--fill', dest='fill_value', type=float,
                        help="fill masked out area with input value. i.e. \n"
                             "np.nan, 0, 1000, ... \n"
                             "By default, it's np.ma.masked for int16 type and np.nan for all the others.")
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2,
                        help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2,
                        help='subset range in y/along-track/row direction')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def mask_matrix(data, mask, fill_value=None):
    """mask a 2D matrxi data with mask
    Parameters: data : 2D / 3D np.array
                mask : 2D np.array of bool
                fill_value : number
    Returns:    data : same shape of array as input
    """
    # Masked Value
    if fill_value is None:
        if data.dtype == np.float32:
            fill_value = np.nan
        else:
            raise ValueError('specify fill_value for input data type: {}'.format(data.dtype))

    if len(data.shape) == 2:
        data[mask == 0] = fill_value
    elif len(data.shape) == 3:
        data[:, mask == 0] = fill_value
    return data


def update_mask_with_inps(mask, inps=None, print_msg=True):
    """Update mask matrix from input options: subset_x/y and threshold"""
    if not inps:
        inps = cmd_line_parse()

    if inps.subset_x:
        mask[:, 0:inps.subset_x[0]] = 0
        mask[:, inps.subset_x[1]:] = 0
        if print_msg:
            print('mask out area not in x: {}'.format(inps_dict['subset_x']))

    if inps.subset_y:
        mask[0:inps.subset_y[0], :] = 0
        mask[inps.subset_y[1]:, :] = 0
        if print_msg:
            print('mask out area not in y: {}'.format(inps_dict['subset_y']))

    if inps.threshold:
        mask[mask < inps.threshold] = 0
        if print_msg:
            print('mask out pixels with value < {} in mask file'.format(inps.threshold))
    return mask


def mask_file(fname, mask_file, out_file, inps=None):
    """ Mask input fname with mask_file
    Inputs:
        fname/mask_file - string, 
        inps_dict - dictionary including the following options:
                    subset_x/y - list of 2 ints, subset in x/y direction
                    threshold - float, threshold/minValue to generate mask
    Output:
        out_file - string
    """
    if not inps:
        inps = cmd_line_parse()

    # read mask_file
    mask = readfile.read(mask_file)[0]
    mask = update_mask_with_inps(mask, inps)

    # masking input file
    dsNames = readfile.get_dataset_list(fname)
    maxDigit = max([len(i) for i in dsNames])
    dsDict = {}
    for dsName in dsNames:
        if dsName not in ['coherence']:
            print('masking {d:<{w}} from {f} ...'.format(d=dsName, w=maxDigit, f=fname))
            data = readfile.read(fname, datasetName=dsName, print_msg=False)[0]
            data = mask_matrix(data, mask, fill_value=inps.fill_value)
        dsDict[dsName] = data

    # default output filename
    if not out_file:
        fbase, fext = os.path.splitext(fname)
        out_file = '{}_msk{}'.format(fbase, fext)

    writefile.write(dsDict, out_file=out_file, ref_file=fname)
    return out_file


def mask_isce_file(in_file, mask_file, out_file=None):
    if not in_file:
        return    

    # read mask_file
    print('read mask from {}'.format(mask_file))
    mask = readfile.read(mask_file)[0]

    # mask isce file
    atr = readfile.read_attribute(in_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    interleave = atr['scheme'].upper()
    num_band = int(atr['number_bands'])

    # default short name for data type from ISCE
    dataTypeDict = {
        'byte': 'bool_',
        'float': 'float32',
        'double': 'float64',
        'cfloat': 'complex64',
    }
    data_type = atr['DATA_TYPE'].lower()
    if data_type in dataTypeDict.keys():
        data_type = dataTypeDict[data_type]

    print('read {}'.format(in_file))
    print('setting the (phase) value on the masked out pixels to zero')
    fbase, ext = os.path.splitext(in_file)
    if ext == '.unw':
        amp = readfile.read_binary(in_file, data_type=data_type, num_band=num_band, band_interleave=interleave, band=1)[0]
        pha = readfile.read_binary(in_file, data_type=data_type, num_band=num_band, band_interleave=interleave, band=2)[0]
        pha[mask == 0] = 0
        data = np.hstack((amp, pha)).flatten()
    elif ext == '.int':
        data = np.fromfile(in_file, dtype=data_type, count=length*width).reshape(-1, width)
        #data[mask == 0] = np.abs(data[mask == 0])  #set the angle of complex data to zero
        data[mask == 0] = 0
    elif ext in ['.cor','.conncomp']:
        data = readfile.read(in_file)[0] #, data_type=data_type, num_band=num_band, band_interleave=interleave, band=1)[0]
        data[mask == 0] = 0
    else:
        raise ValueError('unsupported ISCE file: {}'.format(in_file))

    # output filename
    if not out_file:
        if ext in ['.int', '.cor', '.unw']:
            out_file = '{}_msk{}'.format(fbase, ext)
        elif in_file.endswith('.unw.conncomp'):
            out_file = '{}_msk.unw.conncomp'.format(in_file.split('.unw.conncomp')[0])
        else:
            raise ValueError('unrecognized input file type: {}'.format(in_file))

    data.tofile(out_file)
    print('finished writing to file {}'.format(out_file))

    # prepare ISCE metadata file by
    # 1. copy and rename metadata files
    # 2. update file path inside files
    for ext in ['xml', 'vrt']:
        # copy
        cmd = 'cp {i}.{e} {o}.{e}'.format(i=in_file, o=out_file, e=ext)
        os.system(cmd)

        msg = cmd
        msg += ' and update the corresponding filename'
        print(msg)

        # update file path
        meta_file = '{o}.{e}'.format(o=out_file, e=ext)
        with open(meta_file, 'r') as f:
            s = f.read()
        s = s.replace(os.path.basename(in_file),
                      os.path.basename(out_file))
        with open(meta_file, 'w') as f:
            f.write(s)
    return out_file


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    if os.path.isfile(inps.file+'.xml'):
        mask_isce_file(inps.file, inps.mask_file, inps.outfile)
    else:
        mask_file(inps.file, inps.mask_file, inps.outfile, inps)

    print('Done.')
    return


############################################################
if __name__ == '__main__':
    main()
