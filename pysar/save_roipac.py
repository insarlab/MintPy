#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import argparse
import numpy as np
from pysar.objects import timeseries
from pysar.utils import readfile, writefile, ptime


##############################################################################
EXAMPLE = """example:
  save_roipac.py  velocity.h5
  save_roipac.py  timeseries.h5    20050601
  save_roipac.py  timeseries.h5    20050601    --ref-date 20040728
  save_roipac.py  INPUTS/ifgramStack.h5  unwrapPhase-20091225_20100723
  save_roipac.py  INPUTS/ifgramStack.h5  unwrapPhase-20091225_20100723  --ref-yx 640 810
  save_roipac.py  INPUTS/ifgramStack.h5    coherence-20091225_20100723
  save_roipac.py  temporal_coherence.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Convert PySAR HDF5 file to ROI_PAC format.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='HDF5 file to be converted.\n' +
                        'for velocity  : the ouput will be a one year interferogram.\n' +
                        'for timeseries: if date is not specified, the last date will be used.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name.')
    parser.add_argument('-r', '--ref-date', dest='ref_date',
                        help='Reference date for timeseries file')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2,
                        help='custom reference pixel in y/x')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if inps.ref_date:
        inps.ref_date = ptime.yyyymmdd(inps.ref_date)
    return inps


def read_data(inps):
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']

    if inps.ref_yx:
        atr['REF_Y'] = inps.ref_yx[0]
        atr['REF_X'] = inps.ref_yx[1]
        print('change reference point to y/x: {}'.format(inps.ref_yx))

    if not inps.dset:
        if k == 'timeseries':
            print('No input date specified >>> continue with the last date')
            inps.dset = timeseries(inps.file).get_date_list()[-1]
        elif k == 'ifgramStack':
            raise Exception("No input dataset! It's required for {} file".format(k))

    print('read {} from file {}'.format(inps.dset, inps.file))
    data = readfile.read(inps.file, datasetName=inps.dset, print_msg=False)[0]
    range2phase = -4 * np.pi / float(atr['WAVELENGTH'])
    if k == 'velocity':
        print("converting velocity to a 1 year interferogram.")
        data *= range2phase
        if inps.ref_yx:
            data -= data[inps.ref_yx[0], inps.ref_yx[1]]
        atr['FILE_TYPE'] = '.unw'
        atr['UNIT'] = 'radian'
        if not inps.outfile:
            inps.outfile = '{}{}'.format(os.path.splitext(inps.file)[0], atr['FILE_TYPE'])

    elif k == 'timeseries':
        if inps.ref_date:
            print('read {} from file {}'.format(inps.ref_date, inps.file))
            data -= readfile.read(inps.file, datasetName=inps.ref_date)[0]
            atr['DATE'] = inps.ref_date[2:8]
            atr['DATE12'] = '{}-{}'.format(inps.ref_date[2:8], inps.dset[2:8])
        else:
            atr['DATE'] = atr['REF_DATE']
            atr['DATE12'] = '{}-{}'.format(atr['REF_DATE'][2:8], inps.dset[2:8])

        data *= range2phase
        if inps.ref_yx:
            data -= data[inps.ref_yx[0], inps.ref_yx[1]]
        atr['FILE_TYPE'] = '.unw'
        atr['UNIT'] = 'radian'
        if not inps.outfile:
            inps.outfile = '{}{}'.format(atr['DATE12'], atr['FILE_TYPE'])

    elif k == 'ifgramStack':
        dsetFamily, atr['DATE12'] = inps.dset.split('-')
        if dsetFamily == 'unwrapPhase':
            if 'REF_X' in atr.keys():
                data -= data[int(atr['REF_Y']), int(atr['REF_X'])]
                print('consider the reference pixel in y/x: ({}, {})'.format(atr['REF_Y'],
                                                                             atr['REF_X']))
            else:
                print('No ref_y/x info found in attributes.')
            atr['FILE_TYPE'] = '.unw'
            atr['UNIT'] = 'radian'
        elif dsetFamily == 'coherence':
            atr['FILE_TYPE'] = '.cor'
            atr['UNIT'] = '1'
        elif dsetFamily == 'wrapPhase':
            atr['FILE_TYPE'] = '.int'
            atr['UNIT'] = 'radian'
        else:
            raise Exception('unrecognized dataset type: {}'.format(inps.dset))

        if not inps.outfile:
            inps.outfile = '{}{}'.format(atr['DATE12'], atr['FILE_TYPE'])

    else:
        if 'coherence' in k.lower():
            atr['FILE_TYPE'] = '.cor'
        elif k in ['mask']:
            atr['FILE_TYPE'] = '.msk'
        elif k in ['geometry'] and inps.dset == 'height':
            if 'Y_FIRST' in atr.keys():
                atr['FILE_TYPE'] = '.dem'
            else:
                atr['FILE_TYPE'] = '.hgt'
        else:
            atr['FILE_TYPE'] = '.unw'

        if not inps.outfile:
            inps.outfile = '{}{}'.format(os.path.splitext(inps.file)[0], atr['FILE_TYPE'])

    if not atr['PROCESSOR'] or atr['PROCESSOR'] == 'pysar':
        atr['PROCESSOR'] = 'roipac'
    return data, atr, inps.outfile


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    data, atr, out_file = read_data(inps)

    print('writing >>> {}'.format(out_file))
    writefile.write(data, out_file=out_file, metadata=atr)
    return inps.outfile


##########################################################################
if __name__ == '__main__':
    main()
