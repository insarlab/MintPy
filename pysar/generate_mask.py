#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys
import argparse
import h5py
import numpy as np
from pysar.utils import readfile, writefile, utils as ut
from pysar.objects import ifgramDatasetNames


################################################################################################
EXAMPLE='''example:
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  generate_mask.py  081018_090118.unw     -m 3 -M 8 -y 100 700 -x 200 800 -o mask_1.h5
  generate_mask.py  srtm1.dem             -m 0.5 -o maskLand.h5
  generate_mask.py  unwrapIfgram.h5 101120-110220 -m 4
  generate_mask.py  unwrapIfgram.h5 --nonzero

  generate_mask.py ifgramStack.h5 unwrapPhase --nonzero -o maskValidPhase.h5
  generate_mask.py ifgramStack.h5 connectComponent --nonzero -o maskConnComp.h5
'''

def createParser():
    parser = argparse.ArgumentParser(description='Generate mask file from input file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='input file')
    parser.add_argument('dset', nargs='?', help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o','--output', dest='outfile', help='output file name.')

    parser.add_argument('-m','--min', dest='vmin', type=float, help='minimum value for selected pixels')
    parser.add_argument('-M','--max', dest='vmax', type=float, help='maximum value for selected pixels')
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2, metavar=('XMIN','XMAX'), \
                              help='selection range in x/cross-track/range direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2, metavar=('YMIN','YMAX'), \
                              help='selection range in y/along-track/azimuth direction')

    parser.add_argument('--nonzero', dest='nonzero', action='store_true',\
                        help='Select all non-zero pixels.\n'+\
                             'i.e. mask.h5 from unwrapIfgram.h5')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps


def create_threshold_mask(inps):
    if inps.dset:
        print('read %s %s' % (inps.file, inps.dset))
    else:
        print('read %s' % (inps.file))
    data, atr = readfile.read(inps.file, datasetName=inps.dset)
    if len(data.shape) > 2:
        print('ERROR: Only 2D dataset is supported for threshold method, input is 3D')
        sys.exit(1)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    print('create initial mask with the same size as the input file and all = 1')
    mask = np.ones((length, width), dtype=np.float32)

    if inps.nonzero:
        print('all pixels with zero value = 0')
        mask[data == 0] = 0

    # min threshold
    if inps.vmin:
        mask[data<inps.vmin] = 0
        print('all pixels with value < %s = 0' % str(inps.vmin))

    # max threshold
    if inps.vmax:
        mask[data>inps.vmax] = 0
        print('all pixels with value > %s = 0' % str(inps.vmax))

    # nan value
    mask[np.isnan(data)] = 0
    print('all pixels with nan value = 0')

    # subset in Y
    if inps.subset_y:
        y0,y1 = sorted(inps.subset_y)
        mask[0:y0,:] = 0
        mask[y1:length,:] = 0
        print('all pixels with y OUT of [%d, %d] = 0' % (y0,y1))

    # subset in x
    if inps.subset_x:
        x0,x1 = sorted(inps.subset_x)
        mask[:,0:x0] = 0
        mask[:,x1:width] = 0
        print('all pixels with x OUT of [%d, %d] = 0' % (x0,x1))
  
    ## Write mask file
    print('writing >>> '+inps.outfile)
    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, atr, inps.outfile)
    return inps.outfile


################################################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print('input {} file: {}'.format(k, inps.file))

    # default output filename
    if not inps.outfile:
        if inps.file.endswith('temporalCoherence.h5'):
            inps.outfile = 'maskTempCoh.h5'
        else:
            inps.outfile = 'mask.h5'
        if inps.file.startswith('geo_'):
            inps.outfile = 'geo_'+inps.outfile

    ##### Mask: Non-zero
    if inps.nonzero and k == 'ifgramStack':
        inps.outfile = ut.nonzero_mask(inps.file, outFile=inps.outfile, datasetName=inps.dset)
        return inps.outfile

    ##### Mask: Threshold
    inps.outfile = create_threshold_mask(inps)
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()

