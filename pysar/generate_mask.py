#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import sys
import argparse
import numpy as np
from pysar.utils import readfile, writefile, utils as ut, plot as pp


################################################################################################
EXAMPLE = """example:
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5

  # exlcude area by min/max value and/or subset in row/col direction
  generate_mask.py  081018_090118.unw -m 3 -M 8 -y 100 700 -x 200 800 -o mask_1.h5

  # exclude / include an circular area
  generate_mask.py  maskTempCoh.h5 -m 0.5 --ex-circle 230 283 100 -o maskTempCoh_nonDef.h5
  generate_mask.py  maskTempCoh.h5 -m 0.5 --in-circle 230 283 100 -o maskTempCoh_Def.h5

  # use an specific dataset from multiple dataset file
  generate_mask.py  geometryRadar.dem height -m 0.5 -o waterMask.h5
  generate_mask.py  ifgramStack.h5 unwrapPhase-20101120_20110220 -m 4

  # common mask file of pixels without zero unwrapped phase
  generate_mask.py  ifgramStack.h5 unwrapPhase --nonzero

  # interative polygon selection of region of interest
  # useful for custom mask generation in unwrap error correction with bridging
  generate_mask.py  waterMask.h5 -m 0.5 --roipoly
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate mask file from input file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='input file')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name.')

    parser.add_argument('-m', '--min', dest='vmin', type=float,
                        help='minimum value for selected pixels')
    parser.add_argument('-M', '--max', dest='vmax', type=float,
                        help='maximum value for selected pixels')
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2, metavar=('XMIN', 'XMAX'),
                        help='selection range in x/cross-track/range direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2, metavar=('YMIN', 'YMAX'),
                        help='selection range in y/along-track/azimuth direction')
    parser.add_argument('--ex-circle', dest='ex_circle', nargs=3, type=int, metavar=('X', 'Y', 'RADIUS'),
                        help='exclude area defined by an circle (x, y, radius) in pixel number')
    parser.add_argument('--in-circle', dest='in_circle', nargs=3, type=int, metavar=('X', 'Y', 'RADIUS'),
                        help='include area defined by an circle (x, y, radius) in pixel number')

    parser.add_argument('--roipoly', action='store_true',
                        help='Interactive polygonal region of interest (ROI) selection.')

    parser.add_argument('--nonzero', dest='nonzero', action='store_true',
                        help='Select all non-zero pixels.\n' +
                             'i.e. mask.h5 from unwrapIfgram.h5')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################################################
def create_threshold_mask(inps):
    if inps.dset:
        print('read %s %s' % (inps.file, inps.dset))
    else:
        print('read %s' % (inps.file))
    data, atr = readfile.read(inps.file, datasetName=inps.dset)
    if len(data.shape) > 2:
        raise Exception('Only 2D dataset is supported for threshold method, input is 3D')
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    nanmask = ~np.isnan(data)

    print('create initial mask with the same size as the input file and all = 1')
    mask = np.ones((length, width), dtype=np.bool_)

    # nan value
    mask *= nanmask
    print('all pixels with nan value = 0')

    if inps.nonzero:
        print('exclude pixels with zero value')
        mask[nanmask] *= ~(data[nanmask] == 0.)

    # min threshold
    if inps.vmin is not None:
        mask[nanmask] *= ~(data[nanmask] < inps.vmin)
        print('exclude pixels with value < %s' % str(inps.vmin))

    # max threshold
    if inps.vmax is not None:
        mask[nanmask] *= ~(data[nanmask] > inps.vmax)
        print('exclude pixels with value > %s' % str(inps.vmax))

    # subset in Y
    if inps.subset_y is not None:
        y0, y1 = sorted(inps.subset_y)
        mask[0:y0, :] = 0
        mask[y1:length, :] = 0
        print('exclude pixels with y OUT of [%d, %d]' % (y0, y1))

    # subset in x
    if inps.subset_x is not None:
        x0, x1 = sorted(inps.subset_x)
        mask[:, 0:x0] = 0
        mask[:, x1:width] = 0
        print('exclude pixels with x OUT of [%d, %d]' % (x0, x1))

    # exclude circular area
    if inps.ex_circle:
        x, y, r = inps.ex_circle
        cmask = ut.get_circular_mask(x, y, r, length, width)
        mask[cmask == 1] = 0
        print('exclude pixels inside of circle defined as (x={}, y={}, r={})'.format(x, y, r))

    # include circular area
    if inps.in_circle:
        x, y, r = inps.ex_circle
        cmask = ut.get_circular_mask(x, y, r, length, width)
        mask[cmask == 0] = 0
        print('exclude pixels outside of circle defined as (x={}, y={}, r={})'.format(x, y, r))

    # interactively select polygonal region of interest (ROI)
    if inps.roipoly:
        mask *= pp.get_poly_mask(data)

    # Write mask file
    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, out_file=inps.outfile, metadata=atr)
    return inps.outfile


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print('input {} file: {}'.format(k, inps.file))

    # default output filename
    if not inps.outfile:
        if 'temporalCoherence' in inps.file:
            suffix = inps.file.split('temporalCoherence')[1]
            inps.outfile = 'maskTempCoh'+suffix
        else:
            if inps.roipoly:
                inps.outfile = 'maskPoly.h5'
            else:
                inps.outfile = 'mask.h5'
        if inps.file.startswith('geo_'):
            inps.outfile = 'geo_'+inps.outfile

    # default vmin for temporal coherence
    if not inps.vmin and inps.file.endswith('temporalCoherence.h5'):
        inps.vmin = 0.7

    ##### Mask: Non-zero
    if inps.nonzero and k == 'ifgramStack':
        inps.outfile = ut.nonzero_mask(inps.file, out_file=inps.outfile, datasetName=inps.dset)
        return inps.outfile

    ##### Mask: Threshold
    inps.outfile = create_threshold_mask(inps)
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()
