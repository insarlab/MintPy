#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import time
import warnings
import argparse
import h5py
import numpy as np
from mintpy.utils import readfile, writefile, utils as ut


################################################################################################
EXAMPLE = """example:
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5 --base inputs/geometryRadar.h5 --base-dset shadow --base-value 1
  generate_mask.py  avgSpatialCoh.h5     -m 0.7 --base waterMask.h5 -o maskSpatialCoh.h5

  # exlcude area by min/max value and/or subset in row/col direction
  generate_mask.py  081018_090118.unw -m 3 -M 8 -y 100 700 -x 200 800 -o mask_1.h5

  # exclude / include an circular area
  generate_mask.py  maskTempCoh.h5 -m 0.5 --ex-circle 230 283 100 -o maskTempCoh_nonDef.h5
  generate_mask.py  maskTempCoh.h5 -m 0.5 --in-circle 230 283 100 -o maskTempCoh_Def.h5
  # maskout an area within a circle AND with height smaller than a threshold
  generate_mask.py  inputs/geometryGeo.h5 height --in-circle 339 370 21 -M 1400 --revert -o maskCrater.h5

  # use an specific dataset from multiple dataset file
  generate_mask.py  geometryRadar.dem height -m 0.5 -o waterMask.h5
  generate_mask.py  ifgramStack.h5 unwrapPhase-20101120_20110220 -m 4

  # common mask file of pixels in all connected components / with non-zero unwrapped phase
  generate_mask.py  ifgramStack.h5  --nonzero  -o maskConnComp.h5  --update

  # interative polygon selection of region of interest
  # useful for custom mask generation in unwrap error correction with bridging
  generate_mask.py  waterMask.h5 -m 0.5 --roipoly
  generate_mask.py  azOff.h5 --roipoly --view-cmd "-v -0.1 0.1"
  generate_mask.py  velocity.h5 --roipoly --view-cmd "--dem ./inputs/geometryGeo.h5 --contour-step 100 --contour-smooth 0.0"
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
    parser.add_argument('--keep-nan', dest='keep_nan', action='store_true',
                        help='Do not exclude pixels with NaN value')
    parser.add_argument('--revert', action='store_true', help='revert 0 and 1 value of output mask file')

    parser.add_argument('-m', '--min', dest='vmin', type=float,
                        help='minimum value for selected pixels')
    parser.add_argument('-M', '--max', dest='vmax', type=float,
                        help='maximum value for selected pixels')

    aoi = parser.add_argument_group('AOI', 'define secondary area of interest')
    # AOI defined by parameters in command line
    aoi.add_argument('-x', dest='subset_x', type=int, nargs=2, metavar=('XMIN', 'XMAX'),
                     help='selection range in x/cross-track/range direction')
    aoi.add_argument('-y', dest='subset_y', type=int, nargs=2, metavar=('YMIN', 'YMAX'),
                     help='selection range in y/along-track/azimuth direction')
    aoi.add_argument('--ex-circle', dest='ex_circle', nargs=3, type=int, metavar=('X', 'Y', 'RADIUS'),
                     help='exclude area defined by an circle (x, y, radius) in pixel number')
    aoi.add_argument('--in-circle', dest='in_circle', nargs=3, type=int, metavar=('X', 'Y', 'RADIUS'),
                     help='include area defined by an circle (x, y, radius) in pixel number')

    # AOI defined by file
    aoi.add_argument('--base', dest='base_file', type=str,
                     help='exclude pixels == base_value\n'
                          'output_mask[base_data == base_value] = 0')
    aoi.add_argument('--base-dset','--base-dataset', dest='base_dataset', type=str,
                     help='dataset in base_file to be used, for file with multiple datasets.\n'
                          'i.e.: --base inputs/geometryRadar.h5 --base-dset shadow --base-value 1')
    aoi.add_argument('--base-value', dest='base_value', type=float, default=0,
                     help='value of pixels in base_file to be excluded.\nDefault: 0')

    # AOI manual selected
    aoi.add_argument('--roipoly', action='store_true',
                     help='Interactive polygonal region of interest (ROI) selection.')
    aoi.add_argument('--view-cmd', dest='view_cmd', type=str,
                     help='view.py command to facilitate the AOI selection.'
                          'E.g. "-v -0.1 0.1"')

    # special type of mask
    parser.add_argument('--nonzero', dest='nonzero', action='store_true',
                        help='Select all non-zero pixels.\n' +
                             'i.e. maskConnComp.h5 from ifgramStack.h5')

    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update checking for --nonzero option.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check the optional --base and --base-dset options
    if inps.base_file and inps.base_dataset:
        base_dataset_list = readfile.get_dataset_list(inps.base_file)
        if inps.base_dataset not in base_dataset_list:
            msg = 'dataset {} NOT found in input base file {}'.format(inps.base_dataset, inps.base_file)
            msg += '\navailable datasets:\n{}'.format(base_dataset_list)
            warnings.warn(msg)
            print('ignore --base --base-dset option and continue')
            inps.base_file = None
            inps.base_dataset = None

    return inps


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file vs input dataset
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print('1) output file {} NOT exist.'.format(inps.outfile))
    else:
        print('1) output file {} already exists.'.format(inps.outfile))
        with h5py.File(inps.file, 'r') as f:
            ti = float(f[inps.dset].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.file)))
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output file is NOT newer than input dataset: {}.'.format(inps.dset))
        else:
            print('2) output file is newer than input dataset: {}.'.format(inps.dset))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


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
    if not inps.keep_nan:
        mask *= nanmask
        print('all pixels with nan value = 0')

    if inps.nonzero:
        mask[nanmask] *= ~(data[nanmask] == 0.)
        print('exclude pixels with zero value')

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
        cmask = ut.get_circular_mask(x, y, r, (length, width))
        mask[cmask == 1] = 0
        print('exclude pixels inside of circle defined as (x={}, y={}, r={})'.format(x, y, r))

    # include circular area
    if inps.in_circle:
        x, y, r = inps.in_circle
        cmask = ut.get_circular_mask(x, y, r, (length, width))
        mask[cmask == 0] = 0
        print('exclude pixels outside of circle defined as (x={}, y={}, r={})'.format(x, y, r))

    # interactively select polygonal region of interest (ROI)
    if inps.roipoly:
        from mintpy.utils import plot_ext
        poly_mask = plot_ext.get_poly_mask(inps.file, datasetName=inps.dset, view_cmd=inps.view_cmd)
        if poly_mask is not None:
            mask *= poly_mask

    # base mask
    if inps.base_file:
        # read base mask file
        base_data = readfile.read(inps.base_file, datasetName=inps.base_dataset)[0]
        if len(base_data.shape) == 3:
            base_data = np.sum(base_data, axis=0)

        # apply base mask
        mask[base_data == float(inps.base_value)] = 0

        # message
        msg = 'exclude pixels in base file {} '.format(os.path.basename(inps.base_file))
        if inps.base_dataset:
            msg += 'dataset {} '.format(inps.base_dataset)
        msg += 'with value == {}'.format(inps.base_value)
        print(msg)

    # revert
    if inps.revert:
        temp = np.array(mask, dtype=np.bool_)
        mask[temp == True] = False
        mask[temp == False] = True
        del temp

    # Write mask file
    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, out_file=inps.outfile, metadata=atr)
    return inps.outfile


################################################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print('input {} file: {}'.format(k, inps.file))

    # default output filename
    if not inps.outfile:
        if inps.roipoly:
            inps.outfile = 'maskPoly.h5'
        elif 'temporalCoherence' in inps.file:
            suffix = inps.file.split('temporalCoherence')[1]
            inps.outfile = 'maskTempCoh'+suffix
        else:
            inps.outfile = 'mask.h5'
        if inps.file.startswith('geo_'):
            inps.outfile = 'geo_'+inps.outfile

    # default vmin for temporal coherence
    if inps.vmin is None and inps.file.endswith('temporalCoherence.h5'):
        inps.vmin = 0.7

    ##### Mask: Non-zero
    if inps.nonzero and k == 'ifgramStack':
        # get dataset name
        if not inps.dset:
            with h5py.File(inps.file, 'r') as f:
                inps.dset = [i for i in ['connectComponent', 'unwrapPhase'] if i in f.keys()][0]

        # update mode
        if inps.update_mode and inps.outfile and run_or_skip(inps) == 'skip':
            return inps.outfile

        # run
        inps.outfile = ut.nonzero_mask(inps.file, out_file=inps.outfile, datasetName=inps.dset)
        return inps.outfile

    ##### Mask: Threshold
    inps.outfile = create_threshold_mask(inps)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
