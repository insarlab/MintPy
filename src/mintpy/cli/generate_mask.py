#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys
import warnings

from mintpy.utils.arg_utils import create_argument_parser

################################################################################################
EXAMPLE = """example:
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  generate_mask.py  temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5 --base inputs/geometryRadar.h5 --base-dset shadow --base-value 1
  generate_mask.py  avgSpatialCoh.h5     -m 0.7 --base waterMask.h5 -o maskSpatialCoh.h5

  # exclude area by min/max value and/or subset in row/col direction
  generate_mask.py  081018_090118.unw -m 3 -M 8 -y 100 700 -x 200 800 -o mask_1.h5

  # exclude pixel cluster based on minimum number of pixels
  generate_mask.py  maskTempCoh.h5 -p 10 mask_1.h5

  # exclude pixels with large velocity STD: |velocity| > cutoff (2 by default) * velocityStd
  generate_mask.py  velocity.h5 --vstd
  generate_mask.py  velocity.h5 --vstd --vstd-num 3

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

  # interactive polygon selection of region of interest
  # useful for custom mask generation in unwrap error correction with bridging
  generate_mask.py  waterMask.h5 -m 0.5 --roipoly
  generate_mask.py  azOff.h5 --roipoly --view-cmd "-v -0.1 0.1"
  generate_mask.py  velocity.h5 --roipoly --view-cmd "--dem ./inputs/geometryGeo.h5 --contour-step 100 --contour-smooth 0.0"
"""


def create_parser(subparsers=None):
    synopsis = 'Generate mask file from input file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

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
    parser.add_argument('-p','--mp','--minpixels', dest='minpixels', type=int,
                        help='minimum cluster size in pixels, to remove small pixel clusters.')

    # vmask
    vmask = parser.add_argument_group('AOI for threshold', 'Define the AOI for thresholding operations.')
    vmask.add_argument('--vx', dest='v_subset_x', nargs=2, type=int, metavar=('XMIN', 'XMAX'),
                       help='AOI range in X for threshold operation (and keep the rest untouched.)')
    vmask.add_argument('--vy', dest='v_subset_y', nargs=2, type=int, metavar=('YMIN', 'YMAX'),
                       help='AOI range in Y for threshold operation (and keep the rest untouched.)')
    vmask.add_argument('--vroipoly', action='store_true',
                       help='AOI via interactive polygonal region of interest (ROI) selection.')

    # velocity masking by velocityStd
    parser.add_argument('--vstd', action='store_true',
                        help='mask according to the formula: |velocity| > a * velocityStd')
    parser.add_argument('--vstd-num', dest='vstd_num', type=int, default=2,
                        help='multiple of velocityStd (a) to use for cutoff')

    aoi = parser.add_argument_group('AOI', 'define secondary area of interest')
    # AOI defined by parameters in command line
    aoi.add_argument('-x','--sub-x', dest='subset_x', type=int, nargs=2, metavar=('XMIN', 'XMAX'),
                     help='selection range in x/cross-track/range direction')
    aoi.add_argument('-y','--sub-y', dest='subset_y', type=int, nargs=2, metavar=('YMIN', 'YMAX'),
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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: --base and --base-dset options
    if inps.base_file and inps.base_dataset:
        base_dataset_list = readfile.get_dataset_list(inps.base_file)
        if inps.base_dataset not in base_dataset_list:
            msg = f'dataset {inps.base_dataset} NOT found in input base file {inps.base_file}'
            msg += f'\navailable datasets:\n{base_dataset_list}'
            warnings.warn(msg)
            print('ignore --base --base-dset option and continue')
            inps.base_file = None
            inps.base_dataset = None

    # default: --output
    if not inps.outfile:
        if inps.roipoly:
            inps.outfile = 'maskPoly.h5'
        elif 'temporalCoherence' in inps.file:
            inps.outfile = 'maskTempCoh' + inps.file.split('temporalCoherence')[1]
        else:
            inps.outfile = 'mask.h5'

        # "geo_" prefix
        if inps.file.startswith('geo_'):
            inps.outfile = 'geo_' + inps.outfile

    # default: --min (for temporal coherence)
    if inps.vmin is None and inps.file.endswith('temporalCoherence.h5'):
        inps.vmin = 0.7

    # default: dset (dataset name for non-zero mask from ifgram stack)
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']
    if not inps.dset and inps.nonzero and ftype == 'ifgramStack':
        dset_list = readfile.get_dataset_list(inps.file)
        inps.dset = [i for i in ['connectComponent', 'unwrapPhase'] if i in dset_list][0]

    return inps


################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.generate_mask import create_mask

    # run
    create_mask(inps)


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
