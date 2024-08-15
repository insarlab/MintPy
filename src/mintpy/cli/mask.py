#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import math
import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################
EXAMPLE = """example:
  mask.py  velocity.h5     -m maskTempCoh.h5
  mask.py  timeseries.h5   -m temporalCoherence.h5  --vmin 0.7
  mask.py  ifgramStack.h5  -m 100102_101120.cor     --vmin 0.9  -y  200 300  -x 300 400

  mask.py  filt_20060924_20090214.int -m waterMask.h5 -o filt_20060924_20090214_msk.int
  mask.py  filt_20060924_20090214.cor -m waterMask.h5 -o filt_20060924_20090214_msk.cor
"""

def create_parser(subparsers=None):
    synopsis = 'Mask file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='File to be masked')
    parser.add_argument('-m', '--mask', dest='mask_file', required=True,
                        help='mask out pixels with mask value == 0.')
    parser.add_argument('-o', '--outfile', help='Output file name.')

    # modify input mask
    parser.add_argument('--vmin','--mask-vmin', dest='mask_vmin', type=float,
                        help='mask out pixels with mask value < vmin.')
    parser.add_argument('--vmax','--mask-vmax', dest='mask_vmax', type=float,
                        help='mask out pixels with mask value > vmax.')
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2,
                        help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2,
                        help='subset range in y/along-track/row direction')
    parser.add_argument('--fill', dest='fill_value', type=float, default=math.nan,
                        help='fill masked out area with input value. i.e. \n'
                             'np.nan (default), 0, 1000, ... \n'
                             'If np.nan and input data matrix is not float/complex, '
                             'convert matrix data type to np.float32.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.mask import mask_file

    # run
    mask_file(
        inps.file,
        mask_file=inps.mask_file,
        out_file=inps.outfile,
        fill_value=inps.fill_value,
        inps=inps,
    )


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
