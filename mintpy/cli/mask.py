############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import math
from mintpy.utils.arg_utils import create_argument_parser


############################################################
EXAMPLE = """example:
  mask.py  velocity.h5     -m Mask.h5
  mask.py  timeseries.h5   -m temporalCoherence.h5  -t 0.7
  mask.py  ifgramStack.h5  -m 100102_101120.cor     -t 0.9  -y  200 300  -x 300 400

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
                        help='mask for pixels used in ramp estimation')
    parser.add_argument('-o', '--outfile', help='Output file name.')

    # modify input mask
    parser.add_argument('-t', dest='threshold', type=float,
                        help='threshold value used for masking.\n' +
                        'if not specified, only pixels with mask value equal to zero is masked out.')
    parser.add_argument('--fill', dest='fill_value', type=float, default=math.nan,
                        help="fill masked out area with input value. i.e. \n"
                             "np.nan (default), 0, 1000, ... \n"
                             "If np.nan and input data matrix is not float/complex, convert matrix data type to np.float32.")
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
def main(iargs=None):
    from ..mask import mask_file
    
    inps = cmd_line_parse(iargs)

    ext = os.path.splitext(inps.file)[1]
    #if os.path.isfile(inps.file+'.xml') and ext in ['.unw','.int','.cor','.conncomp']:
    #    mask_isce_file(inps.file, inps.mask_file, inps.outfile)
    #else:
    mask_file(inps.file, inps.mask_file, inps.outfile, inps)

    print('Done.')


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
