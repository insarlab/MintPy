############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
import warnings
from mintpy.utils.arg_utils import create_argument_parser


# from mintpy.objects import RAMP_LIST
RAMP_LIST = [
    'linear', 
    'linear_range',
    'linear_azimuth',
    'quadratic',
    'quadratic_range',
    'quadratic_azimuth',
]


###########################################################################################
EXAMPLE = """example:
  remove_ramp.py  timeseries.h5      -m maskTempCoh.h5
  remove_ramp.py  ifgramStack.h5     -m maskTempCoh.h5  -d unwrapPhase_bridging
  remove_ramp.py  090214_101120.unw  -m maskTempCoh.h5  -s quadratic
"""


def create_parser(subparsers=None):
    synopsis = 'Remove 2D ramp(s) from the input file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='File for ramp removal')
    parser.add_argument('-m', '--mask', dest='mask_file', default='maskTempCoh.h5',
                        help='mask for pixels used in ramp estimation\n'
                             'default - maskTempCoh.h5\n'
                             'no - use the whole area')
    parser.add_argument('-s', dest='surface_type', default='linear', choices=RAMP_LIST,
                        help='type of surface/ramp to remove, linear by default')
    parser.add_argument('-d','--dset', dest='dset',
                        help='dataset name to be derampped in ifgramStack file\n'
                             'e.g.: unwrapPhase\n'
                             '      unwrapPhase_bridging')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('--save-ramp-coeff', dest='save_ramp_coeff', action='store_true',
                        help='Save the estimated ramp coefficients into text file.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if:\n'
                             '1) output file already exists, readable '
                             'and newer than input file\n'
                             '2) all configuration parameters are the same.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # --update requires --outfile
    if inps.update_mode and not inps.outfile:
        inps.update_mode = False
        warnings.warn('update_mode is chosen but NOT turned on because the required --outfile is missing.')
    return inps


###########################################################################################
def main(iargs=None):
    from ..utils import utils as ut
    from ..remove_ramp import run_or_skip, config_keys

    inps = cmd_line_parse(iargs)

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    out_file = ut.run_deramp(
        inps.file,
        ramp_type=inps.surface_type,
        mask_file=inps.mask_file,
        out_file=inps.outfile,
        datasetName=inps.dset,
        save_ramp_coeff=inps.save_ramp_coeff)

    # config parameter
    print('add/update the following configuration metadata to file:\n{}'.format(config_keys))
    atr_new = {}
    atr_new['mintpy.deramp'] = inps.surface_type
    atr_new['mintpy.deramp.maskFile'] = inps.mask_file
    ut.add_attribute(out_file, atr_new)


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
