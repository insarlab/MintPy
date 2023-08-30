#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Yujie Zheng, Aug 2022         #
############################################################


import os
import sys

from mintpy.utils import arg_utils

################################################################################
REFERENCE = """reference:
  Y. Zheng, H. Fattahi, P. Agram, M. Simons and P. Rosen, (2022). On Closure Phase
    and Systematic Bias in Multi-looked SAR Interferometry, in IEEE Trans. Geosci.
    Remote Sens., doi:10.1109/TGRS.2022.3167648.
"""

EXAMPLE = """example:
  # Note: ONLY sequential network is supported in this implementation.
  # Notebook tutorial:
  #   https://nbviewer.org/github/insarlab/MintPy-tutorial/blob/main/applications/closure_phase_bias.ipynb

  # create mask for areas susceptible to biases
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 5  -a mask
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 -a mask --num-sigma 2.5

  # estimate non-closure phase bias time-series [quick and approximate solution]
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 5  --bw 3  -a quick_estimate --num-worker 6

  # estimate non-closure phase bias time-series
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 5  --bw 3  -a estimate --num-worker 6 -c local
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 10 -a estimate --num-worker 6 -c local
"""

def create_parser(subparsers=None):
    synopsis = 'Phase non-closure related biases correction'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-i','--ifgramstack', type=str, dest='stack_file',
                        help='interferogram stack file that contains the unwrapped phases')
    parser.add_argument('--wm','--water-mask', dest='water_mask_file',
                        help='Water mask to skip pixels on water body.\n'
                             'Default: waterMask.h5 if exists, otherwise None.')

    # bandwidth and bias free connection level
    parser.add_argument('--nl','--conn-level', dest='nl', type=int, default=20,
                        help='connection level that we are correcting to (or consider as no bias)\n'
                             '(default: %(default)s)')
    parser.add_argument('--bw', dest='bw', type=int, default=10,
                        help='bandwidth of time-series analysis that you want to correct')

    parser.add_argument('-a','--action', dest='action', type=str, default='mask',
                        choices={'mask', 'quick_estimate', 'estimate'},
                        help='action to take (default: %(default)s):\n'
                             'mask           - create a mask of areas susceptible to closure phase errors\n'
                             'quick_estimate - quick and approximate estimation on how bias decays with time\n'
                             '                 output sequential closure phase files\n'
                             'estimate       - estimate how bias decays with time\n'
                             '                 processed for each pixel on a pixel by pixel basis [slow]')

    # mask configuration
    mask = parser.add_argument_group('Mask', 'Configuration for closure phase bias mask')
    mask.add_argument('--num-sigma', dest='num_sigma', type=float, default=3,
                      help='Threashold for phase, in number of sigmas (default: %(default)s).\n'
                           'Assuming a Gaussian distribution for the cumulative closure phase'
                           ' with sigma = pi / sqrt(3*num_cp)')
    mask.add_argument('--eps','--epsilon', dest='epsilon', type=float, default=0.3,
                      help='Threashold for the normalized amplitude in [0-1] (default: %(default)s).')

    # compute
    parser = arg_utils.add_parallel_argument(parser)
    parser = arg_utils.add_memory_argument(parser)

    # output
    parser.add_argument('-o', dest='outdir', type=str, default='./', help='output file directory')

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: --water-mask option
    if not inps.water_mask_file and os.path.isfile('./waterMask.h5'):
        inps.water_mask_file = os.path.abspath('./waterMask.h5')

    return inps


################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.closure_phase_bias import run_closure_phase_bias

    # run
    run_closure_phase_bias(inps)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
