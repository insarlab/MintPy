############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Yujie Zheng, Aug 2022         #
############################################################


import os
import sys
import time

from mintpy.utils import arg_utils


################################################################################
REFERENCE = """reference:
  Y. Zheng, H. Fattahi, P. Agram, M. Simons and P. Rosen, (2022). On Closure Phase
    and Systematic Bias in Multi-looked SAR Interferometry, in IEEE Trans. Geosci.
    Remote Sens., doi:10.1109/TGRS.2022.3167648.
"""

EXAMPLE = """example:
  # Note: ONLY sequential network is supported in this implementation.

  # create mask for areas suseptible to biases
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

    # check --water-mask option
    if not inps.water_mask_file and os.path.isfile('./waterMask.h5'):
        inps.water_mask_file = os.path.abspath('./waterMask.h5')

    return inps


################################################################################
def main(iargs=None):
    # parse args
    inps = cmd_line_parse(iargs)

    # import
    import numpy as np
    from ..closure_phase_bias import (
        calc_closure_phase_mask,
        compute_unwrap_closure_phase,
        estimate_bias_timeseries_approx,
        estimate_bias_timeseries,
    )

    # run
    start_time = time.time()

    # common inputs
    kwargs = dict(outdir=inps.outdir, max_memory=inps.maxMemory)

    if inps.action == 'mask':
        calc_closure_phase_mask(
            stack_file=inps.stack_file,
            bias_free_conn=inps.nl,
            num_sigma=inps.num_sigma,
            threshold_amp=inps.epsilon,
            **kwargs)

    elif inps.action.endswith('estimate'):
        # compute the unwrapped closure phase bias time-series
        # and re-unwrap to mitigate the impact of phase unwrapping errors
        # which can dominate the true non-closure phase.
        # max(2, inps.bw) is used to ensure we have conn-2 closure phase processed
        conn_list = np.arange(2, max(2, inps.bw) + 1).tolist() + [inps.nl]
        conn_list = sorted(list(set(conn_list)))
        for conn in conn_list:
            print('\n'+'-'*80)
            print('calculating the unwrapped closure phase for '
                  f'connection level = {conn} out of {conn_list} ...')
            compute_unwrap_closure_phase(
                stack_file=inps.stack_file,
                conn=conn,
                num_worker=int(inps.numWorker),
                **kwargs)

        if inps.action == 'quick_estimate':
            estimate_bias_timeseries_approx(
                stack_file=inps.stack_file,
                bias_free_conn=inps.nl,
                bw=inps.bw,
                water_mask_file=inps.water_mask_file,
                **kwargs)

        elif inps.action == 'estimate':
            cluster_kwargs = {
                "cluster_type" : inps.cluster,
                "num_worker"   : inps.numWorker,
                "config_name"  : inps.config}
            estimate_bias_timeseries(
                stack_file=inps.stack_file,
                bias_free_conn=inps.nl,
                bw=inps.bw,
                cluster_kwargs=cluster_kwargs,
                water_mask_file=inps.water_mask_file,
                **kwargs)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
