#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_utils

################################################################################################
TEMPLATE = get_template_content('invert_network')

REFERENCE = """references:
  Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002). A new algorithm for surface
    deformation monitoring based on small baseline differential SAR interferograms. IEEE TGRS,
    40(11), 2375-2383. doi:10.1109/TGRS.2002.803792
  Pepe, A., and Lanari, R. (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    of multitemporal differential SAR interferograms, IEEE-TGRS, 44(9), 2374-2383.
  Perissin, D., and Wang, T. (2012), Repeat-pass SAR interferometry with partially coherent targets, IEEE TGRS,
    50(1), 271-280, doi:10.1109/tgrs.2011.2160644.
  Samiei-Esfahany, S., Martins, J. E., Van Leijen, F., and Hanssen, R. F. (2016), Phase Estimation for Distributed
    Scatterers in InSAR Stacks Using Integer Least Squares Estimation, IEEE TGRS, 54(10), 5671-5687.
  Seymour, M. S., and Cumming, I. G. (1994), Maximum likelihood estimation for SAR interferometry, 1994.
    IGARSS '94., 8-12 Aug 1994.
  Yunjun, Z., Fattahi, H., and Amelung, F. (2019), Small baseline InSAR time series analysis: Unwrapping error
    correction and noise reduction, Computers & Geosciences, 133, 104331, doi:10.1016/j.cageo.2019.104331.
  Yunjun, Z., Fattahi, H., Brancato, V., Rosen, P., Simons, M. (2021), Oral: Tectonic displacement mapping from SAR
    offset time series: noise reduction and uncertainty quantification, ID 590, FRINGE 2021, 31 May – 4 Jun, 2021, Virtual.
"""

EXAMPLE = """example:
  ifgram_inversion.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg --update
  ifgram_inversion.py inputs/ifgramStack.h5 -w no  # turn off weight for fast processing
  ifgram_inversion.py inputs/ifgramStack.h5 -c no  # turn off parallel processing
  # offset
  ifgram_inversion.py inputs/ifgramStack.h5 -i rangeOffset   -w no -m waterMask.h5 --md offsetSNR --mt 5
  ifgram_inversion.py inputs/ifgramStack.h5 -i azimuthOffset -w no -m waterMask.h5 --md offsetSNR --mt 5
"""

def create_parser(subparsers=None):
    synopsis = 'Invert network of interferograms into time-series.'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # input dataset
    parser.add_argument('ifgramStackFile', help='interferograms stack file to be inverted')
    parser.add_argument('-t','--template', dest='templateFile', help='template text file with options')

    parser.add_argument('-i','-d', '--dset', dest='obsDatasetName', type=str,
                        help='dataset name of unwrap phase / offset to be used for inversion'
                             '\ne.g.: unwrapPhase, unwrapPhase_bridging, ...')
    parser.add_argument('-m','--water-mask', dest='waterMaskFile',
                        help='Skip inversion on the masked out region, i.e. water.')

    # options rarely used or changed
    parser.add_argument('-o', '--output', dest='outfile', nargs=3,
                        metavar=('TS_FILE', 'TCOH_FILE', 'NUM_INV_FILE'),
                        help='Output file name. (default: %(default)s).')
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--skip-reference','--skip-ref', dest='skip_ref', action='store_true',
                        help='[for offset and testing] do not apply spatial referencing.')

    # solver
    solver = parser.add_argument_group('solver', 'solver for the network inversion problem')
    solver.add_argument('-w', '--weight-func', dest='weightFunc', default='var',
                        choices={'var', 'fim', 'coh', 'no'},
                        help='function used to convert coherence to weight for inversion:\n' +
                             'var - inverse of phase variance due to temporal decorrelation (default)\n' +
                             'fim - Fisher Information Matrix as weight' +
                             'coh - spatial coherence\n' +
                             'no  - no/uniform weight')
    solver.add_argument('--min-norm-phase', dest='minNormVelocity', action='store_false',
                        help=('Enable inversion with minimum-norm deformation phase,'
                              ' instead of the default minimum-norm deformation velocity.'))
    #solver.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1', 'L2'],
    #                    help='Optimization method, L1 or L2 norm. (default: %(default)s).')

    # uncertainty propagation
    parser.add_argument('--calc-cov', dest='calcCov', action='store_true',
                        help='Calculate time-series STD via linear propagation '
                             'from the network of interferograms or offset pairs.')

    # mask
    mask = parser.add_argument_group('mask', 'mask observation data before inversion')
    mask.add_argument('--mask-dset','--mask-dataset','--md', dest='maskDataset',
                      help='dataset used to mask unwrapPhase, e.g. coherence, connectComponent')
    mask.add_argument('--mask-thres','--mask-threshold','--mt', dest='maskThreshold', metavar='NUM', type=float, default=0.4,
                      help='threshold to generate mask when mask is coherence (default: %(default)s).')
    mask.add_argument('--min-redun','--min-redundancy','--mr', dest='minRedundancy', metavar='NUM', type=float, default=1.0,
                      help='minimum redundancy of interferograms for every SAR acquisition. (default: %(default)s).')
    # for offset ONLY
    #mask.add_argument('--mask-min-snr', dest='maskMinSNR', type=float, default=10.0,
    #                  help='minimum SNR to disable/ignore the threshold-based masking [for offset only].')
    #mask.add_argument('--mask-min-area-size', dest='maskMinAreaSize', type=float, default=16.0,
    #                  help='minimum area size to disable/ignore the threshold-based masking [for offset only]')

    # computing
    parser = arg_utils.add_memory_argument(parser)
    parser = arg_utils.add_parallel_argument(parser)

    # update / skip
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n' +
                        'readable and newer than input interferograms file')

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.objects import cluster, ifgramStack
    from mintpy.utils import readfile

    # check
    atr = readfile.read_attribute(inps.ifgramStackFile)

    # check: input file type
    if atr['FILE_TYPE'] not in ['ifgramStack']:
        raise ValueError('input is {} file, support ifgramStack file only.'.format(atr['FILE_TYPE']))

    # check: template file
    if inps.templateFile:
        inps, template = read_template2inps(inps.templateFile, inps)
    else:
        template = dict()

    # check: --cluster and --num-worker option
    inps.numWorker = str(cluster.DaskCluster.format_num_worker(inps.cluster, inps.numWorker))
    if inps.cluster and inps.numWorker == '1':
        print('WARNING: number of workers is 1, turn OFF parallel processing and continue')
        inps.cluster = None

    # default: --dset option
    if not inps.obsDatasetName:
        inps.obsDatasetName = 'unwrapPhase'

        # determine suffix based on unwrapping error correction method
        obs_suffix_map = {
            'bridging'               : '_bridging',
            'phase_closure'          : '_phaseClosure',
            'bridging+phase_closure' : '_bridging_phaseClosure',
        }
        key = 'mintpy.unwrapError.method'
        if key in template.keys() and template[key]:
            unw_err_method = template[key].lower().replace(' ','')   # fix potential typo
            inps.obsDatasetName += obs_suffix_map[unw_err_method]
            print(f'phase unwrapping error correction "{unw_err_method}" is turned ON')
        print(f'use dataset "{inps.obsDatasetName}" by default')

    # check: --dset option (if input observation dataset exists)
    stack_obj = ifgramStack(inps.ifgramStackFile)
    stack_obj.open(print_msg=False)
    if inps.obsDatasetName not in stack_obj.datasetNames:
        msg = f'input dataset name "{inps.obsDatasetName}" not found in file: {inps.ifgramStackFile}'
        raise ValueError(msg)

    # default: --skip-ref option
    if ('offset' in inps.obsDatasetName.lower()
            and 'REF_X' not in atr.keys()
            and 'REF_Y' not in atr.keys()):
        inps.skip_ref = True

    # default: --output option
    if not inps.outfile:
        if inps.obsDatasetName.startswith('unwrapPhase'):
            if os.path.basename(inps.ifgramStackFile).startswith('ion'):
                inps.outfile = ['ion.h5', 'temporalCoherenceIon.h5', 'numInvIon.h5']
            else:
                inps.outfile = ['timeseries.h5', 'temporalCoherence.h5', 'numInvIfgram.h5']

        elif inps.obsDatasetName.startswith('azimuthOffset'):
            inps.outfile = ['timeseriesAz.h5', 'residualInvAz.h5', 'numInvOffAz.h5']

        elif inps.obsDatasetName.startswith('rangeOffset'):
            inps.outfile = ['timeseriesRg.h5', 'residualInvRg.h5', 'numInvOffRg.h5']

        else:
            raise ValueError(f'un-recognized input observation dataset name: {inps.obsDatasetName}')

    # default: --output (split for easy reference)
    inps.tsFile, inps.invQualityFile, inps.numInvFile = inps.outfile

    # default: --water-mask option
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file:', template_file)

    from mintpy.ifgram_inversion import key_prefix
    from mintpy.utils import readfile, utils1 as ut

    iDict = vars(inps)
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(iDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if key in ['weightFunc', 'maskDataset', 'minNormVelocity']:
            iDict[key] = value
        elif value:
            if key in ['maskThreshold', 'minRedundancy']:
                iDict[key] = float(value)
            elif key in ['residualNorm', 'waterMaskFile']:
                iDict[key] = value

    # computing configurations
    dask_key_prefix = 'mintpy.compute.'
    key_list = [i for i in list(iDict.keys()) if dask_key_prefix+i in template.keys()]
    for key in key_list:
        value = template[dask_key_prefix+key]
        if key in ['cluster', 'config']:
            iDict[key] = value
        elif value:
            if key in ['numWorker']:
                iDict[key] = str(value)
            elif key in ['maxMemory']:
                iDict[key] = float(value)

    # False/None --> 'no'
    for key in ['weightFunc']:
        if not iDict[key]:
            iDict[key] = 'no'

    return inps, template


################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.ifgram_inversion import run_ifgram_inversion, run_or_skip

    # run or skip
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return

    # run
    run_ifgram_inversion(inps)


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
