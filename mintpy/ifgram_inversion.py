#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
# Parallel support added by David Grossman, Joshua Zahner  #
############################################################
# Recommend import:
#     from mintpy import ifgram_inversion as ifginv


import os
import sys
import time
import argparse
import warnings
import h5py
import numpy as np
from scipy import linalg   # more effieint than numpy.linalg
from mintpy.objects import ifgramStack, timeseries, cluster
from mintpy.simulation import decorrelation as decor
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, writefile, ptime, utils as ut


# key configuration parameter name
key_prefix = 'mintpy.networkInversion.'
configKeys = ['obsDatasetName',
              'numIfgram',
              'weightFunc',
              'maskDataset',
              'maskThreshold',
              'minRedundancy',
              'minNormVelocity']


################################################################################################
EXAMPLE = """example:
  ifgram_inversion.py  inputs/ifgramStack.h5 -t smallbaselineApp.cfg --update
  ifgram_inversion.py  inputs/ifgramStack.h5 -w no  # turn off weight for fast processing
  ifgram_inversion.py. inputs/ifgramStack.h5 -c no  # turn off parallel processing
  # invert offset stack
  ifgram_inversion.py  inputs/ifgramStack.h5 -i azimuthOffset --water-mask waterMask.h5 --mask-dset offsetSNR --mask-threshold 5
"""

TEMPLATE = get_template_content('invert_network')

REFERENCE = """references:
  Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002). A new algorithm for surface
    deformation monitoring based on small baseline differential SAR interferograms. IEEE TGRS,
    40(11), 2375-2383. doi:10.1109/TGRS.2002.803792
  Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    of multitemporal differential SAR interferograms, IEEE-TGRS, 44(9), 2374-2383.
  Perissin, D., and T. Wang (2012), Repeat-pass SAR interferometry with partially coherent targets, IEEE TGRS,
    50(1), 271-280, doi:10.1109/tgrs.2011.2160644.
  Samiei-Esfahany, S., J. E. Martins, F. v. Leijen, and R. F. Hanssen (2016), Phase Estimation for Distributed
    Scatterers in InSAR Stacks Using Integer Least Squares Estimation, IEEE TGRS, 54(10), 5671-5687.
  Seymour, M. S., and I. G. Cumming (1994), Maximum likelihood estimation for SAR interferometry, 1994.
    IGARSS '94., 8-12 Aug 1994.
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis: Unwrapping error
    correction and noise reduction, Computers & Geosciences, 133, 104331, doi:10.1016/j.cageo.2019.104331.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Invert network of interferograms into time-series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+TEMPLATE+'\n'+EXAMPLE)
    # input dataset
    parser.add_argument('ifgramStackFile', help='interferograms stack file to be inverted')
    parser.add_argument('--template', '-t', dest='templateFile', help='template text file with options')

    parser.add_argument('-i','-d', '--dset', dest='obsDatasetName', type=str,
                        help='dataset name of unwrap phase / offset to be used for inversion'
                             '\ne.g.: unwrapPhase, unwrapPhase_bridging, ...')
    parser.add_argument('--water-mask', '-m', dest='waterMaskFile',
                        help='Skip inversion on the masked out region, i.e. water.')

    # options rarely used or changed
    parser.add_argument('-o', '--output', dest='outfile', nargs=2, metavar=('TS_FILE', 'TCOH_FILE'),
                        default=['timeseries.h5', 'temporalCoherence.h5', 'numInvIfgram.h5'],
                        help='Output file name. (default: %(default)s).')
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',
                        help='Skip checking reference pixel value, for simulation testing.')

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
    solver.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1', 'L2'],
                        help='Optimization mehtod, L1 or L2 norm. (default: %(default)s).')

    # mask
    mask = parser.add_argument_group('mask', 'mask observation data before inversion')
    mask.add_argument('--mask-dset', dest='maskDataset',
                      help='dataset used to mask unwrapPhase, e.g. coherence, connectComponent')
    mask.add_argument('--mask-threshold', dest='maskThreshold', metavar='NUM', type=float, default=0.4,
                      help='threshold to generate mask when mask is coherence (default: %(default)s).')
    mask.add_argument('--min-redundancy', dest='minRedundancy', metavar='NUM', type=float, default=1.0,
                      help='minimum redundancy of interferograms for every SAR acquisition. (default: %(default)s).')

    # computing
    parser.add_argument('-r', '--ram', '--memory', dest='memorySize', type=float, default=2,
                        help='Max amount of memory in GB to use (default: %(default)s).\n' +
                             'Adjust according to your computer memory.')

    par = parser.add_argument_group('parallel', 'parallel processing using dask')
    par.add_argument('-c', '--cluster', '--cluster-type', dest='cluster', type=str,
                     default='local', choices=cluster.CLUSTER_LIST + ['no'],
                     help='Cluster to use for parallel computing, no to turn OFF. (default: %(default)s).')
    par.add_argument('--num-worker', dest='numWorker', type=str, default='4',
                     help='Number of workers to use (default: %(default)s).')
    par.add_argument('--config', '--config-name', dest='config', type=str, default=None,
                     help='Configuration name to use in dask.yaml (default: %(default)s).')

    # update / skip
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n' +
                        'readable and newer than input interferograms file')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.ifgramStackFile)
    if atr['FILE_TYPE'] not in ['ifgramStack']:
        raise ValueError('input is {} file, support ifgramStack file only.'.format(atr['FILE_TYPE']))

    if inps.templateFile:
        inps, template = read_template2inps(inps.templateFile, inps)
    else:
        template = dict()

    # --cluster and --num-worker option
    inps.numWorker = str(cluster.DaskCluster.format_num_worker(inps.cluster, inps.numWorker))
    if inps.cluster != 'no' and inps.numWorker == '1':
        print('WARNING: number of workers is 1, turn OFF parallel processing and continue')
        inps.cluster = 'no'

    # --water-mask option
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    # --dset option
    if not inps.obsDatasetName:
        inps.obsDatasetName = 'unwrapPhase'

        # determine suffix based on unwrapping error correction method
        obs_suffix_map = {'bridging'               : '_bridging',
                          'phase_closure'          : '_phaseClosure',
                          'bridging+phase_closure' : '_bridging_phaseClosure'}
        key = 'mintpy.unwrapError.method'
        if key in template.keys() and template[key]:
            unw_err_method = template[key].lower().replace(' ','')   # fix potential typo
            inps.obsDatasetName += obs_suffix_map[unw_err_method]
            print('phase unwrapping error correction "{}" is turned ON'.format(unw_err_method))
        print('use dataset "{}" by default'.format(inps.obsDatasetName))

        # check if input observation dataset exists.
        stack_obj = ifgramStack(inps.ifgramStackFile)
        stack_obj.open(print_msg=False)
        if inps.obsDatasetName not in stack_obj.datasetNames:
            msg = 'input dataset name "{}" not found in file: {}'.format(inps.obsDatasetName, inps.ifgramStackFile)
            raise ValueError(msg)

    # --skip-ref option
    if 'offset' in inps.obsDatasetName.lower():
        inps.skip_ref = True

    # --output option
    if not inps.outfile:
        if inps.obsDatasetName.startswith('unwrapPhase'):
            inps.outfile = ['timeseries.h5', 'temporalCoherence.h5', 'numInvIfgram.h5']

        elif inps.obsDatasetName.startswith('azimuthOffset'):
            inps.outfile = ['timeseriesAz.h5', 'temporalCoherenceAz.h5', 'numInvOffset.h5']

        elif inps.obsDatasetName.startswith('rangeOffset'):
            inps.outfile = ['timeseriesRg.h5', 'temporalCoherenceRg.h5', 'numInvOffset.h5']

        else:
            raise ValueError('un-recognized input observation dataset name: {}'.format(inps.obsDatasetName))

    inps.tsFile, inps.tempCohFile, inps.numInvFile = inps.outfile

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    iDict = vars(inps)

    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)
    keyList = [i for i in list(iDict.keys()) if key_prefix+i in template.keys()]
    for key in keyList:
        value = template[key_prefix+key]
        if key in ['weightFunc', 'maskDataset', 'minNormVelocity']:
            iDict[key] = value
        elif value:
            if key in ['maskThreshold', 'minRedundancy', 'memorySize']:
                iDict[key] = float(value)
            elif key in ['residualNorm', 'waterMaskFile']:
                iDict[key] = value

    # computing configurations
    dask_key_prefix = 'mintpy.compute.'
    keyList = [i for i in list(iDict.keys()) if dask_key_prefix+i in template.keys()]
    for key in keyList:
        value = template[dask_key_prefix+key]
        if key in ['cluster', 'config']:
            iDict[key] = value
        elif value:
            if key in ['numWorker']:
                iDict[key] = str(value)
            elif key in ['memorySize']:
                iDict[key] = float(value)

    # False/None --> 'no'
    for key in ['weightFunc', 'cluster']:
        if not iDict[key]:
            iDict[key] = 'no'

    return inps, template


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output files vs input dataset
    if not all(os.path.isfile(i) for i in inps.outfile):
        flag = 'run'
        print('1) NOT ALL output files found: {}.'.format(inps.outfile))
    else:
        print('1) output files already exist: {}.'.format(inps.outfile))
        with h5py.File(inps.ifgramStackFile, 'r') as f:
            ti = float(f[inps.obsDatasetName].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgramStackFile)))
        to = min(os.path.getmtime(i) for i in inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output files are NOT newer than input dataset: {}.'.format(inps.obsDatasetName))
        else:
            print('2) output dataset is newer than input dataset: {}.'.format(inps.obsDatasetName))

    # check configuration
    if flag == 'skip':
        meta_keys = ['REF_Y', 'REF_X']
        atr_ifg = readfile.read_attribute(inps.ifgramStackFile)
        atr_ts = readfile.read_attribute(inps.tsFile)
        inps.numIfgram = len(ifgramStack(inps.ifgramStackFile).get_date12_list(dropIfgram=True))

        if any(str(vars(inps)[key]) != atr_ts.get(key_prefix+key, 'None') for key in configKeys):
            flag = 'run'
            print('3) NOT all key configration parameters are the same: {}'.format(configKeys))
        elif any(atr_ts[key] != atr_ifg[key] for key in meta_keys):
            flag = 'run'
            print('3) NOT all the metadata are the same: {}'.format(meta_keys))
        else:
            print('3) all key configuration parameters are the same: {}.'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


################################# Time-series Estimator ###################################
def estimate_timeseries(A, B, tbase_diff, ifgram, weight_sqrt=None, min_norm_velocity=True,
                        rcond=1e-5, min_redundancy=1., skip_zero_value=True):
    """Estimate time-series from a stack/network of interferograms with
    Least Square minimization on deformation phase / velocity.

    opt 1: X = np.dot(np.dot(numpy.linalg.inv(np.dot(B.T, B)), B.T), ifgram)
    opt 2: X = np.dot(numpy.linalg.pinv(B), ifgram)
    opt 3: X = np.dot(scipy.linalg.pinv2(B), ifgram)
    opt 4: X = scipy.linalg.lstsq(B, ifgram)[0] [recommend and used]

    opt 4 supports weight.
    scipy.linalg provides more advanced and slighted faster performance than numpy.linalg.
    This function relies on the LAPACK routine gelsd. It computes the minimum-norm
    solution to a linear least squares problem using the singular value decomposition
    of A and a divide and conquer method.

    opt 4 is faster than opt 1/2/3 because it estimates X directly without calculating
    the A_inv matrix.

    opt 2/3 is better than opt 1 because numpy.linalg.inv() can not handle rank defiency of
    design matrix B

    Traditional Small BAseline Subsets (SBAS) algorithm (Berardino et al., 2002, IEEE-TGRS)
    is equivalent to the setting of:
        min_norm_velocity=True
        weight_sqrt=None

    Parameters: A                 - 2D np.array in size of (num_ifgram, num_date-1)
                B                 - 2D np.array in size of (num_ifgram, num_date-1),
                                    design matrix B, each row represents differential temporal
                                    baseline history between reference and secondary date of one interferogram
                tbase_diff        - 2D np.array in size of (num_date-1, 1),
                                    differential temporal baseline history
                ifgram            - 2D np.array in size of (num_ifgram, num_pixel),
                                    phase of all interferograms
                weight_sqrt       - 2D np.array in size of (num_ifgram, num_pixel),
                                    square root of weight of all interferograms
                min_norm_velocity - bool, assume minimum-norm deformation velocity, or not
                rcond             - cut-off ratio of small singular values of A or B, to maintain robustness.
                                    It's recommend to >= 1e-5 by experience, to generate reasonable result.
                min_redundancy    - min redundancy defined as min num_ifgram for every SAR acquisition
    Returns:    ts                - 2D np.array in size of (num_date, num_pixel), phase time-series
                temp_coh          - 1D np.array in size of (num_pixel), temporal coherence
                num_inv_obs       - 1D np.array in size of (num_pixel), number of observatons (ifgrams / offsets)
                                    used during the inversion
    """
    ifgram = ifgram.reshape(A.shape[0], -1)
    if weight_sqrt is not None:
        weight_sqrt = weight_sqrt.reshape(A.shape[0], -1)
    num_date = A.shape[1] + 1
    num_pixel = ifgram.shape[1]

    # initial output value
    ts = np.zeros((num_date, num_pixel), dtype=np.float32)
    temp_coh = 0.
    num_inv_obs = 0

    # skip zero phase value
    if skip_zero_value and not np.all(ifgram):
        idx = (ifgram[:, 0] != 0.).flatten()
        A = A[idx, :]
        B = B[idx, :]

        # skip the pixel if its redundancy < threshold
        if np.min(np.sum(A != 0., axis=0)) < min_redundancy:
            return ts, temp_coh, num_inv_obs

        # check matrix invertability
        if weight_sqrt is not None:  #for WLS only because OLS contains it already
            try:
                linalg.inv(np.dot(B.T, B))
            except linalg.LinAlgError:
                return ts, temp_coh, num_inv_obs

        ifgram = ifgram[idx, :]
        if weight_sqrt is not None:
            weight_sqrt = weight_sqrt[idx, :]

    # update number of observations used for inversion
    num_inv_obs = A.shape[0]

    # invert time-series
    try:
        # assume minimum-norm deformation velocity
        if min_norm_velocity:
            if weight_sqrt is not None:
                X = linalg.lstsq(np.multiply(B, weight_sqrt),
                                 np.multiply(ifgram, weight_sqrt),
                                 cond=rcond)[0]
            else:
                X = linalg.lstsq(B, ifgram, cond=rcond)[0]

            # calc temporal coherence
            temp_coh = calc_temporal_coherence(ifgram, B, X)

            # assemble time-series
            ts_diff = X * np.tile(tbase_diff, (1, num_pixel))
            ts[1:, :] = np.cumsum(ts_diff, axis=0)

        # assume minimum-norm deformation phase
        else:
            if weight_sqrt is not None:
                X = linalg.lstsq(np.multiply(A, weight_sqrt),
                                 np.multiply(ifgram, weight_sqrt),
                                 cond=rcond)[0]
            else:
                X = linalg.lstsq(A, ifgram, cond=rcond)[0]

            # calc temporal coherence
            temp_coh = calc_temporal_coherence(ifgram, A, X)

            # assemble time-series
            ts[1: ,:] = X

    except linalg.LinAlgError:
        pass

    return ts, temp_coh, num_inv_obs


def calc_temporal_coherence(ifgram, G, X):
    """Calculate the temporal coherence from the network inversion results

    Parameters: ifgram   - 2D np.array in size of (num_ifgram, num_pixel), phase or offset
                G        - 2D np.array in size of (num_ifgram, num_date-1), design matrix A or B
                X        - 2D np.array in size of (num_date-1, num_pixel), solution
    Returns:    temp_coh - 1D np.array in size of (num_pixel), temporal coherence
    """

    num_ifgram, num_pixel = ifgram.shape
    temp_coh = np.zeros(num_pixel, dtype=np.float32)

    # chunk_size as the number of pixels
    chunk_size = int(ut.round_to_1(2e5 / num_ifgram))
    if num_pixel > chunk_size:
        num_chunk = int(np.ceil(num_pixel / chunk_size))
        num_chunk_step = int(ut.round_to_1(num_chunk / 5))
        print(('calcualting temporal coherence in chunks of {} pixels'
               ': {} chunks in total ...').format(chunk_size, num_chunk))

        for i in range(num_chunk):
            c0 = i * chunk_size
            c1 = min((i + 1) * chunk_size, num_pixel)

            # calc residual
            ifgram_diff = ifgram[:, c0:c1] - np.dot(G, X[:, c0:c1])

            # calc temporal coherence
            temp_coh[c0:c1] = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / num_ifgram

            # print out message
            if (i+1) % num_chunk_step == 0:
                print('chunk {} / {}'.format(i+1, num_chunk))

    else:
        # calc residual
        ifgram_diff = ifgram - np.dot(G, X)

        # calc temporal coherence
        temp_coh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / num_ifgram

    return temp_coh



###################################### File IO ############################################
def write2hdf5_file(ifgram_file, metadata, ts, temp_coh, num_inv_ifg=None,
                    suffix='', inps=None):
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    date_list = stack_obj.get_date_list(dropIfgram=True)

    # File 1 - timeseries.h5
    ts_file = '{}{}.h5'.format(suffix, os.path.splitext(inps.outfile[0])[0])
    metadata['REF_DATE'] = date_list[0]
    metadata['FILE_TYPE'] = 'timeseries'
    metadata['UNIT'] = 'm'

    print('-'*50)
    print('calculating perpendicular baseline timeseries')
    pbase = stack_obj.get_perp_baseline_timeseries(dropIfgram=True)

    ts_obj = timeseries(ts_file)
    ts_obj.write2hdf5(data=ts, dates=date_list, bperp=pbase, metadata=metadata)

    # File 2 - temporalCoherence.h5
    out_file = '{}{}.h5'.format(suffix, os.path.splitext(inps.outfile[1])[0])
    metadata['FILE_TYPE'] = 'temporalCoherence'
    metadata['UNIT'] = '1'
    print('-'*50)
    writefile.write(temp_coh, out_file=out_file, metadata=metadata)

    ## File 3 - timeseriesDecorStd.h5
    #if not np.all(ts_std == 0.):
    #    out_file = 'timeseriesDecorStd{}.h5'.format(suffix)
    #    metadata['FILE_TYPE'] = 'timeseries'
    #    metadata['UNIT'] = 'm'
    #    phase2range = -1*float(stack_obj.metadata['WAVELENGTH'])/(4.*np.pi)
    #    ts_std *= abs(phase2range)
    #    print('-'*50)
    #    writefile.write(ts_std, out_file=out_file, metadata=metadata, ref_file=ts_file)

    # File 3 - numInvIfgram.h5
    out_file = 'numInvIfgram{}.h5'.format(suffix)
    metadata['FILE_TYPE'] = 'mask'
    metadata['UNIT'] = '1'
    print('-'*50)
    writefile.write(num_inv_ifg, out_file=out_file, metadata=metadata)
    return


def split2boxes(ifgram_file, memory_size=4, print_msg=True):
    """Split into chunks in rows to reduce memory usage
    Parameters: dataset_shape - tuple of 3 int
                memory_size   - float, max memory to use in GB
                print_msg     - bool
    Returns:    box_list      - list of tuple of 4 int
                num_box       - int, number of boxes
    """
    ifg_obj = ifgramStack(ifgram_file)
    ifg_obj.open(print_msg=False)

    # 1st dimension size: defo obs (phase / offset) + weight + time-series
    num_epoch = ifg_obj.numIfgram * 2 + ifg_obj.numDate + 5
    length, width = ifg_obj.length, ifg_obj.width

    # split in lines based on the input memory limit
    y_step = (memory_size * (1e3**3)) / (num_epoch * width * 4)

    # calibrate based on experience
    y_step = int(ut.round_to_1(y_step * 0.6))

    num_box = int((length - 1) / y_step) + 1
    if print_msg and num_box > 1:
        print('maximum memory size: %.1E GB' % memory_size)
        print('split %d lines into %d patches for processing' % (length, num_box))
        print('    with each patch up to %d lines' % y_step)

    # y_step / num_box --> box_list
    box_list = []
    for i in range(num_box):
        y0 = i * y_step
        y1 = min([length, y0 + y_step])
        box = (0, y0, width, y1)
        box_list.append(box)

    return box_list, num_box


def check_design_matrix(ifgram_file, weight_func='var'):
    """
    Check Rank of Design matrix for weighted inversion
    """

    date12_list = ifgramStack(ifgram_file).get_date12_list(dropIfgram=True)
    A = ifgramStack.get_design_matrix4timeseries(date12_list)[0]
    if weight_func == 'no':
        if np.linalg.matrix_rank(A) < A.shape[1]:
            print('WARNING: singular design matrix! Inversion result can be biased!')
            print('continue using its SVD solution on all pixels')
    else:
        if np.linalg.matrix_rank(A) < A.shape[1]:
            print('ERROR: singular design matrix!')
            print('    Input network of interferograms is not fully connected!')
            print('    Can not invert the weighted least square solution.')
            print('You could try:')
            print('    1) Add more interferograms to make the network fully connected:')
            print('       a.k.a., no multiple subsets nor network islands')
            print("    2) Use '-w no' option for non-weighted SVD solution.")
            raise Exception()
    return A


def read_unwrap_phase(stack_obj, box, ref_phase, obs_ds_name='unwrapPhase', dropIfgram=True,
                      print_msg=True):
    """Read unwrapPhase from ifgramStack file

    Parameters: stack_obj - ifgramStack object
                box       - tuple of 4 int
                ref_phase - 1D array or None
    Returns:    pha_data  - 2D array of unwrapPhase in size of (num_ifgram, num_pixel)
    """
    # Read unwrapPhase
    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if print_msg:
        print('reading {} in {} * {} ...'.format(obs_ds_name, box, num_ifgram))
    pha_data = stack_obj.read(datasetName=obs_ds_name,
                              box=box,
                              dropIfgram=dropIfgram,
                              print_msg=False).reshape(num_ifgram, -1)
    pha_data[np.isnan(pha_data)] = 0.

    # read ref_phase
    if ref_phase is not None:
        # use input ref_phase array
        if print_msg:
            print('use input reference phase')

    elif 'refPhase' in stack_obj.datasetNames:
        # read refPhase from file itself
        if print_msg:
            print('read reference phase from file')
        with h5py.File(stack_obj.file, 'r') as f:
            ref_phase = f['refPhase'][:]

    else:
        raise Exception('No reference phase input/found on file!'+
                        ' unwrapped phase is not referenced!')

    # reference unwrapPhase
    for i in range(num_ifgram):
        mask = pha_data[i, :] != 0.
        pha_data[i, :][mask] -= ref_phase[i]
    return pha_data


def mask_unwrap_phase(pha_data, stack_obj, box, mask_ds_name=None, mask_threshold=0.4,
                      dropIfgram=True, print_msg=True):
    """
    Mask input unwrapped phase.
    """

    # Read/Generate Mask
    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if mask_ds_name and mask_ds_name in stack_obj.datasetNames:
        if print_msg:
            print('reading {} in {} * {} ...'.format(mask_ds_name, box, num_ifgram))
        msk_data = stack_obj.read(datasetName=mask_ds_name,
                                  box=box,
                                  dropIfgram=dropIfgram,
                                  print_msg=False).reshape(num_ifgram, -1)
        msk_data[np.isnan(msk_data)] = 0
        if mask_ds_name in ['coherence', 'offsetSNR']:
            msk_data = msk_data >= mask_threshold
            if print_msg:
                print('mask out pixels with {} < {}'.format(mask_ds_name, mask_threshold))
        elif mask_ds_name in ['connectComponent']:
            if print_msg:
                print('mask out pixels with {} == 0'.format(mask_ds_name))
        pha_data[msk_data == 0.] = 0.
        del msk_data
    return pha_data


def read_coherence(stack_obj, box, dropIfgram=True, print_msg=True):
    """
    Read spatial coherence
    """

    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if print_msg:
        print('reading coherence in {} * {} ...'.format(box, num_ifgram))
    coh_data = stack_obj.read(datasetName='coherence',
                              box=box,
                              dropIfgram=dropIfgram,
                              print_msg=False).reshape(num_ifgram, -1)
    coh_data[np.isnan(coh_data)] = 0.
    return coh_data


def calc_weight(stack_obj, box, weight_func='var', dropIfgram=True, chunk_size=100000):
    """Read coherence and calculate weight from it, chunk by chunk to save memory
    """

    print('calculating weight from spatial coherence ...')

    # read coherence
    weight = read_coherence(stack_obj, box=box, dropIfgram=dropIfgram)
    num_pixel = weight.shape[1]

    if 'NCORRLOOKS' in stack_obj.metadata.keys():
        L = np.rint(float(stack_obj.metadata['NCORRLOOKS'])).astype(np.int16)
    else:
        # use the typical ratio of resolution vs pixel size of Sentinel-1 IW mode
        L = int(stack_obj.metadata['ALOOKS']) * int(stack_obj.metadata['RLOOKS'])
        L = np.rint(L / 1.94)
    # make sure L >= 1
    L = max(L, 1)

    # convert coherence to weight chunk-by-chunk to save memory
    num_chunk = int(np.ceil(num_pixel / chunk_size))
    print(('convert coherence to weight in chunks of {c} pixels'
           ': {n} chunks in total ...').format(c=chunk_size, n=num_chunk))

    for i in range(num_chunk):
        c0 = i * chunk_size
        c1 = min((i + 1) * chunk_size, num_pixel)
        if i == 0:
            print_msg = True
        else:
            print_msg = False

        # calc weight from coherence
        weight[:, c0:c1] = decor.coherence2weight(weight[:, c0:c1],
                                                  weight_func,
                                                  L=L,
                                                  epsilon=5e-2,
                                                  print_msg=print_msg)
        weight[:, c0:c1] = np.sqrt(weight[:, c0:c1])

        # print out message
        if (i+1) % 1 == 0:
            print('chunk {} / {}'.format(i+1, num_chunk))

    return weight


def ifgram_inversion_patch(ifgram_file, box=None, ref_phase=None, obs_ds_name='unwrapPhase',
                           weight_func='var', water_mask_file=None, min_norm_velocity=True,
                           mask_ds_name=None, mask_threshold=0.4, min_redundancy=1.0):
    """Invert one patch of an ifgram stack into timeseries.

    Parameters: box               - tuple of 4 int, indicating (x0, y0, x1, y1) of the area of interest
                                    or None for the whole image
                ifgram_file       - str, interferograms stack HDF5 file, e.g. ./inputs/ifgramStack.h5
                ref_phase         - 1D array in size of (num_ifgram), or None
                obs_ds_name       - str, dataset to feed the inversion.
                weight_func       - str, weight function, choose in ['no', 'fim', 'var', 'coh']
                water_mask_file   - str, water mask filename if available, to skip inversion on water
                min_norm_velocity - bool, minimize the residual phase or phase velocity
                mask_ds_name      - str, dataset name in ifgram_file used to mask unwrapPhase pixelwisely
                mask_threshold    - float, min coherence of pixels if mask_dataset_name='coherence'
                min_redundancy    - float, the min number of ifgrams for every acquisition.
    Returns:    ts                - 3D array in size of (num_date, num_row, num_col)
                temp_coh          - 2D array in size of (num_row, num_col)
                num_inv_ifg       - 2D array in size of (num_row, num_col)
                box               - tuple of 4 int
    Example:    ifgram_inversion_patch('ifgramStack.h5', box=(0,200,1316,400))
    """

    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)

    # debug
    #y, x = 258, 454
    #box = (x, y, x+1, y+1)

    ## 1. input info

    # size
    if box:
        num_row = box[3] - box[1]
        num_col = box[2] - box[0]
    else:
        num_row = stack_obj.length
        num_col = stack_obj.width
    num_pixel = num_row * num_col

    # get tbase_diff in the unit of year
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)
    tbase = np.array(ptime.date_list2tbase(date_list)[0], np.float32) / 365.25
    tbase_diff = np.diff(tbase).reshape(-1, 1)

    # design matrix
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    A, B = stack_obj.get_design_matrix4timeseries(date12_list=date12_list)[0:2]

    # prep for decor std time-series
    #if os.path.isfile('reference_date.txt'):
    #    ref_date = str(np.loadtxt('reference_date.txt', dtype=bytes).astype(str))
    #else:
    #    ref_date = date_list[0]
    #Astd = stack_obj.get_design_matrix4timeseries(date12_list=date12_list, refDate=ref_date)[0]
    #ref_idx = date_list.index(ref_date)
    #time_idx = [i for i in range(num_date)]
    #time_idx.remove(ref_idx)

    # skip zero value in the network inversion for phase
    if 'phase' in obs_ds_name.lower():
        skip_zero_value = True
    else:
        skip_zero_value = False

    # 1.1 read / calcualte weight
    if weight_func in ['no', 'sbas']:
        weight = None
    else:
        weight = calc_weight(stack_obj,
                             box,
                             weight_func=weight_func,
                             dropIfgram=True,
                             chunk_size=100000)

    # 1.2 read / mask unwrapPhase / offset
    pha_data = read_unwrap_phase(stack_obj,
                                 box,
                                 ref_phase,
                                 obs_ds_name=obs_ds_name,
                                 dropIfgram=True)

    pha_data = mask_unwrap_phase(pha_data,
                                 stack_obj,
                                 box,
                                 dropIfgram=True,
                                 mask_ds_name=mask_ds_name,
                                 mask_threshold=mask_threshold)

    # 1.3 mask of pixels to invert
    mask = np.ones(num_pixel, np.bool_)

    # 1.3.1 - Water Mask
    if water_mask_file:
        print('skip pixels on water with mask from file: {}'.format(os.path.basename(water_mask_file)))
        atr_msk = readfile.read_attribute(water_mask_file)
        len_msk, wid_msk = int(atr_msk['LENGTH']), int(atr_msk['WIDTH'])
        if (len_msk, wid_msk) != (stack_obj.length, stack_obj.width):
            raise ValueError('Input water mask file has different size from ifgramStack file.')

        dsName = [i for i in readfile.get_dataset_list(water_mask_file) if i in ['waterMask', 'mask']][0]
        waterMask = readfile.read(water_mask_file, datasetName=dsName, box=box)[0].flatten()
        mask *= np.array(waterMask, dtype=np.bool_)
        del waterMask

    # 1.3.2 - Mask for Zero Phase in ALL ifgrams
    if 'phase' in obs_ds_name.lower():
        print('skip pixels with zero/nan value in all interferograms')
        with warnings.catch_warnings():
            # ignore warning message for all-NaN slices
            warnings.simplefilter("ignore", category=RuntimeWarning)
            phase_stack = np.nanmean(pha_data, axis=0)
        mask *= np.multiply(~np.isnan(phase_stack), phase_stack != 0.)
        del phase_stack

    # 1.3.3 invert pixels on mask 1+2
    num_pixel2inv = int(np.sum(mask))
    idx_pixel2inv = np.where(mask)[0]
    print('number of pixels to invert: {} out of {} ({:.1f}%)'.format(
        num_pixel2inv, num_pixel, num_pixel2inv/num_pixel*100))

    ## 2. inversion

    # 2.1 initiale the output matrices
    ts = np.zeros((num_date, num_pixel), np.float32)
    #ts_std = np.zeros((num_date, num_pixel), np.float32)
    temp_coh    = np.zeros(num_pixel, np.float32)
    num_inv_ifg = np.zeros(num_pixel, np.int16)

    # return directly if there is nothing to invert
    if num_pixel2inv < 1:
        ts = ts.reshape(num_date, num_row, num_col)
        #ts_std = ts_std.reshape(num_date, num_row, num_col)
        temp_coh = temp_coh.reshape(num_row, num_col)
        num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)
        return ts, temp_coh, num_inv_ifg, box

    # 2.2 un-weighted inversion (classic SBAS)
    if weight_func in ['no', 'sbas']:

        # a. mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        if 'phase' in obs_ds_name.lower():
            mask_all_net = np.all(pha_data, axis=0)
            mask_all_net *= mask
        else:
            mask_all_net = np.array(mask)
        mask_part_net = mask ^ mask_all_net
        del mask

        # b. invert once for all pixels with obs in all ifgrams
        if np.sum(mask_all_net) > 0:
            print(('inverting pixels with valid phase in all  ifgrams'
                   ' ({:.0f} pixels) ...').format(np.sum(mask_all_net)))
            tsi, tcohi, num_ifgi = estimate_timeseries(A, B, tbase_diff,
                                                       ifgram=pha_data[:, mask_all_net],
                                                       weight_sqrt=None,
                                                       min_norm_velocity=min_norm_velocity,
                                                       min_redundancy=min_redundancy,
                                                       skip_zero_value=skip_zero_value)
            ts[:, mask_all_net] = tsi
            temp_coh[mask_all_net] = tcohi
            num_inv_ifg[mask_all_net] = num_ifgi

        # c. pixel-by-pixel for pixels with obs not in all ifgrams
        if np.sum(mask_part_net) > 0:
            print(('inverting pixels with valid phase in some ifgrams'
                   ' ({:.0f} pixels) ...').format(np.sum(mask_part_net)))
            num_pixel2inv = int(np.sum(mask_part_net))
            idx_pixel2inv = np.where(mask_part_net)[0]
            prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
            for i in range(num_pixel2inv):
                idx = idx_pixel2inv[i]
                tsi, tcohi, num_ifgi = estimate_timeseries(A, B, tbase_diff,
                                                           ifgram=pha_data[:, idx],
                                                           weight_sqrt=None,
                                                           min_norm_velocity=min_norm_velocity,
                                                           min_redundancy=min_redundancy,
                                                           skip_zero_value=skip_zero_value)
                ts[:, idx] = tsi.flatten()
                temp_coh[idx] = tcohi
                num_inv_ifg[idx] = num_ifgi
                prog_bar.update(i+1, every=2000, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
            prog_bar.close()

    # 2.3 weighted inversion - pixel-by-pixel
    else:
        print('inverting network of interferograms into time-series ...')
        prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
        for i in range(num_pixel2inv):
            idx = idx_pixel2inv[i]
            tsi, tcohi, num_ifgi = estimate_timeseries(A, B, tbase_diff,
                                                       ifgram=pha_data[:, idx],
                                                       weight_sqrt=weight[:, idx],
                                                       min_norm_velocity=min_norm_velocity,
                                                       min_redundancy=min_redundancy,
                                                       skip_zero_value=skip_zero_value)
            ts[:, idx] = tsi.flatten()
            temp_coh[idx] = tcohi
            num_inv_ifg[idx] = num_ifgi

            prog_bar.update(i+1, every=2000, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
        prog_bar.close()
        del weight
    del pha_data

    ## 3. prepare output

    # 3.1 reshape
    ts = ts.reshape(num_date, num_row, num_col)
    #ts_std = ts_std.reshape(num_date, num_row, num_col)
    temp_coh = temp_coh.reshape(num_row, num_col)
    num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)

    # 3.2 convert displacement unit to meter
    if obs_ds_name.startswith('unwrapPhase'):
        phase2range = -1 * float(stack_obj.metadata['WAVELENGTH']) / (4.*np.pi)
        ts *= phase2range
        print('converting LOS phase unit from radian to meter')

    elif obs_ds_name == 'azimuthOffset':
        az_pixel_size = ut.azimuth_ground_resolution(stack_obj.metadata)
        ts *= az_pixel_size
        print('converting azimuth offset unit from pixel ({:.2f} m) to meter'.format(az_pixel_size))

    elif obs_ds_name == 'rangeOffset':
        rg_pixel_size = float(stack_obj.metadata['RANGE_PIXEL_SIZE'])
        ts *= rg_pixel_size
        print('converting range offset unit from pixel ({:.2f} m) to meter'.format(rg_pixel_size))

    return ts, temp_coh, num_inv_ifg, box


def ifgram_inversion(inps=None):
    """Phase triangulatino of small baseline interferograms

    Parameters: inps - namespace
    Example:    inps = cmd_line_parse()
                ifgram_inversion(inps)
    """

    if not inps:
        inps = cmd_line_parse()
    start_time = time.time()

    ## 1. input info

    stack_obj = ifgramStack(inps.ifgramStackFile)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date_list = stack_obj.get_date_list(dropIfgram=True)
    length, width = stack_obj.length, stack_obj.width

    # 1.1 read values on the reference pixel
    inps.refPhase = stack_obj.get_reference_phase(unwDatasetName=inps.obsDatasetName,
                                                  skip_reference=inps.skip_ref,
                                                  dropIfgram=True)

    # 1.2 design matrix
    A = stack_obj.get_design_matrix4timeseries(date12_list)[0]
    num_ifgram, num_date = A.shape[0], A.shape[1]+1
    inps.numIfgram = num_ifgram

    # 1.3 print key setup info
    msg = '-------------------------------------------------------------------------------\n'
    if inps.minNormVelocity:
        suffix = 'deformation velocity'
    else:
        suffix = 'deformation phase'
    msg += 'least-squares solution with L2 min-norm on: {}\n'.format(suffix)
    msg += 'minimum redundancy: {}\n'.format(inps.minRedundancy)
    msg += 'weight function: {}\n'.format(inps.weightFunc)

    if inps.maskDataset:
        if inps.maskDataset in ['coherence', 'offsetSNR']:
            suffix = '{} < {}'.format(inps.maskDataset, inps.maskThreshold)
        else:
            suffix = '{} == 0'.format(inps.maskDataset)
        msg += 'mask out pixels with: {}\n'.format(suffix)
    else:
        msg += 'mask: no\n'

    if np.linalg.matrix_rank(A) < A.shape[1]:
        msg += '***WARNING: the network is NOT fully connected.\n'
        msg += '\tInversion result can be biased!\n'
        msg += '\tContinue to use SVD to resolve the offset between different subsets.\n'
    msg += '-------------------------------------------------------------------------------'
    print(msg)

    print('number of interferograms: {}'.format(num_ifgram))
    print('number of acquisitions  : {}'.format(num_date))
    print('number of lines   : {}'.format(length))
    print('number of columns : {}'.format(width))

    ## 2. prepare output

    # 2.1 metadata
    meta = dict(stack_obj.metadata)
    for key in configKeys:
        meta[key_prefix+key] = str(vars(inps)[key])

    # 2.2 instantiate time-series
    date_dtype = np.dtype('S{}'.format(len(date_list[0])))
    dsNameDict = {
        "date"       : (date_dtype, (num_date,)),
        "bperp"      : (np.float32, (num_date,)),
        "timeseries" : (np.float32, (num_date, length, width)),
    }

    meta['FILE_TYPE'] = 'timeseries'
    meta['UNIT'] = 'm'
    meta['REF_DATE'] = date_list[0]

    ts_obj = timeseries(inps.tsFile)
    ts_obj.layout_hdf5(dsNameDict, meta)

    # write date time-series
    date_list_utf8 = [dt.encode('utf-8') for dt in date_list]
    writefile.write_hdf5_block(inps.tsFile, date_list_utf8, datasetName='date')

    # write bperp time-series
    pbase = stack_obj.get_perp_baseline_timeseries(dropIfgram=True)
    writefile.write_hdf5_block(inps.tsFile, pbase, datasetName='bperp')

    # 2.3 instantiate temporal coherence
    dsNameDict = {"temporalCoherence" : (np.float32, (length, width))}
    meta['FILE_TYPE'] = 'temporalCoherence'
    meta['UNIT'] = '1'
    meta.pop('REF_DATE')
    writefile.layout_hdf5(inps.tempCohFile, dsNameDict, metadata=meta)

    # 2.4 instantiate number of inverted observations
    dsNameDict = {"mask" : (np.float32, (length, width))}
    meta['FILE_TYPE'] = 'mask'
    meta['UNIT'] = '1'
    writefile.layout_hdf5(inps.numInvFile, dsNameDict, metadata=meta)

    ## 3. run the inversion / estimation and write to disk

    # 3.1 split ifgram_file into blocks to save memory
    box_list, num_box = split2boxes(inps.ifgramStackFile, memory_size=inps.memorySize)

    # 3.2 prepare the input arguments for *_patch()
    data_kwargs = {
        "ifgram_file"       : inps.ifgramStackFile,
        "ref_phase"         : inps.refPhase,
        "obs_ds_name"       : inps.obsDatasetName,
        "weight_func"       : inps.weightFunc,
        "min_norm_velocity" : inps.minNormVelocity,
        "water_mask_file"   : inps.waterMaskFile,
        "mask_ds_name"      : inps.maskDataset,
        "mask_threshold"    : inps.maskThreshold,
        "min_redundancy"    : inps.minRedundancy
    }

    # 3.3 invert / write block-by-block
    for i, box in enumerate(box_list):
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_width))
            print('box length: {}'.format(box_length))

        # update box argument in the input data
        data_kwargs['box'] = box

        if inps.cluster == 'no':
            # non-parallel
            ts, temp_coh, num_inv_ifg = ifgram_inversion_patch(**data_kwargs)[:-1]

        else:
            # parallel
            print('\n\n------- start parallel processing using Dask -------')

            # initiate the output data
            ts = np.zeros((num_date, box_length, box_width), np.float32)
            temp_coh     = np.zeros((box_length, box_width), np.float32)
            num_inv_ifg  = np.zeros((box_length, box_width), np.float32)

            # initiate dask cluster and client
            cluster_obj = cluster.DaskCluster(inps.cluster, inps.numWorker, config_name=inps.config)
            cluster_obj.open()

            # run dask
            ts, temp_coh, num_inv_ifg = cluster_obj.run(func=ifgram_inversion_patch,
                                                        func_data=data_kwargs,
                                                        results=[ts, temp_coh, num_inv_ifg])

            # close dask cluster and client
            cluster_obj.close()

            print('------- finished parallel processing -------\n\n')

        # write the block to disk
        # with 3D block in [z0, z1, y0, y1, x0, x1]
        # and  2D block in         [y0, y1, x0, x1]
        # time-series - 3D
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(inps.tsFile,
                                   data=ts,
                                   datasetName='timeseries',
                                   block=block)

        # temporal coherence - 2D
        block = [box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(inps.tempCohFile,
                                   data=temp_coh,
                                   datasetName='temporalCoherence',
                                   block=block)

        # number of inverted obs - 2D
        writefile.write_hdf5_block(inps.numInvFile,
                                   data=num_inv_ifg,
                                   datasetName='mask',
                                   block=block)

        m, s = divmod(time.time() - start_time, 60)
        print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    # 3.4 update output data on the reference pixel
    if not inps.skip_ref:
        # grab ref_y/x
        ref_y = int(stack_obj.metadata['REF_Y'])
        ref_x = int(stack_obj.metadata['REF_X'])
        print('-'*50)
        print('update values on the reference pixel: ({}, {})'.format(ref_y, ref_x))

        print('set temporal coherence on the reference pixel to 1.')
        with h5py.File(inps.tempCohFile, 'r+') as f:
            f['temporalCoherence'][ref_y, ref_x] = 1.

        print('set  # of observations on the reference pixel as {}'.format(num_ifgram))
        with h5py.File(inps.numInvFile, 'r+') as f:
            f['mask'][ref_y, ref_x] = num_ifgram

    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    # Network Inversion
    if inps.residualNorm == 'L2':
        ifgram_inversion(inps)

    else:
        raise NotImplementedError('L1 norm minimization is not fully tested.')
        #ut.timeseries_inversion_L1(inps.ifgramStackFile, inps.tsFile)

    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()
