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
import multiprocessing
import h5py
import numpy as np
from scipy import linalg   # more effieint than numpy.linalg
from mintpy.objects import ifgramStack, timeseries
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
  ifgram_inversion.py  inputs/ifgramStack.h5 -w var

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
                        default=['timeseries.h5', 'temporalCoherence.h5'],
                        help='Output file name. (default: %(default)s).')
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',
                        help='Skip checking reference pixel value, for simulation testing.')

    # solver
    solver = parser.add_argument_group('solver', 'solver for the network inversion problem')
    solver.add_argument('-w', '--weight-func', dest='weightFunc', default='no',
                        choices={'var', 'fim', 'coh', 'no'},
                        help='function used to convert coherence to weight for inversion:\n' +
                             'var - inverse of phase variance due to temporal decorrelation\n' +
                             'fim - Fisher Information Matrix as weight' +
                             'coh - spatial coherence\n' +
                             'no  - no/uniform weight (default)')
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
    parser.add_argument('-r', '--ram', '--memory', dest='memorySize', type=float, default=4,
                        help='Max amount of memory in GB to use (default: %(default)s).\n' +
                             'Adjust according to your computer memory.')

    par = parser.add_argument_group('parallel', 'parallel processing using dask')
    par.add_argument('--cluster', '--cluster-type', dest='cluster', type=str,
                     default='local', choices={'local', 'lsf', 'pbs', 'slurm', 'no'},
                     help='Cluster to use for parallel computing, no to turn OFF. (default: %(default)s).')
    par.add_argument('--num-worker', dest='numWorker', type=str, default='4',
                     help='Number of workers to use (default: %(default)s).')

    par.add_argument('--config', '--config-name', dest='config', type=str, default='no', 
                     help='Configuration name to use in dask.yaml (default: %(default)s).')
    par.add_argument('--walltime', dest='walltime', type=str, default='00:40',
                     help='Walltime for each dask worker (default: %(default)s).')

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

    # --num-worker option
    if inps.cluster == 'local':
        # translate numWorker = all
        num_core = multiprocessing.cpu_count()
        if inps.numWorker == 'all':
            inps.numWorker = num_core

        inps.numWorker = int(inps.numWorker)
        if inps.numWorker > num_core:
            msg = '\nWARNING: input number of worker: {} > available cores: {}'.format(inps.numWorker, num_core)
            msg += '\nchange number of worker to {} and continue\n'.format(int(num_core/2))
            print(msg)
            inps.numWorker = int(num_core / 2)

    else:
        if inps.numWorker == 'all':
            msg = 'numWorker = all is NOT supported for cluster type: {}'.format(inps.cluster)
            raise ValueError(msg)
        inps.numWorker = int(inps.numWorker)

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
            inps.outfile = ['timeseries.h5', 'temporalCoherence.h5']
        elif inps.obsDatasetName.startswith('azimuthOffset'):
            inps.outfile = ['timeseriesAz.h5', 'temporalCoherenceAz.h5']
        elif inps.obsDatasetName.startswith('rangeOffset'):
            inps.outfile = ['timeseriesRg.h5', 'temporalCoherenceRg.h5']
        else:
            raise ValueError('un-recognized input observation dataset name: {}'.format(inps.obsDatasetName))

    inps.timeseriesFile, inps.tempCohFile = inps.outfile

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
        if key in ['maskDataset', 'minNormVelocity']:
            iDict[key] = value
        elif value:
            if key in ['maskThreshold', 'minRedundancy', 'memorySize']:
                iDict[key] = float(value)
            elif key in ['weightFunc', 'residualNorm', 'waterMaskFile']:
                iDict[key] = value

    # computing configurations
    dask_key_prefix = 'mintpy.compute.'
    keyList = [i for i in list(iDict.keys()) if dask_key_prefix+i in template.keys()]
    for key in keyList:
        value = template[dask_key_prefix+key]
        if key in ['cluster', 'config']:
            iDict[key] = value
        elif value:
            if key in ['walltime', 'numWorker']:
                iDict[key] = str(value)
            elif key in ['memorySize']:
                iDict[key] = float(value)

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
        atr_ts = readfile.read_attribute(inps.timeseriesFile)
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
                                    baseline history between master and slave date of one interferogram
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
                num_inv_ifg       - 1D np.array in size of (num_pixel), number of ifgrams
                                    used during the inversion
    """
    ifgram = ifgram.reshape(A.shape[0], -1)
    if weight_sqrt is not None:
        weight_sqrt = weight_sqrt.reshape(A.shape[0], -1)
    num_date = A.shape[1] + 1
    num_pixel = ifgram.shape[1]

    # Initial output value
    ts = np.zeros((num_date, num_pixel), np.float32)
    temp_coh = 0.
    num_inv_ifg = 0

    # Skip Zero Phase Value
    if skip_zero_value and not np.all(ifgram):
        idx = (ifgram[:, 0] != 0.).flatten()
        A = A[idx, :]
        B = B[idx, :]

        # Skip the pixel if its redundancy < threshold
        if np.min(np.sum(A != 0., axis=0)) < min_redundancy:
            return ts, temp_coh, num_inv_ifg

        # check matrix invertability
        if weight_sqrt is not None:  #for WLS only because OLS contains it already
            try:
                linalg.inv(np.dot(B.T, B))
            except linalg.LinAlgError:
                return ts, temp_coh, num_inv_ifg

        ifgram = ifgram[idx, :]
        if weight_sqrt is not None:
            weight_sqrt = weight_sqrt[idx, :]

    # invert time-series
    try:
        # assume minimum-norm deformation velocity
        if min_norm_velocity:
            if weight_sqrt is not None:
                B_w = np.multiply(B, weight_sqrt)
                ifgram_w = np.multiply(ifgram, weight_sqrt)
                X = linalg.lstsq(B_w, ifgram_w, cond=rcond)[0]
            else:
                X = linalg.lstsq(B, ifgram, cond=rcond)[0]

            ts_diff = X * np.tile(tbase_diff, (1, num_pixel))
            ts[1:, :] = np.cumsum(ts_diff, axis=0)
            ifgram_diff = ifgram - np.dot(B, X)

        # assume minimum-norm deformation phase
        else:
            if weight_sqrt is not None:
                A_w = np.multiply(A, weight_sqrt)
                ifgram_w = np.multiply(ifgram, weight_sqrt)
                X = linalg.lstsq(A_w, ifgram_w, cond=rcond)[0]
            else:
                X = linalg.lstsq(A, ifgram, cond=rcond)[0]
            ts[1: ,:] = X
            ifgram_diff = ifgram - np.dot(A, X)

        # calculate temporal coherence
        num_inv_ifg = A.shape[0]
        temp_coh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / num_inv_ifg

    except linalg.LinAlgError:
        pass

    return ts, temp_coh, num_inv_ifg


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


def write_aux2hdf5_file(metadata, temp_coh, num_inv_ifg=None, suffix='', inps=None):

    # File 2 - temporalCoherence.h5
    out_file = '{}{}.h5'.format(suffix, os.path.splitext(inps.outfile[1])[0])
    metadata['FILE_TYPE'] = 'temporalCoherence'
    metadata['UNIT'] = '1'
    print('-'*50)
    writefile.write(temp_coh, out_file=out_file, metadata=metadata)

    # File 3 - numInvIfgram.h5
    out_file = 'numInvIfgram{}.h5'.format(suffix)
    metadata['FILE_TYPE'] = 'mask'
    metadata['UNIT'] = '1'
    print('-'*50)
    writefile.write(num_inv_ifg, out_file=out_file, metadata=metadata)

    return None


def split2boxes(dataset_shape, memory_size=4, print_msg=True):
    """Split into chunks in rows to reduce memory usage
    Parameters: dataset_shape - tuple of 3 int
                memory_size   - float, max memory to use in GB
                print_msg     - bool
    Returns:    box_list      - list of tuple of 4 int
    """
    # memory_size --> chunk_size
    # 10 is from phase (4 bytes), weight (4 bytes)
    # and time-series (4 bytes but half the size on the 1st dimension on average)
    chunk_size = memory_size * (1024**3) / 10

    # chunk_size --> y_step / chunk_num
    length, width = dataset_shape[1:3]
    y_step = chunk_size / (dataset_shape[0] * width)         # split in lines
    y_step = int(ut.round_to_1(y_step))
    chunk_num = int((length - 1) / y_step) + 1

    if print_msg and chunk_num > 1:
        print('maximum memory size: %.1E GB' % memory_size)
        print('maximum chunk  size: %.1E' % chunk_size)
        print('split %d lines into %d patches for processing' % (length, chunk_num))
        print('    with each patch up to %d lines' % y_step)

    # y_step / chunk_num --> box_list
    box_list = []
    for i in range(chunk_num):
        y0 = i * y_step
        y1 = min([length, y0 + y_step])
        box = (0, y0, width, y1)
        box_list.append(box)

    return box_list


def split_box2sub_boxes(box, num_split, dimension='x'):
    """Further divides the box size into `num_split` different sub_boxes.
    Note that this is different from `split2boxes()`, whic splits based on chunk_size (memory-based).

    :param box: [x0, y0, x1, y1]: list[int] of size 4
    :param num_split: int, the number of sub_boxes to split a box into
    :param dimension: str = 'y' or 'x', the dimension along which to split the boxes
    """
    x0, y0, x1, y1 = box
    length, width = y1 - y0, x1 - x0

    sub_boxes = []
    if dimension == 'y':
        for i in range(num_split):
            start = (i * length) // num_split + y0
            end = ((i + 1) * length) // num_split + y0
            sub_boxes.append([x0, start, x1, end])

    else:
        for i in range(num_split):
            start = (i * width) // num_split + x0
            end = ((i + 1) * width) // num_split + x0
            sub_boxes.append([start, y0, end, y1])

    return sub_boxes


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

    ## debug
    #y, x = 258, 454
    #box = (x, y, x+1, y+1)

    # Size Info - Patch
    if box:
        num_row = box[3] - box[1]
        num_col = box[2] - box[0]
    else:
        num_row = stack_obj.length
        num_col = stack_obj.width
    num_pixel = num_row * num_col

    # get tbase_diff
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)
    tbase = np.array(ptime.date_list2tbase(date_list)[0], np.float32) / 365.25
    tbase_diff = np.diff(tbase).reshape(-1, 1)

    # Design matrix
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

    # Initialization of output matrix
    ts     = np.zeros((num_date, num_pixel), np.float32)
    #ts_std = np.zeros((num_date, num_pixel), np.float32)
    temp_coh    = np.zeros(num_pixel, np.float32)
    num_inv_ifg = np.zeros(num_pixel, np.int16)

    # Read/Mask unwrapPhase / offset
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

    # Mask for pixels to invert
    mask = np.ones(num_pixel, np.bool_)

    # 1 - Water Mask
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

    # 2 - Mask for Zero Phase in ALL ifgrams
    if 'phase' in obs_ds_name.lower():
        print('skip pixels with zero/nan value in all interferograms')
        with warnings.catch_warnings():
            # ignore warning message for all-NaN slices
            warnings.simplefilter("ignore", category=RuntimeWarning)
            phase_stack = np.nanmean(pha_data, axis=0)
        mask *= np.multiply(~np.isnan(phase_stack), phase_stack != 0.)
        del phase_stack

    # Invert pixels on mask 1+2
    num_pixel2inv = int(np.sum(mask))
    idx_pixel2inv = np.where(mask)[0]
    print('number of pixels to invert: {} out of {} ({:.1f}%)'.format(num_pixel2inv,
                                                                      num_pixel,
                                                                      num_pixel2inv/num_pixel*100))

    if num_pixel2inv < 1:
        ts = ts.reshape(num_date, num_row, num_col)
        #ts_std = ts_std.reshape(num_date, num_row, num_col)
        temp_coh = temp_coh.reshape(num_row, num_col)
        num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)
        return ts, temp_coh, num_inv_ifg

    # skip zero value in the network inversion for phase
    if 'phase' in obs_ds_name.lower():
        skip_zero_value = True
    else:
        skip_zero_value = False

    # Inversion - SBAS
    if weight_func in ['no', 'sbas']:
        # Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        if 'phase' in obs_ds_name.lower():
            mask_all_net = np.all(pha_data, axis=0)
            mask_all_net *= mask
        else:
            mask_all_net = np.array(mask)
        mask_part_net = mask ^ mask_all_net

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

    # Inversion - WLS
    else:
        # calculate weight from coherence
        L = int(stack_obj.metadata['ALOOKS']) * int(stack_obj.metadata['RLOOKS'])
        weight = read_coherence(stack_obj, box=box, dropIfgram=True)
        weight = decor.coherence2weight(weight, weight_func, L=L, epsilon=5e-2)
        weight = np.sqrt(weight)

        # weighted inversion pixel by pixel
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

    ts = ts.reshape(num_date, num_row, num_col)
    #ts_std = ts_std.reshape(num_date, num_row, num_col)
    temp_coh = temp_coh.reshape(num_row, num_col)
    num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)

    # convert displacement unit to meter
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
    start_time = time.time()

    # Check Inputs
    if not inps:
        inps = cmd_line_parse()

    # basic info
    stack_obj = ifgramStack(inps.ifgramStackFile)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date_list = stack_obj.get_date_list(dropIfgram=True)
    length, width = stack_obj.length, stack_obj.width

    # design matrix
    A = stack_obj.get_design_matrix4timeseries(date12_list)[0]
    num_ifgram, num_date = A.shape[0], A.shape[1]+1
    inps.numIfgram = num_ifgram

    # print key setup info
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

    # split ifgram_file into blocks to save memory
    box_list = split2boxes(dataset_shape=stack_obj.get_size(), memory_size=inps.memorySize)
    num_box = len(box_list)

    # read ifgram_file in small patches and write them together
    inps.refPhase = stack_obj.get_reference_phase(unwDatasetName=inps.obsDatasetName,
                                                  skip_reference=inps.skip_ref,
                                                  dropIfgram=True)

    # Initialization of output matrix
    temp_coh    = np.zeros((length, width), np.float32)
    num_inv_ifg = np.zeros((length, width), np.int16)

    # metadata
    metadata = dict(stack_obj.metadata)
    for key in configKeys:
        metadata[key_prefix+key] = str(vars(inps)[key])

    metadata['REF_DATE'] = date_list[0]
    metadata['FILE_TYPE'] = 'timeseries'
    metadata['UNIT'] = 'm'

    # instantiate a timeseries object
    ts_file = '{}.h5'.format(os.path.splitext(inps.outfile[0])[0])
    ts_obj = timeseries(ts_file)

    # layout the HDF5 file for the datasets and the metadata
    dsNameDict = {"date": ((np.dtype('S8'), (num_date,))),
                  "bperp": (np.float32, (num_date,)),
                  "timeseries": (np.float32, (num_date, length, width))}
    ts_obj.layout_hdf5(dsNameDict, metadata)

    # prepare the input arguments for *_patch()
    kwargs = {
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

    # invert & write block by block
    for box_i, box in enumerate(box_list):
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(box_i+1, num_box))

        # initiate the output matrice for the box
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        tsi = np.zeros((num_date, box_length, box_width), np.float32)
        temp_cohi    = np.zeros((box_length, box_width), np.float32)
        num_inv_ifgi = np.zeros((box_length, box_width), np.float32)
        print('box width:  {}'.format(box_width))
        print('box length: {}'.format(box_length))

        if inps.cluster.lower() == 'no':
            kwargs['box'] = box
            tsi, temp_cohi, num_inv_ifgi, box = ifgram_inversion_patch(**kwargs)

        else:
            print('\n\n'+'------- start parallel processing using dask -------'+'\n\n')

            try:
                from dask.distributed import Client, as_completed
            except ImportError:
                raise ImportError('Cannot import dask.distributed!')
            from mintpy.objects.cluster import DaskCluster

            # intiate the cluster client
            # Look at the ~/.config/dask/mintpy.yaml file for changing the Dask configuration defaults
            print('initiate dask cluster')
            cluster_obj = DaskCluster(cluster_type=inps.cluster, walltime=inps.walltime, config_name=inps.config)
            cluster = cluster_obj.cluster
            # This line submits NUM_WORKERS jobs to the cluster to start a bunch of workers
            # In tests on Pegasus `general` queue in Jan 2019, no more than 40 workers could RUN
            # at once (other user's jobs gained higher priority in the general at that point)
            print('scale the cluster to {} workers'.format(inps.numWorker))
            NUM_WORKERS = inps.numWorker
            cluster.scale(NUM_WORKERS)

            # This line needs to be in a function or in a `if __name__ == "__main__":` block. If it is in no function
            # or "main" block, each worker will try to create its own client (which is bad) when loading the module
            print('initiate dask client')
            client = Client(cluster)

            # split the primary box into sub boxes for each worker
            sub_boxes = split_box2sub_boxes(box, num_split=1*NUM_WORKERS, dimension='x')
            print('split the patch to {} sub boxes in x direction for workers to process'.format(len(sub_boxes)))

            # submit jobs for each worker
            start_time_sub = time.time()
            futures = []
            for i, sub_box in enumerate(sub_boxes):
                print('submit job to workers for sub box {}: {}'.format(i, sub_box))
                kwargs['box'] = sub_box

                # David: I haven't played with fussing with `retries`, however sometimes a future fails
                # on a worker for an unknown reason. retrying will save the whole process from failing.
                # TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
                # thrown in that case, therefore I don't know how to recognize when this happens.
                future = client.submit(ifgram_inversion_patch, **kwargs, retries=3)
                futures.append(future)

            # assemble results from all workers
            i_future = 0
            for future, result in as_completed(futures, with_results=True):
                # catch result
                sub_tsi, sub_temp_cohi, sub_num_inv_ifgi, sub_box = result

                # message
                i_future += 1
                sub_t = time.time() - start_time_sub
                print("FUTURE #{} box {} complete. Time used: {:.0f} seconds".format(i_future, sub_box, sub_t))

                # convert the abosulte sub_box into local col/row start/end relative to the primary box
                # to assemble the result from each worker
                x0, y0, x1, y1 = sub_box
                x0 -= box[0]
                x1 -= box[0]
                y0 -= box[1]
                y1 -= box[1]

                tsi[:, y0:y1, x0:x1] = sub_tsi
                temp_cohi[y0:y1, x0:x1] = sub_temp_cohi
                num_inv_ifgi[y0:y1, x0:x1] = sub_num_inv_ifgi

            # close dask cluster and client
            cluster.close()
            client.close()
            print('close dask cluster')
            print('close dask client')

            # move *.o/.e files produced by dask in stdout/stderr
            ut.move_dask_stdout_stderr_files()

            print('\n\n------- finished parallel processing -------\n\n')

        # write the block of timeseries to disk
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        ts_obj.write2hdf5_block(tsi, datasetName='timeseries', block=block)

        # save the block of aux datasets
        temp_coh[box[1]:box[3], box[0]:box[2]] = temp_cohi
        num_inv_ifg[box[1]:box[3], box[0]:box[2]] = num_inv_ifgi

    # write date and bperp to disk
    print('-'*50)
    date_list_utf8 = [dt.encode('utf-8') for dt in date_list]
    ts_obj.write2hdf5_block(date_list_utf8, datasetName='date')

    pbase = stack_obj.get_perp_baseline_timeseries(dropIfgram=True)
    ts_obj.write2hdf5_block(pbase, datasetName='bperp')

    # reference pixel
    if not inps.skip_ref:
        ref_y = int(stack_obj.metadata['REF_Y'])
        ref_x = int(stack_obj.metadata['REF_X'])
        num_inv_ifg[ref_y, ref_x] = num_ifgram
        temp_coh[ref_y, ref_x] = 1.

    # write auxliary data to files: temporal coherence, number of inv ifgrams, etc.
    write_aux2hdf5_file(metadata, temp_coh, num_inv_ifg, suffix='', inps=inps)

    m, s = divmod(time.time()-start_time, 60)
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
        #ut.timeseries_inversion_L1(inps.ifgramStackFile, inps.timeseriesFile)

    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()
