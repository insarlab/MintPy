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
  ifgram_inversion.py  inputs/ifgramStack.h5 -t smallbaselineApp.cfg --fast
  ifgram_inversion.py  inputs/ifgramStack.h5 -w var
  ifgram_inversion.py  inputs/ifgramStack.h5 -w fim
  ifgram_inversion.py  inputs/ifgramStack.h5 -w coh

  # parallel processing for HPC
  ifgram_inversion.py  inputs/ifgramStack.h5 -w var --parallel
  ifgram_inversion.py  inputs/ifgramStack.h5 -w var --parallel --num-worker 25

  # invert offset stack
  ifgram_inversion.py  inputs/ifgramStack.h5 -i azimuthOffset --water-mask waterMask.h5 --mask-dset offsetSNR --mask-threshold 5
"""

TEMPLATE = get_template_content('invert_network')

REFERENCE = """references:
  Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002). A new algorithm for surface
    deformation monitoring based on small baseline differential SAR interferograms. IEEE TGRS,
    40(11), 2375-2383. doi:10.1109/TGRS.2002.803792
  Guarnieri, A. M., and S. Tebaldini (2008), On the exploitation of target statistics for SAR
    interferometry applications, Geoscience and Remote Sensing, IEEE Transactions on, 46(11), 3436-3443.
  Just, D., & Bamler, R. (1994). Phase statistics of interferograms with applications to synthetic
    aperture radar. Applied optics, 33(20), 4361-4368.
  Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    of multitemporal differential SAR interferograms, IEEE-TGRS, 44(9), 2374-2383.
  Perissin, D., and T. Wang (2012), Repeat-pass SAR interferometry with partially coherent targets, IEEE TGRS,
    50(1), 271-280, doi:10.1109/tgrs.2011.2160644.
  Samiei-Esfahany, S., J. E. Martins, F. v. Leijen, and R. F. Hanssen (2016), Phase Estimation for Distributed
    Scatterers in InSAR Stacks Using Integer Least Squares Estimation, IEEE TGRS, 54(10), 5671-5687.
  Seymour, M. S., and I. G. Cumming (1994), Maximum likelihood estimation for SAR interferometry, 1994.
    IGARSS '94., 8-12 Aug 1994.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Invert network of interferograms into time-series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+TEMPLATE+'\n'+EXAMPLE)
    # input dataset
    parser.add_argument('ifgramStackFile', help='interferograms stack file to be inverted')
    parser.add_argument('-i','-d', '--dset', dest='obsDatasetName', type=str,
                        help='dataset name of unwrap phase / offset to be used for inversion\n'
                             'e.g.: unwrapPhase, unwrapPhase_bridging, ...')
    parser.add_argument('--water-mask', '-m', dest='waterMaskFile',
                        help='Skip inversion on the masked out region, i.e. water.')
    parser.add_argument('--template', '-t', dest='templateFile',
                        help='template text file with options')
    parser.add_argument('-o', '--output', dest='outfile', nargs=2,
                        metavar=('TS_FILE', 'TCOH_FILE'), default=['timeseries.h5', 'temporalCoherence.h5'],
                        help='Output file name for timeseries and temporal coherence, default:\n' +
                        'timeseries.h5 temporalCoherence.h5')

    # options rarely used or changed
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--chunk-size', dest='chunk_size', type=float, default=100e6,
                        help='max number of data (= ifgram_num * num_row * num_col) to read per loop\n' +
                        'default: 0.2 G; adjust it according to your computer memory.')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',
                        help='Skip checking reference pixel value, for simulation testing.')

    # solver
    solver = parser.add_argument_group('solver', 'solver for the network inversion problem')
    solver.add_argument('--weight-function', '-w', dest='weightFunc', default='no', choices={'var', 'fim', 'coh', 'no'},
                        help='function used to convert coherence to weight for inversion:\n' +
                        'var - inverse of phase variance due to temporal decorrelation\n' +
                        'fim - Fisher Information Matrix as weight' +
                        'coh - spatial coherence\n' +
                        'no  - no/uniform weight')
    solver.add_argument('--min-norm-phase', dest='minNormVelocity', action='store_false',
                        help=('Enable inversion with minimum-norm deformation phase,'
                              ' instead of the default minimum-norm deformation velocity.'))
    solver.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1', 'L2'],
                        help='Inverse method used to residual optimization, L1 or L2 norm minimization. Default: L2')

    # mask
    mask = parser.add_argument_group('mask', 'mask observation data before inversion')
    mask.add_argument('--mask-dset', dest='maskDataset',
                      help='dataset used to mask unwrapPhase, e.g. coherence, connectComponent')
    mask.add_argument('--mask-threshold', dest='maskThreshold', metavar='NUM', type=float, default=0.4,
                      help='threshold to generate mask when mask is coherence')
    mask.add_argument('--min-redundancy', dest='minRedundancy', metavar='NUM', type=float, default=1.0,
                      help='minimum redundancy of interferograms for every SAR acquisition.')

    # parallel computing
    par = parser.add_argument_group('parallel', 'parallel processing configuration for Dask')
    par.add_argument('--parallel', dest='parallel', action='store_true',
                     help='Enable parallel processing for the pixelwise weighted inversion.')
    par.add_argument('--cluster', '--cluster-type', dest='cluster', type=str,
                     default='SLURM', choices={'LSF', 'PBS', 'SLURM'},
                     help='Type of HPC cluster you are running on (default: %(default)s).')
    par.add_argument('--config', '--config-name', dest='config', type=str, default='no', 
                     help='Configuration name to use in dask.yaml (default: %(default)s).')
    par.add_argument('--num-worker', dest='numWorker', type=int, default=40,
                     help='Number of workers the Dask cluster should use (default: %(default)s).')
    par.add_argument('--walltime', dest='walltime', type=str, default='00:40',
                     help='Walltime for each dask worker (default: %(default)s).')

    # efficiency
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n' +
                        'readable and newer than input interferograms file')
    parser.add_argument('--fast','--sbas', action='store_true',
                        help='Fast network invertion by forcing the following options:\n'+
                             '\t--weight-function = no\n'+
                             '\t--mask-dset = no\n'+
                             'This is equivalent to SBAS algorithm (Berardino et al., 2002)')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.ifgramStackFile)
    if atr['FILE_TYPE'] not in ['ifgramStack']:
        raise ValueError('input is {} file, support ifgramStack file only.'.format(atr['FILE_TYPE']))

    if inps.templateFile:
        inps = read_template2inps(inps.templateFile, inps)

    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    # --fast option
    if inps.fast:
        print("Enable fast network inversion.")
        if inps.weightFunc != 'no':
            inps.weightFunc = 'no'
            print("\tforcing weightFunc = 'no'")
        if inps.maskDataset is not None:
            inps.maskDataset = None
            print("\tforcing maskDataset = None")

    # --dset option
    if not inps.obsDatasetName:
        stack_obj = ifgramStack(inps.ifgramStackFile)
        stack_obj.open(print_msg=False)
        inps.obsDatasetName = [i for i in ['unwrapPhase_bridging_phaseClosure',
                                           'unwrapPhase_bridging',
                                           'unwrapPhase_phaseClosure',
                                           'unwrapPhase']
                               if i in stack_obj.datasetNames][0]

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
    
    # --config option
    if inps.config == 'no':
        inps.config = inps.cluster.lower()

    if inps.parallel and inps.config != inps.cluster.lower():
        import dask
        try:
            dask.config.get('jobqueue.{}'.format(inps.config))
        except KeyError:
            msg = 'Dask configuration "{}" was not found in ~/.config/dask/mintpy.yaml'.format(inps.config)
            msg += '\nFall back to default config name: "{}"'.format(inps.cluster.lower())
            print(msg)
            inps.config = inps.cluster.lower()

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
        if key in ['maskDataset', 'minNormVelocity', 'parallel', 'cluster']:
            iDict[key] = value
        elif value:
            if key in ['numWorker']:
                iDict[key] = int(value)
            elif key in ['walltime', 'config']:
                iDict[key] = str(value)
            elif key in ['maskThreshold', 'minRedundancy']:
                iDict[key] = float(value)
            elif key in ['weightFunc', 'residualNorm', 'waterMaskFile']:
                iDict[key] = value
    return inps


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


#################################### Weight Functions #####################################
def coherence2phase_variance(coherence, L=32, epsilon=1e-3, print_msg=False):
    """Convert coherence to phase variance based on DS phase PDF (Tough et al., 1995)"""
    lineStr = '    number of looks L={}'.format(L)
    if L > 80:
        L = 80
        lineStr += ', use L=80 to avoid dividing by 0 in calculation with negligible effect'
    if print_msg:
        print(lineStr)

    coh_num = 1000
    coh_min = 0.0 + epsilon
    coh_max = 1.0 - epsilon
    coh_lut = np.linspace(coh_min, coh_max, coh_num)
    coh_min = np.min(coh_lut)
    coh_max = np.max(coh_lut)
    coh_step = (coh_max - coh_min) / (coh_num - 1)

    coherence = np.array(coherence)
    coherence[coherence < coh_min] = coh_min
    coherence[coherence > coh_max] = coh_max
    coherence_idx = np.array((coherence - coh_min) / coh_step, np.int16)

    var_lut = decor.phase_variance_ds(int(L), coh_lut)[0]
    variance = var_lut[coherence_idx]
    return variance


def coherence2fisher_info_index(data, L=32, epsilon=1e-3):
    """Convert coherence to Fisher information index (Seymour & Cumming, 1994, IGARSS)"""
    if data.dtype != np.float64:
        data = np.array(data, np.float64)
    data[data > 1-epsilon] = 1-epsilon
    data = 2.0 * int(L) * np.square(data) / (1 - np.square(data))
    return data


def coherence2weight(coh_data, weight_func='var', L=20, epsilon=5e-2, print_msg=True):
    coh_data[np.isnan(coh_data)] = epsilon
    coh_data[coh_data < epsilon] = epsilon
    coh_data = np.array(coh_data, np.float64)

    # Calculate Weight matrix
    weight_func = weight_func.lower()
    if 'var' in weight_func:
        if print_msg:
            print('convert coherence to weight using inverse of phase variance')
            print('    with phase PDF for distributed scatterers from Tough et al. (1995)')
        weight = 1.0 / coherence2phase_variance(coh_data, L, print_msg=print_msg)

    elif any(i in weight_func for i in ['coh', 'lin']):
        if print_msg:
            print('use coherence as weight directly (Perissin & Wang, 2012; Tong et al., 2016)')
        weight = coh_data

    elif any(i in weight_func for i in ['fim', 'fisher']):
        if print_msg:
            print('convert coherence to weight using Fisher Information Index (Seymour & Cumming, 1994)')
        weight = coherence2fisher_info_index(coh_data, L)

    elif weight_func in ['no', 'sbas', 'uniform']:
        weight = None

    else:
        raise Exception('Un-recognized weight function: %s' % weight_func)

    if weight is not None:
        weight = np.array(weight, np.float32)
    del coh_data
    return weight


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

    Parameters: A - 2D np.array in size of (num_ifgram, num_date-1)
                B - 2D np.array in size of (num_ifgram, num_date-1),
                    design matrix B, each row represents differential temporal
                    baseline history between master and slave date of one interferogram
                tbase_diff - 2D np.array in size of (num_date-1, 1),
                    differential temporal baseline history
                ifgram - 2D np.array in size of (num_ifgram, num_pixel),
                    phase of all interferograms
                weight_sqrt - 2D np.array in size of (num_ifgram, num_pixel),
                    square root of weight of all interferograms
                min_norm_velocity - bool, assume minimum-norm deformation velocity, or not
                rcond - cut-off ratio of small singular values of A or B, to maintain robustness.
                    It's recommend to >= 1e-5 by experience, to generate reasonable result.
                min_redundancy - min redundancy defined as min num_ifgram for every SAR acquisition
    Returns:    ts - 2D np.array in size of (num_date, num_pixel), phase time-series
                temp_coh - 1D np.array in size of (num_pixel), temporal coherence
                num_inv_ifg - 1D np.array in size of (num_pixel), number of ifgrams
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


def write2hdf5_auxFiles(metadata, temp_coh, num_inv_ifg=None, suffix='', inps=None):

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


def split_ifgram_file(ifgram_file, chunk_size=100e6):
    """Split ifgramStack file into several smaller files."""
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    metadata = dict(stack_obj.metadata)

    # get reference phase
    ref_phase = stack_obj.get_reference_phase(dropIfgram=False)

    # get list of boxes
    box_list = split2boxes(dataset_shape=stack_obj.get_size(),
                           chunk_size=chunk_size,
                           print_msg=True)
    num_box = len(box_list)

    # read/write each patch file
    outfile_list = []
    for i in range(num_box):
        box = box_list[i]
        outfile = '{}_{:03d}{}'.format(os.path.splitext(ifgram_file)[0], i+1,
                                       os.path.splitext(ifgram_file)[1])

        # datasets
        print('-'*50)
        print('reading all datasets in {} from file: {} ...'.format(box, ifgram_file))
        dsNames = readfile.get_dataset_list(ifgram_file)
        dsDict = {}
        dsDict['refPhase'] = ref_phase
        for dsName in dsNames:
            data = stack_obj.read(datasetName=dsName, box=box, print_msg=False)
            dsDict[dsName] = data

        # metadata
        metadata['LENGTH'] = box[3] - box[1]
        metadata['WIDTH'] = box[2] - box[0]
        writefile.write(dsDict, out_file=outfile, metadata=metadata, ref_file=ifgram_file)
        outfile_list.append(outfile)
    return outfile_list


def split2boxes(dataset_shape, chunk_size=100e6, print_msg=True):
    """Split into chunks in rows to reduce memory usage
    Parameters:
    """
    # Get r_step / chunk_num
    r_step = chunk_size / (dataset_shape[0] * dataset_shape[2])         # split in lines
    r_step = int(ut.round_to_1(r_step))
    chunk_num = int((dataset_shape[1]-1)/r_step) + 1

    if print_msg and chunk_num > 1:
        print('maximum chunk size: %.1E' % chunk_size)
        print('split %d lines into %d patches for processing' % (dataset_shape[1], chunk_num))
        print('    with each patch up to %d lines' % r_step)

    # Computing the inversion
    box_list = []
    for i in range(chunk_num):
        r0 = i * r_step
        r1 = min([dataset_shape[1], r0+r_step])
        box = (0, r0, dataset_shape[2], r1)
        box_list.append(box)
    return box_list


def subsplit_boxes4_workers(box, num_split, dimension='y'):
    """ This is a bit hacky, but after creating the patches,
    this function further divides the box size into `num_split` different subboxes.
    Note that `split2boxes`  splits based on chunk_size (memory-based).

    :param box: [x0, y0, x1, y1]: list[int] of size 4
    :param num_split: int, the number of subboxes to split a box into
    :param dimension: str = 'y' or 'x', the dimension along which to split the boxes
    """

    # Flip x and y coordinates if splitting along 'x' dimension
    x0, y0, x1, y1 = box
    subboxes = []

    if dimension == 'y':
        y_diff = y1 - y0
        # `start` and `end` are the new bounds of the subdivided box
        for i in range(num_split):
            start = (i * y_diff) // num_split
            end = ((i + 1) * y_diff) // num_split
            subboxes.append([x0, start, x1, end])
    elif dimension == 'x':
        x_diff = x1 - x0
        for i in range(num_split):
            start = (i * x_diff) // num_split
            end = ((i + 1) * x_diff) // num_split
            subboxes.append([start, y0, end, y1])
    else:
        raise Exception("Unknown value for dimension parameter:", dimension)

    return subboxes


def check_design_matrix(ifgram_file, weight_func='var'):
    """Check Rank of Design matrix for weighted inversion"""
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


def read_unwrap_phase(stack_obj, box, ref_phase, obsDatasetName='unwrapPhase', dropIfgram=True,
                      print_msg=True):
    """Read unwrapPhase from ifgramStack file
    Parameters: stack_obj : ifgramStack object
                box : tuple of 4 int
                ref_phase : 1D array or None
    Returns:    pha_data : 2D array of unwrapPhase in size of (num_ifgram, num_pixel)
    """
    # Read unwrapPhase
    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if print_msg:
        print('reading {} in {} * {} ...'.format(obsDatasetName, box, num_ifgram))
    pha_data = stack_obj.read(datasetName=obsDatasetName,
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
        # read ref_phase from file itself
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


def mask_unwrap_phase(pha_data, stack_obj, box, mask_ds_name=None, mask_threshold=0.4, dropIfgram=True,
                      print_msg=True):
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
    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if print_msg:
        print('reading coherence in {} * {} ...'.format(box, num_ifgram))
    coh_data = stack_obj.read(datasetName='coherence',
                              box=box,
                              dropIfgram=dropIfgram,
                              print_msg=False).reshape(num_ifgram, -1)
    coh_data[np.isnan(coh_data)] = 0.
    return coh_data


def ifgram_inversion_patch(ifgram_file, box=None, ref_phase=None, obsDatasetName='unwrapPhase',
                           weight_func='var', min_norm_velocity=True,
                           mask_dataset_name=None, mask_threshold=0.4, min_redundancy=1.0,
                           water_mask_file=None):
    """Invert one patch of an ifgram stack into timeseries.
    Parameters: ifgram_file       : str, interferograms stack HDF5 file, e.g. ./inputs/ifgramStack.h5
                box               : tuple of 4 int, indicating (x0, y0, x1, y1) pixel coordinate of area of interest
                                    or None, to process the whole file and write output file
                ref_phase         : 1D array in size of (num_ifgram)
                                    or None
                weight_func       : str, weight function, choose in ['no', 'fim', 'var', 'coh']
                mask_dataset_name : str, dataset name in ifgram_file used to mask unwrapPhase pixelwisely
                mask_threshold    : float, min coherence of pixels if mask_dataset_name='coherence'
                water_mask_file   : str, water mask filename if available,
                                    skip inversion on water to speed up the process
    Returns:    ts          : 3D array in size of (num_date, num_row, num_col)
                temp_coh    : 2D array in size of (num_row, num_col)
                num_inv_ifg : 2D array in size of (num_row, num_col)
    Example:    ifgram_inversion_patch('ifgramStack.h5', box=(0,200,1316,400), ref_phase=np.array(),
                                       weight_func='var', min_norm_velocity=True, mask_dataset_name='coherence')
                ifgram_inversion_patch('ifgramStack_001.h5', box=None, ref_phase=None,
                                       weight_func='var', min_norm_velocity=True, mask_dataset_name='coherence')
    """

    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)

    ## debug
    #y, x = 258, 454
    #box = (x, y, x+1, y+1)

    # Size Info - Patch
    if box:
        #print('processing \t %d-%d / %d lines ...' % (box[1], box[3], stack_obj.length))
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
                                 obsDatasetName=obsDatasetName,
                                 dropIfgram=True)

    pha_data = mask_unwrap_phase(pha_data,
                                 stack_obj,
                                 box,
                                 dropIfgram=True,
                                 mask_ds_name=mask_dataset_name,
                                 mask_threshold=mask_threshold)

    # Mask for pixels to invert
    mask = np.ones(num_pixel, np.bool_)
    # 1 - Water Mask
    if water_mask_file:
        print(('skip pixels on water with mask from'
               ' file: {}').format(os.path.basename(water_mask_file)))
        atr_msk = readfile.read_attribute(water_mask_file)
        if (int(atr_msk['LENGTH']), int(atr_msk['WIDTH'])) != (stack_obj.length, stack_obj.width):
            raise ValueError('Input water mask file has different size from ifgramStack file.')
        del atr_msk
        dsName = [i for i in readfile.get_dataset_list(water_mask_file)
                  if i in ['waterMask', 'mask']][0]
        waterMask = readfile.read(water_mask_file,
                                  datasetName=dsName,
                                  box=box)[0].flatten()
        mask *= np.array(waterMask, np.bool_)
        del waterMask

    # 2 - Mask for Zero Phase in ALL ifgrams
    if 'phase' in obsDatasetName.lower():
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
    print(('number of pixels to invert: {} out of {}'
           ' ({:.1f}%)').format(num_pixel2inv, num_pixel,
                                num_pixel2inv/num_pixel*100))
    if num_pixel2inv < 1:
        ts = ts.reshape(num_date, num_row, num_col)
        #ts_std = ts_std.reshape(num_date, num_row, num_col)
        temp_coh = temp_coh.reshape(num_row, num_col)
        num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)
        return ts, temp_coh, num_inv_ifg

    # skip zero value in the network inversion for phase
    if 'phase' in obsDatasetName.lower():
        skip_zero_value = True
    else:
        skip_zero_value = False

    # Inversion - SBAS
    if weight_func in ['no', 'sbas']:
        # Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        if 'phase' in obsDatasetName.lower():
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
        L = int(stack_obj.metadata['ALOOKS']) * int(stack_obj.metadata['RLOOKS'])
        weight = read_coherence(stack_obj, box=box, dropIfgram=True)
        weight = coherence2weight(weight, weight_func=weight_func, L=L, epsilon=5e-2)
        weight = np.sqrt(weight)

        # Weighted Inversion pixel by pixel
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
    if obsDatasetName.startswith('unwrapPhase'):
        phase2range = -1 * float(stack_obj.metadata['WAVELENGTH']) / (4.*np.pi)
        ts *= phase2range
        print('converting LOS phase displacement unit from radian to meter')

    elif obsDatasetName == 'azimuthOffset':
        az_pixel_size = ut.azimuth_ground_resolution(stack_obj.metadata)
        ts *= az_pixel_size
        print('converting azimuth offset displacement unit from pixel ({:.2f} m) to meter'.format(az_pixel_size))

    elif obsDatasetName == 'rangeOffset':
        rg_pixel_size = float(stack_obj.metadata['RANGE_PIXEL_SIZE'])
        ts *= rg_pixel_size
        print('converting range offset displacement unit from pixel ({:.2f} m) to meter'.format(rg_pixel_size))

    return ts, temp_coh, num_inv_ifg


def ifgram_inversion(ifgram_file='ifgramStack.h5', inps=None):
    """Implementation of the SBAS algorithm.
    Parameters: ifgram_file : string,
                    HDF5 file name of the interferograms stck
                inps : namespace, including the following options:
    Returns:    timeseriesFile : string
                    HDF5 file name of the output timeseries
                tempCohFile : string
                    HDF5 file name of temporal coherence
    Example:
        inps = cmd_line_parse()
        ifgram_inversion('ifgramStack.h5', inps)
    """
    start_time = time.time()

    # Check Inputs
    if not inps:
        inps = cmd_line_parse()

    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    A = stack_obj.get_design_matrix4timeseries(stack_obj.get_date12_list(dropIfgram=True))[0]
    num_ifgram, num_date = A.shape[0], A.shape[1]+1
    length, width = stack_obj.length, stack_obj.width
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
    box_list = split2boxes(dataset_shape=stack_obj.get_size(), chunk_size=inps.chunk_size)
    num_box = len(box_list)

    # read ifgram_file in small patches and write them together
    ref_phase = stack_obj.get_reference_phase(unwDatasetName=inps.obsDatasetName,
                                              skip_reference=inps.skip_ref,
                                              dropIfgram=True)

    # Initialization of output matrix
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    pbase = stack_obj.get_perp_baseline_timeseries(dropIfgram=True)
    date_list = stack_obj.get_date_list(dropIfgram=True)

    #ts     = np.zeros((num_date, length, width), np.float32)
    #ts_std = np.zeros((num_date, length, width), np.float32)
    temp_coh    = np.zeros((length, width), np.float32)
    num_inv_ifg = np.zeros((length, width), np.int16)

    # metadata
    metadata = dict(stack_obj.metadata)
    for key in configKeys:
        metadata[key_prefix+key] = str(vars(inps)[key])

    metadata['REF_DATE'] = date_list[0]
    metadata['FILE_TYPE'] = 'timeseries'
    metadata['UNIT'] = 'm'

    # Loop
    if not inps.parallel:
        # instantiate a timeseries object
        ts_file = '{}.h5'.format(os.path.splitext(inps.outfile[0])[0])
        ts_obj = timeseries(ts_file)

        # A dictionary of the datasets which we like to have in the timeseries
        dsNameDict = {
                "date": ((np.dtype('S8'), (num_date,))),
                "bperp": (np.float32, (num_date,)),
                "timeseries": (np.float32, (num_date, length, width)),
            }

        # layout the HDF5 file for the datasets and the metadata
        ts_obj.layout_hdf5(dsNameDict, metadata)

        # invert & write block by block
        for i in range(num_box):
            box = box_list[i]
            if num_box > 1:
                print('\n------- Processing Patch {} out of {} --------------'.format(i+1, num_box))

            # invert the network
            (tsi,
             temp_cohi,
             ifg_numi) = ifgram_inversion_patch(ifgram_file,
                                                box=box,
                                                ref_phase=ref_phase,
                                                obsDatasetName=inps.obsDatasetName,
                                                weight_func=inps.weightFunc,
                                                min_norm_velocity=inps.minNormVelocity,
                                                mask_dataset_name=inps.maskDataset,
                                                mask_threshold=inps.maskThreshold,
                                                min_redundancy=inps.minRedundancy,
                                                water_mask_file=inps.waterMaskFile)

            # write the block of timeseries to disk
            block = [0, num_date, box[1], box[3], box[0], box[2]]
            ts_obj.write2hdf5_block(tsi, datasetName='timeseries', block=block)

            # save the block of aux datasets
            temp_coh[box[1]:box[3], box[0]:box[2]] = temp_cohi
            num_inv_ifg[box[1]:box[3], box[0]:box[2]] = ifg_numi

        # write date and bperp to disk
        print('-'*50)
        date_list_utf8 = [dt.encode('utf-8') for dt in date_list]
        ts_obj.write2hdf5_block(date_list_utf8, datasetName='date')
        ts_obj.write2hdf5_block(pbase, datasetName='bperp')

    # Parallel loop
    else:
        try:
            from dask.distributed import Client, as_completed
            import mintpy.objects.cluster as cl
        except ImportError:
            raise ImportError('Cannot import dask.distributed!')

        ts = np.zeros((num_date, length, width), np.float32)

        # Look at the ~/.config/dask/mintpy.yaml file for Changing the Dask configuration defaults

        # This line submits NUM_WORKERS jobs to Pegasus to start a bunch of workers
        # In tests on Pegasus `general` queue in Jan 2019, no more than 40 workers could RUN
        # at once (other user's jobs gained higher priority in the general at that point)
        NUM_WORKERS = inps.numWorker

        # FA: the following command starts the jobs
        cluster = cl.get_cluster(cluster_type=inps.cluster, walltime=inps.walltime, config_name=inps.config)
        cluster.scale(NUM_WORKERS)
        print("JOB COMMAND CALLED FROM PYTHON:", cluster.job_script())
        with open('dask_command_run_from_python.txt', 'w') as f:
              f.write(cluster.job_script() + '\n')

        # This line needs to be in a function or in a `if __name__ == "__main__":` block. If it is in no function
        # or "main" block, each worker will try to create its own client (which is bad) when loading the module
        client = Client(cluster)

        all_boxes = []
        for box in box_list:
            # `box_list` is split into smaller boxes and then each box is processed in parallel
            # With larger jobs, increasing the `num_split` factor may improve runtime
            all_boxes += subsplit_boxes4_workers(box, num_split=1 * NUM_WORKERS, dimension='x')

        futures = []
        start_time_subboxes = time.time()
        for i, subbox in enumerate(all_boxes):
            print(i, subbox)

            data = (ifgram_file,
                    subbox,
                    ref_phase,
                    inps.obsDatasetName,
                    inps.weightFunc,
                    inps.minNormVelocity,
                    inps.maskDataset,
                    inps.maskThreshold,
                    inps.minRedundancy,
                    inps.waterMaskFile)

            # David: I haven't played with fussing with `retries`, however sometimes a future fails
            # on a worker for an unknown reason. retrying will save the whole process from failing.
            # TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
            # thrown in that case, therefore I don't know how to recognize when this happens.
            future = client.submit(parallel_ifgram_inversion_patch, data, retries=3)
            futures.append(future)

        # Some workers are slower than others. When #(futures to complete) < #workers, we should
        # investigate whether `Client.replicate(future)` will speed up work on the final futures.
        # It's a slight speedup, but depending how computationally complex each future is, this could be a
        # decent speedup
        i_future = 0
        for future, result in as_completed(futures, with_results=True):
            i_future += 1
            print("FUTURE #" + str(i_future), "complete in", time.time() - start_time_subboxes,
                  "seconds. Box:", subbox, "Time:", time.time())
            tsi, temp_cohi, ifg_numi, subbox = result

            ts[:, subbox[1]:subbox[3], subbox[0]:subbox[2]] = tsi
            #ts_std[:, subbox[1]:subbox[3], subbox[0]:subbox[2]] = ts_stdi
            temp_coh[subbox[1]:subbox[3], subbox[0]:subbox[2]] = temp_cohi
            num_inv_ifg[subbox[1]:subbox[3], subbox[0]:subbox[2]] = ifg_numi

        # Shut down Dask workers gracefully
        cluster.close()
        client.close()

        ut.move_dask_stdout_stderr_files()

    # reference pixel
    if not inps.skip_ref:
        ref_y = int(stack_obj.metadata['REF_Y'])
        ref_x = int(stack_obj.metadata['REF_X'])
        num_inv_ifg[ref_y, ref_x] = num_ifgram
        temp_coh[ref_y, ref_x] = 1.

    if inps.parallel:
        # for dask still use the old function to write.
        # consider to migrate to block-by-block writing, if HDF5 support multiple 
        write2hdf5_file(ifgram_file, metadata, ts, temp_coh, num_inv_ifg, suffix='', inps=inps)
    else:
        write2hdf5_auxFiles(metadata, temp_coh, num_inv_ifg, suffix='', inps=inps)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return


def parallel_ifgram_inversion_patch(data):
    """
    This is the starting point for Dask futures. Futures start executing code here.
    :param data:
    :return: The box
    """
    (ifgram_file,
     box, 
     ref_phase,
     obsDatasetName,
     weight_func,
     min_norm_velocity,
     mask_dataset_name,
     mask_threshold,
     min_redundancy,
     water_mask_file) = data
    print("BOX DIMS:", box)

    # This line is where all of the processing happens.
    (tsi,
     temp_cohi,
     ifg_numi) = ifgram_inversion_patch(ifgram_file,
                                        box=box,
                                        ref_phase=ref_phase,
                                        obsDatasetName=obsDatasetName,
                                        weight_func=weight_func,
                                        min_norm_velocity=min_norm_velocity,
                                        mask_dataset_name=mask_dataset_name,
                                        mask_threshold=mask_threshold,
                                        min_redundancy=min_redundancy,
                                        water_mask_file=water_mask_file)

    return tsi, temp_cohi, ifg_numi, box


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    # Network Inversion
    if inps.residualNorm == 'L2':
        ifgram_inversion(inps.ifgramStackFile, inps)
    else:
        raise NotImplementedError('L1 norm minimization is not fully tested.')
        #ut.timeseries_inversion_L1(inps.ifgramStackFile, inps.timeseriesFile)
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()
