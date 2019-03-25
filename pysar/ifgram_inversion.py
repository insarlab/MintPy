#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend import:
#     from pysar import ifgram_inversion as ifginv


import os
import sys
import re
import time
import argparse
import h5py
import numpy as np
from scipy import linalg   # more effieint than numpy.linalg
from scipy.special import gamma
from pysar.objects import ifgramStack, timeseries
from pysar.utils import readfile, writefile, ptime, utils as ut


# key configuration parameter name
key_prefix = 'pysar.networkInversion.'
configKeys = ['unwDatasetName',
              'numIfgram',
              'weightFunc',
              'maskDataset',
              'maskThreshold',
              'minRedundancy',
              'minNormVelocity']


################################################################################################
EXAMPLE = """example:
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -t pysarApp_template.txt --update
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -t pysarApp_template.txt --fast
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w var
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w fim
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w coh
"""

TEMPLATE = """
## Invert network of interferograms into time-series using weighted least sqaure (WLS) estimator.
## weighting options for least square inversion:
## 1) var - use inverse of covariance as weight (Guarnieri & Tebaldini, 2008, TGRS) [recommended]
## 2) fim - use Fisher Information Matrix as weight (Seymour & Cumming, 1994, IGARSS)
## 3) coh - use coherence as weight (Perissin & Wang, 2012, IEEE-TGRS)
## 4) no  - uniform weight
## mask options for unwrapPhase of each interferogram before inversion:
## 1) coherence        - mask out pixels with spatial coherence < maskThreshold [recommended]
## 2) connectComponent - mask out pixels with False/0 value
## 3) no               - no masking. [Recommended]
## Temporal coherence is calculated and used to generate final mask (Pepe & Lanari, 2006, IEEE-TGRS)
## SBAS (Berardino et al., 2002) = minNormVelocity (yes) + weightFunc (no)
pysar.networkInversion.weightFunc      = auto #[var / fim / coh / no], auto for var
pysar.networkInversion.maskDataset     = auto #[coherence / connectComponent / no], auto for no
pysar.networkInversion.maskThreshold   = auto #[0-1], auto for 0.4
pysar.networkInversion.minRedundancy   = auto #[1-inf], auto for 1.0, min num_ifgram for every SAR acquisition
pysar.networkInversion.waterMaskFile   = auto #[filename / no], auto for no
pysar.networkInversion.minNormVelocity = auto #[yes / no], auto for yes, min-norm deformation velocity or phase
pysar.networkInversion.residualNorm    = auto #[L2 ], auto for L2, norm minimization solution
pysar.networkInversion.minTempCoh      = auto #[0.0-1.0], auto for 0.7, min temporal coherence for mask
pysar.networkInversion.minNumPixel     = auto #[int > 0], auto for 100, min number of pixels in mask above
"""

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
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgramStackFile',
                        help='interferograms stack file to be inverted')
    parser.add_argument('-i','-d', '--dset', dest='unwDatasetName', type=str,
                        help='dataset name of unwrap phase in ifgram file to be used for inversion\n'
                             'e.g.: unwrapPhase, unwrapPhase_bridging, ...')
    parser.add_argument('--template', '-t', dest='templateFile',
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--ref-date', dest='ref_date',
                        help='Reference date, first date by default.')
    parser.add_argument('--mask-dset', dest='maskDataset',
                        help='dataset used to mask unwrapPhase, e.g. coherence, connectComponent')
    parser.add_argument('--mask-threshold', dest='maskThreshold', metavar='NUM', type=float, default=0.4,
                        help='threshold to generate mask when mask is coherence')
    parser.add_argument('--min-redundancy', dest='minRedundancy', metavar='NUM', type=float, default=1.0, 
                        help='minimum redundancy of interferograms for every SAR acquisition.')

    parser.add_argument('--weight-function', '-w', dest='weightFunc', default='no', choices={'var', 'fim', 'coh', 'no'},
                        help='function used to convert coherence to weight for inversion:\n' +
                        'var - inverse of phase variance due to temporal decorrelation\n' +
                        'fim - Fisher Information Matrix as weight' +
                        'coh - spatial coherence\n' +
                        'no  - no/uniform weight')
    parser.add_argument('--min-norm-velocity', dest='minNormVelocity', action='store_true',
                        help=('Enable inversion with minimum-norm deformation velocity,'
                              ' instead of minimum-norm deformation phase'))
    parser.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1', 'L2'],
                        help='Inverse method used to residual optimization, L1 or L2 norm minimization. Default: L2')

    parser.add_argument('--chunk-size', dest='chunk_size', type=float, default=100e6,
                        help='max number of data (= ifgram_num * num_row * num_col) to read per loop\n' +
                        'default: 0.2 G; adjust it according to your computer memory.')
    parser.add_argument('--parallel', dest='parallel', action='store_true',
                        help='Enable parallel processing for the pixelwise weighted inversion. [not working yet]')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',
                        help='Skip checking reference pixel value, for simulation testing.')
    parser.add_argument('-o', '--output', dest='outfile', nargs=2,
                        metavar=('TS_FILE', 'TCOH_FILE'), default=['timeseries.h5', 'temporalCoherence.h5'],
                        help='Output file name for timeseries and temporal coherence, default:\n' +
                        'timeseries.h5 temporalCoherence.h5')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n' +
                        'readable and newer than input interferograms file')
    parser.add_argument('--noskip-zero-phase', dest='skip_zero_phase', action='store_false',
                        help='Do not skip interferograms with zero phase.')
    parser.add_argument('--water-mask', '-m', dest='waterMaskFile',
                        help='Skip inversion on the masked out region, i.e. water.')
    parser.add_argument('--split-file', dest='split_file', action='store_true',
                        help='Split ifgramStack file into small files and invert them separately')
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
    if atr['FILE_TYPE'] != 'ifgramStack':
        raise ValueError('input is {} file, only support ifgramStack file.'.format(k))
    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    keyList = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in keyList:
        value = template[key_prefix+key]
        if key in ['maskDataset', 'minNormVelocity']:
            inpsDict[key] = value
        elif value:
            if key in ['maskThreshold', 'minRedundancy']:
                inpsDict[key] = float(value)
            elif key in ['weightFunc', 'residualNorm', 'waterMaskFile']:
                inpsDict[key] = value
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
            ti = float(f[inps.unwDatasetName].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgramStackFile)))
        to = min(os.path.getmtime(i) for i in inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output files are NOT newer than input dataset: {}.'.format(inps.unwDatasetName))
        else:
            print('2) output dataset is newer than input dataset: {}.'.format(inps.unwDatasetName))

    # check configuration
    if flag == 'skip':
        atr = readfile.read_attribute(inps.timeseriesFile)
        if any(str(vars(inps)[key]) != atr.get(key_prefix+key, 'None') for key in configKeys):
            flag = 'run'
            print('3) NOT all key configration parameters are the same: {}'.format(configKeys))
        else:
            print('3) all key configuration parameters are the same: {}.'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


################################################################################################
def phase_pdf_ds(L, coherence=None, phi_num=1000, epsilon=1e-3):
    """Marginal PDF of interferometric phase for distributed scatterers (DS)
    Eq. 66 (Tough et al., 1995) and Eq. 4.2.23 (Hanssen, 2001)
    Inputs:
        L         - int, number of independent looks
        coherence - 1D np.array for the range of coherence, with value < 1.0 for valid operation
        phi_num    - int, number of phase sample for the numerical calculation
    Output:
        pdf       - 2D np.array, phase pdf in size of (phi_num, len(coherence))
        coherence - 1D np.array for the range of coherence
    Example:
        epsilon = 1e-4
        coh = np.linspace(0., 1-epsilon, 1000)
        pdf, coh = phase_pdf_ds(1, coherence=coh)
    """
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000)
    coherence = np.array(coherence, np.float64).reshape(1, -1)
    phi = np.linspace(-np.pi, np.pi, phi_num, dtype=np.float64).reshape(-1, 1)

    # Phase PDF - Eq. 4.2.32 (Hanssen, 2001)
    A = np.power((1-np.square(coherence)), L) / (2*np.pi)
    A = np.tile(A, (phi_num, 1))
    B = gamma(2*L - 1) / ((gamma(L))**2 * 2**(2*(L-1)))

    beta = np.multiply(np.abs(coherence), np.cos(phi))
    C = np.divide((2*L - 1) * beta, np.power((1 - np.square(beta)), L+0.5))
    C = np.multiply(C, (np.pi/2 + np.arcsin(beta)))
    C += 1 / np.power((1 - np.square(beta)), L)

    sumD = 0
    if L > 1:
        for r in range(int(L)-1):
            D = gamma(L-0.5) / gamma(L-0.5-r)
            D *= gamma(L-1-r) / gamma(L-1)
            D *= (1 + (2*r+1)*np.square(beta)) / np.power((1 - np.square(beta)), r+2)
            sumD += D
        sumD /= (2*(L-1))

    pdf = B*C + sumD
    pdf = np.multiply(A, pdf)
    return pdf, coherence.flatten()


def phase_variance_ds(L,  coherence=None, epsilon=1e-3):
    """Interferometric phase variance for distributed scatterers (DS)
    Eq. 2.1.2 (Box et al., 2015) and Eq. 4.2.27 (Hanssen, 2001)
    Inputs:
        L         - int, number of independent looks
        coherence - 1D np.array for the range of coherence, with value < 1.0 for valid operation
        phiNum    - int, number of phase sample for the numerical calculation
    Output:
        var       - 1D np.array, phase variance in size of (len(coherence))
        coherence - 1D np.array for the range of coherence
    Example:
        epsilon = 1e-4
        coh = np.linspace(0., 1-epsilon, 1000)
        var, coh = phase_variance_ds(1, coherence=coh)
    """
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000, dtype=np.float64)
    phiNum = len(coherence)

    phi = np.linspace(-np.pi, np.pi, phiNum, dtype=np.float64).reshape(-1, 1)
    phi_step = 2*np.pi/phiNum

    pdf, coherence = phase_pdf_ds(L, coherence=coherence)
    var = np.sum(np.multiply(np.square(np.tile(phi, (1, len(coherence)))), pdf)*phi_step, axis=0)
    return var, coherence


def phase_variance_ps(L, coherence=None, epsilon=1e-3):
    """the Cramer-Rao bound (CRB) of phase variance
    Given by Eq. 25 (Rodriguez and Martin, 1992)and Eq 4.2.32 (Hanssen, 2001)
    Valid when coherence is close to 1.
    """
    if coherence is None:
        coherence = np.linspace(0.9, 1.-epsilon, 1000, dtype=np.float64)
    var = (1 - coherence**2) / (2 * int(L) * coherence**2)
    return var, coherence


def coherence2phase_variance_ds(coherence, L=32, epsilon=1e-3, print_msg=False):
    """Convert coherence to phase variance based on DS phase PDF (Tough et al., 1995)"""
    lineStr = '    number of looks L={}'.format(L)
    if L > 80:
        L = 80
        lineStr += ', use L=80 to avoid dividing by 0 in calculation with Negligible effect'
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

    var_lut = phase_variance_ds(int(L), coh_lut)[0]
    variance = var_lut[coherence_idx]
    return variance


def coherence2fisher_info_index(data, L=32, epsilon=1e-3):
    """Convert coherence to Fisher information index (Seymour & Cumming, 1994, IGARSS)"""
    if data.dtype != np.float64:
        data = np.array(data, np.float64)
    data[data > 1-epsilon] = 1-epsilon
    data = 2.0 * int(L) * np.square(data) / (1 - np.square(data))
    return data


def estimate_timeseries(A, B, tbase_diff, ifgram, weight_sqrt=None, min_norm_velocity=True,
                        skip_zero_phase=True, rcond=1e-5, min_redundancy=1.):
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
                skip_zero_phase - bool, skip ifgram with zero phase value
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
    if skip_zero_phase and not np.all(ifgram):
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


###########################################################################################
def write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std=None, num_inv_ifg=None,
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
    print('converting phase to range')
    phase2range = -1*float(stack_obj.metadata['WAVELENGTH'])/(4.*np.pi)
    ts *= phase2range

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

    # File 3 - timeseriesDecorStd.h5
    if not np.all(ts_std == 0.):
        out_file = 'timeseriesDecorStd{}.h5'.format(suffix)
        metadata['FILE_TYPE'] = 'timeseries'
        metadata['UNIT'] = 'm'
        ts_std *= abs(phase2range)
        print('-'*50)
        writefile.write(ts_std, out_file=out_file, metadata=metadata, ref_file=ts_file)

    # File 4 - numInvIfgram.h5
    out_file = 'numInvIfgram{}.h5'.format(suffix)
    metadata['FILE_TYPE'] = 'mask'
    metadata['UNIT'] = '1'
    print('-'*50)
    writefile.write(num_inv_ifg, out_file=out_file, metadata=metadata)
    return


def split_ifgram_file(ifgram_file, chunk_size=100e6):
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    metadata = dict(stack_obj.metadata)

    # get reference phase
    ref_phase = stack_obj.get_reference_phase(dropIfgram=False)

    # get list of boxes
    box_list = split_into_boxes(dataset_shape=stack_obj.get_size(),
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


def split_into_boxes(dataset_shape, chunk_size=100e6, print_msg=True):
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


def subsplit_boxes_by_num_workers(box, num_workers, dimension='y'):
    """ David: This is a bit hacky, but after creating the patches,
    I wanted to further divide the box size into `num_subboxes` different subboxes.
    Note that `split_into_boxes` only splits based on chunk_size (memory-based).

    :param box: [x0, y0, x1, y1]: list[int] of size 4,
    """

    # Flip x and y coordinates if splitting along 'x' dimension
    x0, y0, x1, y1 = box
    subboxes = []

    if dimension == 'y':
        y_diff = y1 - y0
        for i in range(num_workers):
            start = (i * y_diff) // num_workers
            end = ((i + 1) * y_diff) // num_workers if i + 1 != num_workers else y_diff
            subboxes.append([x0, start, x1, end])
    elif dimension == 'x':
        x_diff = x1 - x0
        for i in range(num_workers):
            start = (i * x_diff) // num_workers
            end = ((i + 1) * x_diff) // num_workers if i + 1 != num_workers else x_diff
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


def read_unwrap_phase(stack_obj, box, ref_phase, unwDatasetName='unwrapPhase', dropIfgram=True,
                      skip_zero_phase=True, print_msg=True):
    """Read unwrapPhase from ifgramStack file
    Parameters: stack_obj : ifgramStack object
                box : tuple of 4 int
                ref_phase : 1D array or None
                skip_zero_phase : bool
    Returns:    pha_data : 2D array of unwrapPhase in size of (num_ifgram, num_pixel)
    """
    # Read unwrapPhase
    num_ifgram = stack_obj.get_size(dropIfgram=dropIfgram)[0]
    if print_msg:
        print('reading {} in {} * {} ...'.format(unwDatasetName, box, num_ifgram))
    pha_data = stack_obj.read(datasetName=unwDatasetName,
                              box=box,
                              dropIfgram=dropIfgram,
                              print_msg=False).reshape(num_ifgram, -1)

    # read ref_phase
    if ref_phase is not None:
        # use input ref_phase array (for split_file=False)
        if print_msg:
            print('use input reference phase')
    elif 'refPhase' in stack_obj.datasetNames:
        # read ref_phase from file itself (for split_file=True)
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
        if mask_ds_name == 'coherence':
            msk_data = msk_data >= mask_threshold
            if print_msg:
                print('mask out pixels with {} < {}'.format(mask_ds_name, mask_threshold))
        else:
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
    return coh_data


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
        weight = 1.0 / coherence2phase_variance_ds(coh_data, L, print_msg=print_msg)

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


def ifgram_inversion_patch(ifgram_file, box=None, ref_phase=None, unwDatasetName='unwrapPhase',
                           weight_func='var', min_norm_velocity=True,
                           mask_dataset_name=None, mask_threshold=0.4, min_redundancy=1.0,
                           water_mask_file=None, skip_zero_phase=True):
    """Invert one patch of an ifgram stack into timeseries.
    Parameters: ifgram_file       : str, interferograms stack HDF5 file, e.g. ./INPUTS/ifgramStack.h5
                box               : tuple of 4 int, indicating (x0, y0, x1, y1) pixel coordinate of area of interest
                                    or None, to process the whole file and write output file
                ref_phase         : 1D array in size of (num_ifgram) 
                                    or None
                weight_func       : str, weight function, choose in ['no', 'fim', 'var', 'coh']
                mask_dataset_name : str, dataset name in ifgram_file used to mask unwrapPhase pixelwisely
                mask_threshold    : float, min coherence of pixels if mask_dataset_name='coherence'
                water_mask_file   : str, water mask filename if available,
                                    skip inversion on water to speed up the process
                skip_zero_phase   : bool, skip zero value of unwrapped phase or not, default yes, for comparison
    Returns:    ts             : 3D array in size of (num_date, num_row, num_col)
                temp_coh       : 2D array in size of (num_row, num_col)
                ts_std         : 3D array in size of (num_date, num_row, num_col)
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
    num_ifgram = len(date12_list)

    # prep for decor std time-series
    try:
        ref_date = str(np.loadtxt('reference_date.txt', dtype=bytes).astype(str))
    except:
        ref_date = date_list[0]
    Astd = stack_obj.get_design_matrix4timeseries(date12_list=date12_list, refDate=ref_date)[0]
    #ref_idx = date_list.index(ref_date)
    #time_idx = [i for i in range(num_date)]
    #time_idx.remove(ref_idx)

    # Initialization of output matrix
    ts     = np.zeros((num_date, num_pixel), np.float32)
    ts_std = np.zeros((num_date, num_pixel), np.float32)
    temp_coh    = np.zeros(num_pixel, np.float32)
    num_inv_ifg = np.zeros(num_pixel, np.int16)

    # Read/Mask unwrapPhase
    pha_data = read_unwrap_phase(stack_obj,
                                 box,
                                 ref_phase,
                                 unwDatasetName=unwDatasetName,
                                 dropIfgram=True,
                                 skip_zero_phase=skip_zero_phase)

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
    print('skip pixels with zero/nan value in all interferograms')
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
        ts_std = ts_std.reshape(num_date, num_row, num_col)
        temp_coh = temp_coh.reshape(num_row, num_col)
        num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)
        return ts, temp_coh, ts_std, num_inv_ifg

    # Inversion - SBAS
    if weight_func in ['no', 'sbas']:
        # Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        mask_all_net = np.all(pha_data, axis=0)
        mask_all_net *= mask
        mask_part_net = mask ^ mask_all_net

        if np.sum(mask_all_net) > 0:
            print(('inverting pixels with valid phase in all  ifgrams'
                   ' ({:.0f} pixels) ...').format(np.sum(mask_all_net)))
            tsi, tcohi, num_ifgi = estimate_timeseries(A, B, tbase_diff, 
                                                       ifgram=pha_data[:, mask_all_net],
                                                       weight_sqrt=None,
                                                       min_norm_velocity=min_norm_velocity,
                                                       skip_zero_phase=skip_zero_phase,
                                                       min_redundancy=min_redundancy)
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
                                                           skip_zero_phase=skip_zero_phase,
                                                           min_redundancy=min_redundancy)
                ts[:, idx] = tsi.flatten()
                temp_coh[idx] = tcohi
                num_inv_ifg[idx] = num_ifgi
                prog_bar.update(i+1, every=1000, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
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
                                                       skip_zero_phase=skip_zero_phase,
                                                       min_redundancy=min_redundancy)
            ts[:, idx] = tsi.flatten()
            temp_coh[idx] = tcohi
            num_inv_ifg[idx] = num_ifgi
            prog_bar.update(i+1, every=1000, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
        prog_bar.close()

    ts = ts.reshape(num_date, num_row, num_col)
    ts_std = ts_std.reshape(num_date, num_row, num_col)
    temp_coh = temp_coh.reshape(num_row, num_col)
    num_inv_ifg = num_inv_ifg.reshape(num_row, num_col)

    # write output files if input file is splitted (box == None)
    if box is None:
        # metadata
        metadata = dict(stack_obj.metadata)
        metadata[key_prefix+'weightFunc'] = weight_func
        suffix = re.findall('_\d{3}', ifgram_file)[0]
        write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std, num_inv_ifg, suffix, inps=inps)
        return
    else:
        return ts, temp_coh, ts_std, num_inv_ifg


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
    #msg += '\tLS  for pixels with full rank      network\n'
    #msg += '\tSVD for pixels with rank deficient network\n'
    msg += 'minimum redundancy: {}\n'.format(inps.minRedundancy)
    msg += 'weight function: {}\n'.format(inps.weightFunc)

    if inps.maskDataset:
        if inps.maskDataset == 'coherence':
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
    box_list = split_into_boxes(dataset_shape=stack_obj.get_size(), chunk_size=inps.chunk_size)
    num_box = len(box_list)
    if inps.split_file:
        # split ifgram_file into small files and write each of them
        print('\n---------------------------- Splitting Input File -----------------------------')
        ifgram_files = split_ifgram_file(ifgram_file, chunk_size=inps.chunk_size)
        num_file = len(ifgram_files)

        # Loop
        for fname in ifgram_files:
            if num_file > 1:
                print('\n------- Processing {} ({} in total) --------------'.format(fname, num_file))
            ifgram_inversion_patch(fname,
                                   box=None,
                                   ref_phase=None,
                                   unwDatasetName=inps.unwDatasetName,
                                   weight_func=inps.weightFunc,
                                   min_norm_velocity=inps.minNormVelocity,
                                   mask_dataset_name=inps.maskDataset,
                                   mask_threshold=inps.maskThreshold,
                                   min_redundancy=inps.minRedundancy,
                                   water_mask_file=inps.waterMaskFile,
                                   skip_zero_phase=inps.skip_zero_phase)
    else:
        # read ifgram_file in small patches and write them together
        ref_phase = stack_obj.get_reference_phase(unwDatasetName=inps.unwDatasetName,
                                                  skip_reference=inps.skip_ref,
                                                  dropIfgram=True)

        # Initialization of output matrix
        stack_obj = ifgramStack(ifgram_file)
        stack_obj.open(print_msg=False)
        ts     = np.zeros((num_date, length, width), np.float32)
        ts_std = np.zeros((num_date, length, width), np.float32)
        temp_coh    = np.zeros((length, width), np.float32)
        num_inv_ifg = np.zeros((length, width), np.int16)

        # Loop
        if inps.parallel:
            from dask.distributed import Client, as_completed
            # David: dask_jobqueue is needed for HPC.
            # PBSCluster (similar to LSFCluster) should also work out of the box
            from dask_jobqueue import LSFCluster

            # Initialize Dask Workers.
            # TODO: Should these params be moved into a config file?
            # We could use our own config setup or Dask's config setup
            NUM_WORKERS = 40
            cluster = LSFCluster(project='insarlab',
                                 name='pysar_worker_bee2',
                                 queue='general',
                                 # David: The first parameter is required by Pegasus. This actually changes memory usage.
                                 job_extra=['-R "rusage[mem=6400]"',
                                            # David: This second line will allow you to read your worker's output
                                            "-o WORKER-%J.out"],
                                 # David: This allows workers to write overflow memory to file. Unsure if necessary
                                 local_directory='/scratch/projects/insarlab/dwg11/',
                                 # David: Somehow, cores=2 is essential. When cores=1 it was super slow and
                                 # when cores=3 it took forever for workers to go from PENDING to RUNNING
                                 cores=2,
                                 walltime='00:30',
                                 # David: This line is required by dask_jobqueue
                                 # but ignored by Pegasus as far as I can tell
                                 memory='2GB',
                                 # David: This line prevents the cluster from defaulting to the system Python.
                                 # You need to point it to your dask environment's Python
                                 python='/nethome/dwg11/anaconda2/envs/pysar_parallel/bin/python')
            # David:  This line submits NUM_WORKERS number of jobs to Pegasus to start a bunch of workers
            cluster.scale(NUM_WORKERS)
            print("JOB FILE:", cluster.job_script())

            # David: This line needs to be in an `if __name__ == "__main__":` block I believe. It should not
            # be floating around or in code that workers might execute (note that if it is in no function,
            # workers will execute it when loading the module)
            #
            client = Client(cluster)

            all_boxes = []
            for i in range(num_box):
                # David: All "patches" are processed at once
                # David: You can change the factor of `num_subboxes` to have smaller jobs per worker.
                # This could be helpful with large jobs
                all_boxes += subsplit_boxes_by_num_workers(box_list[i], num_workers=1 * NUM_WORKERS, dimension='x')

            futures = []
            start_time_subboxes = time.time()
            for i, subbox in enumerate(all_boxes):
                print(i, subbox)

                data = (ifgram_file,
                        subbox,
                        ref_phase,
                        inps.unwDatasetName,
                        inps.weightFunc,
                        inps.minNormVelocity,
                        inps.maskDataset,
                        inps.maskThreshold,
                        inps.minRedundancy,
                        inps.waterMaskFile,
                        inps.skip_zero_phase)


                # David: I haven't played with fussing with `retries`, however sometimes a future fails
                # on a worker for an unknown reason. retrying will save the whole process from failing.
                # David: TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
                # thrown in that case, therefore I don't know how to recognize when this happens.
                future = client.submit(parallel_ifgram_inversion_patch, data, retries=3)
                futures.append(future)

            # David: Some workers are slower than others. When #(futures to complete) < #workers, we should
            # investigate whether `Client.replicate(future)` will speed up work on the final futures.
            # It's a slight speedup, but depending how computationally complex each future is, this could be a
            # decent speedup
            i_future = 0
            for future, result in as_completed(futures, with_results=True):
                i_future += 1
                print("FUTURE #" + str(i_future), "complete in", time.time() - start_time_subboxes,
                      "seconds. Box:", subbox, "Time:", time.time())
                tsi, temp_cohi, ts_stdi, ifg_numi, subbox = result

                ts[:, subbox[1]:subbox[3], subbox[0]:subbox[2]] = tsi
                ts_std[:, subbox[1]:subbox[3], subbox[0]:subbox[2]] = ts_stdi
                temp_coh[subbox[1]:subbox[3], subbox[0]:subbox[2]] = temp_cohi
                num_inv_ifg[subbox[1]:subbox[3], subbox[0]:subbox[2]] = ifg_numi

            # Shut down Dask workers gracefully
            client.close()
            cluster.close()
        else:
            for i in range(num_box):
                box = box_list[i]
                if num_box > 1:
                    print('\n------- Processing Patch {} out of {} --------------'.format(i+1, num_box))
                (tsi,
                 temp_cohi,
                 ts_stdi,
                 ifg_numi) = ifgram_inversion_patch(ifgram_file,
                                                    box=box,
                                                    ref_phase=ref_phase,
                                                    unwDatasetName=inps.unwDatasetName,
                                                    weight_func=inps.weightFunc,
                                                    min_norm_velocity=inps.minNormVelocity,
                                                    mask_dataset_name=inps.maskDataset,
                                                    mask_threshold=inps.maskThreshold,
                                                    min_redundancy=inps.minRedundancy,
                                                    water_mask_file=inps.waterMaskFile,
                                                    skip_zero_phase=inps.skip_zero_phase)

                ts[:, box[1]:box[3], box[0]:box[2]] = tsi
                ts_std[:, box[1]:box[3], box[0]:box[2]] = ts_stdi
                temp_coh[box[1]:box[3], box[0]:box[2]] = temp_cohi
                num_inv_ifg[box[1]:box[3], box[0]:box[2]] = ifg_numi

        # reference pixel
        ref_y = int(stack_obj.metadata['REF_Y'])
        ref_x = int(stack_obj.metadata['REF_X'])
        num_inv_ifg[ref_y, ref_x] = num_ifgram
        temp_coh[ref_y, ref_x] = 1.

        # metadata
        metadata = dict(stack_obj.metadata)
        for key in configKeys:
            metadata[key_prefix+key] = str(vars(inps)[key])
        write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std, num_inv_ifg, suffix='', inps=inps)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return


def parallel_ifgram_inversion_patch(data):
    """
    David: This is the starting point for futures. Futures start executing code here.
    :param data:
    :return: The box
    """
    (ifgram_file, box, ref_phase, unwDatasetName,
     weight_func, min_norm_velocity,
     mask_dataset_name, mask_threshold,
     min_redundancy, water_mask_file, skip_zero_phase) = data

    print("BOX DIMS:", box)

    # David: This is the main call. This code was copied from the old
    # `ifgram_inversion` function
    (tsi, temp_cohi, ts_stdi, ifg_numi) = ifgram_inversion_patch(ifgram_file,
                                           box= box,
                                           ref_phase= ref_phase,
                                           unwDatasetName=unwDatasetName,
                                           weight_func=weight_func,
                                           min_norm_velocity=min_norm_velocity,
                                           mask_dataset_name=mask_dataset_name,
                                           mask_threshold=mask_threshold,
                                           min_redundancy=min_redundancy,
                                           water_mask_file=water_mask_file,
                                           skip_zero_phase=skip_zero_phase)

    return tsi, temp_cohi, ts_stdi, ifg_numi, box


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.templateFile:
        inps = read_template2inps(inps.templateFile, inps)
    inps.timeseriesFile, inps.tempCohFile = inps.outfile

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
    if not inps.unwDatasetName:
        stack_obj = ifgramStack(inps.ifgramStackFile)
        stack_obj.open(print_msg=False)
        inps.unwDatasetName = [i for i in ['unwrapPhase_bridging_phaseClosure',
                                           'unwrapPhase_bridging',
                                           'unwrapPhase_phaseClosure',
                                           'unwrapPhase']
                               if i in stack_obj.datasetNames][0]

    # --update option
    inps.numIfgram = len(ifgramStack(inps.ifgramStackFile).get_date12_list(dropIfgram=True))
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
