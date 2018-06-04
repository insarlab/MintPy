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
from scipy.special import gamma
from pysar.objects import ifgramStack, timeseries
from pysar.utils import readfile, writefile, ptime, utils as ut

key_prefix = 'pysar.networkInversion.'


################################################################################################
EXAMPLE = """example:
  ifgram_inversion.py  INPUTS/ifgramStack.h5
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -t pysarApp_template.txt
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w var
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w fim
  ifgram_inversion.py  INPUTS/ifgramStack.h5 -w coh
"""

TEMPLATE = """
## Invert network of interferograms into time series using weighted least sqaure (WLS) estimator.
## mask options for unwrapPhase of each interferogram before inversion:
## 1) coherence        - mask out pixels with spatial coherence < maskThreshold [Recommended]
## 2) connectComponent - mask out pixels with False/0 value
## 3) no               - no masking.
## weighting options for least square inversion:
## 1) fim  - WLS, use Fisher Information Matrix as weight (Seymour & Cumming, 1994, IGARSS). [Recommended]
## 2) var  - WLS, use inverse of covariance as weight (Guarnieri & Tebaldini, 2008, TGRS)
## 3) coh  - WLS, use coherence as weight (Perissin & Wang, 2012, IEEE-TGRS)
## 4) sbas - LS/SVD, uniform weight (Berardino et al., 2002, TGRS)
## Temporal coherence is calculated and used to generate final mask (Pepe & Lanari, 2006, IEEE-TGRS)
pysar.networkInversion.weightFunc    = auto #[fim / var / coh / sbas], auto for fim
pysar.networkInversion.maskDataset   = auto #[coherence / connectComponent / no], auto for coherence
pysar.networkInversion.maskThreshold = auto #[0-1], auto for 0.4
pysar.networkInversion.waterMaskFile = auto #[filename / no], auto for no
pysar.networkInversion.residualNorm  = auto #[L2 ], auto for L2, norm minimization solution
pysar.networkInversion.minTempCoh    = auto #[0.0-1.0], auto for 0.7, min temporal coherence for mask
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
    parser = argparse.ArgumentParser(description='Invert network of interferograms into timeseries.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgramStackFile',
                        help='interferograms stack file to be inverted')
    parser.add_argument('--template', '-t', dest='templateFile',
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--ref-date', dest='ref_date',
                        help='Reference date, first date by default.')
    parser.add_argument('--mask-dset', dest='maskDataset',
                        help='dataset used to mask unwrapPhase, e.g. coherence, connectComponent')
    parser.add_argument('--mask-threshold', dest='maskThreshold', type=float, default=0.4,
                        help='threshold to generate mask when mask is coherence')

    parser.add_argument('--weight-function', '-w', dest='weightFunc', default='sbas', choices={'fim', 'var', 'coh', 'sbas'},
                        help='function used to convert coherence to weight for inversion:\n' +
                        'fim  - Fisher Information Matrix as weight' +
                        'var  - phase variance due to temporal decorrelation\n' +
                        'coh  - uniform distribution CDF function\n' +
                        'sbas - uniform weight')
    parser.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1', 'L2'],
                        help='Inverse method used to residual optimization, L1 or L2 norm minimization. Default: L2')

    parser.add_argument('--chunk-size', dest='chunk_size', type=float, default=100e6,
                        help='max number of data (= ifgram_num * num_row * num_col) to read per loop\n' +
                        'default: 0.2 G; adjust it according to your computer memory.')
    parser.add_argument('--parallel', dest='parallel', action='store_true',
                        help='Enable parallel processing for the pixelwise weighted inversion. [not working yet]')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',
                        help='Skip checking reference pixel value, for simulation testing.')
    parser.add_argument('-o', '--output', dest='outfile', nargs=2, default=['timeseries.h5', 'temporalCoherence.h5'],
                        help='Output file name for timeseries and temporal coherence, default:\n' +
                        'timeseries.h5 temporalCoherence.h5')
    parser.add_argument('--update-mode', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n' +
                        'readable and newer than input interferograms file')
    parser.add_argument('--noskip-zero-phase', dest='skip_zero_phase', action='store_false',
                        help='Do not skip interferograms with zero phase.')
    parser.add_argument('--water-mask', '-m', dest='waterMaskFile',
                        help='Skip inversion on the masked out region, i.e. water.')
    parser.add_argument('--split-file', dest='split_file', action='store_true',
                        help='Split ifgramStack file into small files and invert them separately')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.parallel = False
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
        if key in ['maskDataset']:
            inpsDict[key] = value
        elif value:
            if key in ['maskThreshold']:
                inpsDict[key] = float(value)
            elif key in ['weightFunc', 'residualNorm', 'waterMaskFile']:
                inpsDict[key] = value
    return inps


################################################################################################
def phase_pdf_ds(l, coherence=None, phi_num=1000):
    """Marginal PDF of interferometric phase for distributed scatterers (DS)
    Eq. 66 (Tough et al., 1995) and Eq. 4.2.23 (Hanssen, 2001)
    Inputs:
        l         - int, number of independent looks
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
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000)
    coherence = np.array(coherence, np.float64).reshape(1, -1)
    phi = np.linspace(-np.pi, np.pi, phi_num, dtype=np.float64).reshape(-1, 1)

    # Phase PDF - Eq. 4.2.32 (Hanssen, 2001)
    A = np.power((1-np.square(coherence)), l) / (2*np.pi)
    A = np.tile(A, (phi_num, 1))
    B = gamma(2*l - 1) / ((gamma(l))**2 * 2**(2*(l-1)))

    beta = np.multiply(np.abs(coherence), np.cos(phi))
    # C1 = np.power((1 - np.square(beta)), l+0.5)
    # C1[C1 == 0.] = epsilon
    # C = np.divide((2*l - 1) * beta, C1)
    C = np.divide((2*l - 1) * beta, np.power((1 - np.square(beta)), l+0.5))
    C = np.multiply(C, (np.pi/2 + np.arcsin(beta)))
    # C2 = np.power((1 - np.square(beta)), l)
    # C2[C2 == 0.0] = epsilon
    # C += 1 / C2
    C += 1 / np.power((1 - np.square(beta)), l)

    sumD = 0
    if l > 1:
        for r in range(l-1):
            D = gamma(l-0.5) / gamma(l-0.5-r)
            D *= gamma(l-1-r) / gamma(l-1)
            # D1 = np.power((1 - np.square(beta)), r+2)
            # D1[D1 == 0.] = epsilon
            # D *= (1 + (2*r+1)*np.square(beta)) / D1
            D *= (1 + (2*r+1)*np.square(beta)) / np.power((1 - np.square(beta)), r+2)
            sumD += D
        sumD /= (2*(l-1))

    pdf = B*C + sumD
    pdf = np.multiply(A, pdf)
    return pdf, coherence.flatten()


def phase_variance_ds(l,  coherence=None):
    """Interferometric phase variance for distributed scatterers (DS)
    Eq. 2.1.2 (Box et al., 2015) and Eq. 4.2.27 (Hanssen, 2001)
    Inputs:
        l         - int, number of independent looks
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
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000, dtype=np.float64)
    phiNum = len(coherence)

    phi = np.linspace(-np.pi, np.pi, phiNum, dtype=np.float64).reshape(-1, 1)
    phi_step = 2*np.pi/phiNum

    pdf, coherence = phase_pdf_ds(l, coherence=coherence)
    var = np.sum(np.multiply(np.square(np.tile(phi, (1, len(coherence)))), pdf)*phi_step, axis=0)
    return var, coherence


def phase_variance_ps(L, coherence=None):
    """the Cramer-Rao bound (CRB) of phase variance
    Given by Eq. 25 (Rodriguez and Martin, 1992)and Eq 4.2.32 (Hanssen, 2001)
    Valid when coherence is close to 1.
    """
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0.9, 1.-epsilon, 1000, dtype=np.float64)
    var = (1-coherence**2) / (2*L*coherence**2)
    return var, coherence


def coherence2phase_variance_ds(coherence, L=32, print_msg=False):
    """Convert coherence to phase variance based on DS phase PDF (Tough et al., 1995)"""
    lineStr = '    number of multilooks L=%d' % L
    if L > 80:
        L = 80
        lineStr += ', use L=80 to avoid dividing by 0 in calculation with Negligible effect'
    if print_msg:
        print(lineStr)

    epsilon = 1e-4
    coh_num = 1000
    coh_min = 0.0
    coh_max = 1.0 - epsilon
    coh_lut = np.linspace(coh_min, coh_max, coh_num)
    coh_min = np.min(coh_lut)
    coh_max = np.max(coh_lut)
    coh_step = (coh_max - coh_min) / (coh_num - 1)

    coherence = np.array(coherence)
    coherence[coherence < coh_min] = coh_min
    coherence[coherence > coh_max] = coh_max
    coherence_idx = np.array((coherence - coh_min) / coh_step, np.int16)

    var_lut = phase_variance_ds(L, coh_lut)[0]
    variance = var_lut[coherence_idx]
    return variance


def coherence2fisher_info_index(data, L=32, epsilon=1e-4):
    """Convert coherence to Fisher information index (Seymour & Cumming, 1994, IGARSS)"""
    if data.dtype != np.float64:
        data = np.array(data, np.float64)
    data[data > 1-epsilon] = 1-epsilon
    data = 2.0 * L * np.square(data) / (1 - np.square(data))
    return data


def round_to_1(x):
    """Return the most significant digit of input number"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)


def ceil_to_1(x):
    """Return the most significant digit of input number and ceiling it"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)+10**digit


def network_inversion_sbas(B, ifgram, tbase_diff, skip_zero_phase=True):
    """ Network inversion based on Small BAseline Subsets (SBAS) algorithm (Berardino et al.,
        2002, IEEE-TGRS). For full rank design matrix, a.k.a., fully connected network, ordinary
        least square (OLS) inversion is applied; otherwise, Singular Value Decomposition (SVD).

    Inputs:
        B          - 2D np.array in size of (ifgram_num, num_date-1)
                     design matrix B, which represents temporal baseline timeseries between
                     master and slave date for each interferogram
        ifgram     - 2D np.array in size of (ifgram_num, num_pixel)
                     phase of all interferograms
        tbase_diff - 2D np.array in size of (num_date-1, 1)
                     differential temporal baseline of time-series
        skip_zero_phase - bool, skip ifgram with zero phase value
    Output:
        ts      - 2D np.array in size of (num_date-1, num_pixel), phase time series
        temp_coh - 1D np.array in size of (num_pixel), temporal coherence
    """
    ifgram = ifgram.reshape(B.shape[0], -1)
    dateNum1 = B.shape[1]
    ts = np.zeros(dateNum1, np.float32)
    temp_coh = 0.
    ifg_num = 0

    # Skip Zero Phase Value
    if skip_zero_phase and not np.all(ifgram):
        idx = (ifgram != 0.).flatten()
        B = B[idx, :]
        if B.shape[0] < dateNum1*2:
            return ts, temp_coh, ifg_num
        ifgram = ifgram[idx, :]

    try:
        # Invert time-series
        B_inv = np.array(np.linalg.pinv(B), np.float32)
        ts_rate = np.dot(B_inv, ifgram)
        ts_diff = ts_rate * np.tile(tbase_diff, (1, ifgram.shape[1]))
        ts = np.cumsum(ts_diff, axis=0)

        # Temporal Coherence
        ifgram_diff = ifgram - np.dot(B, ts_rate)
        temp_coh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / B.shape[0]
        ifg_num = B.shape[0]
    except:
        pass

    return ts, temp_coh, ifg_num


def network_inversion_wls(A, ifgram, weight, skip_zero_phase=True, Astd=None):
    """Network inversion based on Weighted Least Square (WLS) solution.
    Inputs:
        A      - 2D np.array in size of (ifgram_num, num_date-1)
                 representing date configuration for each interferogram
                 (-1 for master, 1 for slave, 0 for others)
        ifgram - np.array in size of (ifgram_num,) or (ifgram_num, 1)
                 phase of all interferograms
        weight - np.array in size of (ifgram_num,) or (ifgram_num, 1)
                 weight of ifgram
        skip_zero_phase - bool, skip ifgram with zero phase value
        Astd   - 2D np.array in size of (ifgram_num, num_date-1)
                 design matrix for STD calculation excluding the reference date
    Output:
        ts      - 1D np.array in size of (num_date-1,), phase time series
        temp_coh - float32, temporal coherence
        ts_std   - 1D np.array in size of (num_date-1,), decor noise std time series
    """
    ifgram = np.reshape(ifgram, (-1, 1))
    if Astd is None:
        Astd = A

    dateNum1 = A.shape[1]
    ts = np.zeros(dateNum1, np.float32)
    ts_std = np.zeros(dateNum1, np.float32)
    temp_coh = 0.
    ifg_num = 0

    # Skip Zero Phase Value
    if skip_zero_phase and not np.all(ifgram):
        idx = (ifgram != 0.).flatten()
        A = A[idx, :]
        if A.shape[0] < dateNum1*2:
            return ts, temp_coh, ts_std, ifg_num
        ifgram = ifgram[idx]
        weight = weight[idx]
        Astd = Astd[idx, :]

    W = np.diag(weight.flatten())
    try:
        # WLS Inversion
        ATW = np.array(A.T.dot(W), np.float32)
        ts = np.array(np.linalg.inv(ATW.dot(A)).dot(ATW), np.float32).dot(ifgram)
        # A_inv_wls = np.linalg.inv(A.T.dot(W).dot(A))
        # ts = A_inv_wls.dot(A.T).dot(W).dot(ifgram.reshape(-1,1))

        # Temporal Coherence
        ifgram_diff = ifgram - np.dot(A, ts)
        temp_coh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / A.shape[0]
        ifg_num = A.shape[0]

        # Decorrelation Noise Std
        # ts_std = np.sqrt(np.diag(np.linalg.inv(Astd.T.dot(W).dot(Astd))))
    except:
        pass

    return ts, temp_coh, ts_std, ifg_num


def temporal_coherence(A, ts, ifgram, weight=None, chunk_size=500):
    """Calculate temporal coherence based on Tizzani et al. (2007, RSE)
    Inputs:
        A      - 2D np.array in size of (ifgram_num, num_date-1)
                 representing date configuration for each interferogram
                 (-1 for master, 1 for slave, 0 for others)
        ts     - 2D np.array in size of (num_date-1, num_pixel), phase time series
        ifgram - 2D np.array in size of (ifgram_num, num_pixel), observed interferometric phase
        weight - 2D np.array in size of (ifgram_num, num_pixel), weight of ifgram
        chunk_size - int, max number of pixels per loop during the calculation
    Output:
        temp_coh - 1D np.array in size of (num_pixel), temporal coherence
    """
    # Default: uniform weight
    if weight is None:
        weight = np.ones(ifgram.shape, np.float32)

    # Calculate weighted temporal coherence
    if ifgram.ndim == 1 or ifgram.shape[1] <= chunk_size:
        ifgram_diff = ifgram - np.dot(A, ts)
        temp_coh = np.abs(np.sum(np.multiply(weight, np.exp(1j*ifgram_diff)),
                                 axis=0)) / np.sum(weight, axis=0)

    else:
        # Loop chunk by chunk to reduce memory usage
        num_pixel = ifgram.shape[1]
        temp_coh = np.zeros(num_pixel, np.float32)

        chunk_num = int((num_pixel-1)/chunk_size) + 1
        for i in range(chunk_num):
            sys.stdout.write('\rcalculating chunk %s/%s ...' % (i+1, chunk_num))
            sys.stdout.flush()
            p0 = i*chunk_size
            p1 = min([p0+chunk_size, num_pixel])
            ifgram_diff = ifgram[:, p0:p1] - np.dot(A, ts[:, p0:p1])
            temp_coh[p0:p1] = np.abs(np.sum(np.exp(1j*ifgram_diff),
                                            axis=0)) / np.sum(weight[:, p0:p1], axis=0)
        print('')
    return temp_coh


###########################################################################################
def write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std=None, num_inv_ifgram=None,
                    suffix=''):
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    date_list = stack_obj.get_date_list(dropIfgram=True)

    # File 1 - timeseries.h5
    ts_file = 'timeseries{}.h5'.format(suffix)
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
    out_file = 'temporalCoherence{}.h5'.format(suffix)
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
    writefile.write(num_inv_ifgram, out_file=out_file, metadata=metadata)
    return


def split_ifgram_file(ifgram_file, chunk_size=100e6):
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    metadata = dict(stack_obj.metadata)

    # get reference phase
    ref_phase = get_ifgram_reference_phase(ifgram_file, drop_ifgram=False)

    # get list of boxes
    box_list = split_into_boxes(ifgram_file,
                                chunk_size=chunk_size,
                                print_msg=True)
    num_box = len(box_list)

    # read/write each patch file
    outfile_list = []
    for i in range(num_box):
        box = box_list[i]
        outfile = '{}_{:03d}{}'.format(os.path.splitext(ifgram_file)[0],
                                       i+1,
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


def split_into_boxes(ifgram_file, chunk_size=100e6, print_msg=True):
    """Split into chunks in rows to reduce memory usage
    Parameters:
    """
    shape = ifgramStack(ifgram_file).get_size()
    # Get r_step / chunk_num
    r_step = chunk_size / (shape[0] * shape[2])         # split in lines
    #if inps.weightFunc != 'sbas':  # more memory usage (coherence) for WLS
    #    r_step /= 2.0
    r_step = int(ceil_to_1(r_step))
    chunk_num = int((shape[1]-1)/r_step) + 1

    if print_msg and chunk_num > 1:
        print('maximum chunk size: %.1E' % chunk_size)
        print('split %d lines into %d patches for processing' % (shape[1], chunk_num))
        print('    with each patch up to %d lines' % r_step)

    # Computing the inversion
    box_list = []
    for i in range(chunk_num):
        r0 = i * r_step
        r1 = min([shape[1], r0+r_step])
        box = (0, r0, shape[2], r1)
        box_list.append(box)
    return box_list


def check_design_matrix(ifgram_file, weight_func='fim'):
    """Check Rank of Design matrix for weighted inversion"""
    A = ifgramStack(ifgram_file).get_design_matrix(dropIfgram=True)[0]
    print('-------------------------------------------------------------------------------')
    if weight_func == 'sbas':
        print('generic least square inversion with min-norm phase velocity')
        print('    based on Berardino et al. (2002, IEEE-TGRS)')
        print('    OLS for pixels with full rank      network')
        print('    SVD for pixels with rank deficient network')
        if np.linalg.matrix_rank(A) < A.shape[1]:
            print('WARNING: singular design matrix! Inversion result can be biased!')
            print('continue using its SVD solution on all pixels')
    else:
        print('weighted least square (WLS) inversion with min-norm phase, pixelwise')
        if np.linalg.matrix_rank(A) < A.shape[1]:
            print('ERROR: singular design matrix!')
            print('    Input network of interferograms is not fully connected!')
            print('    Can not invert the weighted least square solution.')
            print('You could try:')
            print('    1) Add more interferograms to make the network fully connected:')
            print('       a.k.a., no multiple subsets nor network islands')
            print("    2) Use '-w no' option for non-weighted SVD solution.")
            raise Exception()
    print('-------------------------------------------------------------------------------')
    return A


def get_ifgram_reference_phase(ifgram_file, skip_reference=False, drop_ifgram=True):
    """Read refPhase"""
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.get_size()
    stack_obj.get_metadata()
    try:
        ref_y = int(stack_obj.metadata['REF_Y'])
        ref_x = int(stack_obj.metadata['REF_X'])
        ref_phase = np.squeeze(stack_obj.read(datasetName='unwrapPhase',
                                              box=(ref_x, ref_y, ref_x+1, ref_y+1),
                                              dropIfgram=drop_ifgram,
                                              print_msg=False))
        print('reference pixel in y/x: {}'.format((ref_y, ref_x)))
    except:
        if skip_reference:
            ref_phase = np.zeros((stack_obj.numIfgram,), np.float32)
            print('skip checking reference pixel info - This is for SIMULATION ONLY.')
        else:
            msg = 'ERROR: No REF_X/Y found! Can not invert interferograms without reference in space.'
            msg += '\nrun reference_point.py {} for a quick referencing.'.format(stack_obj.file)
            raise Exception(msg)
    return ref_phase


def read_unwrap_phase(stack_obj, box, ref_phase, skip_zero_phase=True):
    """Read unwrapPhase from ifgramStack file
    Parameters: stack_obj : ifgramStack object
                box : tuple of 4 int
                ref_phase : 1D array or None
                skip_zero_phase : bool
    Returns:    pha_data : 3D array of unwrapPhase
    """
    # Read unwrapPhase
    num_ifgram = np.sum(stack_obj.dropIfgram)
    print('reading unwrapPhase in {} * {} ...'.format(box, num_ifgram))
    pha_data = stack_obj.read(datasetName='unwrapPhase',
                              box=box,
                              dropIfgram=True,
                              print_msg=False).reshape(num_ifgram, -1)

    # read ref_phase
    if ref_phase is not None:
        # use input ref_phase array (for split_file=False)
        print('use input reference phase')
    elif 'refPhase' in stack_obj.datasetNames:
        # read ref_phase from file itself (for split_file=True)
        print('read reference phase from file')
        with h5py.File(stack_obj.file, 'r') as f:
            ref_phase = f['refPhase'][:]
    else:
        raise Exception('No reference phase input/found on file!'+
                        ' unwrapped phase is not referenced!')

    # determine msk_value
    if skip_zero_phase:
        # print('skip zero value in unwrapped phase')
        msk_value = 0.
    else:
        msk_value = np.nan

    # reference unwrapPhase
    for i in range(num_ifgram):
        mask = pha_data[i, :] != msk_value
        pha_data[i, :][mask] -= ref_phase[i]
    return pha_data


def mask_unwrap_phase(pha_data, stack_obj, box, mask_ds_name=None, mask_threshold=0.4):
    # Read/Generate Mask
    num_ifgram = np.sum(stack_obj.dropIfgram)
    if mask_ds_name and mask_ds_name in stack_obj.datasetNames:
        print('reading {} in {} * {} ...'.format(mask_ds_name, box, num_ifgram))
        msk_data = stack_obj.read(datasetName=mask_ds_name,
                                  box=box,
                                  dropIfgram=True,
                                  print_msg=False).reshape(num_ifgram, -1)
        if mask_ds_name == 'coherence':
            msk_data = msk_data >= mask_threshold
            print('mask out pixels with {} < {}'.format(mask_ds_name, mask_threshold))
        else:
            print('mask out pixels with {} == 0'.format(mask_ds_name))
        pha_data[msk_data == 0.] = 0.
        del msk_data
    return pha_data


def read_coherence2weight(stack_obj, box, weight_func='fim'):
    epsilon = 1e-4
    num_ifgram = np.sum(stack_obj.dropIfgram)
    print('reading coherence in {} * {} ...'.format(box, num_ifgram))
    coh_data = stack_obj.read(datasetName='coherence',
                              box=box,
                              dropIfgram=True,
                              print_msg=False).reshape(num_ifgram, -1)
    coh_data[np.isnan(coh_data)] = epsilon

    # Calculate Weight matrix
    weight = np.array(coh_data, np.float64)
    del coh_data
    L = int(stack_obj.metadata['ALOOKS']) * int(stack_obj.metadata['RLOOKS'])
    if weight_func == 'var':
        print('convert coherence to weight using inverse of phase variance')
        print('    with phase PDF for distributed scatterers from Tough et al. (1995)')
        weight = 1.0 / coherence2phase_variance_ds(weight, L, print_msg=True)

    elif weight_func == 'coh':
        print('use coherence as weight directly (Perissin & Wang, 2012; Tong et al., 2016)')
        weight[weight < epsilon] = epsilon

    elif weight_func == 'fim':
        print('convert coherence to weight using Fisher Information Index (Seymour & Cumming, 1994)')
        weight = coherence2fisher_info_index(weight, L)

    else:
        raise Exception('Un-recognized weight function: %s' % weight_func)

    weight = np.array(weight, np.float32)
    return weight


def ifgram_inversion_patch(ifgram_file, box=None, ref_phase=None, weight_func='fim',
                           mask_dataset_name=None, mask_threshold=0.4,
                           water_mask_file=None, skip_zero_phase=True):
    """Invert one patch of an ifgram stack into timeseries.
    Parameters: ifgram_file       : str, interferograms stack HDF5 file, e.g. ./INPUTS/ifgramStack.h5
                box               : tuple of 4 int, indicating (x0, y0, x1, y1) pixel coordinate of area of interest
                                    or None, to process the whole file and write output file
                ref_phase         : 1D array in size of (num_ifgram) 
                                    or None
                weight_func       : str, weight function, choose in ['sbas', 'fim', 'var', 'coh']
                mask_dataset_name : str, dataset name in ifgram_file used to mask unwrapPhase pixelwisely
                mask_threshold    : float, min coherence of pixels if mask_dataset_name='coherence'
                water_mask_file   : str, water mask filename if available,
                                    skip inversion on water to speed up the process
                skip_zero_phase   : bool, skip zero value of unwrapped phase or not, default yes, for comparison
    Returns:    ts             : 3D array in size of (num_date, num_row, num_col)
                temp_coh       : 2D array in size of (num_row, num_col)
                ts_std         : 3D array in size of (num_date, num_row, num_col)
                num_inv_ifgram : 2D array in size of (num_row, num_col)
    Example:    ifgram_inversion_patch('ifgramStack.h5', box=(0,200,1316,400), ref_phase=np.array(),
                                       weight_func='fim', mask_dataset_name='coherence')
                ifgram_inversion_patch('ifgramStack_001.h5', box=None, ref_phase=None,
                                       weight_func='fim', mask_dataset_name='coherence')
    """

    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)

    # Size Info - Patch
    if box:
        print('processing %8d/%d lines ...' % (box[3], stack_obj.length))
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
    A, B = stack_obj.get_design_matrix(date12_list=date12_list)
    num_ifgram = len(date12_list)
    try:
        ref_date = str(np.loadtxt('reference_date.txt', dtype=bytes).astype(str))
    except:
        ref_date = date_list[0]
    ref_idx = date_list.index(ref_date)
    time_idx = [i for i in range(num_date)]
    time_idx.remove(ref_idx)
    Astd = stack_obj.get_design_matrix(refDate=ref_date, dropIfgram=True)[0]

    # Initialization of output matrix
    print('number of interferograms: {}'.format(num_ifgram))
    print('number of acquisitions  : {}'.format(num_date))
    print('number of lines  : {}'.format(stack_obj.length))
    print('number of columns: {}'.format(stack_obj.width))
    ts = np.zeros((num_date, num_pixel), np.float32)
    ts_std = np.zeros((num_date, num_pixel), np.float32)
    temp_coh = np.zeros(num_pixel, np.float32)
    num_inv_ifgram = np.zeros(num_pixel, np.int16)

    # Read/Mask unwrapPhase
    pha_data = read_unwrap_phase(stack_obj,
                                 box,
                                 ref_phase,
                                 skip_zero_phase=skip_zero_phase)

    pha_data = mask_unwrap_phase(pha_data,
                                 stack_obj,
                                 box,
                                 mask_ds_name=mask_dataset_name,
                                 mask_threshold=mask_threshold)

    # Mask for pixels to invert
    mask = np.ones(num_pixel, np.bool_)
    # 1 - Water Mask
    if water_mask_file:
        print(('skip pixels on water with mask from'
               ' file: {}').format(os.path.basename(water_mask_file)))
        dsNames = readfile.get_dataset_list(water_mask_file)
        dsName = [i for i in dsNames
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
           ' ({:.1f}%)').format(num_pixel2inv,
                                num_pixel,
                                num_pixel2inv/num_pixel*100))
    if num_pixel2inv < 1:
        ts = ts.reshape(num_date, num_row, num_col)
        temp_coh = temp_coh.reshape(num_row, num_col)
        ts_std = ts_std.reshape(num_date, num_row, num_col)
        num_inv_ifgram = num_inv_ifgram.reshape(num_row, num_col)
        return ts, temp_coh, ts_std, num_inv_ifgram

    # Inversion - SBAS
    if weight_func == 'sbas':
        # Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        mask_all_net = np.all(pha_data, axis=0)
        mask_all_net *= mask
        mask_part_net = mask ^ mask_all_net

        if np.sum(mask_all_net) > 0:
            print(('inverting pixels with valid phase in all  ifgrams'
                   ' ({:.0f} pixels) ...').format(np.sum(mask_all_net)))
            # num_all_net = int(np.sum(mask_all_net))
            # pha_data_temp = pha_data[:, mask_all_net]
            # ts1 = np.zeros((num_date-1, num_all_net))
            # temp_coh1 = np.zeros(num_all_net)
            # step = 1000
            # loop_num = int(np.floor(num_all_net/step))
            # prog_bar = ptime.progressBar(maxValue=loop_num)
            # for i in range(loop_num):
            #     [i0, i1] = [i * step, min((i + 1) * step, num_all_net)]
            #     ts1[:, i0:i1], temp_coh1[i0:i1] = network_inversion_sbas(B,
            #                                                              ifgram=pha_data_temp[:, i0:i1],
            #                                                              tbase_diff=tbase_diff,
            #                                                              skip_zero_phase=False)
            #     prog_bar.update(i+1, suffix=i0)
            # prog_bar.close()
            ts1, temp_coh1, ifg_num1 = network_inversion_sbas(B, ifgram=pha_data[:, mask_all_net], 
                                                              tbase_diff=tbase_diff,
                                                              skip_zero_phase=False)
            ts[1:, mask_all_net] = ts1
            temp_coh[mask_all_net] = temp_coh1
            num_inv_ifgram[mask_all_net] = ifg_num1

        if np.sum(mask_part_net) > 0:
            print(('inverting pixels with valid phase in some ifgrams'
                   ' ({:.0f} pixels) ...').format(np.sum(mask_part_net)))
            num_pixel2inv = int(np.sum(mask_part_net))
            idx_pixel2inv = np.where(mask_part_net)[0]
            prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
            for i in range(num_pixel2inv):
                idx = idx_pixel2inv[i]
                ts1, temp_coh1, ifg_num1 = network_inversion_sbas(B, ifgram=pha_data[:, idx],
                                                                  tbase_diff=tbase_diff,
                                                                  skip_zero_phase=skip_zero_phase)
                ts[1:, idx] = ts1.flatten()
                temp_coh[idx] = temp_coh1
                num_inv_ifgram[idx] = ifg_num1
                prog_bar.update(i+1, every=100, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
            prog_bar.close()

    # Inversion - WLS
    else:
        weight = read_coherence2weight(stack_obj, box=box, weight_func=weight_func)

        # Converting to 32 bit floats leads to 2X speedup
        # (comment it out as we now convert it beforehand)
        # A = np.array(A, np.float32)
        # pha_data = np.array(pha_data, np.float32)
        # weight = np.array(weight, np.float32)
        # Astd = np.array(Astd, np.float32)

        # Weighted Inversion pixel by pixel
        print('inverting network of interferograms into time series ...')
        prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
        for i in range(num_pixel2inv):
            idx = idx_pixel2inv[i]
            ts1, temp_coh1, ts_std1, ifg_numi = network_inversion_wls(A, ifgram=pha_data[:, idx],
                                                                      weight=weight[:, idx],
                                                                      skip_zero_phase=skip_zero_phase,
                                                                      Astd=Astd)
            ts[1:, idx] = ts1.flatten()
            temp_coh[idx] = temp_coh1
            ts_std[time_idx, idx] = ts_std1.flatten()
            num_inv_ifgram[idx] = ifg_numi
            prog_bar.update(i+1, every=100, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
        prog_bar.close()

    ts = ts.reshape(num_date, num_row, num_col)
    ts_std = ts_std.reshape(num_date, num_row, num_col)
    temp_coh = temp_coh.reshape(num_row, num_col)
    num_inv_ifgram = num_inv_ifgram.reshape(num_row, num_col)

    # write output files if input file is splitted (box == None)
    if box is None:
        # metadata
        metadata = dict(stack_obj.metadata)
        metadata[key_prefix+'weightFunc'] = weight_func
        suffix = re.findall('_\d{3}', ifgram_file)[0]
        write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std, num_inv_ifgram, suffix)
        return
    else:
        return ts, temp_coh, ts_std, num_inv_ifgram


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

    if inps.update_mode and not ut.update_file(inps.timeseriesFile, ifgram_file):
        return inps.timeseriesFile, inps.tempCohFile

    A = check_design_matrix(ifgram_file, weight_func=inps.weightFunc)
    num_date = A.shape[1] + 1

    # split ifgram_file into blocks to save memory
    box_list = split_into_boxes(ifgram_file, chunk_size=inps.chunk_size)
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
                                   weight_func=inps.weightFunc,
                                   mask_dataset_name=inps.maskDataset,
                                   mask_threshold=inps.maskThreshold,
                                   water_mask_file=inps.waterMaskFile,
                                   skip_zero_phase=inps.skip_zero_phase)
    else:
        # read ifgram_file in small patches and write them together
        ref_phase = get_ifgram_reference_phase(ifgram_file,
                                               skip_reference=inps.skip_ref,
                                               drop_ifgram=True)

        # Initialization of output matrix
        stack_obj = ifgramStack(ifgram_file)
        stack_obj.open(print_msg=False)
        ts = np.zeros((num_date, stack_obj.length, stack_obj.width), np.float32)
        temp_coh = np.zeros((stack_obj.length, stack_obj.width), np.float32)
        ts_std = np.zeros(ts.shape, np.float32)
        num_inv_ifgram = np.zeros(temp_coh.shape, np.int16)

        # Loop
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
                                                weight_func=inps.weightFunc,
                                                mask_dataset_name=inps.maskDataset,
                                                mask_threshold=inps.maskThreshold,
                                                water_mask_file=inps.waterMaskFile,
                                                skip_zero_phase=inps.skip_zero_phase)

            temp_coh[box[1]:box[3], box[0]:box[2]] = temp_cohi
            ts[:, box[1]:box[3], box[0]:box[2]] = tsi
            ts_std[:, box[1]:box[3], box[0]:box[2]] = ts_stdi
            num_inv_ifgram[box[1]:box[3], box[0]:box[2]] = ifg_numi

        # metadata
        metadata = dict(stack_obj.metadata)
        metadata[key_prefix+'weightFunc'] = inps.weightFunc
        write2hdf5_file(ifgram_file, metadata, ts, temp_coh, ts_std, num_inv_ifgram, suffix='')

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.templateFile:
        inps = read_template2inps(inps.templateFile, inps)
    inps.timeseriesFile, inps.tempCohFile = inps.outfile

    # Input file info
    atr = readfile.read_attribute(inps.ifgramStackFile)
    if atr['FILE_TYPE'] != 'ifgramStack':
        print('ERROR: only ifgramStack file supported, input is {} file!'.format(atr['FILE_TYPE']))
        sys.exit(1)

    # Network Inversion
    if inps.residualNorm == 'L2':
        print('inverse time-series using L2 norm minimization')
        ifgram_inversion(inps.ifgramStackFile, inps)
    else:
        print('inverse time-series using L1 norm minimization')
        ut.timeseries_inversion_L1(inps.ifgramStackFile, inps.timeseriesFile)
    return


################################################################################################
if __name__ == '__main__':
    main()
