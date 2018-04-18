#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys
import time
import argparse
import string
import h5py
import numpy as np
from scipy.special import gamma
from pysar.utils import readfile, writefile, datetime as ptime, utils as ut
from pysar.objects import ifgramStack, timeseries


################################################################################################
EXAMPLE='''example:
  ifgram_inversion.py  ifgramStack.h5
  ifgram_inversion.py  ifgramStack.h5 -t pysarApp_template.txt
  ifgram_inversion.py  ifgramStack.h5 -w var
  ifgram_inversion.py  ifgramStack.h5 -w fim
  ifgram_inversion.py  ifgramStack.h5 -w coh
'''

TEMPLATE='''
## Invert network of interferograms into time series using weighted least sqaure (WLS) estimator.
## There are 4 weighting options:
## 1) fim  - WLS, use Fisher Information Matrix as weight (Seymour & Cumming, 1994, IGARSS). [Recommended]
## 2) var  - WLS, use inverse of covariance as weight (Guarnieri & Tebaldini, 2008, TGRS)
## 3) coh  - WLS, use coherence as weight (Perissin & Wang, 2012, IEEE-TGRS)
## 4) sbas - LS/SVD, uniform weight (Berardino et al., 2002, TGRS)
## There are 3 mask options to mask unwrapPhase for each interferogram:
## 1) connectComponent - mask out pixels with False/0 value [Recommended]
## 2) coherence        - mask out pixels with value < maskThreshold
## 3) no               - no masking.
## Temporal coherence is calculated and used to generate final mask (Pepe & Lanari, 2006, IEEE-TGRS)
pysar.networkInversion.weightFunc    = auto #[fim / var / coh / sbas], auto for sbas
pysar.networkInversion.maskDataset   = auto #[connectComponent / coherence / no], auto for connectComponent
pysar.networkInversion.maskThreshold = auto #[0-1], auto for 0.4
pysar.networkInversion.waterMaskFile = auto #[filename / no], auto for no
pysar.networkInversion.residualNorm  = auto #[L2 ], auto for L2, norm minimization solution
pysar.networkInversion.minTempCoh    = auto #[0.0-1.0], auto for 0.7, min temporal coherence for mask
'''

REFERENCE='''references:
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
'''

def createParser():
    parser = argparse.ArgumentParser(description='Invert network of interferograms into timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgramStackFile', help='interferograms stack file to be inverted')
    parser.add_argument('--template','-t', dest='templateFile',\
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--maskDataset', dest='maskDataset', default='connComp',\
                        help='dataset used to mask unwrapPhase')
    parser.add_argument('--maskThreshold', type=float, default=0.4,\
                        help='threshold to generate mask when mask is coherence')

    parser.add_argument('--weight-function','-w', dest='weightFunc', default='no',\
                        help='function used to convert coherence to weight for inversion:\n'+\
                             'fim  - Fisher Information Matrix as weight'+\
                             'var  - phase variance due to temporal decorrelation\n'+\
                             'coh  - uniform distribution CDF function\n'+\
                             'sbas - uniform weight')
    parser.add_argument('--norm', dest='residualNorm', default='L2', choices=['L1','L2'],\
                        help='Inverse method used to residual optimization, L1 or L2 norm minimization. Default: L2')

    parser.add_argument('--chunk-size', dest='chunk_size', type=float, default=0.5e9,\
                        help='max number of data (= ifgram_num * row_num * col_num) to read per loop\n'+\
                             'default: 0.5G; adjust it according to your computer memory.')
    parser.add_argument('--parallel', dest='parallel', action='store_true',\
                        help='Enable parallel processing for the pixelwise weighted inversion. [not working yet]')
    parser.add_argument('--skip-reference', dest='skip_ref', action='store_true',\
                        help='Skip checking reference pixel value, for simulation testing.')
    parser.add_argument('-o','--output', dest='outfile', nargs=2, default=['timeseries.h5','temporalCoherence.h5'],\
                        help='Output file name for timeseries and temporal coherence, default:\n'+\
                             'timeseries.h5 temporalCoherence.h5')
    parser.add_argument('--update-mode', dest='update_mode', action='store_true',\
                        help='Enable update mode, and skip inversion if output timeseries file already exists,\n'+\
                             'readable and newer than input interferograms file')
    parser.add_argument('--noskip-zero-phase', dest='skip_zero_phase', action='store_false',\
                        help='Do not skip interferograms with zero phase.')
    parser.add_argument('--water-mask','-m', dest='waterMaskFile',\
                        help='Skip inversion on the masked out region, i.e. water.')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    inps.parallel = False
    return inps


def read_template2inps(template_file, inps):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()
    inpsDict = vars(inps)
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    prefix = 'pysar.networkInversion.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if key in ['maskDataset']:
            inpsDict[key] = value
        elif value:
            if key in ['maskThreshold']:
                inpsDict[key] = float(value)
            elif key in ['weightFunc','residualNorm','waterMaskFile']:
                inpsDict[key] = value
    return inps


################################################################################################
def phase_pdf_ds(L, coherence=None, phiNum=1000):
    '''Marginal PDF of interferometric phase for distributed scatterers (DS)
    Eq. 66 (Tough et al., 1995) and Eq. 4.2.23 (Hanssen, 2001)
    Inputs:
        L         - int, number of independent looks
        coherence - 1D np.array for the range of coherence, with value < 1.0 for valid operation
        phiNum    - int, number of phase sample for the numerical calculation
    Output:
        pdf       - 2D np.array, phase pdf in size of (phiNum, len(coherence))
        coherence - 1D np.array for the range of coherence
    Example:
        epsilon = 1e-4
        coh = np.linspace(0., 1-epsilon, 1000)
        pdf, coh = phase_pdf_ds(1, coherence=coh)
    '''
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000)
    coherence = np.array(coherence, np.float64).reshape(1,-1)
    phi = np.linspace(-np.pi, np.pi, phiNum, dtype=np.float64).reshape(-1,1)

    ### Phase PDF - Eq. 4.2.32 (Hanssen, 2001)
    A = np.power((1-np.square(coherence)), L) / (2*np.pi)
    A = np.tile(A, (phiNum, 1))
    B = gamma(2*L - 1) / ((gamma(L))**2 * 2**(2*(L-1)))

    beta = np.multiply(np.abs(coherence), np.cos(phi), dtype=np.float64)
    #C1 = np.power((1 - np.square(beta)), L+0.5)
    #C1[C1 == 0.] = epsilon
    #C = np.divide((2*L - 1) * beta, C1)
    C = np.divide((2*L - 1) * beta, np.power((1 - np.square(beta)), L+0.5))
    C = np.multiply(C, (np.pi/2 + np.arcsin(beta)))
    #C2 = np.power((1 - np.square(beta)), L)
    #C2[C2 == 0.0] = epsilon
    #C += 1 / C2
    C += 1 / np.power((1 - np.square(beta)), L)

    sumD = 0
    if L > 1:
        for r in range(L-1):
            D = gamma(L-0.5) / gamma(L-0.5-r)
            D *= gamma(L-1-r) / gamma(L-1)
            #D1 = np.power((1 - np.square(beta)), r+2)
            #D1[D1 == 0.] = epsilon
            #D *= (1 + (2*r+1)*np.square(beta)) / D1
            D *= (1 + (2*r+1)*np.square(beta)) / np.power((1 - np.square(beta)), r+2)
            sumD += D
        sumD /= (2*(L-1))

    pdf = B*C + sumD
    pdf = np.multiply(A, pdf)
    return pdf, coherence.flatten()


def phase_variance_ds(L,  coherence=None):
    '''Interferometric phase variance for distributed scatterers (DS)
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
    '''
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000, np.float64)
    phiNum = len(coherence)

    phi = np.linspace(-np.pi, np.pi, phiNum, np.float64).reshape(-1,1)
    phi_step = 2*np.pi/phiNum

    pdf, coherence = phase_pdf_ds(L, coherence=coherence)
    var = np.sum(np.multiply(np.square(np.tile(phi, (1, len(coherence)))), pdf)*phi_step, axis=0)
    return var, coherence


def phase_variance_ps(L, coherence=None):
    '''the Cramer-Rao bound (CRB) of phase variance
    Given by Eq. 25 (Rodriguez and Martin, 1992)and Eq 4.2.32 (Hanssen, 2001)
    Valid when coherence is close to 1.
    '''
    epsilon = 1e-4
    if coherence is None:
        coherence = np.linspace(0.9, 1.-epsilon, 1000, np.float64)
    var = (1-coherence**2) / (2*L*coherence**2)
    return var, coherence


def coherence2phase_variance_ds(coherence, L=32, print_msg=False):
    '''Convert coherence to phase variance based on DS phase PDF (Tough et al., 1995)'''
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

def coherence2fisher_info_index(coherence, L=32, epsilon=1e-4):
    '''Convert coherence to Fisher information index (Seymour & Cumming, 1994, IGARSS)'''
    coherence = np.array(coherence, np.float64)
    coherence[coherence > 1-epsilon] = 1-epsilon
    weight = 2.0 * L * np.square(coherence) / (1 - np.square(coherence))
    return weight


def round_to_1(x):
    '''Return the most significant digit of input number'''
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)

def ceil_to_1(x):
    '''Return the most significant digit of input number and ceiling it'''
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)+10**digit


def network_inversion_sbas(B, ifgram, tbase_diff, skipZeroPhase=True):
    ''' Network inversion based on Small BAseline Subsets (SBAS) algorithm (Berardino et al.,
        2002, IEEE-TGRS). For full rank design matrix, a.k.a., fully connected network, ordinary
        least square (OLS) inversion is applied; otherwise, Singular Value Decomposition (SVD).

    Inputs:
        B          - 2D np.array in size of (ifgram_num, date_num-1)
                     design matrix B, which represents temporal baseline timeseries between
                     master and slave date for each interferogram
        ifgram     - 2D np.array in size of (ifgram_num, pixel_num)
                     phase of all interferograms
        tbase_diff - 2D np.array in size of (date_num-1, 1)
                     differential temporal baseline of time-series
        skipZeroPhase - bool, skip ifgram with zero phase value
    Output:
        ts      - 2D np.array in size of (date_num-1, pixel_num), phase time series
        tempCoh - 1D np.array in size of (pixel_num), temporal coherence
    '''
    ifgram = ifgram.reshape(B.shape[0],-1)
    dateNum1 = B.shape[1]
    ts = np.zeros(dateNum1, np.float32)
    tempCoh = 0.

    ## Skip Zero Phase Value
    if skipZeroPhase and not np.all(ifgram):
        idx = (ifgram != 0.).flatten()
        B = B[idx, :]
        if B.shape[0] < dateNum1*2:
            return ts, tempCoh
        ifgram = ifgram[idx, :]

    try:
        ## Invert time-series
        B_inv = np.array(np.linalg.pinv(B), np.float32)
        ts_rate = np.dot(B_inv, ifgram)
        ts_diff = ts_rate * np.tile(tbase_diff, (1, ifgram.shape[1]))
        ts = np.cumsum(ts_diff, axis=0)

        ## Temporal Coherence
        ifgram_diff = ifgram - np.dot(B, ts_rate)
        tempCoh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / B.shape[0]
    except:
        pass

    return ts, tempCoh


def network_inversion_wls(A, ifgram, weight, skipZeroPhase=True, Astd=None):
    '''Network inversion based on Weighted Least Square (WLS) solution.
    Inputs:
        A      - 2D np.array in size of (ifgram_num, date_num-1)
                 representing date configuration for each interferogram
                 (-1 for master, 1 for slave, 0 for others)
        ifgram - np.array in size of (ifgram_num,) or (ifgram_num, 1)
                 phase of all interferograms
        weight - np.array in size of (ifgram_num,) or (ifgram_num, 1)
                 weight of ifgram
        skipZeroPhase - bool, skip ifgram with zero phase value
        Astd   - 2D np.array in size of (ifgram_num, date_num-1)
                 design matrix for STD calculation excluding the reference date
    Output:
        ts      - 1D np.array in size of (date_num-1,), phase time series
        tempCoh - float32, temporal coherence
        tsStd   - 1D np.array in size of (date_num-1,), decor noise std time series
    '''
    if Astd is None:
        Astd = A

    dateNum1 = A.shape[1]
    ts = np.zeros(dateNum1, np.float32)
    tsStd = np.zeros(dateNum1, np.float32)
    tempCoh = 0.
    ifgram = np.reshape(ifgram, (-1,1))

    ## Skip Zero Phase Value
    if skipZeroPhase and not np.all(ifgram):
        idx = ifgram != 0.
        A = A[idx,:]
        if A.shape[0] < dateNum1:
            return ts, tempCoh, tsStd
        ifgram = ifgram[idx]
        weight = weight[idx]
        Astd = Astd[idx,:]

    W = np.diag(weight.flatten())
    try:
        ## WLS Inversion
        ATW = np.array(A.T.dot(W), np.float32)
        ts = np.array(np.linalg.inv(ATW.dot(A)).dot(ATW), np.float32).dot(ifgram)
        #A_inv_wls = np.linalg.inv(A.T.dot(W).dot(A))
        #ts = A_inv_wls.dot(A.T).dot(W).dot(ifgram.reshape(-1,1))

        ## Temporal Coherence
        ifgram_diff = ifgram - np.dot(A, ts)
        tempCoh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / A.shape[0]

        ## Decorrelation Noise Std
        #tsStd = np.sqrt(np.diag(np.linalg.inv(Astd.T.dot(W).dot(Astd))))
        connNet = True
    except:
        connNet = False

    return ts, tempCoh, tsStd


def temporal_coherence(A, ts, ifgram, weight=None, chunk_size=500):
    '''Calculate temporal coherence based on Tizzani et al. (2007, RSE)
    Inputs:
        A      - 2D np.array in size of (ifgram_num, date_num-1)
                 representing date configuration for each interferogram
                 (-1 for master, 1 for slave, 0 for others)
        ts     - 2D np.array in size of (date_num-1, pixel_num), phase time series
        ifgram - 2D np.array in size of (ifgram_num, pixel_num), observed interferometric phase
        weight - 2D np.array in size of (ifgram_num, pixel_num), weight of ifgram
        chunk_size - int, max number of pixels per loop during the calculation
    Output:
        temp_coh - 1D np.array in size of (pixel_num), temporal coherence
    '''
    # Default: uniform weight
    if weight is None:
        weight = np.ones(ifgram.shape, np.float32)

    # Calculate weighted temporal coherence
    if ifgram.ndim == 1 or ifgram.shape[1] <= chunk_size:
        ifgram_diff = ifgram - np.dot(A, ts)
        temp_coh = np.abs(np.sum(np.multiply(weight, np.exp(1j*ifgram_diff)), axis=0)) / np.sum(weight, axis=0)

    else:
        #Loop chunk by chunk to reduce memory usage
        pixel_num = ifgram.shape[1]
        temp_coh = np.zeros(pixel_num, np.float32)

        chunk_num = int((pixel_num-1)/chunk_size) + 1
        for i in range(chunk_num):
            sys.stdout.write('\rcalculating chunk %s/%s ...' % (i+1, chunk_num))
            sys.stdout.flush()
            p0 = i*chunk_size
            p1 = min([p0+chunk_size, pixel_num])
            ifgram_diff = ifgram[:,p0:p1] - np.dot(A, ts[:,p0:p1])
            temp_coh[p0:p1] = np.abs(np.sum(np.multiply(weight[:,p0:p1], np.exp(1j*ifgram_diff)), axis=0)) /\
                              np.sum(weight[:,p0:p1], axis=0)
        print('')
    return temp_coh


def read_unwrap_phase(stackobj, inps, box):
    ##### Read unwrapPhase
    print('reading unwrapPhase in {} * {} ...'.format(box, inps.numIfgram))
    pha_data = None
    pha_data = stackobj.read(datasetName='unwrapPhase', box=box, dropIfgram=True).reshape(inps.numIfgram,-1)
    if inps.skip_zero_phase:
        #print('skip zero phase value (masked out and filled during phase unwrapping)')
        for i in range(inps.numIfgram):
            pha_data[i,:][pha_data[i,:] != 0.] -= inps.refPhase[i]            

    ##### Read/Generate Mask
    if inps.maskDataset and inps.maskDataset in stackobj.datasetNames:
        print('reading {} in {} * {} ...'.format(inps.maskDataset, box, inps.numIfgram))
        msk_data = stackobj.read(datasetName=inps.maskDataset, box=box, dropIfgram=True).reshape(inps.numIfgram,-1)
        if inps.maskDataset == 'coherence':
            msk_data = msk_data >= inps.maskThreshold
            print('mask out pixels with {} < {}'.format(inps.maskDataset, inps.maskThreshold))
        else:
            print('mask out pixels with {} == 0'.format(inps.maskDataset))
        pha_data[msk_data==0.] = 0.
    return pha_data


def read_coherence2weight(stackobj, inps, box):
    epsilon = 1e-4
    print('reading coherence in {} * {} ...'.format(box, inps.numIfgram))
    coh_data = stackobj.read(datasetName='coherence', box=box, dropIfgram=True).reshape(inps.numIfgram,-1)
    coh_data[np.isnan(coh_data)] = epsilon

    ##### Calculate Weight matrix
    weight = np.array(coh_data, np.float64)
    L = int(stackobj.metadata['ALOOKS']) * int(stackobj.metadata['RLOOKS'])
    if inps.weightFunc == 'var':
        print('convert coherence to weight using inverse of phase variance')
        print('    with phase PDF for distributed scatterers from Tough et al. (1995)')
        weight = 1.0 / coherence2phase_variance_ds(weight, L, print_msg=True)

    elif inps.weightFunc == 'coh':
        print('use coherence as weight directly (Perissin & Wang, 2012; Tong et al., 2016)')
        weight[weight < epsilon] = epsilon

    elif inps.weightFunc == 'fim':
        print('convert coherence to weight using Fisher Information Index (Seymour & Cumming, 1994)')
        weight = coherence2fisher_info_index(weight, L)

    else:
        print('Un-recognized weight function: %s' % inps.weightFunc)
        sys.exit(-1)
    return weight


def ifgram_inversion_patch(stackobj, inps, box=None):
    '''
    Inputs:
        ifgramStackFile    - string, interferograms hdf5 file
        coherenceFile - string, coherence hdf5 file
        box           - 4-tuple, left, upper, right, and lower pixel coordinate of area of interest
        meta          - dict, including the following attributes:

                        #Interferograms
                        length/width - int, file size for each interferogram
                        ifgram_list  - list of string, interferogram dataset name
                        date12_list  - list of string, YYMMDD-YYMMDD
                        ref_value    - np.array in size of (ifgram_num, 1)
                                       reference pixel coordinate in row/column number
                        ref_y/x      - int, reference pixel coordinate in row/column number

                        #Time-series
                        date8_list   - list of string in YYYYMMDD
                        tbase_diff   - np.array in size of (date_num-1, 1), differential temporal baseline

                        #Inversion
                        weightFunc   - sbas, fim, var, coh
    Outputs:
        ts       - 3D np.array in size of (date_num, row_num, col_num)
        temp_coh - 2D np.array in size of (row_num, col_num)
        tsStd    - 3D np.array in size of (date_num, row_num, col_num)
    '''

    ##### Get patch size/index
    if not box:
        box = (0,0,stackobj.width,stackobj.length)
    print('processing %8d/%d lines ...' % (box[3], stackobj.length))
    row_num = box[3] - box[1]
    col_num = box[2] - box[0]
    pixel_num = row_num * col_num
    date_num = stackobj.numDate

    ts    = np.zeros((date_num, pixel_num), np.float32)
    tsStd = np.zeros((date_num, pixel_num), np.float32)
    temp_coh = np.zeros(pixel_num, np.float32)

    ##### Read unwrapPhase
    pha_data = read_unwrap_phase(stackobj, inps, box)


    ##### Mask for pixels to invert
    mask = np.ones(pixel_num, np.bool_)
    ## 1 - Water Mask
    if inps.waterMaskFile:
        print('skip pixels on water with mask from file: %s' % (os.path.basename(inps.waterMaskFile)))
        try:    waterMask = readfile.read(inps.waterMaskFile, datasetName='waterMask', box=box)[0].flatten()
        except: waterMask = readfile.read(inps.waterMaskFile, datasetName='mask', box=box)[0].flatten()
        mask *= np.array(waterMask, np.bool_)

    ## 2 - Mask for Zero Phase in ALL ifgrams
    print('skip pixels with zero/nan value in all interferograms')
    phaseStack = np.nanmean(pha_data, axis=0)
    mask *= np.multiply(~np.isnan(phaseStack), phaseStack != 0.)

    ## Invert pixels on mask 1+2
    pixel_num2inv = np.sum(mask)
    pixel_idx2inv = np.where(mask)[0]
    print('number of pixels to invert: %s out of %s' % (pixel_num2inv, pixel_num))
    if pixel_num2inv < 1:
        ts = ts.reshape(date_num, row_num, col_num)
        temp_coh = temp_coh.reshape(row_num, col_num)
        tsStd = tsStd.reshape(date_num, row_num, col_num)
        return ts, temp_coh, tsStd


    ##### Design matrix
    A, B = stackobj.get_design_matrix(dropIfgram=True)
    try:    ref_date = str(np.loadtxt('reference_date.txt', dtype=bytes).astype(str))
    except: ref_date = stackobj.dateList[0]
    refIdx = stackobj.dateList.index(ref_date)
    timeIdx = [i for i in range(date_num)]
    timeIdx.remove(refIdx)
    Astd = stackobj.get_design_matrix(refDate=ref_date, dropIfgram=True)[0]


    ##### Inversion - SBAS
    if inps.weightFunc == 'sbas':
        ## Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
        maskAllNet = np.all(pha_data, axis=0)
        maskAllNet *= mask
        maskPartNet = mask ^ maskAllNet

        if np.sum(maskAllNet) > 0:
            print('inverting pixels with valid phase in all  ifgrams: {} pixels ...'.format(np.sum(maskAllNet)))
            numAllNet = np.sum(maskAllNet)
            phaDataTemp = pha_data[:,maskAllNet]
            ts1 = np.zeros((date_num-1,numAllNet))
            tempCoh1 = np.zeros((numAllNet))
            step = 1000
            loopNum = int(np.floor(numAllNet/step))
            progBar = ptime.progress_bar(maxValue=loopNum)
            for i in range(loopNum):
                [i0, i1] = [i*step, min((i+1)*step, numAllNet)]
                ts1[:,i0:i1], tempCoh1[i0:i1] = network_inversion_sbas(B, phaDataTemp[:,i0:i1], inps.tbaseDiff,\
                                                                       skipZeroPhase=False)
                progBar.update(i+1, suffix=i0)
            progBar.close()
            #ts1, tempCoh1 = network_inversion_sbas(B, pha_data[:,maskAllNet], inps.tbaseDiff, skipZeroPhase=False)
            ts[1:,maskAllNet] = ts1
            temp_coh[maskAllNet] = tempCoh1

        if np.sum(maskPartNet) > 0:
            print('inverting pixels with valid phase in some ifgrams (%.0f pixels) ...' % (np.sum(maskPartNet)))
            pixel_num2inv = np.sum(maskPartNet)
            pixel_idx2inv = np.where(maskPartNet)[0]
            progBar = ptime.progress_bar(maxValue=pixel_num2inv)
            for i in range(pixel_num2inv):
                idx = pixel_idx2inv[i]
                ts1, tempCoh1 = network_inversion_sbas(B, pha_data[:,idx], inps.tbaseDiff, inps.skip_zero_phase)
                ts[1:, idx] = ts1.flatten()
                temp_coh[idx] = tempCoh1
                progBar.update(i+1, every=100, suffix=str(i+1)+'/'+str(pixel_num2inv)+' pixels')
            progBar.close()


    ##### Inversion - WLS
    else:
        weight = read_coherence2weight(stackobj, inps, box)

        ##### Weighted Inversion pixel by pixel
        print('inverting time series ...')
        progBar = ptime.progress_bar(maxValue=pixel_num2inv)
        for i in range(pixel_num2inv):
            idx = pixel_idx2inv[i]
            ts1, tempCoh1, tsStd1 = network_inversion_wls(A, pha_data[:,idx], weight[:,idx], Astd=Astd,\
                                                          skipZeroPhase=inps.skip_zero_phase)
            ts[1:, idx] = ts1.flatten()
            temp_coh[idx] = tempCoh1
            tsStd[timeIdx, idx] = tsStd1.flatten()
            progBar.update(i+1, every=100, suffix=str(i+1)+'/'+str(pixel_num2inv)+' pixels')
        progBar.close()


    ts = ts.reshape(date_num, row_num, col_num)
    tsStd = tsStd.reshape(date_num, row_num, col_num)
    temp_coh = temp_coh.reshape(row_num, col_num)
    return ts, temp_coh, tsStd


def ifgram_inversion(ifgramStackFile='ifgramStack.h5', inps=None):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 

    Parameters: ifgramStackFile : string,
                    HDF5 file name of the interferograms stck
                inps : namespace, including the following options:
    Returns:    timeseriesFile : string
                    HDF5 file name of the output timeseries
                tempCohFile : string
                    HDF5 file name of temporal coherence
    Example:
        inps = cmdLineParse()
        ifgram_inversion('ifgramStack.h5', inps)
    '''
    total = time.time()

    ##Check Inputs
    if not inps:
        inps = cmdLineParse()
    inps.timeseriesStdFile = os.path.join(os.path.dirname(inps.timeseriesFile),'timeseriesDecorStd.h5')

    if inps.update_mode and not ut.update_file(inps.timeseriesFile, ifgramStackFile):
        return inps.timeseriesFile, inps.tempCohFile

    ##### IfgramStack Info
    stackobj = ifgramStack(ifgramStackFile)
    stackobj.open()
    inps.numIfgram = np.sum(stackobj.dropIfgram)
    print('number of interferograms: {}'.format(inps.numIfgram))
    print('number of acquisitions  : {}'.format(stackobj.numDate))
    print('number of lines  : {}'.format(stackobj.length))
    print('number of columns: {}'.format(stackobj.width))

    inps.dateList = stackobj.get_date_list(dropIfgram=True)
    inps.tbase = ptime.date_list2tbase(inps.dateList)[0]
    inps.tbaseDiff = np.diff(inps.tbase).reshape(-1,1)

    inps.refPhase = check_ifgram_reference_phase(stackobj, inps)
    check_design_matrix(stackobj, inps)

    ##### Invert time-series phase
    boxList = split_into_boxes(inps, inps.numIfgram, stackobj.length, stackobj.width)
    numBox = len(boxList)
    ts    = np.zeros((stackobj.numDate, stackobj.length, stackobj.width), np.float32)
    tsStd = np.zeros((stackobj.numDate, stackobj.length, stackobj.width), np.float32)
    tempCoh = np.zeros((stackobj.length, stackobj.width), np.float32)
    for i in range(numBox):
        if numBox > 1:
            print('\n------- Processing Patch %d out of %d --------------' % (i+1, numBox))
        box = boxList[i]
        tsi, tcohi, tsStdi = ifgram_inversion_patch(stackobj, inps, box)
        tempCoh[box[1]:box[3],box[0]:box[2]] = tcohi
        ts[:,box[1]:box[3],box[0]:box[2]] = tsi
        tsStd[:,box[1]:box[3],box[0]:box[2]] = tsStdi

    print('converting phase to range')
    phase2range = -1*float(stackobj.metadata['WAVELENGTH'])/(4.*np.pi)
    ts *= phase2range
    tsStd *= abs(phase2range)

    ##### Output 
    print('calculating perpendicular baseline timeseries')
    pbase = stackobj.get_perp_baseline_timeseries(dropIfgram=True)
    atr = dict(stackobj.metadata)
    stackobj.close()

    atr['REF_DATE'] = inps.dateList[0]
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'

    print('-'*50)
    tsobj = timeseries(inps.timeseriesFile)
    tsobj.write2hdf5(data=ts, dates=inps.dateList, bperp=pbase, metadata=atr)

    if not np.all(tsStd == 0.):
        print('-'*50)
        tsobj = timeseries(inps.timeseriesStdFile)
        tsobj.write2hdf5(data=tsStd, refFile=inps.timeseriesFile)

    print('-'*50)
    print('writing >>> '+inps.tempCohFile)
    atr['FILE_TYPE'] = 'coherence'
    atr['UNIT'] = '1'
    writefile.write(tempCoh, atr, inps.tempCohFile)

    print('network inversion took {:.1f} seconds\nDone.'.format(time.time()-total))
    return inps.timeseriesFile, inps.tempCohFile


def split_into_boxes(inps, numIfgram, length, width, print_msg=True):
    '''Split into chunks in rows to reduce memory usage'''
    #Get r_step / chunk_num
    r_step = inps.chunk_size / (numIfgram * width)         #split in lines
    if inps.weightFunc != 'sbas':  #more memory usage (coherence) for WLS
        r_step /= 2.0
    r_step = int(ceil_to_1(r_step))
    chunk_num = int((length-1)/r_step) + 1

    if print_msg and chunk_num > 1:
        print('maximum chunk size: %.1E' % (inps.chunk_size))
        print('split %d lines into %d patches for processing' % (length, chunk_num))
        print('    with each patch up to %d lines' % (r_step))

    ##Computing the inversion
    boxList = []
    for i in range(chunk_num):
        r0 = i*r_step
        r1 = min([length, r0+r_step])
        box = (0,r0,width,r1)
        boxList.append(box)
    return boxList


def check_design_matrix(stackobj, inps):
    '''Check Rank of Design matrix for weighted inversion'''
    A = stackobj.get_design_matrix(dropIfgram=True)[0]
    print('-------------------------------------------------------------------------------')
    if inps.weightFunc == 'sbas':
        print('generic least square inversion with min-norm phase velocity')
        print('    based on Berardino et al. (2002, IEEE-TGRS)')
        print('    OLS for pixels with full rank      network')
        print('    SVD for pixels with rank deficient network')
        if np.linalg.matrix_rank(A) < stackobj.numDate-1:
            print('WARNING: singular design matrix! Inversion result can be biased!')
            print('continue using its SVD solution on all pixels')
    else:
        print('weighted least square (WLS) inversion with min-norm phase, pixelwise')
        if np.linalg.matrix_rank(A) < stackobj.numDate-1:
            print('ERROR: singular design matrix!')
            print('    Input network of interferograms is not fully connected!')
            print('    Can not invert the weighted least square solution.')
            print('You could try:')
            print('    1) Add more interferograms to make the network fully connected:')
            print('       a.k.a., no multiple subsets nor network islands')
            print("    2) Use '-w no' option for non-weighted SVD solution.")
            sys.exit(-1)
    print('-------------------------------------------------------------------------------')
    return


def check_ifgram_reference_phase(stackobj, inps):
    '''Read refPhase'''
    try:
        ref_y = int(stackobj.metadata['REF_Y'])
        ref_x = int(stackobj.metadata['REF_X'])
        box = [ref_x, ref_y, ref_x+1, ref_y+1]
        inps.refPhase = np.squeeze(stackobj.read(datasetName='unwrapPhase', box=box, dropIfgram=True))
        print('reference pixel in y/x: {}'.format((ref_y, ref_x)))
    except:
        if inps.skip_ref:
            inps.refPhase = 0.0
            print('skip checking reference pixel info - This is for SIMULATION ONLY.')
        else:
            print('ERROR: No REF_X/Y found! Can not invert interferograms without reference in space.')
            print('run reference_point.py '+inps.ifgramStackFile+' --mark-attribute for a quick referencing.')
            sys.exit(1)
    return inps.refPhase


################################################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    if inps.templateFile:
        inps = read_template2inps(inps.templateFile, inps)
    inps.timeseriesFile, inps.tempCohFile = inps.outfile

    # Input file info
    atr = readfile.read_attribute(inps.ifgramStackFile)
    if atr['FILE_TYPE'] != 'ifgramStack':
        print('ERROR: only ifgramStack file supported, input is '+k+' file!')
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

