#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################
# timeseries_inversion are modified from a software originally
# written by Scott Baker with the following licence:
###############################################################################
#  Copyright (c) 2011, Scott Baker 
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
############################################################################### 


import os
import sys
import time
import argparse
import string

import h5py
import numpy as np
from scipy.special import gamma

import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile
import pysar.utils.utils as ut


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
        if B.shape[0] < dateNum1:
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
        ATW = A.T.dot(W)
        ts = np.linalg.inv(ATW.dot(A)).dot(ATW).dot(ifgram.reshape(-1,1))
        #A_inv_wls = np.linalg.inv(A.T.dot(W).dot(A))
        #ts = A_inv_wls.dot(A.T).dot(W).dot(ifgram.reshape(-1,1))

        ## Temporal Coherence
        ifgram_diff = ifgram - np.dot(A, ts).flatten()
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


def ifgram_inversion_patch(ifgramFile, coherenceFile, meta, box=None):
    '''
    Inputs:
        ifgramFile    - string, interferograms hdf5 file
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
                        weight_function   - no, fim, var, coh
    Outputs:
        ts       - 3D np.array in size of (date_num, row_num, col_num)
        temp_coh - 2D np.array in size of (row_num, col_num)
        tsStd    - 3D np.array in size of (date_num, row_num, col_num)
    '''

    ##### Get patch size/index
    if not box:
        box = (0,0,meta['width'],meta['length'])
    c0,r0,c1,r1 = box
    print('processing %8d/%d lines ...' % (r1, meta['length']))

    ## Initiate output data matrixs
    row_num = r1-r0
    col_num = c1-c0
    pixel_num = row_num * col_num
    date_num = len(meta['date8_list'])
    ts = np.zeros((date_num, pixel_num), np.float32)
    tsStd = np.zeros((date_num, pixel_num), np.float32)
    temp_coh = np.zeros(pixel_num, np.float32)

    ##### Mask for pixels to invert
    mask = np.ones(pixel_num, np.bool_)
    ## 1 - Water Mask
    if meta['water_mask_file']:
        print('skip pixels on water with mask from file: %s' % (os.path.basename(meta['water_mask_file'])))
        try:    waterMask = readfile.read(meta['water_mask_file'], epoch='waterMask')[0][r0:r1,c0:c1].flatten()
        except: waterMask = readfile.read(meta['water_mask_file'], epoch='mask')[0][r0:r1,c0:c1].flatten()
        mask *= np.array(waterMask, np.bool_)

    ## 2 - Mask for Zero Phase in ALL ifgrams
    print('skip pixels with zero/nan value in all interferograms')
    ifgram_stack = ut.get_file_stack(ifgramFile)[r0:r1,c0:c1].flatten()
    mask *= ~np.isnan(ifgram_stack)
    mask *= ifgram_stack != 0.

    ## Invert pixels on mask 1+2
    pixel_num2inv = np.sum(mask)
    pixel_idx2inv = np.where(mask)[0]
    print('number of pixels to invert: %s out of %s' % (pixel_num2inv, pixel_num))
    if pixel_num2inv < 1:
        ts = ts.reshape(date_num, row_num, col_num)
        temp_coh = temp_coh.reshape(row_num, col_num)
        tsStd = tsStd.reshape(date_num, row_num, col_num)
        return ts, temp_coh, tsStd

    ##### Read interferograms
    ifgram_num = len(meta['ifgram_list'])
    ifgram_data = np.zeros((ifgram_num, pixel_num), np.float32)
    date12_list = meta['date12_list']

    if meta['skip_zero_phase']:
        print('skip zero phase value (masked out and filled during phase unwrapping)')
    atr = readfile.read_attribute(ifgramFile)
    h5ifgram = h5py.File(ifgramFile,'r')
    for j in range(ifgram_num):
        ifgram = meta['ifgram_list'][j]
        d = h5ifgram['interferograms'][ifgram].get(ifgram)[r0:r1,c0:c1].flatten()
        if meta['skip_zero_phase']:
            d[d != 0.] -= meta['ref_value'][j]
        else:
            d -= meta['ref_value'][j]
        ifgram_data[j] = d
        sys.stdout.write('\rreading interferograms %s/%s ...' % (j+1, ifgram_num))
        sys.stdout.flush()
    print(' ')
    h5ifgram.close()
    #ifgram_data -= meta['ref_value']

    ## 3 - Mask for Non-Zero Phase in ALL ifgrams (share one B in sbas inversion)
    maskAllNet = np.all(ifgram_data, axis=0)
    maskAllNet *= mask
    maskPartNet = mask ^ maskAllNet

    ##### Design matrix
    A,B = ut.design_matrix(ifgramFile, date12_list)
    try:    ref_date = str(np.loadtxt('reference_date.txt', dtype=bytes).astype(str))
    except: ref_date = meta['date8_list'][0]
    #print 'calculate decorrelation noise covariance with reference date = %s' % (ref_date)
    refIdx = meta['date8_list'].index(ref_date)
    timeIdx = [i for i in range(date_num)]
    timeIdx.remove(refIdx)
    Astd = ut.design_matrix(ifgramFile, date12_list, referenceDate=ref_date)[0]

    ##### Inversion
    if meta['weight_function'] in ['no','uniform']:
        if np.sum(maskAllNet) > 0:
            print('inverting pixels with valid phase in all  ifgrams (%.0f pixels) ...' % (np.sum(maskAllNet)))
            ts1, tempCoh1 = network_inversion_sbas(B, ifgram_data[:,maskAllNet], meta['tbase_diff'], skipZeroPhase=False)
            ts[1:,maskAllNet] = ts1
            temp_coh[maskAllNet] = tempCoh1

        if np.sum(maskPartNet) > 0:
            print('inverting pixels with valid phase in some ifgrams ...')
            pixel_num2inv = np.sum(maskPartNet)
            pixel_idx2inv = np.where(maskPartNet)[0]
            prog_bar = ptime.progress_bar(maxValue=pixel_num2inv)
            for i in range(pixel_num2inv):
                idx = pixel_idx2inv[i]
                ts1, tempCoh1 = network_inversion_sbas(B, ifgram_data[:,idx], meta['tbase_diff'], meta['skip_zero_phase'])
                ts[1:, idx] = ts1.flatten()
                temp_coh[idx] = tempCoh1
                prog_bar.update(i+1, every=100, suffix=str(i+1)+'/'+str(pixel_num2inv)+' pixels')
            prog_bar.close()

    else:
        ##### Read coherence
        coh_data = np.zeros((ifgram_num, pixel_num), np.float32)
        h5coh = h5py.File(coherenceFile,'r')
        coh_list = sorted(h5coh['coherence'].keys())
        coh_list = ut.check_drop_ifgram(h5coh)
        for j in range(ifgram_num):
            ifgram = coh_list[j]
            d = h5coh['coherence'][ifgram].get(ifgram)[r0:r1,c0:c1]
            d[np.isnan(d)] = 0.
            coh_data[j] = d.flatten()
            sys.stdout.write('\rreading coherence %s/%s ...' % (j+1, ifgram_num))
            sys.stdout.flush()
        print(' ')
        h5coh.close()

        ##### Calculate Weight matrix
        weight = np.array(coh_data, np.float64)
        L = int(atr['ALOOKS']) * int(atr['RLOOKS'])
        epsilon = 1e-4
        if meta['weight_function'].startswith('var'):
            print('convert coherence to weight using inverse of phase variance')
            print('    with phase PDF for distributed scatterers from Tough et al. (1995)')
            weight = 1.0 / coherence2phase_variance_ds(weight, L, print_msg=True)

        elif meta['weight_function'].startswith(('lin','coh','cor')):
            print('use coherence as weight directly (Perissin & Wang, 2012; Tong et al., 2016)')
            weight[weight < epsilon] = epsilon

        elif meta['weight_function'].startswith(('fim','fisher')):
            print('convert coherence to weight using Fisher Information Index (Seymour & Cumming, 1994)')
            weight = coherence2fisher_info_index(weight, L)

        else:
            print('Un-recognized weight function: %s' % meta['weight_function'])
            sys.exit(-1)

        ##### Weighted Inversion pixel by pixel
        print('inverting time series ...')
        prog_bar = ptime.progress_bar(maxValue=pixel_num2inv)
        for i in range(pixel_num2inv):
            idx = pixel_idx2inv[i]
            ts1, tempCoh1, tsStd1 = network_inversion_wls(A, ifgram_data[:,idx], weight[:,idx], Astd=Astd,\
                                                          skipZeroPhase=meta['skip_zero_phase'])
            ts[1:, idx] = ts1.flatten()
            temp_coh[idx] = tempCoh1
            tsStd[timeIdx, idx] = tsStd1.flatten()
            prog_bar.update(i+1, every=100, suffix=str(i+1)+'/'+str(pixel_num2inv)+' pixels')
        prog_bar.close()

    ts = ts.reshape(date_num, row_num, col_num)
    temp_coh = temp_coh.reshape(row_num, col_num)
    tsStd = tsStd.reshape(date_num, row_num, col_num)


    ##Write to temp hdf5 files for parallel processing
    if meta['parallel']:
        fname = meta['ftemp_base']+str(int(r0/meta['row_step']))+'.h5'
        print('writing >>> '+fname)
        h5temp = h5py.File(fname, 'w')
        group = h5temp.create_group('timeseries')
        dset = group.create_dataset('timeseries', shape=(date_num+1, row_num, col_num), dtype=np.float32)
        dset[0:-1,:,:] = ts
        dset[1,:,:] = temp_coh
        h5temp.close()
        return
    else:
        return ts, temp_coh, tsStd


def ifgram_inversion(ifgramFile='unwrapIfgram.h5', coherenceFile='coherence.h5', meta=None):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 

    Inputs:
        ifgramFile    - string, HDF5 file name of the interferograms
        coherenceFile - string, HDF5 file name of the coherence
        meta          - dict, including the following options:
                        weight_function
                        chunk_size - float, max number of data (ifgram_num*row_num*col_num)
                                     to read per loop; to control the memory
    Output:
        timeseriesFile - string, HDF5 file name of the output timeseries
        tempCohFile    - string, HDF5 file name of temporal coherence
    Example:
        meta = dict()
        meta['weight_function'] = 'variance'
        meta['chunk_size'] = 0.5e9
        meta['timeseriesFile'] = 'timeseries_var.h5'
        meta['tempCohFile'] = 'temporalCoherence_var.h5'
        ifgram_inversion('unwrapIfgram.h5', 'coherence.h5', meta)
    '''
    if 'tempCohFile' not in meta.keys():
        meta['tempCohFile'] = 'temporalCoherence.h5'
    meta['timeseriesStdFile'] = 'timeseriesDecorStd.h5'
    total = time.time()

    if not meta:
        meta = vars(cmdLineParse())

    if meta['update_mode'] and not ut.update_file(meta['timeseriesFile'], ifgramFile):
        return meta['timeseriesFile'], meta['tempCohFile']

    ##### Basic Info
    # length/width
    atr = readfile.read_attribute(ifgramFile)
    length = int(atr['LENGTH'])
    width  = int(atr['WIDTH'])
    meta['length'] = length
    meta['width']  = width

    # ifgram_list
    h5ifgram = h5py.File(ifgramFile,'r')
    ifgram_list = sorted(h5ifgram['interferograms'].keys())
    #if meta['weight_function'] in ['no','uniform']:
    #    ifgram_list = ut.check_drop_ifgram(h5ifgram)
    ifgram_list = ut.check_drop_ifgram(h5ifgram)
    meta['ifgram_list'] = ifgram_list
    ifgram_num = len(ifgram_list)

    # date12_list/date8_list/tbase_diff
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    m_dates = [i.split('-')[0] for i in date12_list]
    s_dates = [i.split('-')[1] for i in date12_list]
    date8_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    date_num = len(date8_list)
    meta['date8_list'] = date8_list
    meta['date12_list'] = date12_list

    tbase_list = ptime.date_list2tbase(date8_list)[0]
    tbase_diff = np.diff(tbase_list).reshape((-1,1))
    meta['tbase_diff'] = tbase_diff

    print('number of interferograms: %d' % (ifgram_num))
    print('number of acquisitions  : %d' % (date_num))
    print('number of columns: %d' % (width))
    print('number of lines  : %d' % (length))

    ##### ref_y/x/value
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print('reference pixel in y/x: [%d, %d]' % (ref_y, ref_x))
        ref_value = np.zeros((ifgram_num,1), np.float32)
        for j in range(ifgram_num):
            ifgram = ifgram_list[j]
            dset = h5ifgram['interferograms'][ifgram].get(ifgram)
            ref_value[j] = dset[ref_y,ref_x]
        meta['ref_y'] = ref_y
        meta['ref_x'] = ref_x
        meta['ref_value'] = ref_value
    except:
        if meta['skip_ref']:
            meta['ref_value'] = 0.0
            print('skip checking reference pixel info - This is for SIMULATION ONLY.')
        else:
            print('ERROR: No ref_x/y found! Can not invert interferograms without reference in space.')
            print('run reference_point.py '+ifgramFile+' --mark-attribute for a quick referencing.')
            sys.exit(1)
    h5ifgram.close()

    ##### Rank of Design matrix for weighted inversion
    A, B = ut.design_matrix(ifgramFile, date12_list)
    print('-------------------------------------------------------------------------------')
    if meta['weight_function'] in ['no','uniform']:
        print('generic least square inversion with min-norm phase velocity')
        print('    based on Berardino et al. (2002, IEEE-TGRS)')
        print('    OLS for pixels with full rank      network')
        print('    SVD for pixels with rank deficient network')
        if np.linalg.matrix_rank(A) < date_num-1:
            print('WARNING: singular design matrix! Inversion result can be biased!')
            print('continue using its SVD solution on all pixels')
    else:
        print('weighted least square (WLS) inversion with min-norm phase, pixelwise')
        if np.linalg.matrix_rank(A) < date_num-1:
            print('ERROR: singular design matrix!')
            print('    Input network of interferograms is not fully connected!')
            print('    Can not invert the weighted least square solution.')
            print('You could try:')
            print('    1) Add more interferograms to make the network fully connected:')
            print('       a.k.a., no multiple subsets nor network islands')
            print("    2) Use '-w no' option for non-weighted SVD solution.")
            sys.exit(-1)
    print('-------------------------------------------------------------------------------')


    ##### Invert time-series phase
    ##Check parallel environment
    if meta['weight_function'] in ['no','uniform']:
        meta['parallel'] = False
    if meta['parallel']:
        num_cores, meta['parallel'], Parallel, delayed = ut.check_parallel(1000, print_msg=False)

    ##Split into chunks to reduce memory usage
    r_step = meta['chunk_size']/ifgram_num/width         #split in lines
    if meta['weight_function'] not in ['no','uniform']:  #more memory usage (coherence) for WLS
        r_step /= 2.0
        if meta['parallel']:
            r_step /= num_cores
    r_step = int(ceil_to_1(r_step))
    meta['row_step'] = r_step
    chunk_num = int((length-1)/r_step)+1

    if chunk_num > 1:
        print('maximum chunk size: %.1E' % (meta['chunk_size']))
        print('split %d lines into %d patches for processing' % (length, chunk_num))
        print('    with each patch up to %d lines' % (r_step))
        if meta['parallel']:
            print('parallel processing using %d cores ...' % (min([num_cores,chunk_num])))

    ##Computing the inversion
    box_list = []
    for i in range(chunk_num):
        r0 = i*r_step
        r1 = min([length, r0+r_step])
        box = (0,r0,width,r1)
        box_list.append(box)
    box_num = len(box_list)

    if not meta['parallel']:
        timeseries = np.zeros((date_num, length, width), np.float32)
        timeseriesStd = np.zeros((date_num, length, width), np.float32)
        tempCoh = np.zeros((length, width), np.float32)
        for i in range(box_num):
            if box_num > 1:
                print('\n------- Processing Patch %d out of %d --------------' % (i+1, box_num))
            box = box_list[i]
            ts, tcoh, tsStd = ifgram_inversion_patch(ifgramFile, coherenceFile, meta, box)
            tempCoh[box[1]:box[3],box[0]:box[2]] = tcoh
            timeseries[:,box[1]:box[3],box[0]:box[2]] = ts
            timeseriesStd[:,box[1]:box[3],box[0]:box[2]] = tsStd

    else:
        ##Temp file list
        meta['ftemp_base'] = 'timeseries_temp_'
        temp_file_list = [meta['ftemp_base']+str(i)+'.h5' for i in range(chunk_num)]

        ##Computation
        Parallel(n_jobs=num_cores)(delayed(ifgram_inversion_patch)\
                                   (ifgramFile, coherenceFile, meta, box) for box in box_list)

        ##Concatenate temp files
        print('concatenating temporary timeseries files ...')
        timeseries = np.zeros((date_num, length, width), np.float32)
        tempCoh = np.zeros((length, width), np.float32)
        rmCmd = 'rm'
        for i in range(chunk_num):
            fname = temp_file_list[i]
            box = box_list[i]
            print('reading '+fname)
            h5temp = h5py.File(fname, 'r')
            dset = h5temp['timeseries'].get('timeseries')
            timeseries[:,box[1]:box[3],box[0]:box[2]] = dset[0:-1,:,:]
            tempCoh[box[1]:box[3],box[0]:box[2]] = dset[-1,:,:]
            h5temp.close()
            rmCmd += ' '+fname
        print(rmCmd)
        os.system(rmCmd)

    print('converting phase to range')
    phase2range = -1*float(atr['WAVELENGTH'])/(4.*np.pi)
    timeseries *= phase2range
    timeseriesStd *= abs(phase2range)

    ##### Calculate time-series attributes
    print('calculating perpendicular baseline timeseries')
    pbase, pbase_top, pbase_bottom = ut.perp_baseline_ifgram2timeseries(ifgramFile, ifgram_list)
    # convert np.array into string separated by white space
    pbase = str(pbase.tolist()).translate(str.maketrans('[],','   ')).strip()
    pbase_top = str(pbase_top.tolist()).translate(str.maketrans('[],','   ')).strip()
    pbase_bottom = str(pbase_bottom.tolist()).translate(str.maketrans('[],','   ')).strip()
    atr['P_BASELINE_TIMESERIES'] = pbase
    atr['P_BASELINE_TOP_TIMESERIES'] = pbase_top
    atr['P_BASELINE_BOTTOM_TIMESERIES'] = pbase_bottom
    atr['ref_date'] = date8_list[0]
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'

    ##### Output
    ## 1. Write time-series file
    meta['timeseriesFile'] = write_timeseries_hdf5_file(timeseries, date8_list, atr,\
                                                        timeseriesFile=meta['timeseriesFile'])
    if not np.all(timeseriesStd == 0.):
        meta['timeseriesStdFile'] = write_timeseries_hdf5_file(timeseriesStd, date8_list, atr,\
                                                               timeseriesFile=meta['timeseriesStdFile'])

    ## 2. Write Temporal Coherence File
    print('writing >>> '+meta['tempCohFile'])
    atr['FILE_TYPE'] = 'temporal_coherence'
    atr['UNIT'] = '1'
    meta['tempCohFile'] = writefile.write(tempCoh, atr, meta['tempCohFile'])

    print('Time series inversion took ' + str(time.time()-total) +' secs\nDone.')
    return meta['timeseriesFile'], meta['tempCohFile']


def write_timeseries_hdf5_file(timeseries, date8_list, atr, timeseriesFile=None):
    ''' Write to timeseries HDF5 file
    Inputs:
        timeseries - 3D np.array in size of (date_num, length, width)
                     cumulative time series phase
        date8_list - list of string in YYYYMMDD format
        atr        - dict, attributes of time-series file, including two parts:
                     1) attributes inherited from interferograms
                     2) attributes of time-series inverted from network of interferograms:
                         P_BASELINE_TIMESERIES
                         P_BASELINE_TOP_TIMESERIES
                         P_BASELINE_BOTTOM_TIMESERIES
                         ref_date
        timeseriesFile - string, file name of output time-series file
    Output:
        timeseriesFile - string, file name of output time-series file
    '''
    if not timeseriesFile:
        timeseriesFile = 'timeseries.h5'
    print('writing >>> '+timeseriesFile)

    date_num = len(date8_list)
    print('number of acquisitions: '+str(date_num))
    h5timeseries = h5py.File(timeseriesFile,'w')
    group = h5timeseries.create_group('timeseries')
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date8_list[i]
        dset = group.create_dataset(date, data=timeseries[i], compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    for key,value in iter(atr.items()):
        group.attrs[key] = value
    h5timeseries.close()

    return timeseriesFile


def read_template2inps(template_file, inps):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(template_file)

    # Coherence-based network modification
    prefix = 'pysar.networkInversion.'

    key = prefix+'residualNorm'
    if key in template.keys() and template[key] in ['L1']:
        inps.resid_norm = 'L1'
    else:
        inps.resid_norm = 'L2'

    key = prefix+'coherenceFile'
    if key in template.keys():
        value = template[key]
        if value in ['auto']:
            inps.coherence_file = 'coherence.h5'
        elif value in ['no']:
            inps.coherence_file = None
        else:
            inps.coherence_file = value

    key = prefix+'weightFunc'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.weight_function = 'no'
        elif value.startswith(('lin','coh','cor')):
            inps.weight_function = 'linear'
        elif value.startswith('var'):
            inps.weight_function = 'variance'
        elif value.startswith(('fim','fisher')):
            inps.weight_function = 'fim'
        else:
            print('Un-recognized input for %s = %s' % (key, value))
            sys.exit(-1)

    key = prefix+'waterMaskFile'
    if key in template.keys():
        value = template[key]
        if value in ['auto', 'no']:
            maskFile = None
            #atr = readfile.read_attribute(inps.ifgram_file)
            #if 'Y_FIRST' in atr.keys():
            #    maskFile = 'geometryGeo.h5'
            #else:
            #    maskFile = 'geometryRadar.h5'
        else:
            maskFile = value
            try:
                data = readfile.read(maskFile, epoch='mask')[0]
                inps.water_mask_file = maskFile
            except:
                print('Can not found mask dataset in file: %s' % (maskFile))
                print('Ignore this input water mask file option and continue.')

    return inps


################################################################################################
EXAMPLE='''example:
  ifgram_inversion.py  unwrapIfgram.h5
  ifgram_inversion.py  unwrapIfgram.h5 -t pysarApp_template.txt
  ifgram_inversion.py  unwrapIfgram.h5 -w var
  ifgram_inversion.py  unwrapIfgram.h5 -w fim
  ifgram_inversion.py  unwrapIfgram.h5 -w coh
'''

TEMPLATE='''
## Invert network of interferograms into time series using weighted least sqaure (WLS) estimator.
## Temporal coherence is calculated using Tazzani et al. (Tizzani et al., 2007, IEEE-TGRS)
## Singular-Value Decomposition (SVD) is applied if network are not fully connected for no weight scenario.
## There are 4 weighting options:
## a. fim       - WLS, use Fisher Information Matrix as weight (Seymour & Cumming, 1994, IGARSS). [Recommended]
## b. variance  - WLS, use inverse of covariance as weight (Guarnieri & Tebaldini, 2008, TGRS)
## c. coherence - WLS, use coherence as weight (Perissin & Wang, 2012, IEEE-TGRS)
## d. no        - LS, no/uniform weight (Berardino et al., 2002, TGRS)
pysar.networkInversion.weightFunc    = auto #[fim / variance / coherence / no], auto for no
pysar.networkInversion.coherenceFile = auto #[filename / no], auto for coherence.h5, file to read weight data
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
Perissin, D., and T. Wang (2012), Repeat-pass SAR interferometry with partially coherent targets, IEEE TGRS,
    50(1), 271-280, doi:10.1109/tgrs.2011.2160644.
Samiei-Esfahany, S., J. E. Martins, F. v. Leijen, and R. F. Hanssen (2016), Phase Estimation for Distributed
    Scatterers in InSAR Stacks Using Integer Least Squares Estimation, IEEE TGRS, 54(10), 5671-5687.
Seymour, M. S., and I. G. Cumming (1994), Maximum likelihood estimation for SAR interferometry, 1994. 
    IGARSS '94., 8-12 Aug 1994.
Tizzani, P., Berardino, P., Casu, F., Euillades, P., Manzo, M., Ricciardi, G. P., Lanari, R.
    (2007). Surface deformation of Long Valley caldera and Mono Basin, California, investigated
    with the SBAS-InSAR approach. Remote Sensing of Environment, 108(3), 277-289.
    doi:http://dx.doi.org/10.1016/j.rse.2006.11.015
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Invert network of interferograms into timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be inverted')
    parser.add_argument('--template','-t', dest='template_file',\
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--ref-date', dest='ref_date', help='Reference date, first date by default.')
    parser.add_argument('--coherence','-c', dest='coherence_file', default='coherence.h5', help='coherence file')
    parser.add_argument('--weight-function','-w', dest='weight_function', default='no',\
                        help='function used to convert coherence to weight for inversion:\n'+\
                             'variance - phase variance due to temporal decorrelation\n'+\
                             'linear   - uniform distribution CDF function\n'+\
                             'no       - no weight, or ordinal inversion with uniform weight')
    parser.add_argument('--norm', dest='resid_norm', default='L2', choices=['L1','L2'],\
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
    parser.add_argument('--water-mask','-m', dest='water_mask_file', help='Skip inversion on the masked out region, i.e. water.')
    inps = parser.parse_args()
    inps.parallel = False
    return inps


################################################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    inps.timeseriesFile = inps.outfile[0]
    inps.tempCohFile = inps.outfile[1]

    # Input file info
    atr = readfile.read_attribute(inps.ifgram_file)
    k = atr['FILE_TYPE']
    if not k == 'interferograms':
        sys.exit('ERROR: only interferograms file supported, input is '+k+' file!')

    # Network Inversion
    if inps.resid_norm == 'L2':
        print('inverse time-series using L2 norm minimization')
        ifgram_inversion(inps.ifgram_file, inps.coherence_file, meta=vars(inps))
    else:
        print('inverse time-series using L1 norm minimization')
        ut.timeseries_inversion_L1(inps.ifgram_file, inps.timeseriesFile)

    return


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

