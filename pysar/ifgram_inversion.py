#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
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

import h5py
import numpy as np
from scipy.stats import norm

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


################################################################################################
def round_to_1(x):
    '''Return the most significant digit of input number'''
    return round(x, -int(np.floor(np.log10(abs(x)))))


def network_inversion_sbas(B, ifgram, tbase_diff, B_inv=None):
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
        B_inv      - 2D np.array in size of (date_num-1, ifgram_num)
                     (Moore-Penrose) pseudo-inverse of design matrix B
    Output:
        ts - 2D np.array in size of (date_num-1, pixel_num), phase time series
    '''
    if B_inv is None:
        B_inv = np.array(np.linalg.pinv(B), np.float32)

    ts_rate = np.dot(B_inv, ifgram)
    ts_diff = ts_rate * np.tile(tbase_diff, (1, ifgram.shape[1]))
    ts = np.cumsum(ts_diff, axis=0)
    return ts


def network_inversion_wls(A, ifgram, weight, ts=None):
    '''Network inversion based on Weighted Least Square (WLS) solution.
    Inputs:
        A      - 2D np.array in size of (ifgram_num, date_num-1)
                 representing date configuration for each interferogram
                 (-1 for master, 1 for slave, 0 for others)
        ifgram - 1D np.array in size of (ifgram_num,)
                 phase of all interferograms
        weight - 1D np.array in size of (ifgram_num,), weight of ifgram
    Output:
        ts - 1D np.array in size of (date_num-1,), phase time series
    '''
    W = np.diag(weight)
    ts = np.linalg.inv(A.T.dot(W).dot(A)).dot(A.T).dot(W).dot(ifgram)
    return ts


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
        prog_bar = ptime.progress_bar(maxValue=chunk_num)
        for i in range(chunk_num):
            p0 = i*chunk_size
            p1 = min([p0+chunk_size, pixel_num])
            ifgram_diff = ifgram[:,p0:p1] - np.dot(A, ts[:,p0:p1])
            temp_coh[p0:p1] = np.abs(np.sum(np.multiply(weight[:,p0:p1], np.exp(1j*ifgram_diff)), axis=0)) /\
                              np.sum(weight[:,p0:p1], axis=0)
            prog_bar.update(i+1, every=10, suffix=str(p1)+'/'+str(pixel_num))
        prog_bar.close()
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
                        weight_function   - 
                        max/min_coherence - 
    '''

    ##### Get patch size/index
    if not box:
        box = (0,0,meta['width'],meta['length'])
    c0,r0,c1,r1 = box

    ## Initiate output data matrixs
    row_num = r1-r0
    col_num = c1-c0
    pixel_num = row_num * col_num
    date_num = len(meta['date8_list'])
    ts = np.zeros((date_num, pixel_num), np.float32)
    temp_coh = np.zeros(pixel_num, np.float32)

    ##### Get mask of non-zero pixels
    print 'skip pixels with zero/nan value in all interferograms'
    ifgram_stack = ut.get_file_stack(ifgramFile)[r0:r1,c0:c1].flatten()
    mask = ~np.isnan(ifgram_stack)
    mask[ifgram_stack == 0.] = 0
    pixel_num2inv = np.sum(mask)
    pixel_idx2inv = np.where(mask)[0]
    print 'number of pixels to inverse: %d' % (pixel_num2inv)
    if pixel_num2inv < 1:
        ts = ts.reshape(date_num, row_num, col_num)
        temp_coh = temp_coh.reshape(row_num, col_num)
        return ts, temp_coh

    ##### Read interferograms
    ifgram_num = len(meta['ifgram_list'])
    ifgram_data = np.zeros((ifgram_num, pixel_num2inv), np.float32)
    date12_list = meta['date12_list']

    print 'reading interferograms ...'
    h5ifgram = h5py.File(ifgramFile,'r')
    prog_bar = ptime.progress_bar(maxValue=ifgram_num)
    for j in range(ifgram_num):
        ifgram = meta['ifgram_list'][j]
        d = h5ifgram['interferograms'][ifgram].get(ifgram)[r0:r1,c0:c1]
        ifgram_data[j] = d.flatten()[mask]
        prog_bar.update(j+1, suffix=date12_list[j])
    h5ifgram.close()
    prog_bar.close()
    ifgram_data -= meta['ref_value']

    ##### Design matrix
    A,B = ut.design_matrix(ifgramFile, date12_list)
    B_inv = np.array(np.linalg.pinv(B), np.float32)

    ##### Inverse
    if meta['weight_function'] in ['no','uniform']:
        print 'inversing time-series ...'
        ts[1:,mask] = network_inversion_sbas(B, ifgram_data, meta['tbase_diff'], B_inv=B_inv)

        print 'calculating temporal coherence ...'
        temp_coh[mask] = temporal_coherence(A, ts[1:,mask], ifgram_data)

    else:
        ##### Read coherence
        coh_data = np.zeros((ifgram_num, pixel_num2inv), np.float32)
        print 'reading coherence ...'
        h5coh = h5py.File(coherenceFile,'r')
        coh_list = sorted(h5coh['coherence'].keys())
        prog_bar = ptime.progress_bar(maxValue=ifgram_num)
        for j in range(ifgram_num):
            ifgram = coh_list[j]
            d = h5coh['coherence'][ifgram].get(ifgram)[r0:r1,c0:c1]
            d[np.isnan(d)] = 0.
            coh_data[j] = d.flatten()[mask]
            prog_bar.update(j+1, suffix=date12_list[j])
        h5coh.close()
        prog_bar.close()

        ##### Calculate Weight matrix
        weight = coh_data
        if meta['weight_function'].startswith('var'):
            print 'convert coherence to weight using inverse of variance: x**2/(1-x**2)'
            weight[weight > 0.999] = 0.999
            weight[weight < 0.001] = 0.001
            #if meta['weight_function'] == 'variance-max-coherence':
            #    print 'constrain the max coherence to %f' % meta['max_coherence']
            #    weight[weight > meta['max_coherence']] = meta['max_coherence']
            weight = np.square(weight)
            weight *= 1. / (1. - weight)
            #if meta['weight_function'] == 'variance-log':
            #    print 'use log(1/variance)+1 as weight'
            #    weight = np.log(weight+1)
        elif meta['weight_function'].startswith('lin'):
            print 'use coherence as weight directly (Tong et al., 2016, RSE)'
            weight[weight < 1e-6] = 1e-6
        elif meta['weight_function'].startswith('norm'):
            mu  = (meta['min_coherence'] + meta['max_coherence']) / 2.0
            std = (meta['max_coherence'] - meta['min_coherence']) / 6.0
            print 'convert coherence to weight using CDF of normal distribution: N(%f, %f)' % (mu, std)
            chunk_size = 1000
            chunk_num = int(pixel_num2inv-1/chunk_size)+1
            prog_bar = ptime.progress_bar(maxValue=chunk_num)
            for k in range(chunk_num):
                p0 = (k-1)*chunk_size
                p1 = min([pixel_num2inv, p0+chunk_size])
                weight[:,p0:p1] = norm.cdf(weight[:,p0:p1], mu, std)
                prog_bar.update(k+1, every=10)
            prog_bar.close()
            #weight = norm.cdf(weight, mu, std)
        else:
            print 'Un-recognized weight function: %s' % meta['weight_function']
            sys.exit(-1)

        ##### Weighted Inversion pixel by pixel
        print 'inversing time series ...'
        prog_bar = ptime.progress_bar(maxValue=pixel_num2inv)
        for i in range(pixel_num2inv):
            ts[1:,pixel_idx2inv[i]] = network_inversion_wls(A, ifgram_data[:,i], weight[:,i])
            prog_bar.update(i+1, every=1000, suffix=str(i+1)+'/'+str(pixel_num2inv)+' pixels')
        prog_bar.close()

        print 'calculating temporal coherence ...'
        temp_coh[mask] = temporal_coherence(A, ts[1:,mask], ifgram_data, weight)

    ts = ts.reshape(date_num, row_num, col_num)
    temp_coh = temp_coh.reshape(row_num, col_num)
    return ts, temp_coh


def ifgram_inversion(ifgramFile='unwrapIfgram.h5', coherenceFile='coherence.h5', meta=None):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 

    Inputs:
        ifgramFile    - string, HDF5 file name of the interferograms
        coherenceFile - string, HDF5 file name of the coherence
        meta          - dict, including the following options:
                        weight_function
                        min_coherence
                        max_coherence
                        chunk_size - float, max number of data (ifgram_num*row_num*col_num)
                                     to read per loop; to control the memory
    Output:
        timeseriesFile - string, HDF5 file name of the output timeseries
        tempCohFile    - string, HDF5 file name of temporal coherence
    '''
    total = time.time()

    if not meta:
        meta = vars(cmdLineParse())

    ##### Basic Info
    # length/width
    atr = readfile.read_attribute(ifgramFile)
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    meta['length'] = length
    meta['width']  = width

    # ifgram_list
    h5ifgram = h5py.File(ifgramFile,'r')
    ifgram_list = sorted(h5ifgram['interferograms'].keys())
    if meta['weight_function'] in ['no','uniform']:
        ifgram_list = ut.check_drop_ifgram(h5ifgram, atr, ifgram_list)
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
    tbase_diff = np.diff(tbase_list).reshape((date_num-1,1))
    meta['tbase_diff'] = tbase_diff

    print 'number of interferograms: %d' % (ifgram_num)
    print 'number of acquisitions  : %d' % (date_num)
    print 'number of columns: %d' % (width)
    print 'number of lines  : %d' % (length)

    ##### ref_y/x/value
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print 'reference pixel in y/x: [%d, %d]' % (ref_y, ref_x)
        ref_value = np.zeros((ifgram_num,1), np.float32)
        for j in range(ifgram_num):
            ifgram = ifgram_list[j]
            dset = h5ifgram['interferograms'][ifgram].get(ifgram)
            ref_value[j] = dset[ref_y,ref_x]
        meta['ref_y'] = ref_y
        meta['ref_x'] = ref_x
        meta['ref_value'] = ref_value
    except:
        print 'ERROR: No ref_x/y found! Can not inverse interferograms without reference in space.'
        print 'run seed_data.py '+ifgramFile+' --mark-attribute for a quick referencing.'
        sys.exit(1)
    h5ifgram.close()

    ##### Rank of Design matrix for weighted inversion
    A, B = ut.design_matrix(ifgramFile, date12_list)
    print '-------------------------------------------------------------------------------'
    if meta['weight_function'] in ['no','uniform']:
        print 'ordinary least square (OLS) inversion with min-norm phase velocity'
        print '    based on Berardino et al. (2002, IEEE-TGRS)'
        if np.linalg.matrix_rank(A) < date_num-1:
            print 'WARNING: singular design matrix! Inversion result can be biased!'
            print 'continue using its SVD solution'
    else:
        print 'weighted least square (WLS) inversion with min-norm phase, pixelwise'
        if np.linalg.matrix_rank(A) < date_num-1:
            print 'ERROR: singular design matrix!'
            print '    Input network of interferograms is not fully connected!'
            print '    Can not inverse the weighted least square solution.'
            print 'You could try:'
            print '    1) Add more interferograms to make the network fully connected:'
            print '       a.k.a., no multiple subsets nor network islands'
            print "    2) Use '-w no' option for non-weighted SVD solution."
            sys.exit(-1)
    print '-------------------------------------------------------------------------------'


    ##### Inverse time-series phase
    timeseries = np.zeros((date_num, length, width), np.float32)
    temp_coh = np.zeros((length, width), np.float32)

    r_step = meta['chunk_size']/ifgram_num/width
    if meta['weight_function'] not in ['no','uniform']:
        r_step /= 2
    r_step = int(round_to_1(r_step))
    chunk_num = int((length-1)/r_step)+1
    print 'maximum chunk size: %.1E' % (meta['chunk_size'])
    print 'divide the whole dataset into %d patches for processing' % (chunk_num)
    print 'each patch has up to %d*%d pixels' % (r_step, width)
    for i in range(chunk_num):
        r0 = i*r_step
        r1 = min([length, r0+r_step])
        print 'processing %8d/%d lines ...' % (r1, length)
        box = (0,r0,width,r1)
        timeseries[:,r0:r1,:], temp_coh[r0:r1,:] = ifgram_inversion_patch(ifgramFile, coherenceFile, meta, box)

    ##### Calculate time-series attributes
    print 'calculating perpendicular baseline timeseries'
    pbase, pbase_top, pbase_bottom = ut.perp_baseline_ifgram2timeseries(ifgramFile, ifgram_list)
    pbase = str(pbase.tolist()).translate(None,'[],')  # convert np.array into string separated by white space
    pbase_top = str(pbase_top.tolist()).translate(None,'[],')
    pbase_bottom = str(pbase_bottom.tolist()).translate(None,'[],')
    atr['P_BASELINE_TIMESERIES'] = pbase
    atr['P_BASELINE_TOP_TIMESERIES'] = pbase_top
    atr['P_BASELINE_BOTTOM_TIMESERIES'] = pbase_bottom
    atr['ref_date'] = date8_list[0]
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'


    ##### Output
    ## 1. Write time-series file
    timeseriesFile = write_timeseries_hdf5_file(timeseries, date8_list, atr)

    ## 2. Write Temporal Coherence File
    tempCohFile = 'temporalCoherence.h5'
    print 'writing >>> '+tempCohFile
    atr['FILE_TYPE'] = 'temporal_coherence'
    atr['UNIT'] = '1'
    writefile.write(temp_coh, atr, tempCohFile)

    print 'Time series inversion took ' + str(time.time()-total) +' secs\nDone.'
    return timeseriesFile, tempCohFile


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

    ## 1 Convert time-series phase to displacement
    print 'converting phase to range'
    phase2range = -1*float(atr['WAVELENGTH'])/(4.*np.pi)
    timeseries *= phase2range

    ## 2 Write time-series data matrix
    if not timeseriesFile:
        timeseriesFile = 'timeseries.h5'
    print 'writing >>> '+timeseriesFile

    date_num = len(date8_list)
    print 'number of acquisitions: '+str(date_num)
    h5timeseries = h5py.File(timeseriesFile,'w')
    group = h5timeseries.create_group('timeseries')
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date8_list[i]
        dset = group.create_dataset(date, data=timeseries[i], compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5timeseries.close()

    return timeseriesFile


def read_template2inps(template_file, inps):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(template_file)
    key_list = template.keys()

    # Coherence-based network modification
    prefix = 'pysar.timeseriesInv.'

    key = prefix+'residualNorm'
    if key in key_list and template[key] in ['L1']:
        inps.resid_norm = 'L1'
    else:
        inps.resid_norm = 'L2'

    key = prefix+'coherenceFile'
    if key in key_list:
        value = template[key]
        if value in ['auto']:
            inps.coherence_file = 'coherence.h5'
        elif value in ['no']:
            inps.coherence_file = None
        else:
            inps.coherence_file = value

    key = prefix+'minCoherence'
    if key in key_list:
        value = template[key]
        if value in ['auto']:
            inps.min_coherence = 0.2
        else:
            inps.min_coherence = float(value)

    key = prefix+'maxCoherence'
    if key in key_list:
        value = template[key]
        if value in ['auto']:
            inps.max_coherence = 0.85
        else:
            inps.max_coherence = float(value)

    key = prefix+'weightFunc'
    if key in key_list:
        value = template[key]
        if value in ['auto','no']:
            inps.weight_function = 'no'
        elif value.startswith('norm'):
            inps.weight_function = 'normal'
        elif value.startswith('lin'):
            inps.weight_function = 'linear'
        elif value.startswith('var'):
            inps.weight_function = 'variance'
        else:
            print 'Un-recognized input for %s = %s' % (key, value)
            sys.exit(-1)

    return inps


################################################################################################
EXAMPLE='''example:
  ifgram_inversion.py  unwrapIfgram.h5
  ifgram_inversion.py  unwrapIfgram.h5 -t pysarApp_template.txt
  ifgram_inversion.py  unwrapIfgram.h5 -c coherence.h5 -w variance
'''

TEMPLATE='''
pysar.timeseriesInv.residualNorm  = auto #[L2 ], auto for L2, norm minimization solution
pysar.timeseriesInv.coherenceFile = auto #[fname / no], auto for coherence.h5, file to read weight data
pysar.timeseriesInv.minCoherence  = auto #[0.0-1.0], auto for 0.20, put 0 weight for pixels with coherence < input
pysar.timeseriesInv.maxCoherence  = auto #[0.0-1.0], auto for 0.85, put 1 weight for pixels with coherence > input
pysar.timeseriesInv.weightFunc    = auto #[variance / no / linear / normal], auto for no, coherence to weight
                                         #variance - phase variance due to temporal decorrelation
                                         #no - no weight, or ordinal inversion with uniform weight
                                         #linear - uniform distribution CDF function
                                         #normal - normal  distribution CDF function
'''

REFERENCE='''references:
Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002). A new algorithm for surface 
    deformation monitoring based on small baseline differential SAR interferograms. IEEE TGRS,
    40(11), 2375-2383. doi:10.1109/TGRS.2002.803792
Tizzani, P., Berardino, P., Casu, F., Euillades, P., Manzo, M., Ricciardi, G. P., Lanari, R.
    (2007). Surface deformation of Long Valley caldera and Mono Basin, California, investigated
    with the SBAS-InSAR approach. Remote Sensing of Environment, 108(3), 277-289.
    doi:http://dx.doi.org/10.1016/j.rse.2006.11.015
Tong, X., & Schmidt, D. (2016). Active movement of the Cascade landslide complex in Washington
    from a coherence-based InSAR time series method. Remote Sensing of Environment, 186, 405-415.
    doi:http://dx.doi.org/10.1016/j.rse.2016.09.008
Just, D., & Bamler, R. (1994). Phase statistics of interferograms with applications to synthetic
    aperture radar. Applied optics, 33(20), 4361-4368. 
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Inverse network of interferograms into timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be inversed')
    parser.add_argument('--template','-t', dest='template_file',\
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--coherence','-c', dest='coherence_file', default='coherence.h5', help='coherence file')
    parser.add_argument('--weight-function','-w', dest='weight_function', default='no',\
                        help='function used to convert coherence to weight for inversion:\n'+\
                             'variance - phase variance due to temporal decorrelation\n'+\
                             'no       - no weight, or ordinal inversion with uniform weight\n'+\
                             'linear   - uniform distribution CDF function\n'+\
                             'normal   - normal  distribution CDF function')
    parser.add_argument('--min-coherence','--min-coh', dest='min_coherence', type=float, default=0.2,\
                        help='set min weight for pixels with coherence < minimum coherence')
    parser.add_argument('--max-coherence','--max-coh', dest='max_coherence', type=float, default=0.85,\
                        help='set max weight for pixels with coherence < maximum coherence')
    parser.add_argument('--norm', dest='resid_norm', default='L2', choices=['L1','L2'],\
                        help='Inverse method used to residual optimization, L1 or L2 norm minimization. Default: L2')

    parser.add_argument('--chunk-size', dest='chunk_size', type=float, default=0.5e9,\
                        help='max number of data (= ifgram_num * row_num * col_num) to read per loop\n'+\
                             'default: 0.5G; adjust it according to your computer memory.')
    #parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
    #                    help='Disable parallel processing.')

    inps = parser.parse_args()
    return inps


################################################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # Input file info
    atr = readfile.read_attribute(inps.ifgram_file)
    k = atr['FILE_TYPE']
    if not k == 'interferograms':
        sys.exit('ERROR: only interferograms file supported, input is '+k+' file!')

    # Network Inversion
    if inps.resid_norm == 'L2':
        print 'inverse time-series using L2 norm minimization'
        ifgram_inversion(inps.ifgram_file, inps.coherence_file, meta=vars(inps))
    else:
        print 'inverse time-series using L1 norm minimization'
        ut.timeseries_inversion_L1(inps.ifgram_file, inps.timeseries_file)

    return


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

