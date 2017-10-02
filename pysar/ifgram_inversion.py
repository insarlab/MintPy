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
def timeseries_inversion(ifgramFile='unwrapIfgram.h5', coherenceFile='coherence.h5', inps_dict=None):
    '''Implementation of the SBAS algorithm.
    modified from sbas.py written by scott baker, 2012 

    Inputs:
        ifgramFile    - string, HDF5 file name of the interferograms
        coherenceFile - string, HDF5 file name of the coherence
        inps_dict     - dict, including the following options:
                        weight_function
                        min_coherence
                        max_coherence
    Output:
        timeseriesFile - string, HDF5 file name of the output timeseries
        tempCohFile    - string, HDF5 file name of temporal coherence
    '''
    total = time.time()

    if not inps_dict:
        inps_dict = vars(cmdLineParse())
    weight_func = inps_dict['weight_function']
    min_coh = inps_dict['min_coherence']
    max_coh = inps_dict['max_coherence']

    # Basic Info
    atr = readfile.read_attribute(ifgramFile)
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    pixel_num = length * width

    h5ifgram = h5py.File(ifgramFile,'r')
    ifgram_list = sorted(h5ifgram['interferograms'].keys())
    if inps_dict['weight_function'] == 'no':
        ifgram_list = ut.check_drop_ifgram(h5ifgram, atr, ifgram_list)
    ifgram_num = len(ifgram_list)

    # Convert ifgram_list to date12/8_list
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    m_dates = [i.split('-')[0] for i in date12_list]
    s_dates = [i.split('-')[1] for i in date12_list]
    date8_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    date_num = len(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    tbase_diff = np.diff(tbase_list).reshape((date_num-1,1))

    print 'number of interferograms: '+str(ifgram_num)
    print 'number of acquisitions  : '+str(date_num)
    print 'number of pixels: '+str(pixel_num)

    # Reference pixel in space
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print 'reference pixel in y/x: [%d, %d]'%(ref_y, ref_x)
    except:
        print 'ERROR: No ref_x/y found! Can not inverse interferograms without reference in space.'
        print 'run seed_data.py '+ifgramFile+' --mark-attribute for a quick referencing.'
        sys.exit(1)

    ##### Read Interferograms
    print 'reading interferograms ...'
    ifgram_data = np.zeros((ifgram_num, pixel_num), np.float32)
    prog_bar = ptime.progress_bar(maxValue=ifgram_num)
    for j in range(ifgram_num):
        ifgram = ifgram_list[j]
        d = h5ifgram['interferograms'][ifgram].get(ifgram)[:]
        #d[d != 0.] -= d[ref_y, ref_x]
        d -= d[ref_y, ref_x]
        ifgram_data[j] = d.flatten()
        prog_bar.update(j+1, suffix=date12_list[j])
    h5ifgram.close()
    prog_bar.close()


    #####---------------------- Inversion ----------------------#####
    # Design matrix
    A,B = ut.design_matrix(ifgramFile, date12_list)

    if weight_func == 'no':
        print 'generalized inversion using SVD (Berardino et al., 2002, IEEE-TGRS)'
        print 'inversing time series ...'
        B_inv = np.array(np.linalg.pinv(B), np.float32)
        ts_rate = np.dot(B_inv, ifgram_data)
        ts1 = ts_rate * np.tile(tbase_diff, (1, pixel_num))
        ts0 = np.array([0.]*pixel_num, np.float32)
        ts_data  = np.vstack((ts0, np.cumsum(ts1, axis=0)))
        del ts_rate, ts0, ts1

        # Temporal coherence
        print 'calculating temporal coherence (Tizzani et al., 2007, RSE)'
        temp_coh = np.zeros((1, pixel_num), np.float32)+0j
        prog_bar = ptime.progress_bar(maxValue=ifgram_num)
        for i in range(ifgram_num):
            ifgram_est = np.dot(A[i,:], ts_data[1:,:])
            ifgram_diff = ifgram_data[i,:] - ifgram_est
            temp_coh += np.exp(1j*ifgram_diff)
            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()
        del ifgram_data, ifgram_est, ifgram_diff
        temp_coh = np.array((np.absolute(temp_coh)/ifgram_num).reshape((length,width)), dtype=np.float32)

    else:
        print 'weighted least square (WLS) inversion using coherence pixel by pixel'
        if np.linalg.matrix_rank(A) < date_num-1:
            print 'ERROR: singular design matrix!'
            print '    Input network of interferograms is not fully connected!'
            print '    Can not inverse the weighted least square solution.'
            print 'You could try:'
            print '    1) Add more interferograms to make the network fully connected:'
            print '       a.k.a., no multiple subsets nor network islands'
            print "    2) Use '-w no' option for non-weighted SVD solution."
            sys.exit(-1)

        pixel_mask = np.ones(pixel_num, np.bool_)
        print 'reading coherence: '+os.path.basename(coherenceFile)
        h5coh = h5py.File(coherenceFile,'r')
        coh_list = sorted(h5coh['coherence'].keys())
        coh_data = np.zeros((ifgram_num, pixel_num), np.float32)
        prog_bar = ptime.progress_bar(maxValue=ifgram_num)
        for j in range(ifgram_num):
            ifgram = coh_list[j]
            d = h5coh['coherence'][ifgram].get(ifgram)[:].flatten()
            d[np.isnan(d)] = 0.
            pixel_mask[d == 0.] = 0
            coh_data[j] = d
            prog_bar.update(j+1, suffix=date12_list[j])
        h5coh.close()
        prog_bar.close()

        # Get mask of valid pixels to inverse
        print 'skip pixels with zero coherence in at least one interferogram'
        print 'skip pixels with zero phase     in all          interferograms'
        ifgram_stack = ut.get_file_stack(ifgramFile).flatten()
        pixel_mask[ifgram_stack == 0.] = 0

        pixel_num2inv = np.sum(pixel_mask)
        pixel_idx2inv = np.where(pixel_mask)[0]
        ifgram_data = ifgram_data[:,pixel_mask]
        coh_data = coh_data[:,pixel_mask]
        print 'number of pixels to inverse: %d' % (pixel_num2inv)

        ##### Calculate Weight matrix
        weight = coh_data
        if weight_func.startswith('var'):
            print 'convert coherence to weight using inverse of variance: x**2/(1-x**2) from Hanssen (2001, for 4.2.32)'
            weight[weight > 0.999] = 0.999
            if weight_func == 'variance-max-coherence':
                print 'constrain the max coherence to %f' % max_coh
                weight[weight > max_coh] = max_coh
            weight = np.square(weight)
            weight *= 1. / (1. - weight)
            if weight_func == 'variance-log':
                print 'use log(1/variance)+1 as weight'
                weight = np.log(weight+1)
        elif weight_func.startswith('lin'):
            print 'use coherence as weight directly (Tong et al., 2016, RSE)'
        elif weight_func.startswith('norm'):
            print 'convert coherence to weight using CDF of normal distribution: N(%f, %f)' % (mu, std)
            mu  = (min_coh + max_coh) / 2.0
            std = (max_coh - min_coh) / 6.0
            chunk_size = 1000
            chunk_num = int(pixel_num2inv/chunk_size)+1
            prog_bar = ptime.progress_bar(maxValue=chunk_num)
            for i in range(chunk_num):
                i0 = (i-1)*chunk_size
                i1 = min([pixel_num2inv, i0+chunk_size])
                weight[:,i0:i1] = norm.cdf(weight[:,i0:i1], mu, std)
                prog_bar.update(i+1, every=10)
            prog_bar.close()
            #weight = norm.cdf(weight, mu, std)
        else:
            print 'Un-recognized weight function: %s' % weight_func
            sys.exit(-1)

        ##### Weighted Inversion pixel by pixel
        print 'inversing time series ...'
        ts_data = np.zeros((date_num, pixel_num), np.float32)
        temp_coh = np.zeros(pixel_num, np.float32)
        prog_bar = ptime.progress_bar(maxValue=pixel_num2inv)
        for i in range(pixel_num2inv):
            # Inverse timeseries
            ifgram_pixel = ifgram_data[:,i]
            weight_pixel = weight[:,i]
            W = np.diag(weight_pixel)
            ts = np.linalg.inv(A.T.dot(W).dot(A)).dot(A.T).dot(W).dot(ifgram_pixel)
            ts_data[1:,pixel_idx2inv[i]] = ts

            # Calculate weighted temporal coherence
            ifgram_diff = ifgram_pixel - np.dot(A, ts)
            temp_coh_pixel = np.abs(np.sum(np.multiply(weight_pixel, np.exp(1j*ifgram_diff)), axis=0)) / np.sum(weight_pixel)
            temp_coh[pixel_idx2inv[i]] = temp_coh_pixel

            prog_bar.update(i+1, every=2000, suffix=str(i+1)+' pixels')
        prog_bar.close()
        del ifgram_data, weight


    #####---------------------- Outputs ----------------------#####
    ## 1.1 Convert time-series phase to displacement
    print 'converting phase to range'
    phase2range = -1*float(atr['WAVELENGTH'])/(4.*np.pi)
    ts_data *= phase2range

    ## 1.2 Write time-series data matrix
    timeseriesFile = 'timeseries.h5'
    print 'writing >>> '+timeseriesFile
    print 'number of acquisitions: '+str(date_num)
    h5timeseries = h5py.File(timeseriesFile,'w')
    group = h5timeseries.create_group('timeseries')
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date8_list[i]
        dset = group.create_dataset(date, data=ts_data[i].reshape(length,width), compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    ## 1.3 Write time-series attributes
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
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5timeseries.close()
    del ts_data

    ## 2. Write Temporal Coherence File
    tempCohFile = 'temporalCoherence.h5'
    print 'writing >>> '+tempCohFile
    atr['FILE_TYPE'] = 'temporal_coherence'
    atr['UNIT'] = '1'
    writefile.write(temp_coh.reshape(length,width), atr, tempCohFile)

    print 'Time series inversion took ' + str(time.time()-total) +' secs\nDone.'
    return timeseriesFile, tempCohFile


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
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Inverse network of interferograms into timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be inversed')
    parser.add_argument('--template','-t', dest='template_file',\
                        help='template text file with the following options:\n'+TEMPLATE)
    parser.add_argument('--coherence','-c', dest='coherence_file', default='coherence.h5', help='coherence file')
    parser.add_argument('--weight-func','-w', dest='weight_function', default='no',\
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
        timeseries_inversion(inps.ifgram_file, inps.coherence_file, vars(inps))
    else:
        print 'inverse time-series using L1 norm minimization'
        ut.timeseries_inversion_L1(inps.ifgram_file, inps.timeseries_file)

    return


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

