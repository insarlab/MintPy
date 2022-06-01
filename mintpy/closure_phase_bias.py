#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yujie Zheng, Feb 2022                            #
############################################################
# Compute average con-nl closure phase and output mask identifying areas suseptible to closure phase errors.

import os
import sys
import time
import argparse
import numpy as np
import glob
from datetime import datetime as dt

from mintpy.objects import ifgramStack, cluster
from mintpy.utils import readfile, writefile, ptime, isce_utils, arg_group
from mintpy import ifgram_inversion as ifginv

################################################################################
REFERENCE = """reference:
  Y. Zheng, H. Fattahi, P. Agram, M. Simons and P. Rosen, "On Closure Phase and Systematic Bias in Multi-looked SAR Interferometry," in IEEE Transactions on Geoscience and Remote Sensing, doi: 10.1109/TGRS.2022.3167648.
"""
EXAMPLE = """example:
    closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --action create_mask
    closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --num_sigma 2.5 --action create_mask
    closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 10 -a quick_bias_estimate
    closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 10 -a bias_estimate --noupdate -c local
"""

def create_parser():
    parser = argparse.ArgumentParser(description = 'Mask / estimate / correct for phase non-closure related biases.',
                                    formatter_class = argparse.RawTextHelpFormatter,
                                    epilog=REFERENCE+'\n'+EXAMPLE)
    parser.add_argument('-i','--ifgramstack',type = str, dest = 'ifgram_stack',help = 'interferogram stack file that contains the unwrapped phases')
    parser.add_argument('--nl', dest = 'nl', type = int, default = 20, help = 'connection level that we are correcting to (or consider as no bias)')
    parser.add_argument('--bw', dest = 'bw', type = int, default = 10, help = 'bandwidth of time-series analysis that you want to correct')
    parser.add_argument('--num_sigma',dest = 'num_sigma', type = float, default = 3, help = 'Threashold for phase (number of sigmas,0-infty), default to be 3 sigma of a Gaussian distribution (assumed distribution for the cumulative closure phase) with sigma = pi/sqrt(3*num_cp)')
    parser.add_argument('--epi',dest = 'episilon', type = float, default = 0.3, help = 'Threashold for amplitude (0-1), default 0.3')
    parser.add_argument('--noupdate',dest = 'update_cp', action = 'store_false', help = 'Use when no need to compute closure phases')
    parser.add_argument('-o', dest = 'outdir', type = str, default = '.', help = 'output file directory')
    parser.add_argument('-a','--action', dest='action', type=str, default='create_mask',
                        choices={'create_mask', 'quick_bias_estimate', 'bias_estimate'},
                        help='action to take (default: %(default)s):\n'+
                             'create_mask -  create a mask of areas susceptible to closure phase errors\n'+
                             'quick_bias_estimate - estimate how bias decays with time, will output sequential closure phase files, and gives a quick and appximate bias estimation\n'
                             'bias_estimate - estimate how bias decays with time, processed for each pixel on a pixel by pixel basis')
    parser = arg_group.add_parallel_argument(parser)
    parser = arg_group.add_memory_argument(parser)
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args = iargs)
    return inps

def seq_closure_phase(slc_list, date12_list_all, ifgram_stack, ref_phase, n, box):
    """Computes wrapped sequential closure phases of conneciton-n

    Parameters : slc_list           - list of string, SLC dates
                 date12_list_all    - list of string, date12 of all the interferograms stored in the ifgramstack file
                 ifgram_stack       - string, file path of ifgramStack.h5
                 ref_phase          - 1D array in size of (num_ifgram,) in float, unwrapped phase of the reference pixel
                 n                  - integer, connection level of the closure phase (e.g., triplets are connection-2)
                 box                - list in size of (4,) in integer, bounding box coordinates
    Returns :    cp_w               - 3D array in size of (num_ifgram, box_length, box_width) in float,  wrapped sequential closure phases of connection n
    """
    cp_idx = [] # initiate a list storing the index of interferograms in each closure phase computation
    nslc = len(slc_list)
    for i in range(nslc-n):
        ifgram = []
        flag = True
        for j in range(n):
            ifgram.append('{}_{}'.format(slc_list[i+j],slc_list[i+j+1]))
        ifgram.append('{}_{}'.format(slc_list[i],slc_list[i+n]))
        for ifgram_name in ifgram:
            if ifgram_name not in date12_list_all:
                flag = False # if missing an interferogram, we won't make the corresponding closure phase
        if flag:
            cp_idx.append([date12_list_all.index(ifgram[j]) for j in range(n+1)])

    cp_idx = np.array(cp_idx, np.int16)
    cp_idx = np.unique(cp_idx, axis = 0)

    num_cp = len(cp_idx)
    print('Number of con-', n, 'closure measurements expected, ', nslc-n)
    print('Number of con-', n, 'closure measurements found, ', num_cp)

    if num_cp < nslc-n:
        print('Missing interferograms, abort')
        raise Exception("Some interferograms are missing")

    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    phase = readfile.read(ifgram_stack, box=box,print_msg=False)[0]
    cp_w = np.zeros((num_cp, box_length, box_width), np.float32)
    for i in range(num_cp):
        cp0_w = np.zeros ((box_length, box_width), np.float32)
        for j in range(n):
                    idx = cp_idx[i,j]
                    cp0_w = cp0_w + phase[idx,:,:] - ref_phase[idx]
        idx = cp_idx[i,n]
        cp0_w = cp0_w - (phase[idx,:,:]-ref_phase[idx])
        cp_w[i,:,:] = np.angle(np.exp(1j*cp0_w))

    return cp_w


def sum_seq_closure_phase(slc_list, date12_list_all, ifgram_stack, ref_phase, n, box):
    """Computes the sum of consecutive complex sequential closure phase of connection n

    Parameters: slc_list            - list of string, SLC dates
                date12_list_all     - list of string, date12 of all the interferograms stored in the ifgramstack file
                ifgram_stack        - string, file path of ifgramStack.h5
                ref_phase           - 1D array in size of (num_ifgram,) in float, unwrapped phase of the reference pixel
                n                   - integer, connection level of the closure phase
                box                 - list in size of (4,) in integer, bounding box coordinates
    Returns:    cum_cp              - 2D array in size of (box_length, box_width) in complex64, sum of consecutive complex sequential closure phase of connection n
                num_cp              - integer, number of closure phases in the sum
    """
    cp_idx = []
    nslc = len(slc_list)
    for i in range(nslc-n):
        ifgram = []
        flag = True
        for j in range(n):
            ifgram.append('{}_{}'.format(slc_list[i+j],slc_list[i+j+1]))
        ifgram.append('{}_{}'.format(slc_list[i],slc_list[i+n]))
        for ifgram_name in ifgram:
            if ifgram_name not in date12_list_all:
                flag = False # if missing an interferogram, we won't make the corresponding closure phase
        if flag:
            cp_idx.append([date12_list_all.index(ifgram[j]) for j in range(n+1)])

    cp_idx = np.array(cp_idx, np.int16)
    cp_idx = np.unique(cp_idx, axis = 0)

    num_cp = len(cp_idx)
    print('Number of closure measurements expected, ', len(slc_list)-n)
    print('Number of closure measurements found, ', num_cp)

    if num_cp <1:
        print('No closure phase measurements found, abort')
        raise Exception("No triplets found!")

    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    phase = readfile.read(ifgram_stack, box=box,print_msg=False)[0]
    cum_cp = np.zeros((box_length, box_width), np.complex64)
    for i in range(num_cp):
        cp0_w = np.zeros ((box_length, box_width), np.float32)
        for j in range(n):
                    idx = cp_idx[i,j]
                    cp0_w = cp0_w + phase[idx,:,:] - ref_phase[idx]
        idx = cp_idx[i,n]
        cp0_w = cp0_w - (phase[idx,:,:]-ref_phase[idx])
        cum_cp = cum_cp + (np.exp(1j*cp0_w))

    return cum_cp, num_cp

def cum_seq_unw_closure_phase(n,filepath,length, width, refY, refX, slc_list, meta):
    '''output cumulative con-n sequential closure phase in time-series format (Eq. 25 in Zheng et al., 2022, but divided by n)

    Parameters: n               - integer, connection level of closure phases
                filepath        - string, filepath of sequential closure phases of connection - n
                width, length   - integer, width and length of the interferograms
                refY, refX      - ingeger, reference point coordinates
                slc_list        - list of string, SLC dates
                meta            - dict, metadata of ifgramStack.h5
    Returns:    conn_seqcumclosurephase.h5              - output hdf5 file of 3D array in size of (num_ifgram, length, width) in float, cumulative sequential closure phase
                conn_seqcumclosurephase_maskconcp.h5    - output hdf5 file of 2D array in size of (length, width) in bool, mask based on connected component
    '''
    outfiledir = os.path.join(filepath, 'conn'+str(n)+'_seqcumclosurephase.h5')
    outmaskdir = os.path.join(filepath, 'conn'+str(n)+'_seqcumclosurephase_maskconcp.h5')
    print('Creating '+outfiledir)
    if not os.path.isfile(outfiledir) or not os.path.isfile(outmaskdir):
        filelist = glob.glob(os.path.join(filepath, '*.unw'))
        numfile = len(filelist)
        filelist_st = sorted(filelist)
        cp_phase_all = np.zeros((numfile, length, width),np.float32)
        mask_all = np.zeros((numfile, length, width),np.float32)

        for i in range(numfile):
            file = filelist_st[i]
            concpfile = file.replace('.unw','.unw.conncomp')
            corfile = file.replace('.unw','.cor')
            cp =np.fromfile(file, dtype='float32')
            cp = cp.reshape([length,width*2])
            cor = np.fromfile(corfile, dtype = 'float32')
            concp = np.fromfile(concpfile, dtype = 'byte')
            cor = cor.reshape([length, width])
            concp = concp.reshape([length, width])
            cp_phase = cp[:,width:]
            cp_phase = cp_phase - cp_phase[refY,refX]
            mask_concp = np.where(concp >= 1, 1, np.nan)
            cp_phase_all[i,:,:] = cp_phase
            mask_all[i,:,:] = mask_concp

        # compute sequential closure phase
        N = len(slc_list)
        biasts = np.zeros([N,length, width], np.float32)
        biasts[0,:,:] = 0
        biasts[1:N-n+1,:,:]= np.cumsum(cp_phase_all,0)
        for i in range(N-n+1,N):
            biasts[i,:,:] = (i-N+n)*cp_phase_all[-1,:,:]+biasts[N-n,:,:]

        slc_list = np.array(slc_list, np.string_)
        dsDict = dict()
        dsDict = {'seqcum_closurephase': [np.float32, (N, length, width), biasts/n],
                          'date': [np.dtype('S8'),np.shape(slc_list), slc_list],}

        meta['FILE_TYPE'] = 'timeseries'
        writefile.layout_hdf5(outfiledir, dsDict, meta)

        mask = np.where(np.isnan(np.sum(mask_all,0)), False, True)
        dsDict = dict()
        dsDict = {'mask_concp': [np.bool_, (length, width),mask ],
                      'date': [np.dtype('S8'),np.shape(slc_list), slc_list],}
        meta['FILE_TYPE'] = 'mask'
        writefile.layout_hdf5(outmaskdir, dsDict, meta)

def seq2cum_closure_phase(conn, outdir, box):
    ''' this script read in cumulative sequential closure phase from individual closure phase directory (Eq. 25) in Zheng et al., 2022

    Parameters: conn    - integer, connection level of sequential closure phases
                outdir  - string, directory of conn_seqcumclosurephase.h5
                box     - list in size of (4,) in integer, coordinates of bounding box
    Returns:    biasts  - 3D array in size of (nslc, box_lengh, box_width) in float, cumulative sequential closure phases
    '''
    filepath = 'conn'+str(conn)+'_cp'
    filename = 'conn'+str(conn)+'_seqcumclosurephase.h5'
    seqcpfile = os.path.join(outdir, 'closurePhase', filepath, filename)
    biasts = readfile.read(seqcpfile, box=box,print_msg=False)[0]
    return biasts

def estimate_ratioX(tbase, n, nl, wvl, box, outdir, mask=False):
    ''' This script estimates w(n\delta_t)/w(delta_t), Eq.(29) in Zheng et al., 2022

    Parameters: tbase           - list in size of (nslc,) in float, time in accumulated years
                n               - integer, connection-level
                nl              - integer, minimum connection-level that we think is bias-free
                wvl             - float, wavelength
                box             - list in size of (4,) in integer, coordinates of bounding box
                outdir          - string, the working directory
                mask            - bool, whether to mask out areas with average bias velocity less than 1 mm/year
    Returns:    wratio          - 2D array of size (box_length, box_width) in float, w(n\delta_t)/w(delta_t), Eq.(29) in Zheng et al., 2022
                wratio_velocity - 2D array of size (box_length, box_width) in float, bias-velocity at n*delta_t temporal baseline
    '''
    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    cum_bias_conn_1 = seq2cum_closure_phase(nl, outdir, box)[-1,:,:]
    coef = -4*np.pi/wvl
    delta_T = tbase[-1]-tbase[0]
    vel_bias_conn1 = cum_bias_conn_1/coef/delta_T

    if n==1:
        wratio = np.ones([box_length, box_width])
        wratio_velocity = np.multiply(wratio,vel_bias_conn1)
        if mask:
            # if average velocity smaller than 1 mm/year (hardcoded here), mask out for better visual
            # this option is only turned on while outputint Wratio.h5 file.
            wratio[abs(vel_bias_conn1)<0.1]=np.nan
    else:
        cum_bias_conn_n =  seq2cum_closure_phase(n, outdir,box)[-1,:,:]
        wratio = np.divide(cum_bias_conn_n,cum_bias_conn_1)
        wratio = 1-wratio
        wratio[wratio>1]=1
        wratio[wratio<0]=0
        wratio_velocity = np.multiply(wratio,vel_bias_conn1)
        if mask:
            wratio[abs(vel_bias_conn1)<0.1]=np.nan # if average velocity smaller than 1 mm/year (hardcoded here), mask out for better visual
    return wratio,wratio_velocity # wratio is a length by width 2D matrix

def estimate_ratioX_all(bw,nl,outdir,box):
    ''' Estimate w(n\delta_t)/w(delta_t) for n=1:bw

    Parameters: nl              - integer, minimum connection-level that we think is bias-free
                bw              - integer, bandwidth of given time-sereis analysis
                box             - list in size of (4,) in integer, coordinates of bounding box
                outdir          - string, the working directory
    Returns:    wratio          - 3D array in size of (bw+1, length, width) in float, the first slice (w[0,:,:]) is a padding to ensure that wratio[n,:,:]=w(n\delta_t)/w(delta_t).
    '''
    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    cum_bias_conn_1 = seq2cum_closure_phase(nl, outdir, box)[-1,:,:]
    wratio = np.zeros([bw+1,box_length, box_width], dtype = np.float32)
    for n in np.arange(2,bw+1):
        cum_bias_conn_n =  seq2cum_closure_phase(n, outdir, box)[-1,:,:]
        wratio[n,:,:] = np.divide(cum_bias_conn_n,cum_bias_conn_1)

    wratio = 1-wratio
    wratio[wratio>1]=1
    wratio[wratio<0]=0
    return wratio

def get_design_matrix_W(M, A, bw, box, tbase, nl, outdir):
    ''' computes the matrix W (Eq. 16 in Zheng et al., 2022) for a bounding box.

    Parameters  : M         - integer, number of interferograms
                  A         - 2D array in size of (M, nslc) in integer, design matrix specifying SAR acquisitions used
                  bw        - integer, bandwidth of time-series analysis
                  tbase     - list in size of (nslc,) in float, time in accumulated years
                  nl        - integer, minimum connection-level that we think is bias-free
                  box       - list in size of (4,) in integer, coordinates of bounding box
                  outdir    - string, the working directory
    Returns:      W         - 2D array in size of (numpix, M) in float, each row stores the diagnal component of W (Eq. 16 in Zheng et al., 2022) for one pixel.
    '''
    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    numpix = box_width * box_length
    # intial output value
    W = np.zeros((numpix,M),dtype = np.float32)
    wratioall = estimate_ratioX_all(bw, nl, outdir, box)
    for i in range(M):
        Aline = list(A[i,:])
        idx1 = Aline.index(-1)
        idx2 = Aline.index(1)
        conn = idx2 - idx1
        if conn > bw:
            print('Interferograms with maximum connection level larger than input bandwidth exists in ifgramStack.h5, use modify_network.py to adjust the maximum connection level')
        wratio = wratioall[conn,:,:]
        wratio = wratio.reshape(-1)
        W[:,i] = wratio

    return W

def average_temporal_span(date_ordinal,bw):
    '''compute average temporal span (days) for interferogram subsets chosen for limited bandwidth analysis

    Parameters:     date_ordinal        - list of size (nslc,) in integer, time in days
                    bw                  - integer, bandwidth of time-series analysis
    Return：        avgtime             - float, average time-span in days for interferograms subsets of bandwith bw.
    '''
    avgtime = 0
    numigram = 0
    for level in range(1, bw+1):
        slcdate_firstn = date_ordinal[0:level]
        slcdate_lastn = date_ordinal[-level:]
        for i in range(level):
            avgtime = avgtime + slcdate_lastn[i] - slcdate_firstn[i]
        numigram = numigram + len(date_ordinal)-level

    avgtime = avgtime/numigram

    return avgtime

def average_connN_igrams(date_ordinal,conn):
    ''' compute average temporal span (days) for connection-n interferograms

    Parameters:     date_ordinal        - list of size (nslc,) in integer, time in days
                    conn                - integer, connection level of interferograms
    Return：        avgtime             - float, average time-span in days for connnection conn interferograms .
    '''
    slcdate_firstn = date_ordinal[0:conn]
    slcdate_lastn = date_ordinal[-conn:]
    avgtime = 0
    for i in range(conn):
        avgtime = avgtime + slcdate_lastn[i] - slcdate_firstn[i]

    numigram = len(date_ordinal)-conn
    avgtime = avgtime/numigram

    return avgtime

def estimate_tsbias_approx(nl, bw, tbase, date_ordinal, wvl, box, outdir):
    ''' This script gives a quick approximate estimate of bias of a time-series of a certain bandwidth (bw) for a bounding box
        This estimate is not exact, but often close enough.
        It is good for a quick estimate to see how big the biases are.

    Parameters: nl              - integer, connection level that we assume bias-free
                bw              - integer, bandwidth of the given time-series analysis
                tbase           - list in size of (nslc,) in float, time in accumulated years
                date_ordinal    - list of size (nslc,) in integer, time in days
                wvl             - float, wavelength of the SAR system
                box             - list in size of (4,) in integer, coordinates of bounding box
                outdir          - string, directory for outputing files
    Returns:    biasts          - 3D array in size of (nslc, box_length, box_width) in float, bias timeseries
    '''
    deltat_n = [average_connN_igrams(date_ordinal,n) for n in range(1,bw+1)] # average temporal span for ifgrams of connection-1 to connection-bw
    avgtimespan = average_temporal_span(date_ordinal,bw)
    p = (np.abs(np.asarray(deltat_n) - avgtimespan)).argmin()+1 # the bias in a bandwidth-bw analysis is similar to bias in connectoin-p interferograms
    print('p = ',p)
    coef = -4*np.pi/wvl
    m1 = 2
    m2 = nl
    wratio_p = estimate_ratioX(tbase, p, nl, wvl, box, outdir)[0]
    wratio_m1 = estimate_ratioX(tbase, m1, nl, wvl, box, outdir)[0]
    wratio_p[np.isnan(wratio_p)] = 0
    wratio_m1[np.isnan(wratio_m1)] = 0
    wratio_m1[abs(wratio_m1-1)<0.1] = np.nan
    ratio1 = np.divide(wratio_p,(1-wratio_m1))
    biasts1 = seq2cum_closure_phase(m1, outdir, box)
    biasts2 = seq2cum_closure_phase(m2, outdir, box)
    for i in range(biasts1.shape[0]):
        biasts1[i,:,:] = np.multiply(biasts1[i,:,:]/coef,ratio1)
        biasts2[i,:,:] = np.multiply(biasts2[i,:,:]/coef,wratio_p)
    biasts = biasts1
    biasts[np.isnan(biasts)]=biasts2[np.isnan(biasts1)]
    return biasts

def quick_bias_correction(ifgram_stack, nl, bw, max_memory, outdir):
    '''Output Wr (eq.20 in Zheng et al., 2022) and a quick approximate solution to bias time-series

    Parameters: ifgram_stack                - string, path for ifgramStack.h5
                nl                          - integer, connection level at which we assume is bias-free
                bw                          - integer, bandwidth of the given time-series.
                wvl                         - float, wavelength of the SAR System
                max_mermory                 - float, maximum memory in GB for each patch processed
                outdir                      - string, directory for output files
    Returns:    Wratio.h5                   - output hdf5 file storing two 3D array of size (bw, length, width) of float, wratios and bias_velocity
                bias_timeseries_approx.h5   - output hdf5 file storing a 3D array of size (nslc, length, width) of float, approximate bias time-series.
    '''
    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)

    date1s = [i.split('_')[0] for i in date12_list]
    date2s = [i.split('_')[1] for i in date12_list]
    slc_list = sorted(list(set(date1s + date2s)))
    # tbase in the unit of years
    date_format = ptime.get_date_str_format(slc_list[0])
    dates = np.array([dt.strptime(i, date_format) for i in slc_list])
    tbase = [i.days + i.seconds / (24 * 60 * 60) for i in (dates - dates[0])]
    tbase = np.array(tbase, dtype=np.float32) / 365.25
    date_ordinal = []
    for date_str in slc_list:
        format_str = '%Y%m%d'
        datetime_obj = dt.strptime(date_str, format_str)
        date_ordinal.append(datetime_obj.toordinal())

    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH']) *100
    slc_list = np.array(slc_list, np.string_)
    connlist = list(np.arange(1,bw+1))
    connlist.append(nl)
    Wr_filedir = os.path.join(outdir, 'Wratio.h5')
    meta['FILE_TYPE'] = None
    ds_name_dict = {'wratio': [np.float32, (len(connlist)-1, length, width), None],
                'bias_velocity': [np.float32, (len(connlist)-1, length, width), None],
                'date': [np.dtype('S8'),np.shape(slc_list), slc_list],}
    writefile.layout_hdf5(Wr_filedir, ds_name_dict, meta)

    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(ifgram_stack, max_memory)

    #process block-by-block
    for i, box in enumerate(box_list):
            box_width  = box[2] - box[0]
            box_length = box[3] - box[1]
            print(box)
            if num_box > 1:
                print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
                print('box width:  {}'.format(box_width))
                print('box length: {}'.format(box_length))

            w_ratios = np.zeros([len(connlist)-1,box_length,box_width])
            w_ratios_velocity = np.zeros([len(connlist)-1,box_length, box_width])
            for idx in range(len(connlist)-1):
                conn = connlist[idx]
                w,wv = estimate_ratioX(tbase, conn, nl, wvl, box, outdir,mask=True)
                w_ratios[idx,:,:] = w
                w_ratios_velocity[idx,:,:] = wv

            # write the block to disk
            block = [0, len(connlist)-1,box[1], box[3], box[0], box[2]]

            writefile.write_hdf5_block(Wr_filedir,
                                       data=w_ratios,
                                       datasetName='wratio',
                                       block=block)

            writefile.write_hdf5_block(Wr_filedir,
                                       data=w_ratios_velocity,
                                       datasetName='bias_velocity',
                                       block=block)

    # a quick/approximate estimate for bias time-series
    biasfile = os.path.join(outdir, 'bias_timeseries_approx.h5')
    meta = dict(stack_obj.metadata)
    ds_name_dict = {'timeseries': [np.float32, (len(slc_list), length, width), None],
            'date': [np.dtype('S8'),np.shape(slc_list), slc_list],}
    writefile.layout_hdf5(biasfile, ds_name_dict, meta)
    for i, box in enumerate(box_list):
        tsbias = estimate_tsbias_approx(nl, bw, tbase, date_ordinal, wvl, box, outdir)
        block = [0, len(slc_list),box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(biasfile,
                                   data=tsbias/100,
                                   datasetName='timeseries',
                                   block=block)

    return

def estimate_bias(ifgram_stack, nl, bw, wvl, box, outdir):
    '''Output bias time-series of a certain bandwidth (bw) for a bounding box using the algorithm provided in Zheng et al., 2022

    Parameters:     ifgram_stack                - string, path for ifgramStack.h5
                    nl                          - integer, connection level at which we assume is bias-free
                    bw                          - integer, bandwidth of the given time-series.
                    wvl                         - float, wavelength of the SAR System
                    box                         - list in size of (4,) in integer, coordinates of bounding box
                    outdir                      - string, directory for output files
    Returns:        biasts_bwn                  - 3D array of size (bw, box_length, box_width) of float, estimated bias time-series
                    box                         - list in size of (4,) in integer, coordinates of bounding box, output for parallel computing
    '''
    coef = -4*np.pi/wvl
    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    numpix = box_width * box_length
    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    A,B = stack_obj.get_design_matrix4timeseries(date12_list = date12_list, refDate = 'no')[0:2]
    B = B[:,:-1]

    # We first need to have the bias time-series for bw-1 analysis
    biasts_bw1_rough = seq2cum_closure_phase(nl, outdir, box)
    m = 2
    biasts_bw1_fine  = seq2cum_closure_phase(m, outdir, box)
    date1s = [i.split('_')[0] for i in date12_list]
    date2s = [i.split('_')[1] for i in date12_list]
    slc_list = sorted(list(set(date1s + date2s)))
    nslc = len(slc_list)
    # tbase in the unit of years
    date_format = ptime.get_date_str_format(slc_list[0])
    dates = np.array([dt.strptime(i, date_format) for i in slc_list])
    tbase = [i.days + i.seconds / (24 * 60 * 60) for i in (dates - dates[0])]
    tbase = np.array(tbase, dtype=np.float32) / 365.25
    tbase_diff = np.diff(tbase).reshape(-1, 1)
    delta_T = tbase[-1]-tbase[0]
    velocity_m = biasts_bw1_fine[-1,:,:]/coef/delta_T
    mask = np.where(np.abs(velocity_m)<0.1, 0,1)

    for i in range(nslc):
        biasts_bw1_fine [i,:,:]  = np.multiply(np.divide(biasts_bw1_rough[-1,:,:],biasts_bw1_fine[-1,:,:]),biasts_bw1_fine[i,:,:])

    biasts_bw1_rough = biasts_bw1_rough.reshape(nslc,-1)
    biasts_bw1_fine = biasts_bw1_fine.reshape(nslc,-1)
    mask = mask.reshape(-1)

    # Then We construct ifgram_bias (W A \Phi^X, or Wr A w(\delta_t)\Phi^X in Eq.(19) in Zheng et al., 2022) , same structure with ifgram_stack
    biasts_bwn = np.zeros((nslc, numpix),dtype = np.float32)
    num_ifgram = np.shape(A)[0]
    if num_ifgram != int(bw*(nslc*2-bw-1)/2): # check the dimensions
        print('Number of interferograms expected: ',int(bw*(nslc*2-bw-1)/2))
        print('Number of interferograms found: ', num_ifgram)
        raise Exception("Modify maximum connection in ifgramStack.h5 to be consistent with input bandwidth!")
    W = get_design_matrix_W(num_ifgram, A, bw, box, tbase, nl, outdir) # this matrix is a numpix by num_ifgram matrix, each row stores the diagnal component of the Wr matrix for that pixel
    prog_bar = ptime.progressBar(maxValue=numpix)
    for i in range(numpix):
        Wr = np.diag(W[i,:])
        WrA = np.matmul(Wr,A)
        Dphi_rough = biasts_bw1_rough[:,i]
        Dphi_fine  = biasts_bw1_fine [:,i]
        if mask[i] == 0 :
            Dphi_bias = np.matmul(WrA,Dphi_rough)
        else:
            Dphi_bias  = np.matmul(WrA,Dphi_fine)
        B_inv  = np.linalg.pinv(B) # here we perform phase velocity inversion as per the original SBAS paper rather doing direct phase inversion.
        biasvel = np.matmul(B_inv,Dphi_bias)
        biasts = np.cumsum(biasvel.reshape(-1)*tbase_diff.reshape(-1))
        biasts_bwn[1:,i] = biasts/coef
        prog_bar.update(i+1, every=200, suffix='{}/{} pixels'.format(i+1, numpix))
    prog_bar.close()
    biasts_bwn = biasts_bwn.reshape(nslc, box_length, box_width)

    return biasts_bwn,box

def bias_correction(ifgram_stack, nl, bw, max_memory, outdir, parallel):
    '''Output a solution to bias time-series

    Parameters: ifgram_stack                - string, path for ifgramStack.h5
                nl                          - integer, connection level at which we assume is bias-free
                bw                          - integer, bandwidth of the given time-series.
                max_memory                 - float, maximum memory in GB for each patch processed
                outdir                      - string, directory for output files
                parallel                    - dictonary containing settings of parallel computing. To turn off, set parallel['clustertype']=''
    Returns:    bias_timeseries.h5          - output hdf5 file storing a 3D array of size (nslc, length, width) of float, estimated bias time-series.
    '''
    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date1s = [i.split('_')[0] for i in date12_list]
    date2s = [i.split('_')[1] for i in date12_list]
    slc_list = sorted(list(set(date1s + date2s)))
    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(ifgram_stack, max_memory)

    # estimate for bias time-series
    biasfile = os.path.join(outdir, 'bias_timeseries.h5')
    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH']) *100 # convert to cm
    slc_list = np.array(slc_list, np.string_)
    nslc = len(slc_list)
    ds_name_dict = {'timeseries': [np.float32, (len(slc_list), length, width), None],
            'date': [np.dtype('S8'),np.shape(slc_list), slc_list],}
    writefile.layout_hdf5(biasfile, ds_name_dict, meta)

    data_kwargs = {
        "ifgram_stack"      : ifgram_stack,
        "nl"                : nl,
        "bw"                : bw,
        "wvl"               : wvl,
        "outdir"            : outdir,
    }
    num_threads_dict = cluster.set_num_threads("1")
    start_time = time.time()
    for i, box in enumerate(box_list):
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        print(box)
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_width))
            print('box length: {}'.format(box_length))
        #update box argument in the input data
        data_kwargs['box'] = box
        if not parallel['clustertype']:
            # non-parallel
            tsbias = estimate_bias(ifgram_stack, nl, bw, wvl, box, outdir)[:-1]
        else:
            # parallel
            print('\n\n------- start parallel processing using Dask -------')
            # initiate the output data
            tsbias = np.zeros((nslc, box_length, box_width), np.float32)
            # initiate dask cluster and client
            cluster_obj = cluster.DaskCluster(parallel['clustertype'], parallel['numWorker'], config_name=parallel['config_name'])
            cluster_obj.open()
            # run dask
            tsbias, box = cluster_obj.run(func=estimate_bias, func_data = data_kwargs, results=[tsbias, box])
            # close dask cluster and client
            cluster_obj.close()
            print('------- finished parallel processing -------\n\n')

        block = [0, len(slc_list),box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(biasfile,
                                   data=tsbias/100,
                                   datasetName='timeseries',
                                   block=block)

    # roll back to the original number of threads
    cluster.roll_back_num_threads(num_threads_dict)
    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return


def create_cp_mask(ifgram_stack, nl, max_memory, num_sigma, threshold_amp, outdir):
    """ create a mask identifying areas most suseptible to bias

    Parameters: ifgram_stack                - string, path for ifgramStack.h5
                nl                          - integer, connection level at which we assume is bias-free
                max_mermory                 - float, maximum memory in GB for each patch processed
                num_sigma                   - float, number of sigmas for computing phase threshold
                threshold_amp               - float, threshold of ampliutde of the cumulative sequential closure phase
                outdir                      - string, directory of output files
    Returns:    maskClosurePhase.h5         - output hdf5 file storing a 2D array of size (length, width) of boolean, areas suseptible to bias is 0.
                avgCpxClosurePhase.h5       - output hdf5 file storing two 2D array of size (length, width) of float, phase and amplitude of average cumulative sequential closure phase of connection-nl
    """
    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date12_list_all = stack_obj.get_date12_list(dropIfgram=False)
    print('scene length, width', length, width)
    ref_phase = stack_obj.get_reference_phase(unwDatasetName = 'unwrapPhase')

    # retrieve the list of SLC dates from ifgramStack.h5
    ifgram0 = date12_list[0]
    date1, date2 = ifgram0.split('_')
    slc_list = [date1, date2]
    for ifgram in date12_list:
        date1, date2 = ifgram.split('_')
        if date1 not in slc_list:
            slc_list.append(date1)
        if date2 not in slc_list:
            slc_list.append(date2)
    slc_list.sort()
    print('number of SLC found : ', len(slc_list))
    print('first SLC: ', slc_list[0])
    print('last  SLC: ', slc_list[-1])

    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(ifgram_stack,max_memory)
    closurephase =  np.zeros([length,width],np.complex64)
    #process block-by-block
    for i, box in enumerate(box_list):
            box_width  = box[2] - box[0]
            box_length = box[3] - box[1]
            print(box)
            if num_box > 1:
                print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
                print('box width:  {}'.format(box_width))
                print('box length: {}'.format(box_length))

            closurephase[box[1]:box[3],box[0]:box[2]], numcp = sum_seq_closure_phase(slc_list, date12_list_all, ifgram_stack, ref_phase,nl,box)
    # What is a good thredshold?
    # Assume that it's pure noise so that the phase is uniform distributed from -pi to pi.
    # The standard deviation of phase in each loop is pi/sqrt(3) (technically should be smaller because when forming loops there should be a reduction in phase variance)
    # The standard deviation of phase in cumulative wrapped closure phase is pi/sqrt(3)/sqrt(numcp) -- again another simplification assuming no correlation.
    # We use 3\delta as threshold -- 99.7% confidence

    threshold_pha = np.pi/np.sqrt(3)/np.sqrt(numcp)*num_sigma

    mask = np.ones([length,width],dtype=bool)
    mask[np.abs(np.angle(closurephase))>threshold_pha] = 0 # this masks areas with potential bias
    mask[np.abs(np.abs(closurephase)/numcp < threshold_amp)] = 1 # this unmasks areas with low correlation (where it's hard to know wheter there is bias either)

    # save mask
    meta = dict(stack_obj.metadata)
    meta['FILE_TYPE'] = 'mask'
    ds_name_dict = {'mask': [np.bool_, (length, width), mask],}
    writefile.layout_hdf5(os.path.join(outdir,'maskClosurePhase.h5'), ds_name_dict, meta)

    # also save the average closure phase
    ds_name_dict2 = {'phase': [np.float32, (length, width), np.angle(closurephase)],
                    'amplitude':[np.float32,(length,width),np.abs(closurephase)/numcp],}
    writefile.layout_hdf5(os.path.join(outdir,'avgCpxClosurePhase.h5'), ds_name_dict2, meta)

    return

def compute_unwrap_closure_phase(ifgram_stack, conn, max_memory, outdir):
    ''' Ouput wrapped, and unwrapped sequential closure phases, and cumulative closure phase time-series of connection-conn in directory outdir/closurePhase/conn{conn}_cp

    Parameters: ifgram_stack                - string, path for ifgramStack.h5
                conn                        - integer, connection level
                max_mermory                 - float, maximum memory in GB for each patch processed
                outdir                      - string, path for output files
    Returns:    various wrapped, unwrapped and cumulative closure phase time-series 

    '''
    stack_obj = ifgramStack(ifgram_stack)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    meta = dict(stack_obj.metadata)
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date12_list_all = stack_obj.get_date12_list(dropIfgram=False)
    print('scene length, width', length, width)
    ref_phase = stack_obj.get_reference_phase(unwDatasetName = 'unwrapPhase')
    refX = stack_obj.refX
    refY = stack_obj.refY
    # retrieve the list of SLC dates from ifgramStack.h5
    ifgram0 = date12_list[0]
    date1, date2 = ifgram0.split('_')
    slc_list = [date1, date2]
    for ifgram in date12_list:
        date1, date2 = ifgram.split('_')
        if date1 not in slc_list:
            slc_list.append(date1)
        if date2 not in slc_list:
            slc_list.append(date2)
    slc_list.sort()
    print('number of SLC found : ', len(slc_list))
    print('first SLC: ', slc_list[0])
    print('last  SLC: ', slc_list[-1])

    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(ifgram_stack,max_memory)

    closurephase =  np.zeros([len(slc_list)-conn, length,width],np.float32)
    #process block-by-block
    for i, box in enumerate(box_list):
            box_width  = box[2] - box[0]
            box_length = box[3] - box[1]
            print(box)
            if num_box > 1:
                print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
                print('box width:  {}'.format(box_width))
                print('box length: {}'.format(box_length))

            closurephase[:,box[1]:box[3],box[0]:box[2]] = seq_closure_phase(slc_list, date12_list_all, ifgram_stack, ref_phase, conn, box)

    # directory
    cpdir = os.path.join(outdir, 'closurePhase')
    if not os.path.isdir(cpdir):
        os.mkdir(cpdir)

    cpdir_conn = os.path.join(cpdir,'conn'+str(conn)+'_cp')
    if not os.path.isdir(cpdir_conn):
        os.mkdir(cpdir_conn)

    # filter and output
    for i in range(len(slc_list)-conn):
        concpname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.int' # some day we will need to make this 4 digits.
        concpdir = os.path.join(cpdir_conn,concpname)
        if not os.path.isfile(concpdir):
            kernel = isce_utils.gaussian_kernel(5,5,1,1)
            closurephase_filt = isce_utils.convolve(np.exp(1j*closurephase[i,:,:]), kernel)
            fid = open(concpdir,mode = 'wb')
            closurephase_filt.tofile(fid)
            fid.close()
            meta['FILE_TYPE']='.int'
            meta['INTERLEAVE']='BIP'
            meta['DATA_TYPE']='complex64'
            meta['BANDS']=1
            writefile.write_isce_xml(meta, concpdir)

    # compute phase sigma and output
    for i in range(len(slc_list)-conn):
        concpname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.int'
        concpdir = os.path.join(cpdir_conn,concpname)
        concpcorname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.cor'
        concpcordir = os.path.join(cpdir_conn, concpcorname)
        if not os.path.isfile(concpcordir):
            isce_utils.estimate_coherence(concpdir, concpcordir)

  #  unwrap
    for i in range(len(slc_list)-conn):
        concpname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.int'
        concpdir = os.path.join(cpdir_conn, concpname)
        concpcorname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.cor'
        concpcordir = os.path.join(cpdir_conn, concpcorname)
        concpunwname = 'conn'+str(conn)+'_filt_'+'{:03}'.format(i)+'.unw'
        concpunwdir = os.path.join(cpdir_conn, concpunwname)
        if not os.path.isfile(concpunwdir):
            isce_utils.unwrap_snaphu(concpdir,concpcordir,concpunwdir,meta)

  # output accumulated unwrapped closure phase time-series
    cum_seq_unw_closure_phase(conn,cpdir_conn,length,width,refY,refX, slc_list, meta)


def main(iargs = None):
    inps = cmd_line_parse(iargs)
    if inps.action == 'create_mask':
        create_cp_mask(inps.ifgram_stack, inps.nl, inps.maxMemory, inps.num_sigma, inps.episilon, inps.outdir)

    if inps.action.endswith ('bias_estimate'):
        maxconn = np.maximum(2,inps.bw) # to make sure we have con-2 closure phase processed
        if inps.update_cp:
            for conn in np.arange(2,maxconn+1):
                compute_unwrap_closure_phase(inps.ifgram_stack, conn, inps.maxMemory, inps.outdir)
            compute_unwrap_closure_phase(inps.ifgram_stack, inps.nl, inps.maxMemory, inps.outdir)
        if inps.action == 'quick_bias_estimate':
            # a quick solution to bias-correction and output diagonal component of Wr (how fast the bias-inducing signal decays with temporal baseline)
            quick_bias_correction(inps.ifgram_stack, inps.nl, inps.bw, inps.maxMemory, inps.outdir)

        if inps.action == 'bias_estimate':
            # bias correction
            parallel={
            "clustertype" : inps.cluster,
            "numWorker"   : inps.numWorker,
            "config_name" : inps.config,
            }
            bias_correction(inps.ifgram_stack, inps.nl, inps.bw, inps.maxMemory, inps.outdir, parallel)

if __name__ == '__main__':
    main(sys.argv[1:])
