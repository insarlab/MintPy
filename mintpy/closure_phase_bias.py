#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yujie Zheng, Zhang Yunjun, Feb 2022              #
############################################################
# Recommend import:
#   from mintpy import closure_phase_bias as cpbias


import os
import sys
import time
import argparse
import numpy as np
import glob
from datetime import datetime as dt

from mintpy.objects import ifgramStack, cluster
from mintpy.utils import (
    arg_group,
    ptime,
    readfile,
    writefile,
    isce_utils,
    utils as ut,
)


################################################################################
REFERENCE = """reference:
  Y. Zheng, H. Fattahi, P. Agram, M. Simons and P. Rosen, (2022). On Closure Phase
    and Systematic Bias in Multi-looked SAR Interferometry, in IEEE Trans. Geosci.
    Remote Sens., doi:10.1109/TGRS.2022.3167648.
"""

EXAMPLE = """example:
  # create mask for areas suseptible to biases
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 5  -a mask
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 -a mask --num-sigma 2.5

  # estimate and correct for biases
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 5  --bw 3  -a quick_estimate --num-worker 6
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 10 -a       estimate -c local
"""

def create_parser():
    parser = argparse.ArgumentParser(description = 'Phase non-closure related biases correction',
                                    formatter_class = argparse.RawTextHelpFormatter,
                                    epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('-i','--ifgramstack', type=str, dest='stack_file',
                        help='interferogram stack file that contains the unwrapped phases')
    parser.add_argument('--nl','--conn-level', dest='nl', type=int, default=20,
                        help='connection level that we are correcting to (or consider as no bias)\n'
                             '(default: %(default)s)')
    parser.add_argument('--bw', dest='bw', type=int, default=10,
                        help='bandwidth of time-series analysis that you want to correct')
    parser.add_argument('-a','--action', dest='action', type=str, default='mask',
                        choices={'mask', 'quick_estimate', 'estimate'},
                        help='action to take (default: %(default)s):\n'
                             'mask           - create a mask of areas susceptible to closure phase errors\n'
                             'quick_estimate - quick and approximate estimation on how bias decays with time\n'
                             '                 output sequential closure phase files\n'
                             'estimate       - estimate how bias decays with time\n'
                             '                 processed for each pixel on a pixel by pixel basis [slow]')

    # mask configuration
    mask = parser.add_argument_group('Mask', 'Configuration for closure phase bias mask')
    mask.add_argument('--num-sigma', dest='num_sigma', type=float, default=3,
                      help='Threashold for phase, in number of sigmas (default: %(default)s).\n'
                           'Assuming a Gaussian distribution for the cumulative closure phase'
                           ' with sigma = pi / sqrt(3*num_cp)')
    mask.add_argument('--eps','--epsilon', dest='epsilon', type=float, default=0.3,
                      help='Threashold for the normalized amplitude in [0-1] (default: %(default)s).')

    # compute
    parser = arg_group.add_parallel_argument(parser)
    parser = arg_group.add_memory_argument(parser)

    # output
    parser.add_argument('-o', dest='outdir', type=str, default='./', help='output file directory')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################  Closure Phase  ##################################
def seq_closure_phase(stack_obj, box, conn):
    """Computes wrapped sequential closure phases for a given conneciton level.

    For conn = 5, seq_closure_phase = p12 + p23 + p34 + p45 + p56 - p16.

    Parameters: stack_obj - ifgramStack object
                box       - tuple of 4 int, bounding box in (x0, y0, x1, y1)
                conn      - int, connection level of the closure phase
                normalize - bool, normalize the output complex magnitude by num_cp
    Returns:    cp_w      - 3D np.ndarray in float32 in size of (num_cp, box_len, box_wid)
                            wrapped sequential closure phases for the given connection level.
    """
    # basic info
    num_date = len(stack_obj.get_date_list(dropIfgram=True))
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]

    ## get the closure index
    cp_idx = stack_obj.get_closure_phase_index(conn=conn, dropIfgram=True)
    num_cp = cp_idx.shape[0]
    print(f'number of closure measurements expected: {num_date - conn}')
    print(f'number of closure measurements found   : {num_cp}')
    if num_cp < num_date - conn:
        msg = f'num_cp ({num_cp}) < num_date - conn ({num_date - conn})'
        msg += ' --> some interferograms are missing!'
        raise Exception(msg)

    ## read data
    phase = readfile.read(stack_obj.file, box=box, print_msg=False)[0]
    ref_phase = stack_obj.get_reference_phase(dropIfgram=False)
    for i in range(phase.shape[0]):
        mask = phase[i] != 0.
        phase[i][mask] -= ref_phase[i]

    ## calculate cum seq closure phase
    cp_w = np.zeros((num_cp, box_len, box_wid), dtype=np.float32)
    for i in range(num_cp):

        # calculate closure phase - cp0_w
        idx_plus, idx_minor = cp_idx[i, :-1], cp_idx[i, -1]
        cp0_w = np.sum(phase[idx_plus], axis=0) - phase[idx_minor]

        # get the wrapped closure phase
        cp_w[i] = np.angle(np.exp(1j * cp0_w))

    return cp_w


def sum_seq_closure_phase(stack_obj, box, conn, normalize=False):
    """Computes the sum of complex sequential closure phase for a given connection level.

    For conn = 5, seq_closure_phase = p12 + p23 + p34 + p45 + p56 - p16.

    Parameters: stack_obj  - ifgramStack object
                box        - tuple of 4 int, bounding box in (x0, y0, x1, y1)
                conn       - int, connection level of the closure phase
                normalize  - bool, normalize the output complex magnitude by num_cp
    Returns:    sum_cp     - 2D np.ndarray in complex64 in (box_len, box_wid)
                             sum of sequential closure phase for the given connection level
                num_cp     - integer, number of closure phases used in the sum
    """

    # basic info
    num_date = len(stack_obj.get_date_list(dropIfgram=True))
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]

    ## get the closure index
    cp_idx = stack_obj.get_closure_phase_index(conn=conn, dropIfgram=True)
    num_cp = cp_idx.shape[0]
    print(f'Number of closure measurements expected: {num_date - conn}')
    print(f'Number of closure measurements found   : {num_cp}')
    if num_cp < 1:
        raise Exception(f"No triplets found at connection level: {conn}!")

    ## read unwrapPhase
    phase = readfile.read(stack_obj.file, box=box, print_msg=False)[0]
    ref_phase = stack_obj.get_reference_phase(dropIfgram=False)
    for i in range(phase.shape[0]):
        mask = phase[i] != 0.
        phase[i][mask] -= ref_phase[i]

    ## calculate cum seq closure phase
    sum_cp = np.zeros((box_len, box_wid), dtype=np.complex64)
    for i in range(num_cp):

        # calculate closure phase - cp0_w
        idx_plus, idx_minor = cp_idx[i, :-1], cp_idx[i, -1]
        cp0_w = np.sum(phase[idx_plus], axis=0) - phase[idx_minor]

        # cumulative
        sum_cp += np.exp(1j * cp0_w)

    # normalize
    if normalize:
        sum_cp /= num_cp

    return sum_cp, num_cp



#################################  Mask  #######################################
def calc_closure_phase_mask(stack_file, bias_free_conn, num_sigma=3, threshold_amp=0.3,
                            outdir='./', max_memory=4.0):
    """Calculate a mask for areas most suseptible to biases.

    Parameters: stack_file        - str, path for ifgramStack.h5 file
                bias_free_conn    - int, connection level at which we assume is bias-free
                num_sigma         - float, number of sigmas for computing phase threshold
                threshold_amp     - float, threshold of ampliutde of the cumulative sequential closure phase
                outdir            - str, directory of output files
                max_mermory       - float, maximum memory in GB for each patch processed
    Returns:    mask              - 2D np.ndarray of size (length, width) in boolean, 0 for areas suseptible to biases.
                                    Saved to file: maskClosurePhase.h5
                avg_closure_phase - 2D np.ndarray of size (length, width) in complex64, average cum. seq. closure phase
                                    Saved to file: avgCpxClosurePhase.h5
    """
    print('calculating the mask to flag areas suseptible to non-closure-phase related biases (as zero) ...')

    # basic info
    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date_list = stack_obj.get_date_list(dropIfgram=True)
    print(f'number of valid acquisitions: {len(date_list)}')
    print(f'start / end date: {date_list[0]} / {date_list[-1]}')
    print(f'length / width: {length} / {width}')

    # calculate the average complex closure phase
    # process block-by-block to save memory
    num_cp = stack_obj.get_closure_phase_index(bias_free_conn).shape[0]
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory,
                                              dim0_size=stack_obj.numIfgram+num_cp*2)
    avg_closure_phase = np.zeros([length,width], dtype=np.complex64)
    for i, box in enumerate(box_list):
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print(f'box: {box}')
            print('box width:  {}'.format(box[2] - box[0]))
            print('box length: {}'.format(box[3] - box[1]))

        (avg_closure_phase[box[1]:box[3], box[0]:box[2]],
         num_cp) = sum_seq_closure_phase(
            stack_obj=stack_obj,
            box=box,
            conn=bias_free_conn,
            normalize=True)

    # create mask
    mask = np.ones([length,width], dtype=bool)

    ## What is a good thredshold?
    # Assume that it's pure noise so that the phase is uniform distributed from -pi to pi.
    # The standard deviation of phase in each loop is:
    #     pi/sqrt(3)
    # (technically should be smaller because when forming loops there should be a reduction in phase variance)
    # The standard deviation of phase in cumulative wrapped closure phase is:
    #     pi/sqrt(3)/sqrt(num_cp) -- again another simplification assuming no correlation.
    # We use 3\delta as threshold -- 99.7% confidence
    threshold_pha = np.pi / np.sqrt(3 * num_cp) * num_sigma

    # mask areas with potential bias
    mask[np.abs(np.angle(avg_closure_phase)) > threshold_pha] = 0

    # unmask areas with low correlation
    # where it's hard to know wheter there is bias or not
    mask[np.abs(np.abs(avg_closure_phase) < threshold_amp)] = 1

    # write file 1 - mask
    mask_file = os.path.join(outdir, 'maskClosurePhase.h5')
    meta = dict(stack_obj.metadata)
    meta['FILE_TYPE'] = 'mask'
    meta['DATA_TYPE'] = 'bool'
    writefile.write(mask, out_file=mask_file, metadata=meta)

    # write file 2 - average closure phase
    avg_cp_file = os.path.join(outdir, 'avgCpxClosurePhase.h5')
    meta['FILE_TYPE'] = 'mask'
    meta['DATA_TYPE'] = 'float32'
    ds_dict = {
        'amplitude' : [np.float32, (length, width), np.abs(avg_closure_phase)],
        'phase'     : [np.float32, (length, width), np.angle(avg_closure_phase)],
    }
    writefile.layout_hdf5(avg_cp_file, ds_dict, metadata=meta)

    return mask, avg_closure_phase


################################################################################
def unwrap_closure_phase(int_file, cor_file, unw_file):
    """Unwrap the input wrapped sequential closure phase interferogram.
    """

    if os.path.isfile(cor_file) and os.path.isfile(unw_file):
        print(f'unwrapped interferogram file exists: {unw_file}, skip re-unwrapping.')
    else:
        isce_utils.estimate_coherence(int_file, cor_file)
        isce_utils.unwrap_snaphu(int_file, cor_file, unw_file)

    return unw_file


def cum_seq_unw_closure_phase_timeseries(conn, conn_dir, date_list, meta):
    '''Outpu cumulative conn-n sequential closure phase in time-series format.

    Reference: Eq. 25 in Zheng et al., 2022, but divided by conn.

    Parameters: conn        - int, connection level of closure phases
                conn_dir    - str, path of sequential closure phases file for connection - n
                date_list   - list of str, SLC dates
                meta        - dict, metadata of ifgramStack.h5
    Returns:    bias_ts     - 3D np.ndarray in size of (num_ifgram, length, width) in float32,
                              cumulative sequential closure phase time series,
                              saved to file: cumSeqClosurePhase.h5
                common_mask - 2D np.ndarray in size of (length, width) in bool,
                              mask based on common connected components,
                              saved to file: maskConnComp.h5
    '''

    # output file
    cum_cp_file = os.path.join(conn_dir, 'cumSeqClosurePhase.h5')
    mask_file = os.path.join(conn_dir, 'maskConnComp.h5')

    # update mode checking
    if os.path.isfile(cum_cp_file) and os.path.isfile(mask_file):
        msg = 'cumulative seq closure phase time series and mask files exist, skip re-generating.'
        msg += f'\n{cum_cp_file}\n{mask_file}'
        print(msg)

        # read data
        bias_ts = readfile.read(cum_cp_file)[0]
        common_mask = readfile.read(mask_file)[0]

        return bias_ts, common_mask

    # basic info
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    ref_y, ref_x = int(meta['REF_Y']), int(meta['REF_X'])

    unw_files = sorted(glob.glob(os.path.join(conn_dir, '*.unw')))
    num_file = len(unw_files)

    print('calculate the cumulative seq closure phase time series ...')
    cp_phase = np.zeros((num_file, length, width), dtype=np.float32)
    mask = np.zeros((num_file, length, width), dtype=np.float32)

    prog_bar = ptime.progressBar(maxValue=num_file)
    for i, unw_file in enumerate(unw_files):
        prog_bar.update(i+1, suffix=f'{i+1}/{num_file} {os.path.basename(unw_file)}')

        unw = readfile.read(unw_file, datasetName='phase')[0]
        unw -= unw[ref_y, ref_x]
        cp_phase[i] = unw

        conn_comp_file = unw_file + '.conncomp'
        conn_comp = readfile.read(conn_comp_file)[0]
        mask[i] = np.where(conn_comp >= 1, 1, np.nan)

    prog_bar.close()

    # compute cumulative sequential closure phase
    num_date = len(date_list)
    bias_ts = np.zeros((num_date, length, width), dtype=np.float32)
    bias_ts[1:num_date-conn+1, :, :] = np.cumsum(cp_phase, 0)
    for i in range(num_date-conn+1, num_date):
        bias_ts[i] = (i - num_date + conn) * cp_phase[-1] + bias_ts[num_date - conn]
    bias_ts /= conn

    # write bias time series to HDF5 file
    ds_dict = {
        'timeseries' : [np.float32,     (num_date, length, width), bias_ts],
        'date'       : [np.dtype('S8'), (num_date,), np.array(date_list, np.string_)],
    }
    meta['FILE_TYPE'] = 'timeseries'
    writefile.layout_hdf5(cum_cp_file, ds_dict, metadata=meta)

    # write mask to HDF5 file
    common_mask = np.where(np.isnan(np.sum(mask,0)), False, True)
    meta['FILE_TYPE'] = 'mask'
    writefile.write(common_mask, out_file=mask_file, metadata=meta)

    return bias_ts, common_mask


def compute_unwrap_closure_phase(stack_file, conn, num_worker=1, outdir='./', max_memory=4.0):
    '''Compute the following phase time-series of connection-conn:

    +   wrapped sequential closure phase
    + unwrapped sequential closure phase
    + cumulative unwrapped sequential closure phase
    at directory: outdir/closurePhase/conn{conn}

    Parameters: stack_file  - str, path for ifgramStack.h5
                conn        - int, connection level
                max_mermory - float, maximum memory in GB for each patch processed
                outdir      - str, path for output files
    Returns:    various wrapped, unwrapped and cumulative closure phase time-series

    '''
    print('-'*60)
    print('step 1/3: calculate and filter the wrapped sequential closure phase stack ...')

    # basic info
    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    meta = dict(stack_obj.metadata)
    print(f'scene size: {length} x {width}')

    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)
    print(f'number of acquisitions found: {num_date}')
    print(f'start / end date: {date_list[0]} / {date_list[-1]}')
    # number of expected closure phase
    num_cp = num_date - conn
    num_digit = len(str(num_cp))

    ## default binary filenames
    # output directory
    conn_dir = os.path.join(outdir, f'closurePhase/conn{conn}')
    os.makedirs(conn_dir, exist_ok=True)
    # output file names
    fbases = [os.path.join(conn_dir, f'filt_{x+1:0{num_digit}}') for x in range(num_cp)]
    int_files = [f'{x}.int' for x in fbases]
    cor_files = [f'{x}.cor' for x in fbases]
    unw_files = [f'{x}.unw' for x in fbases]

    if all(os.path.isfile(x) for x in int_files):
        print('ALL the filtered closure phase file exist, skip re-generation.')

    else:
        # process block-by-block
        # split igram_file into blocks to save memory
        box_list, num_box = stack_obj.split2boxes(max_memory=max_memory,
                                                  dim0_size=stack_obj.numIfgram+num_cp*2)
        closure_phase = np.zeros([num_cp, length, width],np.float32)
        for i, box in enumerate(box_list):
            print(box)
            if num_box > 1:
                print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
                print('box length: {}'.format(box[3] - box[1]))
                print('box width : {}'.format(box[2] - box[0]))

            closure_phase[:,
                          box[1]:box[3],
                          box[0]:box[2]] = seq_closure_phase(stack_obj, box=box, conn=conn)

        ## filter the closure phase and re-unwrap
        print('-'*80)
        print('filter the wrapped closure phase stack with a Gaussian kernel of 5 x 5 ...')
        print(f'number of wrapped closure phase: {num_cp}')

        kernel = isce_utils.gaussian_kernel(5, 5, 1, 1)
        for i, int_file in enumerate(int_files):
            if not os.path.isfile(int_file):
                # filter the closure phase interferogram
                closure_phase_filt = isce_utils.convolve(
                    data=np.exp(1j*closure_phase[i]),
                    kernel=kernel,
                ).astype(np.complex64)

                # write to binary file in isce2 format
                print(f'write file: {int_file}')
                with open(int_file, mode='wb') as fid:
                    closure_phase_filt.tofile(fid)

                # write metadata in isce2/roipac format
                meta['FILE_TYPE'] = '.int'
                meta['INTERLEAVE'] = 'BIP'
                meta['DATA_TYPE'] = 'complex64'
                meta['BANDS'] = 1
                writefile.write_isce_xml(meta, int_file)
                writefile.write_roipac_rsc(meta, int_file+'.rsc')
        del closure_phase

    print('-'*60)
    print('step 2/3: unwrap the filtered wrapped closure phase stack ...')
    print(f'number of closure phase: {num_cp}')
    if all(os.path.isfile(x) for x in unw_files):
        print('ALL the unwrapped closure phase file exist, skip re-generation.')

    else:
        num_core, run_parallel, Parallel, delayed = ut.check_parallel(
            num_cp,
            print_msg=False,
            maxParallelNum=num_worker)

        if run_parallel and num_core > 1:
            print(f'parallel processing using {num_core} cores')
            Parallel(n_jobs=num_core)(delayed(unwrap_closure_phase)(x, y, z)
                                      for x, y, z in zip(int_files, cor_files, unw_files))

        else:
            for x, y, z in zip(int_files, cor_files, unw_files):
                unwrap_closure_phase(x, y, z)

    ## calc the cumulativev unwrapped closure phase time-series
    print('-'*60)
    print('step 3/3: calculate the unwrapped cumulative sequential closure phase time-series ...')
    cum_seq_unw_closure_phase_timeseries(conn, conn_dir, date_list, meta)

    return


################################################################################
def read_cum_seq_closure_phase4conn(conn, outdir, box):
    '''Read cumulative sequential closure phase from individual closure phase directory.

    Reference: Eq. (25) in Zheng et al. (2022).

    Parameters: conn    - integer, connection level of sequential closure phases
                outdir  - string, directory of conn{n}_cumSeqClosurePhase.h5
                box     - list in size of (4,) in integer, coordinates of bounding box
    Returns:    biasts  - 3D array in size of (num_date, box_lengh, box_wid) in float,
                          cumulative sequential closure phases
    '''
    seq_cp_file = os.path.join(outdir, f'closurePhase/conn{conn}/cumSeqClosurePhase.h5')
    bias_ts = readfile.read(seq_cp_file, box=box, print_msg=False)[0]
    return bias_ts


def estimate_wratio(tbase, conn, bias_free_conn, wvl, box, outdir='./', mask=False):
    '''Estimate w(n\delta_t)/w(delta_t) - Eq.(29) in Zheng et al., 2022.

    Parameters: tbase          - list(float) or array, in size of (num_date), time in accumulated years
                conn           - integer, connection-level
                bias_free_conn - integer, minimum connection-level that we think is bias-free
                wvl            - float, wavelength
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - str, the working directory
                mask           - bool, whether to mask out areas with average bias velocity less than 1 mm/year
    Returns:    wratio         - 2D array of size (box_len, box_wid) in float, w(n\delta_t)/w(delta_t)
                vel_bias_connN - 2D array of size (box_len, box_wid) in float,
                                 bias-velocity at n*delta_t temporal baseline
    '''
    delta_t = tbase[-1] - tbase[0]
    phase2range = -1 * wvl / (4 * np.pi)

    cum_bias_conn1 = read_cum_seq_closure_phase4conn(bias_free_conn, outdir, box)[-1,:,:]
    flag = np.multiply(~np.isnan(cum_bias_conn1), cum_bias_conn1 != 0)
    vel_bias_conn1 = cum_bias_conn1 / delta_t * phase2range

    # calc wratio
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]
    wratio = np.ones([box_len, box_wid], dtype=np.float32)

    if conn > 1:
        cum_bias_connN = read_cum_seq_closure_phase4conn(conn, outdir, box)[-1,:,:]
        wratio[flag] = 1 - cum_bias_connN[flag] / cum_bias_conn1[flag]

        # bound within [0, 1]
        wratio[wratio > 1] = 1
        wratio[wratio < 0] = 0

    else:
        cum_bias_connN = None

    # vel_bias_connN
    vel_bias_connN = np.multiply(wratio, vel_bias_conn1)
    if mask:
        # if average velocity smaller than 1 mm/year (hardcoded here), mask out for better visual
        # this option is only turned on while outputing Wratio.h5 file.
        wratio[abs(vel_bias_conn1) < 0.1] = np.nan

    # debug mode
    debug_mode = False
    if debug_mode:
        from matplotlib import pyplot as plt

        data_list = [wratio, vel_bias_connN, cum_bias_connN,
                     None,   vel_bias_conn1, cum_bias_conn1]
        titles = ['w_ratio', f'bias_vel_{conn}', f'bias_{conn}',
                  None,       'bias_vel_1',       'bias_1']

        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=[12, 6])
        for ax, data, title in zip(axs.flatten(), data_list, titles):
            if data is not None:
                im = ax.imshow(data, cmap='jet', interpolation='nearest')
                ax.set_title(title)
                fig.colorbar(im, ax=ax)
            else:
                ax.axis('off')
        fig.tight_layout()
        plt.show()

    return wratio, vel_bias_connN


def estimate_wratio_all(bw, bias_free_conn, outdir, box):
    ''' Estimate w(n\delta_t)/w(delta_t) for n=1:bw

    Parameters: bias_free_conn - integer, minimum connection-level that we think is bias-free
                bw             - integer, bandwidth of given time-sereis analysis
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, the working directory
    Returns:    wratio         - 3D array in size of (bw+1, length, width) in float32,
                                 the first slice (w[0,:,:]) is a padding 
                                 to ensure that wratio[n,:,:] = w(n\delta_t)/w(delta_t).
    '''
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]
    cum_bias_conn1 = read_cum_seq_closure_phase4conn(
        bias_free_conn,
        outdir,
        box)[-1,:,:]

    wratio = np.ones([bw+1, box_len, box_wid], dtype=np.float32)
    for n in np.arange(2, bw+1):
        cum_bias_connN = read_cum_seq_closure_phase4conn(n, outdir, box)[-1,:,:]
        wratio[n,:,:] = 1 - cum_bias_connN / cum_bias_conn1

    wratio[wratio > 1] = 1
    wratio[wratio < 0] = 0

    return wratio


def get_design_matrix_W(num_ifgram, A, bw, box, bias_free_conn, outdir):
    '''Computes the matrix W for a given bounding box. (Eq. 16 in Zheng et al., 2022).

    Parameters: num_ifgram     - integer, number of interferograms
                A              - 2D array in size of (num_ifgram, num_date) in int, design matrix specifying SAR acquisitions used
                bw             - integer, bandwidth of time-series analysis
                bias_free_conn - integer, minimum connection-level that we think is bias-free
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, the working directory
    Returns:    W              - 2D array in size of (num_pix, num_ifgram) in float, 
                                 each row stores the diagnal component of W (Eq. 16 in Zheng et al., 2022) for one pixel.
    '''
    wratio_all = estimate_wratio_all(bw, bias_free_conn, outdir, box)

    # intial output value
    num_pix = (box[2] - box[0]) * (box[3] - box[1])
    W = np.zeros((num_pix, num_ifgram),dtype = np.float32)
    for i in range(num_ifgram):
        Aline = list(A[i,:])
        idx1 = Aline.index(-1)
        idx2 = Aline.index(1)
        conn = idx2 - idx1
        if conn > bw:
            print('Existing max-conn-level > the input bandwidth, '
                  'use modify_network.py inputs/ifgramStack.h5 to adjust the max-conn-level.')
        wratio = wratio_all[conn,:,:].reshape(-1)
        W[:,i] = wratio

    return W


def get_avg_time_span_within_bandwidth(date_ordinal, bw):
    '''Compute average temporal span (days) for all interferogram within the given bandwidth

    Parameters: date_ordinal - list of size (num_date,) in integer, time in days
                bw           - integer, bandwidth of time-series analysis
    Return:     avg_time     - float, average time-span in days
    '''
    avg_time = 0
    num_ifgram = 0
    for level in range(1, bw+1):
        slc_date_firstN = date_ordinal[0:level]
        slc_date_lastN  = date_ordinal[-level:]
        for i in range(level):
            avg_time += slc_date_lastN[i] - slc_date_firstN[i]
        num_ifgram += len(date_ordinal) - level

    avg_time /= num_ifgram

    return avg_time


def get_avg_time_span4conn(date_ordinal, conn):
    '''Compute the average temporal span (days) for connection-n interferograms

    Parameters: date_ordinal - list of size (num_date,) in integer, time in days
                conn         - int, connection level of interferograms
    Return:     avg_time     - float, average time-span in days
    '''
    slc_date_firstN = date_ordinal[0:conn]
    slc_date_lastN = date_ordinal[-conn:]
    avg_time = 0
    for i in range(conn):
        avg_time += slc_date_lastN[i] - slc_date_firstN[i]

    num_ifgram = len(date_ordinal) - conn
    avg_time /= num_ifgram

    return avg_time


def estimate_bias_timeseries_approx(bias_free_conn, bw, tbase, date_ordinal, wvl, box, outdir):
    '''Quick and approximate estimate of the bias time-series of a certain bandwidth (bw) for a bounding box

    Note: This estimate is not exact, but often close enough.
      It is good for a quick estimate to see how big the biases are.

    Parameters: bias_free_conn - integer, connection level that we assume bias-free
                bw             - integer, bandwidth of the given time-series analysis
                tbase          - list in size of (num_date,) in float, time in accumulated years
                date_ordinal   - list of size (num_date,) in integer, time in days
                wvl            - float, wavelength of the SAR system
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, directory for outputing files
    Returns:    bias_ts        - 3D array in size of (num_date, box_len, box_wid) in float, bias timeseries
    '''
    print('\n'+'-'*60)
    print(f'quick and approximate estimation of bias time series for bandwidth = {bw}')
    # basic info
    phase2range = -1 * wvl / (4 * np.pi)
    num_date = len(tbase)

    # average temporal span for ifgrams of connection-1 to connection-bw
    deltat_n = np.asarray([get_avg_time_span4conn(date_ordinal, n) for n in range(1, bw+1)])
    avg_time_span = get_avg_time_span_within_bandwidth(date_ordinal, bw)

    # the bias in a bandwidth-bw analysis is similar to bias in connectoin-p interferograms
    p = (np.abs(deltat_n - avg_time_span)).argmin() + 1
    print(f'average connection level within the bandwidth = {p}')

    # get wratio
    kwargs = dict(
        bias_free_conn=bias_free_conn,
        wvl=wvl,
        box=box,
        outdir=outdir,
    )
    wratio_p = estimate_wratio(tbase, conn=p, **kwargs)[0]
    wratio_2 = estimate_wratio(tbase, conn=2, **kwargs)[0]
    wratio_p[np.isnan(wratio_p)] = 0
    wratio_2[np.isnan(wratio_2)] = 0

    wratio_2[abs(wratio_2 - 1) < 0.1] = np.nan
    ratio_p2 = wratio_p / (1 - wratio_2)

    # get bias_ts
    bias_ts = read_cum_seq_closure_phase4conn(2, outdir, box) * phase2range
    bias_ts_bf = read_cum_seq_closure_phase4conn(bias_free_conn, outdir, box) * phase2range
    for i in range(num_date):
        bias_ts[i,:,:] *= ratio_p2
        bias_ts_bf[i,:,:] *= wratio_p

    flag = np.isnan(bias_ts)
    bias_ts[flag] = bias_ts_bf[flag]

    return bias_ts


def quick_bias_estimation(stack_file, bias_free_conn, bw, outdir, max_memory=4.0):
    '''Quick & approximate estimation of the bias time series and Wr.

    Reference: Eq. (20) in Zheng et al. (2022, TGRS).

    Parameters: stack_file     - string, path for ifgramStack.h5
                bias_free_conn - integer, connection level at which we assume is bias-free
                bw             - integer, bandwidth of the given time-series.
                wvl            - float, wavelength of the SAR System
                outdir         - str, directory for output files
                max_mermory    - float, maximum memory in GB for each patch processed
    Returns:    bias_ts_file   - str, path to the HDF5 file for the approximate bias time series in (num_date, length, width)
                wratio_file    - str, path to the HDF5 file for the wratio and bias velocity in (bw, length, width)
    '''
    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH'])

    # time info
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)

    tbase = [x / 365.25 for x in ptime.date_list2tbase(date_list)[0]]
    date_str_fmt = ptime.get_date_str_format(date_list[0])
    date_ordinal = [dt.strptime(x, date_str_fmt).toordinal() for x in date_list]

    # split igram_file into blocks to save memory
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory,
                                              dim0_size=num_date)

    # initiate output files
    # 1 - wratio file
    wratio_file = os.path.join(outdir, 'wratio.h5')
    meta['FILE_TYPE'] = 'wratio'
    meta['DATA_TYPE'] = 'float32'
    meta['UNIT'] = '1'
    ds_name_dict = {
        'wratio'       : [np.float32,     (bw, length, width), None],
        'velocityBias' : [np.float32,     (bw, length, width), None],
        #'date'     : [np.dtype('S8'), (num_date,),  np.array(date_list, np.string_)],
    }
    writefile.layout_hdf5(wratio_file, ds_name_dict, meta)

    # 2 - time series file
    bias_ts_file = os.path.join(outdir, 'timeseriesBiasApprox.h5')
    meta['FILE_TYPE'] = 'timeseries'
    meta['UNIT'] = 'm'
    ds_name_dict = {
        'timeseries' : [np.float32,     (num_date, length, width), None],
        'date'       : [np.dtype('S8'), (num_date,),  np.array(date_list, np.string_)],
    }
    writefile.layout_hdf5(bias_ts_file, ds_name_dict, meta)

    # process block-by-block
    for i, box in enumerate(box_list):
        box_wid = box[2] - box[0]
        box_len = box[3] - box[1]
        print(box)
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_wid))
            print('box length: {}'.format(box_len))

        # 1 - estimate the wratio(_velocity)
        wratio = np.zeros([bw, box_len, box_wid], dtype=np.float32)
        bias_vel = np.zeros([bw, box_len, box_wid], dtype=np.float32)
        for j in range(bw):
            print(f'estimation W_ratio for bandwidth = {j+1}')
            wratio[j, :, :], bias_vel[j, :, :] = estimate_wratio(
                tbase,
                conn=j+1,
                bias_free_conn=bias_free_conn,
                wvl=wvl,
                box=box,
                outdir=outdir,
                mask=True)

        # write the block to disk
        block = [0, bw, box[1], box[3], box[0], box[2]]

        writefile.write_hdf5_block(wratio_file,
                                   data=wratio,
                                   datasetName='wratio',
                                   block=block)

        writefile.write_hdf5_block(wratio_file,
                                   data=bias_vel,
                                   datasetName='velocityBias',
                                   block=block)

        # 2 - estimate the bias time series
        bias_ts = estimate_bias_timeseries_approx(
            bias_free_conn=bias_free_conn,
            bw=bw,
            tbase=tbase,
            date_ordinal=date_ordinal,
            wvl=wvl,
            box=box,
            outdir=outdir)

        # write the block to disk
        block = [0, len(date_list), box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(bias_ts_file,
                                   data=bias_ts,
                                   datasetName='timeseries',
                                   block=block)

    return bias_ts_file, wratio_file



################################################################################
def estimate_bias_timeseries(stack_file, nl, bw, wvl, box, outdir):
    '''Output bias time-series of a certain bandwidth (bw) for a bounding box using the algorithm provided in Zheng et al., 2022

    Parameters: stack_file - string, path for ifgramStack.h5
                nl           - integer, connection level at which we assume is bias-free
                bw           - integer, bandwidth of the given time-series.
                wvl          - float, wavelength of the SAR System
                box          - list in size of (4,) in integer, coordinates of bounding box
                outdir       - string, directory for output files
    Returns:    biasts_bwn   - 3D array of size (bw, box_len, box_wid) of float, estimated bias time-series
                box          - list in size of (4,) in integer, coordinates of bounding box, output for parallel computing
    '''
    phase2range = -1 * wvl / (4 * np.pi)
    box_wid  = box[2] - box[0]
    box_len = box[3] - box[1]
    num_pix = box_wid * box_len
    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    A,B = stack_obj.get_design_matrix4timeseries(date12_list = date12_list, refDate = 'no')[0:2]
    B = B[:,:-1]

    # We first need to have the bias time-series for bw-1 analysis
    biasts_bw1_rough = read_cum_seq_closure_phase4conn(nl, outdir, box)
    m = 2
    biasts_bw1_fine  = read_cum_seq_closure_phase4conn(m, outdir, box)
    date1s = [i.split('_')[0] for i in date12_list]
    date2s = [i.split('_')[1] for i in date12_list]
    date_list = sorted(list(set(date1s + date2s)))
    num_date = len(date_list)
    # tbase in the unit of years
    date_format = ptime.get_date_str_format(date_list[0])
    dates = np.array([dt.strptime(i, date_format) for i in date_list])
    tbase = [i.days + i.seconds / (24 * 60 * 60) for i in (dates - dates[0])]
    tbase = np.array(tbase, dtype=np.float32) / 365.25
    tbase_diff = np.diff(tbase).reshape(-1, 1)
    delta_T = tbase[-1]-tbase[0]
    velocity_m = biasts_bw1_fine[-1,:,:] / delta_T * phase2range
    mask = np.where(np.abs(velocity_m)<0.1, 0,1)

    for i in range(num_date):
        biasts_bw1_fine[i,:,:] = np.multiply(np.divide(biasts_bw1_rough[-1,:,:],
                                                       biasts_bw1_fine[-1,:,:]),
                                             biasts_bw1_fine[i,:,:])

    biasts_bw1_rough = biasts_bw1_rough.reshape(num_date,-1)
    biasts_bw1_fine = biasts_bw1_fine.reshape(num_date,-1)
    mask = mask.reshape(-1)

    # Then We construct ifgram_bias (W A \Phi^X, or Wr A w(\delta_t)\Phi^X
    # Eq.(19) in Zheng et al., 2022), same structure with stack_file
    biasts_bwn = np.zeros((num_date, num_pix),dtype = np.float32)
    num_ifgram = np.shape(A)[0]
    if num_ifgram != int(bw*(num_date*2-bw-1)/2): # check the dimensions
        print('Number of interferograms expected: ',int(bw*(num_date*2-bw-1)/2))
        print('Number of interferograms found: ', num_ifgram)
        raise Exception("Modify maximum connection in ifgramStack.h5 to be consistent with input bandwidth!")

    # this matrix is a num_pix by num_ifgram matrix, each row stores the diagnal component of the Wr matrix for that pixel
    W = get_design_matrix_W(num_ifgram, A, bw, box, nl, outdir)
    prog_bar = ptime.progressBar(maxValue=num_pix)
    for i in range(num_pix):
        Wr = np.diag(W[i,:])
        WrA = np.matmul(Wr,A)
        Dphi_rough = biasts_bw1_rough[:,i]
        Dphi_fine  = biasts_bw1_fine [:,i]
        if mask[i] == 0 :
            Dphi_bias = np.matmul(WrA,Dphi_rough)
        else:
            Dphi_bias  = np.matmul(WrA,Dphi_fine)

        # here we perform phase velocity inversion as per the original SBAS paper rather doing direct phase inversion.
        B_inv  = np.linalg.pinv(B)
        biasvel = np.matmul(B_inv,Dphi_bias)
        biasts = np.cumsum(biasvel.reshape(-1)*tbase_diff.reshape(-1))
        biasts_bwn[1:,i] = biasts * phase2range
        prog_bar.update(i+1, every=200, suffix='{}/{} pixels'.format(i+1, num_pix))
    prog_bar.close()
    biasts_bwn = biasts_bwn.reshape(num_date, box_len, box_wid)

    return biasts_bwn,box


def bias_estimation(stack_file, nl, bw, parallel, outdir, max_memory=4.0):
    '''Output a solution to bias time-series

    Parameters: stack_file       - string, path for ifgramStack.h5
                nl                 - integer, connection level at which we assume is bias-free
                bw                 - integer, bandwidth of the given time-series.
                max_memory         - float, maximum memory in GB for each patch processed
                outdir             - string, directory for output files
                parallel           - dictonary containing settings of parallel computing. To turn off, set parallel['clustertype']=''
    Returns:    bias_timeseries.h5 - output hdf5 file storing a 3D array of size (num_date, length, width) of float, estimated bias time-series.
    '''
    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date1s = [i.split('_')[0] for i in date12_list]
    date2s = [i.split('_')[1] for i in date12_list]
    date_list = sorted(list(set(date1s + date2s)))
    # split igram_file into blocks to save memory
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory,
                                              dim0_size=stack_obj.numIfgram*2)

    # estimate for bias time-series
    biasfile = os.path.join(outdir, 'bias_timeseries.h5')
    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH']) *100 # convert to cm
    date_list = np.array(date_list, np.string_)
    num_date = len(date_list)
    ds_name_dict = {
        'timeseries' : [np.float32,     (len(date_list), length, width), None],
        'date'       : [np.dtype('S8'), np.shape(date_list),             date_list],}
    writefile.layout_hdf5(biasfile, ds_name_dict, meta)

    data_kwargs = {
        "stack_file" : stack_file,
        "nl"         : nl,
        "bw"         : bw,
        "wvl"        : wvl,
        "outdir"     : outdir,
    }
    num_threads_dict = cluster.set_num_threads("1")
    start_time = time.time()
    for i, box in enumerate(box_list):
        box_wid  = box[2] - box[0]
        box_len = box[3] - box[1]
        print(box)
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_wid))
            print('box length: {}'.format(box_len))
        #update box argument in the input data
        data_kwargs['box'] = box
        if not parallel['clustertype']:
            # non-parallel
            tsbias = estimate_bias_timeseries(stack_file, nl, bw, wvl, box, outdir)[:-1]
        else:
            # parallel
            print('\n\n------- start parallel processing using Dask -------')
            # initiate the output data
            tsbias = np.zeros((num_date, box_len, box_wid), np.float32)
            # initiate dask cluster and client
            cluster_obj = cluster.DaskCluster(
                cluster_type=parallel['clustertype'],
                num_worker=parallel['numWorker'],
                config_name=parallel['config_name'],
            )
            cluster_obj.open()
            # run dask
            tsbias, box = cluster_obj.run(
                func=estimate_bias_timeseries,
                func_data=data_kwargs,
                results=[tsbias, box],
            )
            # close dask cluster and client
            cluster_obj.close()
            print('------- finished parallel processing -------\n\n')

        block = [0, len(date_list),box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(biasfile,
                                   data=tsbias/100,
                                   datasetName='timeseries',
                                   block=block)

    # roll back to the original number of threads
    cluster.roll_back_num_threads(num_threads_dict)
    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return






################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # common inputs
    kwargs = dict(outdir=inps.outdir, max_memory=inps.maxMemory)    

    if inps.action == 'mask':
        calc_closure_phase_mask(
            stack_file=inps.stack_file,
            bias_free_conn=inps.nl,
            num_sigma=inps.num_sigma,
            threshold_amp=inps.epsilon,
            **kwargs)

    elif inps.action.endswith ('estimate'):
        # compute the unwrapped closure phase bias time-series
        # to make sure we have conn-2 closure phase processed
        max_conn = np.maximum(2, inps.bw)
        conn_list = np.arange(2, max_conn + 1).tolist() + [inps.nl]
        conn_list = sorted(list(set(conn_list)))
        for conn in conn_list:
            print('\n'+'-'*80)
            print(f'calculating the unwrapped closure phase for connection level = {conn} out of {conn_list} ...')
            compute_unwrap_closure_phase(
                stack_file=inps.stack_file,
                conn=conn,
                num_worker=int(inps.numWorker),
                **kwargs)

        if inps.action == 'quick_estimate':
            # a quick solution to bias-correction
            # output diagonal component of Wr (how fast the bias-inducing signal decays with temporal baseline)
            print('\n'+'-'*80)
            print('quick estimation of the bias time-series ...')
            quick_bias_estimation(
                stack_file=inps.stack_file,
                bias_free_conn=inps.nl,
                bw=inps.bw,
                **kwargs)

        elif inps.action == 'estimate':
            # bias correction
            parallel={
                "clustertype" : inps.cluster,
                "numWorker"   : inps.numWorker,
                "config_name" : inps.config,
            }
            bias_estimation(inps.stack_file, inps.nl, inps.bw, parallel, **kwargs)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
