############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yujie Zheng, Zhang Yunjun, Feb 2022              #
############################################################
# Recommend import:
#   from mintpy import closure_phase_bias as cpbias


import glob
import os
import time
from datetime import datetime as dt

import numpy as np

from mintpy.ifgram_inversion import estimate_timeseries
from mintpy.objects import cluster, ifgramStack
from mintpy.utils import isce_utils, ptime, readfile, utils as ut, writefile


#################################  Mask  #######################################
def calc_closure_phase_mask(stack_file, bias_free_conn, num_sigma=3, threshold_amp=0.3,
                            outdir='./', max_memory=4.0):
    """Calculate a mask for areas susceptible to biases, based on the average closure phase tau.

    Equation: tau = 1 / K * Sigma_{k=1}^K (np.exp(j * Phi_k^{nl}))
      where K is the number of closure phase for connection nl, Phi_k^{nl} is the k-th sequential
      closure phase for connection nl, as defined in equation (21).
    Reference: Section VI in Zheng et al. (2022, TGRS).

    Parameters: stack_file        - str, path for ifgramStack.h5 file
                bias_free_conn    - int, connection level at which we assume is bias-free
                num_sigma         - float, number of sigmas for computing phase threshold
                threshold_amp     - float, threshold of ampliutde of the cumulative sequential closure phase
                outdir            - str, directory of output files
                max_mermory       - float, maximum memory in GB for each patch processed
    Returns:    mask              - 2D np.ndarray of size (length, width) in boolean, 0 for areas susceptible to biases.
                                    Saved to file: maskClosurePhase.h5
                avg_cp - 2D np.ndarray of size (length, width) in complex64, average cum. seq. closure phase
                                    Saved to file: avgCpxClosurePhase.h5
    """

    # basic info
    stack_obj = ifgramStack(stack_file)
    stack_obj.open(print_msg=False)
    meta = dict(stack_obj.metadata)
    length, width = stack_obj.length, stack_obj.width
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_cp = stack_obj.get_closure_phase_index(bias_free_conn).shape[0]

    ## What is a good thredshold?
    # Assume that it's pure noise so that the phase is uniform distributed from -pi to pi.
    # The standard deviation of phase in each loop is:
    #     pi/sqrt(3)
    # (technically should be smaller because when forming loops there should be a reduction in phase variance)
    # The standard deviation of phase in cumulative wrapped closure phase is:
    #     pi/sqrt(3)/sqrt(num_cp) -- again another simplification assuming no correlation.
    # We use 3\delta as threshold -- 99.7% confidence
    threshold_pha = np.pi / np.sqrt(3 * num_cp) * num_sigma

    # key info
    print('\n'+'-'*80)
    print('calculating the mask to flag areas susceptible to non-closure-phase related biases (as zero) ...')
    print(f'number of valid acquisitions: {len(date_list)} ({date_list[0]} - {date_list[-1]})')
    print(f'average complex closure phase threshold in amplitude/correlation: {threshold_amp}')
    print(f'average complex closure phase threshold in phase: {num_sigma} sigma ({threshold_pha:.1f} rad)')

    # calculate the average complex closure phase
    # process block-by-block to save memory
    print('\ncalculating the average complex closure phase')
    print(f'length / width: {length} / {width}')
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory, dim0_size=stack_obj.numIfgram+num_cp*2)

    avg_cp = np.zeros([length,width], dtype=np.complex64)
    for i, box in enumerate(box_list):
        if num_box > 1:
            print(f'\n------- processing patch {i+1} out of {num_box} --------------')
            print(f'box: {box}')
            print(f'box width:  {box[2] - box[0]}')
            print(f'box length: {box[3] - box[1]}')

        avg_cp[box[1]:box[3], box[0]:box[2]], num_cp = stack_obj.get_sequential_closure_phase(
            box=box,
            conn=bias_free_conn,
            post_proc='mean',
        )[1:]

    # mask out no-data pixels
    geom_file = 'geometryGeo.h5' if 'Y_FIRST' in meta.keys() else 'geometryRadar.h5'
    geom_file = os.path.join(os.path.dirname(stack_file), geom_file)
    if os.path.isfile(geom_file):
        geom_ds_names = readfile.get_dataset_list(geom_file)
        ds_names = [x for x in ['incidenceAngle', 'waterMask'] if x in geom_ds_names]
        if len(ds_names) > 0:
            print(f'mask out pixels with no-data-value (zero {ds_names[0]} from file: {os.path.basename(geom_file)})')
            no_data_mask = readfile.read(geom_file, datasetName=ds_names[0])[0] == 0
            avg_cp[no_data_mask] = np.nan

    # create mask
    print('\ncreate mask for areas susceptible to non-closure phase biases')
    mask = np.ones([length,width], dtype=bool)

    # mask areas with potential bias
    print(f'set pixels with average complex closure phase angle > {num_sigma} sigma ({threshold_pha:.1f} rad) to 0.')
    mask[np.abs(np.angle(avg_cp)) > threshold_pha] = 0

    # unmask areas with low correlation
    # where it's hard to know whether there is bias or not
    print(f'set pixels with average complex closure phase amplitude (correlation) < {threshold_amp} to 1.')
    mask[np.abs(np.abs(avg_cp) < threshold_amp)] = 1

    # write file 1 - mask
    mask_file = os.path.join(outdir, 'maskClosurePhase.h5')
    meta['FILE_TYPE'] = 'mask'
    meta['DATA_TYPE'] = 'bool'
    writefile.write(mask, out_file=mask_file, metadata=meta)

    # write file 2 - average closure phase
    avg_cp_file = os.path.join(outdir, 'avgCpxClosurePhase.h5')
    meta['FILE_TYPE'] = 'mask'
    meta['DATA_TYPE'] = 'float32'
    ds_dict = {
        'amplitude' : np.abs(avg_cp),
        'phase'     : np.angle(avg_cp),
    }
    writefile.write(ds_dict, out_file=avg_cp_file, metadata=meta)

    return mask, avg_cp


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
    '''Output cumulative conn-n sequential closure phase in time-series format,
    which is the weighted phase history of the temporally inconsistent process (f^n / n_l).

    Reference: Equation (25) and (28) in Zheng et al. (2022, TGRS).

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

    # compute cumulative sequential closure phase - f^n
    # equation (25) in Zheng et al. (2022, TGRS)
    num_date = len(date_list)
    bias_ts = np.zeros((num_date, length, width), dtype=np.float32)
    bias_ts[1:num_date-conn+1, :, :] = np.cumsum(cp_phase, 0)
    for i in range(num_date-conn+1, num_date):
        bias_ts[i] = (i - num_date + conn) * cp_phase[-1] + bias_ts[num_date - conn]

    # equation (28)
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
    '''Compute the following phase stack & time-series of connection-conn:

    +   wrapped seq closure phase stack
    + unwrapped seq closure phase stack
    + cumulative unwrapped seq closure phase time-series
    at directory: outdir/closurePhase/conn{conn}

    Parameters: stack_file  - str, path for ifgramStack.h5
                conn        - int, connection level
                max_mermory - float, maximum memory in GB for each patch processed
                outdir      - str, path for output files
    '''
    # output directory
    conn_dir = os.path.join(outdir, f'closurePhase/conn{conn}')
    os.makedirs(conn_dir, exist_ok=True)

    # update mode checking
    cum_cp_file = os.path.join(conn_dir, 'cumSeqClosurePhase.h5')
    if os.path.isfile(cum_cp_file):
        print(f'cumulative unwrapped seq closure phase time-series exists at: {cum_cp_file}, skip re-generating.')
        return

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

    ## default output binary filenames
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
                print(f'\n------- processing patch {i+1} out of {num_box} --------------')
                print(f'box length: {box[3] - box[1]}')
                print(f'box width : {box[2] - box[0]}')

            closure_phase[:, box[1]:box[3], box[0]:box[2]] = stack_obj.get_sequential_closure_phase(
                box=box,
                conn=conn,
            )[0]

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
    print('  Note that a reference point in the ifgramStack.h5 (as attributes "REF_Y/X") is needed to continue. ')
    print('  A good reference point should be a pixel that has good temporal coherence and no bias.')
    cum_seq_unw_closure_phase_timeseries(conn, conn_dir, date_list, meta)

    return


################################################################################
def read_cum_seq_closure_phase4conn(conn, outdir='./', box=None, print_msg=False):
    '''Read cumulative sequential closure phase from individual closure phase directory.

    Reference: Eq. (25) in Zheng et al. (2022).

    Parameters: conn    - integer, connection level of sequential closure phases
                outdir  - string, directory of conn{n}_cumSeqClosurePhase.h5
                box     - list in size of (4,) in integer, coordinates of bounding box
    Returns:    bias_ts - 3D array in size of (num_date, box_lengh, box_wid) in float,
                          cumulative sequential closure phases
    '''
    cum_cp_file = os.path.join(outdir, f'closurePhase/conn{conn}/cumSeqClosurePhase.h5')
    if print_msg:
        print(f'read timeseries from file: {cum_cp_file}')

    bias_ts = readfile.read(cum_cp_file, box=box, print_msg=False)[0]
    return bias_ts


def estimate_wratio(tbase, conn, bias_free_conn, wvl, box, outdir='./', mask=False):
    '''Estimate W_r & velocity bias for the given connection level.

    W_r is the M x M diagonal matrix with W_r(ii) = w(i*delta_t) / w(delta_t), i = 1,2,...,M.
    This is defined in Equation (20), to be used for bias estimation; and can be calculated from
    the weighted phase history (cumulative seq closure phase) based on Equation (29).

    Parameters: tbase          - list(float) or array, in size of (num_date), time in accumulated years
                conn           - integer, connection-level
                bias_free_conn - integer, minimum connection-level that we think is bias-free
                wvl            - float, wavelength
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - str, the working directory
                mask           - bool, whether to mask out areas with average bias velocity less than 1 mm/year
    Returns:    wratio_connN   - 2D np.ndarray of size (box_len, box_wid) in float32, W_r of the input connection level
                vel_bias_connN - 2D np.ndarray of size (box_len, box_wid) in float32, velocity bias of the input connection level
    '''
    delta_t = tbase[-1] - tbase[0]
    phase2range = -1 * wvl / (4 * np.pi)

    # cum/vel_bias at bias free connection level
    cum_bias_connF = read_cum_seq_closure_phase4conn(bias_free_conn, outdir, box)[-1,:,:]
    vel_bias_connF = cum_bias_connF / delta_t * phase2range

    # calc wratio at input connection level
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]
    wratio_connN = np.ones([box_len, box_wid], dtype=np.float32)

    if conn > 1:
        cum_bias_connN = read_cum_seq_closure_phase4conn(conn, outdir, box)[-1,:,:]
        # skip invalid pixels
        flag = np.multiply(~np.isnan(cum_bias_connF), cum_bias_connF != 0)
        # Equation (29)
        wratio_connN[flag] = 1 - cum_bias_connN[flag] / cum_bias_connF[flag]

        # bound within [0, 1]
        wratio_connN[wratio_connN > 1] = 1
        wratio_connN[wratio_connN < 0] = 0

    else:
        cum_bias_connN = None

    # vel_bias at input connection level
    vel_bias_connN = np.multiply(wratio_connN, vel_bias_connF)
    if mask:
        # if average velocity smaller than 1 mm/year (hardcoded here), mask out for better visual
        # this option is only turned on while outputting wratio.h5 file.
        wratio_connN[abs(vel_bias_connF) < 0.001] = np.nan

    # debug mode
    debug_mode = False
    if debug_mode:
        from matplotlib import pyplot as plt

        data_list = [wratio_connN, vel_bias_connN, cum_bias_connN,
                     None,         vel_bias_connF, cum_bias_connF]
        titles = ['w_ratio', f'bias_vel_{conn}', f'bias_{conn}',
                  None,       'bias_vel_F',       'bias_F']

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

    return wratio_connN, vel_bias_connN


def estimate_wratio_all(bw, bias_free_conn, outdir, box):
    '''Estimate diaginal matrix W_r for all connections levels within the given bandwidth.

    Parameters: bias_free_conn - integer, minimum connection-level that we think is bias-free
                bw             - integer, bandwidth of given time-sereis analysis
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, the working directory
    Returns:    wratio         - 3D array in size of (bw+1, length, width) in float32,
                                 the first slice (w[0,:,:]) is a padding
                                 to ensure that wratio[n,:,:] = w(n * delta_t) / w(delta_t).
    '''
    box_wid = box[2] - box[0]
    box_len = box[3] - box[1]
    cum_bias_connF = read_cum_seq_closure_phase4conn(bias_free_conn, outdir, box)[-1,:,:]
    # skip invalid pixels
    flag = np.multiply(~np.isnan(cum_bias_connF), cum_bias_connF != 0)

    wratio = np.ones([bw+1, box_len, box_wid], dtype=np.float32)
    for n in np.arange(2, bw+1):
        cum_bias_connN = read_cum_seq_closure_phase4conn(n, outdir, box)[-1,:,:]
        wratio[n,flag] = 1 - cum_bias_connN[flag] / cum_bias_connF[flag]

    # bound into [0, 1]
    wratio[wratio > 1] = 1
    wratio[wratio < 0] = 0

    return wratio


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


def estimate_bias_timeseries_approx_patch(bias_free_conn, bw, tbase, date_ordinal, wvl, box, outdir):
    '''Quick and approximate estimate of the bias time-series of a certain bandwidth (bw) for a bounding box

    Note: This estimate is not exact, but often close enough.
      It is good for a quick estimate to see how big the biases are.

    Parameters: bias_free_conn - integer, connection level that we assume bias-free
                bw             - integer, bandwidth of the given time-series analysis
                tbase          - 1D np.ndarray in size of (num_date,) in float, time in accumulated years
                date_ordinal   - list of size (num_date,) in integer, time in days
                wvl            - float, wavelength of the SAR system
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, directory for outputting files
    Returns:    bias_ts        - 3D array in size of (num_date, box_len, box_wid) in float, bias timeseries
    '''
    print('\n'+'-'*60)
    print(f'quick and approximate estimation of bias time series for bandwidth = {bw}')
    # basic info
    phase2range = -1 * wvl / (4 * np.pi)
    num_date = tbase.size

    # average temporal span for ifgrams of connection-1 to connection-bw
    deltat_n = np.asarray([get_avg_time_span4conn(date_ordinal, n) for n in range(1, bw+1)])
    avg_time_span = get_avg_time_span_within_bandwidth(date_ordinal, bw)

    # the bias in a bandwidth-bw analysis is similar to bias in connectoin-p interferograms
    p = (np.abs(deltat_n - avg_time_span)).argmin() + 1
    print(f'average connection level within the bandwidth = {p}')

    # get wratio
    kwargs1 = dict(
        bias_free_conn=bias_free_conn,
        wvl=wvl,
        box=box,
        outdir=outdir,
    )
    wratio_p = estimate_wratio(tbase, conn=p, **kwargs1)[0]
    wratio_2 = estimate_wratio(tbase, conn=2, **kwargs1)[0]
    wratio_p[np.isnan(wratio_p)] = 0
    wratio_2[np.isnan(wratio_2)] = 0

    wratio_2[abs(wratio_2 - 1) < 0.1] = np.nan
    ratio_p2 = wratio_p / (1 - wratio_2)

    # get bias_ts
    kwargs2 = dict(outdir=outdir, box=box, print_msg=True)
    bias_ts = read_cum_seq_closure_phase4conn(2, **kwargs2) * phase2range
    bias_ts_bf = read_cum_seq_closure_phase4conn(bias_free_conn, **kwargs2) * phase2range
    for i in range(num_date):
        bias_ts[i,:,:] *= ratio_p2
        bias_ts_bf[i,:,:] *= wratio_p

    flag = np.isnan(bias_ts)
    bias_ts[flag] = bias_ts_bf[flag]

    return bias_ts


def estimate_bias_timeseries_approx(stack_file, bias_free_conn, bw, water_mask_file=None, outdir='./', max_memory=4.0):
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
                                 Shows how fast the bias-inducing signal decays with temporal baseline.
    '''
    print('\n'+'-'*80)
    print(f'quick estimation of the non-closure phase bias time-series for bandwidth={bw} (Zheng et al., 2022) ...')

    stack_obj = ifgramStack(stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH'])

    # time info
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)

    tbase = np.array(ptime.date_list2tbase(date_list)[0], dtype=np.float32) / 365.25
    date_str_fmt = ptime.get_date_str_format(date_list[0])
    date_ordinal = [dt.strptime(x, date_str_fmt).toordinal() for x in date_list]

    # split igram_file into blocks to save memory
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory, dim0_size=num_date)

    # initiate output files
    # 1 - wratio file
    wratio_file = os.path.join(outdir, 'wratio.h5')
    meta['FILE_TYPE'] = 'wratio'
    meta['DATA_TYPE'] = 'float32'
    meta['UNIT'] = '1'
    ds_name_dict = {
        'wratio'       : [np.float32,     (bw, length, width), None],
        'velocityBias' : [np.float32,     (bw, length, width), None],
    }
    ds_unit_dict = {
        'wratio'       : '1',
        'velocityBias' : 'm/year',
    }
    writefile.layout_hdf5(wratio_file, ds_name_dict, metadata=meta, ds_unit_dict=ds_unit_dict)

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
            print(f'\n------- processing patch {i+1} out of {num_box} --------------')
            print(f'box width:  {box_wid}')
            print(f'box length: {box_len}')

        # read water mask
        if water_mask_file:
            print(f'skip pixels on water (zero value from file: {os.path.basename(water_mask_file)})')
            water_mask = readfile.read(water_mask_file, box=box)[0]
        else:
            water_mask_file = None

        # 1 - estimate the wratio(_velocity)
        wratio = np.zeros([bw, box_len, box_wid], dtype=np.float32)
        bias_vel = np.zeros([bw, box_len, box_wid], dtype=np.float32)
        for j in range(bw):
            print(f'estimation W_ratio for connection level: {j+1}')
            wratio[j, :, :], bias_vel[j, :, :] = estimate_wratio(
                tbase,
                conn=j+1,
                bias_free_conn=bias_free_conn,
                wvl=wvl,
                box=box,
                outdir=outdir,
                mask=True)

        if water_mask_file:
            wratio[:, water_mask==0] = np.nan
            bias_vel[:, water_mask==0] = np.nan

        # write the block to disk
        block = [0, bw, box[1], box[3], box[0], box[2]]

        writefile.write_hdf5_block(
            wratio_file,
            data=wratio,
            datasetName='wratio',
            block=block)

        writefile.write_hdf5_block(
            wratio_file,
            data=bias_vel,
            datasetName='velocityBias',
            block=block)

        # 2 - estimate the bias time series
        bias_ts = estimate_bias_timeseries_approx_patch(
            bias_free_conn=bias_free_conn,
            bw=bw,
            tbase=tbase,
            date_ordinal=date_ordinal,
            wvl=wvl,
            box=box,
            outdir=outdir)

        if water_mask_file:
            bias_ts[:, water_mask==0] = np.nan

        # write the block to disk
        block = [0, len(date_list), box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(
            bias_ts_file,
            data=bias_ts,
            datasetName='timeseries',
            block=block)

    return bias_ts_file, wratio_file



################################################################################
def bandwidth2num_ifgram(bw, num_date):
    '''Get the number of interferograms for the network with the given bandwidth.

    Reference: Equation (15) in Zheng et al. (2022)

    Parameters: bw         - int, bandwidth
                num_date   - int, number of acquisitions
    Returns:    num_ifgram - int, number of interferograms
    '''
    return int(bw * (num_date * 2 - bw - 1) / 2)


def get_design_matrix_Wr(date12_list, bw, box, bias_free_conn, outdir='./'):
    '''Computes the matrix W_r for a given bounding box, following Equation (20).

    Parameters: date12_list    - list(str), interfeorgram pairs list in YYYYMMDD_YYYYMMDD
                bw             - integer, bandwidth of time-series analysis
                bias_free_conn - integer, minimum connection-level that we think is bias-free
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, the working directory
    Returns:    Wr             - 2D np.ndarray in size of (num_ifgram, num_pix) in float32,
                                 each row stores the diagnal component of W (Eq. 16 in Zheng et al., 2022) for one pixel.
                A              - 2D np.ndarray in size of (num_ifgram, num_date) in float32
    '''
    # get design matrix A
    A = ifgramStack.get_design_matrix4timeseries(date12_list=date12_list, refDate='no')[0]
    num_ifgram = A.shape[0]

    # get w(delta_t) * phi^x - section VI-A
    wratio_all = estimate_wratio_all(bw, bias_free_conn, outdir, box)

    # initial output value
    num_pix = (box[2] - box[0]) * (box[3] - box[1])
    Wr = np.zeros((num_ifgram, num_pix), dtype=np.float32)
    for i in range(num_ifgram):
        # get the connection level
        Aline = list(A[i,:])
        idx1 = Aline.index(-1)
        idx2 = Aline.index(1)
        conn = idx2 - idx1
        if conn > bw:
            print('Existing max-conn-level > the input bandwidth, '
                  'use modify_network.py inputs/ifgramStack.h5 to adjust the max-conn-level.')

        # assign to Wr matrix
        Wr[i, :] = wratio_all[conn, :, :].reshape(-1)

    return Wr, A


def estimate_bias_timeseries_patch(stack_file, bias_free_conn, bw, wvl, box, water_mask_file=None, outdir='./'):
    '''Estimate the bias time-series of a certain bandwidth (bw) for a bounding box.

    Reference: Zheng et al. (2022, TGRS).

    Parameters: stack_file     - string, path for ifgramStack.h5
                bias_free_conn - integer, connection level at which we assume is bias-free
                bw             - integer, bandwidth of the given time-series.
                wvl            - float, wavelength of the SAR System
                box            - list in size of (4,) in integer, coordinates of bounding box
                outdir         - string, directory for output files
    Returns:    bias_ts_bwn     - 3D array of size (bw, box_len, box_wid) of float, estimated bias time-series
                box            - list in size of (4,) in integer, coordinates of bounding box, output for parallel computing
    '''
    phase2range = -1 * wvl / (4 * np.pi)
    box_wid, box_len = box[2] - box[0], box[3] - box[1]
    num_pix = box_wid * box_len

    stack_obj = ifgramStack(stack_file)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    num_ifgram = len(date12_list)

    # time info
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_date = len(date_list)
    # tbase in the unit of years
    tbase = np.array(ptime.date_list2tbase(date_list)[0], dtype=np.float32) / 365.25
    tbase_diff = np.diff(tbase).reshape(-1,1)


    ## mask of pixels to invert
    mask = np.ones(num_pix, dtype=np.bool_)

    # water mask
    if water_mask_file and os.path.isfile(water_mask_file):
        print(f'skip pixels (on the water) with zero value in file: {os.path.basename(water_mask_file)}')
        water_mask = readfile.read(water_mask_file, box=box)[0].flatten()
        mask *= np.array(water_mask, dtype=np.bool_)
        del water_mask
    else:
        water_mask_file = None

    # invert pixels on mask(s)
    num_pix2inv = int(np.sum(mask))
    idx_pix2inv = np.where(mask)[0]
    print('number of pixels to invert: {} out of {} ({:.1f}%)'.format(
        num_pix2inv, num_pix, num_pix2inv/num_pix*100))


    ## 1. get bias time-series for bw-1 analysis
    kwargs = dict(outdir=outdir, box=box, print_msg=True)
    bias_ts_bw1_rough = read_cum_seq_closure_phase4conn(bias_free_conn, **kwargs).reshape(num_date, -1)
    bias_ts_bw1_fine  = read_cum_seq_closure_phase4conn(2, **kwargs).reshape(num_date, -1)

    bias_vel_bw1 = bias_ts_bw1_fine[-1,:] * phase2range / (tbase[-1] - tbase[0])
    flag = np.where(np.abs(bias_vel_bw1) < 0.001, 0, 1).astype(np.bool_)
    num_pix_less, num_pix_more = np.sum(flag[mask]), np.sum(~flag[mask])
    digit = len(str(num_pix))
    msg = 'number of pixels with bandwidth=1 velocity bias '
    msg += f'< | >= 0.1 mm/yr: {num_pix_less:{digit}d} | {num_pix_more:{digit}d} out of {num_pix} '
    msg += f'({num_pix_less/num_pix*100:.0f}% | {num_pix_less/num_pix*100:.0f}%)'
    print(msg)

    # scale bias_ts_bw1_fine based on bias_ts_bw1_rough
    r2f_flag = np.multiply(~np.isnan(bias_ts_bw1_fine[-1,:]), bias_ts_bw1_fine[-1,:] != 0)
    r2f_scale = np.ones((num_pix), dtype=np.float32)
    r2f_scale[r2f_flag] = bias_ts_bw1_rough[-1,r2f_flag] / bias_ts_bw1_fine[-1,r2f_flag]
    for i in range(num_date):
        bias_ts_bw1_fine[i,:] *= r2f_scale
    del r2f_flag, r2f_scale


    # 2. construct bias_stack = W * A * Phi^X = Wr * A * w(delta_t) * Phi^X
    # Equation (20) in Zheng et al. (2022, TGRS)
    # this matrix is a num_pix by num_ifgram matrix, each row stores the diagnal component of the Wr matrix for that pixel
    print('estimating bias_stack = Wr * A * w(delta_t) * Phi^X (Zheng et al., 2022, TGRS) ...')
    Wr, A = get_design_matrix_Wr(date12_list, bw, box, bias_free_conn, outdir)
    wPhi_x = np.array(bias_ts_bw1_rough, dtype=np.float32)
    wPhi_x[:, flag] = bias_ts_bw1_fine[:, flag]

    bias_stack = np.zeros((num_ifgram, num_pix), dtype=np.float32)
    prog_bar = ptime.progressBar(maxValue=num_pix2inv)
    for i in range(num_pix2inv):
        idx = idx_pix2inv[i]

        # calculate the bias_stack = W * A * phi^x = W^r * A * w(delta_t) * phi^x
        bias_stack[:, idx] = np.linalg.multi_dot([np.diag(Wr[:, idx]), A, wPhi_x[:, idx]]).flatten()

        prog_bar.update(i+1, every=3000, suffix=f'{i+1}/{num_pix2inv} pixels')
    prog_bar.close()
    del bias_ts_bw1_rough, bias_ts_bw1_fine, wPhi_x, Wr


    # 3. estimate bias time-series from bias stack: bias_ts = A+ * bias_stack
    # Equation (20) in Zheng et al. (2022, TGRS)
    # perform phase velocity inversion as per the original SBAS paper rather doing direct phase inversion.
    print('estimating bias time-series from bias stack via the SBAS approach ...')
    A1, B1 = stack_obj.get_design_matrix4timeseries(date12_list=date12_list)
    kwargs = {
        'A'                 : A1,
        'B'                 : B1,
        'tbase_diff'        : tbase_diff,
        'weight_sqrt'       : None,
        'min_norm_velocity' : True,
        'inv_quality_name'  : 'no',
    }

    # a. split mask into mask_all/par_net
    # mask for valid (~NaN) observations in ALL ifgrams (share one B in sbas inversion)
    mask_all_net = np.all(~np.isnan(bias_stack), axis=0) * np.all(bias_stack != 0, axis=0)
    mask_all_net *= mask
    mask_par_net = mask ^ mask_all_net
    msg = 'estimating time-series for pixels with valid stack values'

    # b. invert once for all pixels with obs in ALL ifgrams
    bias_ts = np.zeros((num_date, num_pix), dtype=np.float32)
    if np.sum(mask_all_net) > 0:
        num_pix_all = int(np.sum(mask_all_net))
        print(f'{msg} in all  ifgrams ({num_pix_all} pixels; {num_pix_all/num_pix2inv*100:.0f}%) ...')

        # invert
        bias_ts[:, mask_all_net] = estimate_timeseries(y=bias_stack[:, mask_all_net], **kwargs)[0]

    # c. invert pixel-by-pixel for pixels with obs NOT in all ifgrams
    if np.sum(mask_par_net) > 0:
        num_pix_par = int(np.sum(mask_par_net))
        idx_pix_par = np.where(mask_par_net)[0]
        print(f'{msg} in some ifgrams ({num_pix_par} pixels; {num_pix_par/num_pix2inv*100:.0f}%) ...')

        prog_bar = ptime.progressBar(maxValue=num_pix_par)
        for i in range(num_pix_par):
            idx = idx_pix_par[i]
            # invert
            bias_ts[:, idx] = estimate_timeseries(y=bias_stack[:, idx], **kwargs)[0].flatten()

            prog_bar.update(i+1, every=200, suffix=f'{i+1}/{num_pix_par} pixels')
        prog_bar.close()
    del bias_stack

    bias_ts = bias_ts.reshape(num_date, box_len, box_wid) * phase2range

    return bias_ts, box


def estimate_bias_timeseries(stack_file, bias_free_conn, bw, cluster_kwargs, water_mask_file=None, outdir='./', max_memory=4.0):
    '''Run the bias time-series estimation.

    Parameters: stack_file     - string, path for ifgramStack.h5
                bias_free_conn - integer, connection level at which we assume is bias-free
                bw                - integer, bandwidth of the given time-series.
                cluster_kwargs - dictionary containing settings of parallel computing. To turn off, set parallel['clustertype']=''
                outdir            - string, directory for output files
                max_memory        - float, maximum memory in GB for each patch processed
    Returns:    bias_ts_file   - str, path to the bias time series file: timeseriesBias.h5
    '''
    print('\n'+'-'*80)
    print(f'estimating the non-closure phase bias time-series for bandwidth={bw} (Zheng et al., 2022) ...')

    stack_obj = ifgramStack(stack_file)
    stack_obj.open(print_msg=False)
    length, width = stack_obj.length, stack_obj.width
    meta = dict(stack_obj.metadata)
    wvl = float(meta['WAVELENGTH'])

    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date_list = stack_obj.get_date_list(dropIfgram=True)
    num_ifgram = len(date12_list)
    num_date = len(date_list)

    # check the bandwidth of ifgramStack
    num_ifgram_exp = bandwidth2num_ifgram(bw, num_date)
    print(f'number of interferograms expected for bandwidth={bw}: {num_ifgram_exp}')
    print(f'number of interferograms kept in ifgramStack.h5  : {num_ifgram}')
    if num_ifgram != num_ifgram_exp:
        msg = f'number of the kept interferograms ({num_ifgram}) is NOT the same as expected ({num_ifgram_exp})!'
        msg += f'\n  This indicates the bandwidth between ifgramStack.h5 and the user input ({bw}) are NOT consistent!'
        msg +=  '\n  Modify the network of interferograms to be consistent via modify_network.py by:'
        msg += f'\n  1) set mintpy.network.connNumMax={bw} and 2) re-run modify_network.py -t smallbaselineApp.cfg'
        raise Exception(msg)

    # estimate for bias time-series
    bias_ts_file = os.path.join(outdir, 'timeseriesBias.h5')
    ds_name_dict = {
        'timeseries' : [np.float32,     (num_date, length, width), None],
        'date'       : [np.dtype('S8'), (num_date,),  np.array(date_list, np.string_)],
    }
    meta['FILE_TYPE'] = 'timeseries'
    meta['DATA_TYPE'] = 'float32'
    meta['REF_DATE'] = date_list[0]
    meta['UNIT'] = 'm'
    meta['BANDS'] = '1'
    meta['DATE12'] = f'{date_list[0][2:]}-{date_list[-1][2:]}'
    writefile.layout_hdf5(bias_ts_file, ds_name_dict, meta)

    data_kwargs = {
        "stack_file"      : stack_file,
        "bias_free_conn"  : bias_free_conn,
        "bw"              : bw,
        "wvl"             : wvl,
        "water_mask_file" : water_mask_file,
        "outdir"          : outdir,
    }

    # split igram_file into blocks to save memory
    box_list, num_box = stack_obj.split2boxes(max_memory=max_memory, dim0_size=num_ifgram*2+num_date*3)
    num_threads_dict = cluster.set_num_threads("1")

    for i, box in enumerate(box_list):
        box_wid, box_len = box[2] - box[0], box[3] - box[1]
        print(box)
        if num_box > 1:
            print(f'\n------- processing patch {i+1} out of {num_box} --------------')
            print(f'box width : {box_wid}')
            print(f'box length: {box_len}')

        #update box argument in the input data
        data_kwargs['box'] = box

        if not cluster_kwargs['cluster_type']:
            # non-parallel
            bias_ts = estimate_bias_timeseries_patch(**data_kwargs)[:-1]

        else:
            # parallel
            print('\n\n------- start parallel processing using Dask -------')

            # initiate the output data
            bias_ts = np.zeros((num_date, box_len, box_wid), dtype=np.float32)

            # initiate dask cluster and client
            cluster_obj = cluster.DaskCluster(**cluster_kwargs)
            cluster_obj.open()

            # run dask
            bias_ts = cluster_obj.run(
                func=estimate_bias_timeseries_patch,
                func_data=data_kwargs,
                results=[bias_ts],
            )

            # close dask cluster and client
            cluster_obj.close()
            print('------- finished parallel processing -------\n\n')

        writefile.write_hdf5_block(
            bias_ts_file,
            data=bias_ts,
            datasetName='timeseries',
            block=[0, num_date, box[1], box[3], box[0], box[2]],
        )

    # roll back to the original number of threads
    cluster.roll_back_num_threads(num_threads_dict)

    return bias_ts_file


################################################################################
def run_closure_phase_bias(inps):
    """Run non-closure phase related bias masking or estimation."""
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

    elif inps.action.endswith('estimate'):
        # compute the unwrapped closure phase bias time-series
        # and re-unwrap to mitigate the impact of phase unwrapping errors
        # which can dominate the true non-closure phase.
        # max(2, inps.bw) is used to ensure we have conn-2 closure phase processed
        conn_list = np.arange(2, max(2, inps.bw) + 1).tolist() + [inps.nl]
        conn_list = sorted(list(set(conn_list)))
        for conn in conn_list:
            print('\n'+'-'*80)
            print('calculating the unwrapped closure phase for '
                  f'connection level = {conn} out of {conn_list} ...')
            compute_unwrap_closure_phase(
                stack_file=inps.stack_file,
                conn=conn,
                num_worker=int(inps.numWorker),
                **kwargs)

        if inps.action == 'quick_estimate':
            estimate_bias_timeseries_approx(
                stack_file=inps.stack_file,
                bias_free_conn=inps.nl,
                bw=inps.bw,
                water_mask_file=inps.water_mask_file,
                **kwargs)

        elif inps.action == 'estimate':
            cluster_kwargs = {
                "cluster_type" : inps.cluster,
                "num_worker"   : inps.numWorker,
                "config_name"  : inps.config}
            estimate_bias_timeseries(
                stack_file=inps.stack_file,
                bias_free_conn=inps.nl,
                bw=inps.bw,
                cluster_kwargs=cluster_kwargs,
                water_mask_file=inps.water_mask_file,
                **kwargs)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return
