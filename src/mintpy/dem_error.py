############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import time

import h5py
import numpy as np
from scipy import linalg

from mintpy.objects import cluster, geometry, timeseries
from mintpy.utils import ptime, readfile, time_func, utils as ut, writefile

# key configuration parameter name
key_prefix = 'mintpy.topographicResidual.'
config_keys = [
    'polyOrder',
    'phaseVelocity',
    'stepFuncDate',
    'excludeDate',
]

# debug mode
debug_mode = False


############################################################################
def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file
    if not os.path.isfile(inps.ts_cor_file):
        flag = 'run'
        print(f'1) output file {inps.ts_cor_file} NOT found.')
    else:
        # check if time-series file is partly written using file size
        # since time-series file is not compressed
        with h5py.File(inps.ts_cor_file, 'r') as f:
            fsize_ref = f['timeseries'].size * 4
        fsize = os.path.getsize(inps.ts_cor_file)
        if fsize <= fsize_ref:
            flag = 'run'
            print(f'1) output file {inps.ts_cor_file} is NOT fully written.')

        else:
            print(f'1) output file {inps.ts_cor_file} already exists.')

            # check modification time
            infiles = [inps.ts_file]
            if inps.geom_file:
                infiles.append(inps.geom_file)
            ti = max(os.path.getmtime(i) for i in infiles)
            to = os.path.getmtime(inps.ts_cor_file)
            if ti > to:
                flag = 'run'
                print(f'2) output file is NOT newer than input file: {infiles}.')
            else:
                print(f'2) output file is newer than input file: {infiles}.')

    # check configuration
    if flag == 'skip':
        date_list_all = timeseries(inps.ts_file).get_date_list()
        inps.excludeDate = read_exclude_date(inps.excludeDate, date_list_all, print_msg=False)[1]
        meta = readfile.read_attribute(inps.ts_cor_file)
        if any(str(vars(inps)[key]) != meta.get(key_prefix+key, 'None') for key in config_keys):
            flag = 'run'
            print(f'3) NOT all key configuration parameters are the same:{config_keys}')
        else:
            print(f'3) all key configuration parameters are the same:{config_keys}')

    # result
    print(f'run or skip: {flag}.')
    return flag


############################################################################
def read_template2inps(template_file, inps):
    """Read input template file into inps.excludeDate"""
    iDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    # Read template option
    keyList = [i for i in list(iDict.keys()) if key_prefix+i in template.keys()]
    for key in keyList:
        value = template[key_prefix+key]
        if key in ['phaseVelocity']:
            iDict[key] = value
        elif value:
            if key in ['polyOrder']:
                iDict[key] = int(value)
            elif key in ['excludeDate','stepFuncDate']:
                iDict[key] = ptime.yyyymmdd(value.split(','))

    # computing configurations
    dask_key_prefix = 'mintpy.compute.'
    keyList = [i for i in list(iDict.keys()) if dask_key_prefix+i in template.keys()]
    for key in keyList:
        value = template[dask_key_prefix+key]
        if key in ['cluster', 'config']:
            iDict[key] = value
        elif value:
            if key in ['numWorker']:
                iDict[key] = str(value)
            elif key in ['maxMemory']:
                iDict[key] = float(value)

    return inps


def read_exclude_date(ex_date_list, date_list_all, print_msg=True):
    """Read exclude dates info
    Parameters: ex_date_list  : list of string, date in YYMMDD or YYYYMMDD format,
                                or text file with date in it
                date_list_all : list of string, date in YYYYMMDD format
    Returns:    date_flag     : 1D array of bool in size of (num_date,)
    """
    # Read exclude date input
    ex_date_list = ptime.read_date_list(ex_date_list)
    if ex_date_list and print_msg:
        print(('exclude the following dates for DEM error estimation:'
               ' ({})\n{}').format(len(ex_date_list), ex_date_list))

    # convert to mark array
    date_flag = np.array([i not in ex_date_list for i in date_list_all],
                         dtype=np.bool_)
    return date_flag, ex_date_list


def get_design_matrix4defo(inps):
    """Get the design matrix for ground surface deformation
    Parameters: inps   - namespace
    Returns:    G_defo - 2D np.ndarray in float32 in size of [num_date, num_param]
    """

    # key msg
    msg = '-'*80
    msg += '\ncorrect topographic phase residual (DEM error) (Fattahi & Amelung, 2013, IEEE-TGRS)'
    msg += '\nordinal least squares (OLS) inversion with L2-norm minimization on: phase'
    if inps.phaseVelocity:
        msg += ' velocity'
    msg += f"\ntemporal deformation model: polynomial order = {inps.polyOrder}"
    if inps.stepFuncDate:
        msg += f"\ntemporal deformation model: step functions at {inps.stepFuncDate}"
    if inps.periodic:
        msg += f"\ntemporal deformation model: periodic functions of {inps.periodic} yr"
    msg += '\n'+'-'*80
    print(msg)

    # prepare temporal deformation model
    model = dict()
    model['polynomial'] = inps.polyOrder
    model['step'] = inps.stepFuncDate
    model['periodic'] = inps.periodic

    # prepare SAR info
    ts_obj = timeseries(inps.ts_file)
    date_list = ts_obj.get_date_list()
    seconds = ts_obj.get_metadata().get('CENTER_LINE_UTC', 0)

    # compose design matrix
    G_defo = time_func.get_design_matrix4time_func(date_list, model, seconds=seconds)

    return G_defo


def read_geometry(ts_file, geom_file=None, box=None):
    """Read the following geometry info in 0/2/3D
    Parameters: ts_file       - str, path of time-series file
                geom_file     - str, path of geometry file
                box           - tuple of 4 int for (x0, y0, x1, y1) of the area of interest
    Returns:    sin_inc_angle - 0/2D array, sin(inc_angle)
                range_dist    - 0/2D array, slant range distance in meter
                pbase         - 0/3D array, perp baseline in meter
    """
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)

    # 0/2/3D geometry
    if geom_file:
        geom_obj = geometry(geom_file)
        geom_obj.open()

        # 0/2D incidence angle / slant range distance
        if 'incidenceAngle' not in geom_obj.datasetNames:
            inc_angle = ut.incidence_angle(ts_obj.metadata, dimension=0)
            range_dist = ut.range_distance(ts_obj.metadata, dimension=0)
        else:
            print('read 2D incidenceAngle, slantRangeDistance from {} file: {}'.format(
                geom_obj.name, os.path.basename(geom_obj.file)))
            inc_angle  = geom_obj.read(datasetName='incidenceAngle', box=box, print_msg=False).flatten()
            range_dist = geom_obj.read(datasetName='slantRangeDistance', box=box, print_msg=False).flatten()

        # 0/3D perp baseline
        if 'bperp' in geom_obj.datasetNames:
            print(f'read 3D bperp from {geom_obj.name} file: {os.path.basename(geom_obj.file)} ...')
            dset_list = [f'bperp-{d}' for d in ts_obj.dateList]
            pbase = geom_obj.read(datasetName=dset_list, box=box, print_msg=False).reshape((ts_obj.numDate, -1))
            pbase -= np.tile(pbase[ts_obj.refIndex, :].reshape(1, -1), (ts_obj.numDate, 1))
        else:
            print(f'read mean bperp from {ts_obj.name} file')
            pbase = ts_obj.pbase.reshape((-1, 1))

    # 0D geometry
    else:
        print(f'read mean incidenceAngle, slantRangeDistance, bperp value from {ts_obj.name} file')
        inc_angle = ut.incidence_angle(ts_obj.metadata, dimension=0)
        range_dist = ut.range_distance(ts_obj.metadata, dimension=0)
        pbase = ts_obj.pbase.reshape((-1, 1))

    sin_inc_angle = np.sin(inc_angle * np.pi / 180.)
    return sin_inc_angle, range_dist, pbase


def estimate_dem_error(ts0, G0, tbase, date_flag=None, phase_velocity=False,
                       cond=1e-8, display=False, sharey=True):
    """Estimate DEM error with least square optimization.
    Parameters: ts0            - 2D np.array in size of (numDate, numPixel), original displacement time-series
                G0             - 2D np.array in size of (numDate, numParam), design matrix in [G_geom, G_defo]
                tbase          - 2D np.array in size of (numDate, 1), temporal baseline
                date_flag      - 1D np.array in bool data type, mark the date used in the estimation
                phase_velocity - bool, use phase history or phase velocity for minimization
                cond           - float, cutoff for ‘small’ singular values in least squares solution.
    Returns:    delta_z        - 2D np.array in size of (1,       numPixel) estimated DEM residual
                ts_cor         - 2D np.array in size of (numDate, numPixel),
                                    corrected timeseries = tsOrig - delta_z_phase
                ts_res         - 2D np.array in size of (numDate, numPixel),
                                    residual timeseries = tsOrig - delta_z_phase - defModel
    Example:    delta_z, ts_cor, ts_res = estimate_dem_error(ts, G, tbase, date_flag)
    """
    if len(ts0.shape) == 1:
        ts0 = ts0.reshape(-1, 1)
    if date_flag is None:
        date_flag = np.ones(ts0.shape[0], np.bool_)

    # skip noisy acquisitions
    G = G0[date_flag, :]
    ts = ts0[date_flag, :]

    if phase_velocity:
        # adjust from phase to phase velocity
        tbase_diff = np.diff(tbase[date_flag], axis=0).reshape(-1,1)
        ts = np.diff(ts, axis=0) / np.repeat(tbase_diff, ts.shape[1], axis=1)
        G = np.diff(G, axis=0) / np.repeat(tbase_diff, G.shape[1], axis=1)
        # remove the all-zero column in G and all-one column in G0
        # as the unknown constant term of the polynomial func is not estimated,
        # resulting in a larger time-series residual ts_res, thus, not suitable
        # for RMS based noise evaluation with timeseries_rms.py
        G = np.hstack((G[:,:1], G[:,2:]))
        G0 = np.hstack((G0[:,:1], G0[:,2:]))

    # Inverse using L-2 norm to get unknown parameters X
    # X = [delta_z, constC, vel, acc, deltaAcc, ..., step1, step2, ...]
    # equivalent to X = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T), ts)
    #               X = np.dot(np.linalg.pinv(G), ts)
    X = linalg.lstsq(G, ts, cond=cond)[0]

    # prepare outputs
    delta_z = X[0, :]
    ts_cor = ts0 - np.dot(G0[:, 0].reshape(-1, 1), delta_z.reshape(1, -1))
    ts_res = ts0 - np.dot(G0, X)

    # for debug
    if debug_mode or display:
        from matplotlib import pyplot as plt
        _, axs = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharex=True, sharey=sharey)
        titles = ['Original TS', 'Corrected TS', 'Fitting residual', 'Fitted defo model']
        for ax, data, title in zip(axs, [ts0, ts_cor, ts_res, ts_cor - ts_res], titles):
            ax.plot(data, '.')
            ax.set_title(title)
        plt.show()

    return delta_z, ts_cor, ts_res


def correct_dem_error_patch(G_defo, ts_file, geom_file=None, box=None,
                            date_flag=None, phase_velocity=False):
    """
    Correct one path of a time-series for DEM error.

    Parameters: G_defo         - 2D np.ndarray in float32 in size of (num_date, num_param)
                ts_file        - str, path of time-series file
                geom_file      - str, path of geometry file
                box            - tuple of 4 int in (x0, y0, x1, y1) for the area of interest
                date_flag      - 1D np.ndarray in bool in size of (num_date), dates used for the estimation
                phase_velocity - bool, minimize the resdiual phase or phase velocity
    Returns:    delta_z        - 2D np.ndarray in size of (num_row, num_col)
                ts_cor         - 3D np.ndarray in size of (num_date, num_row, num_col)
                ts_res         - 3D np.ndarray in size of (num_date, num_row, num_col)
                box            - tuple of 4 int in (x0, y0, x1, y1) for the area of interest
    """

    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)

    ## 1. input info
    if debug_mode:
        debug_y, debug_x = 611, 713
        print(f'DEBUG at Y/X = {debug_y}/{debug_x}')
        box = [debug_x, debug_y, debug_x+1, debug_y+1]

    # size
    if box:
        num_row = box[3] - box[1]
        num_col = box[2] - box[0]
    else:
        num_row = ts_obj.length
        num_col = ts_obj.width
    num_pixel = num_row * num_col

    # get date info
    tbase = np.array(ts_obj.tbase, np.float32) / 365.25
    num_date = ts_obj.numDate

    # 1.1 read time-series
    ts_data = readfile.read(ts_file, box=box)[0].reshape(num_date, -1)

    # 1.2 read geometry
    sin_inc_angle, range_dist, pbase = read_geometry(ts_file, geom_file, box=box)

    # 1.3 mask of pixels to invert
    print('skip pixels with ZERO in ALL acquisitions')
    mask = np.nanmean(ts_data, axis=0) != 0.

    print('skip pixels with NaN  in ANY acquisitions')
    mask *= np.sum(np.isnan(ts_data), axis=0) == 0

    tcoh_file = os.path.join(os.path.dirname(ts_file), 'temporalCoherence.h5')
    if os.path.isfile(tcoh_file):
        print('skip pixels with ZERO temporal coherence')
        tcoh = readfile.read(tcoh_file, box=box)[0].flatten()
        mask *= tcoh != 0.
        del tcoh

    if range_dist.size != 1:
        print('skip pixels with ZERO / NaN value in incidenceAngle / slantRangeDistance')
        for geom_data in [sin_inc_angle, range_dist]:
            mask *= geom_data != 0.
            mask *= ~np.isnan(geom_data)

    num_pixel2inv = int(np.sum(mask))
    idx_pixel2inv = np.where(mask)[0]
    perc = num_pixel2inv/num_pixel*100
    print(f'number of pixels to invert: {num_pixel2inv} out of {num_pixel} ({perc:.1f}%)')


    ## 2. estimation

    # 2.1 initiate the output matrices
    delta_z = np.zeros(num_pixel, dtype=np.float32)
    ts_cor = np.zeros((num_date, num_pixel), dtype=np.float32)
    ts_res = np.zeros((num_date, num_pixel), dtype=np.float32)

    # return directly if there is nothing to invert
    if num_pixel2inv < 1:
        delta_z = delta_z.reshape((num_row, num_col))
        ts_cor = ts_cor.reshape((num_date, num_row, num_col))
        ts_res = ts_res.reshape((num_date, num_row, num_col))
        return delta_z, ts_cor, ts_res, box

    # 2.2 estimate
    if range_dist.size == 1:
        print('estimating DEM error ...')
        # compose design matrix
        G_geom = pbase / (range_dist * sin_inc_angle)
        G = np.hstack((G_geom, G_defo))

        # run
        delta_z_i, ts_cor_i, ts_res_i = estimate_dem_error(
            ts0=ts_data[:, mask],
            G0=G,
            tbase=tbase,
            date_flag=date_flag,
            phase_velocity=phase_velocity,
        )

        # assemble
        delta_z[mask] = delta_z_i
        ts_cor[:, mask] = ts_cor_i
        ts_res[:, mask] = ts_res_i

    else:
        print('estimating DEM error pixel-wisely ...')
        prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
        for i in range(num_pixel2inv):
            idx = idx_pixel2inv[i]

            # compose design matrix
            if pbase.shape[1] == 1:
                pbase_i = pbase
            else:
                pbase_i = pbase[:, idx].reshape(-1, 1)

            G_geom = pbase_i / (range_dist[idx] * sin_inc_angle[idx])
            G = np.hstack((G_geom, G_defo))

            # run
            delta_z_i, ts_cor_i, ts_res_i = estimate_dem_error(
                ts0=ts_data[:, idx],
                G0=G,
                tbase=tbase,
                date_flag=date_flag,
                phase_velocity=phase_velocity,
            )

            # assemble
            delta_z[idx] = delta_z_i
            ts_cor[:, idx] = ts_cor_i.flatten()
            ts_res[:, idx] = ts_res_i.flatten()

            prog_bar.update(i+1, every=2000, suffix=f'{i+1}/{num_pixel2inv}')
        prog_bar.close()
    del ts_data, pbase

    ## 3. prepare output
    delta_z = delta_z.reshape((num_row, num_col))
    ts_cor = ts_cor.reshape((num_date, num_row, num_col))
    ts_res = ts_res.reshape((num_date, num_row, num_col))

    return delta_z, ts_cor, ts_res, box


def correct_dem_error(inps):
    """Correct DEM error of input timeseries file"""

    start_time = time.time()

    # limit the number of threads to 1
    # for slight speedup and big CPU usage save
    num_threads_dict = cluster.set_num_threads("1")

    ## 1. input info

    # 1.1 read date info
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    num_date = ts_obj.numDate
    length, width = ts_obj.length, ts_obj.width

    num_step = len(inps.stepFuncDate)

    # exclude dates
    date_flag = read_exclude_date(inps.excludeDate, ts_obj.dateList)[0]
    if inps.polyOrder > np.sum(date_flag):
        raise ValueError("input poly order {} > number of acquisition {}! Reduce it!".format(
            inps.polyOrder, np.sum(date_flag)))

    # 1.2 design matrix part 1 - time func for surface deformation
    G_defo = get_design_matrix4defo(inps)


    ## 2. prepare output

    # 2.1 metadata
    meta = dict(ts_obj.metadata)
    print(f'add/update the following configuration metadata to file:\n{config_keys}')
    for key in config_keys:
        meta[key_prefix+key] = str(vars(inps)[key])

    # 2.2 instantiate est. DEM error
    meta['FILE_TYPE'] = 'dem'
    meta['UNIT'] = 'm'
    ds_name_dict = {'dem' : [np.float32, (length, width), None]}
    writefile.layout_hdf5(inps.dem_err_file, ds_name_dict, metadata=meta)

    # 2.3 instantiate corrected time-series
    meta['FILE_TYPE'] = 'timeseries'
    writefile.layout_hdf5(inps.ts_cor_file, metadata=meta, ref_file=inps.ts_file)

    # 2.4 instantiate residual phase time-series
    ts_res_file = os.path.join(os.path.dirname(inps.ts_cor_file), 'timeseriesResidual.h5')
    writefile.layout_hdf5(ts_res_file, metadata=meta, ref_file=inps.ts_file)


    ## 3. run the estimation and write to disk

    # 3.1 split ts_file into blocks to save memory
    # 1st dimension size: ts (obs / cor / res / step) + dem_err/inc_angle/rg_dist (+pbase)
    num_epoch = num_date * 3 + num_step + 3
    if inps.geom_file:
        geom_obj = geometry(inps.geom_file)
        geom_obj.open(print_msg=False)
        if 'bperp' in geom_obj.datasetNames:
            num_epoch += num_date

    # split in row/line direction based on the input memory limit
    num_box = int(np.ceil((num_epoch * length * width * 4) * 2.5 / (inps.maxMemory * 1024**3)))
    box_list, num_box = cluster.split_box2sub_boxes(
        box=(0, 0, width, length),
        num_split=num_box,
        dimension='y',
    )

    # 3.2 prepare the input arguments for *_patch()
    data_kwargs = {
        'G_defo'         : G_defo,
        'ts_file'        : inps.ts_file,
        'geom_file'      : inps.geom_file,
        'date_flag'      : date_flag,
        'phase_velocity' : inps.phaseVelocity,
    }

    # 3.3 invert / write block-by-block
    for i, box in enumerate(box_list):
        box_wid = box[2] - box[0]
        box_len = box[3] - box[1]
        if num_box > 1:
            print(f'\n------- processing patch {i+1} out of {num_box} --------------')
            print(f'box width:  {box_wid}')
            print(f'box length: {box_len}')

        # update box argument in the input data
        data_kwargs['box'] = box

        # invert
        if not inps.cluster:
            # non-parallel
            delta_z, ts_cor, ts_res = correct_dem_error_patch(**data_kwargs)[:-1]

        else:
            # parallel
            print('\n\n------- start parallel processing using Dask -------')

            # initiate the output data
            delta_z = np.zeros((box_len, box_wid), dtype=np.float32)
            ts_cor = np.zeros((num_date, box_len, box_wid), dtype=np.float32)
            ts_res = np.zeros((num_date, box_len, box_wid), dtype=np.float32)

            # initiate dask cluster and client
            cluster_obj = cluster.DaskCluster(inps.cluster, inps.numWorker, config_name=inps.config)
            cluster_obj.open()

            # run dask
            delta_z, ts_cor, ts_res = cluster_obj.run(
                func=correct_dem_error_patch,
                func_data=data_kwargs,
                results=[delta_z, ts_cor, ts_res],
            )

            # close dask cluster and client
            cluster_obj.close()

            print('------- finished parallel processing -------\n\n')

        # write the block to disk
        # with 3D block in [z0, z1, y0, y1, x0, x1]
        # and  2D block in         [y0, y1, x0, x1]

        # DEM error - 2D
        block = [box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(
            inps.dem_err_file,
            data=delta_z,
            datasetName='dem',
            block=block)

        # corrected time-series - 3D
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(
            inps.ts_cor_file,
            data=ts_cor,
            datasetName='timeseries',
            block=block)

        # residual time-series - 3D
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(
            ts_res_file,
            data=ts_res,
            datasetName='timeseries',
            block=block)

    # roll back to the original number of threads
    cluster.roll_back_num_threads(num_threads_dict)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return inps.dem_err_file, inps.ts_cor_file, ts_res_file
