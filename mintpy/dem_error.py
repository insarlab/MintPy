#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import time
import argparse
import h5py
import numpy as np
from scipy import linalg
from mintpy.objects import timeseries, geometry, cluster
from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_group, ptime, readfile, writefile, utils as ut


# key configuration parameter name
key_prefix = 'mintpy.topographicResidual.'
configKeys = [
    'polyOrder',
    'phaseVelocity',
    'stepFuncDate',
    'excludeDate',
]


############################################################################
TEMPLATE = get_template_content('correct_topography')

EXAMPLE = """example:
  # correct DEM error with pixel-wise geometry parameters [slow]
  dem_error.py  timeseries_ERA5_ramp.h5 -g inputs/geometryRadar.h5 -t smallbaselineApp.cfg

  # correct DEM error with mean geometry parameters [fast]
  dem_error.py  timeseries_ERA5_ramp.h5 -t smallbaselineApp.cfg

  # get updated/corrected DEM
  save_roipac.py inputs/geometryGeo.h5 -o dem.h5   #for dataset in geo coordinates
  mask.py demErr.h5 -m maskTempCoh.h5 -o demErr_msk.h5
  add.py demErr_msk.h5 dem.h5 -o demNew.h5
"""

REFERENCE = """reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
  IEEE TGRS, 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='DEM Error (Topographic Residual) Correction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='{}\n{}\n{}'.format(REFERENCE, TEMPLATE, EXAMPLE))

    parser.add_argument('timeseries_file',
                        help='Timeseries file to be corrrected')
    parser.add_argument('-g', '--geometry', dest='geom_file',
                        help='geometry file including datasets:\n'+
                             'incidence angle\n'+
                             'slant range distance\n' +
                             'and/or 3D perpendicular baseline')
    parser.add_argument('-o', '--outfile',
                        help='Output file name for corrected time-series')

    defo_model = parser.add_argument_group('temporal deformation model')
    defo_model.add_argument('-t', '--template', dest='template_file',
                            help='template file with the options')
    defo_model.add_argument('--ex', '--exclude', dest='excludeDate', nargs='*', default=[],
                            help='Exclude date(s) for DEM error estimation.\n' +
                                 'All dates will be corrected for DEM residual phase still.')
    defo_model.add_argument('-p', '--poly-order', dest='polyOrder', type=int, default=2,
                            help='polynomial order number of temporal deformation model (default: %(default)s).')
    defo_model.add_argument('-s', '--step-date', dest='stepFuncDate', nargs='*', default=[],
                            help='Date of step jump for temporal deformation model (default: %(default)s).'+
                                 ' i.e. date of earthquake/volcanic eruption')
    defo_model.add_argument('--periodic', '--period', '--peri', dest='periodic', type=float, nargs='+', default=[],
                            help='periodic functinos of temporal deformation model (default: %(default)s).')

    parser.add_argument('--phase-velocity', dest='phaseVelocity', action='store_true',
                        help='Use phase velocity instead of phase for inversion constrain.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if:\n'+
                             '1) output timeseries file already exists, readable '+
                             'and newer than input interferograms file\n' +
                             '2) all configuration parameters are the same.')
    # computing
    parser = arg_group.add_memory_argument(parser)
    parser = arg_group.add_parallel_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # --cluster and --num-worker option
    inps.numWorker = str(cluster.DaskCluster.format_num_worker(inps.cluster, inps.numWorker))
    if inps.cluster and inps.numWorker == '1':
        print('WARNING: number of workers is 1, turn OFF parallel processing and continue')
        inps.cluster = None

    # ignore non-existed exclude_date.txt
    if inps.excludeDate == 'exclude_date.txt' and not os.path.isfile(inps.excludeDate):
        inps.excludeDate = []

    if inps.polyOrder < 1:
        raise argparse.ArgumentTypeError("Minimum polynomial order is 1")

    if not inps.outfile:
        inps.outfile = '{}_demErr.h5'.format(os.path.splitext(inps.timeseries_file)[0])
    return inps


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print('1) output file {} NOT found.'.format(inps.outfile))
    else:
        # check if time-series file is partly written using file size
        # since time-series file is not compressed
        with h5py.File(inps.outfile, 'r') as f:
            fsize_ref = f['timeseries'].size * 4
        fsize = os.path.getsize(inps.outfile)
        if fsize <= fsize_ref:
            flag = 'run'
            print('1) output file {} is NOT fully written.'.format(inps.outfile))

        else:
            print('1) output file {} already exists.'.format(inps.outfile))

            # check modification time
            infiles = [inps.timeseries_file]
            if inps.geom_file:
                infiles.append(inps.geom_file)
            ti = max(os.path.getmtime(i) for i in infiles)
            to = os.path.getmtime(inps.outfile)
            if ti > to:
                flag = 'run'
                print('2) output file is NOT newer than input file: {}.'.format(infiles))
            else:
                print('2) output file is newer than input file: {}.'.format(infiles))

    # check configuration
    if flag == 'skip':
        date_list_all = timeseries(inps.timeseries_file).get_date_list()
        inps.excludeDate = read_exclude_date(inps.excludeDate, date_list_all, print_msg=False)[1]
        meta = readfile.read_attribute(inps.outfile)
        if any(str(vars(inps)[key]) != meta.get(key_prefix+key, 'None') for key in configKeys):
            flag = 'run'
            print('3) NOT all key configuration parameters are the same:{}'.format(configKeys))
        else:
            print('3) all key configuration parameters are the same:{}'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


############################################################################
def read_template2inps(template_file, inps=None):
    """Read input template file into inps.excludeDate"""
    if not inps:
        inps = cmd_line_parse()
    iDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(template_file)
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
                value = value.replace('[','').replace(']','').replace(',', ' ')
                iDict[key] = ptime.yyyymmdd(value.split())

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
    msg += "\ntemporal deformation model: polynomial order = {}".format(inps.polyOrder)
    if inps.stepFuncDate:
        msg += "\ntemporal deformation model: step functions at {}".format(inps.stepFuncDate)
    if inps.periodic:
        msg += "\ntemporal deformation model: periodic functions of {} yr".format(inps.periodic)
    msg += '\n'+'-'*80
    print(msg)

    # get design matrix for temporal deformation model
    model = dict()
    model['polynomial'] = inps.polyOrder
    model['step'] = inps.stepFuncDate
    model['periodic'] = inps.periodic
    date_list = timeseries(inps.timeseries_file).get_date_list()
    G_defo = timeseries.get_design_matrix4time_func(date_list, model)

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

    # size
    if box:
        num_row = box[3] - box[1]
        num_col = box[2] - box[0]
    else:
        num_row = ts_obj.length
        num_col = ts_obj.width

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
            print('read 3D bperp from {} file: {} ...'.format(geom_obj.name, os.path.basename(geom_obj.file)))
            dset_list = ['bperp-{}'.format(d) for d in ts_obj.dateList]
            pbase = geom_obj.read(datasetName=dset_list, box=box, print_msg=False).reshape((ts_obj.numDate, -1))
            pbase -= np.tile(pbase[ts_obj.refIndex, :].reshape(1, -1), (ts_obj.numDate, 1))
        else:
            print('read mean bperp from {} file'.format(ts_obj.name))
            pbase = ts_obj.pbase.reshape((-1, 1))

    # 0D geometry
    else:
        print('read mean incidenceAngle, slantRangeDistance, bperp value from {} file'.format(ts_obj.name))
        inc_angle = ut.incidence_angle(ts_obj.metadata, dimension=0)
        range_dist = ut.range_distance(ts_obj.metadata, dimension=0)
        pbase = ts_obj.pbase.reshape((-1, 1))

    sin_inc_angle = np.sin(inc_angle * np.pi / 180.)
    return sin_inc_angle, range_dist, pbase


def estimate_dem_error(ts0, G0, tbase, date_flag=None, phase_velocity=False):
    """Estimate DEM error with least square optimization.
    Parameters: ts0            - 2D np.array in size of (numDate, numPixel), original displacement time-series
                G0             - 2D np.array in size of (numDate, numParam), design matrix in [G_geom, G_defo]
                tbase          - 2D np.array in size of (numDate, 1), temporal baseline
                date_flag      - 1D np.array in bool data type, mark the date used in the estimation
                phase_velocity - bool, use phase history or phase velocity for minimization
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

    # Prepare Design matrix G and observations ts for inversion
    G = G0[date_flag, :]
    ts = ts0[date_flag, :]
    if phase_velocity:
        tbase = tbase[date_flag, :]
        G = np.diff(G, axis=0) / np.diff(tbase, axis=0)
        ts = np.diff(ts, axis=0) / np.diff(tbase, axis=0)

    # Inverse using L-2 norm to get unknown parameters X
    # X = [delta_z, constC, vel, acc, deltaAcc, ..., step1, step2, ...]
    # equivalent to X = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T), ts)
    #               X = np.dot(np.linalg.pinv(G), ts)
    X = linalg.lstsq(G, ts, cond=1e-15)[0]

    # Prepare Outputs
    delta_z = X[0, :]
    ts_cor = ts0 - np.dot(G0[:, 0].reshape(-1, 1), delta_z.reshape(1, -1))
    ts_res = ts0 - np.dot(G0, X)

    # for debug
    debug_mode = False
    if debug_mode:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8))
        ts_all = np.hstack((ts0, ts_res, ts_cor))
        ymin = np.min(ts_all)
        ymax = np.max(ts_all)
        ax1.plot(ts0, '.');           ax1.set_ylim((ymin, ymax)); ax1.set_title('Original  Timeseries')
        ax2.plot(ts_cor, '.');        ax2.set_ylim((ymin, ymax)); ax2.set_title('Corrected Timeseries')
        ax3.plot(ts_res, '.');        ax3.set_ylim((ymin, ymax)); ax3.set_title('Fitting Residual')
        ax4.plot(ts_cor-ts_res, '.'); ax4.set_ylim((ymin, ymax)); ax4.set_title('Fitted Deformation Model')
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
    print(('number of pixels to invert: {} out of {}'
           ' ({:.1f}%)').format(num_pixel2inv,
                                num_pixel,
                                num_pixel2inv/num_pixel*100))


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
        (delta_z_i,
         ts_cor_i,
         ts_res_i) = estimate_dem_error(ts_data[:, mask], G,
                                        tbase=tbase,
                                        date_flag=date_flag,
                                        phase_velocity=phase_velocity)

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
            (delta_z_i,
             ts_cor_i,
             ts_res_i) = estimate_dem_error(ts_data[:, idx], G,
                                            tbase=tbase,
                                            date_flag=date_flag,
                                            phase_velocity=phase_velocity)

            # assemble
            delta_z[idx] = delta_z_i
            ts_cor[:, idx] = ts_cor_i.flatten()
            ts_res[:, idx] = ts_res_i.flatten()

            prog_bar.update(i+1, every=2000, suffix='{}/{}'.format(i+1, num_pixel2inv))
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
    ts_obj = timeseries(inps.timeseries_file)
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
    print('add/update the following configuration metadata to file:\n{}'.format(configKeys))
    for key in configKeys:
        meta[key_prefix+key] = str(vars(inps)[key])

    # 2.2 instantiate est. DEM error
    dem_err_file = 'demErr.h5'
    meta['FILE_TYPE'] = 'dem'
    meta['UNIT'] = 'm'
    ds_name_dict = {'dem' : [np.float32, (length, width), None]}
    writefile.layout_hdf5(dem_err_file, ds_name_dict, metadata=meta)

    # 2.3 instantiate corrected time-series
    ts_cor_file = inps.outfile
    meta['FILE_TYPE'] = 'timeseries'
    writefile.layout_hdf5(ts_cor_file, metadata=meta, ref_file=inps.timeseries_file)

    # 2.4 instantiate residual phase time-series
    ts_res_file = os.path.join(os.path.dirname(inps.outfile), 'timeseriesResidual.h5')
    writefile.layout_hdf5(ts_res_file, metadata=meta, ref_file=inps.timeseries_file)


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
    box_list = cluster.split_box2sub_boxes(box=(0, 0, width, length),
                                           num_split=num_box,
                                           dimension='y')

    # 3.2 prepare the input arguments for *_patch()
    data_kwargs = {
        'G_defo'         : G_defo,
        'ts_file'        : inps.timeseries_file,
        'geom_file'      : inps.geom_file,
        'date_flag'      : date_flag,
        'phase_velocity' : inps.phaseVelocity,
    }

    # 3.3 invert / write block-by-block
    for i, box in enumerate(box_list):
        box_wid = box[2] - box[0]
        box_len = box[3] - box[1]
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_wid))
            print('box length: {}'.format(box_len))

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
            delta_z, ts_cor, ts_res = cluster_obj.run(func=correct_dem_error_patch,
                                                      func_data=data_kwargs,
                                                      results=[delta_z, ts_cor, ts_res])

            # close dask cluster and client
            cluster_obj.close()

            print('------- finished parallel processing -------\n\n')

        # write the block to disk
        # with 3D block in [z0, z1, y0, y1, x0, x1]
        # and  2D block in         [y0, y1, x0, x1]

        # DEM error - 2D
        block = [box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(dem_err_file,
                                   data=delta_z,
                                   datasetName='dem',
                                   block=block)

        # corrected time-series - 3D
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(ts_cor_file,
                                   data=ts_cor,
                                   datasetName='timeseries',
                                   block=block)

        # residual time-series - 3D
        block = [0, num_date, box[1], box[3], box[0], box[2]]
        writefile.write_hdf5_block(ts_res_file,
                                   data=ts_res,
                                   datasetName='timeseries',
                                   block=block)

    # roll back to the origial number of threads
    cluster.roll_back_num_threads(num_threads_dict)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    return dem_err_file, ts_cor_file, ts_res_file


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    # run
    correct_dem_error(inps)

    return inps.outfile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
