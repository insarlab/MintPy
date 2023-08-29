############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2020                           #
############################################################
# Recommend import:
#   from mintpy import solid_earth_tides as SET


import datetime as dt
import os
import time

import h5py
import numpy as np
import pysolid
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 12})

import mintpy.cli.diff
from mintpy.objects import timeseries
from mintpy.objects.resample import resample
from mintpy.utils import ptime, readfile, utils as ut, writefile


###############################################################
def get_datetime_list(ts_file, date_wise_acq_time=False):
    """Prepare exact datetime for each acquisition in the time-series file.

    Parameters: ts_file            - str, path of the time-series HDF5 or .unw file
                date_wise_acq_time - bool, use the exact date-wise acquisition time
    Returns:    sensingMid         - list of datetime.datetime objects
                date_list          - list(str), dates in YYYYMMDD
    """
    print('\nprepare datetime info for each acquisition')

    ts_file = os.path.abspath(ts_file)
    date_list = timeseries(ts_file).get_date_list()

    proj_dir = os.path.dirname(os.path.dirname(ts_file))
    xml_dirs = [os.path.join(proj_dir, i) for i in ['reference', 'secondarys']]

    # list of existing dataset names
    with h5py.File(ts_file, 'r') as f:
        ds_names = [i for i in f.keys() if isinstance(f[i], h5py.Dataset)]

    dt_name = 'sensingMid'
    if dt_name in ds_names:
        # opt 1. read sensingMid if exists
        print(f'read exact datetime info from /{dt_name} in file: {os.path.basename(ts_file)}')
        with h5py.File(ts_file, 'r') as f:
            sensingMidStr = [i.decode('utf-8') for i in f[dt_name][:]]

        # convert string to datetime object
        date_str_format = ptime.get_date_str_format(sensingMidStr[0])
        sensingMid = [dt.datetime.strptime(i, date_str_format) for i in sensingMidStr]

    elif date_wise_acq_time and all(os.path.isdir(i) for i in xml_dirs):
        # opt 2. read sensingMid in xml files [for Sentinel-1 with topsStack]
        print('read exact datetime info in XML files from ISCE-2/topsStack results in directory:', proj_dir)
        from mintpy.utils import isce_utils
        sensingMid = isce_utils.get_sensing_datetime_list(proj_dir, date_list=date_list)[0]

        # plot
        plot_sensingMid_variation(sensingMid)

    elif "T" in date_list[0]:
        # opt 3. use the time info in the `date` dataset [as provided by UAVSAR stack]
        date_format = ptime.get_date_str_format(date_list[0])
        sensingMid = [dt.datetime.strptime(i, date_format) for i in date_list]

    else:
        # opt 4. use constant time of the day for all acquisitions
        atr = readfile.read_attribute(ts_file)
        utc_sec = dt.timedelta(seconds=float(atr['CENTER_LINE_UTC']))
        sensingMid = [dt.datetime.strptime(i, '%Y%m%d') + utc_sec for i in date_list]

        msg =  'Use the same time of the day for all acquisitions from CENTER_LINE_UTC\n'
        if atr.get('PLATFORM', 'Unknow').lower().startswith('sen'):
            msg += 'With <= 1 min variation for Sentinel-1A/B for example, this simplication has negligible impact on SET calculation.'
        print(msg)


    return sensingMid, date_list


def plot_sensingMid_variation(sensingMid, save_fig=True, disp_fig=False, figsize=[8, 3]):
    # calc diff in secs
    dt0 = sensingMid[0]
    sensingMidTime = [i.replace(year=dt0.year, month=dt0.month, day=dt0.day, microsecond=0) for i in sensingMid]

    # plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    ax.plot(sensingMid, sensingMidTime, '.')
    ax.set_ylabel('time of the day\nsensingMid')
    fig.tight_layout()

    # output
    if save_fig:
        out_fig = os.path.abspath('sensingMid_variation.png')
        print('save figure to file', out_fig)
        plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)
    if disp_fig:
        plt.show()
    else:
        plt.close()
    return


###############################################################
def calc_solid_earth_tides_timeseries(ts_file, geom_file, set_comp='enu2los',
        set_file=None, date_wise_acq_time=False, update_mode=True, verbose=False):
    """Calculate the time-series of solid Earth tides (SET) in LOS direction.

    Parameters: ts_file   - str, path of the time-series HDF5 file, or the date list text file
                geom_file - str, path of the geometry HDF5 file
                set_file  - str, output SET time-sereis file
                date_wise_acq_time - bool, use the exact date-wise acquisition time
    Returns:    ts_set    - 3D np.ndarray, SET component time series in meter
    """

    if update_mode and set_file and os.path.isfile(set_file):
        print('update mode: ON')
        print(f'skip re-calculating and use existing file: {set_file}')
        return set_file

    # prepare LOS geometry: geocoding if in radar-coordinates
    los_inc_angle, los_az_angle, atr_geo = ut.prepare_geo_los_geometry(geom_file, unit='deg')

    # get unit vector for the component of interest
    unit_vec = ut.get_unit_vector4component_of_interest(los_inc_angle, los_az_angle, comp=set_comp)
    msg = f'get the unit vector for {set_comp} projection with mean value: '
    msg += f'E = {np.nanmean(unit_vec[0]):.3f}, '
    msg += f'N = {np.nanmean(unit_vec[1]):.3f}, '
    msg += f'U = {np.nanmean(unit_vec[2]):.3f}.'
    print(msg)

    # prepare datetime
    # prefer ts_file to avoid potential missing CENTER_LINE_UTC attributes in geom_file from alosStack
    if ts_file.endswith('.h5'):
        atr = readfile.read_attribute(ts_file)
        dt_objs, date_list = get_datetime_list(ts_file, date_wise_acq_time=date_wise_acq_time)
    else:
        if ts_file.endswith(('.unw', '.geo', '.rdr', '.bip')):
            atr = readfile.read_attribute(ts_file)
            date_list = ptime.yyyymmdd(atr['DATE12'].split('-'))
        else:
            # text list file
            atr = readfile.read_attribute(geom_file)
            date_list = ptime.read_date_txt(ts_file)
        utc_sec = dt.timedelta(seconds=float(atr['CENTER_LINE_UTC']))
        dt_objs = [dt.datetime.strptime(x, '%Y%m%d') + utc_sec for x in date_list]

    # initiate data matrix
    num_date = len(dt_objs)
    length = int(atr_geo['LENGTH'])
    width = int(atr_geo['WIDTH'])
    ts_set = np.zeros((num_date, length, width), dtype=np.float32)
    # default step size in meter: ~30 pixels
    step_size = ut.round_to_1(abs(float(atr_geo['Y_STEP'])) * 108e3 * 30)

    # loop for calc
    print('\n'+'-'*50)
    print('calculating solid Earth tides using PySolid (Milbert, 2018; Yunjun et al., 2022) ...')
    prog_bar = ptime.progressBar(maxValue=num_date, print_msg=not verbose)
    for i, dt_obj in enumerate(dt_objs):
        # calculate tide in ENU direction
        set_e, set_n, set_u = pysolid.calc_solid_earth_tides_grid(
            dt_obj,
            atr_geo,
            step_size=step_size,
            display=False,
            verbose=verbose,
        )

        # convert ENU to LOS direction
        # sign convention: positive for motion towards satellite
        ts_set[i,:,:] = (  set_e * unit_vec[0]
                         + set_n * unit_vec[1]
                         + set_u * unit_vec[2])

        prog_bar.update(i+1, suffix=f'{dt_obj.isoformat()} ({i+1}/{num_date})')
    prog_bar.close()

    # radar-coding if input in radar-coordinates
    if 'Y_FIRST' not in atr.keys():
        print('radar-coding the tides time-series ...')
        res_obj = resample(lut_file=geom_file, max_memory=0)
        res_obj.open()
        res_obj.src_meta = atr_geo
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        ts_set = res_obj.run_resample(
            src_data=ts_set[:, box[1]:box[3], box[0]:box[2]],
        )

    ## output
    if set_file:
        # attribute
        atr['FILE_TYPE'] = 'timeseries'
        atr['DATA_TYPE'] = 'float32'
        atr['UNIT'] = 'm'
        for key in ['REF_Y', 'REF_X', 'REF_DATE']:
            if key in atr.keys():
                atr.pop(key)

        # dataset
        ds_dict = {
            'timeseries' : ts_set,
            'date'       : np.array(date_list, dtype=np.string_),
            'sensingMid' : np.array([i.strftime('%Y%m%dT%H%M%S') for i in dt_objs],
                                    dtype=np.string_),
        }

        # write
        writefile.write(ds_dict, out_file=set_file, metadata=atr)

    return ts_set


###############################################################
def run_solid_earth_tides(inps):
    start_time = time.time()

    # calc SET
    calc_solid_earth_tides_timeseries(
        ts_file=inps.dis_file,
        geom_file=inps.geom_file,
        set_comp=inps.set_comp,
        set_file=inps.set_file,
        date_wise_acq_time=inps.date_wise_acq_time,
        update_mode=inps.update_mode,
        verbose=inps.verbose,
    )

    # correct SET (using diff.py)
    # diff.py can handle different reference in space and time
    # e.g. the absolute phase and the double referenced time-series
    print('correcting tide for using diff.py')
    iargs = [inps.dis_file, inps.set_file, '-o', inps.cor_dis_file, '--force']
    print('diff.py', ' '.join(iargs))
    mintpy.cli.diff.main(iargs)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return
