############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2020                           #
############################################################
# Recomend import:
#   from mintpy import solid_earth_tides as SET


import os
import datetime as dt
import time
import warnings

import h5py
import numpy as np
import pysolid
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})

from mintpy.objects import timeseries
from mintpy.objects.resample import resample
from mintpy.utils import (
    ptime,
    readfile,
    writefile,
    utils as ut,
)


###############################################################
def get_datetime_list(ts_file, date_wise_acq_time=False):
    """Prepare exact datetime for each acquisition in the time-series file.

    Parameters: ts_file            - str, path of the time-series HDF5 file
                date_wise_acq_time - bool, use the exact date-wise acquisition time
    Returns:    sensingMid         - list of datetime.datetime objects
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
        print('read exact datetime info from /{} in file: {}'.format(dt_name, os.path.basename(ts_file)))
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


    return sensingMid


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
def calc_solid_earth_tides_timeseries(ts_file, geom_file, set_file, date_wise_acq_time=False,
                                      update_mode=True, verbose=False):
    """Calculate the time-series of solid Earth tides (SET) in LOS direction.
    Parameters: ts_file   - str, path of the time-series HDF5 file
                geom_file - str, path of the geometry HDF5 file
                set_file  - str, output SET time-sereis file
                date_wise_acq_time - bool, use the exact date-wise acquisition time
    Returns:    set_file  - str, output SET time-sereis file
    """

    if update_mode and os.path.isfile(set_file):
        print('update mode: ON')
        print('skip re-calculating and use existing file: {}'.format(set_file))
        return set_file

    # prepare LOS geometry: geocoding if in radar-coordinates
    inc_angle, head_angle, atr_geo = ut.prepare_geo_los_geometry(geom_file, unit='rad')

    # get LOS unit vector
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        unit_vec = [
            np.sin(inc_angle) * np.cos(head_angle) * -1,
            np.sin(inc_angle) * np.sin(head_angle),
            np.cos(inc_angle),
        ]

    # prepare datetime
    dt_objs = get_datetime_list(ts_file, date_wise_acq_time=date_wise_acq_time)

    # initiate data matrix
    num_date = len(dt_objs)
    length = int(atr_geo['LENGTH'])
    width = int(atr_geo['WIDTH'])
    ts_tide = np.zeros((num_date, length, width), dtype=np.float32)
    # default step size in meter: ~30 pixels
    step_size = ut.round_to_1(abs(float(atr_geo['Y_STEP'])) * 108e3 * 30)

    # loop for calc
    print('\n'+'-'*50)
    print('calculating solid Earth tides using PySolid (Milbert, 2018; Yunjun et al., 2022) ...')
    prog_bar = ptime.progressBar(maxValue=num_date, print_msg=not verbose)
    for i, dt_obj in enumerate(dt_objs):
        # calculate tide in ENU direction
        (tide_e,
         tide_n,
         tide_u) = pysolid.calc_solid_earth_tides_grid(dt_obj, atr_geo,
                                                       step_size=step_size,
                                                       display=False,
                                                       verbose=verbose)

        # convert ENU to LOS direction
        # sign convention: positive for motion towards satellite
        ts_tide[i,:,:] = (  tide_e * unit_vec[0]
                          + tide_n * unit_vec[1]
                          + tide_u * unit_vec[2])

        prog_bar.update(i+1, suffix='{} ({}/{})'.format(dt_obj.isoformat(), i+1, num_date))
    prog_bar.close()

    # radar-coding if input in radar-coordinates
    # use ts_file to avoid potential missing CENTER_LINE_UTC attributes in geom_file from alosStack
    atr = readfile.read_attribute(ts_file)
    if 'Y_FIRST' not in atr.keys():
        print('radar-coding the LOS tides time-series ...')
        res_obj = resample(lut_file=geom_file)
        res_obj.open()
        res_obj.src_meta = atr_geo
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        ts_tide = res_obj.run_resample(src_data=ts_tide[:,
                                                        box[1]:box[3],
                                                        box[0]:box[2]])

    ## output
    # attribute
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'
    for key in ['REF_Y', 'REF_X', 'REF_DATE']:
        if key in atr.keys():
            atr.pop(key)

    # write
    ds_dict = {}
    ds_dict['timeseries'] = ts_tide
    ds_dict['sensingMid'] = np.array([i.strftime('%Y%m%dT%H%M%S') for i in dt_objs], dtype=np.string_)
    writefile.write(ds_dict, out_file=set_file, metadata=atr, ref_file=ts_file)

    return set_file


def correct_timeseries(dis_file, set_file, cor_dis_file):
    """Correct time-series for the solid Earth tides."""
    # diff.py can handle different reference in space and time
    # between the absolute solid Earth tides and the double referenced time-series
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input time-series using diff.py')

    iargs = [dis_file, set_file, '-o', cor_dis_file]
    print('diff.py', ' '.join(iargs))

    import mintpy.cli.diff
    mintpy.cli.diff.main(iargs)

    return cor_dis_file


###############################################################
def run_solid_earth_tides(inps):
    start_time = time.time()

    # calc SET
    calc_solid_earth_tides_timeseries(
        ts_file=inps.dis_file,
        geom_file=inps.geom_file,
        set_file=inps.set_file,
        date_wise_acq_time=inps.date_wise_acq_time,
        update_mode=inps.update_mode,
        verbose=inps.verbose)

    # correct SET
    correct_timeseries(
        dis_file=inps.dis_file,
        set_file=inps.set_file,
        cor_dis_file=inps.cor_dis_file)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return
