############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Feb 2021                           #
############################################################


import datetime as dt
import os
import time

import h5py
import numpy as np

import mintpy.cli.diff
from mintpy.constants import SPEED_OF_LIGHT
from mintpy.objects import ionex, timeseries
from mintpy.simulation import iono
from mintpy.utils import ptime, readfile, utils as ut, writefile


#####################################################################################
def get_dataset_size(fname):
    atr = readfile.read_attribute(fname)
    shape = (int(atr['LENGTH']), int(atr['WIDTH']))
    return shape


def run_or_skip(iono_file, grib_files, dis_file, geom_file):
    print('update mode: ON')
    print(f'output file: {iono_file}')
    flag = 'skip'

    # check existence and modification time
    if ut.run_or_skip(out_file=iono_file, in_file=grib_files, print_msg=False) == 'run':
        flag = 'run'
        print('1) output file either do NOT exist or is NOT newer than all IONEX files.')

    else:
        print('1) output file exists and is newer than all IONEX files.')

        # check dataset size in space / time
        ds_size_dis = get_dataset_size(dis_file)
        ds_size_ion = get_dataset_size(geom_file)
        date_list_dis = timeseries(dis_file).get_date_list()
        date_list_ion = timeseries(iono_file).get_date_list()
        if ds_size_ion != ds_size_dis or any (x not in date_list_ion for x in date_list_dis):
            flag = 'run'
            print(f'2) output file does NOT have the same len/wid as the geometry file {geom_file}'
                  ' or does NOT contain all dates')
        else:
            print('2) output file has the same len/wid as the geometry file and contains all dates')

            # check if output file is fully written
            with h5py.File(iono_file, 'r') as f:
                if np.all(f['timeseries'][-1,:,:] == 0):
                    flag = 'run'
                    print('3) output file is NOT fully written.')
                else:
                    print('3) output file is fully written.')

    # result
    print(f'run or skip: {flag}')
    return flag


#####################################################################################
def download_ionex_files(date_list, tec_dir, sol_code='jpl'):
    """Download IGS TEC products in IONEX format for the input list of dates.

    Parameters: date_list - list of str, in YYYYMMDD
                tec_dir   - str, path to IGS_TEC directory, e.g. ~/data/aux/IGS_TEC
                sol_code   - str, TEC solution center, e.g. jpl, cod, igs
    Returns:    fnames    - list of str, path of the downloaded TEC files
    """
    print("\n------------------------------------------------------------------------------")
    print("downloading GNSS-based TEC products in IONEX format from NASA/CDDIS ...")
    print('https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html')
    num_date = len(date_list)
    n = len(str(num_date))
    print(f'number of TEC files to download: {num_date}')
    print(f'local TEC file directory: {tec_dir}')

    # output file names/sizes
    fnames = []
    for date_str in date_list:
        fnames.append(ionex.get_ionex_filename(date_str, tec_dir=tec_dir, sol_code=sol_code))

    # remove all existing files
    debug_mode = False
    if debug_mode:
        for fname in fnames:
            for x in [fname, fname+'.Z']:
                if os.path.isfile(x):
                    os.remove(x)

    fsizes = [os.path.getsize(i) / 1024 if os.path.isfile(i) else 0 for i in fnames]

    # download: skip existing ones
    fsizec = ut.most_common(fsizes)
    if fsizec < 400:
        # too small, does not seem right --> download them all
        date_list2dload = list(date_list)

    else:
        # download missing ones
        date_list2dload = [d for d, s in zip(date_list, fsizes) if s < fsizec * 0.9]

    num_date2dload = len(date_list2dload)
    if num_date2dload == 0:
        print(f'ALL files exists with consistent file size (~{fsizec:.0f} KB)'
              ' --> skip re-downloading.\n')

    else:
        for i, date_str in enumerate(date_list2dload):
            print('-'*20)
            print(f'DATE {i+1}/{num_date2dload}: {date_str}')
            ionex.dload_ionex(date_str, tec_dir=tec_dir, sol_code=sol_code, print_msg=True)

        # print file size info, after downloading
        fsizes = [os.path.getsize(i) / 1024 if os.path.isfile(i) else 0 for i in fnames]
        for i in range(num_date):
            print(f'[{i+1:0{n}d}/{num_date}] {fnames[i]}: {fsizes[i]:.2f} KB')

    return fnames


def calc_iono_ramp_timeseries_igs(tec_dir, sol_code, interp_method, ts_file, geom_file, iono_file,
                                  rotate_tec_map=True, sub_tec_ratio=None, update_mode=True):
    """Calculate the time-series of 2D ionospheric delay from IGS TEC data.
    Considering the variation of the incidence angle along range direction.

    Parameters: tec_dir   - str, path of the local TEC directory
                ts_file   - str, path of the time-series file
                geom_file - str, path of the geometry file including incidenceAngle data
                iono_file - str, path of output iono ramp time-series file
    Returns:    iono_file - str, path of output iono ramp time-series file
    """
    print("\n------------------------------------------------------------------------------")
    # prepare geometry
    iono_lat, iono_lon = iono.prep_geometry_iono(geom_file, print_msg=True)[1:3]

    # prepare date/time
    date_list = timeseries(ts_file).get_date_list()
    meta = readfile.read_attribute(ts_file)
    utc_sec = float(meta['CENTER_LINE_UTC'])
    print(f'CENTER_LINE_TUC: {utc_sec}')

    # UTC time & local solar time
    # use an arbitrary date to construct the datetime object
    lon_c = (float(meta['LON_REF1']) + float(meta['LON_REF2'])) / 2
    utc_dt = dt.datetime(2020, 1, 1) + dt.timedelta(seconds=utc_sec)
    local_dt = ptime.utc2solar_time(utc_dt, lon_c)
    print('UTC time:', utc_dt.strftime("%H:%M:%S"))
    print('Local solar time:', local_dt.strftime("%I:%M %p"))

    # read IGS TEC
    print('read IGS TEC file ...')
    print(f'interpolation method: {interp_method}')
    if interp_method == 'linear3d':
        print(f'rotate TEC maps: {rotate_tec_map}')

    vtec_list = []
    prog_bar = ptime.progressBar(maxValue=len(date_list))
    for i, date_str in enumerate(date_list):
        # read zenith TEC
        tec_file = ionex.get_ionex_filename(date_str, tec_dir=tec_dir, sol_code=sol_code)
        vtec = ionex.get_ionex_value(tec_file, utc_sec,
                                     lat=iono_lat, lon=iono_lon,
                                     interp_method=interp_method,
                                     rotate_tec_map=rotate_tec_map)
        vtec_list.append(vtec)
        prog_bar.update(i+1, suffix=date_str)
    prog_bar.close()

    # TEC --> iono ramp
    vtec2iono_ramp_timeseries(
        date_list=date_list,
        vtec_list=vtec_list,
        geom_file=geom_file,
        iono_file=iono_file,
        sub_tec_ratio=sub_tec_ratio,
        update_mode=update_mode,
    )

    return iono_file


def vtec2iono_ramp_timeseries(date_list, vtec_list, geom_file, iono_file, sub_tec_ratio=None,
                              ds_dict_ext=None, update_mode=True):
    """Convert zenith TEC to 2D slant range delay (ramp due to the incidence angle variation)
    and write to HDF5 time-series file.

    Parameters: date_list   - list of str, dates in YYYYMMDD format
                vtec_list   - list of float32, zenith TEC in TECU
                geom_file   - str, path of the geometry file including incidenceAngle data
                iono_file   - str, path of output iono ramp time-series file
                update_mode - bool,
                ds_dict_ext - dict, extra dictionary of dataset to be saved into the HDF5 file.
    Returns:    iono_file   - str, path of output iono ramp time-series file
    """
    top_perc_file = os.path.join(os.path.dirname(mintpy.__file__), 'data', 'top_tec_perc_s1.txt')

    # prepare geometry
    (iono_inc_angle,
     iono_lat,
     iono_lon,
     iono_height) = iono.prep_geometry_iono(geom_file, print_msg=True)

    # prepare date/time
    num_date = len(date_list)
    if len(vtec_list) != num_date:
        msg = 'Input tec_list and date_list have different size!'
        msg += '\nFor acquisitions without TEC data, set it to NaN.'
        raise ValueError(msg)

    meta = readfile.read_attribute(geom_file)
    length = int(meta['LENGTH'])
    width = int(meta['WIDTH'])
    freq = SPEED_OF_LIGHT / float(meta['WAVELENGTH'])

    # Note: Scaling gives slightly better RMSE for SenD but much worse RMSE for SenA and Alos2
    # thus this is not used by default.
    if sub_tec_ratio is not None:
        if ut.is_number(sub_tec_ratio):
            print(f'multiply VTEC by {sub_tec_ratio}')
            vtec_list = (np.array(vtec_list).flatten() * float(sub_tec_ratio)).tolist()

        elif sub_tec_ratio.startswith('adap'):
            dates = ptime.date_list2vector(date_list)[0]
            ydays = np.array([x.timetuple().tm_yday for x in dates])
            fc = np.loadtxt(top_perc_file, dtype=bytes).astype(np.float32)
            print(f'multiply VTEC adaptively based on the day of the year from: {top_perc_file}')
            sub_perc = fc[:,2][np.array(ydays)]
            vtec_list = (np.array(vtec_list).flatten() * sub_perc).tolist()

    # loop to calculate the range delay (ramp)
    print('calculating ionospheric phase ramp time-series from TEC ...')
    ts_ramp = np.zeros((num_date, length, width), dtype=np.float32)
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i, date_str in enumerate(date_list):
        ts_ramp[i,:,:] = iono.vtec2range_delay(
            vtec_list[i],
            inc_angle=iono_inc_angle,
            freq=freq,
        )
        prog_bar.update(i+1, suffix=date_str)
    prog_bar.close()

    ## output
    # prepare metadata
    meta['FILE_TYPE'] = 'timeseries'
    meta['UNIT'] = 'm'
    meta['IONO_LAT'] = iono_lat
    meta['IONO_LON'] = iono_lon
    meta['IONO_HEIGHT'] = iono_height
    meta['IONO_INCIDENCE_ANGLE'] = np.nanmean(iono_inc_angle)
    # absolute delay without double reference
    for key in ['REF_X','REF_Y','REF_LAT','REF_LON','REF_DATE']:
        if key in meta.keys():
            meta.pop(key)

    # prepare data matrix
    ds_dict = {}
    ds_dict['date'] = np.array(date_list, dtype=np.string_)
    ds_dict['vtec'] = np.array(vtec_list, dtype=np.float32)
    ds_dict['timeseries'] = ts_ramp

    # add the extra dataset if specified, e.g. vtec_gim, vtec_top, vtec_sub
    if ds_dict_ext is not None:
        ds_names = ds_dict.keys()
        for ds_name, ds_val in ds_dict_ext.items():
            if ds_name not in ds_names:
                ds_dict[ds_name] = ds_val

    # write to disk
    writefile.write(ds_dict, iono_file, metadata=meta)

    return iono_file


def run_iono_tec(inps):
    """Calculate and/or correct for the ionospheric delay using TEC from GIM."""

    start_time = time.time()

    # download
    date_list = timeseries(inps.dis_file).get_date_list()
    tec_files = download_ionex_files(date_list, tec_dir=inps.tec_dir, sol_code=inps.sol_code)

    # calculate
    if run_or_skip(inps.iono_file, tec_files, inps.dis_file, inps.geom_file) == 'run':
        calc_iono_ramp_timeseries_igs(
            tec_dir=inps.tec_dir,
            sol_code=inps.sol_code,
            interp_method=inps.interp_method,
            ts_file=inps.dis_file,
            geom_file=inps.geom_file,
            iono_file=inps.iono_file,
            rotate_tec_map=inps.rotate_tec_map,
            sub_tec_ratio=inps.sub_tec_ratio,
            update_mode=inps.update_mode,
        )

    ## correct (using diff.py)
    ## diff.py can handle different reference in space and time
    ## e.g. the absolute delay and the double referenced time-series
    #print('correcting delay for using diff.py')
    #iargs = [inps.dis_file, inps.iono_file, '-o', inps.cor_dis_file, '--force']
    #print('diff.py', ' '.join(iargs))
    #mintpy.cli.diff.main(iargs)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return
