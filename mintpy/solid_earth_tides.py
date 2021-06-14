#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2020                           #
############################################################
# Recomend import:
#   from mintpy import solid_earth_tides as SET


import os
import sys
import time
import datetime as dt
import argparse
import warnings
import h5py
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})

try:
    import pysolid
except ImportError:
    raise ImportError('Can not import pysolid! Check https://github.com/insarlab/PySolid.')

from mintpy.objects import timeseries
from mintpy.objects.resample import resample
from mintpy.defaults.template import get_template_content
from mintpy.utils import (
    ptime,
    readfile,
    writefile,
    utils as ut,
    attribute as attr,
)


###############################################################
TEMPLATE = get_template_content('correct_SET')

EXAMPLE = """example:
  solid_earth_tides.py timeseries.h5 -g inputs/geometryRadar.h5
  solid_earth_tides.py timeseries.h5 -g inputs/geometryGeo.h5
  solid_earth_tides.py geo/geo_timeseries_ERA5_demErr.h5 -g geo/geo_geometryRadar.h5
"""

REFERENCE = """reference:
  Milbert, D., Solid Earth Tide, http://geodesyworld.github.io/SOFTS/solid.htm, Accessd 2020 September 6.
  Fattahi, H., Z. Yunjun, X. Pi, P. S. Agram, P. Rosen, and Y. Aoki (2020), Absolute geolocation of SAR 
    Big-Data: The first step for operational InSAR time-series analysis, AGU Fall Meeting 2020, 1-17 Dec 2020.
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Solid Earth tides (SET) correction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='{}\n{}\n{}'.format(REFERENCE, TEMPLATE, EXAMPLE))

    parser.add_argument('dis_file', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g','--geomtry', dest='geom_file', type=str, required=True,
                        help='geometry file including incidence/azimuthAngle.')
    parser.add_argument('--date-wise-acq-time', dest='date_wise_acq_time', action='store_true',
                        help='Use the exact date-wise acquisition time instead of the common one for tides calculation.\n' +
                             'For ISCE-2/topsStack products only, and requires ../reference and ../secondarys folder.\n' +
                             'There is <1 min difference btw. S1A/B -> Negligible impact for InSAR.')

    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose message.')
    parser.add_argument('--update', dest='update_mode', action='store_true', help='Enable update mode.')

    # output
    parser.add_argument('--set-file', dest='set_file', help='line-of-sight solid earth tide file name')
    parser.add_argument('-o', dest='cor_dis_file', help='Output file name for the corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input files - processors & coordinates
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    coord1 = 'geo' if 'Y_FIRST' in atr1.keys() else 'radar'
    coord2 = 'geo' if 'Y_FIRST' in atr2.keys() else 'radar'
    proc = atr1.get('PROCESSOR', 'isce')

    if coord1 == 'radar' and proc in ['gamma', 'roipac']:
        msg = 'Radar-coded file from {} is NOT supported!'.format(proc)
        msg += '\n    Try to geocode the time-series and geometry files and re-run with them instead.'
        raise ValueError(msg)

    if coord1 != coord2:
        n = max(len(os.path.basename(i)) for i in [inps.dis_file, inps.geom_file])
        msg = 'Input time-series and geometry file are NOT in the same coordinate!'
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.dis_file),  n=n, c=coord1)
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.geom_file), n=n, c=coord2)
        raise ValueError(msg)

    # default SET filename
    if not inps.set_file:
        geom_dir = os.path.dirname(inps.geom_file)
        inps.set_file = os.path.join(geom_dir, 'SET.h5')

    # default corrected time-series filename
    if not inps.cor_dis_file:
        dis_dir = os.path.dirname(inps.dis_file)
        fbase, fext = os.path.splitext(os.path.basename(inps.dis_file))
        inps.cor_dis_file = os.path.join(dis_dir, '{}_SET{}'.format(fbase, fext))

    return inps


###############################################################
def prepare_los_geometry(geom_file):
    """Prepare LOS geometry data/info in geo-coordinates
    Parameters: geom_file  - str, path of geometry file
    Returns:    inc_angle  - 2D np.ndarray, incidence angle in radians
                head_angle - 2D np.ndarray, heading   angle in radians
                atr        - dict, metadata in geo-coordinate
    """

    print('prepare LOS geometry in geo-coordinates from file: {}'.format(geom_file))
    atr = readfile.read_attribute(geom_file)

    print('read incidenceAngle from file: {}'.format(geom_file))
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]

    if 'azimuthAngle' in readfile.get_dataset_list(geom_file):
        print('read azimuthAngle   from file: {}'.format(geom_file))
        print('convert azimuth angle to heading angle')
        az_angle  = readfile.read(geom_file, datasetName='azimuthAngle')[0]
        head_angle = ut.azimuth2heading_angle(az_angle)
    else:
        print('use the HEADING attribute as the mean heading angle')
        head_angle = np.ones(inc_angle.shape, dtype=np.float32) * float(atr['HEADING'])

    # geocode inc/az angle data if in radar-coord
    if 'Y_FIRST' not in atr.keys():
        print('-'*50)
        print('geocoding the incidence / heading angles ...')
        res_obj = resample(lut_file=geom_file, src_file=geom_file)
        res_obj.open()
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        inc_angle  = res_obj.run_resample(src_data=inc_angle[box[1]:box[3], box[0]:box[2]])
        head_angle = res_obj.run_resample(src_data=head_angle[box[1]:box[3], box[0]:box[2]])

        # update attribute
        atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)

    # for 'Y_FIRST' not in 'degree'
    # e.g. meters for UTM projection from ASF HyP3
    if not atr['Y_UNIT'].lower().startswith('deg'):
        # get SNWE in meter
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])
        N = float(atr['Y_FIRST'])
        W = float(atr['X_FIRST'])
        y_step = float(atr['Y_STEP'])
        x_step = float(atr['X_STEP'])
        S = N + y_step * length
        E = W + x_step * width

        # SNWE in meter --> degree
        lat0, lon0 = ut.to_latlon(atr['OG_FILE_PATH'], W, N)
        lat1, lon1 = ut.to_latlon(atr['OG_FILE_PATH'], E, S)
        lat_step = (lat1 - lat0) / length
        lon_step = (lon1 - lon0) / width

        # update Y/X_FIRST/STEP/UNIT
        atr['Y_FIRST'] = lat0
        atr['X_FIRST'] = lon0
        atr['Y_STEP'] = lat_step
        atr['X_STEP'] = lon_step
        atr['Y_UNIT'] = 'degrees'
        atr['X_UNIT'] = 'degrees'

    # unit: degree to radian
    inc_angle *= np.pi / 180.
    head_angle *= np.pi / 180.

    return inc_angle, head_angle, atr


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
        # opt 2. read sensingMid in xml files
        print('read exact datetime info in XML files from ISCE-2/topsStack results in directory:', proj_dir)
        from mintpy.utils import isce_utils
        sensingMid = isce_utils.get_sensing_datetime_list(proj_dir, date_list=date_list)[0]

        # plot
        plot_sensingMid_variation(sensingMid)

    else:
        # opt 3. use constant time of the day for all acquisitions
        msg =  'Use the same time of the day for all acquisitions from CENTER_LINE_UTC\n'
        msg += 'With <= 1 min variation for Sentinel-1A/B for example, this simplication has negligible impact on SET calculation.'
        print(msg)
        atr = readfile.read_attribute(ts_file)
        utc_sec = dt.timedelta(seconds=float(atr['CENTER_LINE_UTC']))
        sensingMid = [dt.datetime.strptime(i, '%Y%m%d') + utc_sec for i in date_list]

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
    inc_angle, head_angle, atr_geo = prepare_los_geometry(geom_file)

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

    # loop for calc
    print('\n'+'-'*50)
    print('calculating solid Earth tides using solid.for (D. Milbert, 2018) ...')
    prog_bar = ptime.progressBar(maxValue=num_date, print_msg=not verbose)
    for i, dt_obj in enumerate(dt_objs):
        # calculate tide in ENU direction
        (tide_e,
         tide_n,
         tide_u) = pysolid.calc_solid_earth_tides_grid(dt_obj, atr_geo,
                                                       display=False,
                                                       verbose=verbose)

        # convert ENU to LOS direction
        # sign convention: positive for motion towards satellite
        ts_tide[i,:,:] = (tide_e * unit_vec[0]
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
    from mintpy import diff

    iargs = [dis_file, set_file, '-o', cor_dis_file]
    print('diff.py', ' '.join(iargs))
    diff.main(iargs)
    return cor_dis_file


###############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
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

    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return inps.cor_dis_file

###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
