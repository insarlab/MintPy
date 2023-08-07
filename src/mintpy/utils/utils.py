"""Miscellaneous utilities - dependent on utils0/1."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.utils import utils as ut


import errno
import os

import numpy as np
from scipy.ndimage import map_coordinates

from mintpy.objects import (
    GEOMETRY_DSET_NAMES,
    geometry,
    ifgramStack,
    timeseries,
)
from mintpy.objects.coord import coordinate
from mintpy.objects.resample import resample
from mintpy.utils import attribute as attr, ptime, readfile
from mintpy.utils.utils0 import *
from mintpy.utils.utils1 import *


#################################################################################
def check_loaded_dataset(work_dir='./', print_msg=True, relpath=False):
    """Check the loaded input files, following two rules:
        1. file existence
        2. file attribute readability

    Parameters: work_dir    - str, MintPy working directory
                print_msg   - bool, print out message
    Returns:    stack_file  - str, path to the interferogram stack file
                geom_file   - str, path to the geometry file
                lookup_file - str, path to the look up table file, for radar-coord dataset only.
                ion_file    - str, path to the ionosphere stack file
    Example:    work_dir = os.path.expandvars('./FernandinaSenDT128/mintpy')
                stack_file, geom_file, lookup_file = ut.check_loaded_dataset(work_dir)[:3]
    """
    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    # tips for prep_aria
    template_file = os.path.join(work_dir, 'smallbaselineApp.cfg')
    proc = readfile.read_template(template_file)['mintpy.load.processor']
    if proc == 'aria':
        msg_aria = '. Re-run "prep_aria.py" as printed out in "load_data" step for more information!'
    else:
        msg_aria = ''

    # 1. [required] interferograms stack file: unwrapPhase, coherence
    stack_file = os.path.join(work_dir, 'inputs/ifgramStack.h5')
    dnames = ['unwrapPhase', 'rangeOffset', 'azimuthOffset']

    if is_file_exist(stack_file, abspath=True):
        obj = ifgramStack(stack_file)
        obj.open(print_msg=False)

        if all(x not in obj.datasetNames for x in dnames):
            msg = f'required dataset is missing in file {stack_file}:\n'
            msg += ' OR '.join(dnames)
            raise ValueError(msg)

        # check coherence for phase stack
        if 'unwrapPhase' in obj.datasetNames and 'coherence' not in obj.datasetNames:
            print(f'WARNING: "coherence" is missing in file {stack_file}')
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), stack_file + msg_aria)

    # get coordinate type of the loaded dataset
    atr = readfile.read_attribute(stack_file)
    coord_type = 'GEO' if 'Y_FIRST' in atr.keys() else 'RADAR'
    processor = atr['PROCESSOR']

    # 2. [required] geom_file: height
    geom_file = os.path.join(work_dir, 'inputs', f'geometry{coord_type.capitalize()}.h5')
    dname = GEOMETRY_DSET_NAMES[0]

    if is_file_exist(geom_file, abspath=True):
        obj = geometry(geom_file)
        obj.open(print_msg=False)
        if dname not in obj.datasetNames:
            raise ValueError(f'required dataset "{dname}" is missing in file {geom_file}')
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), geom_file + msg_aria)

    # 3. [required for radar-coord] lookup_file: latitude,longitude or rangeCoord,azimuthCoord
    # could be different than geometry file in case of roipac and gamma
    lookup_file = os.path.join(work_dir, 'inputs/geometry*.h5')
    lookup_file = get_lookup_file(lookup_file, abspath=True, print_msg=print_msg)
    if coord_type == 'RADAR':
        if lookup_file is not None:
            obj = geometry(lookup_file)
            obj.open(print_msg=False)

            # get the proper lookup table dataset names
            if processor in ['isce', 'doris']:
                dnames = GEOMETRY_DSET_NAMES[1:3]
            elif processor in ['gamma', 'roipac']:
                dnames = GEOMETRY_DSET_NAMES[3:5]
            else:
                msg = f'Unknown InSAR processor: {processor} to locate look up table!'
                raise AttributeError(msg)

            for dname in dnames:
                if dname not in obj.datasetNames:
                    raise Exception(f'required dataset "{dname}" is missing in file {lookup_file}')
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), lookup_file)
    else:
        print("Input data seems to be geocoded. Lookup file not needed.")

    # 4. [optional] ionosphere stack file: unwrapPhase, coherence
    ion_file = os.path.join(work_dir, 'inputs/ionStack.h5')
    dname = 'unwrapPhase'

    if is_file_exist(ion_file, abspath=True):
        obj = ifgramStack(ion_file)
        obj.open(print_msg=False)
        if dname not in obj.datasetNames:
            raise ValueError(f'required dataset "{dname}" is missing in file {ion_file}')

        # check coherence for phase stack
        if 'unwrapPhase' in obj.datasetNames and 'coherence' not in obj.datasetNames:
            print(f'WARNING: "coherence" is missing in file {ion_file}')
    else:
        ion_file = None

    if relpath:
        stack_file  = os.path.relpath(stack_file)  if stack_file  else stack_file
        geom_file   = os.path.relpath(geom_file)   if geom_file   else geom_file
        lookup_file = os.path.relpath(lookup_file) if lookup_file else lookup_file
        ion_file    = os.path.relpath(ion_file)    if ion_file    else ion_file

    # print message
    if print_msg:
        msg  = f'Loaded dataset are processed by InSAR software: {processor}'
        msg += f'\nLoaded dataset are in {coord_type} coordinates'
        msg += f'\nInterferogram Stack: {stack_file}'
        msg += f'\nIonosphere    Stack: {ion_file}' if ion_file else ''
        msg += f'\nGeometry      File : {geom_file}'
        msg += f'\nLookup Table  File : {lookup_file}'
        msg += '\n' + '-' * 50
        print(msg)

    return stack_file, geom_file, lookup_file, ion_file


#################################################################################
def read_timeseries_lalo(lat, lon, ts_file, lookup_file=None, ref_lat=None, ref_lon=None,
                         zero_first=True, win_size=1, unit='m', method='mean', print_msg=True):
    """ Read time-series of one pixel with input lat/lon
    Parameters: lat/lon     - float, latitude/longitude
                ts_file     - string, filename of time-series HDF5 file
                lookup_file - string, filename of lookup table file
                ref_lat/lon - float, latitude/longitude of reference pixel
                zero_first  - bool, shift the time-series so that it starts from zero
                win_size    - int, windows size centered at point of interest
                unit        - str, output displacement unit
                method      - str, method to calculate the output displacement and its dispersity
    Returns:    dates       - 1D np.ndarray of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis         - 1D np.ndarray of float32, displacement
                dis_std     - 1D np.ndarray of float32, displacement dispersity
    """
    atr = readfile.read_attribute(ts_file)
    coord = coordinate(atr, lookup_file=lookup_file)
    y, x = coord.geo2radar(lat, lon)[0:2]
    if print_msg:
        print(f'input lat / lon: {lat} / {lon}')
        print(f'corresponding y / x: {y} / {x}')

    # reference pixel
    ref_y, ref_x = None, None
    if ref_lat is not None:
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]

    # call read_timeseries_yx()
    dates, dis, dis_std = read_timeseries_yx(y, x, ts_file,
                                             ref_y=ref_y,
                                             ref_x=ref_x,
                                             zero_first=zero_first,
                                             win_size=win_size,
                                             unit=unit,
                                             method=method,
                                             print_msg=False)
    return dates, dis, dis_std


def read_timeseries_yx(y, x, ts_file, ref_y=None, ref_x=None, zero_first=True,
                       win_size=1, unit='m', method='mean', print_msg=True):
    """ Read time-series of one pixel with input y/x
    Parameters: y/x        - int, row/column number of interest
                ts_file    - string, filename of time-series HDF5 file
                ref_y/x    - int, row/column number of reference pixel
                zero_first - bool, shift the time-series so that it starts from zero
                win_size   - int, windows size centered at point of interest
                unit       - str, output displacement unit
                method     - str, method to calculate the output displacement and its dispersity
    Returns:    dates      - 1D np.ndarray of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis        - 1D np.ndarray of float32, displacement
                dis_std    - 1D np.ndarray of float32, displacement dispersity
    """
    # read date
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    dates = ptime.date_list2vector(obj.dateList)[0]
    dates = np.array(dates)

    # read displacement
    if print_msg:
        print(f'input y / x: {y} / {x}')
    box = (x, y, x+1, y+1)
    dis = readfile.read(ts_file, box=box)[0]
    dis_std = None

    if win_size != 1:
        buf = int(win_size / 2)
        box_win = (x-buf, y-buf, x+buf+1, y+buf+1)
        dis_win = readfile.read(ts_file, box=box_win)[0].reshape(obj.numDate, -1)

        if method == 'mean':
            dis = np.nanmean(dis_win, axis=1)
            dis_std = np.nanstd(dis_win, axis=1)

        elif method == 'median':
            dis = np.nanmedian(dis_win, axis=1)
            dis_std = median_abs_deviation(dis_win)

        else:
            raise ValueError(f'un-recognized method: {method}')

    # reference pixel
    if ref_y is not None:
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]

    #start at zero
    if zero_first:
        dis -= dis[0]

    # custom output unit
    if unit == 'm':
        pass
    elif unit == 'cm':
        dis *= 100.
        dis_std = None if dis_std is None else dis_std * 100.
    elif unit == 'mm':
        dis *= 1000.
        dis_std = None if dis_std is None else dis_std * 1000.
    else:
        raise ValueError(f'un-supported output unit: {unit}')

    return dates, dis, dis_std


#####################################################################
def transect_yx(z, atr, start_yx, end_yx, interpolation='nearest'):
    """Extract 2D matrix (z) value along the line [x0,y0;x1,y1]
    Link: http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array

    Parameters: z : (np.array) 2D data matrix
                atr : (dict) attribute
                start_yx : (list) y,x coordinate of start point
                end_yx : (list) y,x coordinate of end   point
                interpolation : str, sampling/interpolation method, including:
                    'nearest' - nearest neighbour
                    'linear'  - linear  spline interpolation (order of 1)
                    'cubic'   - cubic   spline interpolation (order of 3)
                    'quintic' - quintic spline interpolation (order of 5)

    Returns:    transect: (dict) containing 1D matrix:
                    'X' - 1D np.array for X/column coordinates in float32
                    'Y' - 1D np.array for Y/row.   coordinates in float32
                    'value' - 1D np.array for z value in float32
                    'distance' - 1D np.array for distance in float32

    Example: transect = transect_yx(dem, demRsc, [10,15], [100,115])
    """
    interpolation = interpolation.lower()
    [y0, x0] = start_yx
    [y1, x1] = end_yx

    # check
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    if not all(0<= i < width and 0<= j < length for i,j in zip([x0,x1], [y0,y1])):
        msg = 'input start/end point is out of data coverage'
        msg += f'\nstart_yx: {start_yx}'
        msg += f'\nend_yx:{end_yx}'
        msg += f'\ndata size: ({length}, {width})'
        raise ValueError(msg)

    # Determine points coordinates along the line
    num_pts = int(np.hypot(x1-x0, y1-y0))
    ys = np.linspace(y0, y1, num_pts, dtype=np.float32)
    xs = np.linspace(x0, x1, num_pts, dtype=np.float32)

    # Extract z value along the line
    # for nearest neighbor sampling, use indexing directly
    # for other interpolation, use scipy.ndimage.map_coordinates
    if interpolation == 'nearest':
        z_line = z[np.rint(ys).astype(int), np.rint(xs).astype(int)]

    else:
        # interpolation name to order
        interpolate_name2order = {
            'linear' : 1,
            'cubic'  : 3,
            'quintic': 5,
        }
        if interpolation not in interpolate_name2order.keys():
            msg = f'un-supported interpolation method: {interpolation}'
            msg += f'\navailable methods: {interpolate_name2order.keys()}'
            raise ValueError(msg)
        interp_order = interpolate_name2order[interpolation.lower()]
        # run interpolation
        z_line = map_coordinates(z, np.vstack((ys, xs)), order=interp_order)

    # Calculate Distance along the line
    earth_radius = 6.3781e6    # in meter
    dist_unit = 'm'
    if 'Y_FIRST' in atr.keys():
        [lat0, lat1] = coordinate(atr).yx2lalo([y0, y1], coord_type='y')
        lat_c = (lat0 + lat1) / 2.
        x_step = float(atr['X_STEP']) * np.pi/180.0 * earth_radius * np.cos(lat_c * np.pi/180)
        y_step = float(atr['Y_STEP']) * np.pi/180.0 * earth_radius
    else:
        try:
            x_step = range_ground_resolution(atr)
            y_step = azimuth_ground_resolution(atr)
        except KeyError:
            x_step = 1
            y_step = 1
            dist_unit = 'pixel'
    dist_line = np.hypot((xs - x0) * x_step,
                         (ys - y0) * y_step)

    # remove points in masked out areas
    mask = ~np.isnan(z_line)
    mask *= z_line != 0.0

    # prepare output
    transect = {}
    transect['Y'] = ys[mask]
    transect['X'] = xs[mask]
    transect['value'] = z_line[mask]
    transect['distance'] = dist_line[mask]
    transect['distance_unit'] = dist_unit

    return transect


def transect_lalo(z, atr, start_lalo, end_lalo, interpolation='nearest'):
    """Extract 2D matrix (z) value along the line [start_lalo, end_lalo]"""
    coord = coordinate(atr)
    [y0, y1] = coord.lalo2yx([start_lalo[0], end_lalo[0]], coord_type='lat')
    [x0, x1] = coord.lalo2yx([start_lalo[1], end_lalo[1]], coord_type='lon')
    transect = transect_yx(z, atr, [y0, x0], [y1, x1], interpolation)
    return transect


def transect_lines(z, atr, lines):
    """Extract 2D matrix (z) value along multiple lines
    Parameters: z     : 2D np.ndarray in size of (l,w)
                atr   : dict, metadata of matrix z
                lines : list of lines with each line is defined as:
                    [[lat0, lon0], [lat1, lon1]] for geo coordinates
                    [[y0, x0], [y1, x1]] for radar coordinates
    Returns: transect : (dict) containing 1D matrix:
                    'X' - 1D np.array for X/column coordinates in float32
                    'Y' - 1D np.array for Y/row.   coordinates in float32
                    'value' - 1D np.array for z value in float32
                    'distance' - 1D np.array for distance in float32
    """
    transect = {}
    start_distance = 0
    transect['start_distance'] = []

    for i, line in enumerate(lines):
        # read segment data
        start_lalo, end_lalo = line[0], line[1]
        if 'Y_FIRST' in atr.keys():
            seg = transect_lalo(z, atr, start_lalo, end_lalo)
        else:
            seg = transect_yx(z, atr, start_lalo, end_lalo)
        seg['distance'] += start_distance

        # connect each segment
        if i == 0:
            # first segment
            for key, value in seg.items():
                transect[key] = np.array(value, dtype=np.float32)
        else:
            for key, value in seg.items():
                transect[key] = np.concatenate((transect[key], value))

        # update start_distance for the next segment
        transect['start_distance'].append(start_distance)
        start_distance = transect['distance'][-1]
    transect['start_distance'] = np.array(transect['start_distance'], dtype=np.float32)
    return transect



#################################################################################

def prepare_geo_los_geometry(geom_file, unit='rad'):
    """Prepare LOS geometry data/info in geo-coordinates.

    Parameters: geom_file - str, path of geometry file
                unit      - str, rad or deg, output angle unit
    Returns:    inc_angle - 2D np.ndarray, incidence angle in radians / degrees
                            measured from the vertical
                az_angle  - 2D np.ndarray, azimuth   angle in radians / degrees
                            measured from the north with anti-clockwise direction as positive
                atr       - dict, metadata in geo-coordinate
    """

    print(f'prepare LOS geometry in geo-coordinates from file: {geom_file}')
    atr = readfile.read_attribute(geom_file)

    print(f'read incidenceAngle from file: {geom_file}')
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]

    if 'azimuthAngle' in readfile.get_dataset_list(geom_file):
        print(f'read azimuthAngle   from file: {geom_file}')
        az_angle  = readfile.read(geom_file, datasetName='azimuthAngle')[0]
    else:
        print('use the HEADING attribute as the mean heading angle')
        print('convert heading angle to azimuth angle')
        head_angle = np.ones(inc_angle.shape, dtype=np.float32) * float(atr['HEADING'])
        az_angle = heading2azimuth_angle(head_angle)

    # geocode inc/az angle data if in radar-coord
    if 'Y_FIRST' not in atr.keys():
        print('-'*50)
        print('geocoding the incidence / heading angles ...')
        res_obj = resample(lut_file=geom_file, src_file=geom_file)
        res_obj.open()
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        inc_angle = res_obj.run_resample(src_data=inc_angle[box[1]:box[3], box[0]:box[2]])
        az_angle = res_obj.run_resample(src_data=az_angle[box[1]:box[3], box[0]:box[2]])

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
        lat0, lon0 = utm2latlon(atr, W, N)
        lat1, lon1 = utm2latlon(atr, E, S)
        lat_step = (lat1 - lat0) / length
        lon_step = (lon1 - lon0) / width

        # update Y/X_FIRST/STEP/UNIT
        atr['Y_FIRST'] = lat0
        atr['X_FIRST'] = lon0
        atr['Y_STEP'] = lat_step
        atr['X_STEP'] = lon_step
        atr['Y_UNIT'] = 'degrees'
        atr['X_UNIT'] = 'degrees'

    # set invalid values to nan
    inc_angle[inc_angle == 0] = np.nan

    # unit: degree to radian
    if unit.startswith('rad'):
        inc_angle *= np.pi / 180.
        az_angle *= np.pi / 180.

    return inc_angle, az_angle, atr
