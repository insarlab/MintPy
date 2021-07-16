############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.utils import utils as ut


import os
import glob
import errno
import numpy as np
from scipy.ndimage import map_coordinates

from mintpy.objects import (
    geometryDatasetNames,
    geometry,
    ifgramStack,
    timeseries,
)

from mintpy.utils import ptime, readfile
from mintpy.utils.utils0 import *
from mintpy.utils.utils1 import *
from mintpy.objects.coord import coordinate


#################################################################################
def check_loaded_dataset(work_dir='./', print_msg=True, relpath=False):
    """Check the result of loading data for the following two rules:
        1. file existance
        2. file attribute readability

    Parameters: work_dir  : string, MintPy working directory
                print_msg : bool, print out message
    Returns:    True, if all required files and dataset exist; otherwise, ERROR
                    If True, PROCESS, SLC folder could be removed.
                stack_file  :
                geom_file   :
                lookup_file :
    Example:    work_dir = os.path.expandvars('./FernandinaSenDT128/mintpy')
                ut.check_loaded_dataset(work_dir)
    """
    load_complete = True

    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    # 1. interferograms stack file: unwrapPhase, coherence
    flist = [os.path.join(work_dir, 'inputs/ifgramStack.h5')]
    stack_file = is_file_exist(flist, abspath=True)
    if stack_file is not None:
        obj = ifgramStack(stack_file)
        obj.open(print_msg=False)
        for dname in ['unwrapPhase', 'coherence']:
            if dname not in obj.datasetNames and 'azimuthOffset' not in obj.datasetNames:
                raise ValueError('required dataset "{}" is missing in file {}'.format(dname, stack_file))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './inputs/ifgramStack.h5')

    atr = readfile.read_attribute(stack_file)

    # 2. geom_file: height
    if 'X_FIRST' in atr.keys():
        flist = [os.path.join(work_dir, 'inputs/geometryGeo.h5')]
    else:
        flist = [os.path.join(work_dir, 'inputs/geometryRadar.h5')]
    geom_file = is_file_exist(flist, abspath=True)
    if geom_file is not None:
        obj = geometry(geom_file)
        obj.open(print_msg=False)
        dname = geometryDatasetNames[0]
        if dname not in obj.datasetNames:
            raise ValueError('required dataset "{}" is missing in file {}'.format(dname, geom_file))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './inputs/geometry*.h5')

    # 3. lookup_file: latitude,longitude or rangeCoord,azimuthCoord
    # could be different than geometry file in case of roipac and gamma
    flist = [os.path.join(work_dir, 'inputs/geometry*.h5')]
    lookup_file = get_lookup_file(flist, abspath=True, print_msg=print_msg)
    if 'X_FIRST' not in atr.keys():
        if lookup_file is not None:
            obj = geometry(lookup_file)
            obj.open(print_msg=False)

            if atr['PROCESSOR'] in ['isce', 'doris']:
                dnames = geometryDatasetNames[1:3]
            elif atr['PROCESSOR'] in ['gamma', 'roipac']:
                dnames = geometryDatasetNames[3:5]
            else:
                raise AttributeError('InSAR processor: {}'.format(atr['PROCESSOR']))

            for dname in dnames:
                if dname not in obj.datasetNames:
                    load_complete = False
                    raise Exception('required dataset "{}" is missing in file {}'.format(dname, lookup_file))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './inputs/geometry*.h5')
    else:
        print("Input data seems to be geocoded. Lookup file not needed.")

    if relpath:
        stack_file  = os.path.relpath(stack_file)  if stack_file  else stack_file
        geom_file   = os.path.relpath(geom_file)   if geom_file   else geom_file
        lookup_file = os.path.relpath(lookup_file) if lookup_file else lookup_file

    # print message
    if print_msg:
        print(('Loaded dataset are processed by '
               'InSAR software: {}'.format(atr['PROCESSOR'])))
        if 'X_FIRST' in atr.keys():
            print('Loaded dataset is in GEO coordinates')
        else:
            print('Loaded dataset is in RADAR coordinates')
        print('Interferograms Stack: {}'.format(stack_file))
        print('Geometry File       : {}'.format(geom_file))
        print('Lookup Table File   : {}'.format(lookup_file))
        if load_complete:
            print('-'*50)
            print('All data needed found/loaded/copied. Processed 2-pass InSAR data can be removed.')
        print('-'*50)

    return load_complete, stack_file, geom_file, lookup_file


#################################################################################
def read_timeseries_lalo(lat, lon, ts_file, lookup_file=None, ref_lat=None, ref_lon=None,
                         win_size=1, unit='m', method='mean', print_msg=True):
    """ Read time-series of one pixel with input lat/lon
    Parameters: lat/lon     - float, latitude/longitude
                ts_file     - string, filename of time-series HDF5 file
                lookup_file - string, filename of lookup table file
                ref_lat/lon - float, latitude/longitude of reference pixel
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
        print('input lat / lon: {} / {}'.format(lat, lon))
        print('corresponding y / x: {} / {}'.format(y, x))

    # reference pixel
    ref_y, ref_x = None, None
    if ref_lat is not None:
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]

    # call read_timeseries_yx()
    dates, dis, dis_std = read_timeseries_yx(y, x, ts_file,
                                             ref_y=ref_y,
                                             ref_x=ref_x,
                                             win_size=win_size,
                                             unit=unit,
                                             method=method,
                                             print_msg=False)
    return dates, dis, dis_std


#################################################################################
def read_timeseries_yx(y, x, ts_file, ref_y=None, ref_x=None,
                       win_size=1, unit='m', method='mean', print_msg=True):
    """ Read time-series of one pixel with input y/x
    Parameters: y/x      - int, row/column number of interest
                ts_file  - string, filename of time-series HDF5 file
                ref_y/x  - int, row/column number of reference pixel
                win_size - int, windows size centered at point of interest
                unit     - str, output displacement unit
                method   - str, method to calculate the output displacement and its dispersity
    Returns:    dates    - 1D np.ndarray of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis      - 1D np.ndarray of float32, displacement
                dis_std  - 1D np.ndarray of float32, displacement dispersity
    """
    # read date
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    dates = ptime.date_list2vector(obj.dateList)[0]
    dates = np.array(dates)

    # read displacement
    if print_msg:
        print('input y / x: {} / {}'.format(y, x))
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
            raise ValueError('un-recognized method: {}'.format(method))

    # reference pixel
    if ref_y is not None:
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]

    #start at zero
    dis -= dis[0]

    # custom output unit
    if unit == 'm':
        pass
    elif unit == 'cm':
        dis *= 100.
        dis_std *= 100.
    elif unit == 'mm':
        dis *= 1000.
        dis_std *= 1000.
    else:
        raise ValueError('un-supported output unit: {}'.format(unit))

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
        msg += '\nstart_yx: {}'.format(start_yx)
        msg += '\nend_yx:{}'.format(end_yx)
        msg += '\ndata size: ({}, {})'.format(length, width)
        raise ValueError(msg)

    # Determine points coordinates along the line
    num_pts = int(np.hypot(x1-x0, y1-y0))
    ys = np.linspace(y0, y1, num_pts, dtype=np.float32)
    xs = np.linspace(x0, x1, num_pts, dtype=np.float32)

    # Extract z value along the line
    # for nearest neighbor sampling, use indexing directly
    # for other interpolation, use scipy.ndimage.map_coordinates
    if interpolation == 'nearest':
        z_line = z[np.rint(ys).astype(np.int), np.rint(xs).astype(np.int)]

    else:
        # interpolation name to order
        interpolate_name2order = {
            'linear' : 1,
            'cubic'  : 3,
            'quintic': 5,
        }
        if interpolation not in interpolate_name2order.keys():
            msg = 'un-supported interpolation method: {}'.format(interpolation)
            msg += '\navailable methods: {}'.format(interpolate_name2order.keys())
            raise ValueError(msg)
        interp_order = interpolate_name2order[interpolation.lower()]
        # run interpolation
        z_line = map_coordinates(z, np.vstack((ys, xs)), order=interp_order)

    # Calculate Distance along the line
    earth_radius = 6.3781e6    # in meter
    if 'Y_FIRST' in atr.keys():
        [lat0, lat1] = coordinate(atr).yx2lalo([y0, y1], coord_type='y')
        lat_c = (lat0 + lat1) / 2.
        x_step = float(atr['X_STEP']) * np.pi/180.0 * earth_radius * np.cos(lat_c * np.pi/180)
        y_step = float(atr['Y_STEP']) * np.pi/180.0 * earth_radius
    else:
        x_step = range_ground_resolution(atr)
        y_step = azimuth_ground_resolution(atr)
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

    for i in range(len(lines)):
        # read segment data
        start_lalo, end_lalo = lines[i][0], lines[i][1]
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

