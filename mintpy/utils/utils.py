############################################################
# Program is part of MintPy                                #
# Copyright(c) 2013-2019, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend import:
#   from mintpy.utils import utils as ut


import os
import shutil
import errno
import numpy as np

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
def check_loaded_dataset(work_dir='./', print_msg=True):
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
    Example:    work_dir = os.path.expandvars('$SCRATCHDIR/SinabungT495F50AlosA/mintpy')
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
            if dname not in obj.datasetNames:
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
                dnames = [geometryDatasetNames[1],
                         geometryDatasetNames[2]]
            elif atr['PROCESSOR'] in ['gamma', 'roipac']:
                dnames = [geometryDatasetNames[3],
                      geometryDatasetNames[4]]
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
                         win_size=1, unit='m', print_msg=True):
    """ Read time-series of one pixel with input lat/lon
    Parameters: lat/lon     : float, latitude/longitude
                ts_file     : string, filename of time-series HDF5 file
                lookup_file : string, filename of lookup table file
                ref_lat/lon : float, latitude/longitude of reference pixel
    Returns:    dates : 1D np.array of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis   : 1D np.array of float in meter
    """
    atr = readfile.read_attribute(ts_file)
    coord = coordinate(atr, lookup_file=lookup_file)
    y, x = coord.geo2radar(lat, lon)[0:2]
    if print_msg:
        print('input lat / lon: {} / {}'.format(lat, lon))

    # reference pixel
    ref_y, ref_x = None, None
    if ref_lat is not None:
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]

    # call read_timeseries_yx()
    dates, dis = read_timeseries_yx(y, x, ts_file,
                                    ref_y=ref_y,
                                    ref_x=ref_x,
                                    win_size=win_size,
                                    unit=unit,
                                    print_msg=False)
    return dates, dis


#################################################################################
def read_timeseries_yx(y, x, ts_file, ref_y=None, ref_x=None,
                       win_size=1, unit='m', print_msg=True):
    """ Read time-series of one pixel with input y/x
    Parameters: y/x         : int, row/column number of interest
                ts_file     : string, filename of time-series HDF5 file
                ref_y/x     : int, row/column number of reference pixel
    Returns:    dates : 1D np.array of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis   : 1D np.array of float in meter
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
    if win_size != 1:
        buf = int(win_size / 2)
        box_win = (x-buf, y-buf, x+buf+1, y+buf+1)
        dis_win = readfile.read(ts_file, box=box_win)[0]
        dis = np.nanmean(dis_win.reshape((obj.numDate, -1)), axis=1)

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
    elif unit == 'mm':
        dis *= 1000.
    else:
        raise ValueError('un-supported output unit: {}'.format(unit))

    return dates, dis

#################################################################################
def move_dask_stdout_stderr_files():
    """ move  *o and *e files produced by dask into stdout and sderr directory """

    stdout_files  = glob.glob('*.o')
    stderr_files  = glob.glob('*.e')
    job_files = glob.glob('dask_command_run_from_python.txt*')

    stdout_folder = 'stdout_ifgram_inversion_dask'
    stderr_folder = 'stderr_ifgram_inversion_dask'
    for std_dir in [stdout_folder, stderr_folder]:
        if os.path.isdir(std_dir):
            shutil.rmtree(std_dir)
        os.mkdir(std_dir)

    for item in stdout_files + job_files:
        shutil.move(item, stdout_folder)
    for item in stderr_files:
        shutil.move(item, stderr_folder)

    return
