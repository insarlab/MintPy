############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend import:
#   from pysar.utils import utils as ut


import os
import sys
import glob
import time
import datetime
import errno
from argparse import Namespace
import h5py
import numpy as np
from scipy import ndimage, linalg
import matplotlib.pyplot as plt
import multiprocessing
from pysar.utils import (ptime,
                         readfile,
                         writefile,
                         network as pnet)
from pysar.objects import (geometryDatasetNames,
                           geometry,
                           ifgramDatasetNames,
                           ifgramStack,
                           timeseries,
                           deramp)


###############################################################################
def read_timeseries_lalo(lat, lon, ts_file, lookup_file=None, ref_lat=None, ref_lon=None):
    """ Read time-series of one pixel with input lat/lon
    Parameters: lat/lon     : float, latitude/longitude
                ts_file     : string, filename of time-series HDF5 file
                lookup_file : string, filename of lookup table file
                ref_lat/lon : float, latitude/longitude of reference pixel
    Returns:    dates : 1D np.array of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis   : 1D np.array of float in meter
    """
    # read date
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    dates = ptime.date_list2vector(obj.dateList)[0]
    dates = np.array(dates)

    # read displacement
    coord = coordinate(obj.metadata, lookup_file=lookup_file)
    y, x = coord.geo2radar(lat, lon)[0:2]
    box = (x, y, x+1, y+1)
    dis = readfile.read(ts_file, box=box)[0]
    # reference pixel
    if ref_lat is not None:
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]
    #start at zero
    dis -= dis[0]
    return dates, dis


def read_timeseries_yx(y, x, ts_file, lookup_file=None, ref_y=None, ref_x=None):
    """ Read time-series of one pixel with input y/x
    Parameters: y/x         : int, row/column number of interest
                ts_file     : string, filename of time-series HDF5 file
                lookup_file : string, filename of lookup table file
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
    box = (x, y, x+1, y+1)
    dis = readfile.read(ts_file, box=box)[0]
    # reference pixel
    if ref_y is not None:
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]
    #start at zero
    dis -= dis[0]
    return dates, dis


###############################################################################
def get_all_conn_components(mask_in, min_num_pixel=1e4):
    """Get all connected component with number of pixels larger than threshold
    Parameters: mask_in  : 2D np.array with zero as background and non-zero as foreground
                min_num_pixel : int/float, minimum number of pixels to be identified and output
    Returns:    mask_out : list of 2D np.array in np.bool_ format
    """
    mask_in = np.array(mask_in)
    mask_out = []  # list of 2D np.array in bool
    mask_cc = get_largest_conn_component(mask_in, min_num_pixel=1e4)
    while not np.all(~mask_cc):
        mask_out.append(mask_cc)
        mask_in ^= mask_cc
        mask_cc = get_largest_conn_component(mask_in, min_num_pixel=1e4)
    return mask_out


def get_largest_conn_component(mask_in, min_num_pixel=1e4, display=False):
    """Extract the largest connected component from an 2D array
       with zero as background value
    Parameters: mask_in  : 2D np.array with zero as background and non-zero as foreground
                min_num_pixel : int/float, minimum number of pixels to be identified and output
                display : bool, display the result or not.
    Returns:    mask_out : 2D np.array in np.bool_ format
    """
    mask_out = np.zeros(mask_in.shape, np.bool_)
    labels, n_features = ndimage.label(mask_in)
    num_pixel = np.max(np.bincount(labels.flatten())[1:])
    if num_pixel < min_num_pixel:
        return mask_out

    max_label = np.argmax(np.bincount(labels.flatten())[1:]) + 1
    mask_out = labels == max_label
    if display:
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 5])
        ax[0].imshow(mask_in)
        ax[1].imshow(mask_out)
        ax[2].imshow(mask_in ^ mask_out)
        plt.show()
    return mask_out


def min_region_distance(mask1, mask2, display=False):
    """Calculate the min distance between two regions of pixels marked by mask1 and mask2
    Parameters: mask1/2 : 2D np.array in size of (length, width) in np.bool_ format
    Returns:    pts1 : tuple of 2 int, bridge point in mask1, in (x, y)
                pts2 : tuple of 2 int, bridge point in mask2, in (x, y)
                min_dist : float, min euclidean distance
    """
    from scipy.spatial import cKDTree
    y, x = np.where(mask1 != 0)
    pts1 = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    tree = cKDTree(pts1)

    y, x = np.where(mask2 != 0)
    pts2 = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    dist, idx = tree.query(pts2)

    idx_min = np.argmin(dist)
    xy2 = pts2[idx_min]
    xy1 = pts1[idx[idx_min]]
    min_dist = dist[idx_min]

    if display:
        plt.figure()
        plt.imshow(mask1 * 1 + mask2 * 2)
        plt.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], '-o')
        plt.show()

    return xy1, xy2, min_dist


def median_abs_deviation_threshold(data, center=None, cutoff=3.):
    """calculate rms_threshold based on the standardised residual
    outlier detection with median absolute deviation.
    https://www.statsmodels.org/dev/generated/statsmodels.robust.scale.mad.html
    """
    data = np.array(data)
    if center is None:
        center = np.median(data)
    #from statsmodels.robust import mad
    #data_mad = mad(data, center=center)
    data_mad = np.median(np.abs(data - center)) / 0.67448975019608171
    threshold = center + cutoff * data_mad
    return threshold


def wrap(data_in, wrap_range=[-1.*np.pi, np.pi]):
    """Wrap data into a range.
    Parameters: data_in    : np.array, array to be wrapped
                wrap_range : list of 2 float, range to be wrapped into
    Returns:    data       : np.array, data after wrapping
    """
    w0, w1 = wrap_range
    data = np.array(data_in)
    data = w0 + np.mod(data - w0, w1 - w0)
    return data


def get_snwe(metadata):
    lat0 = float(metadata['Y_FIRST'])
    lon0 = float(metadata['X_FIRST'])
    lat_step = float(metadata['Y_STEP'])
    lon_step = float(metadata['X_STEP'])
    length = int(metadata['LENGTH'])
    width = int(metadata['WIDTH'])
    lat1 = lat0 + lat_step * length
    lon1 = lon0 + lon_step * width
    SNWE = (lat1, lat0, lon0, lon1)
    return SNWE


def azimuth2heading_angle(az_angle):
    """Convert azimuth angle from ISCE los.rdr band2 into satellite orbit heading angle
    ISCE los.band2 is azimuth angle of LOS vector from ground target to the satellite 
    measured from the north in anti-clockwise as positive
    """
    head_angle = -1 * (180 + az_angle + 90)
    head_angle -= np.round(head_angle / 360.) * 360.
    return head_angle


def enu2los(e, n, u, inc_angle=34., head_angle=-168.):
    """
    Parameters: e : np.array or float, displacement in east-west direction, east as positive
                n : np.array or float, displacement in north-south direction, north as positive
                u : np.array or float, displacement in vertical direction, up as positive
                inc_angle  : np.array or float, local incidence angle from vertical
                head_angle : np.array or float, satellite orbit from the north in clock-wise direction as positive
    For AlosA: inc_angle = 34, head_angle = -12.873
    For AlosD: inc_angle = 34, head_angle = -167.157
    For SenD: inc_angle = 34, head_angle = -168
    """
    # if input angle is azimuth angle
    if (head_angle + 180.) > 45.:
        head_angle = azimuth2heading_angle(head_angle)

    inc_angle *= np.pi/180.
    head_angle *= np.pi/180.
    v_los = (-1 * e * np.cos(head_angle) * np.sin(inc_angle)
             + n * np.sin(head_angle) * np.sin(inc_angle)
             + u * np.cos(inc_angle))
    return v_los


def yes_or_no(question):
    """garrettdreyfus on Github: https://gist.github.com/garrettdreyfus/8153571"""
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    elif reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhhh... please enter ")


def interpolate_data(inData, outShape, interpMethod='linear'):
    """Interpolate input 2D matrix into different shape.
    Used to get full resolution perp baseline from ISCE coarse grid baseline file.
    Parameters: inData : 2D array
                outShape : tuple of 2 int in (length, width)
                interpMethod : string, choose in [nearest, linear, cubic]
    Returns:    outData : 2D array in outShape
    """
    from scipy.interpolate import RegularGridInterpolator as RGI
    inShape = inData.shape
    inPts = (np.arange(inShape[0]), np.arange(inShape[1]))
    xx, yy = np.meshgrid(np.linspace(0, inShape[1]-1, outShape[1], endpoint=False),
                         np.linspace(0, inShape[0]-1, outShape[0], endpoint=False))
    outPts = np.hstack((yy.reshape(-1, 1), xx.reshape(-1, 1)))
    outData = RGI(inPts, inData, method=interpMethod,
                  bounds_error=False)(outPts).reshape(outShape)
    return outData


def standardize_trop_model(tropModel, standardWeatherModelNames):
    tropModel = tropModel.replace('-', '').upper()
    if tropModel in standardWeatherModelNames.keys():
        tropModel = standardWeatherModelNames[tropModel]
    return tropModel


def subset_attribute(atr_dict, subset_box, print_msg=True):
    """Update attributes dictionary due to subset
    Parameters: atr_dict : dict, data attributes to update
                subset_box : 4-tuple of int, subset box defined in (x0, y0, x1, y1)
    Returns:    atr : dict, updated data attributes
    """
    if subset_box is None:
        return atr_dict

    sub_x = [subset_box[0], subset_box[2]]
    sub_y = [subset_box[1], subset_box[3]]

    # Update attribute variable
    atr = dict(atr_dict)
    atr['LENGTH'] = str(sub_y[1]-sub_y[0])
    atr['WIDTH'] = str(sub_x[1]-sub_x[0])
    atr['YMAX'] = str(sub_y[1]-sub_y[0] - 1)
    atr['XMAX'] = str(sub_x[1]-sub_x[0] - 1)
    if print_msg:
        print('update LENGTH, WIDTH, Y/XMAX')

    # Subset atribute
    atr['SUBSET_YMAX'] = str(sub_y[1] + int(atr_dict.get('SUBSET_YMIN', '0')))
    atr['SUBSET_YMIN'] = str(sub_y[0] + int(atr_dict.get('SUBSET_YMIN', '0')))
    atr['SUBSET_XMAX'] = str(sub_x[1] + int(atr_dict.get('SUBSET_XMIN', '0')))
    atr['SUBSET_XMIN'] = str(sub_x[0] + int(atr_dict.get('SUBSET_XMIN', '0')))
    if print_msg:
        print(('update/add SUBSET_XMIN/YMIN/XMAX/YMAX: '
               '{x0}/{y0}/{x1}/{y1}').format(x0=atr['SUBSET_XMIN'],
                                             y0=atr['SUBSET_YMIN'],
                                             x1=atr['SUBSET_XMAX'],
                                             y1=atr['SUBSET_YMAX']))

    # Geo coord
    if 'Y_FIRST' in atr.keys():
        atr['Y_FIRST'] = str(float(atr['Y_FIRST']) + sub_y[0]*float(atr['Y_STEP']))
        atr['X_FIRST'] = str(float(atr['X_FIRST']) + sub_x[0]*float(atr['X_STEP']))
        if print_msg:
            print('update Y/X_FIRST')

    # Reference in space
    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str(int(atr['REF_Y']) - sub_y[0])
        atr['REF_X'] = str(int(atr['REF_X']) - sub_x[0])
        if print_msg:
            print('update REF_Y/X')

    # Starting Range for file in radar coord
    if not 'Y_FIRST' in atr_dict.keys():
        try:
            atr['STARTING_RANGE'] = float(atr['STARTING_RANGE'])
            atr['STARTING_RANGE'] += float(atr['RANGE_PIXEL_SIZE'])*sub_x[0]
            if print_msg:
                print('update STARTING_RANGE')
        except:
            pass

    return atr


def ceil_to_1(x):
    """Return the most significant digit of input number and ceiling it"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)+10**digit


def round_to_1(x):
    """Return the most significant digit of input number"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -1*digit)


def touch(fname_list, times=None):
    """python equivalent function to Unix utily - touch
    It sets the modification and access times of files to the current time of day.
    If the file doesn't exist, it is created with default permissions.
    Inputs/Output:
        fname_list - string / list of string
    """
    if not fname_list:
        return None

    if isinstance(fname_list, str):
        fname_list = [fname_list]

    fname_list = [x for x in fname_list if x != None]
    for fname in fname_list:
        if os.path.isfile(fname):
            with open(fname, 'a'):
                os.utime(fname, times)
                print('touch '+fname)

    if len(fname_list) == 1:
        fname_list = fname_list[0]
    return fname_list


def get_lookup_file(filePattern=None, abspath=False, print_msg=True):
    """Find lookup table file with/without input file pattern"""
    # Search Existing Files
    if not filePattern:
        filePattern = ['geometryRadar.h5',
                       'geometryGeo_tight.h5', 'geometryGeo.h5',
                       'geomap*lks_tight.trans', 'geomap*lks.trans',
                       'sim*_tight.UTM_TO_RDC', 'sim*.UTM_TO_RDC']
        filePattern = [os.path.join('INPUTS', i) for i in filePattern]
    existFiles = []
    try:
        existFiles = get_file_list(filePattern)
    except:
        if print_msg:
            print('ERROR: No geometry / lookup table file found!')
            print('It should be like:')
            print(filePattern)
        return None

    # Check Files Info
    outFile = None
    for fname in existFiles:
        atr = readfile.read_attribute(fname)
        for dsName in ['rangeCoord', 'longitude']:
            try:
                dset = readfile.read(fname, datasetName=dsName, print_msg=False)[0]
                outFile = fname
                break
            except:
                pass

    if not outFile:
        if print_msg:
            print('No lookup table info range/lat found in files.')
        return None

    # Path Format
    if abspath:
        outFile = os.path.abspath(outFile)
    return outFile


def get_geometry_file(dset, coordType=None, filePattern=None, abspath=False, print_msg=True):
    """Find geometry file containing input specific dataset"""
    if dset not in geometryDatasetNames:
        raise ValueError('Unrecognized geometry dataset name: %s' % (dset))

    # Search Existing Files
    if not filePattern:
        filePattern = ['geometryRadar.h5', 'geometryGeo.h5']
        if dset in ['rangeCoord', 'azimuthCoord']:
            filePattern += ['geomap*lks_tight.trans', 'geomap*lks.trans',
                            'sim*_tight.UTM_TO_RDC', 'sim*.UTM_TO_RDC']
        elif dset == 'height':
            filePattern += ['demRadar.h5', 'demGeo.h5',
                            'radar*.hgt', '*.dem',
                            '*.dem.wgs84', '*.hgt_sim']
        elif dset == 'incidenceAngle':
            filePattern += ['*incidenceAngle.h5']
        elif dset == 'slantRangeDistance':
            filePattern += ['*rangeDistance.h5']
        elif dset == 'waterMask':
            filePattern += ['*waterMask.h5']
        elif dset == 'shadowMask':
            filePattern += ['*shadowMask.h5']

    existFiles = []
    try:
        existFiles = get_file_list(filePattern)
    except:
        if print_msg:
            print('ERROR: No %s file found!' % (dset))
            print('It should be like:')
            print(filePattern)
        return None

    # Check Files Info
    outFile = None
    for fname in existFiles:
        # Check coord type
        if coordType:
            atr = readfile.read_attribute(fname)
            if ((coordType == 'radar' and 'Y_FIRST' in atr.keys()) or
                    (coordType == 'geo' and 'Y_FIRST' not in atr.keys())):
                continue
        # Check dataset
        try:
            dset = readfile.read(fname, datasetName=dset, print_msg=False)[0]
            outFile = fname
            break
        except:
            pass
    if not outFile:
        if print_msg:
            print('No %s info found in files.' % (dset))
        return None

    # Path Format
    if abspath:
        outFile = os.path.abspath(outFile)
    return outFile


def check_loaded_dataset(work_dir='./', inps=None, print_msg=True):
    """Check the result of loading data for the following two rules:
        1. file existance
        2. file attribute readability

    If inps is valid/not_empty: return updated inps;
    Otherwise, return True/False if all recommended file are loaded and readably or not

    Parameters: work_dir : string, PySAR working directory
                inps : Namespace, optional, variable for pysarApp.py.
                    Not needed for check loading result.
    Returns:    loadComplete : bool, complete loading or not
                    #if True, PROCESS, SLC folder could be removed.
                or  inps : Namespace, if it's inputed
                    atr : dict, metadata of found ifgramStack file
    Example:    True = ut.check_loaded_dataset($SCRATCHDIR+'/SinabungT495F50AlosA/PYSAR')
                inps, atr = ut.check_loaded_dataset(inps.workDir, inps)
    """
    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    if inps:
        inps.stackFile = None
        inps.geomFile = None
        inps.lookupFile = None

    #--------------------------- Search file ------------------------#
    # 1. interferograms stack
    file_list = [os.path.join(work_dir, 'INPUTS/ifgramStack.h5')]
    stack_file = is_file_exist(file_list, abspath=True)
    if not stack_file:
        if inps:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    './INPUTS/ifgramStack.h5')
            #return inps, None
        else:
            return False
    else:
        atr = readfile.read_attribute(stack_file)

    # 2. geometry
    if 'X_FIRST' in atr.keys():
        geocoded = True
        file_list = [os.path.join(work_dir, 'INPUTS/geometryGeo.h5')]
    else:
        geocoded = False
        file_list = [os.path.join(work_dir, 'INPUTS/geometryRadar.h5')]
    geom_file = is_file_exist(file_list, abspath=True)

    # 3. lookup table
    # could be different than geometry file in case of roipac and gamma
    file_list = [os.path.join(work_dir, 'INPUTS/geometry*.h5')]
    lookup_file = get_lookup_file(file_list,
                                  abspath=True,
                                  print_msg=print_msg)

    #------------------ Check required datasets -----------------------#
    # set loadComplete to False if any required dataset is missing
    loadComplete = True

    # 1. stack_file: unwrapPhase, coherence
    if stack_file is not None:
        stack_obj = ifgramStack(stack_file)
        stack_obj.open(print_msg=False)
        for dsName in ['unwrapPhase', 'coherence']:
            if dsName not in stack_obj.datasetNames:
                loadComplete = False
                raise Exception(('required dataset "{}" is missing'
                                 ' in file {}'.format(dsName, stack_file)))
    else:
        loadComplete = False
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                './INPUTS/ifgramStack.h5')

    # 2. geom_file: height
    if geom_file is not None:
        geom_obj = geometry(geom_file)
        geom_obj.open(print_msg=False)
        dsName = geometryDatasetNames[0]
        if dsName not in geom_obj.datasetNames:
            loadComplete = False
            raise Exception(('required dataset "{}" is missing'
                             ' in file {}'.format(dsName, geom_file)))
    else:
        loadComplete = False
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                './INPUTS/geometry*.h5')

    # 3. lookup_file: latitude,longitude or rangeCoord,azimuthCoord
    if lookup_file is not None:
        lut_obj = geometry(lookup_file)
        lut_obj.open(print_msg=False)

        if atr['PROCESSOR'] in ['isce', 'doris']:
            dsNames = [geometryDatasetNames[1],
                       geometryDatasetNames[2]]
        elif atr['PROCESSOR'] in ['gamma', 'roipac']:
            dsNames = [geometryDatasetNames[3],
                       geometryDatasetNames[4]]
        else:
            raise AttributeError('InSAR processor: {}'.format(atr['PROCESSOR']))

        for dsName in dsNames:
            if dsName not in lut_obj.datasetNames:
                loadComplete = False
                raise Exception(('required dataset "{}" is missing'
                                 ' in file {}'.format(dsName, lookup_file)))
    else:
        loadComplete = False
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                './INPUTS/geometry*.h5')

    # print message
    if print_msg:
        print(('Loaded dataset are processed by '
               'InSAR software: {}'.format(atr['PROCESSOR'])))
        if geocoded:
            print('Loaded dataset is in GEO coordinates')
        else:
            print('Loaded dataset is in RADAR coordinates')
        print('Interferograms Stack: {}'.format(stack_file))
        print('Geometry File       : {}'.format(geom_file))
        print('Lookup Table File   : {}'.format(lookup_file))
        if loadComplete:
            print(('-'*50+'\nAll data needed found/loaded/copied. '
                   'Processed 2-pass InSAR data can be removed.'))
        print('-'*50)

    # Return & Update namespace inps if inputed
    if inps:
        inps.stackFile = stack_file
        inps.geomFile = geom_file
        inps.lookupFile = lookup_file
        inps.geocoded = geocoded
        return inps, atr
    else:
        return loadComplete


def is_file_exist(file_list, abspath=True):
    """Check if any file in the file list 1) exists and 2) readable
    Parameters: file_list : list of string, file name with/without wildcards
                abspath   : bool, return absolute file name/path or not
    Returns:    file_path : string, found file name/path; None if not.
    """
    try:
        file = get_file_list(file_list, abspath=abspath)[0]
        atr = readfile.read_attribute(file)
    except:
        file = None
    return file


def four_corners(atr):
    """Return 4 corners lat/lon"""
    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    west = float(atr['X_FIRST'])
    north = float(atr['Y_FIRST'])
    south = north + lat_step*length
    east = west + lon_step*width

    return west, east, south, north


def get_circular_mask(x, y, radius, shape):
    """Get mask of pixels within circle defined by (x, y, r)"""
    length, width = shape
    yy, xx = np.ogrid[-y:length-y,
                      -x:width-x]
    cmask = (xx**2 + yy**2 <= radius**2)
    return cmask


def circle_index(atr, circle_par):
    """Return Index of Elements within a Circle centered at input pixel
    Parameters: atr : dictionary
                    containging the following attributes:
                    WIDT
                    LENGTH
                circle_par : string in the format of 'y,x,radius'
                    i.e. '200,300,20'          for radar coord
                         '31.0214,130.5699,20' for geo   coord
    Returns:    idx : 2D np.array in bool type
                    mask matrix for those pixel falling into the circle
                    defined by circle_par
    Examples:   idx_mat = ut.circle_index(atr, '200,300,20')
                idx_mat = ut.circle_index(atr, '31.0214,130.5699,20')
    """

    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])

    if type(circle_par) == tuple:
        cir_par = circle_par
    elif type(circle_par) == list:
        cir_par = circle_par
    else:
        cir_par = circle_par.replace(',', ' ').split()
    cir_par = [str(i) for i in cir_par]

    try:
        c_y = int(cir_par[0])
        c_x = int(cir_par[1])
        radius = int(float(cir_par[2]))
        print('Input circle index in y/x coord: %d, %d, %d' % (c_y, c_x, radius))
    except:
        try:
            c_lat = float(cir_par[0])
            c_lon = float(cir_par[1])
            radius = int(float(cir_par[2]))
            c_y = np.rint((c_lat-float(atr['Y_FIRST'])) / float(atr['Y_STEP']))
            c_x = np.rint((c_lon-float(atr['X_FIRST'])) / float(atr['X_STEP']))
            print(('Input circle index in lat/lon coord: '
                   '{:.4f}, {:.4f}, {}'.format(c_lat, c_lon, radius)))
        except:
            print('\nERROR: Unrecognized circle index format: '+circle_par)
            print('Supported format:')
            print('--circle 200,300,20            for radar coord input')
            print('--circle 31.0214,130.5699,20   for geo   coord input\n')
            return 0

    y, x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx


def check_template_auto_value(templateDict, auto_file='../defaults/pysarApp.cfg'):
    """Replace auto value based on $PYSAR_HOME/pysar/defaults/template.cfg file."""
    # Read default template value and turn yes/no to True/False
    templateAutoFile = os.path.join(os.path.dirname(__file__), auto_file)
    templateAutoDict = readfile.read_template(templateAutoFile)

    # Update auto value of input template dict
    for key, value in templateDict.items():
        if value == 'auto' and key in templateAutoDict.keys():
            templateDict[key] = templateAutoDict[key]

    # Change yes --> True and no --> False
    specialValues = {'yes': True,
                     'no': False,
                     'none': None
                     }
    for key, value in templateDict.items():
        if value in specialValues.keys():
            templateDict[key] = specialValues[value]

    return templateDict


def update_template_file(template_file, extra_dict):
    """Update option value in template_file with value from input extra_dict"""
    # Compare and skip updating template_file if no new option value found.
    update = False
    orig_dict = readfile.read_template(template_file)
    for key, value in orig_dict.items():
        if key in extra_dict.keys() and extra_dict[key] != value:
            update = True
    if not update:
        print('No new option value found, skip updating '+template_file)
        return template_file

    # Update template_file with new value from extra_dict
    tmp_file = template_file+'.tmp'
    f_tmp = open(tmp_file, 'w')
    for line in open(template_file, 'r'):
        c = [i.strip() for i in line.strip().split('=', 1)]
        if not line.startswith(('%', '#')) and len(c) > 1:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].strip()
            if key in extra_dict.keys() and extra_dict[key] != value:
                line = line.replace(value, extra_dict[key], 1)
                print('    {}: {} --> {}'.format(key, value, extra_dict[key]))
        f_tmp.write(line)
    f_tmp.close()

    # Overwrite exsting original template file
    mvCmd = 'mv {} {}'.format(tmp_file, template_file)
    os.system(mvCmd)
    return template_file


def get_residual_std(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    """Calculate deramped standard deviation in space for each epoch of input timeseries file.
    Parameters: timeseries_resid_file - string, timeseries HDF5 file,
                    e.g. timeseries_ECMWF_demErrInvResid.h5
                mask_file - string, mask file, e.g. maskTempCoh.h5
                ramp_type - string, ramp type, e.g. linear, quadratic, no for do not remove ramp
    Returns:    std_list  - list of float, standard deviation of deramped input timeseries file
                date_list - list of string in YYYYMMDD format, corresponding dates
    Example:    import pysar.utils.utils as ut
                std_list, date_list = ut.get_residual_std('timeseries_ECMWF_demErrInvResid.h5',
                                                          'maskTempCoh.h5')
    """
    # Intermediate files name
    if ramp_type == 'no':
        print('No ramp removal')
        deramped_file = timeseries_resid_file
    else:
        deramped_file = '{}_ramp.h5'.format(os.path.splitext(timeseries_resid_file)[0])
    std_file = os.path.splitext(deramped_file)[0]+'_std.txt'

    # Get residual std text file
    if run_or_skip(out_file=std_file, in_file=[deramped_file, mask_file], check_readable=False) == 'run':
        if run_or_skip(out_file=deramped_file, in_file=timeseries_resid_file) == 'run':
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print('removing a {} ramp from file: '.format(ramp_type, timeseries_resid_file))
                deramped_file = run_deramp(timeseries_resid_file,
                                           ramp_type=ramp_type,
                                           mask_file=mask_file,
                                           out_file=deramped_file)
        print('calculating residual standard deviation for each epoch from file: '+deramped_file)
        std_file = timeseries(deramped_file).timeseries_std(maskFile=mask_file, outFile=std_file)

    # Read residual std text file
    print('read timeseries RMS from file: '+std_file)
    std_fileContent = np.loadtxt(std_file, dtype=bytes).astype(str)
    std_list = std_fileContent[:, 1].astype(np.float).tolist()
    date_list = list(std_fileContent[:, 0])
    return std_list, date_list


def get_residual_rms(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    """Calculate deramped Root Mean Square in space for each epoch of input timeseries file.
    Parameters: timeseries_resid_file : string, 
                    timeseries HDF5 file, e.g. timeseries_ECMWF_demErrInvResid.h5
                mask_file : string,
                    mask file, e.g. maskTempCoh.h5
                ramp_type : string, 
                    ramp type, e.g. linear, quadratic, no for do not remove ramp
    Returns:    rms_list : list of float,
                    Root Mean Square of deramped input timeseries file
                date_list : list of string in YYYYMMDD format,
                    corresponding dates
                rms_file : string, text file with rms and date info.
    Example:
        import pysar.utils.utils as ut
        rms_list, date_list = ut.get_residual_rms('timeseriesResidual.h5', 'maskTempCoh.h5')
    """
    # Intermediate files name
    if ramp_type == 'no':
        print('No ramp removal')
        deramped_file = timeseries_resid_file
    else:
        deramped_file = '{}_ramp.h5'.format(os.path.splitext(timeseries_resid_file)[0])
    rms_file = os.path.join(os.path.dirname(os.path.abspath(deramped_file)),
                            'rms_{}.txt'.format(os.path.splitext(deramped_file)[0]))

    # Get residual RMS text file
    if run_or_skip(out_file=rms_file, in_file=[deramped_file, mask_file], check_readable=False) == 'run':
        if run_or_skip(out_file=deramped_file, in_file=timeseries_resid_file) == 'run':
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print('removing a {} ramp from file: {}'.format(ramp_type, timeseries_resid_file))
                deramped_file = run_deramp(timeseries_resid_file,
                                           ramp_type=ramp_type,
                                           mask_file=mask_file,
                                           out_file=deramped_file)
        print('Calculating residual RMS for each epoch from file: '+deramped_file)
        rms_file = timeseries(deramped_file).timeseries_rms(maskFile=mask_file, outFile=rms_file)

    # Read residual RMS text file
    print('read timeseries residual RMS from file: '+rms_file)
    rms_fileContent = np.loadtxt(rms_file, dtype=bytes).astype(str)
    rms_list = rms_fileContent[:, 1].astype(np.float).tolist()
    date_list = list(rms_fileContent[:, 0])
    return rms_list, date_list, rms_file


def normalize_timeseries(ts_mat, nanValue=0):
    """Normalize timeseries of 2D matrix in time domain"""
    ts_mat -= np.min(ts_mat, 0)
    ts_mat *= 1/np.max(ts_mat, 0)
    ts_mat[np.isnan(ts_mat)] = 0
    return ts_mat


def normalize_timeseries_old(ts_mat, nanValue=0):
    ts_mat -= np.max(ts_mat, 0)
    ts_mat *= -1
    ts_mat /= np.max(ts_mat, 0)
    ts_mat[np.isnan(ts_mat)] = 1
    return ts_mat


############################################################
def run_or_skip(out_file, in_file=None, check_readable=True, print_msg=True):
    """Check whether to update out_file or not.
    return run if any of the following meets:
        1. out_file is empty, e.g. None, []
        2. out_file is not existed
        3. out_file is not readable by readfile.read_attribute() when check_readable=True
        4. out_file is older than in_file, if in_file is not None
    Otherwise, return skip.

    If in_file=None and out_file exists and readable, return skip

    Parameters: out_file : string or list of string, output file(s)
                in_file  : string or list of string, input file(s)
                check_readable : bool, check if the 1st output file has attribute 'WIDTH'
                print_msg      : bool, print message
    Returns:    run/skip : str, whether to update output file or not
    Example:    if ut.run_or_skip(out_file='timeseries_ECMWF_demErr.h5', in_file='timeseries_ECMWF.h5'):
                if ut.run_or_skip(out_file='exclude_date.txt',
                                  in_file=['timeseries_ECMWF_demErrInvResid.h5',
                                           'maskTempCoh.h5',
                                           'pysar_template.txt'],  
                                  check_readable=False):
    """
    # 1 - check existance of output files
    if not out_file:
        return 'run'
    else:
        if isinstance(out_file, str):
            out_file = [out_file]
        if not all(os.path.isfile(i) for i in out_file):
            return 'run'

    # 2 - check readability of output files
    if check_readable:
        try:
            atr = readfile.read_attribute(out_file[0])
            width = atr['WIDTH']
        except:
            if print_msg:
                print('{} exists, but can not read, remove it.'.format(out_file[0]))
            rmCmd = 'rm {}'.format(out_file[0])
            print(rmCmd)
            os.system(rmCmd)
            return 'run'

    # 3 - check modification time of output and input files
    if in_file:
        in_file = get_file_list(in_file)
        # Check modification time
        if in_file:
            t_in  = max([os.path.getmtime(i) for i in in_file])
            t_out = min([os.path.getmtime(i) for i in out_file])
            if t_in > t_out:
                return 'run'
            elif print_msg:
                print('{} exists and is newer than {}, skip updating.'.format(out_file, in_file))
    return 'skip'


def update_attribute_or_not(atr_new, atr_orig):
    """Compare new attributes with exsiting ones"""
    update = False
    for key in atr_new.keys():
        value = str(atr_new[key])
        if ((key in atr_orig.keys() and value == str(atr_orig[key]) and value != 'None')
                or (key not in atr_orig.keys() and value == 'None')):
            next
        else:
            update = True
    return update


def add_attribute(File, atr_new=dict(), print_msg=False):
    """Add/update input attribute into File
    Parameters: File - string, path/name of file
                atr_new - dict, attributes to be added/updated
                    if value is None, delete the item from input File attributes
    Returns:    File - string, path/name of updated file
    """
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']

    # Compare new attributes with exsiting ones
    update = update_attribute_or_not(atr_new, atr)
    if not update:
        print(('All updated (removed) attributes already exists (do not exists)'
               ' and have the same value, skip update.'))
        return File

    # Update attributes
    f = h5py.File(File, 'r+')
    for key, value in iter(atr_new.items()):
        # delete the item is new value is None
        if value == 'None' or value is None:
            try:
                f.attrs.pop(key)
                if print_msg:
                    print('remove {}'.format(key))
            except:
                pass
        else:
            f.attrs[key] = str(value)
            if print_msg:
                print('{} = {}'.format(key, str(value)))
    f.close()
    return File


def check_parallel(file_num=1, print_msg=True, maxParallelNum=8):
    """Check parallel option based file num and installed module
    Examples:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(1000)
    """
    enable_parallel = True

    # Disable parallel option for one input file
    if file_num <= 1:
        enable_parallel = False
        if print_msg:
            print('parallel processing is diabled for one input file')
        return 1, enable_parallel, None, None

    # Check required python module
    try:
        from joblib import Parallel, delayed
    except:
        print('Can not import joblib')
        print('parallel is disabled.')
        enable_parallel = False
        return 1, enable_parallel, None, None

    # Find proper number of cores for parallel processing
    num_cores = min(multiprocessing.cpu_count(), file_num, maxParallelNum)
    if num_cores <= 1:
        enable_parallel = False
        print('parallel processing is disabled because min of the following two numbers <= 1:')
        print('available cpu number of the computer: {}'.format(multiprocessing.cpu_count()))
    elif print_msg:
        print('parallel processing using %d cores ...' % (num_cores))

    try:
        return num_cores, enable_parallel, Parallel, delayed
    except:
        return num_cores, enable_parallel, None, None


def perp_baseline_timeseries(atr, dimension=1):
    """Calculate perpendicular baseline for each acquisition within timeseries
    Parameters: atr - dict, including the following PySAR attribute
                    LENGTH
                    P_BASELINE_TIMESERIES
                    P_BASELINE_TOP_TIMESERIES (optional)
                    P_BASELINE_BOTTOM_TIMESERIES (optional)
                dimension - int, choices = [0, 1]
                    0 for constant P_BASELINE in azimuth direction
                    1 for linear P_BASELINE in azimuth direction, for radar coord only
    Returns:    pbase - np.array, with shape = [date_num, 1] or [date_num, length]
    """
    if dimension > 0 and 'Y_FIRST' in atr.keys():
        dimension = 0
        print('file is in geo coordinates, return constant P_BASELINE for one interferogram')
    if dimension > 0 and any(i not in atr.keys() for i in ['P_BASELINE_TOP_TIMESERIES',
                                                           'P_BASELINE_BOTTOM_TIMESERIES']):
        dimension = 0
        print(('No P_BASELINE_TOP/BOTTOM_TIMESERIES attributes found,'
               ' return constant P_BASELINE for one interferogram'))

    pbase_center = np.array(
        [float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]).reshape(-1, 1)
    if dimension == 0:
        pbase = pbase_center
    elif dimension == 1:
        pbase_top = [float(i) for i in atr['P_BASELINE_TOP_TIMESERIES'].split()]
        pbase_bottom = [float(i) for i in atr['P_BASELINE_BOTTOM_TIMESERIES'].split()]
        pbase_top = np.array(pbase_top).reshape(-1, 1)
        pbase_bottom = np.array(pbase_bottom).reshape(-1, 1)
        length = int(atr['LENGTH'])
        date_num = pbase_center.shape[0]
        pbase = np.zeros((date_num, length))
        for i in range(date_num):
            pbase[i, :] = np.linspace(pbase_top[i], pbase_bottom[i], num=length, endpoint='FALSE')
    else:
        raise ValueError('Input pbase dimension: %s, only support 0 and 1.' % (str(dimension)))

    return pbase


def range_distance(atr, dimension=2, print_msg=True):
    """Calculate slant range distance from input attribute dict
    Parameters: atr : dict, including the following ROI_PAC attributes:
                    STARTING_RANGE
                    RANGE_PIXEL_SIZE
                    LENGTH
                    WIDTH
                dimension : int, choices = [0,1,2]
                    2 for 2d matrix, vary in range direction, constant in az direction,
                        for radar coord only
                    1 for 1d matrix, in range direction, for radar coord file
                    0 for center value
    Returns:    np.array (0, 1 or 2 D) : range distance between antenna and ground target in meters
    """
    # return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        if print_msg:
            print('input file is geocoded, return center range distance for the whole area')

    range_n, dR = float(atr['STARTING_RANGE']), float(atr['RANGE_PIXEL_SIZE'])
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])

    range_f = range_n + dR*(width-1)
    range_c = (range_f + range_n)/2.0
    if print_msg:
        print('center range : %.2f m' % (range_c))
        print('near   range : %.2f m' % (range_n))
        print('far    range : %.2f m' % (range_f))

    if dimension == 0:
        return np.array(range_c, np.float32)

    range_x = np.linspace(range_n, range_f, num=width)
    if dimension == 1:
        return np.array(range_x, np.float32)
    else:
        range_xy = np.tile(range_x, (length, 1))
        return np.array(range_xy, np.float32)


def incidence_angle(atr, dem=None, dimension=2, print_msg=True):
    """Calculate 2D matrix of incidence angle from ROI_PAC attributes, very accurate.
    Parameters: atr : dict - ROI_PAC attributes including the following items:
                     STARTING_RANGE
                     RANGE_PIXEL_SIZE
                     EARTH_RADIUS
                     HEIGHT
                     LENGTH
                     WIDTH
                dem : 2D array for height to calculate local incidence angle
                dimension : int,
                            2 for 2d matrix
                            1 for 1d array
                            0 for one center value
                print_msg : bool
    Returns:    inc_angle : 2D np.array, incidence angle in degree for each pixel
    Example:    dem = readfile.read('hgt.rdr')[0]
                atr = readfile.read_attribute('filt_fine.unw')
                inc_angle = ut.incidence_angle(atr, dem=dem)
    """
    # Return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        if print_msg:
            print('input file is geocoded, return center incident angle only')

    # Read Attributes
    range_n = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    r = float(atr['EARTH_RADIUS'])
    H = float(atr['HEIGHT'])
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    # Calculation
    range_f = range_n+dR*width
    inc_angle_n = (np.pi - np.arccos((r**2 + range_n**2 - (r+H)**2)/(2*r*range_n))) * 180.0/np.pi
    inc_angle_f = (np.pi - np.arccos((r**2 + range_f**2 - (r+H)**2)/(2*r*range_f))) * 180.0/np.pi
    if print_msg:
        print('near   incidence angle : {:.4f} degree'.format(inc_angle_n))
        print('far    incidence angle : {:.4f} degree'.format(inc_angle_f))

    if dimension == 0:
        inc_angle = (inc_angle_n + inc_angle_f) / 2.0
        if print_msg:
            print('center incidence angle : {:.4f} degree'.format(inc_angle))

    elif dimension == 1:
        inc_angle = np.linspace(inc_angle_n, inc_angle_f, num=width,
                                endpoint='FALSE', dtype=np.float32)

    elif dimension == 2:
        # consider the local variable due to topography
        if dem is not None:
            range_dist = range_distance(atr, dimension=2, print_msg=False)
            inc_angle = (np.pi - np.arccos(((r+dem)**2 + range_dist**2 - (r+H)**2) / 
                                           (2*(r+dem)*range_dist))) * 180.0/np.pi
        else:
            inc_angle = np.tile(np.linspace(inc_angle_n, inc_angle_f, num=width,
                                            endpoint='FALSE', dtype=np.float32), (length, 1))
    else:
        raise Exception('un-supported dimension input: {}'.format(dimension))
    return inc_angle


def which(program):
    """Test if executable exists"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def check_drop_ifgram(h5, print_msg=True):
    """Update ifgram_list based on 'DROP_IFGRAM' attribute
    Parameters: h5 : HDF5 file object
    Returns:    dsListOut : list of string, group name with DROP_IFGRAM = 'yes'
    Example:    h5 = h5py.File('unwrapIfgram.h5','r')
                ifgram_list = ut.check_drop_ifgram(h5)
    """
    # Return all interferogram list if 'DROP_IFGRAM' do not exist
    k = list(h5.keys())[0]
    dsList = sorted(h5[k].keys())
    atr = h5[k][dsList[0]].attrs
    if 'DROP_IFGRAM' not in atr.keys():
        return dsList

    dsListOut = list(dsList)
    for ds in dsList:
        if ('DROP_IFGRAM' in h5[k][ds].attrs.keys()
                and h5[k][ds].attrs['DROP_IFGRAM'] == 'yes'):
            dsListOut.remove(ds)

    if len(dsList) > len(dsListOut) and print_msg:
        print("remove interferograms with 'DROP_IFGRAM'='yes'")
    return dsListOut


def nonzero_mask(File, out_file='mask.h5', datasetName=None):
    """Generate mask file for non-zero value of input multi-group hdf5 file"""
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    if k == 'ifgramStack':
        mask = ifgramStack(File).nonzero_mask(datasetName=datasetName)
    else:
        print('Only ifgramStack file is supported for now, input is '+k)
        return None

    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, out_file=out_file, metadata=atr)
    return out_file


######################################################################################################
def spatial_average(File, datasetName='coherence', maskFile=None, box=None,
                    saveList=False, checkAoi=True):
    """Read/Calculate Spatial Average of input file.

    If input file is text file, read it directly;
    If input file is data matrix file:
        If corresponding text file exists with the same mask file/AOI info, read it directly;
        Otherwise, calculate it from data file.

        Only non-nan pixel is considered.
    Parameters: File     : string, path of input file
                maskFile : string, path of mask file, e.g. maskTempCoh.h5
                box      : 4-tuple defining the left, upper, right, and lower pixel coordinate
                saveList : bool, save (list of) mean value into text file
    Returns:    meanList : list for float, average value in space for each epoch of input file
                dateList : list of string for date info
                    date12_list, e.g. 101120-110220, for interferograms/coherence
                    date8_list, e.g. 20101120, for timeseries
                    file name, e.g. velocity.h5, for all the other file types
    Example:    meanList = spatial_average('coherence.h5')[0]
                meanList, date12_list = spatial_average('coherence.h5',
                                                        'maskTempCoh.h5',
                                                        saveList=True)
    """
    def read_text_file(fname):
        txtContent = np.loadtxt(fname, dtype=bytes).astype(str)
        meanList = [float(i) for i in txtContent[:, 1]]
        dateList = [i for i in txtContent[:, 0]]
        return meanList, dateList

    # Baic File Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    if not box:
        box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))

    # If input is text file
    suffix = ''
    if k == 'ifgramStack':
        suffix += '_'+datasetName
    suffix += '_spatialAvg.txt'
    if File.endswith(suffix):
        print('Input file is spatial average txt already, read it directly')
        meanList, dateList = read_text_file(File)
        return meanList, dateList

    # Read existing txt file only if 1) data file is older AND 2) same AOI
    txtFile = os.path.splitext(os.path.basename(File))[0]+suffix
    file_line = '# Data file: {}\n'.format(os.path.basename(File))
    mask_line = '# Mask file: {}\n'.format(maskFile)
    aoi_line = '# AOI box: {}\n'.format(box)
    try:
        # Read AOI line from existing txt file
        fl = open(txtFile, 'r')
        lines = fl.readlines()
        fl.close()
        if checkAoi:
            try:
                aoi_line_orig = [i for i in lines if '# AOI box:' in i][0]
            except:
                aoi_line_orig = ''
        else:
            aoi_line_orig = aoi_line
        try:
            mask_line_orig = [i for i in lines if '# Mask file:' in i][0]
        except:
            mask_line_orig = ''
        if (aoi_line_orig == aoi_line 
                and mask_line_orig == mask_line
                and run_or_skip(out_file=txtFile,
                                in_file=[File, maskFile],
                                check_readable=False) == 'skip'):
            print(txtFile+' already exists, read it directly')
            meanList, dateList = read_text_file(txtFile)
            return meanList, dateList
    except:
        pass

    # Calculate mean coherence list
    if k == 'ifgramStack':
        obj = ifgramStack(File)
        obj.open(print_msg=False)
        meanList, dateList = obj.spatial_average(datasetName=datasetName,
                                                 maskFile=maskFile,
                                                 box=box)
        pbase = obj.pbaseIfgram
        tbase = obj.tbaseIfgram
        obj.close()
    elif k == 'timeseries':
        meanList, dateList = timeseries(File).spatial_average(maskFile=maskFile,
                                                              box=box)
    else:
        data = readfile.read(File, box=box)[0]
        if maskFile and os.path.isfile(maskFile):
            print('mask from file: '+maskFile)
            mask = readfile.read(maskFile, datasetName='mask', box=box)[0]
            data[mask == 0.] = np.nan
        meanList = np.nanmean(data)
        dateList = [os.path.basename(File)]

    # Write mean coherence list into text file
    if saveList:
        print('write average value in space into text file: '+txtFile)
        fl = open(txtFile, 'w')
        # Write comments
        fl.write(file_line+mask_line+aoi_line)
        # Write data list
        numLine = len(dateList)
        if k == 'ifgramStack':
            fl.write('#\tDATE12\t\tMean\tBtemp/days\tBperp/m\t\tNum\n')
            for i in range(numLine):
                fl.write('%s\t%.4f\t%8.0f\t%8.1f\t%d\n' %
                         (dateList[i], meanList[i], tbase[i], pbase[i], i))
        else:
            fl.write('#\tDATE12\t\tMean\n')
            for i in range(numLine):
                fl.write('%s\t%.4f\n' % (dateList[i], meanList[i]))
        fl.close()

    if len(meanList) == 1:
        meanList = meanList[0]
        dateList = dateList[0]
    return meanList, dateList


def temporal_average(File, datasetName='coherence', updateMode=False, outFile=None):
    """Calculate temporal average of multi-temporal dataset, equivalent to stacking
    For ifgramStakc/unwrapPhase, return average phase velocity

    Parameters: File : string, file to be averaged in time
                outFile : string, output file name
                datasetName : string, dataset to be read from input file, for multiple
                    datasets file - ifgramStack - only
                    e.g.: coherence, unwrapPhase
                ignoreNan: bool, ignore NaNs for calculate or not.
    Returns:    dataMean : 2D array
                outFile : string, output file name
    Examples:   avgPhaseVel = ut.temporal_average('ifgramStack.h5', datasetName='unwrapPhase')[0]
                ut.temporal_average('ifgramStack.h5', datasetName='coherence',
                                    outFile='avgSpatialCoherence.h5', updateMode=True)
    """
    atr = readfile.read_attribute(File, datasetName=datasetName)
    k = atr['FILE_TYPE']
    if k not in ['ifgramStack', 'timeseries']:
        print('WARNING: input file is not multi-temporal file: {}, return itself.'.format(File))
        data = readfile.read(File)[0]
        return data, File

    # Default output filename
    if not outFile:
        ext = os.path.splitext(File)[1]
        if not outFile:
            if k == 'ifgramStack':
                if datasetName == 'coherence':
                    outFile = 'avgSpatialCoherence.h5'
                elif 'unwrapPhase' in datasetName:
                    outFile = 'avgPhaseVelocity.h5'
                else:
                    outFile = 'avg{}.h5'.format(datasetName)
            elif k == 'timeseries':
                if k in File:
                    processMark = os.path.basename(File).split('timeseries')[1].split(ext)[0]
                    outFile = 'avgDisplacement{}.h5'.format(processMark)
            else:
                outFile = 'avg{}.h5'.format(File)

    if updateMode and os.path.isfile(outFile):
        dataMean = readfile.read(outFile)[0]
        return dataMean, outFile

    # Calculate temporal average
    if k == 'ifgramStack':
        dataMean = ifgramStack(File).temporal_average(datasetName=datasetName)
        if 'unwrapPhase' in datasetName:
            atr['FILE_TYPE'] = 'velocity'
            atr['UNIT'] = 'm/year'
        else:
            atr['FILE_TYPE'] = datasetName
    elif k == 'timeseries':
        dataMean = timeseries(File).temporal_average()
        atr['FILE_TYPE'] = 'displacement'

    if outFile:
        writefile.write(dataMean, out_file=outFile, metadata=atr)
    return dataMean, outFile


######################################################################################################
def get_file_list(file_list, abspath=False, coord=None):
    """Get all existed files matching the input list of file pattern
    Parameters: file_list - string or list of string, input file/directory pattern
                abspath - bool, return absolute path or not
                coord - string, return files with specific coordinate type: geo or radar
                    if none, skip the checking and return all files
    Returns:    file_list_out - list of string, existed file path/name, [] if not existed
    Example:    file_list = get_file_list(['*velocity*.h5','timeseries*.h5'])
                file_list = get_file_list('timeseries*.h5')
    """
    if not file_list:
        return []

    if isinstance(file_list, str):
        file_list = [file_list]

    # Get rid of None element
    file_list = [x for x in file_list if x != None]
    file_list_out = []
    for i in range(len(file_list)):
        file0 = file_list[i]
        file_list0 = glob.glob(file0)
        file_list_out += sorted(list(set(file_list0) - set(file_list_out)))

    if abspath:
        file_list_out = [os.path.abspath(i) for i in file_list_out]

    if coord is not None:
        for fname in list(file_list_out):
            atr = readfile.read_attribute(fname)
            if coord in ['geo'] and 'Y_FIRST' not in atr.keys():
                file_list_out.remove(fname)
            elif coord in ['radar', 'rdr', 'rdc'] and 'Y_FIRST' in atr.keys():
                file_list_out.remove(fname)
            else:
                raise ValueError('Input coord type: '+str(coord) +
                                 '\n. Only support geo, radar, rdr, rdc inputs.')

    return file_list_out


##################################################################
def check_file_size(fname_list, mode_width=None, mode_length=None):
    """Check file size in the list of files, and drop those not in the same size with majority."""
    # If input file list is empty
    if not fname_list:
        return fname_list, None, None

    # Read Width/Length list
    width_list = []
    length_list = []
    for fname in fname_list:
        atr = readfile.read_attribute(fname)
        width_list.append(atr['WIDTH'])
        length_list.append(atr['LENGTH'])

    # Mode of Width and Length
    if not mode_width:
        mode_width = most_common(width_list)
    if not mode_length:
        mode_length = most_common(length_list)

    # Update Input List
    fname_list_out = list(fname_list)
    if (width_list.count(mode_width) != len(width_list)
            or length_list.count(mode_length) != len(length_list)):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('WARNING: Some files may have the wrong dimensions!')
        print('All files should have the same size.')
        print('The width and length of the majority of files are: %s, %s' %
              (mode_width, mode_length))
        print('But the following files have different dimensions and thus will not be loaded:')
        for i in range(len(fname_list)):
            if width_list[i] != mode_width or length_list[i] != mode_length:
                print('%s    width: %s  length: %s' %
                      (fname_list[i], width_list[i], length_list[i]))
                fname_list_out.remove(fname_list[i])
        print('\nNumber of files left: '+str(len(fname_list_out)))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    return fname_list_out, mode_width, mode_length


def most_common(L, k=1):
    """Return the k most common item in the list L.
    Examples:
        5, 8 = most_common([4,5,5,5,5,8,8,8,9], k=2)
        'duck' = most_common(['goose','duck','duck','dog'])
        'goose' = most_common(['goose','duck','duck','goose'])
    """
    from collections import Counter
    cnt = Counter(L)
    item_mm = [i[0] for i in cnt.most_common(k)]
    if k == 1:
        item_mm = item_mm[0]
    return item_mm


######################################################################################################
def range_ground_resolution(atr, print_msg=False):
    """Get range resolution on the ground in meters,
        from ROI_PAC attributes, for file in radar coord
    """
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no range resolution info.')
        return
    inc_angle = incidence_angle(atr, dimension=0, print_msg=print_msg)
    rg_step = float(atr['RANGE_PIXEL_SIZE'])/np.sin(inc_angle/180.0*np.pi)
    return rg_step


def azimuth_ground_resolution(atr):
    """Get azimuth resolution on the ground in meters,
        from ROI_PAC attributes, for file in radar coord
    """
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no azimuth resolution info.')
        return
    try:
        proc = atr['PROCESSOR']
    except:
        proc = 'isce'
    if proc in ['roipac', 'isce']:
        Re = float(atr['EARTH_RADIUS'])
        Height = float(atr['HEIGHT'])
        az_step = float(atr['AZIMUTH_PIXEL_SIZE']) * Re/(Re+Height)
    elif proc == 'gamma':
        try:
            atr = readfile.attribute_gamma2roipac(atr)
        except:
            pass
        az_step = float(atr['AZIMUTH_PIXEL_SIZE'])
    return az_step


###################################################
def timeseries_inversion_FGLS(h5flat, h5timeseries):
    """Implementation of the SBAS algorithm.

    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms 
      h5timeseries: hdf5 file with the output from the inversion
    ##################################################
    """

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    B1 = np.array(np.linalg.pinv(B), dtype=np.float32)
    ifgram_list = list(h5flat['interferograms'].keys())
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset = h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0])
    data = dset[0:dset.shape[0], 0:dset.shape[1]]
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print('Reading in the interferograms')
    # print ifgram_num,pixel_num
    print('number of interferograms: '+str(ifgram_num))
    print('number of pixels: '+str(pixel_num))
    pixel_num_step = int(pixel_num/10)

    data = np.zeros((ifgram_num, pixel_num), np.float32)
    for ni in range(ifgram_num):
        dset = h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0], 0:dset.shape[1]]
        # print np.shape(d)

    del d
    dataPoint = np.zeros((ifgram_num, 1), np.float32)
    modelDimension = np.shape(B)[1]
    ts_data = np.zeros((date_num, pixel_num), np.float32)
    for ni in range(pixel_num):
        dataPoint = data[:, ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):
            B1tmp = np.dot(B1, np.diag(fin_ndx))
            tmpe_ratea = np.dot(B1tmp, dataPoint)
            zero = np.array([0.], np.float32)
            defo = np.concatenate((zero, np.cumsum([tmpe_ratea*dt])))
            ts_data[:, ni] = defo
        if not np.remainder(ni, pixel_num_step):
            print('Processing point: %8d of %8d, %3d' %
                  (ni, pixel_num, (10*ni/pixel_num_step))+'%')
    del data
    timeseries = np.zeros((date_num, np.shape(dset)[0], np.shape(dset)[1]), np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1], np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.items():
        timeseriesDict[key] = value

    dateIndex = {}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key, value in iter(timeseriesDict.items()):
            group.attrs[key] = value

    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]])
    print('Time series inversion took ' + str(time.time()-total) + ' secs')
    return


def timeseries_inversion_L1(h5flat, h5timeseries):
    try:
        from .l1 import l1
        from cvxopt import normal, matrix
    except ImportError:
        raise ImportError('cvxopt should be installed to be able to use the L1 norm minimization.')
        # modified from sbas.py written by scott baker, 2012

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    BL1 = matrix(B)
    B1 = np.array(np.linalg.pinv(B), dtype=np.float32)
    ifgram_list = list(h5flat['interferograms'].keys())
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset = h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0])
    data = dset[0:dset.shape[0], 0:dset.shape[1]]
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print('Reading in the interferograms')
    print(ifgram_num, pixel_num)

    #data = np.zeros((ifgram_num,pixel_num),np.float32)
    data = np.zeros((ifgram_num, pixel_num))
    for ni in range(ifgram_num):
        dset = h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0], 0:dset.shape[1]]
        # print np.shape(d)

        data[ni] = d.flatten(1)
    del d
    dataPoint = np.zeros((ifgram_num, 1), np.float32)
    modelDimension = np.shape(B)[1]
    ts_data = np.zeros((date_num, pixel_num), np.float32)
    print(data.shape)
    DataL1 = matrix(data)
    L1ORL2 = np.ones((pixel_num, 1))
    for ni in range(pixel_num):
        print(ni)
        dataPoint = data[:, ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):

            B1tmp = np.dot(B1, np.diag(fin_ndx))
            #tmpe_ratea = np.dot(B1tmp,dataPoint)
            try:
                tmpe_ratea = np.array(l1(BL1, DataL1[:, ni]))
                zero = np.array([0.], np.float32)
                defo = np.concatenate((zero, np.cumsum([tmpe_ratea[:, 0]*dt])))
            except:
                tmpe_ratea = np.dot(B1tmp, dataPoint)
                L1ORL2[ni] = 0
                zero = np.array([0.], np.float32)
                defo = np.concatenate((zero, np.cumsum([tmpe_ratea*dt])))

            ts_data[:, ni] = defo
        if not np.remainder(ni, 10000):
            print('Processing point: %7d of %7d ' % (ni, pixel_num))
    del data
    timeseries = np.zeros((date_num, np.shape(dset)[0], np.shape(dset)[1]), np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1], np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
    L1ORL2 = np.reshape(L1ORL2, (np.shape(dset)[1], np.shape(dset)[0])).T

    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.items():
        timeseriesDict[key] = value

    dateIndex = {}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key, value in iter(timeseriesDict.items()):
            group.attrs[key] = value

    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]])
    print('Time series inversion took ' + str(time.time()-total) + ' secs')
    L1orL2h5 = h5py.File('L1orL2.h5', 'w')
    gr = L1orL2h5.create_group('mask')
    dset = gr.create_dataset('mask', data=L1ORL2, compression='gzip')
    L1orL2h5.close()
    return


#########################################################################################
def run_deramp(fname, ramp_type, mask_file=None, out_file=None, datasetName=None):
    """ Remove ramp from each 2D matrix of input file
    Parameters: fname     : str, data file to be derampped
                ramp_type : str, name of ramp to be estimated.
                mask_file : str, file of mask of pixels used for ramp estimation
                out_file  : str, output file name
                datasetName : str, output dataset name, for ifgramStack file type only
    Returns:    out_file  : str, output file name
    """
    print('remove {} ramp from file: {}'.format(ramp_type, fname))
    if not out_file:
        fbase, fext = os.path.splitext(fname)
        out_file = '{}_ramp{}'.format(fbase, fext)

    start_time = time.time()
    atr = readfile.read_attribute(fname)

    # mask
    if os.path.isfile(mask_file):
        mask = readfile.read(mask_file, datasetName='mask')[0]
        print('read mask file: '+mask_file)
    else:
        mask = np.ones((int(atr['LENGTH']), int(atr['WIDTH'])))
        print('use mask of the whole area')

    # deramping
    k = atr['FILE_TYPE']
    if k == 'timeseries':
        print('reading data ...')
        data = readfile.read(fname)[0]
        print('estimating phase ramp ...')
        data = deramp(data, mask, ramp_type=ramp_type, metadata=atr)[0]
        writefile.write(data, out_file, ref_file=fname)

    elif k == 'ifgramStack':
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        if not datasetName:
            datasetName = 'unwrapPhase'
        with h5py.File(fname, 'a') as f:
            ds = f[datasetName]
            dsNameOut = '{}_ramp'.format(datasetName)
            if dsNameOut in f.keys():
                dsOut = f[dsNameOut]
                print('access HDF5 dataset /{}'.format(dsNameOut))
            else:
                dsOut = f.create_dataset(dsNameOut, shape=(obj.numIfgram, obj.length, obj.width),
                                         dtype=np.float32, chunks=True, compression=None)
                print('create HDF5 dataset /{}'.format(dsNameOut))

            prog_bar = ptime.progressBar(maxValue=obj.numIfgram)
            for i in range(obj.numIfgram):
                data = ds[i, :, :]
                data = deramp(data, mask, ramp_type=ramp_type, metadata=atr)[0]
                dsOut[i, :, :] = data
                prog_bar.update(i+1, suffix='{}/{}'.format(i+1, obj.numIfgram))
            prog_bar.close()
            print('finished writing to file: '.format(fname))

    # Single Dataset File
    else:
        data = readfile.read(fname)[0]
        data = deramp(data, mask, ramp_type, metadata=atr)[0]
        print('writing >>> {}'.format(out_file))
        writefile.write(data, out_file=out_file, ref_file=fname)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return out_file
#########################################################################################


#####################################  coordinate class begin ##############################################
class coordinate:
    """
    Coordinates conversion based lookup file in pixel-level accuracy

    Example:
        atr = readfile.read('velocity.h5')
        coord = ut.coordinate(atr, lookup_file='INPUTS/geometryRadar.h5')
        y, x = coord.geo2radar(lat, lon)
        lat, lon = coord.radar2geo(y, x)
    """

    def __init__(self, metadata, lookup_file=None):
        """Define a coordinate object
        Parameters: metadata : dict, source metadata
                    lookup_file : list of 2 strings, or string, lookup table file(s)
        Example:    from pysar.utils import readfile, utils as ut
                    atr = readfile.read_attribute('./velocity.h5')
                    coord = ut.coordinate(atr, './INPUTS/geometryRadar.h5')
                    coord.geo2radar(33.450, -90.22)
                    coord.radar2geo(50, 200)
        """
        self.src_metadata = metadata
        if lookup_file is None:
            lookup_file = get_lookup_file(lookup_file, print_msg=False)
        if isinstance(lookup_file, str):
            lookup_file = [lookup_file, lookup_file]
        self.lookup_file = lookup_file
        self.lut_y = None
        self.lut_x = None

    def open(self):
        try:
            self.earth_radius = float(self.src_metadata['EARTH_RADIUS'])
        except:
            self.earth_radius = 6371.0e3

        if 'Y_FIRST' in self.src_metadata.keys():
            self.geocoded = True
            self.lat0 = float(self.src_metadata['Y_FIRST'])
            self.lon0 = float(self.src_metadata['X_FIRST'])
            self.lat_step = float(self.src_metadata['Y_STEP'])
            self.lon_step = float(self.src_metadata['X_STEP'])
        else:
            self.geocoded = False
            if self.lookup_file:
                self.lut_metadata = readfile.read_attribute(self.lookup_file[0])

    def lalo2yx(self, coord_in, coord_type):
        """convert geo coordinates into radar coordinates (round to nearest integer)
            for Geocoded file only
        Parameters: geoCoord  : coordinate (list / tuple) in latitude/longitude in float
                    metadata : dictionary of file attributes
                    coord_type : coordinate type: latitude, longitude
        Example:    300 = coordinate.lalo2yx(32.104990,    metadata,'lat')
                    [1000,1500] = coordinate.lalo2yx([130.5,131.4],metadata,'lon')
        """
        self.open()
        if not self.geocoded:
            raise ValueError('Input file is not geocoded.')

        # input format
        if isinstance(coord_in, np.ndarray):
            coord_in = coord_in.tolist()
        if isinstance(coord_in, float):
            coord_in = [coord_in]
        coord_in = list(coord_in)

        # convert coordinates
        coord_type = coord_type.lower()
        coord_out = []
        for i in range(len(coord_in)):
            if coord_type.startswith('lat'):
                coord = int(np.rint((coord_in[i] - self.lat0) / self.lat_step))
            elif coord_type.startswith('lon'):
                coord = int(np.rint((coord_in[i] - self.lon0) / self.lon_step))
            else:
                raise ValueError('Unrecognized coordinate type: '+coord_type)
            coord_out.append(coord)

        # output format
        if len(coord_out) == 1:
            coord_out = coord_out[0]
        elif isinstance(coord_in, tuple):
            coord_out = tuple(coord_out)
        return coord_out


    def yx2lalo(self, coord_in, coord_type):
        """convert radar coordinates into geo coordinates (pixel UL corner)
            for Geocoded file only
        Parameters: coord_in : coordinate (list) in row/col in int
                    metadata : dictionary of file attributes
                    coord_type  : coordinate type: row, col, y, x
        Example:    32.104990 = coord_yx2lalo(300, metadata, 'y')
                    [130.5,131.4] = coord_yx2lalo([1000,1500], metadata, 'x')
        """
        self.open()
        if not self.geocoded:
            raise ValueError('Input file is not geocoded.')

        # Convert to List if input is String
        if isinstance(coord_in, np.ndarray):
            coord_in = coord_in.tolist()
        if isinstance(coord_in, int):
            coord_in = [coord_in]
        coord_in = list(coord_in)

        coord_type = coord_type.lower()
        coord_out = []
        for i in range(len(coord_in)):
            if coord_type.startswith(('row', 'y', 'az', 'azimuth')):
                coord = coord_in[i] * self.lat_step + self.lat0
            elif coord_type.startswith(('col', 'x', 'rg', 'range')):
                coord = coord_in[i] * self.lon_step + self.lon0
            else:
                raise ValueError('Unrecognized coordinate type: '+coord_type)
            coord_out.append(coord)

        if len(coord_out) == 1:
            coord_out = coord_out[0]
        elif isinstance(coord_in, tuple):
            coord_out = tuple(coord_out)
        return coord_out


    def _get_lookup_row_col(self, y, x, y_factor=10, x_factor=10, geo_coord=False, debug_mode=False):
        """Get row/col number in y/x value matrix from input y/x
        Use overlap mean value between y and x buffer;
        To support point outside of value pool/matrix, could use np.polyfit to fit a line
        for y and x value buffer and return the intersection point row/col
        """
        ymin = y - y_factor;  ymax = y + y_factor
        xmin = x - x_factor;  xmax = x + x_factor
        if not geo_coord:
            ymin = max(ymin, 0.5)
            xmin = max(xmin, 0.5)
        mask_y = np.multiply(self.lut_y >= ymin, self.lut_y <= ymax)
        mask_x = np.multiply(self.lut_x >= xmin, self.lut_x <= xmax)
        mask_yx = np.multiply(mask_y, mask_x)
        row, col = np.nanmean(np.where(mask_yx), axis=1)

        # for debugging only
        if debug_mode:
            print('Debug mode is ON.\nShow the row/col number searching result.')
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(nrows=1, ncols=3, figsize=[12, 5])
            axs[0].imshow(mask_y);  axs[0].set_title('Buffer in Y direction')
            axs[1].imshow(mask_x);  axs[1].set_title('Buffer in X direction')
            axs[2].imshow(mask_yx); axs[2].set_title('Y & X overlap (zoom in)')

            idx = np.where(np.sum(mask_yx, axis=0))[0]
            idy = np.where(np.sum(mask_yx, axis=1))[0]
            axs[2].set_xlim(idx[0], idx[-1])
            axs[2].set_ylim(idy[0], idy[-1])
            axs[1].set_yticklabels([])
            plt.show()

        # Error message
        if any(np.isnan(i) for i in [row, col]):
            raise RuntimeError('No coresponding coordinate found for y/x: {}/{}'.format(y, x))

        return row, col


    def read_lookup_table(self, print_msg=True):
        if 'Y_FIRST' in self.lut_metadata.keys():
            self.lut_y = readfile.read(self.lookup_file[0],
                                       datasetName='azimuthCoord',
                                       print_msg=print_msg)[0]
            self.lut_x = readfile.read(self.lookup_file[1],
                                       datasetName='rangeCoord',
                                       print_msg=print_msg)[0]
        else:
            self.lut_y = readfile.read(self.lookup_file[0],
                                       datasetName='latitude',
                                       print_msg=print_msg)[0]
            self.lut_x = readfile.read(self.lookup_file[1],
                                       datasetName='longitude',
                                       print_msg=print_msg)[0]
        return self.lut_y, self.lut_x

    def _read_geo_lut_metadata(self):
        """Read lat/lon0, lat/lon_step_deg, lat/lon_step into a Namespace - lut"""
        lut = Namespace()
        lut.lat0 = float(self.lut_metadata['Y_FIRST'])
        lut.lon0 = float(self.lut_metadata['X_FIRST'])
        lut.lat_step_deg = float(self.lut_metadata['Y_STEP'])
        lut.lon_step_deg = float(self.lut_metadata['X_STEP'])

        # Get lat/lon resolution/step in meter
        length = int(self.lut_metadata['LENGTH'])
        lat_c = lut.lat0 + lut.lat_step_deg * length / 2.
        lut.lat_step = lut.lat_step_deg * np.pi/180.0 * self.earth_radius
        lut.lon_step = lut.lon_step_deg * np.pi/180.0 * self.earth_radius * np.cos(lat_c * np.pi/180)
        return lut

    def geo2radar(self, lat, lon, print_msg=True, debug_mode=False):
        """Convert geo coordinates into radar coordinates.
        Parameters: lat/lon : np.array, float, latitude/longitude
        Returns:    az/rg : np.array, float, range/azimuth pixel number
                    az/rg_res : float, residul/uncertainty of coordinate conversion
        """
        self.open()
        if self.geocoded:
            az = self.lalo2yx(lat, coord_type='lat')
            rg = self.lalo2yx(lon, coord_type='lon')
            return az, rg, 0, 0

        if not isinstance(lat, np.ndarray):
            lat = np.array(lat)
            lon = np.array(lon)

        # read lookup table
        if self.lookup_file is None:
            raise FileNotFoundError('No lookup table file found!')
        if self.lut_y is None or self.lut_x is None:
            self.read_lookup_table(print_msg=print_msg)

        # For lookup table in geo-coord, read value directly (GAMMA and ROI_PAC)
        if 'Y_FIRST' in self.lut_metadata.keys():
            lut = self._read_geo_lut_metadata()

            # if source data file is subsetted before
            az0 = 0
            rg0 = 0
            if 'SUBSET_YMIN' in self.src_metadata.keys():
                az0 = int(self.src_metadata['SUBSET_YMIN'])
            if 'SUBSET_XMIN' in self.src_metadata.keys():
                rg0 = int(self.src_metadata['SUBSET_XMIN'])

            # uncertainty due to different resolution between src and lut file
            try:
                az_step = azimuth_ground_resolution(self.src_metadata)
                rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
                x_factor = np.ceil(abs(lut.lon_step) / rg_step).astype(int)
                y_factor = np.ceil(abs(lut.lat_step) / az_step).astype(int)
            except:
                x_factor = 2
                y_factor = 2

            # read y/x value from lookup table
            row = np.rint((lat - lut.lat0) / lut.lat_step_deg).astype(int)
            col = np.rint((lon - lut.lon0) / lut.lon_step_deg).astype(int)
            rg = np.rint(self.lut_x[row, col]).astype(int) - rg0
            az = np.rint(self.lut_y[row, col]).astype(int) - az0

        # For lookup table in radar-coord, search the buffer and use center pixel (ISCE)
        else:
            # get resolution in degree in range/azimuth direction
            az_step = azimuth_ground_resolution(self.src_metadata)
            rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
            lat_c = (np.nanmax(lat) + np.nanmin(lat)) / 2.
            az_step_deg = 180./np.pi * az_step / (self.earth_radius)
            rg_step_deg = 180./np.pi * rg_step / (self.earth_radius * np.cos(lat_c * np.pi/180.))

            az, rg = np.zeros(lat.shape), np.zeros(lat.shape)
            x_factor = 10
            y_factor = 10

            # search the overlap area of buffer in x/y direction and use the cross center
            if lat.size == 1:
                az, rg = self._get_lookup_row_col(lat, lon,
                                                  y_factor*az_step_deg,
                                                  x_factor*rg_step_deg,
                                                  geo_coord=True,
                                                  debug_mode=debug_mode)
            else:
                for i in range(rg.size):
                    az[i], rg[i] = self._get_lookup_row_col(lat[i], lon[i],
                                                            y_factor*az_step_deg,
                                                            x_factor*rg_step_deg,
                                                            geo_coord=True,
                                                            debug_mode=debug_mode)
            az = np.rint(az).astype(int)
            rg = np.rint(rg).astype(int)

        rg_resid = x_factor
        az_resid = y_factor
        return az, rg, az_resid, rg_resid


    def radar2geo(self, az, rg, print_msg=True, debug_mode=False):
        """Convert radar coordinates into geo coordinates
        Parameters: rg/az : np.array, int, range/azimuth pixel number
        Returns:    lon/lat : np.array, float, longitude/latitude of input point (rg,az);
                        nan if not found.
                    latlon_res : float, residul/uncertainty of coordinate conversion
        """
        self.open()
        if self.geocoded:
            lat = self.yx2lalo(az, coord_type='az')
            lon = self.yx2lalo(rg, coord_type='rg')
            return lat, lon, 0, 0

        if not isinstance(az, np.ndarray):
            az = np.array(az)
            rg = np.array(rg)

        # read lookup table file
        if self.lookup_file is None:
            raise FileNotFoundError('No lookup table file found!')
        if self.lut_y is None or self.lut_x is None:
            self.read_lookup_table(print_msg=print_msg)

        # For lookup table in geo-coord, search the buffer and use center pixel
        if 'Y_FIRST' in self.lut_metadata.keys():
            lut = self._read_geo_lut_metadata()

            # Get buffer ratio from range/azimuth ground resolution/step
            try:
                az_step = azimuth_ground_resolution(self.src_metadata)
                rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
                x_factor = 2 * np.ceil(abs(lut.lon_step) / rg_step)
                y_factor = 2 * np.ceil(abs(lut.lat_step) / az_step)
            except:
                x_factor = 10
                y_factor = 10

            if 'SUBSET_XMIN' in self.src_metadata.keys():
                rg += int(self.src_metadata['SUBSET_XMIN'])
                az += int(self.src_metadata['SUBSET_YMIN'])

            lut_row = np.zeros(rg.shape)
            lut_col = np.zeros(rg.shape)
            if rg.size == 1:
                lut_row, lut_col = self._get_lookup_row_col(az, rg, y_factor, x_factor,
                                                            debug_mode=debug_mode)
            else:
                for i in range(rg.size):
                    (lut_row[i],
                     lut_col[i]) = self._get_lookup_row_col(az[i], rg[i],
                                                            y_factor, x_factor,
                                                            debug_mode=debug_mode)
            lat = lut_row * lut.lat_step_deg + lut.lat0
            lon = lut_col * lut.lon_step_deg + lut.lon0
            lat_resid = abs(y_factor * lut.lat_step_deg)
            lon_resid = abs(x_factor * lut.lon_step_deg)

        # For lookup table in radar-coord, read the value directly.
        else:
            lat = self.lut_y[az, rg]
            lon = self.lut_x[az, rg]
            lat_resid = 0.
            lon_resid = 0.
        return lat, lon, lat_resid, lon_resid

    def box_pixel2geo(self, pixel_box):
        """Convert pixel_box to geo_box
        Parameters: pixel_box : list/tuple of 4 int   in (x0, y0, x1, y1)
        Returns:    geo_box   : tuple      of 4 float in (W, N, E, S)
        """
        try:
            lat = self.yx2lalo([pixel_box[1], pixel_box[3]], coord_type='y')
            lon = self.yx2lalo([pixel_box[0], pixel_box[2]], coord_type='x')
            geo_box = (lon[0], lat[0], lon[1], lat[1])
        except:
            geo_box = None
        return geo_box

    def box_geo2pixel(self, geo_box):
        """Convert geo_box to pixel_box
        Parameters: geo_box   : tuple      of 4 float in (W, N, E, S)
        Returns:    pixel_box : list/tuple of 4 int   in (x0, y0, x1, y1)
        """
        try:
            y = self.lalo2yx([geo_box[1], geo_box[3]], coord_type='latitude')
            x = self.lalo2yx([geo_box[0], geo_box[2]], coord_type='longitude')
            pixel_box = (x[0], y[0], x[1], y[1])
        except:
            pixel_box = None
        return pixel_box

    def bbox_radar2geo(self, pix_box, print_msg=False):
        """Calculate bounding box in lat/lon for file in geo coord, based on input radar/pixel box
        Parameters: pix_box - tuple of 4 int, indicating the UL/LR x/y
        Returns:    geo_box - tuple of 4 float, indicating the UL/LR lon/lat of the bounding box
        """
        x = np.array([pix_box[0], pix_box[2], pix_box[0], pix_box[2]])
        y = np.array([pix_box[1], pix_box[1], pix_box[3], pix_box[3]])
        lat, lon, lat_res, lon_res = self.radar2geo(y, x, print_msg=print_msg)
        buf = 2*(np.max(np.abs([lat_res, lon_res])))
        geo_box = (np.min(lon) - buf, np.max(lat) + buf,
                   np.max(lon) + buf, np.min(lat) - buf)
        return geo_box

    def bbox_geo2radar(self, geo_box, print_msg=False):
        """Calculate bounding box in x/y for file in radar coord, based on input geo box.
        Parameters: geo_box - tuple of 4 float, indicating the UL/LR lon/lat 
        Returns:    pix_box - tuple of 4 int, indicating the UL/LR x/y of the bounding box in radar coord
                          for the corresponding lat/lon coverage.
        """
        lat = np.array([geo_box[3], geo_box[3], geo_box[1], geo_box[1]])
        lon = np.array([geo_box[0], geo_box[2], geo_box[0], geo_box[2]])
        y, x, y_res, x_res = self.geo2radar(lat, lon, print_msg=print_msg)
        buf = 2 * np.max(np.abs([x_res, y_res]))
        pix_box = (np.min(x) - buf, np.min(y) - buf,
                   np.max(x) + buf, np.max(y) + buf)
        return pix_box

    def check_box_within_data_coverage(self, pixel_box):
        """Check the subset box's conflict with data coverage
        Parameters:  pixel_box : 4-tuple of int, indicating y/x coordinates of subset
        Returns:     out_box   : 4-tuple of int
        """
        self.open()
        length, width = int(self.src_metadata['LENGTH']), int(self.src_metadata['WIDTH'])
        sub_x = [pixel_box[0], pixel_box[2]]
        sub_y = [pixel_box[1], pixel_box[3]]

        if sub_y[0] >= length or sub_y[1] <= 0 or sub_x[0] >= width or sub_x[1] <= 0:
            data_box = (0, 0, width, length)
            msg = 'ERROR: input index is out of data range!\n'
            msg += '\tdata   range in x/y: {}\n'.format(data_box)
            msg += '\tsubset range in x/y: {}\n'.format(pixel_box)
            msg += '\tdata   range in lat/lon: {}\n'.format(self.box_pixel2geo(data_box))
            msg += '\tsubset range in lat/lon: {}\n'.format(self.box_pixel2geo(pixel_box))
            raise ValueError(msg)

        # Check Y/Azimuth/Latitude subset range
        if sub_y[0] < 0:
            sub_y[0] = 0
            print('WARNING: input y < min (0)! Set it to min.')
        if sub_y[1] > length:
            sub_y[1] = length
            print('WARNING: input y > max ({})! Set it to max.'.format(length))

        # Check X/Range/Longitude subset range
        if sub_x[0] < 0:
            sub_x[0] = 0
            print('WARNING: input x < min (0)! Set it to min.')
        if sub_x[1] > width:
            sub_x[1] = width
            print('WARNING: input x > max ({})! Set it to max.'.format(width))

        out_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])
        return out_box

#####################################  coordinate class end ##############################################
