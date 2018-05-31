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
import warnings
import h5py
import numpy as np
import multiprocessing
from pysar.utils import readfile, writefile, ptime, network as pnet, deramp
from pysar.objects import timeseries, geometry, ifgramStack, geometryDatasetNames, ifgramDatasetNames


###############################################################################
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


def coord_geo2radar(coord_in, metadata, coord_name):
    """convert geo coordinates into radar coordinates (round to nearest integer)
        for Geocoded file only
    Parameters: geoCoord  : coordinate (list / tuple) in latitude/longitude in float
                metadata : dictionary of file attributes
                coord_name : coordinate type: latitude, longitude
    Example:    300 = coord_geo2radar(32.104990,    metadata,'lat')
                [1000,1500] = coord_geo2radar([130.5,131.4],metadata,'lon')
    """
    try:
        metadata['X_FIRST']
    except:
        raise Exception('Support geocoded file only!')

    # input format
    if isinstance(coord_in, float):
        coord_in = [coord_in]
    coord_in = list(coord_in)

    # convert coordinates
    coord_out = []
    for i in range(len(coord_in)):
        if coord_name.lower().startswith('lat'):
            coord = np.rint((coord_in[i]-float(metadata['Y_FIRST'])) / float(metadata['Y_STEP']))
        elif coord_name.lower().startswith('lon'):
            coord = np.rint((coord_in[i]-float(metadata['X_FIRST'])) / float(metadata['X_STEP']))
        else:
            print('Unrecognized coordinate type: '+coord_name)
        coord_out.append(int(coord))

    # output format
    if len(coord_out) == 1:
        coord_out = coord_out[0]
    elif isinstance(coord_in, tuple):
        coord_out = tuple(coord_out)

    return coord_out


################################################################
def coord_radar2geo(coord_in, metadata, coord_name):
    """convert radar coordinates into geo coordinates (pixel UL corner)
        for Geocoded file only
    Parameters: coord_in : coordinate (list) in row/col in int
                metadata : dictionary of file attributes
                coord_name  : coordinate type: row, col, y, x
    Example:    32.104990 = coord_radar2geo(300, metadata, 'y')
                [130.5,131.4] = coord_radar2geo([1000,1500], metadata, 'x')
    """
    try:
        metadata['X_FIRST']
    except:
        raise Exception('Support geocoded file only!')

    # Convert to List if input is String
    if isinstance(coord_in, int):
        coord_in = [coord_in]
    coord_in = list(coord_in)

    coord_out = []
    for i in range(len(coord_in)):
        if coord_name.lower().startswith(('row', 'y')):
            coord = coord_in[i]*float(metadata['Y_STEP']) + float(metadata['Y_FIRST'])
        elif coord_name.lower().startswith(('col', 'x')):
            coord = coord_in[i]*float(metadata['X_STEP']) + float(metadata['X_FIRST'])
        else:
            print('Unrecognized coordinate type: '+coord_name)
        coord_out.append(coord)
    # coord_out.sort()

    if len(coord_out) == 1:
        coord_out = coord_out[0]
    elif isinstance(coord_in, tuple):
        coord_out = tuple(coord_out)

    return coord_out


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

    atr = dict()
    for key, value in iter(atr_dict.items()):
        atr[key] = str(value)

    # Update attribute variable
    atr['LENGTH'] = str(sub_y[1]-sub_y[0])
    atr['WIDTH'] = str(sub_x[1]-sub_x[0])
    atr['YMAX'] = str(sub_y[1]-sub_y[0] - 1)
    atr['XMAX'] = str(sub_x[1]-sub_x[0] - 1)
    if print_msg:
        print('update LENGTH, WIDTH, Y/XMAX')

    # Subset atribute
    if print_msg:
        print('update/add SUBSET_YMIN/YMAX/XMIN/XMAX')
    try:
        subset_y0_ori = int(atr['SUBSET_YMIN'])
        atr['SUBSET_YMIN'] = str(sub_y[0] + subset_y0_ori)
        atr['SUBSET_YMAX'] = str(sub_y[1] + subset_y0_ori)
    except:
        atr['SUBSET_YMIN'] = str(sub_y[0])
        atr['SUBSET_YMAX'] = str(sub_y[1])
    try:
        subset_x0_ori = int(atr['SUBSET_XMIN'])
        atr['SUBSET_XMIN'] = str(sub_x[0] + subset_x0_ori)
        atr['SUBSET_XMAX'] = str(sub_x[1] + subset_x0_ori)
    except:
        atr['SUBSET_XMIN'] = str(sub_x[0])
        atr['SUBSET_XMAX'] = str(sub_x[1])

    # Geo coord
    try:
        atr['Y_FIRST'] = str(float(atr['Y_FIRST']) + sub_y[0]*float(atr['Y_STEP']))
        atr['X_FIRST'] = str(float(atr['X_FIRST']) + sub_x[0]*float(atr['X_STEP']))
        if print_msg:
            print('update Y/X_FIRST')
    except:
        pass

    # Reference in space
    try:
        atr['REF_Y'] = str(int(atr['REF_Y']) - sub_y[0])
        atr['REF_X'] = str(int(atr['REF_X']) - sub_x[0])
        if print_msg:
            print('update REF_Y/X')
    except:
        pass

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


def round_to_1(x):
    """Return the most significant digit of input number"""
    return round(x, -int(np.floor(np.log10(abs(x)))))


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
        if 'Y_FIRST' in atr.keys():
            dsName2check = 'rangeCoord'
        else:
            dsName2check = 'latitude'
        try:
            dset = readfile.read(fname, datasetName=dsName2check, print_msg=False)[0]
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
        sys.exit('Unrecognized geometry dataset name: %s' % (dset))

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

    Parameters: work_dir : string,
                    PySAR working directory
                inps : Namespace, optional
                    variable for pysarApp.py. Not needed for check loading result.
    Returns:    loadComplete : bool,
                    complete loading or not
            or  inps : Namespace, if it's inputed
                atr : dict,
                    metadata of found ifgramStack file
    Example:
        #if True, PROCESS, SLC folder could be removed.
        True = ut.check_loaded_dataset($SCRATCHDIR+'/SinabungT495F50AlosA/PYSAR')
        inps,atr = ut.check_loaded_dataset(inps.workDir, inps)
    """
    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    if inps:
        inps.stackFile = None
        inps.geomFile = None
        inps.lookupFile = None

    # Required files - interferograms stack
    file_list = [os.path.join(work_dir, 'INPUTS/ifgramStack.h5')]
    stack_file = is_file_exist(file_list, abspath=True)

    if not stack_file:
        if inps:
            return inps, None
        else:
            return False

    atr = readfile.read_attribute(stack_file)
    # Check required dataset - unwrapPhase
    dsList = []
    with h5py.File(stack_file, 'r') as f:
        dsList = list(f.keys())
    if ifgramDatasetNames[0] not in dsList:
        stack_file = None

    # Recommended files - geometry (None if not found)
    if 'X_FIRST' in atr.keys():
        geocoded = True
        file_list = [os.path.join(work_dir, 'INPUTS/geometryGeo.h5')]
    else:
        geocoded = False
        file_list = [os.path.join(work_dir, 'INPUTS/geometryRadar.h5')]
    geom_file = is_file_exist(file_list, abspath=True)
    # Check required dataset - height
    if geom_file is not None:
        geom_obj = geometry(geom_file)
        geom_obj.open(print_msg=False)
        if geometryDatasetNames[0] not in geom_obj.datasetNames:
            geom_file = None

    # Recommended files - lookup table (None if not found)
    # could be different than geometry file in case of roipac and gamma
    file_list = [os.path.join(work_dir, 'INPUTS/geometry*.h5')]
    lookup_file = get_lookup_file(file_list, abspath=True, print_msg=print_msg)
    # Check required dataset
    if lookup_file is not None:
        lut_obj = geometry(lookup_file)
        lut_obj.open(print_msg=False)
        if not (all(i in lut_obj.datasetNames for i in ['latitude', 'longitude'])
                or all(i in lut_obj.datasetNames for i in ['rangeCoord', 'azimuthCoord'])):
            lookup_file = None

    # Set loadComplete to True only if all required datasets exists
    if any(i is None for i in [stack_file, geom_file, lookup_file]):
        loadComplete = False
    else:
        loadComplete = True

    # print message
    if print_msg:
        print('Loaded dataset are processed by InSAR software: {}'.format(atr['PROCESSOR']))
        if geocoded:
            print('Loaded dataset is in GEO coordinates')
        else:
            print('Loaded dataset is in RADAR coordinates')
        print('Interferograms Stack: {}'.format(stack_file))
        print('Geometry File       : {}'.format(geom_file))
        print('Lookup Table File   : {}'.format(lookup_file))
        if loadComplete:
            print('-'*50+'\nAll data needed found/loaded/copied. Processed 2-pass InSAR data can be removed.')
        print('-'*50)

    # Update namespace inps if inputed
    if inps:
        inps.stackFile = stack_file
        inps.geomFile = geom_file
        inps.lookupFile = lookup_file
        inps.geocoded = geocoded
        return inps, atr

    # Check
    else:
        return loadComplete


def is_file_exist(file_list, abspath=True):
    """Check if any file in the file list 1) exists and 2) readable
    Inputs:
        file_list : list of string, file name with/without wildcards
        abspath   : bool, return absolute file name/path or not
    Output:
        file_path : string, found file name/path; None if not.
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


def circle_index(atr, circle_par):
    """Return Index of Elements within a Circle centered at input pixel
    Inputs: atr : dictionary
                containging the following attributes:
                WIDT
                LENGTH
            circle_par : string in the format of 'y,x,radius'
                i.e. '200,300,20'          for radar coord
                     '31.0214,130.5699,20' for geo   coord
    Output: idx : 2D np.array in bool type
                mask matrix for those pixel falling into the circle defined by circle_par
    Examples: idx_mat = ut.circle_index(atr, '200,300,20')
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
            c_y = np.rint((c_lat-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
            c_x = np.rint((c_lon-float(atr['X_FIRST']))/float(atr['X_STEP']))
            print('Input circle index in lat/lon coord: %.4f, %.4f, %d' % (c_lat, c_lon, radius))
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
    Inputs:
        timeseries_resid_file - string, timeseries HDF5 file, e.g. timeseries_ECMWF_demErrInvResid.h5
        mask_file - string, mask file, e.g. maskTempCoh.h5
        ramp_type - string, ramp type, e.g. plane, quadratic, no for do not remove ramp
    outputs:
        std_list  - list of float, standard deviation of deramped input timeseries file
        date_list - list of string in YYYYMMDD format, corresponding dates
    Example:
        import pysar.utils.utils as ut
        std_list, date_list = ut.get_residual_std('timeseries_ECMWF_demErrInvResid.h5', 'maskTempCoh.h5')
    """
    # Intermediate files name
    if ramp_type == 'no':
        print('No ramp removal')
        deramp_file = timeseries_resid_file
    else:
        deramp_file = '{}_{}.h5'.format(os.path.splitext(timeseries_resid_file)[0], ramp_type)
    std_file = os.path.splitext(deramp_file)[0]+'_std.txt'

    # Get residual std text file
    if update_file(std_file, [deramp_file, mask_file], check_readable=False):
        if update_file(deramp_file, timeseries_resid_file):
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print('removing a {} ramp from file: '.format(ramp_type, timeseries_resid_file))
                deramp_file = deramp.remove_surface(timeseries_resid_file,
                                                    ramp_type,
                                                    mask_file,
                                                    deramp_file)
        print('calculating residual standard deviation for each epoch from file: '+deramp_file)
        std_file = timeseries(deramp_file).timeseries_std(maskFile=mask_file,
                                                          outFile=std_file)

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
                    ramp type, e.g. plane, quadratic, no for do not remove ramp
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
        deramp_file = timeseries_resid_file
    else:
        deramp_file = os.path.splitext(timeseries_resid_file)[0]+'_'+ramp_type+'.h5'
    rms_file = os.path.join(os.path.dirname(os.path.abspath(deramp_file)),
                            'rms_{}.txt'.format(os.path.splitext(deramp_file)[0]))

    # Get residual RMS text file
    if update_file(rms_file, [deramp_file, mask_file], check_readable=False):
        if update_file(deramp_file, timeseries_resid_file):
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                print('removing a {} ramp from file: {}'.format(ramp_type, timeseries_resid_file))
                deramp_file = deramp.remove_surface(timeseries_resid_file,
                                                    ramp_type,
                                                    mask_file,
                                                    deramp_file)
        print('Calculating residual RMS for each epoch from file: '+deramp_file)
        rms_file = timeseries(deramp_file).timeseries_rms(maskFile=mask_file,
                                                          outFile=rms_file)

    # Read residual RMS text file
    print('read timeseries residual RMS from file: '+rms_file)
    rms_fileContent = np.loadtxt(rms_file, dtype=bytes).astype(str)
    rms_list = rms_fileContent[:, 1].astype(np.float).tolist()
    date_list = list(rms_fileContent[:, 0])
    return rms_list, date_list, rms_file


def timeseries_coherence(inFile, maskFile='maskTempCoh.h5', outFile=None):
    """Calculate spatial average coherence for each epoch of input time series file
    Inputs:
        inFile   - string, timeseries HDF5 file
        maskFile - string, mask file 
        outFile  - string, output text file 
    Example:
        txtFile = timeseries_coherence('timeseries_ECMWF_demErrInvResid_quadratic.h5')
    """
    try:
        mask = readfile.read(maskFile, datasetName='mask')[0]
        print('read mask from file: '+maskFile)
    except:
        maskFile = None
        print('no mask input, use all pixels')

    if not outFile:
        outFile = os.path.splitext(inFile)[0]+'_coh.txt'

    atr = readfile.read_attribute(inFile)
    k = atr['FILE_TYPE']
    if not k in ['timeseries']:
        raise Exception('Only timeseries file is supported, input file is: '+k)
    range2phase = -4*np.pi/float(atr['WAVELENGTH'])

    h5 = h5py.File(inFile, 'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)

    f = open(outFile, 'w')
    f.write('# Date      spatial_average_coherence\n')
    for i in range(date_num):
        date = date_list[i]
        data = h5[k].get(date)[:]
        data = np.exp(1j*range2phase*data)
        if maskFile:
            data[mask == 0] = np.nan
        coh = np.absolute(np.nanmean(data))
        msg = '%s    %.4f' % (date, coh)
        f.write(msg+'\n')
        print(msg)
    h5.close()
    f.close()
    print('write to '+outFile)

    return outFile


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
def update_file(outFile, inFile=None, overwrite=False, check_readable=True):
    """Check whether to update outFile/outDir or not.
    return True if any of the following meets:
        1. if overwrite option set to True
        2. outFile is empty, e.g. None, []
        3. outFile is not existed
        4. outFile is not readable by readfile.read_attribute() when check_readable=True
        5. outFile is older than inFile, if inFile is not None
    Otherwise, return False.

    If inFile=None and outFile exists and readable, return False

    Parameters: inFile : string or list of string, input file(s)/directories
    Returns:    True/False : bool, whether to update output file or not
    Example:    if ut.update_file('timeseries_ECMWF_demErr.h5', 'timeseries_ECMWF.h5'):
                if ut.update_file('exclude_date.txt', check_readable=False,
                                  inFile=['timeseries_ECMWF_demErrInvResid.h5',
                                          'maskTempCoh.h5',
                                          'pysar_template.txt']):
    """
    if overwrite:
        return True

    if not outFile or (not os.path.isfile(outFile) and not os.path.isdir(outFile)):
        return True

    if check_readable:
        try:
            atr = readfile.read_attribute(outFile)
            width = atr['WIDTH']
        except:
            print(outFile+' exists, but can not read, remove it.')
            rmCmd = 'rm '+outFile
            print(rmCmd)
            os.system(rmCmd)
            return True

    if inFile:
        inFile = get_file_list(inFile)

        # Check modification time
        if inFile:
            if any(os.path.getmtime(outFile) < os.path.getmtime(File) for File in inFile):
                return True
            else:
                print('{} exists and is newer than {}, skip updating.'.format(outFile, inFile))
                return False

    return False


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


def add_attribute(File, atr_new=dict()):
    """Add/update input attribute into File
    Inputs:
        File - string, path/name of file
        atr_new - dict, attributes to be added/updated
                  if value is None, delete the item from input File attributes
    Output:
        File - string, path/name of updated file
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
            except:
                pass
        else:
            f.attrs[key] = str(value)
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
    Inputs:
        atr - dict, including the following PySAR attribute
              LENGTH
              P_BASELINE_TIMESERIES
              P_BASELINE_TOP_TIMESERIES (optional)
              P_BASELINE_BOTTOM_TIMESERIES (optional)
        dimension - int, choices = [0, 1]
                    0 for constant P_BASELINE in azimuth direction
                    1 for linear P_BASELINE in azimuth direction, for radar coord only
    Output:
        pbase - np.array, with shape = [date_num, 1] or [date_num, length]
    """
    if dimension > 0 and 'Y_FIRST' in atr.keys():
        dimension = 0
        print('file is in geo coordinates, return constant P_BASELINE for one interferogram')
    if dimension > 0 and any(i not in atr.keys() for i in ['P_BASELINE_TOP_TIMESERIES',
                                                           'P_BASELINE_BOTTOM_TIMESERIES']):
        dimension = 0
        print('No P_BASELINE_TOP/BOTTOM_TIMESERIES attributes found, return constant P_BASELINE for one interferogram')

    pbase_center = np.array(
        [float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]).reshape(-1, 1)
    if dimension == 0:
        pbase = pbase_center
    elif dimension == 1:
        pbase_top = np.array([float(i) for i in atr['P_BASELINE_TOP_TIMESERIES'].split()]).reshape(-1, 1)
        pbase_bottom = np.array([float(i) for i in atr['P_BASELINE_BOTTOM_TIMESERIES'].split()]).reshape(-1, 1)
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
                    2 for 2d matrix, vary in range direction, constant in az direction, for radar coord only
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
        inc_angle = np.array((inc_angle_n+inc_angle_f)/2.0, np.float32)
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
def spatial_average(File, datasetName=ifgramDatasetNames[1], maskFile=None, box=None, saveList=False, checkAoi=True):
    """Read/Calculate Spatial Average of input file.

    If input file is text file, read it directly;
    If input file is data matrix file:
        If corresponding text file exists with the same mask file/AOI info, read it directly;
        Otherwise, calculate it from data file.

        Only non-nan pixel is considered.
    Input:
        File     : string, path of input file
        maskFile : string, path of mask file, e.g. maskTempCoh.h5
        box      : 4-tuple defining the left, upper, right, and lower pixel coordinate
        saveList : bool, save (list of) mean value into text file
    Output:
        meanList : list for float, average value in space for each epoch of input file
        dateList : list of string for date info
                    date12_list, e.g. 101120-110220, for interferograms/coherence
                    date8_list, e.g. 20101120, for timeseries
                    file name, e.g. velocity.h5, for all the other file types
    Example:
        meanList = spatial_average('coherence.h5')[0]
        meanList, date12_list = spatial_average('coherence.h5', 'maskTempCoh.h5', saveList=True)
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
                and not update_file(txtFile, [File, maskFile], check_readable=False)):
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


def temporal_average(File, datasetName=ifgramDatasetNames[1], updateMode=False, outFile=None):
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
                elif datasetName == 'unwrapPhase':
                    outFile = 'avgPhaseVelocity.h5'
                else:
                    outFile = 'avg{}.h5'.format(datasetName)
            elif k == 'timeseries':
                if k in File:
                    processMark = os.path.basename(File).split('timeseries')[1].split(ext)[0]
                    outFile = 'avgDisplacement{}.h5'.format(processMark)
            else:
                outFile = 'avg{}.h5'.format(File)

    if updateMode and not update_file(outFile, [File]):
        dataMean = readfile.read(outFile)[0]
        return dataMean, outFile

    # Calculate temporal average
    if k == 'ifgramStack':
        dataMean = ifgramStack(File).temporal_average(datasetName=datasetName)
        if datasetName == 'unwrapPhase':
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
    Inputs:
        file_list - string or list of string, input file/directory pattern
        abspath  - bool, return absolute path or not
        coord    - string, return files with specific coordinate type: geo or radar
                   if none, skip the checking and return all files
    Output:
        file_list_out - list of string, existed file path/name, [] if not existed
    Example:
        file_list = get_file_list(['*velocity*.h5','timeseries*.h5'])
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
    if width_list.count(mode_width) != len(width_list) or length_list.count(mode_length) != len(length_list):
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


def most_common(L):
    """Return the most common item in the list L. From Alex Martelli on Stack Overflow.
    If the "most common" items with the same highest count are > 1, return the earliest-occurring one.
    Link: https://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list
    Examples:
        5 = most_common([4,5,5,5,5,8,9])
        'goose' = most_common(['goose','duck','duck','goose'])
    """
    import itertools
    import operator
    # get an iterable of (item, iterable) pairs
    SL = sorted((x, i) for i, x in enumerate(L))
    groups = itertools.groupby(SL, key=operator.itemgetter(0))
    # auxiliary function to get "quality" for an item

    def _auxfun(g):
        item, iterable = g
        count = 0
        min_index = len(L)
        for _, where in iterable:
            count += 1
            min_index = min(min_index, where)
        return count, -min_index
    # pick the highest-count/earliest item
    return max(groups, key=_auxfun)[0]


def mode(thelist):
    """Find Mode (most common) item in the list"""
    if not thelist:
        return None
    if len(thelist) == 1:
        return thelist[0]

    counts = {}
    for item in thelist:
        counts[item] = counts.get(item, 0) + 1
    maxcount = 0
    maxitem = None
    for k, v in iter(counts.items()):
        if v > maxcount:
            maxitem = k
            maxcount = v

    if maxcount == 1:
        print("All values only appear once")
        return None
    elif list(counts.values()).count(maxcount) > 1:
        print("List has multiple modes")
        return None
    else:
        return maxitem


######################################################################################################
def range_ground_resolution(atr, print_msg=False):
    """Get range resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord"""
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no range resolution info.')
        return
    inc_angle = incidence_angle(atr, dimension=0, print_msg=print_msg)
    rg_step = float(atr['RANGE_PIXEL_SIZE'])/np.sin(inc_angle/180.0*np.pi)
    return rg_step


def azimuth_ground_resolution(atr):
    """Get azimuth resolution on the ground in meters, from ROI_PAC attributes, for file in radar coord"""
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no azimuth resolution info.')
        return
    if atr['PROCESSOR'] in ['roipac', 'isce']:
        Re = float(atr['EARTH_RADIUS'])
        Height = float(atr['HEIGHT'])
        az_step = float(atr['AZIMUTH_PIXEL_SIZE']) * Re/(Re+Height)
    elif atr['PROCESSOR'] == 'gamma':
        try:
            atr = readfile.attribute_gamma2roipac(atr)
        except:
            pass
        az_step = float(atr['AZIMUTH_PIXEL_SIZE'])
    return az_step


#########################################################################
# Use geomap*.trans file for precious (pixel-level) coord conversion
def get_lookup_row_col(y, x, lut_y, lut_x, y_factor=10, x_factor=10, geoCoord=False):
    """Get row/col number in y/x value matrix from input y/x
    Use overlap mean value between y and x buffer;
    To support point outside of value pool/matrix, could use np.polyfit to fit a line
    for y and x value buffer and return the intersection point row/col
    """
    ymin = y - y_factor
    xmin = x - x_factor
    if not geoCoord:
        ymin = max(ymin, 0.5)
        xmin = max(xmin, 0.5)
    mask_y = np.multiply(lut_y >= ymin, lut_y <= (y+y_factor))
    mask_x = np.multiply(lut_x >= xmin, lut_x <= (x+x_factor))
    row, col = np.nanmean(np.where(np.multiply(mask_y, mask_x)), axis=1)
    return row, col


def glob2radar(lat, lon, lookupFile=None, atr_rdr=dict(), print_msg=True):
    """Convert geo coordinates into radar coordinates.
    Parameters: lat/lon : np.array, float, latitude/longitude
                lookupFile : string, trans/look up file
                atr_rdr : dict, attributes of file in radar coord, optional but recommended.
    Returns:    az/rg : np.array, float, range/azimuth pixel number
                az/rg_res : float, residul/uncertainty of coordinate conversion
    """
    if lookupFile is None:
        lookupFile = get_lookup_file(lookupFile)
        if lookupFile is None:
            if print_msg:
                print('WARNING: No lookup table found! Can not convert coordinates without it.')
            return None
    if isinstance(lookupFile, str):
        lookupFile = [lookupFile, lookupFile]
    atr_lut = readfile.read_attribute(lookupFile[0])
    if print_msg:
        print('reading file: '+lookupFile[0])

    # For lookup table in geo-coord, read value directly
    if 'Y_FIRST' in atr_lut.keys():
        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3
        lut_x = readfile.read(lookupFile[1], datasetName='rangeCoord', print_msg=False)[0]
        lut_y = readfile.read(lookupFile[0], datasetName='azimuthCoord', print_msg=False)[0]
        lat0 = float(atr_lut['Y_FIRST'])
        lon0 = float(atr_lut['X_FIRST'])
        lat_center = lat0 + float(atr_lut['Y_STEP'])*float(atr_lut['LENGTH'])/2
        lat_step_deg = float(atr_lut['Y_STEP'])
        lon_step_deg = float(atr_lut['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0 * earth_radius
        lon_step = lon_step_deg*np.pi/180.0 * earth_radius*np.cos(lat_center*np.pi/180)

        # Get range/azimuth ground resolution/step in meter
        x_factor = 2
        y_factor = 2
        az0 = 0
        rg0 = 0
        if 'Y_FIRST' not in atr_rdr.keys():
            try:
                az_step = azimuth_ground_resolution(atr_rdr)
                rg_step = range_ground_resolution(atr_rdr, print_msg=False)
                x_factor = np.ceil(abs(lon_step)/rg_step).astype(int)
                y_factor = np.ceil(abs(lat_step)/az_step).astype(int)
                if 'SUBSET_YMIN' in atr_rdr.keys():
                    az0 = int(atr_rdr['SUBSET_YMIN'])
                if 'SUBSET_XMIN' in atr_rdr.keys():
                    rg0 = int(atr_rdr['SUBSET_XMIN'])
            except:
                pass

        width = int(atr_lut['WIDTH'])
        row = np.rint((lat - lat0)/lat_step_deg).astype(int)
        col = np.rint((lon - lon0)/lon_step_deg).astype(int)
        rg = np.rint(lut_x[row, col]).astype(int) - rg0
        az = np.rint(lut_y[row, col]).astype(int) - az0

    # For lookup table in radar-coord, search the buffer and use center pixel
    else:
        lut_x = readfile.read(lookupFile[1], datasetName='longitude', print_msg=print_msg)[0]
        lut_y = readfile.read(lookupFile[0], datasetName='latitude', print_msg=print_msg)[0]
        az = np.zeros(lat.shape)
        rg = np.zeros(lat.shape)
        x_factor = 10
        y_factor = 10
        try:
            earth_radius = float(atr_rdr['EARTH_RADIUS'])
        except:
            earth_radius = 6371.0e3
        az_step = azimuth_ground_resolution(atr_rdr)
        rg_step = range_ground_resolution(atr_rdr, print_msg=False)
        lat0 = np.nanmax(lat)
        lat1 = np.nanmin(lat)
        az_step_deg = 180. / np.pi * az_step / earth_radius
        rg_step_deg = 180. / np.pi * rg_step / (earth_radius*np.cos((lat0+lat1)/2*np.pi/180.))

        if lat.size == 1:
            az, rg = get_lookup_row_col(lat, lon, lut_y, lut_x,
                                        y_factor*az_step_deg,
                                        x_factor*rg_step_deg,
                                        geoCoord=True)
        else:
            for i in range(rg.size):
                az[i], rg[i] = get_lookup_row_col(lat[i], lon[i], lut_y, lut_x,
                                                  y_factor*az_step_deg,
                                                  x_factor*rg_step_deg,
                                                  geoCoord=True)
        az = np.rint(az).astype(int)
        rg = np.rint(rg).astype(int)
    rg_resid = x_factor
    az_resid = y_factor
    return az, rg, az_resid, rg_resid


def radar2glob(az, rg, lookupFile=None, atr_rdr=dict(), print_msg=True):
    """Convert radar coordinates into geo coordinates
    Parameters: rg/az : np.array, int, range/azimuth pixel number
                lookupFile : string, trans/look up file
                atr_rdr : dict, attributes of file in radar coord, optional but recommended.
    Returns:    lon/lat : np.array, float, longitude/latitude of input point (rg,az); nan if not found.
                latlon_res : float, residul/uncertainty of coordinate conversion
    """
    if lookupFile is None:
        lookupFile = get_lookup_file(lookupFile)
        if lookupFile is None:
            if print_msg:
                print('WARNING: No lookup table found! Can not convert coordinates without it.')
            return None
    if isinstance(lookupFile, str):
        lookupFile = [lookupFile, lookupFile]
    atr_lut = readfile.read_attribute(lookupFile[0])
    if print_msg:
        print('reading file: '+lookupFile[0])

    # For lookup table in geo-coord, search the buffer and use center pixel
    if 'Y_FIRST' in atr_lut.keys():
        if 'SUBSET_XMIN' in atr_rdr.keys():
            rg += int(atr_rdr['SUBSET_XMIN'])
            az += int(atr_rdr['SUBSET_YMIN'])

        # Get lat/lon resolution/step in meter
        earth_radius = 6371.0e3    # in meter
        lut_x = readfile.read(lookupFile[1], datasetName='rangeCoord', print_msg=print_msg)[0]
        lut_y = readfile.read(lookupFile[0], datasetName='azimuthCoord', print_msg=print_msg)[0]
        lat0 = float(atr_lut['Y_FIRST'])
        lon0 = float(atr_lut['X_FIRST'])
        lat_center = lat0 + float(atr_lut['Y_STEP'])*float(atr_lut['LENGTH'])/2
        lat_step_deg = float(atr_lut['Y_STEP'])
        lon_step_deg = float(atr_lut['X_STEP'])
        lat_step = lat_step_deg*np.pi/180.0 * earth_radius
        lon_step = lon_step_deg*np.pi/180.0 * earth_radius*np.cos(lat_center*np.pi/180)

        # Get range/azimuth ground resolution/step
        x_factor = 10
        y_factor = 10
        if 'Y_FIRST' not in atr_rdr.keys():
            try:
                az_step = azimuth_ground_resolution(atr_rdr)
                rg_step = range_ground_resolution(atr_rdr, print_msg=False)
                x_factor = 2*np.ceil(abs(lon_step)/rg_step)
                y_factor = 2*np.ceil(abs(lat_step)/az_step)
            except:
                pass

        lut_row = np.zeros(rg.shape)
        lut_col = np.zeros(rg.shape)
        if rg.size == 1:
            lut_row, lut_col = get_lookup_row_col(az, rg, lut_y, lut_x,
                                                  y_factor, x_factor)
        else:
            for i in range(rg.size):
                lut_row[i], lut_col[i] = get_lookup_row_col(az[i], rg[i],
                                                            lut_y, lut_x,
                                                            y_factor, x_factor)
        lat = lut_row*lat_step_deg + lat0
        lon = lut_col*lon_step_deg + lon0
        lat_resid = abs(y_factor*lat_step_deg)
        lon_resid = abs(x_factor*lon_step_deg)

    # For lookup table in radar-coord, read the value directly.
    else:
        lut_x = readfile.read(lookupFile[1], datasetName='longitude', print_msg=print_msg)[0]
        lut_y = readfile.read(lookupFile[0], datasetName='latitude', print_msg=print_msg)[0]
        lat = lut_y[az, rg]
        lon = lut_x[az, rg]

        x_factor = 2
        y_factor = 2
        try:
            earth_radius = float(atr_rdr['EARTH_RADIUS'])
        except:
            earth_radius = 6371.0e3
        az_step = azimuth_ground_resolution(atr_rdr)
        rg_step = range_ground_resolution(atr_rdr, print_msg=False)
        lat0 = np.nanmax(lat)
        lat1 = np.nanmin(lat)
        az_step_deg = 180. / np.pi * az_step / earth_radius
        rg_step_deg = 180. / np.pi * rg_step / \
            (earth_radius*np.cos((lat0+lat1)/2*np.pi/180.))
        lat_resid = abs(y_factor * az_step_deg)
        lon_resid = abs(x_factor * rg_step_deg)

    return lat, lon, lat_resid, lon_resid


#########################################################################
def check_variable_name(path):
    s = path.split("/")[0]
    if len(s) > 0 and s[0] == "$":
        p0 = os.getenv(s[1:])
        path = path.replace(path.split("/")[0], p0)
    return path


#########################################################################
def hillshade(data, scale):
    """from scott baker, ptisk library """
    azdeg = 315.0
    altdeg = 45.0
    az = azdeg * np.pi/180.0
    alt = altdeg * np.pi/180.0
    dx, dy = np.gradient(data/scale)
    slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dx, dy)
    data = np.sin(alt)*np.sin(slope) + np.cos(alt) * \
        np.cos(slope)*np.cos(-az - aspect - 0.5*np.pi)
    return data


#################################################################
def date_list(h5file):
    dateList = []
    tbase = []
    k = list(h5file.keys())
    if 'interferograms' in k:
        k[0] = 'interferograms'
    elif 'coherence' in k:
        k[0] = 'coherence'
    ifgram_list = sorted(h5file[k[0]].keys())
    for ifgram in ifgram_list:
        dates = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
        dates1 = h5file[k[0]][ifgram].attrs['DATE12'].split('-')
        if dates[0][0] == '9':
            dates[0] = '19'+dates[0]
        else:
            dates[0] = '20'+dates[0]
        if dates[1][0] == '9':
            dates[1] = '19'+dates[1]
        else:
            dates[1] = '20'+dates[1]
        if not dates[0] in dateList:
            dateList.append(dates[0])
        if not dates[1] in dateList:
            dateList.append(dates[1])

    dateList.sort()
    dateList1 = []
    for ni in range(len(dateList)):
        dateList1.append(dateList[ni][2:])

    d1 = datetime.datetime(*time.strptime(dateList[0], "%Y%m%d")[0:5])
    for ni in range(len(dateList)):
        d2 = datetime.datetime(*time.strptime(dateList[ni], "%Y%m%d")[0:5])
        diff = d2-d1
        tbase.append(diff.days)
    dateDict = {}
    for i in range(len(dateList)):
        dateDict[dateList[i]] = tbase[i]
    return tbase, dateList, dateDict, dateList1


######################################
def design_matrix(ifgramFile=None, date12_list=[], referenceDate=None):
    """Make the design matrix for the inversion based on date12_list.
    Reference:
        Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002).
        A new algorithm for surface deformation monitoring based on small
        baseline differential SAR interferograms. IEEE TGRS, 40(11), 2375-2383.

    Parameters: ifgramFile : string, 
                    name/path of interferograms stack file
                date12_list : list of string in YYMMDD-YYMMDD format
                    use all date12 from ifgramFile if input is empty
    Returns:    A : 2D np.array in size of (ifgram_num, date_num-1)
                    representing date combination for each interferogram (-1 for master, 1 for slave, 0 for others)
                    used for LS and WLS optimization
                B : 2D np.array in size of (ifgram_num, date_num-1)
                    representing temporal baseline timeseries between master and slave date for each interferogram
                    used for SBAS algorithm
    """
    # Get date12_list from Inputs
    if not date12_list:
        if ifgramFile:
            from pysar.objects import ifgramStack
            date12_list = ifgramStack(ifgramFile).get_date12_list()
        else:
            raise ValueError

    # date12_list to date6_list
    m_dates = [i.split('_')[0] for i in date12_list]
    s_dates = [i.split('_')[1] for i in date12_list]
    date6_list = sorted(list(set(m_dates + s_dates)))
    tbase = np.array(ptime.date_list2tbase(date6_list)[0])
    date_num = len(date6_list)
    ifgram_num = len(date12_list)

    if not referenceDate:
        referenceDate = date6_list[0]
    referenceDate = ptime.yymmdd(referenceDate)
    refIndex = date6_list.index(referenceDate)

    # calculate design matrix
    A = np.zeros((ifgram_num, date_num))
    B = np.zeros(A.shape)
    #t = np.zeros((ifgram_num, 2))
    for i in range(ifgram_num):
        m_idx, s_idx = [date6_list.index(j) for j in date12_list[i].split('_')]
        A[i, m_idx] = -1
        A[i, s_idx] = 1
        B[i, m_idx:s_idx] = tbase[m_idx+1:s_idx+1] - tbase[m_idx:s_idx]
        #t[i,:] = [tbase[m_idx], tbase[s_idx]]

    # Remove reference date as it can not be resolved
    #A = A[:,1:]
    A = np.hstack((A[:, 0:refIndex], A[:, (refIndex+1):]))
    B = B[:, :-1]

    return A, B


###################################################
def timeseries_inversion_FGLS(h5flat, h5timeseries):
    """Implementation of the SBAS algorithm.

    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms 
      h5timeseries: hdf5 file with the output from the inversion
    ##################################################"""

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    B1 = np.linalg.pinv(B)
    B1 = np.array(B1, np.float32)
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
        # if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,pixel_num)
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


def timeseries_inversion_L1(h5flat, h5timeseries):
    try:
        from .l1 import l1
        from cvxopt import normal, matrix
    except:
        print('-----------------------------------------------------------------------')
        print('cvxopt should be installed to be able to use the L1 norm minimization.')
        print('-----------------------------------------------------------------------')
        sys.exit(1)
        # modified from sbas.py written by scott baker, 2012

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    BL1 = matrix(B)
    B1 = np.linalg.pinv(B)
    B1 = np.array(B1, np.float32)
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


def perp_baseline_ifgram2timeseries(ifgramFile, ifgram_list=[]):
    """Calculate perpendicular baseline timeseries from input interferograms file
    Input:
        ifgramFile - string, file name/path of interferograms file
        ifgram_list - list of string, group name that is used for calculation
                      use all if it's empty
    Outputs:
        pbase        - 1D np.array, P_BASELINE_TIMESERIES
        pbase_top    - 1D np.array, P_BASELINE_TOP_TIMESERIES
        pbase_bottom - 1D np.array, P_BASELINE_BOTTOM_TIMESERIES
    """
    k = readfile.read_attribute(ifgramFile)['FILE_TYPE']
    h5file = h5py.File(ifgramFile, 'r')

    if not ifgram_list:
        ifgram_list = sorted(h5file[k].keys())

    # P_BASELINE of all interferograms
    pbase_ifgram = []
    pbase_top_ifgram = []
    pbase_bottom_ifgram = []
    for ifgram in ifgram_list:
        pbase_top = float(h5file[k][ifgram].attrs['P_BASELINE_TOP_HDR'])
        pbase_bottom = float(h5file[k][ifgram].attrs['P_BASELINE_BOTTOM_HDR'])
        pbase_ifgram.append((pbase_bottom+pbase_top)/2.0)
        pbase_top_ifgram.append(pbase_top)
        pbase_bottom_ifgram.append(pbase_bottom)
    h5file.close()

    # Temporal baseline velocity
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    m_dates = [i.split('_')[0] for i in date12_list]
    s_dates = [i.split('_')[1] for i in date12_list]
    date8_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    tbase_v = np.diff(tbase_list)

    A, B = design_matrix(ifgramFile, date12_list)
    B_inv = np.linalg.pinv(B)

    pbase_rate = np.dot(B_inv, pbase_ifgram)
    pbase_top_rate = np.dot(B_inv, pbase_top_ifgram)
    pbase_bottom_rate = np.dot(B_inv, pbase_bottom_ifgram)

    zero = np.array([0.], np.float32)
    pbase = np.concatenate((zero, np.cumsum([pbase_rate*tbase_v])))
    pbase_top = np.concatenate((zero, np.cumsum([pbase_top_rate*tbase_v])))
    pbase_bottom = np.concatenate((zero, np.cumsum([pbase_bottom_rate*tbase_v])))

    return pbase, pbase_top, pbase_bottom


def dBh_dBv_timeseries(ifgramFile):
    h5file = h5py.File(ifgramFile)
    k = list(h5file.keys())
    if 'interferograms' in k:
        k[0] = 'interferograms'
    elif 'coherence' in k:
        k[0] = 'coherence'
    igramList = list(h5file[k[0]].keys())
    dBh_igram = []
    dBv_igram = []
    for igram in igramList:
        dBh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_RATE_HDR']))
        dBv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_RATE_HDR']))

    A, B = design_matrix(ifgramFile)
    tbase, dateList, dateDict, dateList1 = date_list(h5file)
    dt = np.diff(tbase)

    Bh_rate = np.dot(np.linalg.pinv(B), dBh_igram)
    zero = np.array([0.], np.float32)
    dBh = np.concatenate((zero, np.cumsum([Bh_rate*dt])))

    Bv_rate = np.dot(np.linalg.pinv(B), dBv_igram)
    zero = np.array([0.], np.float32)
    dBv = np.concatenate((zero, np.cumsum([Bv_rate*dt])))

    h5file.close()

    return dBh, dBv


def Bh_Bv_timeseries(ifgramFile):
    h5file = h5py.File(ifgramFile)
    k = list(h5file.keys())
    if 'interferograms' in k:
        k[0] = 'interferograms'
    elif 'coherence' in k:
        k[0] = 'coherence'
    igramList = list(h5file[k[0]].keys())
    Bh_igram = []
    Bv_igram = []
    for igram in igramList:
        Bh_igram.append(float(h5file[k[0]][igram].attrs['H_BASELINE_TOP_HDR']))
        Bv_igram.append(float(h5file[k[0]][igram].attrs['V_BASELINE_TOP_HDR']))

    A, B = design_matrix(ifgramFile)
    tbase, dateList, dateDict, dateList1 = date_list(h5file)
    dt = np.diff(tbase)

    Bh_rate = np.dot(np.linalg.pinv(B), Bh_igram)
    zero = np.array([0.], np.float32)
    Bh = np.concatenate((zero, np.cumsum([Bh_rate*dt])))

    Bv_rate = np.dot(np.linalg.pinv(B), Bv_igram)
    zero = np.array([0.], np.float32)
    Bv = np.concatenate((zero, np.cumsum([Bv_rate*dt])))

    h5file.close()

    return Bh, Bv


def yymmdd2YYYYMMDD(date):
    if date[0] == '9':
        date = '19'+date
    else:
        date = '20'+date
    return date


def yyyymmdd(dates):
    datesOut = []
    for date in dates:
        if len(date) == 6:
            if date[0] == '9':
                date = '19'+date
            else:
                date = '20'+date
        datesOut.append(date)
    return datesOut


def yymmdd(dates):
    datesOut = []
    for date in dates:
        if len(date) == 8:
            date = date[2:8]
        datesOut.append(date)
    return datesOut


def make_triangle(dates12, igram1, igram2, igram3):
    dates = []
    dates.append(igram1.split('-')[0])
    dates.append(igram1.split('-')[1])
    dates.append(igram2.split('-')[1])
    datesyy = []
    for d in dates:
        datesyy.append(yymmdd2YYYYMMDD(d))

    datesyy.sort()
    Igramtriangle = []
    Igramtriangle.append(datesyy[0][2:]+'-'+datesyy[1][2:])
    Igramtriangle.append(datesyy[0][2:]+'-'+datesyy[2][2:])
    Igramtriangle.append(datesyy[1][2:]+'-'+datesyy[2][2:])

    IgramtriangleIndexes = [dates12.index(Igramtriangle[0]),
                            dates12.index(Igramtriangle[1]),
                            dates12.index(Igramtriangle[2])]
    return Igramtriangle, IgramtriangleIndexes


def get_triangles(h5file):
    k = list(h5file.keys())
    igramList = list(h5file[k[0]].keys())

    dates12 = []
    for igram in igramList:
        dates12.append(h5file[k[0]][igram].attrs['DATE12'])
    Triangles = []
    Triangles_indexes = []
    for igram1 in dates12:
        igram1_date1 = igram1.split('-')[0]
        igram1_date2 = igram1.split('-')[1]

        igram2 = []
        igram2_date2 = []
        for d in dates12:
            if igram1_date2 == d.split('-')[0]:
                igram2.append(d)
                igram2_date2.append(d.split('-')[1])

        igram3 = []
        igram3_date2 = []
        for d in dates12:
            if igram1_date1 == d.split('-')[0] and d != igram1:
                igram3.append(d)
                igram3_date2.append(d.split('-')[1])

        for date in igram2_date2:
            if date in igram3_date2:
                Igramtriangle, IgramtriangleIndexes = make_triangle(dates12, igram1,
                                                                    igram2[igram2_date2.index(date)],
                                                                    igram3[igram3_date2.index(date)])
                if not Igramtriangle in Triangles:
                    Triangles.append(Igramtriangle)
                    Triangles_indexes.append(IgramtriangleIndexes)

    numTriangles = np.shape(Triangles_indexes)[0]
    curls = np.zeros([numTriangles, 3], dtype=np.int)
    for i in range(numTriangles):
        curls[i][:] = Triangles_indexes[i]

    numIgrams = len(igramList)
    C = np.zeros([numTriangles, numIgrams])
    for ni in range(numTriangles):
        C[ni][curls[ni][0]] = 1
        C[ni][curls[ni][1]] = -1
        C[ni][curls[ni][2]] = 1

    return curls, Triangles, C


def generate_curls(curlfile, h5file, Triangles, curls):
    ifgram_list = list(h5file['interferograms'].keys())
    h5curlfile = h5py.File(curlfile, 'w')
    gg = h5curlfile.create_group('interferograms')

    curl_num = np.shape(curls)[0]
    prog_bar = ptime.progressBar(maxValue=curl_num)
    for i in range(curl_num):
        ifgram1 = ifgram_list[curls[i, 0]]
        ifgram2 = ifgram_list[curls[i, 1]]
        ifgram3 = ifgram_list[curls[i, 2]]
        d1 = h5file['interferograms'][ifgram1].get(ifgram1)[:]
        d2 = h5file['interferograms'][ifgram2].get(ifgram2)[:]
        d3 = h5file['interferograms'][ifgram3].get(ifgram3)[:]

        triangle_date = Triangles[i][0]+'_'+Triangles[i][1]+'_'+Triangles[i][2]
        group = gg.create_group(triangle_date)
        dset = group.create_dataset(triangle_date, data=d1+d3-d2)
        for key, value in h5file['interferograms'][ifgram1].attrs.items():
            group.attrs[key] = value
        prog_bar.update(i+1)

    h5curlfile.close()
    prog_bar.close()
    return curlfile
