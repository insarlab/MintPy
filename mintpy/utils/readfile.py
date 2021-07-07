############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.utils import readfile


import os
import sys
import re
import warnings
import defusedxml.ElementTree as ET

import h5py
import json
import numpy as np

from mintpy.objects import (
    datasetUnitDict,
    geometry,
    giantIfgramStack,
    giantTimeseries,
    ifgramStack,
    timeseriesDatasetNames,
    timeseries,
    HDFEOS
)


SPEED_OF_LIGHT = 299792458.0   # meter/second

standardMetadataKeys = {
    # ROI_PAC/MintPy attributes
    'ALOOKS'             : ['azimuth_looks'],
    'RLOOKS'             : ['range_looks'],
    'AZIMUTH_PIXEL_SIZE' : ['azimuthPixelSize', 'azimuth_pixel_spacing', 'az_pixel_spacing'],
    'RANGE_PIXEL_SIZE'   : ['rangePixelSize', 'range_pixel_spacing', 'rg_pixel_spacing'],
    'CENTER_LINE_UTC'    : ['center_time'],
    'DATA_TYPE'          : ['dataType', 'data_type'],
    'EARTH_RADIUS'       : ['earthRadius', 'earth_radius_below_sensor', 'earth_radius'],
    'HEADING'            : ['HEADING_DEG', 'heading'],
    'HEIGHT'             : ['altitude', 'SC_height'],
    'LENGTH'             : ['length', 'FILE_LENGTH', 'lines', 'azimuth_lines', 'nlines', 'az_samp', 'interferogram_azimuth_lines'],
    'ORBIT_DIRECTION'    : ['passDirection'],
    'PLATFORM'           : ['spacecraftName', 'sensor'],
    'POLARIZATION'       : ['polarization'],
    'PRF'                : ['prf'],
    'STARTING_RANGE'     : ['startingRange', 'near_range_slc', 'near_range'],
    'WAVELENGTH'         : ['wavelength', 'Wavelength', 'radarWavelength', 'radar_wavelength'],
    'WIDTH'              : ['width', 'Width', 'samples', 'range_samp', 'interferogram_width'],
    # from PySAR [MintPy<=1.1.1]
    'REF_DATE'           : ['ref_date'],
    'REF_LAT'            : ['ref_lat'],
    'REF_LON'            : ['ref_lon'],
    'REF_X'              : ['ref_x'],
    'REF_Y'              : ['ref_y'],
    'SUBSET_XMIN'        : ['subset_x0'],
    'SUBSET_XMAX'        : ['subset_x1'],
    'SUBSET_YMIN'        : ['subset_y0'],
    'SUBSET_YMAX'        : ['subset_y1'],
    # from Gamma geo-coordinates - degree / meter
    'X_FIRST'            : ['corner_lon', 'corner_east'],
    'Y_FIRST'            : ['corner_lat', 'corner_north'],
    'X_STEP'             : ['post_lon', 'post_east'],
    'Y_STEP'             : ['post_lat', 'post_north'],

    # HDF-EOS5 attributes
    'beam_swath'     : ['swathNumber'],
    'first_frame'    : ['firstFrameNumber'],
    'last_frame'     : ['lastFrameNumber'],
    'relative_orbit' : ['trackNumber'],

    # ISCE attributes
    'number_bands' : ['bands'],
    'scheme'       : ['interleave'],
}

GDAL2ISCE_DATATYPE = {
    1 : 'BYTE',
    2 : 'uint16',
    3 : 'SHORT',
    4 : 'uint32',
    5 : 'INT',
    6 : 'FLOAT',
    7 : 'DOUBLE',
    10: 'CFLOAT',
    11: 'complex128',
}

GDAL2NUMPY_DATATYPE = {
    1 : 'uint8',
    2 : 'uint16',
    3 : 'int16',
    4 : 'uint32',
    5 : 'int32',
    6 : 'float32',
    7 : 'float64',
    10: 'complex64',
    11: 'complex128',
}

# single file (data + attributes) supported by GDAL
GDAL_FILE_EXTS = ['.tif', '.grd']

# reference: https://subversion.renater.fr/efidir/trunk/efidir_soft/doc/Programming_C_EFIDIR/header_envi.pdf
ENVI2NUMPY_DATATYPE = {
    '1' : 'uint8',
    '2' : 'int16',
    '3' : 'int32',
    '4' : 'float32',
    '5' : 'float64',
    '6' : 'complex64',
    '9' : 'complex128',
    '12': 'uint16',
    '13': 'uint32',
    '14': 'int64',
    '15': 'uint64',
}

ENVI_BAND_INTERLEAVE = {
    'BAND' : 'BSQ',
    'LINE' : 'BIL',
    'PIXEL': 'BIP',
}

ENVI_BYTE_ORDER = {
    '0': 'little-endian',
    '1': 'big-endian',
}


###########################################################################
## Slice-based data identification for data reading:
##
## slice   : np.ndarray in 2D,  with/without '-' in their names
##           each slice is unique within a file
## dataset : np.ndarray in 2D or 3D, without '-' in their names
##
## one 2D dataset can be present as one slice
##   e.g.: temporalCoherence
##         velocity
##         mask
##         ...
## one 3D dataset can be present as multiple slices
##   with '-' to connect dataset name and time info
##   e.g.: unwrapPhase-20150115_20150127
##         unwrapPhase-20150115_20150208
##         ...
##         timeseries-20150115
##         timeseries-20150127
##         ...
## one HDF5 file can be present as a combination of multiple 2D/3D datasets
##                           or as a list of slices
##
## ----------------------- Slice Nameing Examples -------------------------
## for version-1.x files:
##     unwrapPhase-20150115_20150127
##     unwrapPhase-20150115_20150208
##     timeseries-20150115
##     timeseries-20150127
##     temporalCoherence
##
## for version-0.x files: (all in 2D dataset)
##     /interferograms/diff_filt_100814-100918_4rlks.unw/diff_filt_100814-100918_4rlks.unw
##     /coherence/diff_filt_100814-100918_4rlks.cor/diff_filt_100814-100918_4rlks.cor
##     /timeseries/20100918
##     /timeseries/20101103
##     /velocity/velocity
##
## for HDF-EOS5 files:
##     /HDFEOS/GRIDS/timeseries/observation/displacement-20150115
##     /HDFEOS/GRIDS/timeseries/observation/displacement-20150127
##     /HDFEOS/GRIDS/timeseries/quality/temporalCoherence
##
## for GIAnT v1.0 files:
##     figram-20150115_20150127
##     recons-20150115
##     cmask
##
###########################################################################



#########################################################################
def read(fname, box=None, datasetName=None, print_msg=True, xstep=1, ystep=1):
    """Read one dataset and its attributes from input file.
    Parameters: fname       : str, path of file to read
                datasetName : str or list of str, slice names
                box         : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
                x/ystep     : int, number of pixels to pick/multilook for each output pixel
    Returns:    data        : 2/3-D matrix in numpy.array format, return None if failed
                atr         : dictionary, attributes of data, return None if failed
    Examples:
        from mintpy.utils import readfile
        data, atr = readfile.read('velocity.h5')
        data, atr = readfile.read('timeseries.h5')
        data, atr = readfile.read('timeseries.h5', datasetName='timeseries-20161020')
        data, atr = readfile.read('ifgramStack.h5', datasetName='unwrapPhase')
        data, atr = readfile.read('ifgramStack.h5', datasetName='unwrapPhase-20161020_20161026')
        data, atr = readfile.read('ifgramStack.h5', datasetName='coherence', box=(100,1100, 500, 2500))
        data, atr = readfile.read('geometryRadar.h5', datasetName='height')
        data, atr = readfile.read('geometryRadar.h5', datasetName='bperp')
        data, atr = readfile.read('100120-110214.unw', box=(100,1100, 500, 2500))
    """
    # metadata
    dsname4atr = None   #used to determine UNIT
    if isinstance(datasetName, list):
        dsname4atr = datasetName[0].split('-')[0]
    elif isinstance(datasetName, str):
        dsname4atr = datasetName.split('-')[0]
    atr = read_attribute(fname, datasetName=dsname4atr)

    # box
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    if not box:
        box = (0, 0, width, length)

    # Read Data
    fext = os.path.splitext(os.path.basename(fname))[1].lower()
    if fext in ['.h5', '.he5']:
        data = read_hdf5_file(fname,
                              datasetName=datasetName,
                              box=box,
                              xstep=xstep,
                              ystep=ystep,
                              print_msg=print_msg)

    else:
        data, atr = read_binary_file(fname,
                                     datasetName=datasetName,
                                     box=box,
                                     xstep=xstep,
                                     ystep=ystep)
    return data, atr


#########################################################################
def read_hdf5_file(fname, datasetName=None, box=None, xstep=1, ystep=1, print_msg=True):
    """
    Parameters: fname       : str, name of HDF5 file to read
                datasetName : str or list of str, dataset name in root level with/without date info
                    'timeseries'
                    'timeseries-20150215'
                    'unwrapPhase'
                    'unwrapPhase-20150215_20150227'
                    'HDFEOS/GRIDS/timeseries/observation/displacement'
                    'recons'
                    'recons-20150215'
                    ['recons-20150215', 'recons-20150227', ...]
                    '20150215'
                    'cmask'
                    'igram-20150215_20150227'
                    ...
                box         : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
                x/ystep     : int, number of pixels to pick/multilook for each output pixel
    Returns:    data        : 2D/3D array
                atr         : dict, metadata
    """
    # File Info: list of slice / dataset / dataset2d / dataset3d
    slice_list = get_slice_list(fname)
    ds_list = []
    for i in [i.split('-')[0] for i in slice_list]:
        if i not in ds_list:
            ds_list.append(i)
    ds_2d_list = [i for i in slice_list if '-' not in i]
    ds_3d_list = [i for i in ds_list if i not in ds_2d_list]

    # Input Argument: convert input datasetName into list of slice
    if not datasetName:
        datasetName = [ds_list[0]]
    elif isinstance(datasetName, str):
        datasetName = [datasetName]

    # if datasetName is all numerical (or include only T), add dsFamily as prefix
    if all(i.replace('T','').isdigit() for i in datasetName):
        datasetName = ['{}-{}'.format(ds_3d_list[0], i) for i in datasetName]

    # Input Argument: decompose slice list into dsFamily and inputDateList
    dsFamily = datasetName[0].split('-')[0]
    inputDateList = [i.replace(dsFamily,'').replace('-','') for i in datasetName]

    # read hdf5
    with h5py.File(fname, 'r') as f:
        # get dataset object
        dsNames = [i for i in [datasetName[0], dsFamily] if i in f.keys()]
        dsNamesOld = [i for i in slice_list if '/{}'.format(datasetName[0]) in i] # support for old mintpy files
        if len(dsNames) > 0:
            ds = f[dsNames[0]]
        elif len(dsNamesOld) > 0:
            ds = f[dsNamesOld[0]]
        else:
            raise ValueError('input dataset {} not found in file {}'.format(datasetName, fname))

        # 2D dataset
        if ds.ndim == 2:
            # read data
            data = ds[box[1]:box[3],
                      box[0]:box[2]]

            # sampling / nearest interplation in y/xstep
            if xstep * ystep > 1:
                # output size if x/ystep > 1
                xsize = int((box[2] - box[0]) / xstep)
                ysize = int((box[3] - box[1]) / ystep)

                # sampling
                data = data[int(ystep/2)::ystep,
                            int(xstep/2)::xstep]
                data = data[:ysize, :xsize]

        # 3D dataset
        elif ds.ndim == 3:
            # define flag matrix for index in time domain
            slice_flag = np.zeros((ds.shape[0]), dtype=np.bool_)
            if not inputDateList or inputDateList == ['']:
                slice_flag[:] = True
            else:
                date_list = [i.split('-')[1] for i in
                             [j for j in slice_list if j.startswith(dsFamily)]]
                for d in inputDateList:
                    slice_flag[date_list.index(d)] = True

            # read data
            if xstep * ystep == 1:
                data = ds[slice_flag,
                          box[1]:box[3],
                          box[0]:box[2]]

            else:
                # output size if x/ystep > 1
                xsize = int((box[2] - box[0]) / xstep)
                ysize = int((box[3] - box[1]) / ystep)

                # sampling / nearest interplation in y/xstep
                # use for loop to save memory
                num_slice = np.sum(slice_flag)
                data = np.zeros((num_slice, ysize, xsize), ds.dtype)

                inds = np.where(slice_flag)[0]
                for i in range(num_slice):
                    # print out msg
                    if print_msg:
                        sys.stdout.write('\r' + f'reading slice {i+1}/{num_slice}...')
                        sys.stdout.flush()

                    # read and index
                    d2 = ds[inds[i],
                            box[1]:box[3],
                            box[0]:box[2]]
                    d2 = d2[int(ystep/2)::ystep,
                            int(xstep/2)::xstep]
                    data[i, :, :] = d2[:ysize, :xsize]

                if print_msg:
                    print('')

            if any(i == 1 for i in data.shape):
                data = np.squeeze(data)

    return data


def read_binary_file(fname, datasetName=None, box=None, xstep=1, ystep=1):
    """Read data from binary file, such as .unw, .cor, etc.
    Parameters: fname       : str, path/name of binary file
                datasetName : str, dataset name for file with multiple bands of data
                    e.g.: incidenceAngle, azimuthAngle, rangeCoord, azimuthCoord, ...
                box         : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
                x/ystep     : int, number of pixels to pick/multilook for each output pixel
    Returns:    data        : 2D array in size of (length, width) in BYTE / int16 / float32 / complex64 / float64 etc.
                atr         : dict, metadata of binary file
    """
    # Basic Info
    fext = os.path.splitext(os.path.basename(fname))[1].lower()

    # metadata
    atr = read_attribute(fname, datasetName=datasetName)
    processor = atr['PROCESSOR']
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    if not box:
        box = (0, 0, width, length)

    # default data structure
    data_type = atr.get('DATA_TYPE', 'float32').lower()
    byte_order = atr.get('BYTE_ORDER', 'little-endian').lower()
    num_band = int(atr.get('number_bands', '1'))
    band_interleave = atr.get('scheme', 'BIL').upper()

    # default data to read
    band = 1
    cpx_band = 'phase'

    # ISCE
    if processor in ['isce']:
        # convert default short name for data type from ISCE
        dataTypeDict = {
            'byte': 'int8',
            'float': 'float32',
            'double': 'float64',
            'cfloat': 'complex64',
        }
        if data_type in dataTypeDict.keys():
            data_type = dataTypeDict[data_type]

        k = atr['FILE_TYPE'].lower().replace('.', '')
        if k in ['unw', 'cor']:
            band = min(2, num_band)
            if datasetName and datasetName in ['band1','intensity','magnitude']:
                band = 1

        elif k in ['slc']:
            cpx_band = 'magnitude'

        elif k.startswith('los') and datasetName and datasetName.startswith(('band2','az','head')):
            band = min(2, num_band)

        elif k in ['incLocal']:
            band = min(2, num_band)
            if datasetName and 'local' not in datasetName.lower():
                band = 1

        elif datasetName:
            if datasetName.lower() == 'band2':
                band = 2
            elif datasetName.lower() == 'band3':
                band = 3

            elif datasetName.startswith(('mag', 'amp')):
                cpx_band = 'magnitude'
            elif datasetName in ['phase', 'angle']:
                cpx_band = 'phase'
            elif datasetName.lower() == 'real':
                cpx_band = 'real'
            elif datasetName.lower().startswith('imag'):
                cpx_band = 'imag'
            elif datasetName.startswith(('cpx', 'complex')):
                cpx_band = 'complex'

            else:
                # flexible band list
                ds_list = get_slice_list(fname)
                if datasetName in ds_list:
                    band = ds_list.index(datasetName) + 1

        band = min(band, num_band)

    # ROI_PAC
    elif processor in ['roipac']:
        # data structure - auto
        band_interleave = 'BIL'
        byte_order = 'little-endian'

        # data structure - file specific based on file extension
        data_type = 'float32'
        num_band = 1

        if fext in ['.unw', '.cor', '.hgt', '.msk']:
            num_band = 2
            band = 2

        elif fext in ['.int']:
            data_type = 'complex64'

        elif fext in ['.amp']:
            data_type = 'complex64'
            cpx_band = 'magnitude'

        elif fext in ['.dem', '.wgs84']:
            data_type = 'int16'

        elif fext in ['.flg', '.byt']:
            data_type = 'bool_'

        elif fext in ['.trans']:
            num_band = 2
            if datasetName and datasetName.startswith(('az', 'azimuth')):
                band = 2

    # Gamma
    elif processor == 'gamma':
        # data structure - auto
        band_interleave = 'BIL'
        byte_order = atr.get('BYTE_ORDER', 'big-endian')

        data_type = 'float32'
        if fext in ['.unw', '.cor', '.hgt_sim', '.dem', '.amp', '.ramp']:
            pass

        elif fext in ['.int']:
            data_type = 'complex64'

        elif fext in ['.utm_to_rdc']:
            data_type = 'complex64'
            if datasetName and datasetName.startswith(('az', 'azimuth')):
                cpx_band = 'imag'
            else:
                cpx_band = 'real'

        elif fext == '.slc':
            data_type = 'complex32'
            cpx_band = 'magnitude'

        elif fext in ['.mli']:
            byte_order = 'little-endian'

    # SNAP
    # BEAM-DIMAP data format
    # https://www.brockmann-consult.de/beam/doc/help/general/BeamDimapFormat.html
    elif processor == 'snap':
        # data structure - auto
        band_interleave = atr.get('scheme', 'BSQ').upper()

        # byte order
        byte_order = atr.get('BYTE_ORDER', 'big-endian')
        if 'byte order' in atr.keys() and atr['byte order'] == '0':
            byte_order = 'little-endian'

    # GDAL / GMTSAR / ASF HyP3
    elif processor in ['gdal', 'gmtsar', 'hyp3']:
        pass

    else:
        print('Unknown InSAR processor: {}'.format(processor))

    # reading
    if processor in ['gdal', 'gmtsar', 'hyp3']:
        data = read_gdal(
            fname,
            box=box,
            band=band,
            cpx_band=cpx_band,
            xstep=xstep,
            ystep=ystep,
        )
    else:
        data = read_binary(
            fname,
            shape=(length, width),
            box=box,
            data_type=data_type,
            byte_order=byte_order,
            num_band=num_band,
            band_interleave=band_interleave,
            band=band,
            cpx_band=cpx_band,
            xstep=xstep,
            ystep=ystep,
        )

    if 'DATA_TYPE' not in atr:
        atr['DATA_TYPE'] = data_type
    return data, atr


#########################################################################
def get_slice_list(fname):
    """Get list of 2D slice existed in file (for display)"""
    fbase, fext = os.path.splitext(os.path.basename(fname))
    fext = fext.lower()
    atr = read_attribute(fname)
    k = atr['FILE_TYPE']

    global slice_list
    # HDF5 Files
    if fext in ['.h5', '.he5']:
        with h5py.File(fname, 'r') as f:
            d1_list = [i for i in f.keys() if isinstance(f[i], h5py.Dataset)]
        if k == 'timeseries' and k in d1_list:
            obj = timeseries(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        elif k in ['geometry'] and k not in d1_list:
            obj = geometry(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        elif k in ['ifgramStack']:
            obj = ifgramStack(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        elif k in ['HDFEOS']:
            obj = HDFEOS(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        elif k in ['giantTimeseries']:
            obj = giantTimeseries(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        elif k in ['giantIfgramStack']:
            obj = giantIfgramStack(fname)
            obj.open(print_msg=False)
            slice_list = obj.sliceList

        else:
            ## Find slice by walking through the file structure
            length, width = int(atr['LENGTH']), int(atr['WIDTH'])
            def get_hdf5_2d_dataset(name, obj):
                global slice_list
                if isinstance(obj, h5py.Dataset) and obj.shape[-2:] == (length, width):
                    if obj.ndim == 2:
                        slice_list.append(name)
                    elif obj.ndim == 3:
                        slice_list += ['{}-{}'.format(name, i+1) for i in range(obj.shape[0])]
                    else:
                        warnings.warn('file has un-defined {}D dataset: {}'.format(obj.ndim, name))
            slice_list = []
            with h5py.File(fname, 'r') as f:
                f.visititems(get_hdf5_2d_dataset)

    # Binary Files
    else:
        num_band = int(atr.get('number_bands', '1'))
        if fext in ['.trans', '.utm_to_rdc'] and num_band == 2:
            # roipac / gamma lookup table
            slice_list = ['rangeCoord', 'azimuthCoord']

        elif fbase.startswith('los') and num_band == 2:
            # isce los file
            slice_list = ['incidenceAngle', 'azimuthAngle']

        elif fext in ['.int', '.unw']:
            # do not check the actual num_band in order to support
            # mag / pha / cpx reading like "multiple bands"
            slice_list = ['magnitude', 'phase']

        elif fbase.startswith('offset') and fext in ['.bip'] and num_band == 2:
            slice_list = ['azimuthOffset', 'rangeOffset']

        elif fbase.startswith('offset') and fname.endswith('cov.bip') and num_band == 3:
            slice_list = ['azimuthOffsetVar', 'rangeOffsetVar', 'offsetCovar']

        else:
            slice_list = ['band{}'.format(i+1) for i in range(num_band)]

    return slice_list


def get_dataset_list(fname, datasetName=None):
    """Get list of 2D and 3D dataset to facilitate systematic file reading"""
    if datasetName:
        return [datasetName]

    fext = os.path.splitext(fname)[1].lower()

    global ds_list
    if fext in ['.h5', '.he5']:
        atr = read_attribute(fname)
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])

        def get_hdf5_dataset(name, obj):
            global ds_list
            if isinstance(obj, h5py.Dataset) and obj.shape[-2:] == (length, width):
                ds_list.append(name)
        ds_list = []
        with h5py.File(fname, 'r') as f:
            f.visititems(get_hdf5_dataset)

    else:
        ds_list = get_slice_list(fname)

    return ds_list


def get_hdf5_compression(fname):
    """Get the compression type of input HDF5 file"""
    ext = os.path.splitext(fname)[1].lower()
    if ext not in ['.h5','.he5']:
        return None

    compression = None
    ds_name = get_dataset_list(fname)[0]
    with h5py.File(fname, 'r') as f:
        compression = f[ds_name].compression
    return compression


#########################################################################
def read_attribute(fname, datasetName=None, metafile_ext=None):
    """Read attributes of input file into a dictionary
    Parameters: fname : str, path/name of data file
                datasetName : str, name of dataset of interest, for file with multiple datasets
                    e.g. unwrapPhase in ifgramStack.h5
                         coherence   in ifgramStack.h5
                         height      in geometryRadar.h5
                         latitude    in geometryRadar.h5
                         ...
    Returns:    atr : dict, attributes dictionary
    """
    fbase, fext = os.path.splitext(os.path.basename(fname))
    fext = fext.lower()
    if not os.path.isfile(fname):
        msg = 'input file not existed: {}\n'.format(fname)
        msg += 'current directory: '+os.getcwd()
        raise Exception(msg)

    # HDF5 files
    if fext in ['.h5', '.he5']:
        with h5py.File(fname, 'r') as f:
            atr = dict(f.attrs)
            g1_list = [i for i in f.keys() if isinstance(f[i], h5py.Group)]
            d1_list = [i for i in f.keys() if isinstance(f[i], h5py.Dataset) and f[i].ndim >= 2]

        # FILE_TYPE - k
        py2_mintpy_stack_files = ['interferograms', 'coherence', 'wrapped'] #obsolete mintpy format
        if any(i in d1_list for i in ['unwrapPhase', 'azimuthOffset']):
            k = 'ifgramStack'
        elif any(i in d1_list for i in ['height', 'latitude', 'azimuthCoord']):
            k = 'geometry'
        elif any(i in g1_list+d1_list for i in ['timeseries', 'displacement']):
            k = 'timeseries'
        elif any(i in d1_list for i in ['velocity']):
            k = 'velocity'
        elif 'HDFEOS' in g1_list:
            k = 'HDFEOS'
        elif 'recons' in d1_list:
            k = 'giantTimeseries'
        elif any(i in d1_list for i in ['igram', 'figram']):
            k = 'giantIfgramStack'
        elif any(i in g1_list for i in py2_mintpy_stack_files):
            k = list(set(g1_list) & set(py2_mintpy_stack_files))[0]
        elif len(d1_list) > 0:
            k = d1_list[0]
        elif len(g1_list) > 0:
            k = g1_list[0]
        else:
            raise ValueError('unrecognized file type: '+fname)

        # metadata dict
        if k == 'giantTimeseries':
            atr = giantTimeseries(fname).get_metadata()
        elif k == 'giantIfgramStack':
            atr = giantIfgramStack(fname).get_metadata()

        elif len(atr) > 0 and 'WIDTH' in atr.keys():
            # use the attribute at root level, which is already read from the begining

            # grab attribute of dataset if specified, e.g. UNIT, no-data value, etc.
            if datasetName and datasetName in d1_list:
                with h5py.File(fname, 'r') as f:
                    atr.update(dict(f[datasetName].attrs))

        else:
            # otherwise, grab the list of attrs in HDF5 file
            # and use the attrs with most items
            global atr_list
            def get_hdf5_attrs(name, obj):
                global atr_list
                if len(obj.attrs) > 0 and 'WIDTH' in obj.attrs.keys():
                    atr_list.append(dict(obj.attrs))

            atr_list = []
            with h5py.File(fname, 'r') as f:
                f.visititems(get_hdf5_attrs)

            # use the attrs with most items
            if atr_list:
                num_list = [len(i) for i in atr_list]
                atr = atr_list[np.argmax(num_list)]
            else:
                raise ValueError('No attribute WIDTH found in file:', fname)

        # decode string format
        for key, value in atr.items():
            try:
                atr[key] = value.decode('utf8')
            except:
                atr[key] = value

        # attribute identified by MintPy
        # 1. FILE_TYPE
        atr['FILE_TYPE'] = str(k)

        # 2. DATA_TYPE
        ds = None
        with h5py.File(fname, 'r') as f:
            if datasetName and datasetName in f.keys():
                # get the dataset in the root level
                ds = f[datasetName]
            else:
                # get the 1st dataset in deeper levels
                global ds_list
                def get_hdf5_dataset(name, obj):
                    global ds_list
                    if isinstance(obj, h5py.Dataset) and obj.ndim >= 2:
                        ds_list.append(obj)
                ds_list = []
                f.visititems(get_hdf5_dataset)
                if ds_list:
                    ds = ds_list[0]

            if ds is not None:
                atr['DATA_TYPE'] = str(ds.dtype)

        # 3. PROCESSOR
        if 'INSAR_PROCESSOR' in atr.keys():
            atr['PROCESSOR'] = atr['INSAR_PROCESSOR']
        if 'PROCESSOR' not in atr.keys():
            atr['PROCESSOR'] = 'mintpy'

    elif fext == '.dehm':
        # 10 m Digital Ellipsoidal Height Model files from GSI
        atr = {}
        atr['LENGTH'] = 6000           # 40 mins in latitude  per grid
        atr['WIDTH']  = 9000           # 60 mins in longitude per grid
        atr['Y_STEP'] = - 0.4 / 3600.  # degree
        atr['X_STEP'] =   0.4 / 3600.  # degree
        atr['Y_UNIT'] = 'degrees'
        atr['X_UNIT'] = 'degrees'

        # Y/X_FIRST based on the naming convention
        yy, xx = float(fbase[:2]), float(fbase[2:])
        atr['Y_FIRST'] = (yy + 1.) / 1.5
        atr['X_FIRST'] = xx + 100.

        atr['PROCESSOR'] = 'GSI'
        atr['FILE_TYPE'] = fext
        atr['DATA_TYPE'] = 'float32'
        atr['PROJECTION'] = 'LATLON'
        atr['GEODETIC_DATUM'] = 'WGS84'
        atr['UNIT'] = 'm'

        # check file size for potential 5m DEHM data
        if os.path.getsize(fname) != atr['LENGTH'] * atr['WIDTH'] * 4:
            msg = 'input DEHM file size do NOT match with the pre-defined 10m DEHM: '
            msg += '{} * {} in {}!'.format(atr['LENGTH'], atr['WIDTH'], atr['DATA_TYPE'])
            raise ValueError(msg)

    else:
        # grab all existed potential metadata file given the data file in prefered order/priority
        # .aux.xml file does not have geo-coordinates info
        # .vrt file (e.g. incLocal.rdr.vrt from isce) does not have band interleavee info
        metafiles = [
            fname + '.rsc',
            fname + '.xml',
            fname + '.par',
            os.path.splitext(fname)[0] + '.hdr',
            fname + '.vrt',
            fname + '.aux.xml',
        ]
        metafiles = [i for i in metafiles if os.path.isfile(i)]

        # use metadata files with the specified extension if requested
        if metafile_ext:
            metafiles = [i for i in metafiles if i.endswith(metafile_ext)]

        # use the GDAL supported data file is no metadata file found
        if len(metafiles) == 0:
            # for .tif/.grd files, extract metadata from the file itself
            if fext in GDAL_FILE_EXTS:
                metafiles = [fname]
            else:
                raise FileNotFoundError('No metadata file found for data file: {}'.format(fname))

        atr = {}
        # PROCESSOR
        if fname.endswith('.img') and any(i.endswith('.hdr') for i in metafiles):
            atr['PROCESSOR'] = 'snap'

        elif any(i.endswith(('.xml', '.hdr', '.vrt')) for i in metafiles):
            atr['PROCESSOR'] = 'isce'
            xml_files = [i for i in metafiles if i.endswith('.xml')]
            if len(xml_files) > 0:
                atr.update(read_isce_xml(xml_files[0]))

        elif any(i.endswith('.par') for i in metafiles):
            atr['PROCESSOR'] = 'gamma'

        elif any(i.endswith('.rsc') for i in metafiles):
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'roipac'

        elif fext in GDAL_FILE_EXTS:
            atr['PROCESSOR'] = 'gdal'

        if 'PROCESSOR' not in atr.keys():
            atr['PROCESSOR'] = 'mintpy'

        # Read metadata file and FILE_TYPE
        metafile = metafiles[0]
        meta_ext = os.path.splitext(metafile)[1]

        # ignore certain meaningless file extensions
        while fext in ['.geo', '.rdr', '.full', '.wgs84', '.grd']:
            fbase, fext = os.path.splitext(fbase)
        if not fext:
            fext = fbase

        if meta_ext == '.rsc':
            atr.update(read_roipac_rsc(metafile))
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = fext

        elif meta_ext == '.xml':
            atr.update(read_isce_xml(metafile))
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = fext

        elif meta_ext == '.par':
            atr.update(read_gamma_par(metafile))
            atr['FILE_TYPE'] = fext

        elif meta_ext == '.hdr':
            atr.update(read_envi_hdr(metafile))

            # both snap and isce produce .hdr file
            # grab file type based on their different naming conventions
            if atr['PROCESSOR'] == 'snap':
                fbase = os.path.basename(fname).lower()
                if fbase.startswith('unw'):
                    atr['FILE_TYPE'] = '.unw'
                elif fbase.startswith(('coh','cor')):
                    atr['FILE_TYPE'] = '.cor'
                elif fbase.startswith('phase_ifg'):
                    atr['FILE_TYPE'] = '.int'
                elif 'dem' in fbase:
                    atr['FILE_TYPE'] = 'dem'
                else:
                    atr['FILE_TYPE'] = atr['file type']
            else:
                atr['FILE_TYPE'] = fext

        elif meta_ext in ['.vrt'] + GDAL_FILE_EXTS:
            atr.update(read_gdal_vrt(metafile))
            atr['FILE_TYPE'] = fext

        # DATA_TYPE for ISCE products
        dataTypeDict = {
            'byte': 'int8',
            'float': 'float32',
            'double': 'float64',
            'cfloat': 'complex64',
        }
        data_type = atr.get('DATA_TYPE', 'none').lower()
        if data_type != 'none' and data_type in dataTypeDict.keys():
            atr['DATA_TYPE'] = dataTypeDict[data_type]

    # UNIT
    k = atr['FILE_TYPE'].replace('.', '')
    if k == 'ifgramStack':
        if datasetName and datasetName in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[datasetName]
        else:
            atr['UNIT'] = 'radian'

    elif datasetName and datasetName in datasetUnitDict.keys():
        atr['UNIT'] = datasetUnitDict[datasetName]

    elif 'UNIT' not in atr.keys():
        if k in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[k]
        else:
            atr['UNIT'] = '1'

    # FILE_PATH
    if 'FILE_PATH' in atr.keys() and 'OG_FILE_PATH' not in atr.keys():
        # need to check original source file to successfully subset legacy-sensor products
        atr['OG_FILE_PATH'] = atr['FILE_PATH']
    atr['FILE_PATH'] = os.path.abspath(fname)

    atr = standardize_metadata(atr)

    return atr


def standardize_metadata(metaDictIn, standardKeys=None):
    """Convert metadata input ROI_PAC/MintPy format (for metadata with the same values)."""
    if standardKeys is None:
        standardKeys = standardMetadataKeys

    # make a copy
    metaDict = dict()
    for key, value in iter(metaDictIn.items()):
        metaDict[key] = value

    # get potential keys to match
    in_keys = [i for i in metaDict.keys() if i not in standardKeys.keys()]
    std_keys = [i for i in standardKeys.keys() if i not in metaDict.keys()]

    # loop to find match and assign values
    for std_key in std_keys:
        cand_keys = standardKeys[std_key]
        cand_keys = [i for i in cand_keys if i in in_keys]
        if len(cand_keys) > 0:
            metaDict[std_key] = metaDict[cand_keys[0]]

    return metaDict


#########################################################################
def read_template(fname, delimiter='=', print_msg=True):
    """Reads the template file into a python dictionary structure.
    Parameters: fname : str
                    full path to the template file
                delimiter : str
                    string to separate the key and value
                print_msg : bool
                    print message or not
    Returns:    template_dict : dict
                    file content
    Examples:
        tmpl = read_template(KyushuT424F610_640AlosA.template)
        tmpl = read_template(R1_54014_ST5_L0_F898.000.pi, ':')
    """
    template_dict = {}

    # insarmaps: the below logic for plotattributes object can be made much more simple
    # if we assume that any plot attribute coming after a > belongs to the
    # same object. Must Ask Falk and Yunjun if we can assume this to eliminate
    # all these conditionals
    plotAttributeDict = {}
    insidePlotObject = False
    plotAttributes = []

    def is_plot_attribute(attribute):
        tokens = attribute.split(".")
        if tokens is None:
            return False
        return tokens[0] == "plot" and len(tokens) > 1

    # read input text file or string
    lines = None
    if os.path.isfile(fname):
        with open(fname, 'r') as f:
            lines = f.readlines()
    elif isinstance(fname, str):
        lines = fname.split('\n')

    # loop to parser/read each line
    for line in lines:
        line = line.strip()
        # split on the 1st occurrence of delimiter
        c = [i.strip() for i in line.split(delimiter, 1)]

        # ignore commented lines or those without variables
        if len(c) < 2 or line.startswith(('%', '#', '!')):
            # next

            # insarmaps:
            if line.startswith(">"):
                plotAttributeDict = {}
                insidePlotObject = True
            # otherwise, if previously inside attributes object, we are now outside
            # unless the line is a comment
            elif insidePlotObject and not line.startswith('%') and not line.startswith('#'):
                # just came from being inside plot object, but now we are outside
                insidePlotObject = False
                plotAttributes.append(plotAttributeDict)

        else:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].strip()
            value = os.path.expanduser(value)  # translate ~ symbol
            value = os.path.expandvars(value)  # translate env variables

            if value != '':
                template_dict[key] = value

            # insarmaps:
            if insidePlotObject:
                if is_plot_attribute(key):
                    plotAttributeDict[key] = value
                else:
                    # just came from being inside plot object, but now we are outside
                    insidePlotObject = False
                    plotAttributes.append(plotAttributeDict)
                    template_dict[key] = value

    # insarmaps: what if no \n at end of file? write out last plot attributes dict
    if insidePlotObject:
        plotAttributes.append(plotAttributeDict)
    if len(plotAttributes) > 0:
        template_dict["plotAttributes"] = json.dumps(plotAttributes)

    return template_dict


def read_roipac_rsc(fname, delimiter=' '):
    """Read ROI_PAC .rsc file into a python dict structure.
    Parameters: fname : str.
                    File path of .rsc file.
    Returns:    rscDict : dict
                    Dictionary of keys and values in RSC file.
    Examples:
        from mintpy.utils import readfile
        atr = readfile.read_roipac_rsc('filt_101120_110220_c10.unw.rsc')
    """
    # read .rsc file
    with open(fname, 'r') as f:
        lines = f.readlines()

    # convert list of str into dict
    rscDict = {}
    for line in lines:
        c = [i.strip() for i in line.strip().replace('\t',' ').split(delimiter, 1)]
        key = c[0]
        value = c[1].replace('\n', '').strip()
        rscDict[key] = value

    rscDict = standardize_metadata(rscDict)

    return rscDict


def read_gamma_par(fname, delimiter=':', skiprows=3):
    """Read GAMMA .par/.off file into a python dict structure.
    Parameters: fname : str.
                    File path of .par, .off file.
                delimiter : str, optional
                    String used to separate values.
                skiprows : int, optional
                    Skip the first skiprows lines.
    Returns:    parDict : dict
                    Attributes dictionary
    """
    # Read txt file
    with open(fname, 'r') as f:
        lines = f.readlines()[skiprows:]

    # convert list of str into dict
    parDict = {}
    for line in lines:
        c = [i.strip() for i in line.strip().split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%', '#')):
            next
        else:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].split()[0].strip()
            parDict[key] = value

    parDict = attribute_gamma2roipac(parDict)
    parDict = standardize_metadata(parDict)

    return parDict


def attribute_gamma2roipac(par_dict_in):
    """Convert Gamma metadata into ROI_PAC/MintPy format."""
    par_dict = dict()
    for key, value in iter(par_dict_in.items()):
        par_dict[key] = value

    # LENGTH - number of rows
    for key in par_dict_in.keys():
        if any(key.startswith(i) for i in ['azimuth_lines',
                                           'nlines',
                                           'az_samp',
                                           'interferogram_azimuth_lines']):
            par_dict['LENGTH'] = par_dict[key]

    # WIDTH - number of columns
    for key in par_dict_in.keys():
        if any(key.startswith(i) for i in ['width',
                                           'range_samp',
                                           'interferogram_width']):
            par_dict['WIDTH'] = par_dict[key]

    # radar_frequency -> WAVELENGTH
    key = 'radar_frequency'
    if key in par_dict_in.keys():
        value = float(par_dict[key])
        par_dict['WAVELENGTH'] = str(SPEED_OF_LIGHT / value)

    # sar_to_earth_center/earth_radius_below_sensor -> HEIGHT/EARTH_RADIUS
    key = 'earth_radius_below_sensor'
    if key in par_dict_in.keys():
        Re = float(par_dict[key])
        par_dict['EARTH_RADIUS'] = str(Re)

        key2 = 'sar_to_earth_center'
        if key2 in par_dict_in.keys():
            value = float(par_dict[key2])
            par_dict['HEIGHT'] = str(value - Re)

    # sensor -> PLATFORM
    key = 'sensor'
    if key in par_dict_in.keys():
        par_dict['PLATFORM'] = par_dict[key]

    # heading -> ORBIT_DIRECTION
    key = 'heading'
    if key in par_dict_in.keys():
        value = float(par_dict[key])
        if (270 < value < 360) or (-90 < value < 90):
            par_dict['ORBIT_DIRECTION'] = 'ascending'
        else:
            par_dict['ORBIT_DIRECTION'] = 'descending'

    # azimuth_angle -> ANTENNA_SIDE
    key = 'azimuth_angle'
    if key in par_dict_in.keys():
        value = float(par_dict[key])
        if 0 < value < 180:
            par_dict['ANTENNA_SIDE'] = '-1'
        else:
            par_dict['ANTENNA_SIDE'] = '1'

    return par_dict


def read_isce_xml(fname):
    """Read ISCE .xml file into a python dict structure."""
    root = ET.parse(fname).getroot()
    xmlDict = {}

    # imageFile, e.g. filt_fine.unw.xml
    if root.tag.startswith('image'):
        for child in root.findall('property'):
            key = child.get('name')
            value = child.find('value').text
            xmlDict[key] = value

        # Read lat/lon info for geocoded file
        # in form: root/component coordinate*/property name/value
        for coord_name, prefix in zip(['coordinate1', 'coordinate2'], ['X', 'Y']):
            child = root.find("./component[@name='{}']".format(coord_name))
            v_step  = float(child.find("./property[@name='delta']").find('value').text)
            v_first = float(child.find("./property[@name='startingvalue']").find('value').text)
            if abs(v_step) < 1. and abs(v_step) > 1e-7:
                xmlDict['{}_STEP'.format(prefix)] = v_step
                xmlDict['{}_FIRST'.format(prefix)] = v_first - v_step / 2.
                xmlDict['{}_UNIT'.format(prefix)] = 'degrees'

    # PAMDataset, e.g. hgt.rdr.aux.xml
    elif root.tag == 'PAMDataset':
        meta = root.find("./Metadata[@domain='ENVI']")
        for child in meta.findall("MDI"):
            key = child.get('key')
            value = child.text
            xmlDict[key] = value
        xmlDict['data_type'] = ENVI2NUMPY_DATATYPE[xmlDict['data_type']]

    xmlDict = standardize_metadata(xmlDict)

    return xmlDict


def read_envi_hdr(fname):
    """Read ENVI .hdr file into a python dict structure"""
    atr = read_template(fname, delimiter='=')
    atr['DATA_TYPE'] = ENVI2NUMPY_DATATYPE[atr.get('data type', '4')]
    atr['BYTE_ORDER'] = ENVI_BYTE_ORDER[atr.get('byte order', '1')]

    if 'map info' in atr.keys():
        map_info = [i.replace('{','').replace('}','').strip() for i in atr['map info'].split(',')]
        x_step = abs(float(map_info[5]))
        y_step = abs(float(map_info[6])) * -1.
        #unit = map_info[-1].replace('}','').split('=')[1].lower()

        if abs(x_step) < 1. and abs(x_step) > 1e-7:
            atr['X_UNIT'] = 'degrees'
            atr['Y_UNIT'] = 'degrees'
            atr['X_STEP'] = str(x_step)
            atr['Y_STEP'] = str(y_step)
            atr['X_FIRST'] = str(float(map_info[3]) - x_step / 2.)
            atr['Y_FIRST'] = str(float(map_info[4]) - y_step / 2.)

    atr = standardize_metadata(atr)

    return atr


def read_gdal_vrt(fname):
    """Read GDAL .vrt file into a python dict structure using gdal

    Modified from $ISCE_HOME/applications/gdal2isce_xml.gdal2isce_xml() written by David Bekaert.
    """
    try:
        from osgeo import gdal, osr
    except ImportError:
        raise ImportError('Cannot import gdal and osr!')

    # read dataset using gdal
    ds = gdal.Open(fname, gdal.GA_ReadOnly)

    atr = {}
    atr['WIDTH']  = ds.RasterXSize
    atr['LENGTH'] = ds.RasterYSize
    atr['number_bands'] = ds.RasterCount

    # data type
    data_type = ds.GetRasterBand(1).DataType
    atr['DATA_TYPE'] = GDAL2ISCE_DATATYPE[data_type]

    # interleave
    scheme = ds.GetMetadata('IMAGE_STRUCTURE').get('INTERLEAVE', 'PIXEL')
    atr['scheme'] = ENVI_BAND_INTERLEAVE[scheme]

    # transformation contains gridcorners
    # (lines/pixels or lonlat and the spacing 1/-1 or deltalon/deltalat)
    transform = ds.GetGeoTransform()
    x0 = transform[0]
    y0 = transform[3]
    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.

    atr['X_STEP'] = x_step
    atr['Y_STEP'] = y_step
    atr['X_FIRST'] = x0 - x_step / 2.
    atr['Y_FIRST'] = y0 - y_step / 2.

    # projection / coordinate unit
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs_name = srs.GetName()
    if 'UTM' in srs_name:
        atr['UTM_ZONE'] = srs_name.split('UTM zone')[-1].strip()
        atr['X_UNIT'] = 'meters'
        atr['Y_UNIT'] = 'meters'

    elif abs(x_step) < 1. and abs(x_step) > 1e-7:
        atr['X_UNIT'] = 'degrees'
        atr['Y_UNIT'] = 'degrees'
        # constrain longitude within (-180, 180]
        if atr['X_FIRST'] > 180.:
            atr['X_FIRST'] -= 360.

    atr = standardize_metadata(atr)

    return atr


def read_gmtsar_prm(fname, delimiter='='):
    """Read GMTSAR .prm file into a python dict structure.
    Parameters: fname : str.
                    File path of .rsc file.
    Returns:    prmDict : dict
                    Dictionary of keys and values in the PRM file.
    """
    # read .prm file
    with open(fname, 'r') as f:
        lines = f.readlines()

    # convert list of str into dict
    prmDict = {}
    for line in lines:
        c = [i.strip() for i in line.strip().replace('\t',' ').split(delimiter, 1)]
        key = c[0]
        value = c[1].replace('\n', '').strip()
        prmDict[key] = value

    prmDict = attribute_gmtsar2roipac(prmDict)
    prmDict = standardize_metadata(prmDict)

    return prmDict


def attribute_gmtsar2roipac(prm_dict_in):
    """Convert GMTSAR metadata into ROI_PAC/MintPy format (for metadata with different values)."""
    prm_dict = dict()
    for key, value in iter(prm_dict_in.items()):
        prm_dict[key] = value

    # lookdir -> ANTENNA_SIDE
    key = 'lookdir'
    if key in prm_dict_in.keys():
        value = prm_dict[key]
        if value.upper() == 'R':
            prm_dict['ANTENNA_SIDE'] = '-1'
        else:
            prm_dict['ANTENNA_SIDE'] = '1'

    # SC_vel -> AZIMUTH_PIXEL_SIZE (in single look)
    key = 'SC_vel'
    if key in prm_dict_in.keys():
        vel = float(prm_dict[key])
        Re = float(prm_dict['earth_radius'])
        PRF = float(prm_dict['PRF'])
        height = float(prm_dict['SC_height'])
        az_pixel_size = vel / PRF * Re / (Re + height)
        prm_dict['AZIMUTH_PIXEL_SIZE'] = az_pixel_size

    # rng_samp_rate -> RANGE_PIXEL_SIZE (in single look)
    key = 'rng_samp_rate'
    if key in prm_dict_in.keys():
        value = float(prm_dict[key])
        prm_dict['RANGE_PIXEL_SIZE'] = SPEED_OF_LIGHT / value / 2.0

    # SC_clock_start/stop -> CENTER_LINE_TUC
    dt_center = (float(prm_dict['SC_clock_start']) + float(prm_dict['SC_clock_stop'])) / 2.0
    t_center = dt_center - int(dt_center)
    prm_dict['CENTER_LINE_UTC'] = str(t_center * 24. * 60. * 60.)

    return prm_dict



#########################################################################
def read_binary(fname, shape, box=None, data_type='float32', byte_order='l',
                num_band=1, band_interleave='BIL', band=1, cpx_band='phase',
                xstep=1, ystep=1):
    """Read binary file using np.fromfile.

    Parameters: fname : str, path/name of data file to read
                shape : tuple of 2 int in (length, width)
                box   : tuple of 4 int in (x0, y0, x1, y1)
                data_type : str, data type of stored array, e.g.:
                    bool_
                    int8, int16, int32
                    float16, float32, float64
                    complex64, complex128
                byte_order      : str, little/big-endian
                num_band        : int, number of bands
                band_interleave : str, band interleaving scheme, e.g.:
                    BIP
                    BIL
                    BSQ
                band     : int, band of interest, between 1 and num_band.
                cpx_band : str, e.g.:
                    real,
                    imag, imaginary
                    phase,
                    mag, magnitude
                    cpx
                x/ystep  : int, number of pixels to pick/multilook for each output pixel
    Returns:    data     : 2D np.array
    Examples:   # ISCE files
                atr = read_attribute(fname)
                shape = (int(atr['LENGTH']), int(atr['WIDTH']))
                data = read_binary('filt_fine.unw', shape, num_band=2, band=2)
                data = read_binary('filt_fine.cor', shape)
                data = read_binary('filt_fine.int', shape, data_type='complex64', cpx_band='phase')
                data = read_binary('burst_01.slc', shape, data_type='complex64', cpx_band='mag')
                data = read_binary('los.rdr', shape, num_band=2, band=1)
    """
    length, width = shape
    if not box:
        box = (0, 0, width, length)

    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        letter, digit = re.findall('(\d+|\D+)', data_type)
        # convert into short style: float32 --> c4
        if len(letter) > 1:
            letter = letter[0]
            digit = int(int(digit) / 8)
        data_type = '>{}{}'.format(letter, digit)

    # read data
    band_interleave = band_interleave.upper()
    if band_interleave == 'BIL':
        data = np.fromfile(fname,
                           dtype=data_type,
                           count=box[3]*width*num_band).reshape(-1, width*num_band)
        data = data[box[1]:box[3],
                    width*(band-1)+box[0]:width*(band-1)+box[2]]

    elif band_interleave == 'BIP':
        data = np.fromfile(fname,
                           dtype=data_type,
                           count=box[3]*width*num_band).reshape(-1, width*num_band)
        data = data[box[1]:box[3],
                    np.arange(box[0], box[2])*num_band+band-1]

    elif band_interleave == 'BSQ':
        data = np.fromfile(fname,
                           dtype=data_type,
                           count=(box[3]+length*(band-1))*width).reshape(-1, width)
        data = data[length*(band-1)+box[1]:length*(band-1)+box[3],
                    box[0]:box[2]]
    else:
        raise ValueError('unrecognized band interleaving:', band_interleave)

    # adjust output band for complex data
    if data_type.replace('>', '').startswith('c'):
        if cpx_band.startswith('real'):
            data = data.real
        elif cpx_band.startswith('imag'):
            data = data.imag
        elif cpx_band.startswith('pha'):
            data = np.angle(data)
        elif cpx_band.startswith('mag'):
            data = np.absolute(data)
        elif cpx_band.startswith(('cpx', 'complex')):
            pass
        else:
            raise ValueError('unrecognized complex band:', cpx_band)

    # skipping/multilooking
    if xstep * ystep > 1:
        # output size if x/ystep > 1
        xsize = int((box[2] - box[0]) / xstep)
        ysize = int((box[3] - box[1]) / ystep)

        # sampling
        data = data[int(ystep/2)::ystep,
                    int(xstep/2)::xstep]
        data = data[:ysize, :xsize]

    return data


def read_gdal(fname, box=None, band=1, cpx_band='phase', xstep=1, ystep=1):
    """Read binary data file using gdal.

    Parameters: fname    : str, path/name of data file to read
                box      : tuple of 4 int in (x0, y0, x1, y1)
                band     : int, band of interest, between 1 and num_band.
                cpx_band : str, e.g.: real, imag, phase, mag, cpx
                x/ystep  : int, number of pixels to pick/multilook for each output pixel
    Returns:    data     : 2D np.array
    """
    try:
        from osgeo import gdal
    except ImportError:
        raise ImportError('Cannot import gdal!')

    # open data file
    ds = gdal.Open(fname, gdal.GA_ReadOnly)
    bnd = ds.GetRasterBand(band)

    # box
    if not box:
        box = (0, 0, ds.RasterXSize, ds.RasterYSize)

    # read
    # link: https://gdal.org/python/osgeo.gdal.Band-class.html#ReadAsArray
    kwargs = dict(xoff=box[0],
                  yoff=box[1],
                  win_xsize=box[2]-box[0],
                  win_ysize=box[3]-box[1])
    data = bnd.ReadAsArray(**kwargs)

    # adjust output band for complex data
    data_type = GDAL2ISCE_DATATYPE[bnd.DataType]
    if data_type.replace('>', '').startswith('c'):
        if cpx_band.startswith('real'):
            data = data.real
        elif cpx_band.startswith('imag'):
            data = data.imag
        elif cpx_band.startswith('pha'):
            data = np.angle(data)
        elif cpx_band.startswith('mag'):
            data = np.absolute(data)
        elif cpx_band.startswith(('cpx', 'complex')):
            pass
        else:
            raise ValueError('unrecognized complex band:', cpx_band)

    # skipping/multilooking
    if xstep * ystep > 1:
        # output size if x/ystep > 1
        xsize = int((box[2] - box[0]) / xstep)
        ysize = int((box[3] - box[1]) / ystep)

        # sampling
        data = data[int(ystep/2)::ystep,
                    int(xstep/2)::xstep]
        data = data[:ysize, :xsize]

    return data


############################ Obsolete Functions ###############################
def read_float32(fname, box=None, byte_order='l'):
    """Reads roi_pac data (RMG format, interleaved line by line).

    ROI_PAC file: .unw, .cor, .hgt, .trans, .msk

    RMG format (named after JPL radar pionner Richard M. Goldstein): made
    up of real*4 numbers in two arrays side-by-side. The two arrays often
    show the magnitude of the radar image and the phase, although not always
    (sometimes the phase is the correlation). The length and width of each
    array are given as lines in the metadata (.rsc) file. Thus the total
    width width of the binary file is (2*width) and length is (length), data
    are stored as:
    magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    magnitude, magnitude, magnitude, ...,phase, phase, phase, ...
    ......

       box  : 4-tuple defining the left, upper, right, and lower pixel coordinate.
    Example:
       a,p,r = read_float32('100102-100403.unw')
       a,p,r = read_float32('100102-100403.unw',(100,1200,500,1500))
    """

    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'rmg'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f4'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>f4'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*2*width).reshape(box[3], 2*width)
    amplitude = data[box[1]:box[3],
                     box[0]:box[2]]
    phase = data[box[1]:box[3],
                 width+box[0]:width+box[2]]

    return amplitude, phase, atr


def read_real_float64(fname, box=None, byte_order='l'):
    """Read real float64/double data matrix, i.e. isce lat/lon.rdr
    """
    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'float64'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f8'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>f8'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3],
                box[0]:box[2]]
    return data, atr


def read_complex_float32(fname, box=None, byte_order='l', band='phase'):
    """Read complex float 32 data matrix, i.e. roi_pac int or slc data.
    old name: read_complex64()

    ROI_PAC file: .slc, .int, .amp

    Data is sotred as:
    real, imaginary, real, imaginary, ...
    real, imaginary, real, imaginary, ...
    ...

    Parameters: fname : str,
                    input file name
                box : 4-tuple
                    defining (left, upper, right, lower) pixel coordinate.
                byte_order : str, optional
                    order of reading byte in the file
                band : str
                    output format, default = phase
                    phase, amplitude, real, imag, complex
    Returns: data : 2D np.array in complex float32
    Example:
        amp, phase, atr = read_complex_float32('geo_070603-070721_0048_00018.int')
        data, atr       = read_complex_float32('150707.slc', 1)
    """

    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'complex64'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'c8'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>c8'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3],
                box[0]:box[2]]

    if band == 'phase':
        dataOut = np.angle(data)
    elif band == 'amplitude':
        dataOut = np.absolute(data)
    elif band == 'real':
        dataOut = data.real
    elif band == 'imag':
        dataOut = data.imag
    else:
        dataOut = data

    return dataOut, atr


def read_real_float32(fname, box=None, byte_order='l'):
    """Read real float 32 data matrix, i.e. GAMMA .mli file
    Parameters: fname     : str, path, filename to be read
                byte_order : str, optional, order of reading byte in the file
    Returns: data : 2D np.array, data matrix
             atr  : dict, attribute dictionary
    Usage: data, atr = read_real_float32('20070603.mli')
           data, atr = read_real_float32('diff_filt_130118-130129_4rlks.unw')
    """
    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'float32'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'f4'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>f4'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3],
                box[0]:box[2]]
    return data, atr


def read_complex_int16(fname, box=None, byte_order='l', cpx=False):
    """Read complex int 16 data matrix, i.e. GAMMA SCOMPLEX file (.slc)

    Gamma file: .slc

    Inputs:
       file: complex data matrix (cpx_int16)
       box: 4-tuple defining the left, upper, right, and lower pixel coordinate.
    Example:
       data,rsc = read_complex_int16('100102.slc')
       data,rsc = read_complex_int16('100102.slc',(100,1200,500,1500))
    """

    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'complex32'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'i2'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>i2'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*2*width).reshape(box[3], 2*width)
    data = data[box[1]:box[3],
                2*box[0]:2*box[2]].flatten()

    odd_idx = np.arange(1, len(data), 2)
    real = data[odd_idx-1].reshape(box[3]-box[1],
                                   box[2]-box[0])
    imag = data[odd_idx].reshape(box[3]-box[1],
                                 box[2]-box[0])

    if cpx:
        return real, imag, atr
    else:
        amplitude = np.hypot(imag, real)
        phase = np.arctan2(imag, real)
        return amplitude, phase, atr


def read_real_int16(fname, box=None, byte_order='l'):
    atr = read_attribute(fname)
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'int16'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data_type = 'i2'
    if byte_order in ['b', 'big', 'big-endian', 'ieee-be']:
        data_type = '>i2'

    data = np.fromfile(fname,
                       dtype=data_type,
                       count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3],
                box[0]:box[2]]
    return data, atr


def read_bool(fname, box=None):
    """Read binary file with flags, 1-byte values with flags set in bits
    For ROI_PAC .flg, *_snap_connect.byt file.
    """
    # Read attribute
    if fname.endswith('_snap_connect.byt'):
        rscFile = fname.split('_snap_connect.byt')[0]+'.unw.rsc'
    else:
        rscFile = fname+'.rsc'
    atr = read_attribute(rscFile.split('.rsc')[0])
    if 'DATA_TYPE' not in atr.keys():
        atr['DATA_TYPE'] = 'bool'
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    if not box:
        box = [0, 0, width, length]

    data = np.fromfile(fname,
                       dtype=np.bool_,
                       count=box[3]*width).reshape(box[3], width)
    data = data[box[1]:box[3],
                box[0]:box[2]]
    return data, atr


def read_GPS_USGS(fname):
    yyyymmdd = np.loadtxt(fname,
                          dtype=bytes,
                          usecols=(0, 1)).astype(str)[:, 0]
    YYYYMMDD = []
    for y in yyyymmdd:
        YYYYMMDD.append(y)
    data = np.loadtxt(fname, usecols=(1, 2, 3, 4))
    dates = data[:, 0]
    north = np.array(data[:, 1])
    east = np.array(data[:, 2])
    up = np.array(data[:, 3])

    return east, north, up, dates, YYYYMMDD


#########################################################################
