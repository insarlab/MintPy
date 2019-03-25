############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi, 2013              #
############################################################
# Recommend import:
#   from pysar.utils import readfile


from datetime import datetime as dt
import os
import re
import sys
import warnings
import xml.etree.ElementTree as ET

import h5py
import json
import numpy as np

from pysar.objects import (datasetUnitDict,
                           geometry,
                           geometryDatasetNames,
                           giantIfgramStack,
                           giantTimeseries,
                           ifgramDatasetNames,
                           ifgramStack,
                           timeseriesDatasetNames,
                           timeseriesKeyNames,
                           timeseries,
                           HDFEOS)

standardMetadataKeys = {
    # ordered in alphabet by value names
    'swathNumber':'beam_swath',
    'firstFrameNumber':'first_frame',
    'lastFrameNumber':'last_frame',
    'trackNumber':'relative_orbit',
    'bands': 'number_bands',
    'interleave': 'scheme',

    'altitude': 'HEIGHT',
    'azimuth_looks': 'ALOOKS',
    'azimuthPixelSize': 'AZIMUTH_PIXEL_SIZE',
    'azimuth_pixel_spacing': 'AZIMUTH_PIXEL_SIZE', 'az_pixel_spacing': 'AZIMUTH_PIXEL_SIZE',
    'center_time': 'CENTER_LINE_UTC',
    'corner_lon': 'X_FIRST', 'post_lon': 'X_STEP',
    'corner_lat': 'Y_FIRST', 'post_lat': 'Y_STEP',
    'dataType': 'DATA_TYPE', 'data_type': 'DATA_TYPE',
    'drop_ifgram': 'DROP_IFGRAM',
    'earthRadius': 'EARTH_RADIUS', 'earth_radius_below_sensor': 'EARTH_RADIUS',
    'HEADING_DEG': 'HEADING',
    'length': 'LENGTH', 'FILE_LENGTH': 'LENGTH', 'lines': 'LENGTH',
    'passDirection':'ORBIT_DIRECTION',
    'polarization':'POLARIZATION',
    'prf': 'PRF',
    'rangePixelSize': 'RANGE_PIXEL_SIZE',
    'range_pixel_spacing': 'RANGE_PIXEL_SIZE', 'rg_pixel_spacing': 'RANGE_PIXEL_SIZE',
    'range_looks': 'RLOOKS',
    'ref_date': 'REF_DATE',
    'ref_x': 'REF_X', 'ref_y': 'REF_Y', 'ref_lat': 'REF_LAT', 'ref_lon': 'REF_LON',
    'spacecraftName':'PLATFORM',
    'startingRange': 'STARTING_RANGE', 'near_range_slc': 'STARTING_RANGE',
    'subset_x0': 'SUBSET_XMIN', 'subset_x1': 'SUBSET_XMAX',
    'subset_y0': 'SUBSET_YMIN', 'subset_y1': 'SUBSET_YMAX',
    'wavelength': 'WAVELENGTH', 'Wavelength': 'WAVELENGTH', 'radarWavelength': 'WAVELENGTH',
    'width': 'WIDTH', 'Width': 'WIDTH', 'samples': 'WIDTH',
}


GDAL2NUMPY_DATATYPE = {
    '1': 'uint8',
    '2': 'uint16',
    '3': 'int16',
    '4': 'uint32',
    '5': 'int32',
    '6': 'float32',
    '7': 'float64',
    '10': 'complex64',
    '11': 'complex128',
}

# reference: https://subversion.renater.fr/efidir/trunk/efidir_soft/doc/Programming_C_EFIDIR/header_envi.pdf
ENVI2NUMPY_DATATYPE = {
    '1': 'uint8',
    '2': 'int16',
    '3': 'int32',
    '4': 'float32',
    '5': 'float64',
    '6': 'complex64',
    '9': 'complex128',
    '12': 'uint16',
    '13': 'uint32',
    '14': 'int64',
    '15': 'uint64',
}

ENVI_BAND_INTERLEAVE = {
    'BAND': 'BSQ',
    'LINE': 'BIL',
    'PIXEL': 'BIP',
}


###########################################################################
# obsolete variables
multi_group_hdf5_file = ['interferograms',
                         'coherence',
                         'wrapped',
                         'snaphu_connect_component']
multi_dataset_hdf5_file = ['timeseries', 'geometry']
single_dataset_hdf5_file = ['dem', 'mask', 'temporal_coherence', 'velocity']


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
## for PySAR v1.x files:
##     unwrapPhase-20150115_20150127
##     unwrapPhase-20150115_20150208
##     timeseries-20150115
##     timeseries-20150127
##     temporalCoherence
##
## for PySAR v0.x files: (all in 2D dataset)
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
def read(fname, box=None, datasetName=None, print_msg=True):
    """Read one dataset and its attributes from input file.
    Parameters: fname : str, path of file to read
                datasetName : str or list of str, slice names
                box : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
    Returns:    data : 2/3-D matrix in numpy.array format, return None if failed
                atr : dictionary, attributes of data, return None if failed
    Examples:
        from pysar.utils import readfile
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
        data = read_hdf5_file(fname, datasetName=datasetName, box=box)
    else:
        data = read_binary_file(fname, datasetName=datasetName, box=box)[0]
    return data, atr


#########################################################################
def read_hdf5_file(fname, datasetName=None, box=None):
    """
    Parameters: fname : str, name of HDF5 file to read
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
                box : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
    Returns:    data : 2D/3D array
                atr : dict, metadata
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
    if all(i.isdigit() for i in datasetName):
        datasetName = ['{}-{}'.format(ds_3d_list[0], i) for i in datasetName]
    # Input Argument: decompose slice list into dsFamily and inputDateList
    dsFamily = datasetName[0].split('-')[0]
    inputDateList = [i.replace(dsFamily,'').replace('-','') for i in datasetName]

    # read hdf5
    with h5py.File(fname, 'r') as f:
        # get dataset object
        dsNames = [i for i in [datasetName[0], dsFamily] if i in f.keys()]
        dsNamesOld = [i for i in slice_list if '/{}'.format(datasetName[0]) in i] # support for old pysar files
        if len(dsNames) > 0:
            ds = f[dsNames[0]]
        elif len(dsNamesOld) > 0:
            ds = f[dsNamesOld[0]]
        else:
            raise ValueError('input dataset {} not found in file {}'.format(datasetName, fname))

        # 2D dataset
        if ds.ndim == 2:
            data = ds[box[1]:box[3], box[0]:box[2]]

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
            data = ds[slice_flag, box[1]:box[3], box[0]:box[2]]
            data = np.squeeze(data)
    return data


def read_binary_file(fname, datasetName=None, box=None):
    """Read data from binary file, such as .unw, .cor, etc.
    Parameters: fname : str, path/name of binary file
                datasetName : str, dataset name for file with multiple bands of data
                    e.g.: incidenceAngle, azimuthAngle, rangeCoord, azimuthCoord, ...
                box  : 4-tuple of int area to read, defined in (x0, y0, x1, y1) in pixel coordinate
    Returns:    data : 2D array in size of (length, width) in BYTE / int16 / float32 / complex64 / float64 etc.
                atr  : dict, metadata of binary file
    """
    # Basic Info
    fbase, fext = os.path.splitext(os.path.basename(fname))
    fext = fext.lower()

    # metadata
    atr = read_attribute(fname)
    processor = atr['PROCESSOR']
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    if not box:
        box = (0, 0, width, length)

    # ISCE
    if processor in ['isce']:
        # default short name for data type from ISCE
        dataTypeDict = {
            'byte': 'bool_',
            'float': 'float32',
            'double': 'float64',
            'cfloat': 'complex64',
        }

        # data structure - auto
        data_type = atr['DATA_TYPE'].lower()
        if data_type in dataTypeDict.keys():
            data_type = dataTypeDict[data_type]
        num_band = int(atr['number_bands'])
        band_interleave = atr['scheme'].upper()
        byte_order = 'l'

        # data structure - file specific based on FILE_TYPE - k
        band = 1
        cpx_band = 'phase'

        k = atr['FILE_TYPE'].lower().replace('.', '')
        if k in ['unw']:
            band = 2

        elif k in ['slc']:
            cpx_band = 'magnitude'

        elif k in ['los'] and datasetName and datasetName.startswith(('az', 'head')):
            band = 2

        elif k in ['incLocal']:
            band = 2
            if datasetName and 'local' not in datasetName.lower():
                band = 1

        elif datasetName:
            if datasetName.lower() == 'band2':
                band = 2
            elif datasetName.lower() == 'band3':
                band = 3

    # ROI_PAC
    elif processor in ['roipac']:
        # data structure - auto
        band_interleave = 'BIL'
        byte_order = 'l'

        # data structure - file specific based on file extension
        data_type = 'float32'
        num_band = 1
        band = 1
        cpx_band = 'phase'

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
        else:
            raise Exception('unrecognized ROI_PAC file: {}'.format(fname))

    # Gamma
    elif processor == 'gamma':
        # data structure - auto
        band_interleave = 'BIL'
        byte_order = 'big-endian'

        # data structure - file specific based on file extension
        data_type = 'float32'
        num_band = 1
        band = 1
        cpx_band = 'phase'

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

        else:
            raise Exception('unecognized GAMMA file: {}'.format(fname))

    # reading
    data, atr = read_binary(fname, box=box, data_type=data_type, byte_order=byte_order,
                            num_band=num_band, band_interleave=band_interleave,
                            band=band, cpx_band=cpx_band)
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
                    else:
                        warnings.warn('file has un-defined {}D dataset: {}'.format(obj.ndim, name))
            slice_list = []
            with h5py.File(fname, 'r') as f:
                f.visititems(get_hdf5_2d_dataset)

    # Binary Files
    else:
        if fext.lower() in ['.trans', '.utm_to_rdc']:
            slice_list = ['rangeCoord', 'azimuthCoord']
        elif fbase.startswith('los'):
            slice_list = ['incidenceAngle', 'azimuthAngle']
        elif atr.get('number_bands', '1') == '2' and 'unw' not in k:
            slice_list = ['band1', 'band2']
        else:
            slice_list = ['']
    return slice_list


def get_dataset_list(fname, datasetName=None):
    """Get list of 2D and 3D dataset to facilitate systematic file reading"""
    if datasetName:
        return [datasetName]

    fbase, fext = os.path.splitext(os.path.basename(fname))
    fext = fext.lower()

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

    elif fext in ['.trans', '.utm_to_rdc']:
        ds_list = ['rangeCoord', 'azimuthCoord']
    elif fbase.startswith('los'):
        ds_list = ['incidenceAngle', 'azimuthAngle']
    else:
        ds_list = [os.path.splitext(fbase)[0]]
    return ds_list


#########################################################################
def read_attribute(fname, datasetName=None, standardize=True, metafile_ext=None):
    """Read attributes of input file into a dictionary
    Parameters: fname : str, path/name of data file
                datasetName : str, name of dataset of interest, for file with multiple datasets
                    e.g. unwrapPhase in ifgramStack.h5
                         coherence   in ifgramStack.h5
                         height      in geometryRadar.h5
                         latitude    in geometryRadar.h5
                         ...
                standardize : bool, grab standardized metadata key name
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
        f = h5py.File(fname, 'r')
        g1_list = [i for i in f.keys() if isinstance(f[i], h5py.Group)]
        d1_list = [i for i in f.keys() if isinstance(f[i], h5py.Dataset) and f[i].ndim >= 2]

        # FILE_TYPE - k
        if any(i in d1_list for i in ['unwrapPhase']):
            k = 'ifgramStack'
        elif any(i in d1_list for i in ['height', 'latitude', 'azimuthCoord']):
            k = 'geometry'
        elif any(i in g1_list+d1_list for i in ['timeseries', 'displacement']):
            k = 'timeseries'
        elif 'HDFEOS' in g1_list:
            k = 'HDFEOS'
        elif 'recons' in d1_list:
            k = 'giantTimeseries'
        elif any(i in d1_list for i in ['igram', 'figram']):
            k = 'giantIfgramStack'
        elif any(i in g1_list for i in multi_group_hdf5_file):      # old pysar format
            k = list(set(g1_list) & set(multi_group_hdf5_file))[0]
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
        else:
            if len(f.attrs) > 0 and 'WIDTH' in f.attrs.keys():
                atr = dict(f.attrs)
            else:
                # grab the list of attrs in HDF5 file
                global atr_list
                def get_hdf5_attrs(name, obj):
                    global atr_list
                    if len(obj.attrs) > 0 and 'WIDTH' in obj.attrs.keys():
                        atr_list.append(dict(obj.attrs))
                atr_list = []
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

        # attribute identified by PySAR
        # 1. FILE_TYPE
        atr['FILE_TYPE'] = str(k)

        # 2. DATA_TYPE
        ds = None
        if datasetName and datasetName in f.keys():
            ds = f[datasetName]
        else:
            # get the 1st dataset
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
        f.close()

        # 3. PROCESSOR
        if 'INSAR_PROCESSOR' in atr.keys():
            atr['PROCESSOR'] = atr['INSAR_PROCESSOR']
        if 'PROCESSOR' not in atr.keys():
            atr['PROCESSOR'] = 'pysar'

    else:
        # get existing metadata file extensions
        metafile_exts = ['.rsc', '.xml', '.aux.xml', '.par', '.hdr']
        if metafile_ext:
            metafile_exts = [i for i in metafile_exts if i.endswith(metafile_ext)]
        metafile_exts = [i for i in metafile_exts if os.path.isfile(fname+i)]     
        if len(metafile_exts) == 0:
            raise FileNotFoundError('No metadata file found for data file: {}'.format(fname))

        atr = {}
        # PROCESSOR
        if any(i.endswith(('.xml', '.hdr')) for i in metafile_exts):
            atr['PROCESSOR'] = 'isce'
            xml_exts = [i for i in metafile_exts if i.endswith('.xml')]
            if len(xml_exts) > 0:
                atr.update(read_isce_xml(fname+xml_exts[0]))
        elif any(i.endswith('.par') for i in metafile_exts):
            atr['PROCESSOR'] = 'gamma'
        elif any(i.endswith('.rsc') for i in metafile_exts):
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'roipac'
        if 'PROCESSOR' not in atr.keys():
            atr['PROCESSOR'] = 'pysar'

        # Read metadata file and FILE_TYPE
        metafile0 = fname + metafile_exts[0]
        while fext in ['.geo', '.rdr']:
            fbase, fext = os.path.splitext(fbase)
        if not fext:
            fext = fbase
        if metafile0.endswith('.rsc'):
            atr.update(read_roipac_rsc(metafile0))
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = fext

        elif metafile0.endswith('.xml'):
            atr.update(read_isce_xml(metafile0))
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = fext  #atr.get('image_type', fext)

        elif metafile0.endswith('.par'):
            atr.update(read_gamma_par(metafile0))
            atr['FILE_TYPE'] = fext

        elif metafile0.endswith('.hdr'):
            atr.update(read_envi_hdr(metafile0))
            atr['FILE_TYPE'] = atr['file type']

    # UNIT
    k = atr['FILE_TYPE'].replace('.', '')
    if k == 'ifgramStack':
        if datasetName and datasetName in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[datasetName]
        else:
            atr['UNIT'] = 'radian'
    elif 'UNIT' not in atr.keys():
        if datasetName and datasetName in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[datasetName]
        elif k in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[k]
        else:
            atr['UNIT'] = '1'

    # FILE_PATH
    atr['FILE_PATH'] = os.path.abspath(fname)

    if standardize:
        atr = standardize_metadata(atr)
    return atr


def standardize_metadata(metaDict, standardKeys=standardMetadataKeys):
    metaDict_out = {}
    for k in metaDict.keys():
        metaDict_out[k] = metaDict[k]
        if k in standardKeys.keys():
            k2 = standardKeys[k]
            if k2 in metaDict.keys():
                metaDict_out[k2] = metaDict[k2]
            else:
                metaDict_out[k2] = metaDict[k]
    return metaDict_out


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
        from pysar.defaults.auto_path import isceAutoPath
        tmpl = read_template(isceAutoPath, print_msg=False)
    """
    template_dict = {}
    plotAttributeDict = {}
    insidePlotObject = False
    plotAttributes = []
    # the below logic for plotattributes object can be made much more simple
    # if we assume that any plot attribute coming after a > belongs to the
    # same object. Must Ask Falk and Yunjung if we can assume this to eliminate
    # all these conditionals

    if os.path.isfile(fname):
        f = open(fname, 'r')
        lines = f.readlines()
    elif isinstance(fname, str):
        lines = fname.split('\n')

    for line in lines:
        line = line.strip()
        # split on the 1st occurrence of delimiter
        c = [i.strip() for i in line.split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%', '#')):
            if line.startswith(">"):
                plotAttributeDict = {}
                insidePlotObject = True
            # otherwise, if previously inside attributes object, we are now outside
            # unless the line is a comment
            elif insidePlotObject and not line.startswith('%') and not line.startswith('#'):
                # just came from being inside plot object, but now we are outside
                insidePlotObject = False
                plotAttributes.append(plotAttributeDict)
            next  # ignore commented lines or those without variables
        else:
            atrName = c[0]
            atrValue = str.replace(c[1], '\n', '').split("#")[0].strip()
            atrValue = os.path.expanduser(atrValue)
            atrValue = os.path.expandvars(atrValue)

            if insidePlotObject:
                if is_plot_attribute(atrName):
                    plotAttributeDict[atrName] = atrValue
                else:
                    # just came from being inside plot object, but now we are outside
                    insidePlotObject = False
                    plotAttributes.append(plotAttributeDict)
                    template_dict[atrName] = atrValue

            elif atrValue != '':
                template_dict[atrName] = atrValue
    if os.path.isfile(fname):
        f.close()

    # what if no \n at end of file? write out last plot attributes dict
    if insidePlotObject:
        plotAttributes.append(plotAttributeDict)

    if len(plotAttributes) > 0:
        template_dict["plotAttributes"] = json.dumps(plotAttributes)

    return template_dict


def is_plot_attribute(attribute):
    tokens = attribute.split(".")
    if tokens is None:
        return False
    return tokens[0] == "plot" and len(tokens) > 1


def read_roipac_rsc(fname, delimiter=' ', standardize=True):
    """Read ROI_PAC style RSC file.
    Parameters: fname : str.
                    File path of .rsc file.
    Returns:    rscDict : dict
                    Dictionary of keys and values in RSC file.
    Examples:
        from pysar.utils import readfile
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

    if standardize:
        rscDict = standardize_metadata(rscDict)
    return rscDict


def read_gamma_par(fname, delimiter=':', skiprows=3, standardize=True):
    """Read GAMMA .par/.off file into a python dictionary structure.
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
    if standardize:
        parDict = standardize_metadata(parDict)
    return parDict


def read_isce_xml(fname, standardize=True):
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
            if ET.iselement(child):
                v_step  = float(child.find("./property[@name='delta']").find('value').text)
                v_first = float(child.find("./property[@name='startingvalue']").find('value').text)
                if abs(v_step) < 1. and abs(v_step) > 1e-7:
                    xmlDict['{}_STEP'.format(prefix)] = v_step
                    xmlDict['{}_FIRST'.format(prefix)] = v_first - v_step / 2.

    # PAMDataset, e.g. hgt.rdr.aux.xml
    elif root.tag == 'PAMDataset':
        meta = root.find("./Metadata[@domain='ENVI']")
        for child in meta.findall("MDI"):
            key = child.get('key')
            value = child.text
            xmlDict[key] = value
        xmlDict['data_type'] = ENVI2NUMPY_DATATYPE[xmlDict['data_type']]

    if standardize:
        xmlDict = standardize_metadata(xmlDict)
    return xmlDict


def read_envi_hdr(fname, standardize=True):
    """Read ENVI .hdr file into a python dict strcture"""
    atr = read_template(fname)
    atr['DATA_TYPE'] = ENVI2NUMPY_DATATYPE[atr.get('data type', '4')]
    if 'map info' in atr.keys():
        map_info = [i.strip() for i in atr['map info'].split(',')]
        x_step = abs(float(map_info[5]))
        y_step = abs(float(map_info[6])) * -1.
        if abs(x_step) < 1. and abs(x_step) > 1e-7:
            atr['X_FIRST'] = str(float(map_info[3]) - x_step / 2.)
            atr['Y_FIRST'] = str(float(map_info[4]) - y_step / 2.)
            atr['X_STEP'] = str(x_step)
            atr['Y_STEP'] = str(y_step)
    if standardize:
        atr = standardize_metadata(atr)
    return atr


def attribute_gamma2roipac(par_dict_in):
    """Convert Gamma par attribute into ROI_PAC format"""
    par_dict = dict()
    for key, value in iter(par_dict_in.items()):
        par_dict[key] = value

    # Length - number of rows
    for key in par_dict_in.keys():
        if any(key.startswith(i) for i in ['azimuth_lines',
                                           'nlines',
                                           'az_samp',
                                           'interferogram_azimuth_lines']):
            par_dict['LENGTH'] = par_dict[key]

    # Width - number of columns
    for key in par_dict_in.keys():
        if any(key.startswith(i) for i in ['width',
                                           'range_samp',
                                           'interferogram_width']):
            par_dict['WIDTH'] = par_dict[key]

    # WAVELENGTH
    speed_of_light = 299792458.0   # meter/second
    key = 'radar_frequency'
    if key in par_dict_in.keys():
        par_dict['WAVELENGTH'] = str(speed_of_light/float(par_dict[key]))

    # HEIGHT & EARTH_RADIUS
    key = 'earth_radius_below_sensor'
    if key in par_dict_in.keys():
        par_dict['EARTH_RADIUS'] = par_dict[key]

        key2 = 'sar_to_earth_center'
        if key2 in par_dict_in.keys():
            par_dict['HEIGHT'] = str(float(par_dict[key2]) - float(par_dict[key]))

    # PLATFORM
    key = 'sensor'
    if key in par_dict_in.keys():
        par_dict['PLATFORM'] = par_dict[key]

    # ORBIT_DIRECTION
    key = 'heading'
    if key in par_dict_in.keys():
        value = float(par_dict[key])
        if 270 < value < 360 or -90 < value < 90:
            par_dict['ORBIT_DIRECTION'] = 'ascending'
        else:
            par_dict['ORBIT_DIRECTION'] = 'descending'
        par_dict['HEADING'] = str(value)

    # Optional attributes for PySAR from ROI_PAC
    # ANTENNA_SIDE
    key = 'azimuth_angle'
    if key in par_dict_in.keys():
        value = float(par_dict[key])
        if 0 < value < 180:
            par_dict['ANTENNA_SIDE'] = '-1'
        else:
            par_dict['ANTENNA_SIDE'] = '1'

    return par_dict


#########################################################################
def read_binary(fname, box=None, data_type='float32', byte_order='l',
                num_band=1, band_interleave='BIL', band=1, cpx_band='phase'):
    """Read binary file using np.fromfile
    Parameters: fname : str, path/name of data file to read
                box   : tuple of 4 int in (x0, y0, x1, y1)
                data_type : str, data type of stored array, e.g.:
                    bool_
                    int8, int16, int32
                    float16, float32, float64
                    complex64, complex128
                byte_order : str, little/big-endian
                num_band   : int, number of bands
                band_interleave : str, band interleaving scheme, e.g.:
                    BIP
                    BIL
                    BSQ
                band : int, band of interest, between 1 and num_band.
                cpx_band : str, e.g.:
                    real,
                    imag, imaginary
                    phase,
                    mag, magnitude
    Returns:    data : 2D np.array
                atr : dict, metadata
    Examples:   # ISCE files
                data, atr = read_binary('filt_fine.unw', num_band=2, band=2)
                data, atr = read_binary('filt_fine.cor')
                data, atr = read_binary('filt_fine.int', data_type='complex64', cpx_band='phase')
                data, atr = read_binary('burst_01.slc',  data_type='complex64', cpx_band='mag')
                data, atr = read_binary('los.rdr', num_band=2, band=1)
                # ROIPAC files                
    """
    atr = read_attribute(fname)
    length, width = int(float(atr['LENGTH'])), int(float(atr['WIDTH']))
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
    if num_band > 1:
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
    else:
        data = np.fromfile(fname,
                           dtype=data_type,
                           count=box[3]*width).reshape(-1, width)
        data = data[box[1]:box[3],
                    box[0]:box[2]]

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
        else:
            raise ValueError('unrecognized complex band:', cpx_band)

    return data, atr


############################ Obsolete Functions ###############################
def read_float32(fname, box=None, byte_order='l'):
    """Reads roi_pac data (RMG format, interleaved line by line)
    should rename it to read_rmg_float32()

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
