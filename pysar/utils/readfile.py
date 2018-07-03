############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi, 2013              #
############################################################
# Recommend import:
#   from pysar.utils import readfile


import os
import sys
import re
from datetime import datetime as dt
import h5py
import numpy as np
#from PIL import Image
import json

from pysar.objects import (datasetUnitDict,
                           geometry,
                           geometryDatasetNames,
                           giantTimeseries,
                           ifgramDatasetNames,
                           ifgramStack,
                           timeseriesDatasetNames,
                           timeseriesKeyNames,
                           timeseries,
                           HDFEOS)


standardMetadataKeys = {'width': 'WIDTH', 'Width': 'WIDTH', 'samples': 'WIDTH',
                        'length': 'LENGTH', 'FILE_LENGTH': 'LENGTH', 'lines': 'LENGTH',
                        'wavelength': 'WAVELENGTH', 'Wavelength': 'WAVELENGTH', 'radarWavelength': 'WAVELENGTH',
                        'prf': 'PRF',
                        'post_lat': 'Y_STEP',
                        'post_lon': 'X_STEP',
                        'range_looks': 'RLOOKS',
                        'azimuth_looks': 'ALOOKS',
                        'dataType': 'DATA_TYPE', 'data_type': 'DATA_TYPE',
                        'rangePixelSize': 'RANGE_PIXEL_SIZE',
                        'range_pixel_spacing': 'RANGE_PIXEL_SIZE', 'rg_pixel_spacing': 'RANGE_PIXEL_SIZE',
                        'azimuthPixelSize': 'AZIMUTH_PIXEL_SIZE',
                        'azimuth_pixel_spacing': 'AZIMUTH_PIXEL_SIZE', 'az_pixel_spacing': 'AZIMUTH_PIXEL_SIZE',
                        'earthRadius': 'EARTH_RADIUS', 'earth_radius_below_sensor': 'EARTH_RADIUS',
                        'altitude': 'HEIGHT',
                        'startingRange': 'STARTING_RANGE',
                        'center_time': 'CENTER_LINE_UTC',
                        'drop_ifgram': 'DROP_IFGRAM',
                        'ref_date': 'REF_DATE',
                        'ref_x': 'REF_X', 'ref_y': 'REF_Y', 'ref_lat': 'REF_LAT', 'ref_lon': 'REF_LON',
                        'subset_x0': 'SUBSET_XMIN', 'subset_x1': 'SUBSET_XMAX',
                        'subset_y0': 'SUBSET_YMIN', 'subset_y1': 'SUBSET_YMAX',
                        }

GDAL2NUMPY_DATATYPE = {

    1: np.uint8,
    2: np.uint16,
    3: np.int16,
    4: np.uint32,
    5: np.int32,
    6: np.float32,
    7: np.float64,
    10: np.complex64,
    11: np.complex128,

}

#########################################################################
# obsolete variables
multi_group_hdf5_file = ['interferograms',
                         'coherence',
                         'wrapped',
                         'snaphu_connect_component']
multi_dataset_hdf5_file = ['timeseries', 'geometry']
single_dataset_hdf5_file = ['dem', 'mask', 'temporal_coherence', 'velocity']


#########################################################################
def read(fname, box=None, datasetName=None, print_msg=True):
    """Read one dataset and its attributes from input file.
    Parameters: fname : str, path of file to read
                    PySAR   file: interferograms, timeseries, velocity, etc.
                    ROI_PAC file: .unw .cor .hgt .dem .trans
                    Gamma   file: .mli .slc
                    Image   file: .jpeg .jpg .png .ras .bmp
                box : 4-tuple of int
                    area to read, defined in (x0, y0, x1, y1) in pixel coordinate
                datasetName : string
                    dataset to read in the format of
                        datasetName
                        datasetName-dateName
                        datasetName-date12Name

                    for ifgramStack:
                        unwrapPhase
                        coherence
                        ...
                        unwrapPhase-20161020_20161026
                        ...
                    for timeseries:
                        timeseries
                        timeseries-20161020
                        timeseries-20161026
                        ...
                    for geometry:
                        height
                        incidenceAngle
                        bperp
                        ...
                        bperp-20161020
                        bperp-20161026
                        bperp-...
                    for .trans and .utm_to_rdc:
                        rangeCoord
                        azimuthCoord
                    for los.rdr:
                        incidenceAngle
                        headingAngle
                    for GIAnT:
                        giantTimeseries
                        giantTimeseries-20161020
                        giantTimeseries-20161026
                        giantTimeseries-...
                    for the other single dataset file:
                        No need and will be ignored.

    Returns: data : 2/3-D matrix in numpy.array format, return None if failed
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

    # Basic Info
    ext = os.path.splitext(fname)[1].lower()
    fbase = os.path.splitext(os.path.basename(fname))[0]
    if isinstance(datasetName, list):
        atr = read_attribute(fname, datasetName=datasetName[0].split('-')[0])
    elif isinstance(datasetName, str):
        atr = read_attribute(fname, datasetName=datasetName.split('-')[0])
    else:
        atr = read_attribute(fname)
    k = atr['FILE_TYPE']
    processor = atr['PROCESSOR']
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    if not box:
        box = (0, 0, width, length)

    # HDF5
    if ext in ['.h5', '.he5']:
        f = h5py.File(fname, 'r')
        if k in ['timeseries']:
            if isinstance(f[k], h5py.Dataset):
                obj = timeseries(fname)
                data = obj.read(datasetName=datasetName,
                                box=box,
                                print_msg=print_msg)
            else:
                # support for old pysar format
                if not datasetName:
                    datasetName = list(f[k].keys())
                if isinstance(datasetName, str):
                    datasetName = [datasetName]
                data = np.zeros((len(datasetName),
                                 box[3]-box[1],
                                 box[2]-box[0]), np.float32)
                for i in range(len(datasetName)):
                    date = datasetName[i].split('-')[1]
                    data[i, :, :] = f[k][date][box[1]:box[3],
                                               box[0]:box[2]]
                data = np.squeeze(data)

        elif k in ['giantTimeseries']:
            obj = giantTimeseries(fname)
            obj.open(print_msg=False)
            dsName = [i for i in ['recons', 'rawts'] if i in f.keys()][0]
            if print_msg:
                print('reading {} from file {}'.format(dsName, fbase+ext))
            if not datasetName:
                datasetName = obj.datasetList
            if isinstance(datasetName, str):
                datasetName = [datasetName]
            data = np.zeros((len(datasetName),
                                 box[3]-box[1],
                                 box[2]-box[0]), np.float32)
            for i in range(len(datasetName)):
                date = datasetName[i].split('-')[1]
                idx = obj.dateList.index(date)
                data[i, :, :] = f[dsName][idx,
                                          box[1]:box[3],
                                          box[0]:box[2]]
            data = np.squeeze(data)

        elif k in ['ifgramStack']:
            obj = ifgramStack(fname)
            data = obj.read(datasetName=datasetName,
                            box=box,
                            print_msg=print_msg)
            if datasetName in ['unwrapPhase', 'wrapPhase', 'iono']:
                atr['UNIT'] = 'radian'
            else:
                atr['UNIT'] = '1'

        elif k in ['geometry']:
            obj = geometry(fname)
            data = obj.read(datasetName=datasetName,
                            box=box,
                            print_msg=print_msg)

        elif k == 'HDFEOS':
            obj = HDFEOS(fname)
            data = obj.read(datasetName=datasetName,
                            box=box,
                            print_msg=print_msg)

        # old pysar format
        elif k in multi_group_hdf5_file and isinstance(f[k], h5py.Group):
            if not datasetName:
                datasetName = list(f[k].keys())
            if isinstance(datasetName, str):
                datasetName = [datasetName]
            data = np.zeros((len(datasetName),
                                 box[3]-box[1],
                                 box[2]-box[0]), np.float32)
            for i in range(len(datasetName)):
                dsName = datasetName[i]
                data[i, :, :] = f[k][dsName][dsName][box[1]:box[3],
                                                     box[0]:box[2]]
            data = np.squeeze(data)

        else:
            k0 = list(f.keys())[0]
            if isinstance(f[k0], h5py.Dataset):
                if datasetName and datasetName in f.keys():
                    dset = f[datasetName]
                else:
                    dset = f[k0]
            else:
                # support for old pysar format
                k1 = list(f[k0].keys())[0]
                if isinstance(f[k0][k1], h5py.Dataset):
                    if datasetName:
                        dset = f[k0][datasetName]
                    else:
                        dset = f[k0][k1]

            data = dset[box[1]:box[3], box[0]:box[2]]
            atr['LENGTH'] = str(dset.shape[0])
            atr['WIDTH'] = str(dset.shape[1])
        f.close()
        return data, atr

    # Image
    # elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
    #    atr = read_roipac_rsc(fname+'.rsc')
    #    data  = Image.open(fname).crop(box)
    #    return data, atr

    # ISCE
    elif processor in ['isce']:
        if k in ['.unw', 'unw']:
            amp, pha, atr = read_float32(fname, box=box)
            return pha, atr

        elif k in ['.cor', 'cor']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif k in ['.int', 'int', '.flat', 'cpx']:
            data, atr = read_complex_float32(fname, box=box, band='phase')
            return data, atr

        elif k in ['.slc']:
            data, atr = read_complex_float32(fname, box=box, band='amplitude')
            return data, atr

        elif fbase.startswith('los'):
            incAngle, azAngle, atr = read_float32(fname, box=box)
            if not datasetName:
                return incAngle, azAngle, atr
            elif datasetName.startswith('inc'):
                return incAngle, atr
            elif datasetName.startswith(('az', 'head')):
                return azAngle, atr
            else:
                raise Exception('Un-recognized datasetName input: '+datasetName)

        elif atr['DATA_TYPE'].lower() in ['float64', 'double']:
            data, atr = read_real_float64(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['cfloat32']:
            data, atr = read_complex_float32(fname, box=box, band='complex')
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['float32', 'float']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['int16', 'short']:
            data, atr = read_real_int16(fname, box=box)
            return data, atr
        elif atr['DATA_TYPE'].lower() in ['bool', 'byte', 'flag']:
            data, atr = read_bool(fname, box=box)
            return data, atr
        else:
            raise Exception('Un-recognized {} file: {}'.format(processor,
                                                               os.path.basename(fname)))

    # ROI_PAC
    elif processor in ['roipac']:
        if ext in ['.unw', '.cor', '.hgt', '.msk']:
            amp, pha, atr = read_float32(fname, box=box)
            return pha, atr

        elif ext in ['.dem', '.wgs84']:
            dem, atr = read_real_int16(fname, box=box)
            return dem, atr

        elif ext in ['.int']:
            data, atr = read_complex_float32(fname, box=box, band='phase')
            return data, atr

        elif ext in ['.amp']:
            data, atr = read_complex_float32(fname, box=box, band='complex')
            return data.real, data.imag, atr

        elif ext in ['.flt']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif ext in ['.flg', '.byt']:
            flag, atr = read_bool(fname, box=box)
            return flag, atr

        elif ext in ['.trans']:
            rg, az, atr = read_float32(fname, box=box)
            if not datasetName:
                return rg, az, atr
            elif datasetName.startswith(('rg', 'range')):
                return rg, atr
            elif datasetName.startswith(('az', 'azimuth')):
                return az, atr
            else:
                raise Exception('Un-recognized datasetName input: '+datasetName)
        else:
            raise Exception('Un-recognized {} file: {}'.format(processor,
                                                               os.path.basename(fname)))

    # Gamma
    elif processor == 'gamma':
        if ext in ['.unw', '.cor', '.hgt_sim', '.dem']:
            data, atr = read_real_float32(fname, box=box, byte_order='ieee-be')
            return data, atr

        elif ext in ['.UTM_TO_RDC', '.utm_to_rdc']:
            data, atr = read_complex_float32(fname, box=box, byte_order='ieee-be',
                                             band='complex')
            if not datasetName:
                return data.real, data.imag, atr
            elif datasetName.startswith(('rg', 'range')):
                return data.real, atr
            elif datasetName.startswith(('az', 'azimuth')):
                return data.imag, atr
            else:
                raise Exception('Un-recognized datasetName input: '+datasetName)

        elif ext in ['.int']:
            data, atr = read_complex_float32(fname, box=box, byte_order='ieee-be',
                                             band='phase')
            return data, atr

        elif ext in ['.mli']:
            data, atr = read_real_float32(fname, box=box)
            return data, atr

        elif ext in ['.amp', '.ramp']:
            data, atr = read_real_float32(fname, box=box, byte_order='ieee-be')
            return data, atr

        elif ext == '.slc':
            amp, pha, atr = read_complex_int16(fname, box=box, cpx=False)
            return amp, atr

        else:
            raise Exception('Un-recognized {} file: {}'.format(processor,
                                                               os.path.basename(fname)))
    else:
        raise Exception('Unrecognized file format: '+ext)


#########################################################################
def get_dataset_list(fname, datasetName=None):
    """Get list of 2D and 3D dataset to facilitate systematic file reading"""
    if datasetName:
        return [datasetName]

    atr = read_attribute(fname)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])

    datasetList = []
    fbase = os.path.basename(fname)
    ext = os.path.splitext(fname)[1].lower()
    if ext in ['.h5', '.he5']:
        with h5py.File(fname, 'r') as f:
            k0 = list(f.keys())[0]
            if isinstance(f[k0], h5py.Dataset):
                for key in f.keys():
                    if (isinstance(f[key], h5py.Dataset)
                            and f[key].shape[-2:] == (length, width)):
                        datasetList.append(key)

            # support for old pysar format
            else:
                k1 = list(f[k0].keys())[0]
                if isinstance(f[k0][k1], h5py.Dataset):
                    for key in f[k0].keys():
                        if (isinstance(f[k0][key], h5py.Dataset)
                                and f[k0][key].shape[-2:] == (length, width)):
                            datasetList.append(key)

    elif ext in ['.trans', '.utm_to_rdc']:
        datasetList = ['rangeCoord', 'azimuthCoord']
    elif fbase.startswith('los'):
        datasetList = ['incidenceAngle', 'headingAngle']
    else:
        datasetList = [os.path.split(fbase)[0]]
    return datasetList


def get_2d_dataset_list(fname):
    """Get list of 2D dataset existed in file (for display)"""
    file_ext = os.path.splitext(fname)[1]
    file_base = os.path.splitext(os.path.basename(fname))[0]
    file_type = read_attribute(fname)['FILE_TYPE']
    datasetList = []
    # HDF5 Files
    if file_ext in ['.h5', '.he5']:
        with h5py.File(fname, 'r') as f:
            if file_type in ['timeseries']:
                if isinstance(f[file_type], h5py.Dataset):
                    obj = timeseries(fname)
                    obj.open(print_msg=False)
                    datasetList = obj.datasetList
                else:
                    date_list = list(f[file_type].keys())
                    datasetList = ['{}-{}'.format(file_type, i) for i in date_list]

            elif file_type in ['geometry']:
                obj = geometry(fname)
                obj.open(print_msg=False)
                datasetList = obj.datasetList

            elif file_type in ['ifgramStack']:
                obj = ifgramStack(fname)
                obj.open(print_msg=False)
                datasetList = obj.datasetList

            elif file_type in ['HDFEOS']:
                obj = HDFEOS(fname)
                obj.open(print_msg=False)
                datasetList = obj.datasetList

            elif file_type in ['giantTimeseries']:
                obj = giantTimeseries(fname)
                obj.open(print_msg=False)
                datasetList = obj.datasetList

            else:
                k0 = list(f.keys())[0]
                if isinstance(f[k0], h5py.Dataset):
                    datasetList = sorted(list(f.keys()))

                # support for old pysar format
                else:
                    k1 = list(f[k0].keys())[0]
                    if isinstance(f[k0][k1], h5py.Dataset):
                        datasetList = sorted(list(f[k0].keys()))

                    # unwrapIfgram.h5, coherence.h5, etc.
                    elif isinstance(f[k0][k1], h5py.Group):
                        datasetList = sorted(list(f[k0].keys()))

    # Binary Files
    else:
        if file_ext.lower() in ['.trans', '.utm_to_rdc']:
            datasetList = ['rangeCoord', 'azimuthCoord']
        elif file_base.startswith('los'):
            datasetList = ['incidenceAngle', 'headingAngle']
        else:
            datasetList = ['']
    return datasetList


#########################################################################
def read_attribute(fname, datasetName=None, standardize=True):
    """Read attributes of input file into a dictionary
    Input  : string, file name
    Output : dictionary, attributes dictionary
    """
    ext = os.path.splitext(fname)[1].lower()
    if not os.path.isfile(fname):
        msg = 'input file not existed: {}\n'.format(fname)
        msg += 'current directory: '+os.getcwd()
        raise Exception(msg)

    # HDF5 files
    if ext in ['.h5', '.he5']:
        f = h5py.File(fname, 'r')
        g1_list = [i for i in f.keys() if isinstance(f[i], h5py.Group)]
        d1_list = [i for i in f.keys() if isinstance(f[i], h5py.Dataset) and f[i].ndim >= 2]

        # FILE_TYPE - k
        if 'unwrapPhase' in d1_list:
            k = 'ifgramStack'
        elif any(i in d1_list for i in ['height', 'latitude', 'azimuthCoord']):
            k = 'geometry'
        elif any(i in g1_list+d1_list for i in timeseriesDatasetNames):
            k = 'timeseries'
        elif 'HDFEOS' in g1_list:
            k = 'HDFEOS'
        elif 'recons' in d1_list:
            k = 'giantTimeseries'
        elif any(i in g1_list for i in multi_group_hdf5_file):      # old pysar format
            k = list(set(g1_list) & set(multi_group_hdf5_file))[0]
        elif len(d1_list) > 0:
            k = d1_list[0]
        elif len(g1_list) > 0:
            k = g1_list[0]
        else:
            raise ValueError('unrecognized file type: '+fname)

        # search existing metadata dict
        atr = None
        key = 'WIDTH'
        if key in f.attrs.keys():
            atr = dict(f.attrs)
        else:
            for g in f.keys():
                if key in f[g].attrs.keys():
                    atr = dict(f[g].attrs)
                    break

        if atr is None:
            if k in multi_group_hdf5_file:     # old pysar format
                k2 = list(f[k].keys())[0]
                atr = dict(f[k][k2].attrs)
            elif k == 'giantTimeseries':
                atr = giantTimeseries(fname).get_metadata()

        if atr is None:
            raise ValueError('No attribute {} found in 0/1/2 group level!'.format(key))

        for key, value in atr.items():
            try:
                atr[key] = value.decode('utf8')
            except:
                atr[key] = value

        atr['FILE_TYPE'] = str(k)

        # DATA_TYPE
        k0 = list(f.keys())[0]
        if isinstance(f[k0], h5py.Dataset):
            if datasetName and datasetName in f.keys():
                dset = f[datasetName]
            else:
                dset = f[k0]
        else:
            # support for old pysar format
            k1 = list(f[k0].keys())[0]
            if isinstance(f[k0][k1], h5py.Dataset):
                if datasetName and datasetName in f[k0].keys():
                    dset = f[k0][datasetName]
                else:
                    dset = f[k0][k1]
            else:
                k2 = list(f[k0][k1].keys())[0]
                if isinstance(f[k0][k1][k2], h5py.Dataset):
                    dset = f[k0][k1][k2]
        atr['DATA_TYPE'] = str(dset.dtype)

        # PROCESSOR
        if 'INSAR_PROCESSOR' in atr.keys():
            atr['PROCESSOR'] = atr['INSAR_PROCESSOR']
        f.close()

    else:
        # Read metadata file
        if os.path.isfile(fname+'.rsc'):
            atr = read_roipac_rsc(fname+'.rsc')
            atr['FILE_TYPE'] = ext

        elif os.path.isfile(fname+'.xml'):
            atr = read_isce_xml(fname+'.xml')
            if 'FILE_TYPE' not in atr.keys():
                atr['FILE_TYPE'] = ext

        elif os.path.isfile(fname+'.par'):
            atr = read_gamma_par(fname+'.par')
            atr['FILE_TYPE'] = ext

        elif os.path.isfile(fname+'.hdr'):
            atr = read_template(fname+'.hdr')
            atr = attribute_envi2roipac(atr)
            atr['FILE_TYPE'] = atr['file type']
        else:
            raise Exception('Unrecognized file extension: '+ext)

        # Get PROCESSOR
        if os.path.isfile(fname+'.xml'):
            atr['PROCESSOR'] = 'isce'
        elif os.path.isfile(fname+'.hdr'):
            atr['PROCESSOR'] = 'isce'
        elif os.path.isfile(fname+'.par'):
            atr['PROCESSOR'] = 'gamma'
        elif os.path.isfile(fname+'.rsc'):
            if 'PROCESSOR' not in atr.keys():
                atr['PROCESSOR'] = 'roipac'

    # Unit - str
    k = atr['FILE_TYPE']
    if k == 'ifgramStack':
        if datasetName and datasetName in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[datasetName]
        else:
            atr['UNIT'] = 'radian'
    elif 'UNIT' not in atr.keys():
        if k in datasetUnitDict.keys():
            atr['UNIT'] = datasetUnitDict[k]
        else:
            atr['UNIT'] = '1'

    if 'PROCESSOR' not in atr.keys():
        atr['PROCESSOR'] = 'pysar'

    atr['FILE_PATH'] = os.path.abspath(fname)
    if ext == '.wgs84' and atr['PROCESSOR'] == 'isce':
        atr['FILE_TYPE'] = 'dem'

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
def check_variable_name(path, print_msg=True):
    s = path.split("/")[0]
    if len(s) > 0 and s[0] == "$":
        try:
            p0 = os.getenv(s[1:])
            path = path.replace(path.split("/")[0], p0)
        except:
            if print_msg:
                print('WARNING: Un-recognized environmental variable: '+s)
    return path


def is_plot_attribute(attribute):
    tokens = attribute.split(".")
    if tokens is None:
        return False
    return tokens[0] == "plot" and len(tokens) > 1


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
            atrValue = check_variable_name(atrValue, print_msg=print_msg)

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


def read_roipac_rsc(fname, standardize=True):
    """Read ROI_PAC style RSC file.
    Parameters: fname : str.
                    File path of .rsc file.
    Returns:    rscDict : dict
                    Dictionary of keys and values in RSC file.
    Examples:
        from pysar.utils import readfile
        atr = readfile.read_roipac_rsc('filt_101120_110220_c10.unw.rsc')
    """
    #rsc_dict = dict(np.loadtxt(fname, dtype=bytes, usecols=(0,1)).astype(str))
    # return rsc_dict
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    rscDict = {}
    for line in lines:
        key, value = line.strip().split()[0:2]
        rscDict[key] = value

    if standardize:
        rscDict = standardize_metadata(rscDict)
    return rscDict


def read_gamma_par(fname, delimiter=':', skiprows=3, convert2roipac=True, standardize=True):
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
    parDict = {}

    # Read txt file
    f = open(fname, 'r')
    lines = f.readlines()[skiprows:]
    for line in lines:
        line = line.strip()
        c = [i.strip() for i in line.split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%', '#')):
            next
        else:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].split()[0].strip()
            parDict[key] = value
    f.close()

    if convert2roipac:
        parDict = attribute_gamma2roipac(parDict)

    if standardize:
        parDict = standardize_metadata(parDict)

    return parDict


def read_isce_xml(fname, convert2roipac=True, standardize=True):
    """Read ISCE .xml file input a python dictionary structure."""
    from lxml import objectify
    xmlDict = {}
    fObj = objectify.parse(fname)
    root = fObj.getroot()

    for child in root.findall('property'):
        xmlDict[child.attrib['name']] = str(child.value)

    # Read lat/lon info for geocoded file
    try:
        comp1 = root.find("./component[@name='coordinate1']")
        x_step = comp1.find("./property[@name='delta']/value").text
        if x_step not in ['1', '-1']:
            xmlDict['X_STEP'] = x_step
            xmlDict['X_FIRST'] = comp1.find("./property[@name='startingvalue']/value").text
            xmlDict['X_LAST'] = comp1.find("./property[@name='endingvalue']/value").text
    except:
        pass

    try:
        comp2 = root.find("./component[@name='coordinate2']")
        y_step = comp2.find("./property[@name='delta']/value").text
        if y_step not in ['1', '-1']:
            xmlDict['Y_STEP'] = y_step
            xmlDict['Y_FIRST'] = comp2.find("./property[@name='startingvalue']/value").text
            xmlDict['Y_LAST'] = comp2.find("./property[@name='endingvalue']/value").text
    except:
        pass

    if convert2roipac:
        xmlDict = attribute_isce2roipac(xmlDict)
    if standardize:
        xmlDict = standardize_metadata(xmlDict)
    return xmlDict


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

    # STARTING_RANGE
    key = 'near_range_slc'
    if key in par_dict_in.keys():
        par_dict['STARTING_RANGE'] = par_dict[key]

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

    # attributes in geo coordinates
    key = 'corner_lat'
    if key in par_dict_in.keys():
        par_dict['Y_FIRST'] = par_dict[key]

    key = 'corner_lon'
    if key in par_dict_in.keys():
        par_dict['X_FIRST'] = par_dict[key]

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


def attribute_isce2roipac(metaDict, dates=[], baselineDict={}):
    """Convert ISCE xml attribute into ROI_PAC format"""

    rscDict = {}
    for key in metaDict.keys():
        rscDict[key] = str(metaDict[key]).strip().split()[0]

    rscDict['WIDTH'] = rscDict['width']
    rscDict['LENGTH'] = rscDict['length']

    rscDict['PROCESSOR'] = 'isce'
    rscDict['PLATFORM'] = 'Sentinel1'

    rscDict['ANTENNA_SIDE'] = '-1'
    if 'passDirection' in rscDict.keys():
        rscDict['ORBIT_DIRECTION'] = rscDict['passDirection']

    if dates:
        rscDict['DATE12'] = str(dates[0][2:]+'-'+dates[1][2:])
        #rscDict['DATE'] = str(dates[0])

    if dates and baselineDict:
        bperp = baselineDict['bperp'][dates[1]] - baselineDict['bperp'][dates[0]]
        bpar = baselineDict['bpar'][dates[1]] - baselineDict['bpar'][dates[0]]
        rscDict['P_BASELINE_TOP_HDR'] = str(bperp)
        rscDict['P_BASELINE_BOTTOM_HDR'] = str(bperp)
        rscDict['H_BASELINE_TOP_HDR'] = str(bpar)
        rscDict['H_BASELINE_BOTTOM_HDR'] = str(bpar)

    return rscDict


def attribute_envi2roipac(metaDict):
    """Convert ISCE xml attribute into ROI_PAC format"""

    rscDict = {}
    for key in metaDict.keys():
        rscDict[key] = str(metaDict[key]).strip().split()[0]

    enviDataType = rscDict['data type']
    if enviDataType == '4':
        rscDict['DATA_TYPE'] = 'float32'
    return rscDict


#########################################################################
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
