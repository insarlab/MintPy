############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend usage:
#   from pysar.utils import writefile
#


import os
import sys
import h5py
import numpy as np
#from PIL import Image
from pysar.objects import timeseries
from pysar.utils import readfile


def write(datasetDict, out_file, metadata=None, ref_file=None, compression=None):
    """ Write one file.
    Parameters: datasetDict : dict of dataset, with key = datasetName and value = 2D/3D array, e.g.:
                    {'height'        : np.ones((   200,300), dtype=np.int16),
                     'incidenceAngle': np.ones((   200,300), dtype=np.float32),
                     'bperp'         : np.ones((80,200,300), dtype=np.float32),
                     ...}
                out_file : str, output file name
                metadata : dict of attributes
                ref_file : str, reference file to get auxliary info
                compression : str, compression while writing to HDF5 file, None, "lzf", "gzip"
    Returns:    out_file : str
    Examples:   dsDict = dict()
                dsDict['velocity'] = np.ones((200,300), dtype=np.float32)
                write(datasetDict=dsDict, out_file='velocity.h5', metadata=atr)
    """
    ext = os.path.splitext(out_file)[1].lower()
    if ref_file and metadata is None:
        metadata = readfile.read_attribute(ref_file)

    if type(datasetDict) is np.ndarray:
        data = np.array(datasetDict)
        datasetDict = dict()
        datasetDict[metadata['FILE_TYPE']] = data

    # HDF5 File
    if ext in ['.h5', '.he5']:
        k = metadata['FILE_TYPE']
        if k == 'timeseries':
            if ref_file is None:
                raise Exception('Can not write {} file without reference file!'.format(k))
            obj = timeseries(out_file)
            obj.write2hdf5(datasetDict[k],
                           metadata=metadata,
                           refFile=ref_file)

        else:
            if os.path.isfile(out_file):
                print('delete exsited file: {}'.format(out_file))
                os.remove(out_file)
            print('create HDF5 file: {} with w mode'.format(out_file))
            f = h5py.File(out_file, 'w')

            # Write input datasets
            maxDigit = max([len(i) for i in list(datasetDict.keys())])
            for dsName in datasetDict.keys():
                data = datasetDict[dsName]
                print(('create dataset /{d:<{w}} of {t:<10}'
                       ' in size of {s}').format(d=dsName,
                                                 w=maxDigit,
                                                 t=str(data.dtype),
                                                 s=data.shape))
                ds = f.create_dataset(dsName,
                                      data=data,
                                      chunks=True,
                                      compression=compression)

            # Write extra/auxliary datasets from ref_file
            if ref_file:
                fr = h5py.File(ref_file, 'r')
                dsNames = [i for i in fr.keys()
                           if (i not in list(datasetDict.keys())
                               and isinstance(fr[i], h5py.Dataset))]
                for dsName in dsNames:
                    ds = fr[dsName]
                    print(('create dataset /{d:<{w}} of {t:<10}'
                           ' in size of {s}').format(d=dsName,
                                                     w=maxDigit,
                                                     t=str(ds.dtype),
                                                     s=ds.shape))
                    f.create_dataset(dsName,
                                     data=ds[:],
                                     chunks=True,
                                     compression=compression)
                fr.close()

            # metadata
            for key, value in metadata.items():
                f.attrs[key] = str(value)
            f.close()
            print('finished writing to {}'.format(out_file))

    # ISCE / ROI_PAC GAMMA / Image product
    else:
        # Write Data File
        if ext in ['.unw', '.cor', '.hgt']:
            write_float32(data, out_file)
        elif ext == '.dem':
            write_real_int16(data, out_file)
        elif ext in ['.trans']:
            write_float32(rg, az, out_file)
        elif ext in ['.utm_to_rdc', '.UTM_TO_RDC']:
            data = np.zeros(rg.shape, dtype=np.complex64)
            data.real = datasetDict['rangeCoord']
            data.imag = datasetDict['azimuthCoord']
            data.astype('>c8').tofile(out_file)
        # elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        #    data.save(out_file)
        elif ext == '.mli':
            write_real_float32(data, out_file)
        elif ext == '.slc':
            write_complex_int16(data, out_file)
        elif ext == '.int':
            write_complex64(data, out_file)
        elif metadata['DATA_TYPE'].lower() in ['float32', 'float']:
            write_real_float32(data, out_file)
        elif metadata['DATA_TYPE'].lower() in ['int16', 'short']:
            write_real_int16(data, out_file)
        else:
            print('Un-supported file type: '+ext)
            return 0

        # Write .rsc File
        write_roipac_rsc(metadata, out_file+'.rsc')
        return out_file


def write_roipac_rsc(metadata, out_file, sorting=True):
    """Write attribute dict into ROI_PAC .rsc file
    Inputs:
        metadata     - dict, attributes dictionary
        out_file - rsc file name, to which attribute is writen
        sorting - bool, sort attributes in alphabetic order while writing
    Output:
        out_file
    """
    # Convert PYSAR attributes to ROI_PAC attributes
    metadata['FILE_LENGTH'] = metadata['LENGTH']

    # Convert 3.333e-4 to 0.0003333
    if 'X_STEP' in metadata.keys():
        metadata['X_STEP'] = str(float(metadata['X_STEP']))
        metadata['Y_STEP'] = str(float(metadata['Y_STEP']))
        metadata['X_FIRST'] = str(float(metadata['X_FIRST']))
        metadata['Y_FIRST'] = str(float(metadata['Y_FIRST']))

    # sorting by key name
    dictKey = metadata.keys()
    if sorting:
        dictKey = sorted(dictKey)

    # writing .rsc file
    maxDigit = max([len(key) for key in metadata.keys()]+[2])
    f = open(out_file, 'w')
    for key in dictKey:
        f.write('{k:<{d}}    {v}\n'.format(k=str(key),
                                           d=maxDigit,
                                           v=str(metadata[key])))
    f.close()
    return out_file


def write_float32(*args):
    """Write ROI_PAC rmg format with float32 precision
    Format of the binary file is same as roi_pac unw, cor, or hgt data.
          should rename to write_rmg_float32()

    Exmaple:
            write_float32(phase, out_file)
            write_float32(amp, phase, out_file)
    """
    if len(args) == 2:
        amp = args[0]
        pha = args[0]
        out_file = args[1]
    elif len(args) == 3:
        amp = args[0]
        pha = args[1]
        out_file = args[2]
    else:
        print('Error while getting args: support 2/3 args only.')
        return

    data = np.hstack((amp, pha)).flatten()
    data.tofile(out_file)
    return out_file


def write_complex64(data, out_file):
    """Writes roi_pac .int data"""
    num_pixel = data.size
    F = np.zeros([2 * num_pixel, 1], np.float32)
    id1 = list(range(0, 2 * num_pixel, 2))
    id2 = list(range(1, 2 * num_pixel, 2))
    F[id1] = np.reshape(np.cos(data), (num_pixel, 1))
    F[id2] = np.reshape(np.sin(data), (num_pixel, 1))
    F.tofile(out_file)
    return out_file


def write_real_int16(data, out_file):
    data = np.array(data, dtype=np.int16)
    data.tofile(out_file)
    return out_file


def write_dem(data, out_file):
    data = np.array(data, dtype=np.int16)
    data.tofile(out_file)
    return out_file


def write_real_float32(data, out_file):
    """write gamma float data, i.e. .mli file."""
    data = np.array(data, dtype=np.float32)
    data.tofile(out_file)
    return out_file


def write_complex_int16(data, out_file):
    """Write gamma scomplex data, i.e. .slc file.
        data is complex 2-D matrix
        real, imagery, real, ...
    """
    num_pixel = data.size
    id1 = list(range(0, 2 * num_pixel, 2))
    id2 = list(range(1, 2 * num_pixel, 2))

    F = np.zeros([2 * num_pixel, 1], np.int16)
    F[id1] = np.reshape(np.array(data.real, np.int16), (num_pixel, 1))
    F[id2] = np.reshape(np.array(data.imag, np.int16), (num_pixel, 1))
    F.tofile(out_file)
    return out_file
