############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.utils import writefile


import os
import shutil
import h5py
import numpy as np
from mintpy.objects import timeseries
from mintpy.utils import readfile


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
    # copy metadata to meta
    if metadata:
        meta = {key: value for key, value in metadata.items()}
    elif ref_file:
        meta = readfile.read_attribute(ref_file)
    else:
        raise ValueError('No metadata or reference file input.')

    # convert ndarray input into dict type
    if isinstance(datasetDict, np.ndarray):
        data = np.array(datasetDict, datasetDict.dtype)
        datasetDict = dict()
        datasetDict[meta['FILE_TYPE']] = data

    ext = os.path.splitext(out_file)[1].lower()
    # HDF5 File
    if ext in ['.h5', '.he5']:
        # grab info from reference h5 file
        if ref_file and os.path.splitext(ref_file)[1] in ['.h5', '.he5']:
            # compression
            if compression is None:
                compression = readfile.get_hdf5_compression(ref_file)

            # list of auxiliary datasets
            atr_ref = readfile.read_attribute(ref_file)
            shape_ref = (int(atr_ref['LENGTH']), int(atr_ref['WIDTH']))
            with h5py.File(ref_file, 'r') as fr:
                auxDsNames = [i for i in fr.keys()
                              if (i not in list(datasetDict.keys())
                                  and isinstance(fr[i], h5py.Dataset)
                                  and fr[i].shape[-2:] != shape_ref)]
        else:
            auxDsNames = []

        # check required datasets
        dsNames = list(datasetDict.keys()) + auxDsNames
        if meta['FILE_TYPE'] in ['timeseries', 'ifgramStack']:
            if 'date' not in dsNames:
                raise Exception("Can not write {} file without 'date' dataset!".format(meta['FILE_TYPE']))

        # remove existing file
        if os.path.isfile(out_file):
            print('delete exsited file: {}'.format(out_file))
            os.remove(out_file)

        # writing
        print('create HDF5 file: {} with w mode'.format(out_file))
        maxDigit = max([len(i) for i in dsNames])
        with h5py.File(out_file, 'w') as f:
            # 1. write input datasets
            for dsName in datasetDict.keys():
                data = datasetDict[dsName]
                print(('create dataset /{d:<{w}} of {t:<10} in size of {s:<20} '
                       'with compression={c}').format(d=dsName,
                                                      w=maxDigit,
                                                      t=str(data.dtype),
                                                      s=str(data.shape),
                                                      c=compression))
                ds = f.create_dataset(dsName,
                                      data=data,
                                      chunks=True,
                                      compression=compression)

            # 2. Write extra/auxliary datasets from ref_file
            if len(auxDsNames) > 0:
                with h5py.File(ref_file, 'r') as fr:
                    for dsName in auxDsNames:
                        ds = fr[dsName]
                        print(('create dataset /{d:<{w}} of {t:<10} in size of {s:<10} '
                               'with compression={c}').format(d=dsName,
                                                              w=maxDigit,
                                                              t=str(ds.dtype),
                                                              s=str(ds.shape),
                                                              c=compression))
                        f.create_dataset(dsName,
                                         data=ds[:],
                                         chunks=True,
                                         compression=compression)

            # 3. metadata
            for key, value in meta.items():
                try:
                    f.attrs[key] = str(value)
                except:
                    f.attrs[key] = str(value.encode('utf-8'))
            print('finished writing to {}'.format(out_file))

    # ISCE / ROI_PAC GAMMA / Image product
    else:
        key_list = list(datasetDict.keys())
        data_list = []
        for key in key_list:
            data_list.append(datasetDict[key])
        data_type = meta.get('DATA_TYPE', str(data_list[0].dtype)).lower()

        # Write Data File
        print('write {}'.format(out_file))
        # determined by ext
        if ext in ['.unw', '.cor', '.hgt']:
            write_float32(data_list[0], out_file)
            meta['DATA_TYPE'] = 'float32'

        elif ext == '.dem':
            write_real_int16(data_list[0], out_file)
            meta['DATA_TYPE'] = 'int16'

        elif ext in ['.trans']:
            write_float32(data_list[0], data_list[1], out_file)
            meta['DATA_TYPE'] = 'float32'

        elif ext in ['.utm_to_rdc', '.UTM_TO_RDC']:
            data = np.zeros(data_list[0].shape, dtype=np.complex64)
            data.real = datasetDict['rangeCoord']
            data.imag = datasetDict['azimuthCoord']
            data.astype('>c8').tofile(out_file)

        elif ext in ['.mli', '.flt']:
            write_real_float32(data_list[0], out_file)

        elif ext == '.slc':
            write_complex_int16(data_list[0], out_file)

        elif ext == '.int':
            write_complex64(data_list[0], out_file)

        elif ext == '.msk':
            write_byte(data_list[0], out_file)
            meta['DATA_TYPE'] = 'byte'

        # determined by DATA_TYPE
        elif data_type in ['float64']:
            write_real_float64(data_list[0], out_file)

        elif data_type in ['float32', 'float']:
            if len(data_list) == 1:
                write_real_float32(data_list[0], out_file)

            elif len(data_list) == 2 and meta['scheme'] == 'BIL':
                write_float32(data_list[0], data_list[1], out_file)

        elif data_type in ['int16', 'short']:
            write_real_int16(data_list[0], out_file)

        elif data_type in ['int8', 'byte']:
            write_byte(data_list[0], out_file)

        elif data_type in ['bool']:
            write_bool(data_list[0], out_file)

        else:
            print('Un-supported file type: '+ext)
            return 0

        # write metadata file
        write_roipac_rsc(meta, out_file+'.rsc', print_msg=True)
    return out_file


def layout_hdf5(fname, dsNameDict, metadata):
    """Create HDF5 file with defined metadata and (empty) dataset structure

    Parameters: fname      - str, HDF5 file path
                dsNameDict - dict, dataset structure definition, as below:
                metadata   - dict, metadata
    Returns:    fname      - str, HDF5 file path

    Example:

    # structure for ifgramStack
    dsNameDict = {
        "date"             : (np.dtype('S8'), (inps.num_pair, 2)),
        "dropIfgram"       : (np.bool_,       (inps.num_pair,)),
        "bperp"            : (np.float32,     (inps.num_pair,)),
        "unwrapPhase"      : (np.float32,     (inps.num_pair, inps.length, inps.width)),
        "coherence"        : (np.float32,     (inps.num_pair, inps.length, inps.width)),
        "connectComponent" : (np.int16,       (inps.num_pair, inps.length, inps.width)),
    }

    # structure for geometry
    dsNameDict = {
        "height"             : (np.float32, (inps.length, inps.width)),
        "incidenceAngle"     : (np.float32, (inps.length, inps.width)),
        "slantRangeDistance" : (np.float32, (inps.length, inps.width)),
    }

    # structure for timeseries
    dsNameDict = {
        "date"       : (np.dtype("S8"), (numDates,)),
        "bperp"      : (np.float32,     (numDates,)),
        "timeseries" : (np.float32,     (numDates, length, width))
    }
    """

    print('-'*50)
    print('create HDF5 file {} with w mode'.format(fname))
    h5 = h5py.File(fname, "w")

    # initiate dataset
    for key in dsNameDict.keys():
        data_type = dsNameDict[key][0]
        data_shape = dsNameDict[key][1]

        # turn ON compression
        if key in ['connectComponent']:
            compression = 'lzf'
        else:
            compression = None

        # changable dataset shape
        if len(data_shape) == 3:
            maxShape = (None, data_shape[1], data_shape[2])
        else:
            maxShape = data_shape

        print("create dataset: {d:<25} of {t:<25} in size of {s}".format(d=key,
                                                                         t=str(data_type),
                                                                         s=data_shape))
        h5.create_dataset(key,
                          shape=data_shape,
                          maxshape=maxShape,
                          dtype=data_type,
                          chunks=True,
                          compression=compression)

    # write attributes
    for key in metadata.keys():
        h5.attrs[key] = metadata[key]

    h5.close()
    print('close  HDF5 file {}'.format(fname))
    return fname


def write_hdf5_block(fname, data, datasetName, block=None, mode='a'):
    """Write data to existing HDF5 dataset in disk block by block.
    Parameters: data        - np.ndarray 1/2/3D matrix
                datasetName - str, dataset name
                block       - list of 2/4/6 int, for
                              [zStart, zEnd,
                               yStart, yEnd,
                               xStart, xEnd]
                mode        - str, open mode
    Returns:    fname
    """

    # default block value
    if block is None:

        # data shape
        if isinstance(data, list):
            shape=(len(data),)
        else:
            shape = data.shape

        # set default block as the entire data
        if len(shape) ==1:
            block = [0, shape[0]]
        elif len(shape) == 2:
            block = [0, shape[0],
                     0, shape[1]]
        elif len(shape) == 3:
            block = [0, shape[0],
                     0, shape[1],
                     0, shape[2]]

    # write
    print('-'*50)
    print('open  HDF5 file {} in {} mode'.format(fname, mode))
    with h5py.File(fname, mode) as f:
        print("writing dataset /{:<25} block: {}".format(datasetName, block))

        if len(block) == 6:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3],
                           block[4]:block[5]] = data

        elif len(block) == 4:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3]] = data

        elif len(block) == 2:
            f[datasetName][block[0]:block[1]] = data

    print('close HDF5 file {}.'.format(fname))
    return fname


def remove_hdf5_dataset(fname, datasetNames, print_msg=True):
    """Remove an existing dataset from an HDF5 file.
    Parameters: fname : str, HDF5 file name/path
                datasetName : (list of) str, dataset name(s)
    Returns:    fname : str,
    Example:    remove_hdf5_dataset('./inputs/ifgramStack.h5', 'unwrapPhase_phaseClosure')
                remove_hdf5_dataset('./inputs/ifgramStack.h5', ['unwrapPhase_phaseClosure',
                                                                'unwrapPhase_bridging'])
    """
    if isinstance(datasetNames, str):
        datasetNames = list(datasetNames)
    if print_msg:
        print('delete {} from file {}'.format(datasetNames, fname))
    # 1. rename the file to a temporary file
    temp_file = os.path.join(os.path.dirname(fname), 'tmp_{}'.format(os.path.basename(fname)))
    print('move {} to {}'.format(fname, temp_file))
    shutil.move(fname, temp_file)

    # 2. write a new file with all data except for the one to be deleted
    if print_msg:
        print('read   HDF5 file: {} with r mode'.format(temp_file))
        print('create HDF5 file: {} with w mode'.format(fname))
    fi = h5py.File(temp_file, 'r')
    fo = h5py.File(fname, 'w')

    # datasets
    compression = None
    maxDigit = max([len(i) for i in list(fi.keys())])
    for dsName in [i for i in fi.keys() if i not in datasetNames]:
        ds = fi[dsName]
        if print_msg:
            print('create dataset /{d:<{w}} of {t:<10} in size of {s:<20} with compression={c}'.format(
                d=dsName, w=maxDigit, t=str(ds.dtype), s=str(ds.shape), c=compression))
        fo.create_dataset(dsName, data=ds[:], chunks=True, compression=compression)

    # metadata
    for key, value in fi.attrs.items():
        fo.attrs[key] = str(value)
    fi.close()
    fo.close()
    if print_msg:
        print('finished writing to {}'.format(fname))
        print('old file is now saved as: {}. Use rm command to delete it.'.format(temp_file))
    return fname


def write_roipac_rsc(metadata, out_file, update_mode=False, print_msg=False):
    """Write attribute dict into ROI_PAC .rsc file
    Inputs:
        metadata : dict, attributes dictionary
        out_file : rsc file name, to which attribute is writen
        update_mode : bool, skip writing if
                      1) output file existed AND
                      2) no new metadata key/value
        print_msg   : bool, print message
    Output:
        out_file
    """
    run = True
    if update_mode:
        rsc_dict = dict()
        if os.path.isfile(out_file):
            rsc_dict = readfile.read_roipac_rsc(out_file)
        # update .rsc file only if there are new metadata key/value
        if set(metadata.items()).issubset(set(rsc_dict.items())):
            run = False

    if run:
        # Convert MintPy attributes to ROI_PAC attributes
        if 'LENGTH' in metadata.keys():
            metadata['FILE_LENGTH'] = metadata['LENGTH']

        # Convert 3.333e-4 to 0.0003333
        if 'X_STEP' in metadata.keys():
            metadata['X_STEP'] = str(float(metadata['X_STEP']))
            metadata['Y_STEP'] = str(float(metadata['Y_STEP']))
            metadata['X_FIRST'] = str(float(metadata['X_FIRST']))
            metadata['Y_FIRST'] = str(float(metadata['Y_FIRST']))

        # writing .rsc file
        if print_msg:
            print('write', out_file)
        maxDigit = max([len(key) for key in metadata.keys()]+[2])
        with open(out_file, 'w') as f:
            for key in sorted(metadata.keys()):
                f.write('{k:<{d}}    {v}\n'.format(k=str(key),
                                                   d=maxDigit,
                                                   v=str(metadata[key])))
    return out_file


def write_float32(*args):
    """Write ROI_PAC rmg format with float32 precision (BIL)
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
    data = np.array(data, dtype=np.float32)
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


def write_real_float64(data, out_file):
    """write isce float data, i.e. hgt.rdr file."""
    data = np.array(data, dtype=np.float64)
    data.tofile(out_file)
    return out_file


def write_real_float32(data, out_file):
    """write gamma float data, i.e. .mli file."""
    data = np.array(data, dtype=np.float32)
    data.tofile(out_file)
    return out_file


def write_real_int16(data, out_file):
    data = np.array(data, dtype=np.int16)
    data.tofile(out_file)
    return out_file


def write_byte(data, out_file):
    data = np.array(data, dtype=np.byte)
    data.tofile(out_file)
    return out_file


def write_bool(data, out_file):
    data = np.array(data, dtype=np.bool_)
    data.tofile(out_file)
    return out_file
