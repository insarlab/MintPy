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

    # output file info
    fbase, fext = os.path.splitext(out_file)
    fext = fext.lower()
    out_dir = os.path.dirname(os.path.abspath(out_file))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        print('create directory: {}'.format(out_dir))

    # HDF5 File
    if fext in ['.h5', '.he5']:
        # grab info from reference h5 file
        if ref_file and os.path.splitext(ref_file)[1] in ['.h5', '.he5']:
            # compression
            if compression is None:
                compression = readfile.get_hdf5_compression(ref_file)

            # list of auxiliary datasets
            shape2d = (int(meta['LENGTH']), int(meta['WIDTH']))
            with h5py.File(ref_file, 'r') as fr:
                auxDsNames = [i for i in fr.keys()
                              if (i not in list(datasetDict.keys())
                                  and isinstance(fr[i], h5py.Dataset)
                                  and fr[i].shape[-2:] != shape2d)]
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

        # ignore certain meaningless file extensions
        while fext in ['.geo', '.rdr', '.full', '.wgs84']:
            fbase, fext = os.path.splitext(fbase)
        if not fext:
            fext = fbase

        # Write Data File
        print('write {}'.format(out_file))
        # determined by fext
        if fext in ['.unw']:
            write_float32(data_list[0], out_file)
            meta['DATA_TYPE'] = 'float32'

        elif fext in ['.cor', '.hgt']:
            if meta.get('PROCESSOR', 'isce') == 'roipac':
                write_float32(data_list[0], out_file)
            else:
                write_real_float32(data_list[0], out_file)
            meta['DATA_TYPE'] = 'float32'

        elif fext == '.dem':
            write_real_int16(data_list[0], out_file)
            meta['DATA_TYPE'] = 'int16'

        elif fext in ['.trans']:
            write_float32(data_list[0], data_list[1], out_file)
            meta['DATA_TYPE'] = 'float32'

        elif fext in ['.utm_to_rdc', '.UTM_TO_RDC']:
            data = np.zeros(data_list[0].shape, dtype=np.complex64)
            data.real = datasetDict['rangeCoord']
            data.imag = datasetDict['azimuthCoord']
            data.astype('>c8').tofile(out_file)

        elif fext in ['.mli', '.flt']:
            write_real_float32(data_list[0], out_file)

        elif fext == '.slc':
            write_complex_int16(data_list[0], out_file)

        elif fext == '.int':
            write_complex64(data_list[0], out_file)

        elif fext == '.msk':
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
            print('Un-supported file type:', fext)
            return 0

        # write metadata file
        write_roipac_rsc(meta, out_file+'.rsc', print_msg=True)
    return out_file


#########################################################################

def layout_hdf5(fname, ds_name_dict=None, metadata=None, ds_unit_dict=None, ref_file=None, compression=None, print_msg=True):
    """Create HDF5 file with defined metadata and (empty) dataset structure

    Parameters: fname        - str, HDF5 file path
                ds_name_dict - dict, dataset structure definition
                               {dname : [dtype, dshape],
                                dname : [dtype, dshape, None],
                                dname : [dtype, dshape, 1/2/3D np.ndarray], #for aux data
                                ...
                               }
                metadata     - dict, metadata
                ds_unit_dict - dict, dataset unit definition
                               {dname : dunit,
                                dname : dunit,
                                ...
                               }
                ref_file     - str, reference file for the data structure
                compression  - str, HDF5 compression type
    Returns:    fname        - str, HDF5 file path

    Example:    layout_hdf5('timeseries_ERA5.h5', ref_file='timeseries.h5')
                layout_hdf5('timeseries_ERA5.5h', ds_name_dict, metadata)

    # structure for ifgramStack
    ds_name_dict = {
        "date"             : [np.dtype('S8'), (num_ifgram, 2)],
        "dropIfgram"       : [np.bool_,       (num_ifgram,)],
        "bperp"            : [np.float32,     (num_ifgram,)],
        "unwrapPhase"      : [np.float32,     (num_ifgram, length, width)],
        "coherence"        : [np.float32,     (num_ifgram, length, width)],
        "connectComponent" : [np.int16,       (num_ifgram, length, width)],
    }

    # structure for geometry
    ds_name_dict = {
        "height"             : [np.float32, (length, width), None],
        "incidenceAngle"     : [np.float32, (length, width), None],
        "slantRangeDistance" : [np.float32, (length, width), None],
    }

    # structure for timeseries
    dates = np.array(date_list, np.string_)
    ds_name_dict = {
        "date"       : [np.dtype("S8"), (num_date,), dates],
        "bperp"      : [np.float32,     (num_date,), pbase],
        "timeseries" : [np.float32,     (num_date, length, width)],
    }
    """

    print('-'*50)

    # get meta from metadata and ref_file
    if metadata:
        meta = {key: value for key, value in metadata.items()}
    elif ref_file:
        with h5py.File(ref_file, 'r') as fr:
            meta = {key: value for key, value in fr.attrs.items()}
        if print_msg:
            print('grab metadata from ref_file: {}'.format(ref_file))
    else:
        raise ValueError('No metadata or ref_file found.')

    # check ds_name_dict
    if ds_name_dict is None:
        ds_name_dict = {}

        if ref_file and os.path.splitext(ref_file)[1] in ['.h5', '.he5']:
            with h5py.File(ref_file, 'r') as fr:
                shape2d = (int(fr.attrs['LENGTH']), int(fr.attrs['WIDTH']))
                shape2d_out = (int(meta['LENGTH']), int(meta['WIDTH']))

                for key in fr.keys():
                    ds = fr[key]
                    if isinstance(ds, h5py.Dataset):

                        # auxliary dataset
                        if ds.shape[-2:] != shape2d:
                            ds_name_dict[key] = [ds.dtype, ds.shape, ds[:]]

                        # dataset
                        else:
                            ds_shape = list(ds.shape)
                            ds_shape[-2:] = shape2d_out
                            ds_name_dict[key] = [ds.dtype, tuple(ds_shape), None]

            if print_msg:
                print('grab dataset structure from ref_file: {}'.format(ref_file))
        else:
            raise ValueError('No ds_name_dict or ref_file found.')

    # directory
    fdir = os.path.dirname(os.path.abspath(fname))
    if not os.path.isdir(fdir):
        os.makedirs(fdir)
        print('crerate directory: {}'.format(fdir))

    # create file
    with h5py.File(fname, "w") as f:
        if print_msg:
            print('create HDF5 file: {} with w mode'.format(fname))

        # initiate dataset
        max_digit = max([len(i) for i in ds_name_dict.keys()])
        for key in ds_name_dict.keys():
            data_type  = ds_name_dict[key][0]
            data_shape = ds_name_dict[key][1]

            # turn ON compression for conn comp
            ds_comp = compression
            if key in ['connectComponent']:
                ds_comp = 'lzf'

            # changable dataset shape
            if len(data_shape) == 3:
                max_shape = (None, data_shape[1], data_shape[2])
            else:
                max_shape = data_shape

            # create empty dataset
            if print_msg:
                print(("create dataset  : {d:<{w}} of {t:<25} in size of {s:<20} with "
                       "compression = {c}").format(d=key,
                                                   w=max_digit,
                                                   t=str(data_type),
                                                   s=str(data_shape),
                                                   c=ds_comp))
            ds = f.create_dataset(key,
                                  shape=data_shape,
                                  maxshape=max_shape,
                                  dtype=data_type,
                                  chunks=True,
                                  compression=ds_comp)

            # write auxliary data
            if len(ds_name_dict[key]) > 2 and ds_name_dict[key][2] is not None:
                ds[:] = np.array(ds_name_dict[key][2])

        # write attributes in root level
        for key, value in meta.items():
            f.attrs[key] = str(value)

        # write attributes in dataset level
        if ds_unit_dict is not None:
            for key, value in ds_unit_dict.items():
                if value is not None:
                    f[key].attrs['UNIT'] = value
                    print(f'add /{key:<{max_digit}} attribute: UNIT = {value}')

    if print_msg:
        print('close  HDF5 file: {}'.format(fname))

    return fname


def write_hdf5_block(fname, data, datasetName, block=None, mode='a', print_msg=True):
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
    if print_msg:
        print('-'*50)
        print('open  HDF5 file {} in {} mode'.format(fname, mode))
        print("writing dataset /{:<25} block: {}".format(datasetName, block))
    with h5py.File(fname, mode) as f:
        if len(block) == 6:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3],
                           block[4]:block[5]] = data

        elif len(block) == 4:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3]] = data

        elif len(block) == 2:
            f[datasetName][block[0]:block[1]] = data

    if print_msg:
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



#########################################################################

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


def write_isce_xml(fname, width, length, bands=1, data_type='FLOAT', scheme='BIP'):
    """Write XML metadata file in ISCE-2 format

    Parameters: fname     - str, path of data file
                width     - int, number of columns
                length    - int, number of rows
                bands     - int, number of band
                data_type - str, data type name in ISCE convention
                            readfile.GDAL2ISCE_DATATYPE
                scheme    - str, band interleave, BIP, BIL, BSQ
    Examples:   write_isce_xml(out_file='filt_fine.cor', width=400, length=500)
    """
    import isce
    import isceobj

    img = isceobj.Image.createImage()
    img.setFilename(fname)
    img.setWidth(width)
    img.setLength(length)
    img.setAccessMode('READ')
    img.bands = bands
    img.dataType = data_type
    img.scheme = scheme
    img.renderHdr()
    img.renderVRT()

    return


def write_isce_file(data, out_file, file_type='isce_unw'):
    """write data to file in ISCE format

    Parameters: data      - 2D np.ndarray, binary data matrix
                out_file  - str, path of output binary data file
                file_type - str, file type
    Returns:    out_file  - str, path of output binary data file
    """
    import isce
    import isceobj

    # fix potential typo
    file_type = file_type.replace('-', '_')

    # write data to binary file
    data.tofile(out_file)

    # write isce xml metadata file
    length, width = data.shape

    if file_type == 'isce_unw':
        width = int(width / 2)
        write_isce_xml(out_file, width, length, bands=2, data_type='FLOAT', scheme='BIL')

    elif file_type == 'isce_int':
        write_isce_xml(out_file, width, length, bands=1, data_type='CFLOAT', scheme='BIL')

    elif file_type == 'isce_cor':
        write_isce_xml(out_file, width, length, bands=1, data_type='FLOAT', scheme='BIL')

    else:
        raise ValueError('un-recognized ISCE file type: {}'.format(file_type))

    return out_file



#########################################################################

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
