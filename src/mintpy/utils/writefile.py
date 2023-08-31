"""Utilities to write files."""
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

from mintpy.utils import readfile


def write(datasetDict, out_file, metadata=None, ref_file=None, compression=None, ds_unit_dict=None, print_msg=True):
    """ Write one file.
    Parameters: datasetDict  - dict of dataset, with key = datasetName and value = 2D/3D array, e.g.:
                    {'height'        : np.ones((   200,300), dtype=np.int16),
                     'incidenceAngle': np.ones((   200,300), dtype=np.float32),
                     'bperp'         : np.ones((80,200,300), dtype=np.float32),
                     ...}
                out_file     - str, output file name
                metadata     - dict of attributes
                ref_file     - str, reference file to get auxliary info
                compression  - str, compression while writing to HDF5 file, None, "lzf", "gzip"
                ds_unit_dict - dict, top-level dataset unit definition
                    {dname : dunit,
                     dname : dunit,
                     ...
                    }
    Returns:    out_file     - strs
    Examples:   dsDict = dict()
                dsDict['velocity'] = np.ones((200,300), dtype=np.float32)
                write(datasetDict=dsDict, out_file='velocity.h5', metadata=atr)
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

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

    # file extension
    fbase, fext = os.path.splitext(out_file)
    # ignore certain meaningless file extensions
    while fext in ['.geo', '.rdr', '.full', '.wgs84']:
        fbase, fext = os.path.splitext(fbase)
    if not fext:
        fext = fbase
    fext = fext.lower()

    # output file info
    out_dir = os.path.dirname(os.path.abspath(out_file))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        vprint(f'create directory: {out_dir}')

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
            os.remove(out_file)
            vprint(f'delete exsited file: {out_file}')

        # writing
        print(f'create HDF5 file: {out_file} with w mode')
        maxDigit = max(len(i) for i in dsNames)
        with h5py.File(out_file, 'w') as f:
            # 1. write input datasets
            for dsName in datasetDict.keys():
                data = datasetDict[dsName]
                vprint(('create dataset /{d:<{w}} of {t:<10} in size of {s:<20} '
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
                        vprint(('create dataset /{d:<{w}} of {t:<10} in size of {s:<10} '
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

            # write attributes in dataset level
            if ds_unit_dict is not None:
                for key, value in ds_unit_dict.items():
                    if key in f.keys() and value is not None:
                        f[key].attrs['UNIT'] = value
                        vprint(f'add /{key:<{maxDigit}} attribute: UNIT = {value}')

        vprint(f'finished writing to {out_file}')

    # ISCE / ROI_PAC GAMMA / Image product
    else:
        # basic info
        key_list = list(datasetDict.keys())
        data_list = list(datasetDict.values())
        meta['BANDS'] = len(key_list)
        meta['INTERLEAVE'] = meta.get('INTERLEAVE', 'BIL').upper()
        # data type
        meta['DATA_TYPE'] = meta.get('DATA_TYPE', 'float32').lower()
        DATA_TYPE_DICT = {'float' : 'float32',
                          'short' : 'int16',
                          'byte'  : 'int8'}
        if meta['DATA_TYPE'] in DATA_TYPE_DICT.keys():
            meta['DATA_TYPE'] = DATA_TYPE_DICT[meta['DATA_TYPE']]

        # adjust for pre-defined files determined by fext
        if fext in ['.unw']:
            meta['DATA_TYPE'] = 'float32'
            meta['INTERLEAVE'] = 'BIL'
            if key_list != ['magnitude', 'phase']:
                data_list = [data_list[0], data_list[0]]
                meta['BANDS'] = 2

        elif fext in ['.cor']:
            # remove .hgt as it can be float64 in isce2.
            meta['DATA_TYPE'] = 'float32'
            meta['INTERLEAVE'] = 'BIL'
            if meta.get('PROCESSOR', 'isce') == 'roipac':
                data_list = [data_list[0], data_list[0]]

        elif fext == '.dem':
            meta['DATA_TYPE'] = 'int16'

        elif fext in ['.trans']:
            # ROI_PAC lookup table
            meta['DATA_TYPE'] = 'float32'
            meta['INTERLEAVE'] = 'BIL'

        elif fext.endswith(('to_rdc', '2_rdc', '2rdc')):
            # Gamma lookup table
            meta['BANDS'] = 2
            meta['DATA_TYPE'] = 'float32'
            meta['INTERLEAVE'] = 'BIP'

        elif fext in ['.mli', '.flt']:
            meta['DATA_TYPE'] = 'float32'

        elif fext == '.slc':
            # SLC: complex 64 or 32
            meta['BANDS'] = 1

        elif fext == '.int':
            # wrapped interferogram: complex 64
            meta['BANDS'] = 1
            meta['DATA_TYPE'] = 'complex64'
            if key_list == ['magnitude', 'phase']:
                data_list[0] = data_list[0] * np.exp(1j * data_list[1])

        elif fext == '.msk':
            meta['DATA_TYPE'] = 'int8'

        data_types = ['bool', 'int8', 'uint8', 'int16', 'float32', 'float64', 'complex32', 'complex64', 'complex128']
        if meta['DATA_TYPE'] not in data_types:
            msg = 'Un-supported file type "{}" with data type "{}"!'.format(fext, meta['DATA_TYPE'])
            msg += f'\nSupported data type list: {data_types}'
            raise ValueError(msg)

        # write binary file
        write_binary(data_list, out_file, data_type=meta['DATA_TYPE'], interleave=meta['INTERLEAVE'])
        vprint(f'write file: {out_file}')

        # write metadata file
        write_roipac_rsc(meta, out_file+'.rsc', print_msg=print_msg)

    return out_file


#########################################################################

def layout_hdf5(fname, ds_name_dict=None, metadata=None, ds_unit_dict=None, ref_file=None, compression=None, print_msg=True):
    """Create HDF5 file with defined metadata and (empty) dataset structure

    Parameters: fname        - str, HDF5 file path
                ds_name_dict - dict, dataset structure definition
                               {dname : [dtype, dshape],
                                dname : [dtype, dshape, None],
                                dname : [dtype, dshape, 1/2/3/4D np.ndarray], #for aux data
                                ...
                               }
                metadata     - dict, metadata
                ds_unit_dict - dict, top-level dataset unit definition
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
    vprint = print if print_msg else lambda *args, **kwargs: None
    vprint('-'*50)

    # get meta from metadata and ref_file
    if metadata:
        meta = {key: value for key, value in metadata.items()}
    elif ref_file:
        with h5py.File(ref_file, 'r') as fr:
            meta = {key: value for key, value in fr.attrs.items()}
        vprint(f'grab metadata from ref_file: {ref_file}')
    else:
        raise ValueError('No metadata or ref_file found.')

    # check ds_name_dict
    if ds_name_dict is None:
        if not ref_file or not os.path.isfile(ref_file):
            raise FileNotFoundError('No ds_name_dict or ref_file found!')
        else:
            vprint(f'grab dataset structure from ref_file: {ref_file}')

        ds_name_dict = {}
        fext = os.path.splitext(ref_file)[1]
        shape2d = (int(meta['LENGTH']), int(meta['WIDTH']))

        if fext in ['.h5', '.he5']:
            # copy dset structure from HDF5 file
            with h5py.File(ref_file, 'r') as fr:
                # in case output mat size is different from the input ref file mat size
                shape2d_orig = (int(fr.attrs['LENGTH']), int(fr.attrs['WIDTH']))

                for key in fr.keys():
                    ds = fr[key]
                    if isinstance(ds, h5py.Dataset):

                        # auxliary dataset
                        if ds.shape[-2:] != shape2d_orig:
                            ds_name_dict[key] = [ds.dtype, ds.shape, ds[:]]

                        # dataset
                        else:
                            ds_shape = list(ds.shape)
                            ds_shape[-2:] = shape2d
                            ds_name_dict[key] = [ds.dtype, tuple(ds_shape), None]

        else:
            # construct dset structure from binary file
            ds_names = readfile.get_slice_list(ref_file)
            ds_dtype = meta['DATA_TYPE']
            for ds_name in ds_names:
                ds_name_dict[ds_name] = [ds_dtype, tuple(shape2d), None]

    # directory
    fdir = os.path.dirname(os.path.abspath(fname))
    if not os.path.isdir(fdir):
        os.makedirs(fdir)
        vprint(f'crerate directory: {fdir}')

    # create file
    with h5py.File(fname, "w") as f:
        vprint(f'create HDF5 file: {fname} with w mode')

        # initiate dataset
        max_digit = max(len(i) for i in ds_name_dict.keys())
        for key in ds_name_dict.keys():
            data_type  = ds_name_dict[key][0]
            data_shape = ds_name_dict[key][1]

            # turn ON compression for conn comp
            ds_comp = compression
            if key in ['connectComponent']:
                ds_comp = 'lzf'

            # changeable dataset shape
            if len(data_shape) == 3:
                max_shape = (None, data_shape[1], data_shape[2])
            else:
                max_shape = data_shape

            # create empty dataset
            vprint(("create dataset  : {d:<{w}} of {t:<25} in size of {s:<20} with "
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
                if key in f.keys() and value is not None:
                    f[key].attrs['UNIT'] = value
                    vprint(f'add /{key:<{max_digit}} attribute: UNIT = {value}')

    vprint(f'close  HDF5 file: {fname}')

    return fname


def write_hdf5_block(fname, data, datasetName, block=None, mode='a', print_msg=True):
    """Write data to existing HDF5 dataset in disk block by block.
    Parameters: data        - np.ndarray 1/2/3/4D matrix
                datasetName - str, dataset name
                block       - list of 2/4/6/8 int, for
                              [d1Start, d1End,
                               d2Start, d2End,
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
        elif len(shape) == 4:
            block = [0, shape[0],
                     0, shape[1],
                     0, shape[2],
                     0, shape[3]]

    # write
    if print_msg:
        print('-'*50)
        print(f'open  HDF5 file {fname} in {mode} mode')
        print(f"writing dataset /{datasetName:<25} block: {block}")
    with h5py.File(fname, mode) as f:
        if len(block) == 8:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3],
                           block[4]:block[5],
                           block[6]:block[7]] = data

        elif len(block) == 6:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3],
                           block[4]:block[5]] = data

        elif len(block) == 4:
            f[datasetName][block[0]:block[1],
                           block[2]:block[3]] = data

        elif len(block) == 2:
            f[datasetName][block[0]:block[1]] = data

    if print_msg:
        print(f'close HDF5 file {fname}.')

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
    vprint = print if print_msg else lambda *args, **kwargs: None

    if isinstance(datasetNames, str):
        datasetNames = list(datasetNames)
    vprint(f'delete {datasetNames} from file {fname}')

    # 1. rename the file to a temporary file
    temp_file = os.path.join(os.path.dirname(fname), f'tmp_{os.path.basename(fname)}')
    shutil.move(fname, temp_file)
    vprint(f'move {fname} to {temp_file}')

    # 2. write a new file with all data except for the one to be deleted
    vprint(f'read   HDF5 file: {temp_file} with r mode')
    vprint(f'create HDF5 file: {fname} with w mode')
    with h5py.File(temp_file, 'r') as fi:
        with h5py.File(fname, 'w') as fo:

            # datasets
            compression = None
            maxDigit = max(len(i) for i in list(fi.keys()))
            for dsName in [i for i in fi.keys() if i not in datasetNames]:
                ds = fi[dsName]
                msg = f'create dataset /{dsName:<{maxDigit}} of {str(ds.dtype):<10}'
                msg += f' in size of {str(ds.shape):<20} with compression={compression}'
                vprint(msg)
                fo.create_dataset(dsName, data=ds[:], chunks=True, compression=compression)

            # metadata
            for key, value in fi.attrs.items():
                fo.attrs[key] = str(value)

    vprint(f'finished writing to {fname}')
    vprint(f'old file is now saved as: {temp_file}. Use rm command to delete it.')

    return fname



#########################################################################

def write_roipac_rsc(metadata, out_file, update_mode=False, print_msg=False):
    """Write attribute dict into ROI_PAC .rsc file
    Inputs:
        metadata : dict, attributes dictionary
        out_file : rsc file name, to which attribute is written
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
            print(f'write file: {out_file}')
        maxDigit = max([len(key) for key in metadata.keys()]+[2])
        with open(out_file, 'w') as f:
            for key in sorted(metadata.keys()):
                f.write('{k:<{d}}    {v}\n'.format(k=str(key),
                                                   d=maxDigit,
                                                   v=str(metadata[key])))
    return out_file


def write_gdal_vrt(meta, out_file):
    """Write GDAL VRT file.

    !!! This function is NOT RIGHT. DO NOT USE IT. Keep here as a placeholder ONLY. !!!
    It needs more work.

    Parameters: meta     - dict, dictionary of metadata
                out_file - str, VRT file name to which attributes are written
    """
    # data type: mintpy to gdal
    dtype_dict = {
        'int8'      : 'Byte',
        'int16'     : 'Int16',
        'float32'   : 'Float32',
        'float64'   : 'Float64',
        'complex64' : 'CFloat32',
        'complex128': 'CFloat64',
    }

    # pixel / line / image offset
    pixel_offset_dict = {
        'int8'      : '2',
        'int16'     : '4',
        'float32'   : '8',
        'float64'   : '16',
        'complex64' : '16',
        'complex128': '32',
    }
    pixel_offset = int(pixel_offset_dict[meta['DATA_TYPE']])
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    num_band = int(meta['BANDS'])

    interleave = meta['INTERLEAVE']
    if interleave == 'BIP':
        line_offset  = pixel_offset * num_band * width
        image_offset = pixel_offset
    elif interleave == 'BIL':
        line_offset  = pixel_offset * width * num_band
        image_offset = pixel_offset * width
    elif interleave == 'BSQ':
        line_offset  = pixel_offset * width
        image_offset = pixel_offset * width * length
    else:
        raise ValueError(f'un-recognized band interleave type: {interleave}')

    # compose VRT file string
    ds_str = '<VRTDataset rasterXSize="{w}" rasterYSize="{l}">\n'.format(w=meta['WIDTH'], l=meta['LENGTH'])
    for band in range(num_band):
        band_str = '''<VRTRasterBand dataType="{d}" band="{b}" subClass="VRTRawRasterBand">
        <SourceFilename relativeToVRT="1">{f}</SourceFilename>
        <ByteOrder>LSB</ByteOrder>
        <ImageOffset>{io}</ImageOffset>
        <PixelOffset>{po}</PixelOffset>
        <LineOffset>{lo}</LineOffset>
    </VRTRasterBand>
        '''.format(
            d=dtype_dict[meta['DATA_TYPE']],
            b=band+1,
            f=os.path.basename(out_file[:-4]),
            io=image_offset * band,
            po=pixel_offset,
            lo=line_offset,
        )
        ds_str += band_str
    ds_str += '</VRTDataset>\n'

    # write VRT file
    with open(out_file, 'w') as f:
        f.write(ds_str)

    return out_file


def write_isce_xml(meta, fname, print_msg=True):
    """Write XML metadata file in ISCE-2 format

    Parameters: meta      - dict, attributes dictionary
                fname     - str, path of data file, not the metadata file
                print_msg - bool, print out message
    Examples:   write_isce_xml(atr, fname='filt_fine.cor')
    """

    import isce
    import isceobj

    # data type
    dtype = readfile.DATA_TYPE_NUMPY2ISCE[meta['DATA_TYPE']]

    # write ISCE XML and GDAL VRT files
    image_type = meta['FILE_TYPE']
    if not image_type:          img = isceobj.Image.createImage()
    elif image_type == '.slc':  img = isceobj.Image.createSlcImage()
    elif image_type == '.unw':  img = isceobj.Image.createUnwImage()
    elif image_type == '.int':  img = isceobj.Image.createIntImage()
    else:                       img = isceobj.Image.createImage()

    img.setFilename(fname)
    img.setWidth(int(meta['WIDTH']))
    img.setLength(int(meta['LENGTH']))
    img.setAccessMode('READ')
    img.bands = int(meta.get('BANDS', '1'))
    img.dataType = dtype
    img.scheme = meta.get('INTERLEAVE', 'BIL')
    img.renderHdr()
    img.renderVRT()
    if print_msg:
        print(f'write file: {fname}.xml')
        print(f'write file: {fname}.vrt')

    return


def write_isce_file(data, out_file, file_type='isce_unw'):
    """write data to file in ISCE format

    Parameters: data      - 2D np.ndarray, binary data matrix
                out_file  - str, path of output binary data file
                file_type - str, file type
    Returns:    out_file  - str, path of output binary data file
    """
    # fix potential typo
    file_type = file_type.replace('-', '_')

    # write data to binary file
    data.tofile(out_file)

    # write isce xml metadata file
    length, width = data.shape
    meta = {
        'LENGTH' : length,
        'WIDTH'  : width,
    }

    if file_type == 'isce_unw':
        meta['FILE_TYPE'] = '.unw'
        meta['BANDS'] = 2
        meta['DATA_TYPE'] = 'float32'
        meta['INTERLEAVE'] = 'BIL'
        meta['WIDTH'] = int(width / 2)

    elif file_type == 'isce_int':
        meta['FILE_TYPE'] = '.int'
        meta['BANDS'] = 1
        meta['DATA_TYPE'] = 'complex64'
        meta['INTERLEAVE'] = 'BIL'

    elif file_type == 'isce_cor':
        meta['FILE_TYPE'] = '.cor'
        meta['BANDS'] = 1
        meta['DATA_TYPE'] = 'float32'
        meta['INTERLEAVE'] = 'BIL'

    elif file_type == 'isce_slc':
        meta['FILE_TYPE'] = '.slc'
        meta['BANDS'] = 1
        meta['DATA_TYPE'] = 'complex64'
        meta['INTERLEAVE'] = 'BIP'

    else:
        raise ValueError(f'un-recognized ISCE file type: {file_type}')

    write_isce_xml(meta, out_file)

    return out_file



#########################################################################

def write_binary(data_list, out_file, data_type=None, interleave='BIL'):
    """Write binary file.
    Parameters: data_list  - list of 2D np.ndarray matrices
                out_file   - str, path of the output binary file
                data_type  - str/np.dtype, numpy data type object
                interleave - str, band interleave type, BSQ, BIL, BIP
    Returns:    out_file   - str, path of the output binary file
    """
    # data type
    if data_type:
        data_type = data_type.lower()
        if data_type != data_list[0].dtype:
            data_list = [np.array(x, dtype=data_type) for x in data_list]

    # stack multi-band - reshape
    interleave = interleave.upper()
    if interleave == 'BIP':
        data_list = [x.reshape(-1,1) for x in data_list]
    elif interleave == 'BSQ':
        data_list = [x.reshape(1,-1) for x in data_list]

    # stack multi-band - stack
    data = data_list[0]
    for datai in data_list[1:]:
        data = np.hstack((data, datai))

    # write to file
    data.tofile(out_file)

    return out_file



######################## Obsolete functions #############################

def write_float32(*args):
    """Write ROI_PAC rmg format with float32 precision (BIL)
    Format of the binary file is same as roi_pac unw, cor, or hgt data.
          should rename to write_rmg_float32()

    Example:
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


def write_complex_float32(data, out_file):
    """write complex float32 data into file"""
    data = np.array(data, dtype=np.complex64)
    data.tofile(out_file)
    return out_file


#def write_complex_float32(data, out_file):
#    """Writes roi_pac .int data"""
#    num_pixel = data.size
#    F = np.zeros([2 * num_pixel, 1], np.float32)
#    id1 = list(range(0, 2 * num_pixel, 2))
#    id2 = list(range(1, 2 * num_pixel, 2))
#    F[id1] = np.reshape(np.cos(data), (num_pixel, 1))
#    F[id2] = np.reshape(np.sin(data), (num_pixel, 1))
#    F.tofile(out_file)
#    return out_file


def write_complex_int16(data, out_file):
    """Write gamma scomplex data, i.e. .slc file.
        data is complex 2-D matrix
        real, imagery, real, ...
    Write in this way, because numpy does not have complex int16 directly.
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

#########################################################################
