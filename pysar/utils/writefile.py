############################################################
# Program is part of PySAR v2.0                            #
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


def write(datasetDict, outFile, metadata=None, refFile=None, compression=None):
    """ Write one file.
    Parameters: datasetDict : dict of dataset, with key = datasetName and value = 2D/3D array, e.g.:
                    {'height'        : np.ones((   200,300), dtype=np.int16),
                     'incidenceAngle': np.ones((   200,300), dtype=np.float32),
                     'bperp'         : np.ones((80,200,300), dtype=np.float32),
                     ...}
                outFile : str, output file name
                metadata : dict of attributes
                refFile : str, reference file to get auxliary info
                compression : str, compression while writing to HDF5 file, None, "lzf", "gzip"
    Returns:    outFile : str
    Examples:   dsDict = dict()
                dsDict['velocity'] = np.ones((200,300), dtype=np.float32)
                write(datasetDict=dsDict, outFile='velocity.h5', metadata=atr)
    """
    ext = os.path.splitext(outFile)[1].lower()
    if refFile and metadata is None:
        metadata = readfile.read_attribute(refFile)

    if type(datasetDict) is np.ndarray:
        data = np.array(datasetDict)
        datasetDict = dict()
        datasetDict[metadata['FILE_TYPE']] = data

    # HDF5 File
    if ext in ['.h5','.he5']:
        k = metadata['FILE_TYPE']
        if k == 'timeseries':
            if refFile is None:
                print('ERROR: can not write {} file without reference file!'.format(k))
                sys.exit(-1)
            obj = timeseries(outFile)
            obj.write2hdf5(datasetDict[k], metadata=metadata, refFile=refFile)

        else:
            print('create HDF5 file: {} with w mode'.format(outFile))
            f = h5py.File(outFile, 'w')

            # Write input datasets
            maxDigit = max([len(i) for i in list(datasetDict.keys())])
            for dsName in datasetDict.keys():
                data = datasetDict[dsName]
                print('create dataset /{d:<{w}} of {t:<10} in size of {s}'.format(d=dsName,
                                                                                  w=maxDigit,
                                                                                  t=str(data.dtype),
                                                                                  s=data.shape))
                ds = f.create_dataset(dsName, data=data, chunks=True, compression=compression)
                if dsName == 'velocity':
                    ds.attrs['MaxValue'] = np.nanmax(data)  #facilitate disp_min/max for mutiple subplots in view.py
                    ds.attrs['MinValue'] = np.nanmin(data)  #facilitate disp_min/max for mutiple subplots in view.py

            # Write extra/auxliary datasets from refFile
            if refFile:
                fr = h5py.File(refFile, 'r')
                dsNames = [i for i in fr.keys() if i not in list(datasetDict.keys())]
                for dsName in dsNames:
                    ds = fr[dsName]
                    print('create dataset /{d:<{w}} of {t:<10} in size of {s}'.format(d=dsName,
                                                                                      w=maxDigit,
                                                                                      t=str(ds.dtype),
                                                                                      s=ds.shape))
                    f.create_dataset(dsName, data=ds[:], chunks=True, compression=compression)
                fr.close()

            # metadata
            for key, value in metadata.items():
                f.attrs[key] = str(value)
            f.close()
            print('finished writing to {}'.format(outFile))

    ##### ISCE / ROI_PAC GAMMA / Image product
    else:
        ##### Write Data File
        if   ext in ['.unw','.cor','.hgt']:
            write_float32(data,outFile)
        elif ext == '.dem':
            write_real_int16(data,outFile)
        elif ext in ['.trans']:
            write_float32(rg,az,outFile)
        elif ext in ['.utm_to_rdc','.UTM_TO_RDC']:
            data = np.zeros(rg.shape, dtype=np.complex64)
            data.real = rg
            data.imag = az
            data.astype('>c8').tofile(outFile)
        #elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        #    data.save(outFile)
        elif ext == '.mli':
            write_real_float32(data,outFile)
        elif ext == '.slc':
            write_complex_int16(data,outFile)
        elif ext == '.int':
            write_complex64(data, outFile)
        elif metadata['DATA_TYPE'].lower() in ['float32','float']:
            write_real_float32(data,outFile)
        elif metadata['DATA_TYPE'].lower() in ['int16','short']:
            write_real_int16(data,outFile)
        else: print('Un-supported file type: '+ext); return 0;

        ##### Write .rsc File
        write_roipac_rsc(metadata, outFile+'.rsc')
        return outFile


def write_roipac_rsc(metadata, outFile, sorting=True):
    '''Write attribute dict into ROI_PAC .rsc file
    Inputs:
        metadata     - dict, attributes dictionary
        outFile - rsc file name, to which attribute is writen
        sorting - bool, sort attributes in alphabetic order while writing
    Output:
        outFile
    '''
    # Convert PYSAR attributes to ROI_PAC attributes
    metadata['FILE_LENGTH'] = metadata['LENGTH']

    # Convert 3.333e-4 to 0.0003333
    if 'X_STEP' in metadata.keys():
        metadata['X_STEP'] = str(float(metadata['X_STEP']))
        metadata['Y_STEP'] = str(float(metadata['Y_STEP']))
        metadata['X_FIRST'] = str(float(metadata['X_FIRST']))
        metadata['Y_FIRST'] = str(float(metadata['Y_FIRST']))

    # max digit for space formating
    digits = max([len(key) for key in metadata.keys()]+[2])
    f = '{0:<%d}    {1}'%(digits)

    # sorting by key name
    dictKey = metadata.keys()
    if sorting:
        dictKey = sorted(dictKey)

    # writing .rsc file
    frsc = open(outFile,'w')
    for key in dictKey:
        frsc.write(f.format(str(key), str(metadata[key]))+'\n')
    frsc.close()
    return outFile


def write_float32(*args):
    '''Write ROI_PAC rmg format with float32 precision
    Format of the binary file is same as roi_pac unw, cor, or hgt data.
          should rename to write_rmg_float32()
    
    Exmaple:
            write_float32(phase, outFile)
            write_float32(amp, phase, outFile)
    '''
 
    if len(args)==2:
        amp     = args[0]
        pha     = args[0]
        outFile = args[1]
    elif len(args)==3:
        amp     = args[0]
        pha     = args[1]
        outFile = args[2]
    else:
        print('Error while getting args: support 2/3 args only.')
        return
 
    nlines = pha.shape[0]
    WIDTH  = pha.shape[1]
    F=np.zeros([2*nlines*WIDTH,1],np.float32)

    for line in range(nlines):
        F[(2*WIDTH)*(line) :       (2*WIDTH)*(line)+WIDTH]=np.reshape(amp[line][:],[WIDTH,1])
        F[(2*WIDTH)*(line)+WIDTH : (2*WIDTH)*(line+1)]    =np.reshape(pha[line][:],[WIDTH,1])
 
    F.tofile(outFile)
    return outFile


def write_complex64(data,outFile):
    '''Writes roi_pac .int data'''
    nlines=data.shape[0]
    WIDTH=data.shape[1]
    R=np.cos(data)
    Im=np.sin(data)
    # F=np.zeros([2*nlines*WIDTH,1],np.complex64) 
    F=np.zeros([2*nlines*WIDTH,1],np.float32)  
    id1=list(range(0,2*nlines*WIDTH,2))
    id2=list(range(1,2*nlines*WIDTH,2))
    F[id1]=np.reshape(R,(nlines*WIDTH,1))
    F[id2]=np.reshape(Im,(nlines*WIDTH,1))
    F.tofile(outFile)
    return outFile


def write_real_int16(data,outFile):
    data=np.array(data,dtype=np.int16)
    data.tofile(outFile)
    return outFile


def write_dem(data,outFile):
    data=np.array(data,dtype=np.int16)
    data.tofile(outFile)
    return outFile


def write_real_float32(data,outFile):
    '''write gamma float data, i.e. .mli file.'''
    data=np.array(data,dtype=np.float32)
    data.tofile(outFile)
    return outFile


def write_complex_int16(data,outFile):
    '''Write gamma scomplex data, i.e. .slc file.
        data is complex 2-D matrix
        real, imagery, real, ...
    '''

    nlines = data.shape[0]
    WIDTH  = data.shape[1]
    id1 = list(range(0,2*nlines*WIDTH,2))
    id2 = list(range(1,2*nlines*WIDTH,2))

    F=np.zeros([2*nlines*WIDTH,1],np.int16)
    F[id1]=np.reshape(np.array(data.real,np.int16),(nlines*WIDTH,1))
    F[id2]=np.reshape(np.array(data.imag,np.int16),(nlines*WIDTH,1))
    F.tofile(outFile)
    return outFile

