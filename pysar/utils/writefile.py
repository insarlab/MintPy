############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################
#
# Recommend usage:
#   import pysar.utils.writefile as writefile
#


import os
import sys
import h5py
import numpy as np
#from PIL import Image
from pysar.objects import timeseries


def write(data, metadata, outfile, infile=None):
    '''Write one dataset, i.e. interferogram, coherence, velocity, dem ...
        Return 0 if failed.

    Usage:
        write(data,metadata,outfile)

    Inputs:
        data : 2D data matrix
        metadata  : attribute object
        outfile : output file name
    
    Output:
        output file name

    Examples:
        write(data,metadata,'velocity.h5')
        write(data,metadata,'temporal_coherence.h5')
        write(data,metadata,'100120-110214.unw')
        write(data,metadata,'strm1.dem')
        write(data,metadata,'100120.mli')
        write(rg,az,metadata,'geomap_4lks.trans')
    '''
    ext = os.path.splitext(outfile)[1].lower()

    ##### PySAR HDF5 product
    if ext in ['.h5','.he5']:
        k = metadata['FILE_TYPE']
        if k == 'timeseries':
            if infile is None:
                print('ERROR: can not write {} file with infile as reference file!'.format(k))
                sys.exit(-1)
            obj = timeseries(outfile)
            obj.write2hdf5(data, metadata=metadata, refFile=infile)
        else:
            print('create HDF5 file: {} with w mode'.format(outfile))
            f = h5py.File(outfile, 'w')

            print('create dataset /{} of {:<10} in size of {}'.format(k, str(data.dtype), data.shape))
            ds = f.create_dataset(k, data=data, chunks=True, compression="gzip")

            for key, value in metadata.items():
                f.attrs[key] = str(value)
            f.close()
            print('finished writing to {}'.format(outfile))

    ##### ISCE / ROI_PAC GAMMA / Image product
    else:
        ##### Write Data File
        if   ext in ['.unw','.cor','.hgt']:
            write_float32(data,outfile)
        elif ext == '.dem':
            write_real_int16(data,outfile)
        elif ext in ['.trans']:
            write_float32(rg,az,outfile)
        elif ext in ['.utm_to_rdc','.UTM_TO_RDC']:
            data = np.zeros(rg.shape, dtype=np.complex64)
            data.real = rg
            data.imag = az
            data.astype('>c8').tofile(outfile)
        #elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
        #    data.save(outfile)
        elif ext == '.mli':
            write_real_float32(data,outfile)
        elif ext == '.slc':
            write_complex_int16(data,outfile)
        elif ext == '.int':
            write_complex64(data, outfile)
        elif metadata['DATA_TYPE'].lower() in ['float32','float']:
            write_real_float32(data,outfile)
        elif metadata['DATA_TYPE'].lower() in ['int16','short']:
            write_real_int16(data,outfile)
        else: print('Un-supported file type: '+ext); return 0;

        ##### Write .rsc File
        write_roipac_rsc(metadata, outfile+'.rsc')
        return outfile


def write_roipac_rsc(metadata, outfile, sorting=True):
    '''Write attribute dict into ROI_PAC .rsc file
    Inputs:
        metadata     - dict, attributes dictionary
        outfile - rsc file name, to which attribute is writen
        sorting - bool, sort attributes in alphabetic order while writing
    Output:
        outfile
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
    frsc = open(outfile,'w')
    for key in dictKey:
        frsc.write(f.format(str(key), str(metadata[key]))+'\n')
    frsc.close()
    return outfile


def write_float32(*args):
    '''Write ROI_PAC rmg format with float32 precision
    Format of the binary file is same as roi_pac unw, cor, or hgt data.
          should rename to write_rmg_float32()
    
    Exmaple:
            write_float32(phase, outfile)
            write_float32(amp, phase, outfile)
    '''
 
    if len(args)==2:
        amp     = args[0]
        pha     = args[0]
        outfile = args[1]
    elif len(args)==3:
        amp     = args[0]
        pha     = args[1]
        outfile = args[2]
    else:
        print('Error while getting args: support 2/3 args only.')
        return
 
    nlines = pha.shape[0]
    WIDTH  = pha.shape[1]
    F=np.zeros([2*nlines*WIDTH,1],np.float32)

    for line in range(nlines):
        F[(2*WIDTH)*(line) :       (2*WIDTH)*(line)+WIDTH]=np.reshape(amp[line][:],[WIDTH,1])
        F[(2*WIDTH)*(line)+WIDTH : (2*WIDTH)*(line+1)]    =np.reshape(pha[line][:],[WIDTH,1])
 
    F.tofile(outfile)
    return outfile


def write_complex64(data,outfile):
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
    F.tofile(outfile)
    return outfile


def write_real_int16(data,outfile):
    data=np.array(data,dtype=np.int16)
    data.tofile(outfile)
    return outfile


def write_dem(data,outfile):
    data=np.array(data,dtype=np.int16)
    data.tofile(outfile)
    return outfile


def write_real_float32(data,outfile):
    '''write gamma float data, i.e. .mli file.'''
    data=np.array(data,dtype=np.float32)
    data.tofile(outfile)
    return outfile


def write_complex_int16(data,outfile):
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
    F.tofile(outfile)
    return outfile

