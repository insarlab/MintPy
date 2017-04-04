############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Sep 2015: Add write_gamma_float() and write_gamma_scomplex()
# Yunjun, Oct 2015: Add support for write_float32(amp, phase, outname)
# Yunjun, Jan 2016: Add write()


import os

import h5py
import numpy as np
from PIL import Image


def write(*args):
    '''Write one dataset, i.e. interferogram, coherence, velocity, dem ...
        Return 0 if failed.
  
    Usage:
        write(data,atr,outname)
        write(rg,az,atr,outname)
    
    Inputs:
        data : 2D data matrix
        atr  : attribute object
        outname : output file name
    
    Output:
        output file name
    
    Examples:
        write(data,atr,'velocity.h5')
        write(data,atr,'temporal_coherence.h5')
        write(data,atr,'100120-110214.unw')
        write(data,atr,'strm1.dem')
        write(data,atr,'100120.mli')
        write(rg,az,atr,'geomap_4lks.trans')
    '''

    ########## Check Inputs ##########
    if len(args) == 4:       ## .trans file
        rg      = args[0]
        az      = args[1]
        atr     = args[2]
        outname = args[3]
    else:
        data    = args[0]
        atr     = args[1]
        outname = args[2]
  
    ext = os.path.splitext(outname)[1].lower()
    ############### Read ###############
    #print 'writing >>> '+outname
    ##### PySAR HDF5 product
    if ext in ['.h5','.he5']:
        k = atr['FILE_TYPE']
        if k in ['interferograms','coherence','wrapped','timeseries']:
            print 'Un-supported file type: '+k
            print 'Only support 1-dataset-1-attribute file, i.e. velocity, mask, ...'
            return 0;
        h5file = h5py.File(outname,'w')
        group = h5file.create_group(k)
        dset = group.create_dataset(k, data=data, compression='gzip')
        for key , value in atr.iteritems():
            group.attrs[key]=value
        h5file.close()
  
        return outname

    ##### ISCE / ROI_PAC GAMMA / Image product
    else:
        ##### Write Data File
        if   ext in ['.unw','.cor','.hgt']:
            write_float32(data,outname)
        elif ext == '.dem':
            write_real_int16(data,outname)
        elif ext == '.trans':
            write_float32(rg,az,outname)
        elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']:
            data.save(outname)
        elif ext == '.mli':
            write_real_float32(data,outname)
        elif ext == '.slc':
            write_complex_int16(data,outname)
        elif ext == '.int':
            write_complex64(data, outname)
        else: print 'Un-supported file type: '+ext; return 0;
  
        ##### Write .rsc File
        write_roipac_rsc(atr, outname+'.rsc')
        
        return outname


def write_roipac_rsc(atr, outname, sorting=True):
    '''Write attribute dict into ROI_PAC .rsc file
    Inputs:
        atr     - dict, attributes dictionary
        outname - rsc file name, to which attribute is writen
        sorting - bool, sort attributes in alphabetic order while writing
    Output:
        outname
    '''
    keyList = atr.iterkeys()
    if sorting:
        keyList = sorted(keyList)
    digits = max([len(key) for key in keyList]+[2])
    f = '{0:<%d}    {1}'%(digits)
    frsc = open(outname,'w')
    for key in keyList:
        frsc.write(f.format(str(key), str(atr[key]))+'\n')
    frsc.close()
    return outname


def write_float32(*args):
    '''Write ROI_PAC rmg format with float32 precision
    Format of the binary file is same as roi_pac unw, cor, or hgt data.
          should rename to write_rmg_float32()
    
    Exmaple:
            write_float32(phase, outname)
            write_float32(amp, phase, outname)
    '''
 
    if len(args)==2:
        amp     = args[0]
        pha     = args[0]
        outname = args[1]
    elif len(args)==3:
        amp     = args[0]
        pha     = args[1]
        outname = args[2]
    else:
        print 'Error while getting args: support 2/3 args only.'
        return
 
    nlines = pha.shape[0]
    WIDTH  = pha.shape[1]
    F=np.zeros([2*nlines*WIDTH,1],np.float32)
 
    for line in range(nlines):
        F[(2*WIDTH)*(line) :       (2*WIDTH)*(line)+WIDTH]=np.reshape(amp[line][:],[WIDTH,1])
        F[(2*WIDTH)*(line)+WIDTH : (2*WIDTH)*(line+1)]    =np.reshape(pha[line][:],[WIDTH,1])
 
    F.tofile(outname)
    return outname


def write_complex64(data,outname):
    '''Writes roi_pac .int data'''
    nlines=data.shape[0]
    WIDTH=data.shape[1]
    R=np.cos(data)
    Im=np.sin(data)
    # F=np.zeros([2*nlines*WIDTH,1],np.complex64) 
    F=np.zeros([2*nlines*WIDTH,1],np.float32)  
    id1=range(0,2*nlines*WIDTH,2)
    id2=range(1,2*nlines*WIDTH,2)
    F[id1]=np.reshape(R,(nlines*WIDTH,1))
    F[id2]=np.reshape(Im,(nlines*WIDTH,1))
    F.tofile(outname)
    return outname


def write_real_int16(data,outname):
    data=np.array(data,dtype=np.int16)
    data.tofile(outname)
    return outname


def write_dem(data,outname):
    data=np.array(data,dtype=np.int16)
    data.tofile(outname)
    return outname


def write_real_float32(data,outname):
    '''write gamma float data, i.e. .mli file.'''
    data=np.array(data,dtype=np.float32)
    data.tofile(outname)
    return outname


def write_complex_int16(data,outname):
    '''Write gamma scomplex data, i.e. .slc file.
        data is complex 2-D matrix
        real, imagery, real, ...
    '''

    nlines = data.shape[0]
    WIDTH  = data.shape[1]
    id1 = range(0,2*nlines*WIDTH,2)
    id2 = range(1,2*nlines*WIDTH,2)

    F=np.zeros([2*nlines*WIDTH,1],np.int16)
    F[id1]=np.reshape(np.array(data.real,np.int16),(nlines*WIDTH,1))
    F[id2]=np.reshape(np.array(data.imag,np.int16),(nlines*WIDTH,1))
    F.tofile(outname)
    return outname

