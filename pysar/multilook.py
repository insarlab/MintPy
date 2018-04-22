#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys, re
import argparse
import warnings
import h5py
import numpy as np
from pysar.utils import readfile, writefile, datetime as ptime, utils as ut


######################################## Sub Functions ############################################
def multilook_matrix(matrix,lks_y,lks_x):
    rows,cols = matrix.shape
    lks_x = int(lks_x)
    lks_y = int(lks_y)
    if lks_x == 1 and lks_y == 1:
        return matrix

    rows_mli=int(np.floor(rows/lks_y))
    cols_mli=int(np.floor(cols/lks_x))
    #thr = np.floor(lks_x*lks_y/2)
    matrix_Cmli=np.zeros((rows,    cols_mli))
    matrix_mli =np.zeros((rows_mli,cols_mli))

    #for c in range(lks_x):   matrix_Cmli = matrix_Cmli + matrix[:,range(c,cols_mli*lks_x,lks_x)]
    #for r in range(lks_y):   matrix_mli  = matrix_mli  + matrix_Cmli[  range(r,rows_mli*lks_y,lks_y),:]
    #for c in range(int(cols_mli)):  matrix_Cmli[:,c]=np.nansum(matrix[:,(c)*lks_x:(c+1)*lks_x],1)
    #for r in range(int(rows_mli)):  matrix_mli[r,:] =np.nansum(matrix_Cmli[(r)*lks_y:(r+1)*lks_y,:],0)
    #matrix_mli=matrix_mli/(lks_y*lks_x)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for c in range(cols_mli):  matrix_Cmli[:,c] = np.nanmean(matrix[:,(c)*lks_x:(c+1)*lks_x],1)
        for r in range(rows_mli):  matrix_mli[r,:]  = np.nanmean(matrix_Cmli[(r)*lks_y:(r+1)*lks_y,:],0)
    del matrix, matrix_Cmli
    return matrix_mli


def multilook_data(data, lksY, lksX):
    """Modified from Praveen on StackOverflow:
    https://stackoverflow.com/questions/34689519/how-to-coarser-the-2-d-array-data-resolution
    Parameters: data : 2D / 3D np.array
                lksY : int, number of multilook in y/azimuth direction
                lksX : int, number of multilook in x/range direction
    Returns:    coarseData : 2D / 3D np.array after multilooking in last two dimension
    """
    shape = np.array(data.shape, dtype=float)
    if len(shape) == 2:
        newShape = np.floor(shape / (lksY, lksX)).astype(int) * (lksY, lksX)
        cropData = data[:newShape[0], :newShape[1]]
        temp = cropData.reshape((newShape[0] // lksY, lksY, newShape[1] // lksX, lksX))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            coarseData = np.nanmean(temp, axis=(1,3))
    elif len(shape) == 3:
        newShape = np.floor(shape / (1, lksY, lksX)).astype(int) * (1, lksY, lksX)
        cropData = data[:newShape[0], :newShape[1], :newShape[2]]
        temp = cropData.reshape((newShape[0], newShape[1] // lksY, lksY, newShape[2] // lksX, lksX))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            coarseData = np.nanmean(temp, axis=(2,4))
    return coarseData


def multilook_attribute(atr_dict,lks_y,lks_x, print_msg=True):
    #####
    atr = dict()
    for key, value in iter(atr_dict.items()):
        atr[key] = str(value)
  
    ##### calculate new data size
    length = int(atr['LENGTH'])
    width  = int(atr['WIDTH'])
    length_mli = int(np.floor(length/lks_y))
    width_mli  = int(np.floor(width/lks_x))
  
    ##### Update attributes
    atr['LENGTH'] = str(length_mli)
    atr['WIDTH'] = str(width_mli)
    atr['XMIN'] = '0'
    atr['YMIN'] = '0'
    atr['XMAX'] = str(width_mli-1)
    atr['YMAX'] = str(length_mli-1)
    if print_msg:
        print('update LENGTH, WIDTH, YMIN, YMAX, XMIN, XMAX')
    
    try:
        atr['Y_STEP'] = str(lks_y*float(atr['Y_STEP']))
        atr['X_STEP'] = str(lks_x*float(atr['X_STEP']))
        if print_msg: print('update Y/X_STEP')
    except: pass
    try:
        atr['AZIMUTH_PIXEL_SIZE'] = str(lks_y*float(atr['AZIMUTH_PIXEL_SIZE']))
        atr['RANGE_PIXEL_SIZE']   = str(lks_x*float(atr['RANGE_PIXEL_SIZE']))
        if print_msg: print('update AZIMUTH/RANGE_PIXEL_SIZE')
    except: pass

    if not 'Y_FIRST' in atr.keys():
        try:
            atr['RLOOKS'] = str(int(atr['RLOOKS'])*lks_x)
            atr['ALOOKS'] = str(int(atr['ALOOKS'])*lks_y)
            if print_msg: print('update R/ALOOKS')
        except: pass

    try:
        atr['REF_Y'] = str(int(int(atr['REF_Y'])/lks_y))
        atr['REF_X'] = str(int(int(atr['REF_X'])/lks_x))
        if print_msg: print('update ref_y/x')
    except: pass
    try:
        atr['SUBSET_YMIN'] = str(int(int(atr['SUBSET_YMIN'])/lks_y))
        atr['SUBSET_YMAX'] = str(int(int(atr['SUBSET_YMAX'])/lks_y))
        atr['SUBSET_XMIN'] = str(int(int(atr['SUBSET_XMIN'])/lks_x))
        atr['SUBSET_XMAX'] = str(int(int(atr['SUBSET_XMAX'])/lks_x))
        if print_msg: print('update subset_y0/y1/x0/x1')
    except: pass

    return atr


def multilook_file(infile,lks_y,lks_x,outfile=None):
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    ## input file info
    atr = readfile.read_attribute(infile)    
    k = atr['FILE_TYPE']
    print('multilooking {} {} file: {}'.format(atr['PROCESSOR'], k, infile))
    print('number of looks in y / azimuth direction: %d' % lks_y)
    print('number of looks in x / range   direction: %d' % lks_x)

    ## output file name
    if not outfile:
        if os.getcwd() == os.path.dirname(os.path.abspath(infile)):
            ext = os.path.splitext(infile)[1]
            outfile = os.path.splitext(infile)[0]+'_'+str(lks_y)+'alks_'+str(lks_x)+'rlks'+ext
        else:
            outfile = os.path.basename(infile)
    print('writing >>> '+outfile)

    ###############################################################################
    ## Read/Write multi-dataset files
    if k in ['interferograms','coherence','wrapped','timeseries']:
        h5 = h5py.File(infile,'r')
        epochList = sorted(h5[k].keys())
        epoch_num = len(epochList)
        prog_bar = ptime.progressBar(maxValue=epoch_num)

        h5out = h5py.File(outfile,'w')
        group = h5out.create_group(k)

        if k in ['interferograms','coherence','wrapped']:
            date12_list = ptime.list_ifgram2date12(epochList)
            print('number of interferograms: '+str(len(epochList)))
            for i in range(epoch_num):
                epoch = epochList[i]
                data = h5[k][epoch].get(epoch)[:]
                atr = h5[k][epoch].attrs

                data_mli = multilook_data(data,lks_y,lks_x)
                atr_mli = multilook_attribute(atr,lks_y,lks_x,print_msg=False)

                gg = group.create_group(epoch)
                dset = gg.create_dataset(epoch, data=data_mli)
                for key, value in iter(atr_mli.items()):
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])

        elif k == 'timeseries':
            print('number of acquisitions: '+str(len(epochList)))
            for i in range(epoch_num):
                epoch = epochList[i]
                data = h5[k].get(epoch)[:]

                data_mli = multilook_data(data,lks_y,lks_x)
                
                dset = group.create_dataset(epoch, data=data_mli)
                prog_bar.update(i+1, suffix=epoch)
            atr = h5[k].attrs
            atr_mli = multilook_attribute(atr,lks_y,lks_x)
            for key, value in iter(atr_mli.items()):
                group.attrs[key] = value

        h5.close()
        h5out.close()
        prog_bar.close()

    ## Read/Write single-dataset files
    elif k in ['.trans','.utm_to_rdc','.UTM_TO_RDC']:        
        rg,az,atr = readfile.read(infile)
        rgmli = multilook_data(rg,lks_y,lks_x); #rgmli *= 1.0/lks_x
        azmli = multilook_data(az,lks_y,lks_x); #azmli *= 1.0/lks_y
        atr = multilook_attribute(atr,lks_y,lks_x)
        dsDict = {}
        dsDict['rangeCoord'] = rgmli
        dsDict['azimuthCoord'] = azmli
        writefile.write(dsDict, out_file=outfile, metadata=atr)
    else:
        data,atr = readfile.read(infile)
        if lks_y < 0 and lks_x < 0:
            lks_x = 1./abs(lks_x)
            lks_y = 1./abs(lks_y)
            outShape = (int(int(atr['LENGTH'])/lks_y), int(int(atr['WIDTH'])/lks_x))
            data_mli = ut.interpolate_data(data, outShape=outShape, interpMethod='linear')
        else:
            data_mli = multilook_data(data, lks_y, lks_x)
        atr = multilook_attribute(atr,lks_y,lks_x)
        writefile.write(data_mli, out_file=outfile, metadata=atr)

    return outfile


##################################################################################################
EXAMPLE="""example:
  multilook.py  velocity.h5  15 15
  multilook.py  srtm30m.dem  10 10  -o srtm30m_300m.dem

  To interpolate input file into larger size file:
  multilook.py  bperp.rdr  -10 -2 -o bperp_full.rdr
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Multilook.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to multilook')
    parser.add_argument('lks_x', type=int, help='number of multilooking in azimuth/y direction')
    parser.add_argument('lks_y', type=int, help='number of multilooking in range  /x direction')
    parser.add_argument('-o','--outfile', help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.file = ut.get_file_list(inps.file)

    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))

    # multilooking
    if len(inps.file) == 1:
        multilook_file(inps.file[0], inps.lks_y, inps.lks_x, inps.outfile)

    elif inps.parallel:
        #num_cores = min(multiprocessing.cpu_count(), len(inps.file))
        #print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(multilook_file)(file, inps.lks_y, inps.lks_x) for file in inps.file)
    else:
        for File in inps.file:
            print('-------------------------------------------')
            multilook_file(File, inps.lks_y, inps.lks_x)

    print('Done.')
    return

###################################################################################################
if __name__ == '__main__':
    main()



