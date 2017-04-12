#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Oct 2015: Merge timeseries/velocity into one
#                   Merge all non-hdf5 into one
# Yunjun, Nov 2015: Support geomap*.trans file
# Yunjun, May 2015: add multilook() and multilook_attribute()
# Yunjun, Dec 2016: add multilook_file(), cmdLineParse() and parallel option
#                   rename multi_looking.py to multilook.py


import sys
import os
import argparse
import warnings

import h5py
import numpy as np
#from joblib import Parallel, delayed
#import multiprocessing

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


######################################## Sub Functions ############################################
def multilook_matrix(matrix,lks_y,lks_x):
    rows,cols = matrix.shape
    lks_x = int(lks_x)
    lks_y = int(lks_y)
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
  
    return matrix_mli


def multilook_attribute(atr_dict,lks_y,lks_x):
    #####
    atr = dict()
    for key, value in atr_dict.iteritems():  atr[key] = str(value)
  
    ##### calculate new data size
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    length_mli = int(np.floor(length/lks_y))
    width_mli  = int(np.floor(width/lks_x))
  
    ##### Update attributes
    atr['FILE_LENGTH'] = str(length_mli)
    atr['WIDTH'] = str(width_mli)
    atr['XMIN'] = '0'
    atr['YMIN'] = '0'
    atr['XMAX'] = str(width_mli-1)
    atr['YMAX'] = str(length_mli-1)
    
    try:
        atr['Y_STEP'] = str(lks_y*float(atr['Y_STEP']))
        atr['X_STEP'] = str(lks_x*float(atr['X_STEP']))
    except: pass
    try:
        atr['AZIMUTH_PIXEL_SIZE'] = str(lks_y*float(atr['AZIMUTH_PIXEL_SIZE']))
        atr['RANGE_PIXEL_SIZE']   = str(lks_x*float(atr['RANGE_PIXEL_SIZE']))
    except: pass
    
    try:
        atr['ref_y'] = str(int(int(atr['ref_y'])/lks_y))
        atr['ref_x'] = str(int(int(atr['ref_x'])/lks_x))
    except: pass
    try:
        atr['subset_y0'] = str(int(int(atr['subset_y0'])/lks_y))
        atr['subset_y1'] = str(int(int(atr['subset_y1'])/lks_y))
        atr['subset_x0'] = str(int(int(atr['subset_x0'])/lks_x))
        atr['subset_x1'] = str(int(int(atr['subset_x1'])/lks_x))
    except: pass
  
    return atr


def multilook_file(infile,lks_y,lks_x,outfile=None):
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    ## input file info
    atr = readfile.read_attribute(infile)
    k = atr['FILE_TYPE']
    print 'multilooking '+k+' file '+infile+' ...'

    ## output file name
    if not outfile:
        ext = os.path.splitext(infile)[1]
        outfile = os.path.splitext(infile)[0]+'_'+str(lks_y)+'alks_'+str(lks_x)+'rlks'+ext
    print 'writing >>> '+outfile

    ###############################################################################
    ## Read/Write multi-dataset files
    if k in ['interferograms','coherence','wrapped','timeseries']:
        h5 = h5py.File(infile,'r')
        epochList = sorted(h5[k].keys())
        print 'number of epochs: '+str(len(epochList))

        h5out = h5py.File(outfile,'w')
        group = h5out.create_group(k)

        if k in ['interferograms','coherence','wrapped']:
            for epoch in epochList:
                print epoch
                data = h5[k][epoch].get(epoch)[:]
                atr = h5[k][epoch].attrs

                data_mli = multilook_matrix(data,lks_y,lks_x)
                atr_mli = multilook_attribute(atr,lks_y,lks_x)

                gg = group.create_group(epoch)
                dset = gg.create_dataset(epoch, data=data_mli, compression='gzip')
                for key, value in atr_mli.iteritems():
                    gg.attrs[key] = value

        elif k == 'timeseries':
            for epoch in epochList:
                print epoch
                data = h5[k].get(epoch)[:]

                data_mli = multilook_matrix(data,lks_y,lks_x)
                
                dset = group.create_dataset(epoch, data=data_mli, compression='gzip')
            atr = h5[k].attrs
            atr_mli = multilook_attribute(atr,lks_y,lks_x)
            for key, value in atr_mli.iteritems():
                group.attrs[key] = value

        h5.close()
        h5out.close()

    ## Read/Write single-dataset files
    elif k == '.trans':        
        rg,az,atr = readfile.read(infile)
        rgmli = multilook_matrix(rg,lks_y,lks_x);
        azmli = multilook_matrix(az,lks_y,lks_x);
        atr = multilook_attribute(atr,lks_y,lks_x)
        writefile.write(rgmli,azmli,atr,outfile)
    else:
        data,atr = readfile.read(infile)
        data_mli = multilook_matrix(data,lks_y,lks_x)
        atr = multilook_attribute(atr,lks_y,lks_x)
        writefile.write(data_mli,atr,outfile)

    return outfile


##################################################################################################
EXAMPLE='''example:
  multilook.py  velocity.h5  15 15
  multilook.py  srtm30m.dem  10 10  -o srtm30m_300m.dem
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Multilook.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to multilook')
    parser.add_argument('lks_x', type=int, help='number of multilooking in azimuth/y direction')
    parser.add_argument('lks_y', type=int, help='number of multilooking in range  /x direction')
    parser.add_argument('-o','--outfile', help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')

    inps = parser.parse_args()
    return inps


##################################################################################################
def main(argv):

    inps = cmdLineParse()
    #print '\n**************** Multilook *********************'
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
            print '-------------------------------------------'
            multilook_file(File, inps.lks_y, inps.lks_x)

    print 'Done.'
    return

###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



