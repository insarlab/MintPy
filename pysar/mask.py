#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Oct 2015: add support for ROI_PAC product
# Yunjun, Jul 2016: add mask_matrix(), mask_file()
#                   add parallel processing using joblib


import os
import sys
import getopt
import glob

import h5py
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar._writefile as writefile


############################################################
def mask_matrix(data_mat,mask_mat):
    ## mask a 2D matrxi data with mask
    try:
        xsub
        ysub
        mask_mat[ysub[0]:ysub[1],xsub[0]:xsub[1]]=0
    except:   pass

    ## Masked Value
    if data_mat.dtype == np.dtype('int16'):
        mask_value = np.ma.masked
    else:
        mask_value = np.nan

    try:     data_mat[mask_mat<thr] = mask_value
    except:  data_mat[mask_mat==0]  = mask_value

    return data_mat


############################################################
def mask_with_multi_masks(in_file,mask_file,out_file=''):

    h5file=h5py.File(in_file,'r')
    h5mask=h5py.File(mask_file,'r')
    kf=h5file.keys()

    if out_file == '':
        ext = os.path.splitext(in_file)[1]
        out_file = in_file.split('.')[0]+'_masked'+ext
    h5out = h5py.File(out_file,'w')

    if kf[0] in ('interferograms','wrapped','coherence') and 'coherence' in h5mask.keys():
        print 'file type: '+kf[0]
        print 'mask each '+kf[0]+' using its coherence file'
        igramList = sorted(h5file[kf[0]].keys())
        cohList   = sorted(h5mask['coherence'].keys())
        gg = h5out.create_group(kf[0])
        for igram in igramList:
            print igram
            date12 = h5file[kf[0]][igram].attrs['DATE12']
            for cohFile in cohList:
                if h5mask['coherence'][cohFile].attrs['DATE12']==date12:
                    igramCoh = cohFile
            print igramCoh

            unwset = h5file[kf[0]][igram].get(igram)
            unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

            cohset=h5mask['coherence'][igramCoh].get(igramCoh)
            coh=cohset[0:cohset.shape[0],0:cohset.shape[1]]

            unw = mask_matrix(unw,coh)

            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=unw, compression='gzip')
            for key, value in h5file[kf[0]][igram].attrs.iteritems():   group.attrs[key] = value
            del unw, coh
        try:
            mask = h5file['mask'].get('mask')
            gm = h5out.create_group('mask')
            dset = gm.create_dataset('mask', data=mask, compression='gzip')
        except: print 'no mask group found.'

    else: print 'Only support multiple dataset file with multiple masks.'; sys.exit(1)


############################################################
def mask_file(in_file,M,out_file=''):
    ## Mask input file with mask matrix M

    atr = readfile.read_attribute(in_file)
    k = atr['FILE_TYPE']
    print 'file type: '+k

    if out_file == '':
        ext      = os.path.splitext(in_file)[1]
        out_file = os.path.basename(in_file).split('.')[0]+'_masked'+ext

    if k in ['timeseries','interferograms','wrapped','coherence']:
        h5file = h5py.File(in_file,'r')
        epochList = sorted(h5file[k].keys())
        print 'number of epochs: '+str(len(epochList))

        h5out = h5py.File(out_file,'w')
        print 'writing >>> '+out_file

    ##### Multiple Dataset File
    if k == 'timeseries':
        group = h5out.create_group(k)
        for d in epochList:
            print d
            unwset = h5file[k].get(d)
            unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

            unw = mask_matrix(unw,M)

            dset = group.create_dataset(d, data=unw, compression='gzip')
        for key,value in atr.iteritems():   group.attrs[key] = value

    elif k in ['interferograms','wrapped','coherence']:
        gg = h5out.create_group(k)
        for igram in epochList:
            print igram
            unwset = h5file[kf[0]][igram].get(igram)
            unw=unwset[0:unwset.shape[0],0:unwset.shape[1]]

            unw = mask_matrix(unw,M)

            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=unw, compression='gzip')
            for key, value in h5file[k][igram].attrs.iteritems():
                group.attrs[key] = value
        try:
            mask = h5file['mask'].get('mask')
            gm = h5out.create_group('mask')
            dset = gm.create_dataset('mask', data=mask, compression='gzip')
        except: print 'no mask group found.'

    ##### Single Dataset File
    else:
        unw,atr = readfile.read(in_file)
        unw     = mask_matrix(unw,M)
        writefile.write(unw,atr,out_file)

    try:
        h5file.close()
        h5out.close()
    except: pass


############################################################
def usage():
    print '''
**************************************************************************
  Mask File(s)

  Usage:
      mask.py file MaskFile
      mask.py -f file -m MaskFile [ -t threshold -o output_name -x x_subset -y y_subset]

      -f : file (list) that needed to be masked
      -m : mask file 
      -o : output file name [works only for one input file]

      -t : threshold value used for masking. if not 
           specified then only pixels with mask value 
           equal to zero is masked out.
      -x : mask subset in range   / column / longtitude direction
      -y : mask subset in azimuth / row    / latitude   direction

      --no-parallel : disable parallel computing [enabled by default for multiple input files]

  Example:
      mask One File
          mask.py velocity Mask.h5
          mask.py -f geo_100102_101120.unw -m Mask.h5
          mask.py -f timeseries.h5         -m temporal_coherence.h5 -t 0.7
          mask.py -f LoadedData.h5         -m 100102_101120.cor     -t 0.9 -y '200:300' -x '300:400'

      mask Multiple Files
          mask.py -f 'timeseries*.h5'              -m Mask.h5
          mask.py -f 'timeseries*.h5'              -m Mask.h5 --no-parallel
          mask.py -f 'timeseries*.h5,*velocity.h5' -m Mask.h5

**************************************************************************
    '''

############################################################
def main(argv):

    global xsub, ysub, thr
    parallel = 'yes'     ## Use parallel by default for multiple input files

    ######################################
    try:    opts, args = getopt.getopt(argv,"h:f:m:t:x:y:o:",['no-parallel'])
    except getopt.GetoptError:    usage() ; sys.exit(1)

    if len(sys.argv) > 3:
        for opt,arg in opts:
            if opt in ("-h","--help"):     usage();  sys.exit()
            elif opt == '-f':        File     = arg.split(',')
            elif opt == '-m':        maskFile = arg
            elif opt == '-t':        thr  = float(arg)
            elif opt == '-y':        ysub = sorted([int(i) for i in arg.split(':')])
            elif opt == '-x':        xsub = sorted([int(i) for i in arg.split(':')])
            elif opt == '-o':        outFile = arg
            elif opt == '--no-parallel':   parallel = 'no'

    elif len(sys.argv)==3:
        if os.path.isfile(argv[0]) and os.path.isfile(argv[1]):
            File     = argv[0].split(',')
            maskFile = argv[1]
        else:  print 'Input file does not existed: '+argv[0]+' / '+argv[1];  sys.exit(1)
    else:   usage();  sys.exit(1)

    try:
        File
        maskFile
    except:    usage() ; sys.exit(1)

    ##### Check Input File List
    print '\n****************** mask *********************'
    fileList = ut.get_file_list(File)
    print 'number of file to mask: '+str(len(fileList))
    print fileList

    if len(fileList) == 1:
        parallel = 'no'
        try: outFile          ## Customized output file name for one input file only
        except:
            ext     = os.path.splitext(fileList[0])[1]
            outFile = os.path.basename(fileList[0]).split('.')[0]+'_masked'+ext
    elif len(fileList) > 1:
        try:
            del outFile
            print 'Disabled customized output name for multiple input files, continue with automatic naming insread.'
        except: pass
    else: print 'ERROR: No input file!'; sys.exit(1)

    ###### Read Mask File
    atr_mask = readfile.read_attribute(maskFile)
    k_mask = atr_mask['FILE_TYPE']
    if not k_mask == 'coherence':    ## Read mask file once 
        M,Matr = readfile.read(maskFile)
        print 'mask file: '+maskFile

    ##### mask - file by file
    if parallel == 'no':
        ##### Single Mask
        if not k_mask == 'coherence':
            for in_file in fileList:
                print '-------------------------------------------'
                print 'masking  : '+in_file
                try:    mask_file(in_file,M,outFile)
                except: mask_file(in_file,M)
        ##### Multiple Mask
        else:
            try:    mask_with_multi_masks(fileList[0],maskFile,outFile)
            except: mask_with_multi_masks(fileList[0],maskFile)

    ##### mask - parallel
    else:
        print '-----------------------'
        print 'parallel masking ...'
        print '-----------------------'
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(mask_file)(in_file,M) for in_file in fileList)


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])


