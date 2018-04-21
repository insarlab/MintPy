#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
from pysar.utils import readfile, writefile, utils as ut
from pysar.utils.readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


############################################################
def mask_matrix(data_mat,mask_mat, fill_value=None):
    '''mask a 2D matrxi data with mask'''
    ## Masked Value
    if fill_value is None:
        if data_mat.dtype == np.dtype('int16'):
            fill_value = np.ma.masked
        else:
            fill_value = np.nan
    #data_mat = data_mat.astype(np.float32)
    #mask_value = np.nan

    data_mat[mask_mat==0]  = fill_value

    return data_mat


############################################################
def update_mask(mask, inps_dict, print_msg=True):
    '''Update mask matrix from input options: subset_x/y and threshold'''
    if inps_dict['subset_x']:
        mask[:,0:inps_dict['subset_x'][0]] = 0
        mask[:,inps_dict['subset_x'][1]:] = 0
        if print_msg:
            print('mask out area not in x: '+str(inps_dict['subset_x']))
    if inps_dict['subset_y']:
        mask[0:inps_dict['subset_y'][0],:] = 0
        mask[inps_dict['subset_y'][1]:,:] = 0
        if print_msg:
            print('mask out area not in y: '+str(inps_dict['subset_y']))
    if inps_dict['thr']:
        mask[mask<inps_dict['thr']] = 0
        if print_msg:
            print('mask out pixels < '+str(inps_dict['thr'])+' in mask file')
    return mask


############################################################
def mask_file(File, maskFile, outFile=None, inps_dict=None):
    ''' Mask input File with maskFile
    Inputs:
        File/maskFile - string, 
        inps_dict - dictionary including the following options:
                    subset_x/y - list of 2 ints, subset in x/y direction
                    thr - float, threshold/minValue to generate mask
    Output:
        outFile - string
    '''
    
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    print('masking '+k+' file: '+File+' ...')

    # Read maskFile
    atrm = readfile.read_attribute(maskFile)
    km = atrm['FILE_TYPE']
    if km not in multi_group_hdf5_file+multi_dataset_hdf5_file:
        print('reading mask file: '+maskFile)
        mask = readfile.read(maskFile, datasetName='mask')[0]
        if inps_dict:
            mask = update_mask(mask, inps_dict)
    
    if not outFile:
        outFile = os.path.splitext(File)[0]+'_masked'+os.path.splitext(File)[1]

    if k in ['timeseries','interferograms','wrapped','coherence']:
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())

        h5out = h5py.File(outFile,'w')
        print('writing >>> '+outFile)

    ##### Multiple Dataset File
    if k == 'timeseries':
        print('number of acquisitions: '+str(len(epochList)))
        group = h5out.create_group(k)
        for d in epochList:
            print(d)
            unw = h5file[k].get(d)[:]

            unw = mask_matrix(unw, mask, inps_dict['fill_value'])

            dset = group.create_dataset(d, data=unw)
        for key,value in iter(atr.items()):
            group.attrs[key] = value

    elif k in ['interferograms','wrapped','coherence']:
        print('number of interferograms: '+str(len(epochList)))
        gg = h5out.create_group(k)
        
        # Mask multi group file with multi group coherence file
        if km == 'coherence':
            h5mask = h5py.File(maskFile, 'r')
            cohList = sorted(h5mask[km].keys())
            if len(cohList) != len(epochList):
                sys.exit('ERROR: cohERROR: erence mask file has different\
                number of interferograms than input file!')

        for i in range(len(epochList)):
            igram = epochList[i]
            if km == 'coherence':
                coh = cohList[i]
                sys.stdout.write('\r%s %s %s/%s ...' % (igram, coh, i+1, len(epochList)))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r%s %s/%s ...' % (igram, i+1, len(epochList)))
                sys.stdout.flush()

            unw = h5file[k][igram].get(igram)[:]

            if km == 'coherence':
                mask = h5mask[km][coh].get(coh)[:]
                if inps_dict:
                    mask = update_mask(mask, inps_dict, print_msg=False)

            unw = mask_matrix(unw, mask, inps_dict['fill_value'])

            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=unw)
            for key, value in h5file[k][igram].attrs.items():
                group.attrs[key] = value

    ##### Single Dataset File
    elif k in ['.trans','.utm_to_rdc','.UTM_TO_RDC']:
        rg, az, atr = readfile.read(File)
        rg = mask_matrix(rg, mask, inps_dict['fill_value'])
        az = mask_matrix(az, mask, inps_dict['fill_value'])
        print('writing >>> '+outFile)
        writefile.write(rg, az, atr, outFile)

    else:
        unw,atr = readfile.read(File)
        unw     = mask_matrix(unw, mask, inps_dict['fill_value'])
        print('writing >>> '+outFile)
        writefile.write(unw,atr,outFile)

    try: h5file.close()
    except: pass
    try: h5out.close()
    except: pass
    try: h5mask.close()
    except: pass
    return outFile
    

############################################################
EXAMPLE='''example:
  mask.py  velocity.h5     -m Mask.h5
  mask.py  timeseries.h5   -m temporal_coherence.h5  -t 0.7
  mask.py  unwrapIfgram.h5 -m 100102_101120.cor      -t 0.9  -y  200 300  -x 300 400
  mask.py  unwrapIfgram.h5 -m coherence.h5           -t 0.1  --fill 0
  mask.py  timeseries*.h5 velocity*.h5  -m temporal_coherence.h5  -t 0.7
'''

def create_parser():
    parser = argparse.ArgumentParser(description='Mask File(s)',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    
    parser.add_argument('file', nargs='+', help='File(s) for ramp removal')
    parser.add_argument('-m','--mask', dest='mask_file', help='mask for pixels used in ramp estimation')
    parser.add_argument('-t', dest='thr', type=float,\
                        help='threshold value used for masking.\n'+\
                        'if not specified, only pixels with mask value equal to zero is masked out.')
    parser.add_argument('--fill', dest='fill_value', type=float,\
                        help="fill masked out area with input value. i.e. \n"
                             "np.nan, 0, 1000, ... \n"
                             "By default, it's np.ma.masked for int16 type and np.nan for all the others.")
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2, help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2, help='subset range in y/along-track/row direction')
    parser.add_argument('-o','--outfile', help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false', default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None): 
    inps = cmd_line_parse(iargs)
    #print '\n****************** mask *********************'
    inps.file = ut.get_file_list(inps.file)
    print('number of file to mask: '+str(len(inps.file)))
    print(inps.file)

    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))

    # masking
    if len(inps.file) == 1:
        mask_file(inps.file[0], inps.mask_file, inps.outfile, vars(inps))
    
    elif inps.parallel:
        #num_cores = min(multiprocessing.cpu_count(), len(inps.file))
        #print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(mask_file)(File, inps.mask_file, inps_dict=vars(inps)) for File in inps.file)
    else:
        for File in inps.file:
            print('-------------------------------------------')
            mask_file(File, inps.mask_file, inps_dict=vars(inps))

    print('Done.')
    return


############################################################
if __name__ == '__main__':
    main()

