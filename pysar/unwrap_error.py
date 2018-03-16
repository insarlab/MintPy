#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse

import h5py
import numpy as np
from scipy.linalg import pinv

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
import pysar._remove_surface as rm


##########################################################################################
def bridging_data(data,mask,x,y):
    '''Phase Jump Correction, using phase continuity on bridge/bonding points in each pair of patches.
    Inputs:
        data : 2D np.array, phase matrix need to be corrected
        mask : mask file marks different patches with different positive integers
        x/y  : list of int, array of bridge points, lied as: x_ref, x, x_ref, x
    Output:
        data : 2D np.array, phase corrected matrix
    '''

    ## loop based on number of bridges
    n_bridge = len(x)/2
    for i in range(1,n_bridge+1):
        p_ref = data[y[2*i-2],x[2*i-2]]
        p     = data[y[2*i-1],x[2*i-1]]
        n_jump = (abs(p-p_ref)+np.pi)//(2*np.pi)
        if not n_jump == 0:
            if p-p_ref >=0:  n_jump *= -1
            id = np.where(mask == mask[y[2*i-1],x[2*i-1]])
            data[id] = data[id] + n_jump*2*np.pi;
  
    return data


def unwrap_error_correction_phase_closure(ifgram_file, mask_file, ifgram_cor_file=None):
    '''Correct unwrapping errors in network of interferograms using phase closure.
    Inputs:
        ifgram_file     - string, name/path of interferograms file
        mask_file       - string, name/path of mask file to mask the pixels to be corrected
        ifgram_cor_file - string, optional, name/path of corrected interferograms file
    Output:
        ifgram_cor_file
    Example:
        'unwrapIfgram_unwCor.h5' = unwrap_error_correction_phase_closure('Seeded_unwrapIfgram.h5','mask.h5')
    '''
    print('read mask from file: '+mask_file)
    mask = readfile.read(mask_file, epoch='mask')[0].flatten(1)

    atr = readfile.read_attribute(ifgram_file)
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])
    k = atr['FILE_TYPE']
    pixel_num = length*width

    # Check reference pixel
    try:
        ref_y = int(atr['ref_y'])
        ref_x = int(atr['ref_x'])
        print('reference pixel in y/x: %d/%d' % (ref_y, ref_x))
    except:
        sys.exit('ERROR: Can not find ref_y/x value, input file is not referenced in space!')

    h5 = h5py.File(ifgram_file,'r')
    ifgram_list = sorted(h5[k].keys())
    ifgram_num = len(ifgram_list)

    ##### Prepare curls
    curls, Triangles, C = ut.get_triangles(h5)
    curl_num = np.shape(curls)[0]
    print('Number of      triangles: '+  str(curl_num))

    curl_file='curls.h5'
    if not os.path.isfile(curl_file):
        print('writing >>> '+curl_file)
        ut.generate_curls(curl_file, h5, Triangles, curls)

    thr=0.50
    curls = np.array(curls);   n1=curls[:,0];   n2=curls[:,1];   n3=curls[:,2]

    print('reading interferograms...')
    print('Number of interferograms: '+ str(ifgram_num))
    data = np.zeros((ifgram_num,pixel_num),np.float32)
    prog_bar = ptime.progress_bar(maxValue=ifgram_num)
    for ni in range(ifgram_num):
        ifgram = ifgram_list[ni]
        d = h5[k][ifgram].get(ifgram)[:].flatten(1)
        data[ni,:] = d
        prog_bar.update(ni+1)
    prog_bar.close()

    print('reading curls ...') 
    print('number of culrs: '+str(curl_num))
    h5curl = h5py.File(curl_file,'r')
    curl_list = sorted(h5curl[k].keys())
    curl_data = np.zeros((curl_num, pixel_num),np.float32)
    prog_bar = ptime.progress_bar(maxValue=curl_num)
    for ni in range(curl_num):
        d = h5curl[k][curl_list[ni]].get(curl_list[ni])[:].flatten(1)
        curl_data[ni,:] = d.flatten(1)
        prog_bar.update(ni+1)
    prog_bar.close()
    h5curl.close() 

    print('estimating unwrapping error pixel by pixel ...')
    EstUnwrap = np.zeros((ifgram_num,pixel_num),np.float32)
    prog_bar = ptime.progress_bar(maxValue=pixel_num)
    for ni in range(pixel_num):
        if mask[ni]==1:
            dU = data[:,ni]
            unwCurl = np.array(curl_data[:,ni])

            ind  = np.abs(unwCurl)>=thr;      N1 =n1[ind];      N2 =n2[ind];      N3 =n3[ind]
            indC = np.abs(unwCurl)< thr;      Nc1=n1[indC];     Nc2=n2[indC];     Nc3=n3[indC]
  
            N =np.hstack([N1, N2, N3]);       UniN =np.unique(N)
            Nc=np.hstack([Nc1,Nc2,Nc3]);      UniNc=np.unique(Nc)

            inter = list(set(UniNc) & set(UniN)) # intersetion
            UniNc = list(UniNc)
            for x in inter:
                UniNc.remove(x)

            D = np.zeros([len(UniNc),ifgram_num])
            for i in range(len(UniNc)):
                D[i,UniNc[i]]=1

            AAA  = np.vstack([-2*np.pi*C,D])
            AAAA = np.vstack([AAA,0.25*np.eye(ifgram_num)])

            ##########
            # with Tikhonov regularization:
            LLL = list(np.dot(C,dU)) + list(np.zeros(np.shape(UniNc)[0])) + list(np.zeros(ifgram_num))
            ind = np.isnan(AAAA)
            M1 = pinv(AAAA)
            M = np.dot(M1,LLL)
            EstUnwrap[:,ni] = np.round(M[0:ifgram_num])*2.0*np.pi
        prog_bar.update(ni+1, suffix='%s/%d' % (ni,pixel_num))
    prog_bar.close()

    dataCor = data + EstUnwrap

    ##### Output
    if not ifgram_cor_file:
        ifgram_cor_file = os.path.splitext(ifgram_file)[0]+'_unwCor.h5'
    print('writing >>> '+ifgram_cor_file)
    h5unwCor = h5py.File(ifgram_cor_file,'w') 
    gg = h5unwCor.create_group(k) 

    prog_bar = ptime.progress_bar(maxValue=ifgram_num)
    for i in range(ifgram_num):
        ifgram = ifgram_list[i]
        group = gg.create_group(ifgram)
        dset = group.create_dataset(ifgram, data=np.reshape(dataCor[i,:],[width,length]).T, compression='gzip')
        for key, value in h5[k][ifgram].attrs.items():
            group.attrs[key] = value
        prog_bar.update(i+1)
    prog_bar.close()
    h5unwCor.close()
    h5.close()
    return ifgram_cor_file


def unwrap_error_correction_bridging(ifgram_file, mask_file, y_list, x_list, ramp_type='plane',\
                                     ifgram_cor_file=None, save_cor_deramp_file=False):
    '''Unwrapping error correction with bridging.
    Inputs:
        ifgram_file : string, name/path of interferogram(s) to be corrected
        mask_file   : string, name/path of mask file to mark different patches 
        y/x_list    : list of int, bonding points in y/x 
        ifgram_cor_file : string, optional, output file name
        save_cor_deramp_file : bool, optional
    Output:
        ifgram_cor_file
    Example:
        y_list = [235, 270, 350, 390]
        x_list = [880, 890, 1200, 1270]
        unwrap_error_correction_bridging('unwrapIfgram.h5', 'mask_all.h5', y_list, x_list, 'quadratic')
    '''
    ##### Mask and Ramp
    mask = readfile.read(mask_file, epoch='mask')[0]
    ramp_mask = mask == 1
    print('estimate phase ramp during the correction')
    print('ramp type: '+ramp_type)

    ##### Bridge Info
    # Check
    for i in range(len(x_list)):
        if mask[y_list[i],x_list[i]] == 0:
            print('\nERROR: Connecting point (%d,%d) is out of masked area! Select them again!\n' % (y_list[i],x_list[i]))
            sys.exit(1)
    print('Number of bridges: '+str(len(x_list)/2))
    print('Bonding points coordinates:\nx: '+str(x_list)+'\ny: '+str(y_list))

    # Plot Connecting Pair of Points
    plot_bonding_points = False
    if plot_bonding_points:
        point_yx = ''
        line_yx  = ''
        n_bridge = len(x)/2
        for i in range(n_bridge):
            pair_yx = str(y[2*i])+','+str(x[2*i])+','+str(y[2*i+1])+','+str(x[2*i+1])
            if not i == n_bridge-1:
                point_yx += pair_yx+','
                line_yx  += pair_yx+';'
            else:
                point_yx += pair_yx
                line_yx  += pair_yx

        try:
            plot_cmd = 'view.py --point-yx="'+point_yx+'" --line-yx="'+line_yx+\
                       '" --nodisplay -o bonding_points.png -f '+maskFile
            print(plot_cmd)
            os.system(plot_cmd)
        except: pass

    # Basic info
    ext = os.path.splitext(ifgram_file)[1]
    atr = readfile.read_attribute(ifgram_file)
    k = atr['FILE_TYPE']

    try:
        ref_y = int(atr['ref_y'])
        ref_x = int(atr['ref_x'])
        print('reference pixel in y/x: %d/%d' % (ref_y, ref_x))
    except:
        sys.exit('ERROR: Can not find ref_y/x value, input file is not referenced in space!')

    # output file name
    if not ifgram_cor_file:
        ifgram_cor_file = os.path.splitext(ifgram_file)[0]+'_unwCor'+ext
    ifgram_cor_deramp_file = os.path.splitext(ifgram_cor_file)[0]+'_'+ramp_type+ext

    ##### HDF5 file
    if ext == '.h5':
        ##### Read
        h5 = h5py.File(ifgram_file,'r')
        ifgram_list = sorted(h5[k].keys())
        ifgram_num = len(ifgram_list)

        h5out = h5py.File(ifgram_cor_file,'w')
        group = h5out.create_group(k)
        print('writing >>> '+ifgram_cor_file)

        if save_cor_deramp_file:
            h5out_deramp = h5py.File(ifgram_cor_deramp_file,'w')
            group_deramp = h5out_deramp.create_group(k)
            print('writing >>> '+ifgram_cor_deramp_file)

        ##### Loop
        print('Number of interferograms: '+str(ifgram_num))
        prog_bar = ptime.progress_bar(maxValue=ifgram_num)
        date12_list = ptime.list_ifgram2date12(ifgram_list)
        for i in range(ifgram_num):
            ifgram = ifgram_list[i]
            data = h5[k][ifgram].get(ifgram)[:]
            data -= data[ref_y, ref_x]

            data_deramp, ramp = rm.remove_data_surface(data, ramp_mask, ramp_type)
            data_derampCor = bridging_data(data_deramp, mask, x_list, y_list)

            ramp[data == 0.] = 0.
            gg = group.create_group(ifgram)
            dset = gg.create_dataset(ifgram, data=data_derampCor+ramp, compression='gzip')
            for key, value in h5[k][ifgram].attrs.items():
                gg.attrs[key]=value

            if save_cor_deramp_file:
                gg_deramp = group_deramp.create_group(ifgram)
                dset = gg_deramp.create_dataset(ifgram, data=data_derampCor, compression='gzip')
                for key, value in h5[k][ifgram].attrs.items():
                    gg_deramp.attrs[key]=value
            prog_bar.update(i+1, suffix=date12_list[i])

        prog_bar.close()
        h5.close()
        h5out.close()
        try: h5out_deramp.close()
        except: pass

    #### .unw file
    elif ext == '.unw':
        print('read '+ifgram_file)
        data = readfile.read(ifgram_file)[0]
        data -= data[ref_y, ref_x]

        data_deramp,ramp = rm.remove_data_surface(data,ramp_mask,ramp_type)
        data_derampCor = bridging_data(data_deramp,mask,x_list,y_list)

        print('writing >>> '+ifgram_cor_file)
        ramp[data == 0.] = 0.
        ifgram_cor_file = writefile.write(data_derampCor+ramp, atr, ifgram_cor_file)
        if save_cor_deramp_file:
            print('writing >>> '+ifgram_cor_deramp_file)
            ifgram_cor_deramp_file = writefile.write(data_derampCor, atr, ifgram_cor_deramp_file)

    else:
        sys.exit('Un-supported file type: '+ext)

    return ifgram_cor_file, ifgram_cor_deramp_file


def read_template2inps(template_file, inps=None):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()

    print('read options from tempalte file: '+os.path.basename(inps.template_file))
    template = readfile.read_template(inps.template_file)
    key_list = list(template.keys())

    # Coherence-based network modification
    prefix = 'pysar.unwrapError.'

    key = prefix+'method'
    if key in key_list:
        value = template[key]
        if value in ['bridging','phase_closure']:
            inps.method = value
        elif value not in ['auto','no']:
            inps.method = None
        else:
            print('Unrecognized input for %s: %s' % (key, value))

    key = prefix+'maskFile'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.mask_file = value

    key = prefix+'yx'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            yx = value.replace(';',' ').replace(',',' ').split()
            yx = [int(i) for i in yx]
            inps.y = yx[0::2]
            inps.x = yx[1::2]

    key = prefix+'ramp'
    if key in key_list:
        value = template[key]
        if value in ['auto']:
            inps.ramp_type = 'plane'
        elif value in ['plane','quadratic']:
            inps.ramp_type = value
        else:
            print('Unrecognized input for %s: %s' % (key, value))

    return inps


####################################################################################################
EXAMPLE='''example:
Phase Closure:
  unwrap_error.py  Seeded_unwrapIfgram.h5  --mask mask.h5
Bridging:
  unwrap_error.py  unwrapIfgram.h5    -t ShikokuT417F650_690AlosA.template
  unwrap_error.py  unwrapIfgram.h5    --mask mask.h5     -x 283 305 -y 1177 1247
  unwrap_error.py  081018_090118.unw  --mask mask_all.h5 -x 283 305 -y 1177 1247 --ramp quadratic
'''

TEMPLATE='''
## 4. Unwrapping Error Correction
## unwrapping error correction based on the following two methods:
## a. phase closure (Fattahi, 2015, PhD Thesis)
## b. connecting bridge
pysar.unwrapError.method   = auto   #[bridging / phase_closure / no], auto for no
pysar.unwrapError.maskFile = auto   #[file name / no], auto for no
pysar.unwrapError.ramp     = auto   #[plane / quadratic], auto for plane
pysar.unwrapError.yx       = auto   #[y1_start,x1_start,y1_end,x1_end;y2_start,...], auto for none
'''

REFERENCE='''reference:
  Fattahi, H. (2015), Geodetic Imaging of Tectonic Deformation with InSAR, 190 pp, University of Miami, Miami, FL.
'''

DESCRIPTION='''
  Two methods: 1) Phase closure, 2) Bridging
  -------------------------------------------------------------------
  1. Phase closure: correct unwrapping errors based on triangular consistency
      Based on phase closure of pairs circle (ab + bc + ca == 0), this method assumes
      a. abundance of network: for interferogram with unwrapping error, there is
         at least of one triangular connection to form a closed circle; with more
         closed circles comes better constrain.
      b. majority rightness: most of interferograms have to be right (no unwrapping
         error) to correct the wrong minority. And if most of interferograms have 
         unwrapping errors, then the minor right interferograms will turn into wrong.

  -------------------------------------------------------------------
  2. Bridging: correct unwrapping errors based on close bonding points
      This method assumes:
      a. no phase unwrapping error within each patch marked by mask file.
      b. the absolute phase difference of bonding points (usually close in space) is 
         smaller than one pi. Considering prevalent ramps in InSAR data might break
         this assumptioin for bonding points that are not very close, across a bay
         for example, we first estimate and remove a linear phase ramp, then applied
         phase continuity constrain, and add the removed ramp back at the end.
 
      Phase unwrapping error is corrected epoch by epoch, following the steps below:
      a. estimate and remove a linear phase ramp from unwrapped phase;
      b. following the pair order of bonding points, correct patch by patch marked
         by point's coordinate and mask file:
         1) use 1st point as reference, calculate integer N, add N*2pi to 2nd point's
            phase to make sure their absolute phase difference is smaller than pi.
         2) add N*2pi to all pixels in 2nd point's patch.
      c. add linear phase ramp estimated in step a back to the corrected phase in step b.
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Unwrapping Error Correction.'+DESCRIPTION,\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be corrected')
    parser.add_argument('--mask', dest='mask_file',\
                        help='mask file used for correction.\n'+\
                             'For phase closure method, to specify those pixels to be corrected for unwrapping errors\n'+\
                             'For bridging method, to mark different patches that want to be corrected.\n'+\
                             '    Masked out area is marked with 0, patches/area needed to be corrected marked with\n'+\
                             '    positive integers, i.e. 1, 2, 3, ...')
    parser.add_argument('--method', dest='method', choices=['bridging','phase_closure'],\
                        help='method used for error correction.')
    parser.add_argument('-o','--outfile', help="output file name. Default is to add suffix '_unwCor.h5'")

    bridging = parser.add_argument_group('Bridging')
    bridging.add_argument('-y', type=int, nargs='*',\
                          help='Y coordinates of bridge bonding points from reference patch to to-be-corrected patch.\n'+\
                               'e.g. 283 305 350 390')
    bridging.add_argument('-x', type=int, nargs='*',\
                          help='X coordinates of bridge bonding points from reference patch to to-be-corrected patch.\n'+\
                               'e.g. 1177 1247 2100 2200\n'+\
                               'Note: choose x/y_ref point in the patch that also have seed point,'+\
                               ' for consistency in multiple images.')
    bridging.add_argument('-t','--template', dest='template_file',\
                          help='template file with bonding point info, e.g.\n'+\
                               'pysar.unwrapError.yx = 283,1177,305,1247;350,2100,390,2200')
    bridging.add_argument('--ramp', dest='ramp_type', choices=['plane','quadratic'], default='plane',\
                          help='type of phase ramp to be removed before correction.')

    inps = parser.parse_args()
    if inps.y and np.mod(len(inps.y),2) != 0:
        raise argparse.ArgumentTypeError('Number of Y coordinates is not even')
    if inps.x and np.mod(len(inps.x),2) != 0:
        raise argparse.ArgumentTypeError('Number of X coordinates is not even')
    return inps


####################################################################################################
def main(argv):
    inps = cmdLineParse()
    # output filename
    ext = os.path.splitext(inps.ifgram_file)[1]
    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.ifgram_file)[0]+'_unwCor'+ext

    # read template file
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # Memthod
    if not inps.method:
        if inps.y and inps.x:
            inps.method = 'bridging'
        else:
            inps.method = 'phase_closure'
    print('unwrapping error correction using method: '+inps.method)

    #####
    if inps.method == 'phase_closure':
        inps.outfile = unwrap_error_correction_phase_closure(inps.ifgram_file, inps.mask_file, inps.outfile)

    elif inps.method == 'bridging':
        inps.outfile = unwrap_error_correction_bridging(inps.ifgram_file, inps.mask_file, inps.y, inps.x,\
                                                        inps.ramp_type, inps.outfile)[0]

    print('Done.')
    return inps.outfile


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

