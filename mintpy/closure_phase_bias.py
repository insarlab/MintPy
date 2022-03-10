#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yujie Zheng, Feb 2022                            #
############################################################
# Compute average con-nl closure phase and output mask identifying areas suseptible to closure phase errors.

import os
import sys
import argparse
import numpy as np

from mintpy.objects import ifgramStack
from mintpy.utils import readfile, writefile
from mintpy import ifgram_inversion as ifginv


################################################################################
REFERENCE = """reference:
  Zheng, Y., et al., (2022) On Closure Phase and Systematic Bias in Multilooked SAR Interferometry, IEEE TGRS, under review (minor revision)
"""
EXAMPLE = """example:
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --numsigma 2.5
"""

def create_parser():
    parser = argparse.ArgumentParser(description = 'Create an indication map for closure phase bias.')
    parser.add_argument('-i','--ifgramstack',type = str, dest = 'ifgram_stack',help = 'interferogram stack file that contains the unwrapped phases')
    parser.add_argument('--nl', dest = 'nl', type = int, default = 20, help = 'connection level that we are correcting to (or consider as no bias)')
    parser.add_argument('--numsigma',dest = 'numsigma', type = float, default = 3, help = 'Threashold for phase (number of sigmas,0-infty), default to be 3 sigma of a Gaussian distribution (assumed distribution for the cumulative closure phase) with sigma = pi/sqrt(3*num_cp)')
    parser.add_argument('--epi',dest = 'episilon', type = float, default = 0.3, help = 'Threashold for amplitude (0-1), default 0.3')
    parser.add_argument('--maxMemory', dest = 'max_memory', type = float, default = 8, help = 'max memory to use in GB')
    parser.add_argument('-o', dest = 'outdir', type = str, default = '.', help = 'output file directory')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args = iargs)
    return inps

# Obtain sum of consecutive complex sequential closure phase of connection n
def cum_seq_closurePhase(SLC_list, date12_list_all, ifgram_stack, ref_phase, n, box):
    """
    Input parameters:
        SLC_list : list of SLC dates
        date12_list_all: date12 of all the interferograms stored in the ifgramstack file
        ifgram_stack: stack file
        refphase : reference phase
        n        : connection level of the closure phase
        box      : bounding box for the patch
    """
    cp_idx = []
    NSLC = len(SLC_list)
    for i in range(NSLC-n):
        ifgram = []
        flag = True
        for j in range(n):
            ifgram.append('{}_{}'.format(SLC_list[i+j],SLC_list[i+j+1]))
        ifgram.append('{}_{}'.format(SLC_list[i],SLC_list[i+n]))
        for ifgram_name in ifgram:
            if ifgram_name not in date12_list_all:
                flag = False # if missing an interferogram, we won't make the corresponding closure phase
        if flag:
            cp_idx.append([date12_list_all.index(ifgram[j]) for j in range(n+1)])

    cp_idx = np.array(cp_idx, np.int16)
    cp_idx = np.unique(cp_idx, axis = 0)

    num_cp = len(cp_idx)
    print('Number of closure measurements expected, ', len(SLC_list)-n)
    print('Number of closure measurements found, ', num_cp)

    if num_cp <1:
        print('No closure phase measurements found, abort')
        raise Exception("No triplets found!")

    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    phase = readfile.read(ifgram_stack, box=box,print_msg=False)[0]
    cum_cp = np.zeros((box_length, box_width), np.complex64)
    for i in range(num_cp):
        cp0_w = np.zeros ((box_length, box_width), np.float32)
        for j in range(n):
                    idx = cp_idx[i,j]
                    cp0_w = cp0_w + phase[idx,:,:] - ref_phase[idx]
        idx = cp_idx[i,n]
        cp0_w = cp0_w - (phase[idx,:,:]-ref_phase[idx])
        cum_cp = cum_cp + (np.exp(1j*cp0_w))

    # cum_cp = np.angle(cum_cp)
    return cum_cp, num_cp

def main(iargs = None):
    inps = cmd_line_parse(iargs)
    stack_obj = ifgramStack(inps.ifgram_stack)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date12_list_all = stack_obj.get_date12_list(dropIfgram=False)
    print('scene length, width', length, width)
    ref_phase = stack_obj.get_reference_phase(unwDatasetName = 'unwrapPhase')
    inps.length = length
    inps.width = width
    # retrieve the list of SLC dates from ifgramStack.h5
    ifgram0 = date12_list[0]
    date1, date2 = ifgram0.split('_')
    SLC_list = [date1, date2]
    for ifgram in date12_list:
        date1, date2 = ifgram.split('_')
        if date1 not in SLC_list:
            SLC_list.append(date1)
        if date2 not in SLC_list:
            SLC_list.append(date2)
    SLC_list.sort()
    print('number of SLC found : ', len(SLC_list))
    print('first SLC: ', SLC_list[0])
    print('last  SLC: ', SLC_list[-1])


    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(inps.ifgram_stack,inps.max_memory)
    closurephase =  np.zeros([length,width],np.complex64)
    #process block-by-block
    for i, box in enumerate(box_list):
            box_width  = box[2] - box[0]
            box_length = box[3] - box[1]
            print(box)
            if num_box > 1:
                print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
                print('box width:  {}'.format(box_width))
                print('box length: {}'.format(box_length))

            closurephase[box[1]:box[3],box[0]:box[2]], numcp = cum_seq_closurePhase(SLC_list, date12_list_all, inps.ifgram_stack, ref_phase,inps.nl,box)



    # What is a good thredshold?
    # Assume that it's pure noise so that the phase is uniform distributed from -pi to pi.
    # The standard deviation of phase in each loop is pi/sqrt(3) (technically should be smaller because when forming loops there should be a reduction in phase variance)
    # The standard deviation of phase in cumulative wrapped closure phase is pi/sqrt(3)/sqrt(num_cp) -- again another simplification assuming no correlation.
    # We use 3\delta as default threshold -- 99.7% confidence

    if inps.numsigma:
        threshold_cp = np.pi/np.sqrt(3)/np.sqrt(numcp)*inps.numsigma
    else:
        threshold_cp = np.pi/np.sqrt(3)/np.sqrt(numcp)*3 # 3/sigma, 99.7% confidence

    mask = np.ones([length,width],np.float32)
    mask[np.abs(np.angle(closurephase))>threshold_cp] = 0 # this masks areas with potential bias
    mask[np.abs(np.abs(closurephase)/numcp < inps.episilon)] = 1 # this unmasks areas with low correlation (where it's hard to know wheter there is bias either)

    # save mask
    meta = dict(stack_obj.metadata)
    meta['FILE_TYPE'] = 'mask'
    ds_name_dict = {'cpmask': [np.float32, (length, width), mask],}
    writefile.layout_hdf5(os.path.join(inps.outdir,'cpmask.h5'), ds_name_dict, meta)

    # also save the average closure phase
    ds_name_dict2 = {'phase': [np.float32, (length, width), np.angle(closurephase)],
                    'amplitude':[np.float32,(length,width),np.abs(closurephase)/numcp],}
    writefile.layout_hdf5(os.path.join(inps.outdir,'avgwcp.h5'), ds_name_dict2, meta)

if __name__ == '__main__':
    main(sys.argv[1:])
