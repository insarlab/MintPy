#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yujie Zheng, Feb 2022                            #
############################################################
# Estimate fading signal in level-n time-series


import os
import sys
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from mintpy.objects import ifgramStack
from mintpy.utils import readfile, writefile, ptime, utils as ut
from mintpy import ifgram_inversion as ifginv



################################################################################
REFERENCE = """reference:
  Zheng, Y., et al., (2022) On the Closure Phase ...
"""

EXAMPLE = """example:
  closure_phase_bias.py -i inputs/ifgramStack.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Create an indication map for closure phase bias',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)
    parser.add_argument('-i','--ifgram-stack', type=str, dest='ifgram_stack_file',
                        help='interferogram stack file that contains the unwrapped phases')
    parser.add_argument('--nl', dest='nl', type=int, default=20,
                        help='connection level that we are correcting to (or consider as no bias)')
    parser.add_argument('--maxMemory', dest='max_memory', type=float, default=8,
                        help='max memory to use in GB')
    parser.add_argument('-o', dest='outdir', type=str, default='.',
                        help='output file directory')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args = iargs)
    return inps


################################################################################
# Obtain sum of consecutive sequential closure phase of connection n
def cum_seq_closure_phase(inps,n,box): # we call triplets 2-level closure phase.
    cp_idx = []
    SLC_list = inps.SLC_list
    for i in range(inps.NSLC-n):
        ifgram = []
        flag = True
        for j in range(n):
            ifgram.append('{}_{}'.format(inps.SLC_list[i+j],inps.SLC_list[i+j+1]))
        ifgram.append('{}_{}'.format(SLC_list[i],SLC_list[i+n]))
        for ifgram_name in ifgram:
            if ifgram_name not in inps.date12_list_all:
                flag = False
        if flag:
            cp_idx.append([inps.date12_list_all.index(ifgram[j]) for j in range(n+1)])

    cp_idx = np.array(cp_idx, np.int16)
    cp_idx = np.unique(cp_idx, axis = 0)
    
    num_cp = len(cp_idx)
    print('Number of closure measurements expected, ', len(SLC_list)-n)
    print('Number of closure measurements found, ', num_cp)
    
    if num_cp <1:
        print('No closure phase measurements found, abort')
        sys.exit()
        
    box_width  = box[2] - box[0]
    box_length = box[3] - box[1]
    phase = readfile.read(inps.ifgram_stack_file, box=box,print_msg=False)[0]
    cum_cp = np.zeros((box_length, box_width), np.complex64)
    for i in range(num_cp):
        cp0_w = np.zeros ((box_length, box_width), np.float32)
        for j in range(n):
                    idx = cp_idx[i,j]
                    cp0_w = cp0_w + phase[idx,:,:] - inps.ref_phase[idx]
        idx = cp_idx[i,n]
        cp0_w = cp0_w - (phase[idx,:,:]-inps.ref_phase[idx])
        cum_cp = cum_cp + (np.exp(1j*cp0_w))
        
    # cum_cp = np.angle(cum_cp)
    return cum_cp, num_cp


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    stack_obj = ifgramStack(inps.ifgram_stack_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    date12_list_all = stack_obj.get_date12_list(dropIfgram=False)
    print('scene length, width', length, width)
    ref_phase = stack_obj.get_reference_phase(unwDatasetName = 'unwrapPhase')
    inps.ref_phase = ref_phase
    inps.length = length
    inps.width = width

    # retrieve the list of SLC dates from ifgramStack.h5
    SLC_list = stack_obj.get_date_list(dropIfgram=True)
    print('number of SLC found : ', len(SLC_list))
    print('first SLC: ', SLC_list[0])
    print('last  SLC: ', SLC_list[-1])
    inps.SLC_list = SLC_list
    inps.NSLC = len(SLC_list)
    inps.date12_list = date12_list
    inps.date12_list_all = date12_list_all

    # split igram_file into blocks to save memory
    box_list, num_box = ifginv.split2boxes(inps.ifgram_stack_file,inps.max_memory)
    closure_phase =  np.zeros([length,width],np.complex64)
    #process block-by-block
    for i, box in enumerate(box_list):
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        print(box)
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_width))
            print('box length: {}'.format(box_length))
        closure_phase[box[1]:box[3],box[0]:box[2]], numcp = cum_seq_closure_phase(inps,inps.nl,box)


    # What is a good threadshod?
    # Assume that it's pure noise so that the phase is uniform distributed from -pi to pi.
    # The standard deviation of phase in each loop is pi/sqrt(3)
    #   (technically should be smaller because when forming loops there should be a reduction in phase variance)
    # The standard deviation of phase in cumulative wrapped closure phase is pi/sqrt(3)/sqrt(numcp)
    #   -- again another simplification assuming no correlation.
    # We use 2\delta as threashold -- 95.4% confidence

    threashold_cp = np.pi/np.sqrt(3)/np.sqrt(numcp)*3   # 3/sigma, 99.7% confidence
    mask = np.ones([length,width],np.float32)
    # this masks areas with potential bias
    mask[np.abs(np.angle(closure_phase))>threashold_cp] = 0
    # this unmasks areas with low correlation (where it's hard to know wheter there is bias either)
    mask[np.abs(np.abs(closure_phase)/numcp < 0.3)] = 1
    
    # save mask
    meta = dict(stack_obj.metadata)
    meta['FILE_TYPE'] = 'mask'
    ds_name_dict = {'cpmask': [np.float32, (length, width), mask],}
    writefile.layout_hdf5(os.path.join(inps.outdir,'cpmask.h5'), ds_name_dict, meta)

    # also save the average closure phase
    ds_name_dict2 = {'phase': [np.float32, (length, width), np.angle(closure_phase)],
                    'amplitude':[np.float32,(length,width),np.abs(closure_phase)/numcp],}
    writefile.layout_hdf5(os.path.join(inps.outdir,'avgwcp.h5'), ds_name_dict2, meta)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
