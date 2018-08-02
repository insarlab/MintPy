#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import time
import argparse
import itertools
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import sparse, ndimage
from pysar.objects import ifgramStack
from pysar.utils import (ptime,
                         readfile,
                         writefile,
                         utils as ut,
                         plot as pp,
                         deramp)
from pysar import ifgram_inversion as ifginv
from pysar.unwrap_error_phase_closure import run_unwrap_error_closure


####################################################################################################
EXAMPLE = """Example:
  unwrap_error_hybrid.py  ./INPUTS/ifgramStack.h5  -m waterMask.h5
"""

REFERENCE = """Reference:
  Yunjun Z., H. Fattahi, F. Amelung (2018), A Hybrid Unwrapping Error Correction Method based on Phase Closure
  and Bridging for a Network of Interferograms. (in prep.)
"""

def create_parser():
    parser = argparse.ArgumentParser(description=('Unwrapping Error Correction '
                                                  'based on Phase Closure and Connected Components'),
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', type=str, help='interferograms file to be corrected')
    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help='name of dataset to be corrected.')
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction')

    parser.add_argument('-m','--mask', dest='maskFile', type=str, required=True,
                        help='name of mask file for pixels to be corrected, e.g. waterMask.h5')
    parser.add_argument('--ramp', dest='ramp', choices=['plane', 'quadratic'], 
                          help='type of phase ramp to be removed before correction.')
    parser.add_argument('-r','--radius', dest='bridge_end_radius', type=int, default=150,
                        help='radius of the end point of bridge to search area to get median representative value')
    parser.add_argument('--cutoff', dest='cutoff', type=float, default=1.,
                        help='cutoff value to detect histogram spikes for coherent connected components')
    parser.add_argument('--update', action='store_true',
                        help='Enable update mode: if unwrapPhase_unwCor dataset exists, skip the correction.')
    parser.add_argument('--no-closure-fowllow-on', dest='run_closure', action='store_false',
                        help=('Do not run unwrap error correction based on phase closure after briding,'
                              ' to correct the remaining unwrap errors.'))
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def get_nonzero_phase_closure(ifgram_file, out_file=None, thres=0.1, unwDatasetName='unwrapPhase'):
    """Calculate/Read number of non-zero phase closure
    Parameters: ifgram_file : string, path of ifgram stack file
                out_file    : string, path of num non-zero phase closure file
    Returns:    num_nonzero_closure : 2D np.array in size of (length, width)
    """
    if not out_file:
        out_file = 'numNonzeroClosure_{}.h5'.format(unwDatasetName)
    if os.path.isfile(out_file) and readfile.read_attribute(out_file):
        print('1. read number of nonzero phase closure from file: {}'.format(out_file))
        num_nonzero_closure = readfile.read(out_file)[0]
    else:
        stack_obj = ifgramStack(ifgram_file)
        stack_obj.open(print_msg=False)
        length, width = stack_obj.length, stack_obj.width

        ref_phase = stack_obj.get_reference_phase(unwDatasetName=unwDatasetName, dropIfgram=False)
        C = stack_obj.get_design_matrix4ifgram_triangle(dropIfgram=False)

        # calculate phase closure line by line to save memory usage
        num_nonzero_closure = np.zeros((length, width), np.float32)
        print('1. calculating phase closure of all pixels from dataset - {} ...'.format(unwDatasetName))
        line_step = 10
        num_loop = int(np.ceil(length / line_step))
        prog_bar = ptime.progressBar(maxValue=num_loop)
        for i in range(num_loop):
            # read phase
            i0, i1 = i*line_step, min(length, (i+1)*line_step)
            box = (0, i0, width, i1)
            pha_data = ifginv.read_unwrap_phase(stack_obj,
                                                box,
                                                ref_phase,
                                                unwDatasetName=unwDatasetName,
                                                dropIfgram=False,
                                                print_msg=False)
            # calculate phase closure
            pha_closure = np.dot(C, pha_data)
            pha_closure = np.abs(pha_closure - ut.wrap(pha_closure))
            # get number of non-zero phase closure
            num_nonzero = np.sum(pha_closure >= thres, axis=0)
            num_nonzero_closure[i0:i1, :] = num_nonzero.reshape(i1-i0, width)
            prog_bar.update(i+1, every=1, suffix='{}/{} lines'.format((i+1)*line_step, length))
        prog_bar.close()

        atr = dict(stack_obj.metadata)
        atr['FILE_TYPE'] = 'mask'
        atr['UNIT'] = 1
        writefile.write(num_nonzero_closure, out_file=out_file, metadata=atr)
    return num_nonzero_closure


def detect_unwrap_error(ifgram_file, mask_file, mask_cc_file='maskConnComp.h5', unwDatasetName='unwrapPhase',
                        cutoff=1., min_num_pixel=1e4):
    """Detect unwrapping error based on phase closure and extract coherent conn comps
    based on its histogram distribution
    Parameters: ifgram_file : string, path of ifgram stack file
                mask_file   : string, path of mask file, e.g. waterMask.h5, mask.h5
                mask_cc_file: string, path of mask file for coherent conn comps
                cutoff : float, cutoff value for the mean number of nonzero phase closure
                    to be selected as coherent conn comps candidate
                min_num_pixel : float, min number of pixels left after morphology operation
                    to be determined as coherent conn comps
    Returns:    mask_cc_file : string, path of mask file for coherent conn comps
    """
    print('-'*50)
    print('detect unwraping error based on phase closure')
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    C = stack_obj.get_design_matrix4ifgram_triangle(dropIfgram=False)

    num_nonzero_closure = get_nonzero_phase_closure(ifgram_file, unwDatasetName=unwDatasetName)

    # get histogram of num_nonzero_phase_closure
    mask = readfile.read(mask_file)[0]
    mask *= num_nonzero_closure != 0.

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 4])
    num4disp = np.array(num_nonzero_closure)
    num4disp[mask == 0] = np.nan
    im = ax[0].imshow(num4disp)
    ax[0] = pp.auto_flip_direction(stack_obj.metadata, ax=ax[0], print_msg=False)
    cbar = fig.colorbar(im, ax=ax[0])
    cbar.set_label('number of non-zero phase closure')

    print('2. extract coherent conn comps with unwrap error based on histogram distribution')
    max_nonzero_closure = int(np.max(num_nonzero_closure[mask]))
    bin_value, bin_edge = ax[1].hist(num_nonzero_closure[mask].flatten(),
                                     range=(0, max_nonzero_closure),
                                     log=True,
                                     bins=max_nonzero_closure)[0:2]
    ax[1].set_xlabel('number of non-zero phase closure')
    ax[1].set_ylabel('number of pixels')

    if 'Closure' not in unwDatasetName:
        print('eliminate pixels with number of nonzero phase closure < 5% of total phase closure number')
        print('\twhich can be corrected using phase closure alone.')
        bin_value[:int(C.shape[0]*0.05)] = 0.
    bin_value_thres = ut.median_abs_deviation_threshold(bin_value, cutoff=cutoff)
    print('median abs deviation cutoff value: {}'.format(cutoff))

    plt.plot([0, max_nonzero_closure], [bin_value_thres, bin_value_thres])
    out_img = 'numUnwErr_stat.png'
    fig.savefig(out_img, bbox_inches='tight', transparent=True, dpi=300)
    print('save unwrap error detection result to {}'.format(out_img))

    # histogram --> candidates of coherence conn comps --> mask_cc
    # find pixel clusters sharing similar number of non-zero phase closure
    print('searching connected components with more than {} pixels'.format(min_num_pixel))
    bin_label, n_bins = ndimage.label(bin_value > bin_value_thres)

    mask_cc = np.zeros(num_nonzero_closure.shape, dtype=np.int16)
    # first conn comp - reference conn comp with zero non-zero phase closure
    num_cc = 1
    mask_cc1 = num_nonzero_closure == 0.
    mask_cc1s = ut.get_all_conn_components(mask_cc1, min_num_pixel=min_num_pixel)
    for mask_cc1 in mask_cc1s:
        mask_cc += mask_cc1

    # other conn comps - target conn comps to be corrected for unwrap error
    for i in range(n_bins):
        idx = np.where(bin_label == i+1)[0]
        mask_cci0 = np.multiply(num_nonzero_closure >= bin_edge[idx[0]],
                                num_nonzero_closure <  bin_edge[idx[-1]+1])
        mask_ccis = ut.get_all_conn_components(mask_cci0, min_num_pixel=min_num_pixel)
        if mask_ccis:
            for mask_cci in mask_ccis:
                num_cc += 1
                mask_cc += mask_cci * num_cc
    
                fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[8, 4])
                im = ax[0].imshow(mask_cci0)
                im = ax[1].imshow(mask_cci)
                fig.savefig('mask_cc{}.png'.format(num_cc),
                            bbox_inches='tight', transparent=True, dpi=300)

    # save to hdf5 file
    num_bridge = num_cc - 1
    atr = dict(stack_obj.metadata)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 1
    writefile.write(mask_cc, out_file=mask_cc_file, metadata=atr)

    # plot and save figure to img file
    out_img = '{}.png'.format(os.path.splitext(mask_cc_file)[0])
    fig, ax = plt.subplots(figsize=[6, 8])
    im = ax.imshow(mask_cc)
    ax = pp.auto_flip_direction(atr, ax=ax, print_msg=False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax, ticks=np.arange(num_bridge+2))
    fig.savefig(out_img, bbox_inches='tight', transparent=True, dpi=300)
    print('save to {}'.format(out_img))

    return mask_cc_file


def search_bridge(mask_cc_file, radius=150, display=False):
    """Search bridges to connect coherent conn comps with min distance
    Parameters: mask_cc_file : string, path of mask file of coherent conn comps
                radius : int, radius of influence for median value calculation to represent the value of cc
    Returns:    bridges : list of dict for bridges connecting pair of conn comps with each dict contains:
                    'x0' : x coordinate of reference point
                    'y0' : y coordinate of reference point
                    'mask0' : 2D np.array in size of (length, width) in np.bool_
                         represent the circular mask used to calculate the representative value in reference conn comp
                    'x1' : x coordinate of target point
                    'y1' : y coordinate of target point
                    'mask1' : 2D np.array in size of (length, width) in np.bool_
                         represent the circular mask used to calculate the representative value in target conn comp
    """
    print('-'*50)
    print('searching bridges to connect coherence conn comps ...')
    mask_cc = readfile.read(mask_cc_file)[0]
    num_bridge = int(np.max(mask_cc) - 1)
    shape = mask_cc.shape

    # calculate min distance between each pair of conn comps and its corresponding bonding points
    print('1. calculating min distance between each pair of coherent conn comps ...')
    connections = {}
    dist_mat = np.zeros((num_bridge+1, num_bridge+1), dtype=np.float32)
    for i, j in itertools.combinations(range(1, num_bridge+2), 2):
        mask1 = mask_cc == i
        mask2 = mask_cc == j
        pts1, pts2, d = ut.min_region_distance(mask1, mask2, display=False)
        dist_mat[i-1, j-1] = dist_mat[j-1, i-1] = d
        conn_dict = dict()
        conn_dict['{}'.format(i)] = pts1
        conn_dict['{}'.format(j)] = pts2
        conn_dict['distance'] = d
        connections['{}{}'.format(i, j)] = conn_dict
        print('conn comp pair: {} {} to {} {} with distance of {}'.format(i, pts1, j, pts2, d))

    # 1. calculate the min-spanning-tree path to connect all conn comps
    # 2. find the order using a breadth-first ordering starting with conn comp 1 
    print('2. search bridging order using breadth-first ordering')
    dist_mat_mst = sparse.csgraph.minimum_spanning_tree(dist_mat)
    nodes, predecessors = sparse.csgraph.breadth_first_order(dist_mat_mst,
                                                             i_start=0,
                                                             directed=False)
    nodes += 1
    predecessors += 1

    # converting bridging order into bridges
    print('3. find circular area around bridge point')
    bridges = []
    for i in range(num_bridge):
        # get x0/y0/x1/y1
        n0, n1 = predecessors[i+1], nodes[i+1]
        conn_dict = connections['{}{}'.format(n0, n1)]
        x0, y0 = conn_dict[str(n0)]
        x1, y1 = conn_dict[str(n1)]

        # get mask0/mask1
        mask0 = (mask_cc == n0) * ut.get_circular_mask(x0, y0, radius, shape)
        mask1 = (mask_cc == n1) * ut.get_circular_mask(x1, y1, radius, shape)

        # save to list of dict
        bridge = dict()
        bridge['x0'] = x0
        bridge['y0'] = y0
        bridge['x1'] = x1
        bridge['y1'] = y1
        bridge['mask0'] = mask0
        bridge['mask1'] = mask1
        bridges.append(bridge)

    plot_bridges(mask_cc_file, bridges, display=display)
    return bridges


def plot_bridges(mask_cc_file, bridges, display=False):
    """Plot mask of connected components with bridges info
    Parameters: mask_cc_file : string, path of mask cc file
                bridges : list of dict
                display : bool
    """
    out_base='maskConnCompBridge'
    mask_cc, metadata = readfile.read(mask_cc_file)

    # check number of bridges
    num_bridge = len(bridges)

    # save to text file
    out_file = '{}.txt'.format(out_base)
    bridge_yx = np.zeros((num_bridge, 4), dtype=np.int16)
    for i in range(num_bridge):
        bridge = bridges[i]
        bridge_yx[i, :] = [bridge['y0'], bridge['x0'], bridge['y1'], bridge['x1']]
    header_info = 'bridge bonding points for unwrap error correction\n'
    header_info += 'y0\tx0\ty1\tx1'
    np.savetxt(out_file, bridge_yx, fmt='%s', delimiter='\t', header=header_info)
    print('save bridge points to file: {}'.format(out_file))

    # plot/save to image file
    out_file = '{}.png'.format(out_base)
    fig, ax = plt.subplots(figsize=[6, 8])

    # plot mask_cc data
    im = ax.imshow(mask_cc)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax, ticks=np.arange(num_bridge+2))

    # plot bridge data
    for i in range(num_bridge):
        bridge = bridges[i]
        ax.imshow(np.ma.masked_where(~bridge['mask0'], np.zeros(mask_cc.shape)),
                  cmap='gray', alpha=0.3, vmin=0, vmax=1)
        ax.imshow(np.ma.masked_where(~bridge['mask1'], np.zeros(mask_cc.shape)),
                  cmap='gray', alpha=0.3, vmin=0, vmax=1)
        ax.plot([bridge['x0'], bridge['x1']],
                [bridge['y0'], bridge['y1']],
                '-', ms=5, mfc='none')

    ax = pp.auto_flip_direction(metadata, ax=ax, print_msg=False)
    fig.savefig(out_file, bbox_inches='tight', transparent=True, dpi=300)
    print('plot/save bridge setting to file: {}'.format(out_file))
    if display:
        print('showing')
        plt.show()
    return


##########################################################################################
def bridge_unwrap_error(data, mask_cc, bridges):
    """Phase Jump Correction, using phase continuity on bridge/bonding points in each pair of conn comps
    Inputs:
        data : 2D np.array in size of (length, width), phase matrix need to be corrected
        mask_cc : 2D np.array in size of (length, width), mask of different connected components
        bridges : list of dict for bridges connecting pair of conn comps with each dict contains:
                    'x0' : x coordinate of reference point
                    'y0' : y coordinate of reference point
                    'mask0' : 2D np.array in size of (length, width) in np.bool_
                         represent the circular mask used to calculate the representative value in reference conn comps
                    'x1' : x coordinate of target point
                    'y1' : y coordinate of target point
                    'mask1' : 2D np.array in size of (length, width) in np.bool_
                         represent the circular mask used to calculate the representative value in target conn comps
    Output:
        data : 2D np.array, phase corrected matrix
    """
    data = np.array(data, dtype=np.float32)
    for i in range(len(bridges)):
        # get phase difference between two ends of the bridge
        bridge = bridges[i]
        value0 = np.nanmedian(data[bridge['mask0']])
        value1 = np.nanmedian(data[bridge['mask1']])
        diff_value = value1 - value0

        #estimate integer number of phase jump
        num_jump = (np.abs(diff_value) + np.pi) // (2.*np.pi)
        if diff_value > 0:
            num_jump *= -1

        # correct unwrap error to target conn comps
        mask1 = mask_cc == mask_cc[bridge['y1'], bridge['x1']]
        data[mask1] += num_jump * 2. * np.pi

    return data


def run_unwrap_error_bridge(inps, mask_cc_file, bridges, dsNameIn='unwrapPhase',
                            dsNameOut='unwrapPhase_hybrid'):
    print('-'*50)
    print('correct unwrapping error in {} with bridging ...'.format(inps.ifgram_file))
    stack_obj = ifgramStack(inps.ifgram_file)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=False)
    num_ifgram = len(date12_list)
    shape_out = (num_ifgram, stack_obj.length, stack_obj.width)

    print('read mask from file: {}'.format(mask_cc_file))
    mask_cc = readfile.read(mask_cc_file)[0]
    if inps.ramp is not None:
        print('estimate and remove phase ramp of {} during the correction'.format(inps.ramp))
        ramp_mask = (mask_cc == mask_cc[stack_obj.refY, stack_obj.refX])

    # prepare output data writing
    print('open {} with r+ mode'.format(inps.ifgram_file))
    f = h5py.File(inps.ifgram_file, 'r+')
    if dsNameOut in f.keys():
        ds = f[dsNameOut]
        print('write to /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))
    else:
        ds = f.create_dataset(dsNameOut, shape_out, maxshape=(None, None, None),
                              chunks=True, compression=None)
        print('create /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))

    # correct unwrap error ifgram by ifgram
    prog_bar = ptime.progressBar(maxValue=num_ifgram)
    for i in range(num_ifgram):
        # read unwrapPhase
        date12 = date12_list[i]
        unw = np.squeeze(f[dsNameIn][i, :, :])
        unw -= unw[stack_obj.refY, stack_obj.refX]

        # remove phase ramp before phase jump estimation
        if inps.ramp is not None:
            unw, unw_ramp = deramp.remove_data_surface(unw, ramp_mask, inps.ramp)

        # estimate/correct phase jump
        unw_cor = bridge_unwrap_error(unw, mask_cc, bridges)
        if inps.ramp is not None:
            unw_cor += unw_ramp

        # write to hdf5 file
        ds[i, :, :] = unw_cor
        prog_bar.update(i+1, suffix=date12)
    prog_bar.close()
    f.close()
    print('close {} file.'.format(inps.ifgram_file))
    return inps.ifgram_file


####################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_hybrid'.format(inps.datasetNameIn)

    # update mode checking
    atr = readfile.read_attribute(inps.ifgram_file)
    if inps.update and atr['FILE_TYPE'] == 'ifgramStack':
        stack_obj = ifgramStack(inps.ifgram_file)
        stack_obj.open(print_msg=False)
        if inps.datasetNameOut in stack_obj.datasetNames:
            print(("update mode is enabled AND {} already exists"
                   " skip this step.").format(inps.datasetNameOut))
            return inps.ifgram_file

    start_time = time.time()

    mask_cc_file = detect_unwrap_error(ifgram_file=inps.ifgram_file,
                                       mask_file=inps.maskFile,
                                       mask_cc_file='maskConnComp.h5',
                                       unwDatasetName=inps.datasetNameIn,
                                       cutoff=inps.cutoff)

    bridges = search_bridge(mask_cc_file, radius=inps.bridge_end_radius)

    run_unwrap_error_bridge(inps, mask_cc_file, bridges,
                            dsNameIn=inps.datasetNameIn,
                            dsNameOut=inps.datasetNameOut)

    if inps.run_closure:
        print('')
        inps.datasetNameIn = inps.datasetNameOut
        inps.datasetNameOut = '{}_closure'.format(inps.datasetNameIn)
        run_unwrap_error_closure(inps,
                                 dsNameIn=inps.datasetNameIn, 
                                 dsNameOut=inps.datasetNameOut,
                                 fast_mode=True)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
