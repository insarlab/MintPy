#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import sparse, ndimage
from pysar.objects import ifgramStack
from pysar.utils import (ptime,
                         readfile,
                         writefile,
                         utils as ut,
                         plot as pp)
from pysar import ifgram_inversion as ifginv
from pysar.unwrap_error_phase_closure import run_unwrap_error_closure
from pysar.unwrap_error_bridging import search_bridge, run_unwrap_error_bridge


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
    parser.add_argument('--ramp', dest='ramp', choices=['linear', 'quadratic'], 
                          help='type of phase ramp to be removed before correction.')
    parser.add_argument('-r','--radius', dest='bridgePtsRadius', type=int, default=150,
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
        out_file = 'numNonzeroPhaseClosure_{}.h5'.format(unwDatasetName)
    if os.path.isfile(out_file) and readfile.read_attribute(out_file):
        print('1. read number of nonzero phase closure from file: {}'.format(out_file))
        num_nonzero_closure = readfile.read(out_file)[0]
    else:
        obj = ifgramStack(ifgram_file)
        obj.open(print_msg=False)
        length, width = obj.length, obj.width

        ref_phase = obj.get_reference_phase(unwDatasetName=unwDatasetName, dropIfgram=False)
        C = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=False))

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
            pha_data = ifginv.read_unwrap_phase(obj,
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

        atr = dict(obj.metadata)
        atr['FILE_TYPE'] = 'mask'
        atr['UNIT'] = 1
        writefile.write(num_nonzero_closure, out_file=out_file, metadata=atr)
    return num_nonzero_closure


def detect_unwrap_error(ifgram_file, mask_file, mask_cc_file='maskConnComp.h5', unwDatasetName='unwrapPhase',
                        cutoff=1., min_num_pixel=1e4):
    """Detect unwrapping error based on phase closure and extract coherent conn comps
    based on its histogram distribution

    Check:
    https://en.wikipedia.org/wiki/Otsu%27s_method
    from skimage.filters import threshold_otsu
    
    Parameters: ifgram_file : string, path of ifgram stack file
                mask_file   : string, path of mask file, e.g. waterMask.h5, maskConnComp.h5
                mask_cc_file: string, path of mask file for coherent conn comps
                cutoff : float, cutoff value for the mean number of nonzero phase closure
                    to be selected as coherent conn comps candidate
                min_num_pixel : float, min number of pixels left after morphology operation
                    to be determined as coherent conn comps
    Returns:    mask_cc_file : string, path of mask file for coherent conn comps
    """
    print('-'*50)
    print('detect unwraping error based on phase closure')
    obj = ifgramStack(ifgram_file)
    obj.open(print_msg=False)
    C = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=False))

    num_nonzero_closure = get_nonzero_phase_closure(ifgram_file, unwDatasetName=unwDatasetName)

    # get histogram of num_nonzero_phase_closure
    mask = readfile.read(mask_file)[0]
    mask *= num_nonzero_closure != 0.

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 4])
    num4disp = np.array(num_nonzero_closure, dtype=np.float32)
    num4disp[mask == 0] = np.nan
    im = ax[0].imshow(num4disp)
    ax[0].set_xlabel('Range [pix.]')
    ax[0].set_ylabel('Azimuth [pix.]')
    ax[0] = pp.auto_flip_direction(obj.metadata, ax=ax[0], print_msg=False)
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
    atr = dict(obj.metadata)
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


####################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_bridging'.format(inps.datasetNameIn)

    # update mode checking
    atr = readfile.read_attribute(inps.ifgram_file)
    if inps.update and atr['FILE_TYPE'] == 'ifgramStack':
        obj = ifgramStack(inps.ifgram_file)
        obj.open(print_msg=False)
        if inps.datasetNameOut in obj.datasetNames:
            print(("update mode is enabled AND {} already exists"
                   " skip this step.").format(inps.datasetNameOut))
            return inps.ifgram_file

    start_time = time.time()
    # get mask_cc from phase closure
    mask_cc_file = detect_unwrap_error(ifgram_file=inps.ifgram_file,
                                       mask_file=inps.maskFile,
                                       mask_cc_file='maskConnComp.h5',
                                       unwDatasetName=inps.datasetNameIn,
                                       cutoff=inps.cutoff)
    # run bridging
    bridges = search_bridge(mask_cc_file, radius=inps.bridgePtsRadius)
    run_unwrap_error_bridge(inps.ifgram_file,
                            mask_cc_file,
                            bridges,
                            dsNameIn=inps.datasetNameIn,
                            dsNameOut=inps.datasetNameOut,
                            ramp_type=inps.ramp)

    if inps.run_closure:
        print('')
        inps.datasetNameIn = inps.datasetNameOut
        inps.datasetNameOut = '{}_phaseClosure'.format(inps.datasetNameIn)
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
