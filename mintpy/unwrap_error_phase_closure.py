#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import argparse
import time
import h5py
import numpy as np
from matplotlib import pyplot as plt, ticker

try:
    from cvxopt import matrix
except ImportError:
    raise ImportError('Cannot import cvxopt')
try:
    from skimage import measure
except ImportError:
    raise ImportError('Could not import skimage!')

from mintpy.objects import ifgramStack, conncomp
from mintpy.defaults.template import get_template_content
from mintpy.utils import ptime, readfile, writefile, utils as ut, plot as pp
from mintpy.utils.solvers import l1regls
from mintpy import ifgram_inversion as ifginv


key_prefix = 'mintpy.unwrapError.'

##########################################################################################
EXAMPLE = """example:
  # correct phase unwrapping error with phase closure
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --cc-mask maskConnComp.h5  -t smallbaselineApp.cfg   --update
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --cc-mask maskConnComp.h5  --water-mask waterMask.h5 --update

  # calculate the number of non-zero closure phase
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --action calculate
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --action calculate  --water-mask waterMask.h5
"""

NOTE = """
  by exploiting the conservertiveness of the integer ambiguity of interferograms triplets.
  This method assumes:
  a. abundance of network: for interferogram with unwrapping error, there is
     at least of one triangular connection to form a closed circle; with more
     closed circles comes better constrain.
  b. majority rightness: most of interferograms have to be right (no unwrapping
     error) to correct the wrong minority. And if most of interferograms have 
     unwrapping errors, then the minor right interferograms will turn into wrong.
"""

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis:
  Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
  doi:10.1016/j.cageo.2019.104331.
"""

TEMPLATE1 = get_template_content('quick_overview')
TEMPLATE2 = get_template_content('correct_unwrap_error')



def create_parser():
    parser = argparse.ArgumentParser(description='Unwrapping Error Correction based on Phase Closure'+NOTE,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+TEMPLATE1+'\n'+TEMPLATE2+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be corrected')
    parser.add_argument('-c','--cc-mask', dest='cc_mask_file', default='maskConnComp.h5',
                        help='common connected components file, required for --action correct')
    parser.add_argument('-n','--num-sample', dest='numSample', type=int, default=100,
                        help='Number of randomly samples/pixels for each common connected component.')

    parser.add_argument('-a','--action', dest='action', type=str, default='correct',
                        choices={'correct', 'calculate'},
                        help='action to take (default: %(default)s):\n'+
                             'correct   - correct phase unwrapping error\n'+
                             'calculate - calculate the number of non-zero closure phase')

    # IO
    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help="name of dataset to be corrected, default: unwrapPhase")
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_phaseClosure')

    # mask
    mask = parser.add_argument_group('mask')
    mask.add_argument('--water-mask','--wm', dest='waterMaskFile', type=str,
                      help='path of water mask file.')
    mask.add_argument('-t', '--template', dest='template_file',
                      help='template file with options for setting.')

    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_phaseClosure dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # backend setting for matplotlib
    plt.switch_backend('Agg')

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check 1 input file type
    k = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if k not in ['ifgramStack']:
        raise ValueError('input file is not ifgramStack: {}'.format(k))

    # check 2 cc_mask_file
    if inps.action == 'correct' and not os.path.isfile(inps.cc_mask_file):
        raise FileNotFoundError(inps.cc_mask_file)

    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_phaseClosure'.format(inps.datasetNameIn)

    # discard water mask file is not found
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from tempalte file: '+os.path.basename(inps.template_file))
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if value:
            if key in ['waterMaskFile']:
                inpsDict[key] = value
            elif key in ['numSample']:
                inpsDict[key] = int(value)

    return inps


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output dataset
    with h5py.File(inps.ifgram_file, 'r') as f:
        if inps.datasetNameOut not in f.keys():
            flag = 'run'
            print('1) output dataset: {} NOT found.'.format(inps.datasetNameOut))
        else:
            print('1) output dataset: {} exists.'.format(inps.datasetNameOut))
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if ti > to:
                flag = 'run'
                print('2) output dataset is NOT newer than input dataset: {}.'.format(inps.datasetNameIn))
            else:
                print('2) output dataset is newer than input dataset: {}.'.format(inps.datasetNameIn))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


##########################################################################################
def calc_num_triplet_with_nonzero_integer_ambiguity(ifgram_file, mask_file=None, dsName='unwrapPhase',
                                                    out_file=None, max_memory=4, update_mode=True):
    """Calculate the number of triplets with non-zero integer ambiguity of closure phase.

    T_int as shown in equation (8-9) and inline in Yunjun et al. (2019, CAGEO).

    Parameters: ifgram_file - str, path of interferogram stack file
                mask_file   - str, path of mask file
                dsName      - str, unwrapped phase dataset name used to calculate the closure phase
                out_file    - str, custom output filename
                update_mode - bool
    Returns:    out_file    - str, custom output filename
    Example:    calc_num_triplet_with_nonzero_integer_ambiguity('inputs/ifgramStack.h5', mask_file='waterMask.h5')
    """

    # default output file path
    out_dir = os.path.dirname(os.path.dirname(ifgram_file))
    if out_file is None:
        if dsName == 'unwrapPhase':
            # skip the default dsName in output filename
            out_file = 'numTriNonzeroIntAmbiguity.h5'
        else:
            out_file = 'numTriNonzeroIntAmbiguity4{}.h5'.format(dsName)
        out_file = os.path.join(out_dir, out_file)

    # update mode
    if update_mode and os.path.isfile(out_file):
        print('update mode: ON')
        print('1) output file "{}" already exists'.format(out_file))
        flag = 'skip'

        # check modification time
        with h5py.File(ifgram_file, 'r') as f:
            ti = float(f[dsName].attrs.get('MODIFICATION_TIME', os.path.getmtime(ifgram_file)))
        to = os.path.getmtime(out_file)
        if ti > to:
            print('2) output file is NOT newer than input dataset')
            flag = 'run'
        else:
            print('2) output file is newer than input dataset')

        # check REF_Y/X
        key_list = ['REF_Y', 'REF_X']
        atri = readfile.read_attribute(ifgram_file)
        atro = readfile.read_attribute(out_file)
        if not all(atri[i] == atro[i] for i in key_list):
            print('3) NOT all key configurations are the same: {}'.format(key_list))
            flag = 'run'
        else:
            print('3) all key configurations are the same: {}'.format(key_list))

        print('run or skip: {}.'.format(flag))
        if flag == 'skip':
            return out_file

    # read ifgramStack file
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    num_ifgram = len(date12_list)

    C = stack_obj.get_design_matrix4triplet(date12_list)
    if C is None:
        msg = 'No triangles found from ifgramStack file: {}!'.format(ifgram_file)
        msg += '\n    Skip calculating the number of triplets with non-zero integer ambiguity.'
        print(msg)
        return None
    else:
        print('get design matrix for the interferogram triplets in size of {}'.format(C.shape))

    # calculate number of nonzero closure phase
    ds_size = (C.shape[0] * 2 + C.shape[1]) * length * width * 4
    num_loop = int(np.ceil(ds_size * 2 / (max_memory * 1024**3)))
    step = int(np.rint(length / num_loop / 10) * 10)
    num_loop = int(np.ceil(length / step))
    num_nonzero_closure = np.zeros((length, width), dtype=np.float32)
    msg = 'calculating the number of triplets with non-zero integer ambiguity of closure phase ...'
    msg += '\n    block by block with size up to {}, {} blocks in total'.format((step, width), num_loop)
    print(msg)

    ref_phase = stack_obj.get_reference_phase(unwDatasetName=dsName, dropIfgram=True).reshape(num_ifgram, -1)
    prog_bar = ptime.progressBar(maxValue=num_loop)
    for i in range(num_loop):
        # box
        r0 = i * step
        r1 = min((r0+step), stack_obj.length)
        box = (0, r0, stack_obj.width, r1)

        # read data
        unw = ifginv.read_unwrap_phase(stack_obj,
                                       box=box,
                                       ref_phase=ref_phase,
                                       obs_ds_name=dsName,
                                       dropIfgram=True,
                                       print_msg=False).reshape(num_ifgram, -1)

        # calculate based on equation (8-9) and T_int equation inline.
        closure_pha = np.dot(C, unw)
        closure_int = np.round((closure_pha - ut.wrap(closure_pha)) / (2.*np.pi))
        num_nonzero_closure[r0:r1, :] = np.sum(closure_int != 0, axis=0).reshape(-1, width)

        prog_bar.update(i+1, every=1, suffix='line {} / {}'.format(r0, length))
    prog_bar.close()

    # mask
    if mask_file is not None:
        mask = readfile.read(mask_file)[0]
        num_nonzero_closure[mask == 0] = np.nan
        print('mask out pixels with zero in file:', mask_file)

    coh_file = os.path.join(out_dir, 'avgSpatialCoh.h5')
    if os.path.isfile(coh_file):
        coh = readfile.read(coh_file)[0]
        num_nonzero_closure[coh == 0] = np.nan
        print('mask out pixels with zero in file:', coh_file)

    # write to disk
    print('write to file', out_file)
    meta = dict(stack_obj.metadata)
    meta['FILE_TYPE'] = 'mask'
    meta['UNIT'] = '1'
    writefile.write(num_nonzero_closure, out_file, meta)

    # plot
    plot_num_triplet_with_nonzero_integer_ambiguity(out_file)

    return out_file


def plot_num_triplet_with_nonzero_integer_ambiguity(fname, display=False, font_size=12, fig_size=[9,3]):
    """Plot the histogram for the number of triplets with non-zero integer ambiguity

    Fig. 3d-e in Yunjun et al. (2019, CAGEO).
    """

    # read data
    data, atr = readfile.read(fname)
    vmax = int(np.nanmax(data))

    # plot
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=fig_size)

    # subplot 1 - map
    ax = axs[0]
    im = ax.imshow(data, cmap='RdBu_r', interpolation='nearest')

    # reference point
    if all(key in atr.keys() for key in ['REF_Y','REF_X']):
        ax.plot(int(atr['REF_X']), int(atr['REF_Y']), 's', color='white', ms=3)

    # format
    pp.auto_flip_direction(atr, ax=ax, print_msg=False)
    fig.colorbar(im, ax=ax)
    ax.set_title(r'$T_{int}$', fontsize=font_size)

    # subplot 2 - histogram
    ax = axs[1]
    ax.hist(data[~np.isnan(data)].flatten(), range=(0, vmax), log=True, bins=vmax)

    # axis format
    ax.set_xlabel(r'# of triplets w non-zero int ambiguity $T_{int}$', fontsize=font_size)
    ax.set_ylabel('# of pixels', fontsize=font_size)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, numticks=15,
                                                 subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    for ax in axs:
        ax.tick_params(which='both', direction='in', labelsize=font_size,
                       bottom=True, top=True, left=True, right=True)

    fig.tight_layout()

    # output
    out_fig = '{}.png'.format(os.path.splitext(fname)[0])
    print('plot and save figure to file', out_fig)
    fig.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)

    if display:
        plt.show()
    else:
        plt.close(fig)

    return


##########################################################################################
def get_common_region_int_ambiguity(ifgram_file, cc_mask_file, water_mask_file=None, num_sample=100,
                                    dsNameIn='unwrapPhase'):
    """Solve the phase unwrapping integer ambiguity for the common regions among all interferograms
    Parameters: ifgram_file     : str, path of interferogram stack file
                cc_mask_file    : str, path of common connected components file
                water_mask_file : str, path of water mask file
                num_sample      : int, number of pixel sampled for each region
                dsNameIn        : str, dataset name of the unwrap phase to be corrected
    Returns:    common_regions  : list of skimage.measure._regionprops._RegionProperties object
                    modified by adding two more variables:
                    sample_coords : 2D np.ndarray in size of (num_sample, 2) in int64 format
                    int_ambiguity : 1D np.ndarray in size of (num_ifgram,) in int format
    """
    print('-'*50)
    print('calculating the integer ambiguity for the common regions defined in', cc_mask_file)
    # stack info
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open()
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    num_ifgram = len(date12_list)
    C = matrix(ifgramStack.get_design_matrix4triplet(date12_list).astype(float))
    ref_phase = stack_obj.get_reference_phase(unwDatasetName=dsNameIn, dropIfgram=True).reshape(num_ifgram, -1)

    # prepare common label
    print('read common mask from', cc_mask_file)
    cc_mask = readfile.read(cc_mask_file)[0]
    if water_mask_file is not None and os.path.isfile(water_mask_file):
        water_mask = readfile.read(water_mask_file)[0]
        print('refine common mask based on water mask file', water_mask_file)
        cc_mask[water_mask == 0] = 0

    label_img, num_label = conncomp.label_conn_comp(cc_mask, min_area=2.5e3, print_msg=True)
    common_regions = measure.regionprops(label_img)
    print('number of common regions:', num_label)

    # add sample_coords / int_ambiguity
    print('number of samples per region:', num_sample)
    print('solving the phase-unwrapping integer ambiguity for {}'.format(dsNameIn))
    print('\tbased on the closure phase of interferograms triplets (Yunjun et al., 2019)')
    print('\tusing the L1-norm regularzed least squares approximation (LASSO) ...')
    for i in range(num_label):
        common_reg = common_regions[i]
        # sample_coords
        idx = sorted(np.random.choice(common_reg.area, num_sample, replace=False))
        common_reg.sample_coords = common_reg.coords[idx, :].astype(int)

        # solve for int_ambiguity
        U = np.zeros((num_ifgram, num_sample))
        if common_reg.label == label_img[stack_obj.refY, stack_obj.refX]:
            print('{}/{} skip calculation for the reference region'.format(i+1, num_label))
        else:
            prog_bar = ptime.progressBar(maxValue=num_sample, prefix='{}/{}'.format(i+1, num_label))
            for j in range(num_sample):
                # read unwrap phase
                y, x = common_reg.sample_coords[j, :]
                unw = ifginv.read_unwrap_phase(stack_obj,
                                               box=(x, y, x+1, y+1),
                                               ref_phase=ref_phase,
                                               obs_ds_name=dsNameIn,
                                               dropIfgram=True,
                                               print_msg=False).reshape(num_ifgram, -1)

                # calculate closure_int
                closure_pha = np.dot(C, unw)
                closure_int = matrix(np.round((closure_pha - ut.wrap(closure_pha)) / (2.*np.pi)))

                # solve for U
                U[:,j] = np.round(l1regls(-C, closure_int, alpha=1e-2, show_progress=0)).flatten()
                prog_bar.update(j+1, every=5)
            prog_bar.close()
        # add int_ambiguity
        common_reg.int_ambiguity = np.median(U, axis=1)
        common_reg.date12_list = date12_list

    #sort regions by size to facilitate the region matching later
    common_regions.sort(key=lambda x: x.area, reverse=True)

    # plot sample result
    fig_size = pp.auto_figure_size(label_img.shape, disp_cbar=False)
    fig, ax = plt.subplots(figsize=fig_size)
    ax.imshow(label_img, cmap='jet')
    for common_reg in common_regions:
        ax.plot(common_reg.sample_coords[:,1],
                common_reg.sample_coords[:,0], 'k.', ms=2)
    pp.auto_flip_direction(stack_obj.metadata, ax, print_msg=False)
    out_img = 'common_region_sample.png'
    fig.savefig(out_img, bbox_inches='tight', transparent=True, dpi=300)
    print('saved common regions and sample pixels to file', out_img)
    plt.close(fig)

    return common_regions


def run_unwrap_error_phase_closure(ifgram_file, common_regions, water_mask_file=None, ccName='connectComponent',
                                   dsNameIn='unwrapPhase', dsNameOut='unwrapPhase_phaseClosure'):
    print('-'*50)
    print('correct unwrapping error in {} with phase closure ...'.format(ifgram_file))
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open()
    length, width = stack_obj.length, stack_obj.width
    ref_y, ref_x = stack_obj.refY, stack_obj.refX
    date12_list = stack_obj.get_date12_list(dropIfgram=False)
    num_ifgram = len(date12_list)
    shape_out = (num_ifgram, length, width)

    # read water mask
    if water_mask_file and os.path.isfile(water_mask_file):
        print('read water mask from file:', water_mask_file)
        water_mask = readfile.read(water_mask_file)[0]
    else:
        water_mask = None

    # prepare output data writing
    print('open {} with r+ mode'.format(ifgram_file))
    with h5py.File(ifgram_file, 'r+') as f:
        print('input  dataset:', dsNameIn)
        print('output dataset:', dsNameOut)
        if dsNameOut in f.keys():
            ds = f[dsNameOut]
            print('access /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))
        else:
            ds = f.create_dataset(dsNameOut,
                                  shape_out,
                                  maxshape=(None, None, None),
                                  chunks=True,
                                  compression=None)
            print('create /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))

        # correct unwrap error ifgram by ifgram
        prog_bar = ptime.progressBar(maxValue=num_ifgram)
        for i in range(num_ifgram):
            date12 = date12_list[i]

            # read unwrap phase to be updated
            unw_cor = np.squeeze(f[dsNameIn][i, :, :]).astype(np.float32)
            unw_cor -= unw_cor[ref_y, ref_x]

            # update kept interferograms only
            if stack_obj.dropIfgram[i]:
                # get local region info from connectComponent
                cc = np.squeeze(f[ccName][i, :, :])
                if water_mask is not None:
                    cc[water_mask == 0] = 0
                cc_obj = conncomp.connectComponent(conncomp=cc, metadata=stack_obj.metadata)
                cc_obj.label()
                local_regions = measure.regionprops(cc_obj.labelImg)

                # matching regions and correct unwrap error
                idx_common = common_regions[0].date12_list.index(date12)
                for local_reg in local_regions:
                    local_mask = cc_obj.labelImg == local_reg.label
                    U = 0
                    for common_reg in common_regions:
                        y = common_reg.sample_coords[:,0]
                        x = common_reg.sample_coords[:,1]
                        if all(local_mask[y, x]):
                            U = common_reg.int_ambiguity[idx_common]
                            break
                    unw_cor[local_mask] += 2. * np.pi * U

            # write to hdf5 file
            ds[i, :, :] = unw_cor
            prog_bar.update(i+1, suffix=date12)
        prog_bar.close()
        ds.attrs['MODIFICATION_TIME'] = str(time.time())
    print('close {} file.'.format(ifgram_file))

    return ifgram_file


####################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    if inps.action == 'correct':
        # update mode
        if inps.update_mode and run_or_skip(inps) == 'skip':
            return inps.ifgram_file

        # solve integer ambiguity for common connected components
        common_regions = get_common_region_int_ambiguity(ifgram_file=inps.ifgram_file,
                                                         cc_mask_file=inps.cc_mask_file,
                                                         water_mask_file=inps.waterMaskFile,
                                                         num_sample=inps.numSample,
                                                         dsNameIn=inps.datasetNameIn)

        # apply the integer ambiguity from common conn comp to the whole ifgram
        run_unwrap_error_phase_closure(inps.ifgram_file, common_regions,
                                       water_mask_file=inps.waterMaskFile,
                                       dsNameIn=inps.datasetNameIn,
                                       dsNameOut=inps.datasetNameOut)

    else:
        # calculate the number of triplets with non-zero integer ambiguity
        out_file = calc_num_triplet_with_nonzero_integer_ambiguity(inps.ifgram_file,
                                                                   mask_file=inps.waterMaskFile,
                                                                   dsName=inps.datasetNameIn,
                                                                   update_mode=inps.update_mode)
        # for debug
        #plot_num_triplet_with_nonzero_integer_ambiguity(out_file)
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
