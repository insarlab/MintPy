############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import time

import h5py
import numpy as np
from cvxopt import matrix
from matplotlib import pyplot as plt
from skimage import measure

from mintpy import ifgram_inversion as ifginv
from mintpy.objects import conncomp, ifgramStack
from mintpy.utils import plot as pp, ptime, readfile, utils as ut, writefile
from mintpy.utils.solvers import l1regls

key_prefix = 'mintpy.unwrapError.'


##########################################################################################
def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output dataset
    with h5py.File(inps.ifgram_file, 'r') as f:
        if inps.datasetNameOut not in f.keys():
            flag = 'run'
            print(f'1) output dataset: {inps.datasetNameOut} NOT found.')
        else:
            print(f'1) output dataset: {inps.datasetNameOut} exists.')
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if ti > to:
                flag = 'run'
                print(f'2) output dataset is NOT newer than input dataset: {inps.datasetNameIn}.')
            else:
                print(f'2) output dataset is newer than input dataset: {inps.datasetNameIn}.')

    # result
    print(f'run or skip: {flag}.')
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
        out_file = 'numTriNonzeroIntAmbiguity'
        out_file += f'4{dsName}' if dsName != 'unwrapPhase' else ''
        out_file += '.h5'
        out_file = os.path.join(out_dir, out_file)

    # update mode
    if update_mode and os.path.isfile(out_file):
        print('update mode: ON')
        print(f'1) output file "{out_file}" already exists')
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
            print(f'3) NOT all key configurations are the same: {key_list}')
            flag = 'run'
        else:
            print(f'3) all key configurations are the same: {key_list}')

        print(f'run or skip: {flag}.')
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
        msg = f'No triangles found from ifgramStack file: {ifgram_file}!'
        msg += '\n    Skip calculating the number of triplets with non-zero integer ambiguity.'
        print(msg)
        return None
    else:
        print(f'number of interferograms: {C.shape[1]}')
        print(f'number of triplets: {C.shape[0]}')

    # calculate number of nonzero closure phase
    ds_size = (C.shape[0] * 2 + C.shape[1]) * length * width * 4
    num_loop = int(np.ceil(ds_size * 2 / (max_memory * 1024**3)))
    # ensure a min step size of 10
    step = int(np.ceil(length / num_loop / 10)) * 10
    num_loop = int(np.ceil(length / step))
    num_nonzero_closure = np.zeros((length, width), dtype=np.float32)
    msg = 'calculating the number of triplets with non-zero integer ambiguity of closure phase ...'
    msg += f'\n    block by block with size up to {(step, width)}, {num_loop} blocks in total'
    print(msg)

    ref_phase = stack_obj.get_reference_phase(
        unwDatasetName=dsName,
        dropIfgram=True,
    ).reshape(num_ifgram, -1)

    prog_bar = ptime.progressBar(maxValue=num_loop)
    for i in range(num_loop):
        # box
        r0 = i * step
        r1 = min((r0 + step), stack_obj.length)
        box = (0, r0, stack_obj.width, r1)

        # read data
        unw = ifginv.read_stack_obs(
            stack_obj,
            box=box,
            ref_phase=ref_phase,
            obs_ds_name=dsName,
            dropIfgram=True,
            print_msg=False,
        ).reshape(num_ifgram, -1)

        # calculate based on equation (8-9) and T_int equation inline.
        closure_pha = np.dot(C, unw)
        closure_int = np.round((closure_pha - ut.wrap(closure_pha)) / (2.*np.pi))
        num_nonzero_closure[r0:r1, :] = np.sum(closure_int != 0, axis=0).reshape(-1, width)

        prog_bar.update(i+1, every=1, suffix=f'line {r0} / {length}')
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
    writefile.write(
        num_nonzero_closure,
        out_file=out_file,
        metadata=meta,
    )

    # plot
    pp.plot_num_triplet_with_nonzero_integer_ambiguity(out_file)

    return out_file


##########################################################################################
def get_common_region_int_ambiguity(ifgram_file, cc_mask_file, water_mask_file=None, num_sample=100,
                                    dsNameIn='unwrapPhase', cc_min_area=2.5e3):
    """Solve the phase unwrapping integer ambiguity for the common regions among all interferograms
    Parameters: ifgram_file     - str, path of interferogram stack file
                cc_mask_file    - str, path of common connected components file
                water_mask_file - str, path of water mask file
                num_sample      - int, number of pixel sampled for each region
                dsNameIn        - str, dataset name of the unwrap phase to be corrected
                cc_min_area     - float, minimum region/area size
    Returns:    common_regions  - list of skimage.measure._regionprops._RegionProperties object
                                  modified by adding two more variables:
                                  sample_coords : 2D np.ndarray in size of (num_sample, 2) in int64 format
                                  int_ambiguity : 1D np.ndarray in size of (num_ifgram,) in int format
    """
    print('-'*50)
    print('calculating the integer ambiguity for the common regions defined in', cc_mask_file)
    # stack info
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=True)
    num_ifgram = len(date12_list)

    C = matrix(ifgramStack.get_design_matrix4triplet(date12_list).astype(float))
    ref_phase = stack_obj.get_reference_phase(
        unwDatasetName=dsNameIn,
        dropIfgram=True,
    ).reshape(num_ifgram, -1)
    print(f'number of interferograms: {num_ifgram}')
    print(f'number of triplets: {int(len(C)/num_ifgram)}')

    # prepare common label
    print('read common mask from', cc_mask_file)
    cc_mask = readfile.read(cc_mask_file)[0]
    if water_mask_file is not None and os.path.isfile(water_mask_file):
        water_mask = readfile.read(water_mask_file)[0]
        print('refine common mask based on water mask file', water_mask_file)
        cc_mask[water_mask == 0] = 0

    label_img, num_label = conncomp.label_conn_comp(cc_mask, min_area=cc_min_area, print_msg=True)
    common_regions = measure.regionprops(label_img)
    print('number of common regions:', num_label)
    if num_label == 0:
        msg = 'WARNING: NO common region found! '
        msg += f'Try a smaller value for the mintpy.unwrapError.connCompMinArea ({cc_min_area}) option.'
        print(msg)

    else:
        # add sample_coords / int_ambiguity
        print('number of samples per region:', num_sample)
        print(f'solving the phase-unwrapping integer ambiguity for {dsNameIn}')
        print('\tbased on the closure phase of interferograms triplets (Yunjun et al., 2019)')
        print('\tusing the L1-norm regularzed least squares approximation (LASSO) ...')
        for i in range(num_label):
            common_reg = common_regions[i]
            # sample_coords
            rng = np.random.default_rng()
            idx = sorted(rng.choice(int(common_reg.area), num_sample, replace=False))
            common_reg.sample_coords = common_reg.coords[idx, :].astype(int)

            # solve for int_ambiguity
            U = np.zeros((num_ifgram, num_sample))
            if common_reg.label == label_img[stack_obj.refY, stack_obj.refX]:
                print(f'{i+1}/{num_label} skip calculation for the reference region')

            else:
                prog_bar = ptime.progressBar(maxValue=num_sample, prefix=f'{i+1}/{num_label}')
                for j in range(num_sample):
                    # read unwrap phase
                    y, x = common_reg.sample_coords[j, :]
                    unw = ifginv.read_stack_obs(
                        stack_obj,
                        box=(x, y, x+1, y+1),
                        ref_phase=ref_phase,
                        obs_ds_name=dsNameIn,
                        dropIfgram=True,
                        print_msg=False,
                    ).reshape(num_ifgram, -1)

                    # calculate closure_int
                    closure_pha = np.dot(C, unw)
                    closure_int = matrix(np.round((closure_pha - ut.wrap(closure_pha)) / (2.*np.pi)))

                    # solve for U
                    U[:,j] = np.round(l1regls(
                        A=-C,
                        y=closure_int,
                        alpha=1e-2,
                        show_progress=0,
                    )).flatten()

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


def correct_unwrap_error_phase_closure(ifgram_file, common_regions, water_mask_file=None,
                                       ccName='connectComponent', dsNameIn='unwrapPhase',
                                       dsNameOut='unwrapPhase_phaseClosure'):
    print('-'*50)
    print(f'correct unwrapping error in {ifgram_file} with phase closure ...')
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
    print(f'open {ifgram_file} with r+ mode')
    with h5py.File(ifgram_file, 'r+') as f:
        print('input  dataset:', dsNameIn)
        print('output dataset:', dsNameOut)
        if dsNameOut in f.keys():
            ds = f[dsNameOut]
            print(f'access /{dsNameOut} of np.float32 in size of {shape_out}')
        else:
            ds = f.create_dataset(
                dsNameOut,
                shape_out,
                maxshape=(None, None, None),
                chunks=True,
                compression=None)
            print(f'create /{dsNameOut} of np.float32 in size of {shape_out}')

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
    print(f'close {ifgram_file} file.')

    return ifgram_file


##########################################################################################
def run_unwrap_error_phase_closure(inps):

    # matplotlib backend setting
    plt.switch_backend('Agg')
    start_time = time.time()

    if inps.action == 'correct':
        # action 1: correct for unwrapping errors

        # update mode
        if inps.update_mode and run_or_skip(inps) == 'skip':
            return

        # solve integer ambiguity for common connected components
        common_regions = get_common_region_int_ambiguity(
            ifgram_file=inps.ifgram_file,
            cc_mask_file=inps.cc_mask_file,
            water_mask_file=inps.waterMaskFile,
            num_sample=inps.numSample,
            dsNameIn=inps.datasetNameIn,
            cc_min_area=inps.connCompMinArea,
        )

        # apply the integer ambiguity from common conn comp to the whole ifgram
        if len(common_regions) == 0:
            print('skip phase closure correction ...')
            return

        correct_unwrap_error_phase_closure(
            inps.ifgram_file, common_regions,
            water_mask_file=inps.waterMaskFile,
            dsNameIn=inps.datasetNameIn,
            dsNameOut=inps.datasetNameOut,
        )

    else:
        # action 2: quick overview
        # calculate the number of triplets with non-zero integer ambiguity
        out_file = calc_num_triplet_with_nonzero_integer_ambiguity(
            inps.ifgram_file,
            mask_file=inps.waterMaskFile,
            dsName=inps.datasetNameIn,
            update_mode=inps.update_mode,
        )

        # for debug
        debug_mode = False
        if debug_mode:
            pp.plot_num_triplet_with_nonzero_integer_ambiguity(out_file)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs\nDone.')

    return
