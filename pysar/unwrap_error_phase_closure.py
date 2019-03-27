#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2019, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import argparse
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt
try:
    from cvxopt import matrix
except ImportError:
    raise ImportError('Cannot import cvxopt')
try:
    from skimage import measure
except ImportError:
    raise ImportError('Could not import skimage!')

from pysar.objects import ifgramStack
from pysar.objects.conncomp import connectComponent
from pysar.utils import ptime, readfile, utils as ut, plot as pp
from pysar.utils.solvers import l1regls


key_prefix = 'pysar.unwrapError.'

##########################################################################################
EXAMPLE = """Example:
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  maskConnComp.h5  -t pysarApp_template.txt  --update
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  maskConnComp.h5  --water-mask waterMask.h5 --update
"""

TEMPLATE = """
## Unwrapping Error Correction based on Phase Closure (Yunjun et al., 2019)
pysar.unwrapError.waterMaskFile   = auto  #[waterMask.h5 / no], auto for no
"""

REFERENCE = """Reference:
  Yunjun, Z., H. Fattahi, F. Amelung, (2019) InSAR time series analysis: error correction
  and noise reduction (submitted).
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


def create_parser():
    parser = argparse.ArgumentParser(description='Unwrapping Error Correction based on Phase Closure'+NOTE,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be corrected')
    parser.add_argument('cc_mask_file', default='maskConnComp.h5', help='common connected components file')
    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help="name of dataset to be corrected, default: unwrapPhase")
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_phaseClosure')

    parser.add_argument('--water-mask','--wm', dest='waterMaskFile', type=str, help='path of water mask file.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options for setting.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_phaseClosure dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check 1 input file type
    k = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if k not in ['ifgramStack']:
        raise ValueError('input file is not ifgramStack: {}'.format(k))
    # check 2 cc_mask_file
    if not os.path.isfile(inps.cc_mask_file):
        raise FileNotFoundError(inps.cc_mask_file)
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
def write_hdf5_file_patch(ifgram_file, data, box=None, dsName='unwrapPhase_phaseClosure'):
    """Write a patch of 3D dataset into an existing h5 file.
    Parameters: ifgram_file : string, name/path of output hdf5 file
                data : 3D np.array to be written
                box  : tuple of 4 int, indicating of (x0, y0, x1, y1) of data in file
                dsName : output dataset name
    Returns:    ifgram_file
    """
    num_ifgram, length, width = ifgramStack(ifgram_file).get_size(dropIfgram=False)
    if not box:
        box = (0, 0, width, length)
    num_row = box[3] - box[1]
    num_col = box[2] - box[0]

    # write to existing HDF5 file
    print('open {} with r+ mode'.format(ifgram_file))
    f = h5py.File(ifgram_file, 'r+')

    # get h5py.Dataset
    msg = 'dataset /{d} of {t:<10} in size of {s}'.format(d=dsName, t=str(data.dtype),
                                                          s=(num_ifgram, box[3], box[2]))
    if dsName in f.keys():
        print('update '+msg)
        ds = f[dsName]
    else:
        print('create '+msg)
        ds = f.create_dataset(dsName, (num_ifgram, num_row, num_col),
                              maxshape=(None, None, None),
                              chunks=True, compression=None)

    # resize h5py.Dataset if current size is not enough
    if ds.shape != (num_ifgram, length, width):
        ds.resize((num_ifgram, box[3], box[2]))

    # write data to file
    ds[:, box[1]:box[3], box[0]:box[2]] = data

    ds.attrs['MODIFICATION_TIME'] = str(time.time())
    f.close()
    print('close {}'.format(ifgram_file))
    return ifgram_file


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

    label_img, num_label = connectComponent.get_large_label(cc_mask, min_area=2.5e3, print_msg=True)
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
                box = (x, y, x+1, y+1)
                unw = readfile.read(ifgram_file, datasetName=dsNameIn, box=box)[0].reshape(num_ifgram, -1)
                unw -= ref_phase

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
    f = h5py.File(ifgram_file, 'r+')
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
        # skip excluded interferograms
        if not stack_obj.dropIfgram[i]:
            next

        date12 = date12_list[i]
        idx_common = common_regions[0].date12_list.index(date12)
        # read unwrap phase to be updated
        unw_cor = np.squeeze(f[dsNameIn][i, :, :]).astype(np.float32)
        unw_cor -= unw_cor[ref_y, ref_x]

        # get local region info from connectComponent
        cc = np.squeeze(f[ccName][i, :, :])
        if water_mask is not None:
            cc[water_mask == 0] = 0
        cc_obj = connectComponent(conncomp=cc, metadata=stack_obj.metadata)
        cc_obj.label()
        local_regions = measure.regionprops(cc_obj.labelImg)

        # matching regions and correct unwrap error
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
    f.close()
    print('close {} file.'.format(ifgram_file))
    return ifgram_file


####################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_phaseClosure'.format(inps.datasetNameIn)

    # update mode
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.ifgram_file

    start_time = time.time()
    common_regions = get_common_region_int_ambiguity(ifgram_file=inps.ifgram_file,
                                                     cc_mask_file=inps.cc_mask_file,
                                                     water_mask_file=inps.waterMaskFile,
                                                     num_sample=100,
                                                     dsNameIn=inps.datasetNameIn)

    run_unwrap_error_phase_closure(inps.ifgram_file, common_regions,
                                   water_mask_file=inps.waterMaskFile,
                                   dsNameIn=inps.datasetNameIn,
                                   dsNameOut=inps.datasetNameOut)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
