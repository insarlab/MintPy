#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2018, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import time
from datetime import datetime as dt
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pysar.objects import ifgramStack
from pysar.utils import (ptime,
                         readfile,
                         writefile,
                         utils as ut,
                         plot as pp,
                         deramp)


key_prefix = 'pysar.unwrapError.'
# key configuration parameter name
config_key_list = ['maskFile',
                   'ramp',
                   'bridgeYX']


####################################################################################################
EXAMPLE = """Example:
  # Check the Jupyter Notebook for a complete tutorial of this procedure:
  https://github.com/yunjunz/pysar/blob/master/examples/run_unwrap_error_bridging.ipynb

  unwrap_error_bridging.py  ./INPUTS/ifgramStack.h5  -t GalapagosSenDT128.template --update
  unwrap_error_bridging.py  20180502_20180619.unw  -m maskConnComp.h5  -t GalapagosSenDT128.template
"""

REFERENCE = """Reference:
  Yunjun Z., H. Fattahi, F. Amelung (2018), PySAR - A Python package for Generic InSAR time-series
  analysis based on Full Rank Network. (in prep.)
"""

NOTE = """
  correct unwrapping errors by manually select bonding points to form bridges
  between patches which are internally correctly unwrapped

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
"""

TEMPLATE = """
## unwrapping error correction with bridging:
pysar.unwrapError.method   = auto   #[bridging / phase_closure / no], auto for no
pysar.unwrapError.maskFile = auto   #[file name / no], auto for no
pysar.unwrapError.ramp     = auto   #[plane / quadratic], auto for plane
pysar.unwrapError.bridgeYX = auto   #[y1_start,x1_start,y1_end,x1_end;y2_start,...], auto for none
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Unwrapping Error Correction with Bridging'+NOTE,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', type=str, help='interferograms file to be corrected')
    parser.add_argument('-t', '--template', dest='template_file', type=str, required=True,
                          help='template file with bonding point info, e.g.\n' +
                               'pysar.unwrapError.yx = 283,1177,305,1247;350,2100,390,2200')

    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help='name of dataset to be corrected, default: unwrapPhase')
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_bridge')

    parser.add_argument('-m','--mask', dest='maskFile', type=str,
                        help='name of mask file to mark different patches that want to be corrected in different value')
    parser.add_argument('--ramp', dest='ramp', choices=['plane', 'quadratic'],
                          help='type of phase ramp to be removed before correction.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_unwCor dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inps.bridgeYX = None
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if key in ['update', 'bridgeYX']:
            inpsDict[key] = value
        elif value:
            if key in ['maskFile', 'ramp']:
                inpsDict[key] = value

    if not inps.maskFile or not inps.bridgeYX:
        msg =  'No mask file or bridgeYX found.\n'
        msg += 'Bridging method is NOT automatic, you need to:\n'
        msg += '  1) prepare the connected components mask file to mark each area with the same unwrapping error\n'
        msg += '  2) select bridging points to connect areas one by one, starting from area where the reference pixel is.\n'
        msg += 'Check the following Jupyter Notebook for an example:\n'
        msg += '  https://github.com/yunjunz/pysar/blob/master/examples/run_unwrap_error_bridging.ipynb'
        raise SystemExit(msg)

    bridge_yx = inps.bridgeYX.replace(';', ' ').replace(',', ' ')  #convert ,/; into whitespace
    inps.bridgeYX = np.array([int(i) for i in bridge_yx.split()]).reshape(-1, 4)

    return inps


def plot_bridge_mask(mask_cc_file, bridge_yx, metadata=dict()):
    if not mask_cc_file:
        raise ValueError('no mask file for bridging found.')
    mask_cc = readfile.read(mask_cc_file)[0]

    # check number of bridges
    num_bridge = bridge_yx.shape[0]
    if np.unique(mask_cc).size < num_bridge + 1:
        raise ValueError('number of marked patches != number of briges.')

    # check the 1st bridge reference point is in the same patch with reference point
    if metadata:
        ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])
        if mask_cc[ref_y, ref_x] != mask_cc[bridge_yx[0, 0], bridge_yx[0, 1]]:
            raise ValueError(('1st bridge reference/start point is not in the same patch'
                              ' as file reference pixel.'))

    # save bridge points into text file
    out_file_base = 'maskConnCompBridge'
    txt_file = '{}.txt'.format(out_file_base)
    header_info = '# bridge bonding points for unwrap error correction\n'
    header_info += '# y0\tx0\ty1\tx1'
    np.savetxt(txt_file, bridge_yx, fmt='%s', delimiter='\t', header=header_info)
    print('save bridge points to file: {}'.format(txt_file))

    # plot/save mask and bridge setting into an image
    out_file = '{}.png'.format(out_file_base)
    fig, ax = plt.subplots(figsize=[6, 8])
    im = ax.imshow(mask_cc)
    ax = pp.auto_flip_direction(metadata, ax=ax, print_msg=False)
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax, ticks=np.arange(num_bridge+2))

    # plot bridges
    for i in range(num_bridge):
        y0, x0, y1, x1 = bridge_yx[i, :]
        ax.plot([x0, x1], [y0, y1], '-o', ms=5)
    plt.savefig(out_file, bbox_inches='tight', transparent=True, dpi=300)
    print('plot/save bridge setting to file: {}'.format(out_file))

    return


##########################################################################################
def bridge_unwrap_error(data, mask, bridge_yx):
    """Phase Jump Correction, using phase continuity on bridge/bonding points in each pair of patches.
    Inputs:
        data : 2D np.array in size of (length, width), phase matrix need to be corrected
        mask : 2D np.array in size of (length, width), mask of different patches
        bridge_yx : 2D np.array in size of (num_bridge, 4), start/reference point in 
    Output:
        data : 2D np.array, phase corrected matrix
    """
    num_bridge = bridge_yx.shape[0]
    for i in range(num_bridge):
        ref_yx = bridge_yx[i, 0:2]
        tgt_yx = bridge_yx[i, 2:4]

        # estimate integer number of phase jump
        ref_value = data[ref_yx[0], ref_yx[1]]
        tgt_value = data[tgt_yx[0], tgt_yx[1]]
        diff_value = tgt_value - ref_value
        num_jump = (abs(diff_value) + np.pi) // (2.*np.pi)
        if diff_value >= 0.:
            num_jump *= -1

        # correct unwrap error to target patch
        tgt_mask = mask == mask[tgt_yx[0], tgt_yx[1]]
        data[tgt_mask] += num_jump * 2. * np.pi

    return data


def run_unwrap_error_bridge(ifgram_file, mask_file, bridge_yx, dsNameIn='unwrapPhase',
                            dsNameOut='unwrapPhase_bridge', ramp_type=None):
    atr = readfile.read_attribute(ifgram_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    ref_y, ref_x = int(atr['REF_Y']), int(atr['REF_X'])
    k = atr['FILE_TYPE']

    print('read mask from file: {}'.format(mask_file))
    mask = readfile.read(mask_file, datasetName='mask')[0]

    if ramp_type is not None:
        print('estimate and remove phase ramp of {} during the correction'.format(ramp_type))
        ramp_mask = (mask == mask[ref_y, ref_x])

    print('correct unwrapping error in {} with bridging ...'.format(ifgram_file))
    if k == 'ifgramStack':
        date12_list = ifgramStack(ifgram_file).get_date12_list(dropIfgram=False)
        num_ifgram = len(date12_list)
        shape_out = (num_ifgram, length, width)

        # prepare output data writing
        print('open {} with r+ mode'.format(ifgram_file))
        f = h5py.File(ifgram_file, 'r+')
        if dsNameOut in f.keys():
            ds = f[dsNameOut]
            print('access /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))
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
            unw -= unw[ref_y, ref_x]

            # remove phase ramp before phase jump estimation
            if ramp_type is not None:
                unw, unw_ramp = deramp.remove_data_surface(unw, ramp_mask, ramp_type)
            else:
                unw_ramp = np.zeros(unw.shape, np.float32)

            # estimate/correct phase jump
            unw_cor = bridge_unwrap_error(unw, mask, bridge_yx)
            unw_cor += unw_ramp

            # write to hdf5 file
            ds[i, :, :] = unw_cor
            prog_bar.update(i+1, suffix=date12)
        prog_bar.close()

        f.close()
        print('close {} file.'.format(ifgram_file))

    if k == '.unw':
        # read data
        unw = readfile.read(ifgram_file)[0]
        unw -= unw[ref_y, ref_x]

        # remove phase ramp before phase jump estimation
        if ramp_type is not None:
            unw, unw_ramp = deramp.remove_data_surface(unw, ramp_mask, ramp_type)
        else:
            unw_ramp = np.zeros(unw.shape, np.float32)

        # estimate/correct phase jump
        unw_cor = bridge_unwrap_error(unw, mask, bridge_yx)
        unw_cor += unw_ramp

        # write to hdf5 file
        out_file = '{}_unwCor{}'.format(os.path.splitext(ifgram_file)[0],
                                        os.path.splitext(ifgram_file)[1])
        print('writing >>> {}'.format(out_file))
        writefile.write(unw_cor, out_file=out_file, ref_file=ifgram_file)

    return ifgram_file


####################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps = read_template2inps(inps.template_file, inps)
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_bridge'.format(inps.datasetNameIn)

    atr = readfile.read_attribute(inps.ifgram_file)

    # update mode
    if inps.update_mode and atr['FILE_TYPE'] == 'ifgramStack':
        print('update mode: {}'.format(inps.update_mode))
        stack_obj = ifgramStack(inps.ifgram_file)
        stack_obj.open(print_msg=False)
        atr = stack_obj.metadata
        if (inps.datasetNameOut in stack_obj.datasetNames 
                and all([str(vars(inps)[key]) == atr.get(key_prefix+key, 'None')
                         for key in config_key_list])):
            print("1) output dataset: {} exists".format(inps.datasetNameOut))
            print("2) all key configuration parameter are the same: \n\t{}".format(config_key_list))
            print("thus, skip this step.")
            return inps.ifgram_file

    # run bridging
    start_time = time.time()
    plot_bridge_mask(inps.maskFile, inps.bridgeYX, atr)
    run_unwrap_error_bridge(inps.ifgram_file,
                            mask_file=inps.maskFile,
                            bridge_yx=inps.bridgeYX,
                            dsNameIn=inps.datasetNameIn,
                            dsNameOut=inps.datasetNameOut,
                            ramp_type=inps.ramp)

    # config parameter
    print('add/update the following configuration metadata to file:')
    config_metadata = dict()
    for key in config_key_list:
        config_metadata[key_prefix+key] = str(vars(inps)[key])
    ut.add_attribute(inps.ifgram_file, config_metadata, print_msg=True)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
