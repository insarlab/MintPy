#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2018, Zhang Yunjun                     #
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
                         plot as pp)


# key configuration parameter name
key_prefix = 'pysar.unwrapError.'
configKeys = ['maskFile',
              'waterMaskFile',
              'ramp',
              'bridgePtsRadius']


####################################################################################################
EXAMPLE = """Example:
  # Check the Jupyter Notebook for a complete tutorial of this procedure:
  https://github.com/yunjunz/pysar/blob/master/examples/run_unwrap_error_bridging.ipynb

  unwrap_error_bridging.py  ./INPUTS/ifgramStack.h5  -t GalapagosSenDT128.template --update
  unwrap_error_bridging.py  ./INPUTS/ifgramStack.h5  --water-mask waterMask.h5
  unwrap_error_bridging.py  20180502_20180619.unw    -m maskConnComp.h5
"""

REFERENCE = """Reference:
  Yunjun Z., H. Fattahi, F. Amelung (2018), PySAR - A Python package for Generic InSAR time-series
  analysis based on Full Rank Network. (in prep.)
"""

NOTE = """
  correct unwrapping errors by building bridges between areas
  that are internally correctly unwrapped based on spatial continuity

  This method assumes:
  a. no phase unwrapping error within each patch marked by mask file.
  b. the absolute phase difference of bonding points (usually close in space) is 
     smaller than one pi. Considering prevalent ramps in InSAR data might break
     this assumptioin, optionally, a ramp can be estimated and removed during the
     phase jump estimation.
"""

TEMPLATE = """
## unwrapping error correction with bridging:
## automatic for islands with waterMaskFile option (unwrapping errors on areas separated by narrow water bodies)
## manual    for all the other scenarios
pysar.unwrapError.method          = auto  #[bridging / phase_closure / no], auto for no
pysar.unwrapError.waterMaskFile   = auto  #[waterMask.h5 / no], auto for no
pysar.unwrapError.maskFile        = auto  #[maskConnComp.h5 / no], auto for no, mask for connected components areas
pysar.unwrapError.bridgePtsRadius = auto  #[1-inf], auto for 150, radius in pixel of circular area around bridge ends
pysar.unwrapError.ramp            = auto  #[linear / quadratic], auto for linear
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Unwrapping Error Correction with Bridging'+NOTE,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', type=str, help='interferograms file to be corrected')
    parser.add_argument('-m','--mask', dest='maskFile', type=str,
                        help='name of mask file to mark different patches that want to be corrected in different value')
    parser.add_argument('-r','--radius', dest='bridgePtsRadius', type=int, default=150,
                        help='radius of the end point of bridge to search area to get median representative value')
    parser.add_argument('--ramp', dest='ramp', choices=['linear', 'quadratic'],
                          help='type of phase ramp to be removed before correction.')
    parser.add_argument('--water-mask', dest='waterMaskFile', type=str,
                        help='path of water mask file for study area with unwrapping errors due to water separation.')
    parser.add_argument('--coh-mask', dest='cohMaskFile', type=str, default='maskSpatialCoh.h5',
                        help='path of mask file from average spatial coherence' + 
                             ' to mask out low coherent pixels for bridge point circular area selection')

    parser.add_argument('-t', '--template', dest='template_file', type=str,
                          help='template file with bonding point info, e.g.\n' +
                               'pysar.unwrapError.yx = 283,1177,305,1247;350,2100,390,2200')

    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help='name of dataset to be corrected, default: unwrapPhase')
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_bridge')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_unwCor dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.ifgram_file)
    if atr['FILE_TYPE'] != 'ifgramStack':
        raise ValueError('input file is not ifgramStack: {}'.format(atr['FILE_TYPE']))

    # default output dataset name
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_bridge'.format(inps.datasetNameIn)
    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if key in ['update']:
            inpsDict[key] = value
        elif value:
            if key in ['maskFile', 'waterMaskFile', 'ramp']:
                inpsDict[key] = value
            elif key in ['bridgePtsRadius']:
                inpsDict[key] = int(value)
    return inps


def run_check(inps):
    print('-'*50)
    print('update mode: ON')
    run = False

    # check output dataset
    with h5py.File(inps.ifgram_file, 'r') as f:
        if inps.datasetNameOut not in f.keys():
            run = True
            print('  1) output dataset: {} not found --> run.'.format(inps.datasetNameOut))
        else:
            print('  1) output dataset: {} exists'.format(inps.datasetNameOut))
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if to <= ti:
                run = True
                print('  2) output dataset is NOT newer than input dataset: {} --> run.'.format(inps.datasetNameIn))
            else:
                print('  2) output dataset is newer than input dataset: {}'.format(inps.datasetNameIn))

    # check configuration
    if not run:
        # convert inps value to common str format
        inps_dict = dict(vars(inps))
        specialValues = {True : 'yes',
                         False: 'no',
                         None : 'no'}
        for key, value in inps_dict.items():
            if value in specialValues.keys():
                inps_dict[key] = specialValues[value]
        atr = readfile.read_attribute(inps.ifgram_file)
        if any(str(inps_dict[key]) != atr.get(key_prefix+key, 'no') for key in configKeys):
            run = True
            print('  3) NOT all key configration parameters are the same --> run.\n\t{}'.format(configKeys))
        else:
            print('  3) all key configuration parameters are the same:\n\t{}'.format(configKeys))

    # result
    if run:
        print('run.')
    else:
        print('skip the run.')
    return run


####################################################################################################
def water_mask2conn_comp_mask(water_mask_file, ref_yx, out_file='maskConnComp.h5',
                              min_num_pixel=5e4, display=False):
    """Generate connected component mask file from water mask file
    Parameters: water_mask_file : str, path of water mask file
                ref_yx          : tuple of 2 int, row/col number of reference point
                out_file        : str, filename of output connected components mask
                min_num_pixel   : float, min number of pixels to be identified as conn comp
                display         : bool, display generated conn comp mask
    Returns:    out_file        : str, filename of output connected components mask
    """
    print('-'*50)
    print('generate connected component mask from water mask file: ', water_mask_file)
    water_mask, atr = readfile.read(water_mask_file)

    mask_cc = np.zeros(water_mask.shape, dtype=np.int16)
    # first conn comp - reference conn comp
    num_cc = 1
    label_mask = ndimage.label(water_mask)[0]
    mask_cc += label_mask == label_mask[ref_yx[0], ref_yx[1]]

    # all the other conn comps
    water_mask ^= mask_cc == mask_cc[ref_yx[0], ref_yx[1]]
    mask_ccs = ut.get_all_conn_components(water_mask, min_num_pixel=min_num_pixel)
    if mask_ccs:
        for mask_cci in mask_ccs:
            num_cc += 1
            mask_cc += mask_cci * num_cc

    # write file
    atr['FILE_TYPE'] = 'mask'
    atr['REF_Y'] = str(ref_yx[0])
    atr['REF_X'] = str(ref_yx[1])
    writefile.write(mask_cc, out_file=out_file, metadata=atr)

    # plot
    out_img = '{}.png'.format(os.path.splitext(out_file)[0])
    fig, ax = plt.subplots(figsize=[6, 8])
    im = ax.imshow(mask_cc)
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax, ticks=np.arange(num_cc+1))
    ax = pp.auto_flip_direction(atr, ax=ax, print_msg=False)
    fig.savefig(out_img, bbox_inches='tight', transparent=True, dpi=300)
    print('save figure to {}'.format(out_img))
    if display:
        plt.show()
    return out_file


def search_bridge(mask_cc_file, radius=150, coh_mask_file='maskSpatialCoh.h5', save_plot=True):
    """Search bridges to connect coherent conn comps with min distance
    Parameters: mask_cc_file : str, path of mask file of coherent conn comps
                radius : int, radius of influence for median value calculation to represent the value of cc
                coh_mask_file : str, path of mask file based on average spatial coherence
                    To mask out low coherent pixel during the circular mask selection
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
    if num_bridge < 1:
        raise RuntimeError('No bridge found. Check the input mask file.')
    print('number of bridges: {}'.format(num_bridge))

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
    print('2. find MST bridges and determine the briding order using breadth-first ordering')
    dist_mat_mst = sparse.csgraph.minimum_spanning_tree(dist_mat)
    nodes, predecessors = sparse.csgraph.breadth_first_order(dist_mat_mst,
                                                             i_start=0,
                                                             directed=False)
    nodes += 1
    predecessors += 1

    # converting bridging order into bridges
    print('3. find circular area around bridge point with radius of {} pixels'.format(radius))
    if os.path.isfile(coh_mask_file):
        print('  with background mask from: '+coh_mask_file)
        coh_mask = readfile.read(coh_mask_file)[0]
    else:
        coh_mask = np.ones(mask_cc.shape, dtype=np.bool_)

    bridges = []
    for i in range(num_bridge):
        # get x0/y0/x1/y1
        n0, n1 = predecessors[i+1], nodes[i+1]
        conn_dict = connections['{}{}'.format(n0, n1)]
        x0, y0 = conn_dict[str(n0)]
        x1, y1 = conn_dict[str(n1)]

        # get mask0/mask1
        mask0 = (mask_cc == n0) * ut.get_circular_mask(x0, y0, radius, mask_cc.shape) * coh_mask
        mask1 = (mask_cc == n1) * ut.get_circular_mask(x1, y1, radius, mask_cc.shape) * coh_mask

        # save to list of dict
        bridge = dict()
        bridge['x0'] = x0
        bridge['y0'] = y0
        bridge['x1'] = x1
        bridge['y1'] = y1
        bridge['mask0'] = mask0
        bridge['mask1'] = mask1
        bridges.append(bridge)

    out_base='maskConnCompBridge'
    # save bridge points to text file
    save_bridge_txt(bridges, out_file='{}.txt'.format(out_base))

    # plot bridge
    if save_plot:
        fig, ax = plt.subplots(figsize=[6, 8])
        ax, im = plot_bridge(ax, mask_cc_file, bridges)
        # colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        cbar = plt.colorbar(im, cax=cax, ticks=np.arange(num_bridge+2))
        out_file='{}.png'.format(out_base)
        plt.savefig(out_file, bbox_inches='tight', transparent=True, dpi=300)
        print('plot bridge setting to file: {}'.format(out_file))
    return bridges


def save_bridge_txt(bridges, out_file='maskConnCompBridge.txt'):
    """save to text file"""
    num_bridge = len(bridges)
    bridge_yx = np.zeros((num_bridge, 4), dtype=np.int16)
    for i in range(num_bridge):
        bridge = bridges[i]
        bridge_yx[i, :] = [bridge['y0'], bridge['x0'], bridge['y1'], bridge['x1']]
    header_info = 'bridge bonding points for unwrap error correction\n'
    header_info += 'y0\tx0\ty1\tx1'
    np.savetxt(out_file, bridge_yx, fmt='%s', delimiter='\t', header=header_info)
    print('save bridge points  to file: {}'.format(out_file))
    return out_file


def plot_bridge(ax, mask_cc_file, bridges):
    """Plot mask of connected components with bridges info
    Parameters: mask_cc_file : string, path of mask cc file
                bridges : list of dict
    """
    mask_cc, metadata = readfile.read(mask_cc_file)
    num_bridge = len(bridges)

    # plot 1. mask_cc data
    im = ax.imshow(mask_cc)

    # plot 2. bridge data
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
    return ax, im


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
    def most_common_median(array):
        """Get the median value of the most common range of input array"""
        bin_value, bin_edge = np.histogram(array)
        idx = np.argmax(bin_value)
        mask = np.multiply(array >= bin_edge[idx], array <= bin_edge[idx+1])
        d_comm = np.nanmedian(array[mask])
        return d_comm

    data = np.array(data, dtype=np.float32)
    mask = data != 0.
    mask *= np.isfinite(data)
    for i in range(len(bridges)):
        # get phase difference between two ends of the bridge
        bridge = bridges[i]
        value0 = np.nanmedian(data[bridge['mask0'] * mask])
        value1 = np.nanmedian(data[bridge['mask1'] * mask])
        diff_value = value1 - value0

        #estimate integer number of phase jump
        num_jump = (np.abs(diff_value) + np.pi) // (2.*np.pi)
        if diff_value > 0:
            num_jump *= -1

        # correct unwrap error to target conn comps
        mask_tgt = mask_cc == mask_cc[bridge['y1'], bridge['x1']]
        data[mask_tgt] += num_jump * 2. * np.pi
    return data


def run_unwrap_error_bridge(ifgram_file, mask_cc_file, bridges, dsNameIn='unwrapPhase',
                            dsNameOut='unwrapPhase_bridge', ramp_type=None):
    """Run unwrapping error correction with bridging
    Parameters: ifgram_file  : str, path of ifgram stack file
                mask_cc_file : str, path of conn comp mask file
                bridges      : list of dicts, check bridge_unwrap_error() for details
                dsNameIn     : str, dataset name of unwrap phase to be corrected
                dsNameOut    : str, dataset name of unwrap phase to be saved after correction
                ramp_type    : str, name of phase ramp to be removed during the phase jump estimation
    Returns:    ifgram_file  : str, path of ifgram stack file
    """
    print('-'*50)
    print('correct unwrapping error in {} with bridging ...'.format(ifgram_file))
    # file info
    atr = readfile.read_attribute(ifgram_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    ref_y, ref_x = int(atr['REF_Y']), int(atr['REF_X'])
    k = atr['FILE_TYPE']

    # read mask
    print('read mask from file: {}'.format(mask_cc_file))
    mask_cc = readfile.read(mask_cc_file, datasetName='mask')[0]
    if ramp_type is not None:
        print('estimate and remove phase ramp of {} during the correction'.format(ramp_type))
        mask4ramp = (mask_cc == mask_cc[ref_y, ref_x])

    # correct unwrap error ifgram by ifgram
    if k == 'ifgramStack':
        date12_list = ifgramStack(ifgram_file).get_date12_list(dropIfgram=False)
        num_ifgram = len(date12_list)
        shape_out = (num_ifgram, length, width)

        # prepare output data writing
        print('open {} with r+ mode'.format(ifgram_file))
        f = h5py.File(ifgram_file, 'r+')
        print('input  dataset:', dsNameIn)
        print('output dataset:', dsNameOut)
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
            unw[unw != 0.] -= unw[ref_y, ref_x]

            # remove phase ramp before phase jump estimation
            if ramp_type is not None:
                unw, unw_ramp = ut.deramp_data(unw, mask4ramp, ramp_type, metadata=atr)

            # estimate/correct phase jump
            unw_cor = bridge_unwrap_error(unw, mask_cc, bridges)
            if ramp_type is not None:
                unw_cor += unw_ramp

            # write to hdf5 file
            ds[i, :, :] = unw_cor
            prog_bar.update(i+1, suffix=date12)
        prog_bar.close()
        ds.attrs['MODIFICATION_TIME'] = str(time.time())
        f.close()
        print('close {} file.'.format(ifgram_file))

    if k == '.unw':
        # read data
        unw = readfile.read(ifgram_file)[0]
        unw[unw != 0.] -= unw[ref_y, ref_x]

        # remove phase ramp before phase jump estimation
        if ramp_type is not None:
            unw, unw_ramp = ut.deramp_data(unw, mask4ramp, ramp_type, metadata=atr)

        # estimate/correct phase jump
        unw_cor = bridge_unwrap_error(unw, mask_cc, bridges)
        if ramp_type is not None:
            unw_cor += unw_ramp

        # write to hdf5 file
        out_file = '{}_unwCor{}'.format(os.path.splitext(ifgram_file)[0],
                                        os.path.splitext(ifgram_file)[1])
        print('writing >>> {}'.format(out_file))
        writefile.write(unw_cor, out_file=out_file, ref_file=ifgram_file)

    return ifgram_file


####################################################################################################
def main(iargs=None):
    # check inputs
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # update mode
    if inps.update_mode and run_check(inps) is False:
        return inps.ifgram_file

    # check maskConnComp.h5
    if os.path.isfile(str(inps.waterMaskFile)) and not inps.maskFile:
        atr = readfile.read_attribute(inps.ifgram_file)
        ref_yx = (int(atr['REF_Y']), int(atr['REF_X']))
        inps.maskFile = water_mask2conn_comp_mask(inps.waterMaskFile,
                                                  ref_yx=ref_yx,
                                                  min_num_pixel=1e4)
    if not inps.maskFile:
        msg =  'No mask of connected components file found. Bridging method is NOT automatic, you need to:\n'
        msg += '  Prepare the connected components mask file to mark each area with the same unwrapping error\n'
        msg += 'Check the following Jupyter Notebook for an example:\n'
        msg += '  https://github.com/yunjunz/pysar/blob/master/examples/run_unwrap_error_bridging.ipynb'
        raise SystemExit(msg)

    # run bridging
    start_time = time.time()
    bridges = search_bridge(inps.maskFile,
                            radius=inps.bridgePtsRadius,
                            coh_mask_file=inps.cohMaskFile)
    run_unwrap_error_bridge(inps.ifgram_file,
                            inps.maskFile,
                            bridges,
                            dsNameIn=inps.datasetNameIn,
                            dsNameOut=inps.datasetNameOut,
                            ramp_type=inps.ramp)

    # config parameter
    print('add/update the following configuration metadata to file:')
    config_metadata = dict()
    for key in configKeys:
        config_metadata[key_prefix+key] = str(vars(inps)[key])
    ut.add_attribute(inps.ifgram_file, config_metadata, print_msg=True)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
