#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import argparse
import time
import h5py
import numpy as np
from scipy import linalg
from pysar.objects import ifgramStack
from pysar.utils import ptime, readfile, writefile, utils as ut
import pysar.ifgram_inversion as ifginv


key_prefix = 'pysar.unwrapError.'

##########################################################################################
EXAMPLE = """Example:
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  -m mask.h5
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  -m waterMask.h5 --update
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  -m waterMask.h5 --fast
"""

TEMPLATE = """Template:
  ## Unwrapping Error Correction based on Phase Closure (Fattahi, 2015, PhD Thesis)
  pysar.unwrapError.maskFile = auto  #[file name / no], auto for no
  pysar.unwrapError.fastMode = auto  #[yes / no], auto for yes, enable fast mode for phase_closure method
"""

REFERENCE = """Reference:
  Fattahi, H. (2015), Geodetic Imaging of Tectonic Deformation with InSAR, 190 pp, 
  University of Miami, Miami, FL. (Chap. 4), Open Access Dissertations. 1456. 
  https://scholarlyrepository.miami.edu/oa_dissertations/1456
"""

NOTE = """
  correct unwrapping errors based on phase closure of pairs circle (ab + bc + ca == 0).
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
    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help="name of dataset to be corrected, default: unwrapPhase")
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_phaseClosure')

    parser.add_argument('-m','--mask', dest='maskFile',
                        help='mask file to specify those pixels to be corrected for unwrapping errors')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options for setting.')
    parser.add_argument('--fast', dest='fastMode', action='store_true',
                        help='Fast (but not the best) unwrap error correction,'
                             ' by diable the extra constraint on ifgrams with no unwrap error.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_phaseClosure dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
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
        if key in ['fastMode']:
            inpsDict[key] = value
        elif value:
            if key in ['maskFile']:
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
            print('  1) output dataset: {} not found --> run.'.format(inps.datasetNameOut))
        else:
            print('  1) output dataset: {} exists'.format(inps.datasetNameOut))
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if ti > to:
                flag = 'run'
                print('  2) output dataset is NOT newer than input dataset: {} --> run.'.format(inps.datasetNameIn))
            else:
                print('  2) output dataset is newer than input dataset: {}'.format(inps.datasetNameIn))

    # result
    print('check result:', flag)
    return flag


####################################################################################################
def correct_unwrap_error(ifgram, C, Dconstraint=True, thres=0.1, alpha=0.25, rcond=1e-5):
    """Estimate unwrapping error from a stack of unwrapped interferometric phase
    Parameters: ifgram : 2D np.array in size of (num_ifgram, num_pixel) of unwrap phase
                C      : 2D np.array in size of (num_triangle, num_ifgram) triangle design matrix
                Dconstraint : bool, apply zero phase jump constraint on ifgrams without unwrapping error.
                    This is disabled in fast option
                thres  : float, threshold of non-zero phase closure beyond which indicates unwrap error
                alpha  : float, Tikhonov factor, regularization parameter
                rcond  : float, cut-off value for least square estimation
    Returns:    ifgram_cor : 2D np.array in size of (num_ifgram, num_pixel) of unwrap phase after correction
                U          : 2D np.array in size of (num_ifgram, num_pixel) of phase jump integer
    Example:    ifgram_cor, U = estimate_unwrap_error(ifgram, C)
    """
    num_tri, num_ifgram = C.shape
    ifgram = ifgram.reshape(num_ifgram, -1)
    U = np.zeros(ifgram.shape, np.float32)

    I = np.eye(num_ifgram, dtype=np.float32)
    A = np.vstack((C, alpha*I))
    # Use Tikhonov regularization (alpha*D0) to solve equation in ill-posed scenario:
    # e.g.: network is not dense enough AND a lot of unwrapping error exist,
    # then A.shape[0] < A.shape[1]: number of unknown > number of observations

    if Dconstraint:
        # get D
        # 1. calculate phase closure for all ifgram triangles --> identy those with unwrap error
        # 2. identy ifgrams involved in at least 1 triangle with unwrap error
        # 3. get D for ifgram 'without any unwrap error' or 'correctly unwrapped'
        pha_closure = np.dot(C, ifgram)
        pha_closure = np.abs(pha_closure - ut.wrap(pha_closure))
        idx_nonzero_closure = (pha_closure >= thres).flatten()
        idx_err_ifgram = np.sum(C[idx_nonzero_closure, :] != 0., axis=0) >= 1.
        D = I[~idx_err_ifgram, :]
        A = np.vstack((A, D))

    # Eq 4.4 (Fattahi, 2015)
    try:
        A_inv = linalg.pinv2(A, rcond=rcond)
        L = np.zeros((A.shape[0], ifgram.shape[1]), np.float32)
        L[0:num_tri, :] = np.dot(C, ifgram) / (-2.*np.pi)
        U = np.round(np.dot(A_inv, L))

    except linalg.LinAlgError:
        pass

    ifgram_cor = ifgram + 2*np.pi*U
    return ifgram_cor, U


def run_unwrap_error_patch(ifgram_file, box=None, mask_file=None, ref_phase=None, fast_mode=False,
                           thres=0.1, dsNameIn='unwrapPhase'):
    """Estimate/Correct unwrapping error in ifgram stack on area defined by box.
    Parameters: ifgram_file : string, ifgramStack file
                box : tuple of 4 int, indicating areas to be read and analyzed
                mask_file : string, file name of mask file for pixels to be analyzed
                ref_pahse : 1D np.array in size of (num_ifgram,) phase value on reference pixel, because:
                    1) phase value stored in pysar is not reference yet
                    2) reference point may be out of box definition
                fast_mode : bool, apply zero jump constraint on ifgrams without unwrapping error.
                thres : float, threshold of non-zero phase closure to be identified as unwrapping error.
    Returns:    pha_data : 3D np.array in size of (num_ifgram_all, box[3]-box[2], box[2]-box[0]),
                    unwrapped phase value after error correction
    """
    # Basic info
    obj = ifgramStack(ifgram_file)
    obj.open(print_msg=False)
    num_ifgram = obj.numIfgram

    # Size Info - Patch
    if box:
        num_row = box[3] - box[1]
        num_col = box[2] - box[0]
    else:
        num_row = obj.length
        num_col = obj.width
    num_pixel = num_row * num_col

    C = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=True))
    print('number of interferograms: {}'.format(C.shape[1]))
    print('number of triangles: {}'.format(C.shape[0]))

    # read unwrapPhase
    pha_data_all = ifginv.read_unwrap_phase(obj, box, ref_phase,
                                            unwDatasetName=dsNameIn,
                                            dropIfgram=False)
    pha_data = np.array(pha_data_all[obj.dropIfgram, :])

    # mask of pixels to analyze
    mask = np.ones((num_pixel), np.bool_)
    print('number of pixels read: {}'.format(num_pixel))
    # mask 1. mask of water or area of interest
    if mask_file:
        dsNames = readfile.get_dataset_list(mask_file)
        dsName = [i for i in dsNames if i in ['waterMask', 'mask']][0]
        waterMask = readfile.read(mask_file, datasetName=dsName, box=box)[0].flatten()
        mask *= np.array(waterMask, np.bool_)
        del waterMask
        print('number of pixels left after mask: {}'.format(np.sum(mask)))

    # mask 2. mask of pixels without unwrap error: : zero phase closure on all triangles
    print('calculating phase closure of all possible triangles ...')
    pha_closure = np.dot(C, pha_data)
    pha_closure = np.abs(pha_closure - ut.wrap(pha_closure))       # Eq 4.2 (Fattahi, 2015)
    num_nonzero_closure = np.sum(pha_closure >= thres, axis=0)
    mask *= (num_nonzero_closure != 0.)
    del pha_closure
    print('number of pixels left after checking phase closure: {}'.format(np.sum(mask)))

    # mask summary
    num_pixel2proc = int(np.sum(mask))
    if num_pixel2proc > 0:
        ifgram = pha_data[:, mask]
        ifgram_cor = np.array(ifgram, np.float32)
        print('number of pixels to process: {} out of {} ({:.2f}%)'.format(num_pixel2proc, num_pixel,
                                                                           num_pixel2proc/num_pixel*100))

        # correcting unwrap error based on phase closure
        print('correcting unwrapping error ...')
        if fast_mode:
            ifgram_cor = correct_unwrap_error(ifgram, C, Dconstraint=False)[0]

        else:
            prog_bar = ptime.progressBar(maxValue=num_pixel2proc)
            for i in range(num_pixel2proc):
                ifgram_cor[:, i] = correct_unwrap_error(ifgram[:, i], C, Dconstraint=True)[0].flatten()
                prog_bar.update(i+1, every=10, suffix='{}/{}'.format(i+1, num_pixel2proc))
            prog_bar.close()

        pha_data[:, mask] = ifgram_cor
        pha_data_all[obj.dropIfgram, :] = pha_data

    pha_data_all = pha_data_all.reshape(num_ifgram, num_row, num_col)
    num_nonzero_closure = num_nonzero_closure.reshape(num_row, num_col)
    return pha_data_all, num_nonzero_closure


def run_unwrap_error_closure(inps, dsNameIn='unwrapPhase', dsNameOut='unwrapPhase_phaseClosure', fast_mode=False):
    """Run unwrapping error correction in network of interferograms using phase closure.
    Parameters: inps : Namespace of input arguments including the following:
                    ifgram_file : string, path of ifgram stack file
                    maskFile    : string, path of mask file mark the pixels to be corrected
                dsNameIn  : string, dataset name to read in
                dsnameOut : string, dataset name to write out
                fast_mode : bool, enable fast processing mode or not
    Returns:    inps.ifgram_file : string, path of corrected ifgram stack file
    """
    print('-'*50)
    print('Unwrapping Error Coorrection based on Phase Closure Consistency (Fattahi, 2015) ...')
    if fast_mode:
        print('fast mode: ON, the following asuumption is ignored for fast processing')
        print('\tzero phase jump constraint on ifgrams without unwrap error')
    else:
        print(('this step can be very slow'
               'you could terminate this and try --fast option for quick correction process.'))
    print('-'*50)

    obj = ifgramStack(inps.ifgram_file)
    obj.open()
    ref_phase = obj.get_reference_phase(unwDatasetName=dsNameIn, dropIfgram=False)

    num_nonzero_closure = np.zeros((obj.length, obj.width), dtype=np.int16)
    # split ifgram_file into blocks to save memory
    num_tri = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=True)).shape[0]
    length, width = obj.get_size()[1:3]
    box_list = ifginv.split_into_boxes(dataset_shape=(num_tri, length, width), chunk_size=200e6)
    num_box = len(box_list)
    for i in range(num_box):
        box = box_list[i]
        if num_box > 1:
            print('\n------- Processing Patch {} out of {} --------------'.format(i+1, num_box))

        # estimate/correct ifgram
        unw_cor, num_closure = run_unwrap_error_patch(inps.ifgram_file,
                                                      box=box,
                                                      mask_file=inps.maskFile,
                                                      ref_phase=ref_phase,
                                                      fast_mode=fast_mode,
                                                      dsNameIn=dsNameIn)
        num_nonzero_closure[box[1]:box[3], box[0]:box[2]] = num_closure

        # write ifgram
        write_hdf5_file_patch(inps.ifgram_file,
                              data=unw_cor,
                              box=box,
                              dsName=dsNameOut)

    # write number of nonzero phase closure into file
    num_file = 'numNonzeroPhaseClosure.h5'
    atr = dict(obj.metadata)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = '1'
    writefile.write(num_nonzero_closure, out_file=num_file, metadata=atr)
    print('writing >>> {}'.format(num_file))

    return inps.ifgram_file


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
    run_unwrap_error_closure(inps,
                             dsNameIn=inps.datasetNameIn,
                             dsNameOut=inps.datasetNameOut,
                             fast_mode=inps.fastMode)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
