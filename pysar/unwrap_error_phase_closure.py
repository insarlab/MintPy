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
try:
    from cvxopt import matrix, sparse
except ImportError:
    raise ImportError('Cannot import cvxopt')
from pysar.objects import ifgramStack
from pysar.utils import ptime, readfile, writefile, utils as ut
from pysar.utils.solvers import l1reg_lstsq
from pysar import ifgram_inversion as ifginv


key_prefix = 'pysar.unwrapError.'

##########################################################################################
EXAMPLE = """Example:
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  -t pysarApp_template.txt  --update
  unwrap_error_phase_closure.py  ./INPUTS/ifgramStack.h5  --water-mask waterMask.h5 --update
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


##########################################################################################
def run_unwrap_error_closure_patch(ifgram_file, box=None, mask_file=None, ref_phase=None, alpha=1e-2,
                                   dsNameIn='unwrapPhase'):
    """Estimate/Correct unwrapping error in ifgram stack on area defined by box.
    Parameters: ifgram_file : string, ifgramStack file
                box : tuple of 4 int, indicating areas to be read and analyzed
                mask_file : string, file name of mask file for pixels to be analyzed
                ref_pahse : 1D np.array in size of (num_ifgram,) phase value on reference pixel, because:
                    1) phase value stored in pysar is not reference yet
                    2) reference point may be out of box definition
                alpha : nonnegative float, tradeoff between L1- and L2-norm terms.
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

    C = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=True)).astype(float)
    print('number of interferograms: {}'.format(C.shape[1]))
    print('number of triangles: {}'.format(C.shape[0]))

    # read unwrapPhase
    pha_data = ifginv.read_unwrap_phase(obj, box, ref_phase,
                                        unwDatasetName=dsNameIn,
                                        dropIfgram=False)
    print('calculating the integer ambiguity of the closure phase of all possible triplets ...')
    closure_pha = np.dot(C, np.array(pha_data[obj.dropIfgram, :]))
    closure_int = np.round((closure_pha - ut.wrap(closure_pha)) / (2.*np.pi))
    num_closure_int = np.sum(closure_int != 0., axis=0)
    del closure_pha

    # mask of pixels to analyze
    mask = np.ones((num_pixel), np.bool_)
    print('number of pixels read: {}'.format(num_pixel))
    # mask 1. mask of water or area of interest
    if mask_file:
        dsNames = readfile.get_dataset_list(mask_file)
        dsName = [i for i in dsNames if i in ['waterMask', 'mask', 'landMask']][0]
        waterMask = readfile.read(mask_file, datasetName=dsName, box=box)[0].flatten()
        mask *= np.array(waterMask, np.bool_)
        del waterMask
        print('number of pixels left after mask: {}'.format(np.sum(mask)))
    # mask 2. mask of pixels without unwrap error: : zero phase closure on all triangles
    mask *= (num_closure_int != 0.)
    # mask summary
    num_pixel2proc = int(np.sum(mask))
    print('number of pixels left after checking phase closure: {}'.format(num_pixel2proc))

    if num_pixel2proc > 0:
        print('number of pixels to process: {} out of {} ({:.2f}%)'.format(num_pixel2proc,
                                                                           num_pixel,
                                                                           num_pixel2proc/num_pixel*100))
        # prepare matrix in cvxopt format
        C = matrix(C)
        closure_int = matrix(closure_int[:, mask])

        # correcting unwrap error based on phase closure
        print('solving the phase-unwrapping integer ambiguity')
        print('\tusing L1-norm regularized least squares approximation (LASSO) ...')
        pha_offset = np.zeros((np.sum(obj.dropIfgram), num_pixel2proc), np.float32)
        prog_bar = ptime.progressBar(maxValue=num_pixel2proc)
        for i in range(num_pixel2proc):
            U = np.round(l1reg_lstsq(-C, closure_int[:, i], lambd=alpha))
            pha_offset[:, i] = 2.*np.pi*U.flatten()
            prog_bar.update(i+1, every=10, suffix='{}/{}'.format(i+1, num_pixel2proc))
        prog_bar.close()

        # update phase value on selected interferograms and pixels
        pha_data[obj.dropIfgram, :][:, mask] += pha_offset

    pha_data = pha_data.reshape(num_ifgram, num_row, num_col)
    num_closure_int = num_closure_int.reshape(num_row, num_col)
    return pha_data, num_closure_int


def run_unwrap_error_closure(inps, dsNameIn='unwrapPhase', dsNameOut='unwrapPhase_phaseClosure'):
    """Run unwrapping error correction in network of interferograms using phase closure.
    Parameters: inps : Namespace of input arguments including the following:
                    ifgram_file : string, path of ifgram stack file
                    maskFile    : string, path of mask file mark the pixels to be corrected
                dsNameIn  : string, dataset name to read in
                dsnameOut : string, dataset name to write out
    Returns:    inps.ifgram_file : string, path of corrected ifgram stack file
    """
    print('-'*50)
    print('Unwrapping Error Coorrection based on Phase Closure Consistency (Fattahi, 2015) ...')
    obj = ifgramStack(inps.ifgram_file)
    obj.open()
    ref_phase = obj.get_reference_phase(unwDatasetName=dsNameIn, dropIfgram=False)

    num_closure_int = np.zeros((obj.length, obj.width), dtype=np.int16)
    # split ifgram_file into blocks to save memory
    num_tri = obj.get_design_matrix4triplet(obj.get_date12_list(dropIfgram=True)).shape[0]
    length, width = obj.get_size()[1:3]
    box_list = ifginv.split_into_boxes(dataset_shape=(num_tri, length, width), chunk_size=20e6)
    num_box = len(box_list)
    for i in range(num_box):
        box = box_list[i]
        if num_box > 1:
            print('\n------- Processing Patch {} out of {} --------------'.format(i+1, num_box))

        # estimate/correct ifgram
        unw_cor, num_int = run_unwrap_error_closure_patch(inps.ifgram_file,
                                                          box=box,
                                                          mask_file=inps.waterMaskFile,
                                                          ref_phase=ref_phase,
                                                          dsNameIn=dsNameIn)
        num_closure_int[box[1]:box[3], box[0]:box[2]] = num_int

        # write ifgram
        write_hdf5_file_patch(inps.ifgram_file,
                              data=unw_cor,
                              box=box,
                              dsName=dsNameOut)

    # write number of nonzero phase closure into file
    num_file = 'numNonzeroClosureInt.h5'
    atr = dict(obj.metadata)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = '1'
    writefile.write(num_closure_int, out_file=num_file, metadata=atr)
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
    run_unwrap_error_closure(inps, dsNameIn=inps.datasetNameIn, dsNameOut=inps.datasetNameOut)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return inps.ifgram_file


####################################################################################################
if __name__ == '__main__':
    main()
