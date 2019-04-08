#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import sys
import argparse
import numpy as np
from pysar.objects import ifgramStack, timeseries
from pysar.utils import readfile, writefile, ptime, utils as ut
from pysar import ifgram_inversion as ifginv


######################################################################################################
REFERENCE = """reference:
  Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for
  phase unwrapping of multitemporal differential SAR interferograms, IEEE-TGRS, 44(9), 2374-2383.
"""

EXAMPLE = """example:
  temporal_coherence.py  INPUTS/ifgramStack.h5  timeseries.h5  -n invertIfgramNum.h5
  temporal_coherence.py  INPUTS/ifgramStack.h5  timeseries.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generates temporal coherence map.\n'+
                                     '[Not recommended, use temporalCoherence.h5 from ifgram_inversion.py.]',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('ifgram_file', help='unwrapped interferograms file')
    parser.add_argument('timeseries_file', help='timeseries file after network inversion')

    parser.add_argument('-n', dest='ifg_num_file',
                        help='file for number of interferograms used for network inversion.')
    parser.add_argument('-o', '--outfile', default='temporalCoherence.h5',
                        help='output file name for temporal coherence')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


######################################################################################################
def calculate_temporal_coherence(ifgram_file, timeseries_file, ifg_num_file=None, chunk_size=100e6):
    """Calculate temporal coherence based on input timeseries file and interferograms file
    Parameters: ifgram_file : str, path of interferograms file
                timeseries_file : str, path of time series file
                ifg_num_file : str, path of file for number of interferograms used in inversion.
                chunk_size : float
    Returns:    temp_coh : 2D np.array, temporal coherence in float32
    """

    # get box list and size info
    box_list = ifginv.split_into_boxes(ifgram_file, chunk_size=chunk_size)
    num_box = len(box_list)

    stack_shape = ifgramStack(ifgram_file).get_size()
    temp_coh = np.zeros(stack_shape[1:3], np.float32)
    for i in range(num_box):
        if num_box > 1:
            print('\n------- Processing Patch %d out of %d --------------' % (i+1, num_box))
        box = box_list[i]
        temp_cohi = calculate_temporal_coherence_patch(ifgram_file, 
                                                       timeseries_file,
                                                       box=box,
                                                       ifg_num_file=ifg_num_file)
        temp_coh[box[1]:box[3], box[0]:box[2]] = temp_cohi
    return temp_coh


def calculate_temporal_coherence_patch(ifgram_file, timeseries_file, box=None, ifg_num_file=None):
    atr = readfile.read_attribute(timeseries_file)
    if not box:
        box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))

    # Read timeseries data
    ts_obj = timeseries(timeseries_file)
    ts_obj.open(print_msg=False)
    print('reading timeseries data from file: {}'.format(timeseries_file))
    ts_data = ts_obj.read(box=box, print_msg=False).reshape(ts_obj.numDate, -1)
    ts_data = ts_data[1:, :]
    ts_data *= -4*np.pi/float(atr['WAVELENGTH'])

    # Read ifgram data
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    A = stack_obj.get_design_matrix4timeseries(stack_obj.get_date12_list(dropIfgram=True))[0]
    print('reading unwrapPhase data from file: {}'.format(ifgram_file))
    ifgram_data = stack_obj.read(datasetName='unwrapPhase', box=box).reshape(A.shape[0], -1)
    ref_value = stack_obj.get_reference_phase(dropIfgram=True).reshape((-1, 1))
    ifgram_data -= np.tile(ref_value, (1, ifgram_data.shape[1]))

    ifgram_diff = ifgram_data - np.dot(A, ts_data)
    del ts_data

    pixel_num = ifgram_data.shape[1]
    temp_coh = np.zeros((pixel_num), np.float32)
    # (fast) nasty solution, which used all phase value including invalid zero phase
    if not ifg_num_file:
        temp_coh = np.abs(np.sum(np.exp(1j*ifgram_diff), axis=0)) / ifgram_diff.shape[0]

    # (slow) same solution as ifgram_inversion.py, considering:
    #   1) invalid zero phase in ifgram
    #   2) design matrix rank deficiency.
    else:
        print('considering different number of interferograms used in network inversion for each pixel')
        ifg_num_map = readfile.read(ifg_num_file, box=box)[0].flatten()
        prog_bar = ptime.progressBar(maxValue=pixel_num)
        for i in range(pixel_num):
            if ifg_num_map[i] > 0:
                idx = ifgram_data[:, i] != 0.
                temp_diff = ifgram_diff[idx, i]
                temp_coh[i] = np.abs(np.sum(np.exp(1j*temp_diff), axis=0)) / temp_diff.shape[0]
            prog_bar.update(i+1, every=1000, suffix='{}/{}'.format(i+1, pixel_num))
        prog_bar.close()

    temp_coh = np.reshape(temp_coh, (box[3]-box[1], box[2]-box[0]))
    return temp_coh


######################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    temp_coh = calculate_temporal_coherence(ifgram_file=inps.ifgram_file,
                                            timeseries_file=inps.timeseries_file,
                                            ifg_num_file=inps.ifg_num_file)

    # write file
    atr = readfile.read_attribute(inps.timeseries_file)
    atr['FILE_TYPE'] = 'temporalCoherence'
    atr['UNIT'] = '1'
    writefile.write(temp_coh, out_file=inps.outfile, metadata=atr)
    return inps.outfile


######################################################################################################
if __name__ == '__main__':
    main()
