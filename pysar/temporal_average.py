#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2018, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import time
import argparse
from pysar.objects import ifgramDatasetNames
from pysar.utils import utils as ut, readfile


#################################  Usage  ####################################
EXAMPLE = """example:
  temporal_average.py ifgramStack.h5 -d coherence -o avgSpatialCoherence.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Calculate temporal average/mean of multi-temporal datasets',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', type=str, help='input file with multi-temporal datasets')
    parser.add_argument('-d', '--ds', '--dataset', dest='datasetName', default=ifgramDatasetNames[1],
                        help='dataset name to be averaged, for file with multiple dataset family,\n'+
                        'e.g. ifgramStack.h5\n' +
                        'default: {}'.format(ifgramDatasetNames[1]))
    parser.add_argument('-o', '--outfile', help='output file name')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def check_output_filename(inps):
    ext = os.path.splitext(inps.file)[1]
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    if not inps.outfile:
        if k == 'ifgramStack':
            if inps.datasetName == 'coherence':
                inps.outfile = 'avgSpatialCoherence.h5'
            elif inps.datasetName == 'unwrapPhase':
                inps.outfile = 'avgPhaseVelocity.h5'
            else:
                inps.outfile = 'avg{}.h5'.format(inps.datasetName)
        elif k == 'timeseries':
            processMark = os.path.basename(inps.file).split('timeseries')[1].split(ext)[0]
            inps.outfile = 'avgDisplacement{}.h5'.format(processMark)
        else:
            inps.outfile = 'avg{}.h5'.format(inps.file)
    print('output file: {}'.format(inps.outfile))
    return inps.outfile


#############################  Main Function  ################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    inps.outfile = check_output_filename(inps)
    inps.outfile = ut.temporal_average(inps.file,
                                       datasetName=inps.datasetName,
                                       outFile=inps.outfile)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return inps.outfile


##############################################################################
if __name__ == '__main__':
    main()
