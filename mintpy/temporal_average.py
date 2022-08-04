#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import sys
import time
import h5py
from mintpy.utils import utils as ut, readfile
from mintpy.utils.arg_utils import create_argument_parser


#################################  Usage  ####################################
EXAMPLE = """example:
  temporal_average.py ./inputs/ifgramStack.h5 -d unwrapPhase -o avgPhaseVelocity.h5
  temporal_average.py ./inputs/ifgramStack.h5 -d coherence   -o avgSpatialCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate temporal average (stacking) of multi-temporal datasets'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='input file with multi-temporal datasets')
    parser.add_argument('-d', '--ds', '--dataset', dest='datasetName', default='coherence',
                        help='dataset name to be averaged, for file with multiple dataset family,\n'+
                        'e.g. ifgramStack.h5\n' +
                        'default: coherence')
    parser.add_argument('-o', '--outfile', help='output file name')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update checking for --nonzero option.')
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
                inps.outfile = 'avgSpatialCoh.h5'
            elif 'unwrapPhase' in inps.datasetName:
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


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file vs input dataset
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print('1) output file {} NOT exist.'.format(inps.outfile))
    else:
        print('1) output file {} already exists.'.format(inps.outfile))
        with h5py.File(inps.file, 'r') as f:
            ti = float(f[inps.datasetName].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.file)))
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output file is NOT newer than input dataset: {}.'.format(inps.datasetName))
        else:
            print('2) output file is newer than input dataset: {}.'.format(inps.datasetName))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


#############################  Main Function  ################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    inps.outfile = check_output_filename(inps)

    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    ut.temporal_average(inps.file, datasetName=inps.datasetName, outFile=inps.outfile)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs\n'.format(m, s))
    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
