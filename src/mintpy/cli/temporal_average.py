#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Antonio Valentino, Aug 2016        #
############################################################


import os
import sys
import time

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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default:
    if not inps.outfile:
        inps.outfile = auto_output_filename(inps)

    return inps


##############################  Sub Function  ################################
def auto_output_filename(inps):
    """Default output file name"""
    from mintpy.utils import readfile

    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']
    fdir = os.path.dirname(inps.file)
    fbase = os.path.basename(inps.file)

    if ftype == 'ifgramStack':
        if inps.datasetName == 'coherence':
            inps.outfile = 'avgSpatialCoh.h5'
        elif 'unwrapPhase' in inps.datasetName:
            inps.outfile = 'avgPhaseVelocity.h5'
        else:
            inps.outfile = f'avg{inps.datasetName}.h5'

    elif ftype == 'timeseries':
        suffix = os.path.splitext(fbase)[0].split('timeseries')[1]
        inps.outfile = f'avgDisp{suffix}.h5'

    else:
        inps.outfile = f'avg{fbase}'

    inps.outfile = os.path.join(fdir, inps.outfile)
    print(f'output file: {inps.outfile}')

    return inps.outfile


def run_or_skip(inps):
    import h5py

    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file vs input dataset
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print(f'1) output file {inps.outfile} NOT exist.')
    else:
        print(f'1) output file {inps.outfile} already exists.')
        with h5py.File(inps.file, 'r') as f:
            atr = f[inps.datasetName].attrs
            ti = float(atr.get('MODIFICATION_TIME', os.path.getmtime(inps.file)))
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print(f'2) output file is NOT newer than input dataset: {inps.datasetName}.')
        else:
            print(f'2) output file is newer than input dataset: {inps.datasetName}.')

    # result
    print(f'run or skip: {flag}.')
    return flag


#############################  Main Function  ################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.utils import utils1 as ut

    # run or skip
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return

    # run
    start_time = time.time()
    ut.temporal_average(
        inps.file,
        datasetName=inps.datasetName,
        outFile=inps.outfile,
    )

    # used time
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs\n')

    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
