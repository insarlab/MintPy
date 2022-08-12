############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import h5py
from mintpy.utils import readfile


#################################  Usage  ####################################
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
