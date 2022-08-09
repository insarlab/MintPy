############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
from mintpy.objects import RAMP_LIST
from mintpy.utils import readfile


# key configuration parameter name
configKeys = [
    'mintpy.deramp',
    'mintpy.deramp.maskFile',
]


###########################################################################################
def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print('1) output file {} NOT found.'.format(inps.outfile))
    else:
        print('1) output file {} already exists.'.format(inps.outfile))
        infiles = [inps.file]
        if inps.mask_file:
            infiles.append(inps.mask_file)
        ti = max(os.path.getmtime(i) for i in infiles)
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output file is NOT newer than input file: {}.'.format(infiles))
        else:
            print('2) output file is newer than input file: {}.'.format(infiles))

    # check configuration
    if flag == 'skip':
        iDict = {}
        iDict['mintpy.deramp'] = inps.surface_type
        iDict['mintpy.deramp.maskFile'] = inps.mask_file
        atr = readfile.read_attribute(inps.outfile)
        if any(str(iDict[key]) != atr.get(key, 'None') for key in configKeys):
            flag = 'run'
            print('3) NOT all key configuration parameters are the same:{}'.format(configKeys))
        else:
            print('3) all key configuration parameters are the same:{}'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag
