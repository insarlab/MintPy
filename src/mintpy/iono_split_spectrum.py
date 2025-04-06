#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2024                           #
############################################################


from mintpy.cli import diff, ifgram_inversion, modify_network, reference_point
from mintpy.utils import utils as ut


###############################################################
def run_iono_split_spectrum(inps):
    """Estimate (and correct) ionospheric delay time-series."""

    # [in the future] load the low resolution iono here [to be independent from load_data.py]
    # 1. modify the network of iono stack [to exclude pairs with failed split spectrum estimations]
    print('\n'+'-'*80)
    print('Modify the network of ionospheric delay stack via modify_network.py ...')
    if inps.excludeDate or inps.excludeDate12:
        cmd = f'modify_network.py {inps.iono_stack_file}'
        cmd += ' --ex-date ' + ' '.join(x for x in inps.excludeDate) if inps.excludeDate else ''
        cmd += ' --ex-date12 ' + ' '.join(x for x in inps.excludeDate12) if inps.excludeDate12 else ''
        print(cmd)
        modify_network.main(cmd.split()[1:])
    else:
        print('No configuration found, skip modifying the network and continue.')

    # 2. select the spatial reference point for the iono stack
    print('\n'+'-'*80)
    print('Apply the spatial reference to the network of ionospheric delay stack via reference_point.py ...')
    cmd = f'reference_point.py {inps.iono_stack_file} -t {inps.template_file}'
    print(cmd)
    reference_point.main(cmd.split()[1:])

    # 3. estimate iono time-series
    # hardwire "--dset unwrapPhase" to ignore dataset name change from unwrapping error correction options
    print('\n'+'-'*80)
    print('Estimate ionospheric delay time-series via ifgram_inversion.py ...')
    cmd = f'ifgram_inversion.py {inps.iono_stack_file} --dset unwrapPhase --weight-func no --update'
    print(cmd)
    ifgram_inversion.main(cmd.split()[1:])

    # 4. correct iono delay from displacement time-series via diff.py
    if inps.dis_file:
        print('\n'+'-'*80)
        print('Apply ionospheric correction to displacement file via diff.py ...')
        if ut.run_or_skip(inps.cor_dis_file, [inps.dis_file, inps.iono_file]) == 'run':
            # [in the future] diff.py should handle different resolutions between
            # the iono time-series and displacement time-series files
            cmd = f'diff.py {inps.dis_file} {inps.iono_file} -o {inps.cor_dis_file} --force'
            print(cmd)
            diff.main(cmd.split()[1:])
        else:
            print(f'Skip re-applying and use existed corrected displacement file: {inps.cor_dis_file}.')
    else:
        print('No input displacement file, skip correcting ionospheric delays.')

    return
