############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os

import numpy as np

from mintpy.utils import plot as pp, readfile, utils as ut


######################################################################################################
def read_template2inps(templateFile, inps):
    """Update inps with mintpy.residualRMS.* option from templateFile"""
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(templateFile))
    template = readfile.read_template(templateFile)
    template = ut.check_template_auto_value(template)

    prefix = 'mintpy.residualRMS.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        # false/none values are valid inputs, thus, should be passed here without an if check
        if key in ['maskFile', 'deramp']:
            inpsDict[key] = value
        elif key in ['cutoff']:
            inpsDict[key] = float(value)
    return inps


def analyze_rms(date_list, rms_list, inps):
    # reference date
    ref_idx = np.argmin(rms_list)
    print('-'*50+'\ndate with min RMS: {} - {:.4f}'.format(date_list[ref_idx],
                                                           rms_list[ref_idx]))
    ref_date_file = 'reference_date.txt'
    if ut.run_or_skip(out_file=ref_date_file,
                      in_file=[inps.timeseries_file, inps.maskFile, inps.template_file],
                      readable=False) == 'run':
        with open(ref_date_file, 'w') as f:
            f.write(date_list[ref_idx]+'\n')
        print('save date to file: '+ref_date_file)

    # exclude date(s) - outliers
    # equivalent calculation using numpy assuming Gaussian distribution as:
    # rms_threshold = np.median(rms_list) / .6745 * inps.cutoff
    rms_threshold = ut.median_abs_deviation_threshold(rms_list, center=0., cutoff=inps.cutoff)

    ex_idx = [rms_list.index(i) for i in rms_list if i > rms_threshold]
    print('-'*50)
    print(f'date(s) with RMS > {inps.cutoff} * median RMS ({rms_threshold:.4f})')

    ex_date_file = 'exclude_date.txt'
    if ex_idx:
        # print
        for i in ex_idx:
            print(f'{date_list[i]} - {rms_list[i]:.4f}')
        # save to text file
        with open(ex_date_file, 'w') as f:
            for i in ex_idx:
                f.write(date_list[i]+'\n')
        print('save date(s) to file: '+ex_date_file)
    else:
        print('None.')
        if os.path.isfile(ex_date_file):
            os.remove(ex_date_file)

    return inps


######################################################################################################
def run_timeseries_rms(inps):

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # calculate timeseries of residual Root Mean Square
    inps.rms_list, inps.date_list, inps.rms_file = ut.get_residual_rms(
        inps.timeseries_file,
        mask_file=inps.maskFile,
        ramp_type=inps.deramp,
    )

    # analyze RMS: generate reference/exclude_date.txt files
    analyze_rms(inps.date_list, inps.rms_list, inps)

    # plot RMS
    if '--figsize' not in inps.argv and len(inps.date_list) > 120:
        fig_wid = min(15, inps.fig_size[0] * len(inps.date_list) / 120)
        inps.fig_size[0] = float(f'{fig_wid:.1f}')

    pp.plot_timeseries_rms(
        rms_file=inps.rms_file,
        cutoff=inps.cutoff,
        out_fig=os.path.splitext(inps.rms_file)[0]+'.pdf',
        disp_fig=False,
        fig_size=inps.fig_size,
        tick_year_num=inps.tick_year_num,
    )

    return
