############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import shutil
from mintpy.utils import readfile, writefile


######################################## Sub Functions ############################################
def extract_metadata(fname):
    """Read/extract attributes from ROI_PAC .unw, .int, .cor file.

    For each unwrapped interferogram or spatial coherence file, there are 2 .rsc files:
        basic metadata file and baseline parameter file. 
        e.g. filt_100901-110117-sim_HDR_4rlks_c10.unw
             filt_100901-110117-sim_HDR_4rlks_c10.unw.rsc
             100901-110117_baseline.rsc
    Inputs:
        fname : string, ROI_PAC interferogram filename or path,
                i.e. /KujuT422F650AlosA/filt_100901-110117-sim_HDR_4rlks_c10.unw
    Outputs:
        atr : dict, Attributes dictionary
    """
    # 1. Read basic metadata file
    basic_rsc_file = fname+'.rsc'
    if not os.path.isfile(basic_rsc_file) and fname.endswith('_snap_connect.byt'):
        unw_rsc_file = '{}.unw.rsc'.format(fname.split('_snap_connect.byt')[0])
        print('copy {} to {}'.format(unw_rsc_file, basic_rsc_file))
        shutil.copy2(unw_rsc_file, basic_rsc_file)
    basic_dict = readfile.read_roipac_rsc(basic_rsc_file)

    # return if baseline attributes are already existed.
    if 'P_BASELINE_TOP_HDR' in basic_dict.keys():
        return basic_rsc_file

    atr = {}
    atr['PROCESSOR'] = 'roipac'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]

    # 2. Read baseline metadata file
    date1, date2 = basic_dict['DATE12'].split('-')
    baseline_rsc_file = os.path.dirname(fname)+'/'+date1+'_'+date2+'_baseline.rsc'
    baseline_dict = readfile.read_roipac_rsc(baseline_rsc_file)

    # 3. Merge
    atr.update(basic_dict)
    atr.update(baseline_dict)

    # Write to rsc file
    basic_rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(basic_rsc_file)
    except:
        atr_orig = dict()
    if not set(atr.items()).issubset(set(atr_orig.items())):
        atr_out = {**atr_orig, **atr}
        print('merging {} into {} '.format(os.path.basename(baseline_rsc_file),
                                           os.path.basename(basic_rsc_file)))
        writefile.write_roipac_rsc(atr_out, out_file=basic_rsc_file)
    return basic_rsc_file
