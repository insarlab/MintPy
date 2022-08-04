#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import sys
import shutil
from mintpy.utils import readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser


##################################################################################################
NOTE = """
  For each binary file (unwrapped/wrapped interferogram, spatial coherence file), there are 2 .rsc files:
  1) basic metadata file and 2) baseline parameter file. This script find those two rsc files based on
  input binary file name, and merge those two metadata files into one.

  For example, if input binary file is filt_100901-110117-sim_HDR_4rlks_c10.unw, this script will find
  1) filt_100901-110117-sim_HDR_4rlks_c10.unw.rsc and 2) 100901-110117_baseline.rsc and merge 1) and 2) into
  one file: filt_100901-110117-sim_HDR_4rlks_c10.unw.rsc
"""

EXAMPLE = """example:
  prep_roipac.py  filt_100901-110117-sim_HDR_4rlks_c10.unw
  prep_roipac.py  ./interferograms/*/filt_*.unw
  prep_roipac.py  ./interferograms/*/filt_*rlks.cor
  prep_roipac.py  ./interferograms/*/filt_*rlks.int
  prep_roipac.py  ./interferograms/*/filt_*_snap_connect.byt
"""


def create_parser(subparsers=None):
    synopsis = 'Prepare attributes file for ROI_PAC products.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='Gamma file(s)')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.file = ut.get_file_list(inps.file, abspath=True)

    # Check input file type
    ext_list = ['.unw', '.cor', '.int', '.byt', '.hgt', '.dem', '.trans']
    ext = os.path.splitext(inps.file[0])[1]
    if ext not in ext_list:
        msg = 'unsupported input file extension: {}'.format(ext)
        msg += '\nsupported file extensions: {}'.format(ext_list)
        raise ValueError(msg)

    return inps


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


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    ext = os.path.splitext(inps.file[0])[1]
    if ext in ['.unw', '.cor', '.int', '.byt']:
        for fname in inps.file:
            extract_metadata(fname)
    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
