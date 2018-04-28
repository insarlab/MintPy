#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2017, Zhang Yunjun, Yunmeng Cao             #
# Author:  Zhang Yunjun, Yunmeng Cao                       #
############################################################


import os
import sys
import re
import argparse
import numpy as np
from pysar.utils import readfile, writefile, datetime as ptime, utils as ut


##################################################################################################
EXAMPLE = """example:
  prep_gamma.py  diff_filt_HDR_130118-130129_4rlks.unw
  prep_gamma.py  IFGRAM*/diff_*rlks.unw
  prep_gamma.py  IFGRAM*/filt_*rlks.cor
  prep_gamma.py  IFGRAM*/diff_*rlks.int
  prep_gamma.py  sim_150911-150922.hgt_sim
  prep_gamma.py  sim_150911-150922.utm.dem
  prep_gamma.py  sim_150911-150922.UTM_TO_RDC
"""

DESCRIPTION = """
  For each interferogram, including unwrapped/wrapped interferograms and coherence, 3 metadata files are required:
  1) master .par file, e.g. 130118.amp.par
  2) slave  .par file, e.g. 130129.amp.par
  3) interferogram .off file, e.g. 130118-130129_4rlks.off

  Other metadata files are recommended and can be generated from the above 3 if not existed, more specifically:
  4) baseline files, e.g. 130118-130129_4rlks.baseline and 130118-130129_4rlks.base_perp
      It can be generated from file 1-3 with Gamma command base_orbit and base_perp)
  5) corner files, e.g. 130118_4rlks.amp.corner_full and 130118_4rlks.amp.corner
      It can be generated from file 1 with Gamma command SLC_corners)

  This script will read all these files (generate them if not existed), merge them into one, convert their name from
  Gamma style to ROI_PAC style, and write to an metadata file, same name as input binary data file with suffix .rsc,
  e.g. diff_filt_HDR_130118-130129_4rlks.unw.rsc


  For DEM file in radar/geo coordinates (.utm.dem/.hgt_sim) and lookup table file for geocoding (.UTM_TO_RDC), 2 metadata
  files are required:
  1) .par file, contains geo coordinates info, for DEM in geo coordinates and lookup table, e.g. sim_150911-150922.utm.dem.par
  2) .diff_par file, contains radar coordinates info, for DEM in radar coordinates, e.g. sim_150911-150922.diff_par


  Here is an example of how your Gamma files should look like, after all interferograms/SLCs are coregistered:
  For each interferogram, 5 files are needed:
      diff_130118-130129_4rlks.unw
      filt_130118-130129_4rlks.cor
      130118-130129_4rlks.off
      130118_4rlks.amp.par
      130129_4rlks.amp.par
  For each dataset, only one sim* folder with 5 files are needed, 
      sim_130118-130129.hgt_sim
      sim_130118-130129.diff_par
      sim_130118-130129.utm.dem
      sim_130118-130129.utm.dem.par
      sim_130118-130129.UTM_TO_RDC

  Notes: both - and _ are supported; 
         both YYMMDD and YYYYMMDD naming are also supported;
         if no multilooking applied, do not "_4rlks" in your file names.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Prepare attributes file for Gamma product for PySAR.\n'+
                                     DESCRIPTION,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='Gamma file(s)')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false', default=True,
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


######################################## Sub Functions ############################################
def get_perp_baseline(m_par_file, s_par_file, off_file, atr_dict={}):
    """Get perpendicular baseline info from master/slave par file and off file.
    Parameters: m_par_file : str, path, master parameter file, i.e. 130118_4rlks.amp.par
                s_par_file : str, path, slave  parameter file, i.e. 130129_4rlks.amp.oar
                off_file   : str, path, interferogram off file, i.e. 130118-130129_4rlks.off
                atr_dict   : dict, optional, attributes dictionary
    Returns:  bperp : str, perpendicular baseline for pixel at [0,0]
    """
    # Get Path Info
    off_file = os.path.abspath(off_file)
    file_dir = os.path.dirname(off_file)
    file_base = os.path.splitext(os.path.basename(off_file))[0]

    base_perp_file = file_dir+'/'+file_base+'.base_perp'
    baseline_file = file_dir+'/'+file_base+'.baseline'

    # Call Gamma command to calculate Bperp
    if not os.path.isfile(base_perp_file) or os.stat(base_perp_file).st_size == 0:
        if not os.path.isfile(baseline_file) or os.stat(baseline_file).st_size == 0:
            baseCmd = 'base_orbit {} {} {}'.format(m_par_file,
                                                   s_par_file,
                                                   baseline_file)
            print(baseCmd)
            os.system(baseCmd)
        bperpCmd = 'base_perp {} {} {} > {}'.format(baseline_file,
                                                    m_par_file,
                                                    off_file,
                                                    base_perp_file)
        print(bperpCmd)
        os.system(bperpCmd)

    # Read bperp txt file
    bperp_list = []

    f = open(base_perp_file, 'r')
    lines = f.readlines()
    f.close()

    start_line_idx = [lines.index(i)+2 for i in lines if 'bpara       bperp       blen' in i][0]
    for i in range(start_line_idx, len(lines)):
        c = lines[i].strip().split()
        if len(c) == 9:
            bperp_list.append(float(c[7]))

    atr_dict['P_BASELINE_TOP_HDR'] = str(bperp_list[0])
    atr_dict['P_BASELINE_BOTTOM_HDR'] = str(bperp_list[-1])

    return atr_dict


def get_lalo_ref(m_par_file, atr_dict={}):
    """Extract LAT/LON_REF1/2/3/4 from corner file, e.g. 130118_4rlks.amp.corner.
    If it's not existed, call Gamma script - SLC_corners - to generate it from SLC par file
        e.g. 130118_4rlks.amp.par

    Parameters: m_par_file : str, path, master date parameter file, i.e. 130118_4rlks.amp.par
                atr_dict   : dict, optional, attributes dictionary
    Returns:  lalo_ref
    """
    m_par_file = os.path.abspath(m_par_file)
    m_corner_file = os.path.splitext(m_par_file)[0]+'.corner'
    m_corner_full_file = m_corner_file+'_full'

    # Call Gamma command to calculate LAT/LON_REF
    if not os.path.isfile(m_corner_file):
        if not os.path.isfile(m_corner_full_file):
            cornerCmd = 'SLC_corners {} > {}'.format(m_par_file,
                                                     m_corner_full_file)
            print(cornerCmd)
            os.system(cornerCmd)
        extractCmd = "awk 'NR==3,NR==6 {print $3,$6} ' {} > {}".format(m_corner_full_file,
                                                                       m_corner_file)
        print(extractCmd)
        os.system(extractCmd)

    # Read corner txt file
    lalo_ref = np.loadtxt(m_corner_file, dtype=bytes).astype(str)
    atr_dict['LAT_REF1'] = lalo_ref[0, 0]
    atr_dict['LAT_REF2'] = lalo_ref[1, 0]
    atr_dict['LAT_REF3'] = lalo_ref[2, 0]
    atr_dict['LAT_REF4'] = lalo_ref[3, 0]
    atr_dict['LON_REF1'] = lalo_ref[0, 1]
    atr_dict['LON_REF2'] = lalo_ref[1, 1]
    atr_dict['LON_REF3'] = lalo_ref[2, 1]
    atr_dict['LON_REF4'] = lalo_ref[3, 1]

    return atr_dict


def extract_metadata4interferogram(fname):
    """Read/extract attributes for PySAR from Gamma .unw, .cor and .int file
    Inputs:
        fname : str, Gamma interferogram filename or path, i.e. /PopoSLT143TsxD/diff_filt_HDR_130118-130129_4rlks.unw
    Output:
        atr : dict, Attributes dictionary
    """
    file_dir = os.path.dirname(fname)
    file_basename = os.path.basename(fname)

    rsc_file = fname+'.rsc'
    # if os.path.isfile(rsc_file):
    #    return rsc_file

    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]

    # Get info: date12, num of loooks
    try:
        date12 = str(re.findall('\d{8}[-_]\d{8}', file_basename)[0])
    except:
        date12 = str(re.findall('\d{6}[-_]\d{6}', file_basename)[0])
    m_date, s_date = date12.replace('-', '_').split('_')
    atr['DATE12'] = ptime.yymmdd(m_date)+'-'+ptime.yymmdd(s_date)
    lks = os.path.splitext(file_basename.split(date12)[1])[0]

    # Read .off and .par file
    off_file = file_dir+'/*'+date12+lks+'.off'
    m_par_file = [file_dir+'/*'+m_date+lks+i for i in ['.amp.par', '.ramp.par']]
    s_par_file = [file_dir+'/*'+s_date+lks+i for i in ['.amp.par', '.ramp.par']]

    try:
        off_file = ut.get_file_list(off_file)[0]
    except:
        print('\nERROR: Can not find .off file, it supposed to be like: '+off_file)
    try:
        m_par_file = ut.get_file_list(m_par_file)[0]
    except:
        print('\nERROR: Can not find master date .par file, it supposed to be like: '+m_par_file)
    try:
        s_par_file = ut.get_file_list(s_par_file)[0]
    except:
        print('\nERROR: Can not find slave date .par file, it supposed to be like: '+s_par_file)

    # print 'read '+m_par_file
    # print 'read '+off_file
    par_dict = readfile.read_gamma_par(m_par_file)
    off_dict = readfile.read_gamma_par(off_file)
    atr.update(par_dict)
    atr.update(off_dict)

    # Perp Baseline Info
    # print 'extract baseline info from %s, %s and %s file' % (m_par_file, s_par_file, off_file)
    atr = get_perp_baseline(m_par_file, s_par_file, off_file, atr)

    # LAT/LON_REF1/2/3/4
    # print 'extract LAT/LON_REF1/2/3/4 from '+m_par_file
    atr = get_lalo_ref(m_par_file, atr)

    # Write to .rsc file
    # print 'writing >>> '+rsc_file
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = None
    keyList = [i for i in atr_orig.keys() if i in atr.keys()]
    if any(atr_orig[i] != atr[i] for i in keyList):
        print('merge %s, %s and %s into %s' % (os.path.basename(m_par_file), os.path.basename(s_par_file),
                                               os.path.basename(off_file), os.path.basename(rsc_file)))
        writefile.write_roipac_rsc(atr, out_file=rsc_file)

    return rsc_file


def extract_metadata4lookup_table(fname):
    """Read/extract attribute for .UTM_TO_RDC file from Gamma to ROI_PAC
    For example, it read input file, sim_150911-150922.UTM_TO_RDC, 
    find its associated par file, sim_150911-150922.utm.dem.par, read it, and
    convert to ROI_PAC style and write it to an rsc file, sim_150911-150922.UTM_TO_RDC.rsc"""

    # Check existed .rsc file
    rsc_file_list = ut.get_file_list(fname+'.rsc')
    if rsc_file_list:
        rsc_file = rsc_file_list[0]
        #print(rsc_file+' is existed, no need to re-extract.')
        return rsc_file

    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]
    atr['Y_UNIT'] = 'degrees'
    atr['X_UNIT'] = 'degrees'

    par_file = os.path.splitext(fname)[0]+'.utm.dem.par'

    #print('read '+os.path.basename(par_file))
    par_dict = readfile.read_gamma_par(par_file)
    atr.update(par_dict)

    # Write to .rsc file
    rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = None
    keyList = [i for i in atr_orig.keys() if i in atr.keys()]
    if any(atr_orig[i] != atr[i] for i in keyList):
        print('writing >>> '+os.path.basename(rsc_file))
        writefile.write_roipac_rsc(atr, out_file=rsc_file)
    return rsc_file


def extract_metadata4dem_geo(fname):
    """Read/extract attribute for .dem file from Gamma to ROI_PAC
    For example, it read input file, sim_150911-150922.utm.dem, 
    find its associated par file, sim_150911-150922.utm.dem.par, read it, and
    convert to ROI_PAC style and write it to an rsc file, sim_150911-150922.utm.dem.rsc
    """
    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]
    atr['Y_UNIT'] = 'degrees'
    atr['X_UNIT'] = 'degrees'

    par_file = fname+'.par'
    #print('read '+os.path.basename(par_file))
    par_dict = readfile.read_gamma_par(par_file)
    atr.update(par_dict)

    # Write to .rsc file
    rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = None
    keyList = [i for i in atr_orig.keys() if i in atr.keys()]
    if any(atr_orig[i] != atr[i] for i in keyList):
        print('writing >>> '+os.path.basename(rsc_file))
        writefile.write_roipac_rsc(atr, out_file=rsc_file)
    return rsc_file


def extract_metadata4dem_radar(fname):
    """Read/extract attribute for .hgt_sim file from Gamma to ROI_PAC
    Input:
        sim_150911-150922.hgt_sim
        sim_150911-150922.rdc.dem
    Search for:
        sim_150911-150922.diff_par
    Output:
        sim_150911-150922.hgt_sim.rsc
        sim_150911-150922.rdc.dem.rsc
    """
    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]

    # Get basename of file
    fname_base = os.path.splitext(fname)[0]
    for i in range(5):
        fname_base = os.path.splitext(fname_base)[0]

    par_file = fname_base+'.diff_par'
    #print('read '+os.path.basename(par_file))
    par_dict = readfile.read_gamma_par(par_file)
    atr.update(par_dict)

    # Write to .rsc file
    rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = None
    keyList = [i for i in atr_orig.keys() if i in atr.keys()]
    if any(atr_orig[i] != atr[i] for i in keyList):
        print('writing >>> '+os.path.basename(rsc_file))
        writefile.write_roipac_rsc(atr, out_file=rsc_file)
    return rsc_file


def prepare_metadata(inps):
    inps.file = ut.get_file_list(inps.file, abspath=True)

    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file), print_msg=False)

    # multiple datasets files
    ext = os.path.splitext(inps.file[0])[1]
    if ext in ['.unw', '.cor', '.int']:
        if len(inps.file) == 1:
            extract_metadata4interferogram(inps.file[0])
        elif inps.parallel:
            Parallel(n_jobs=num_cores)(delayed(extract_metadata4interferogram)(file) for file in inps.file)
        else:
            for File in inps.file:
                extract_metadata4interferogram(File)

    # Single dataset files
    elif inps.file[0].endswith('.utm.dem'):
        for File in inps.file:
            atr_file = extract_metadata4dem_geo(File)
    elif inps.file[0].endswith(('.rdc.dem', '.hgt_sim')):
        for File in inps.file:
            atr_file = extract_metadata4dem_radar(File)
    elif ext in ['.UTM_TO_RDC']:
        for File in inps.file:
            atr_file = extract_metadata4lookup_table(File)
    return


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    prepare_metadata(inps)
    return


###################################################################################################
if __name__ == '__main__':
    main()
