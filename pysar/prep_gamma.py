#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import argparse
import re

import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


######################################## Sub Functions ############################################
def attribute_gamma2roipac(par_dict):
    '''Convert Gamma par attribute into ROI_PAC format'''
    key_list = par_dict.keys()

    # Length - number of rows
    key = 'azimuth_lines'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    key = 'interferogram_azimuth_lines'
    if key in key_list:
        par_dict['FILE_LENGTH'] = par_dict[key]

    # Width - number of columns
    key = 'range_samples'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    key = 'interferogram_width'
    if key in key_list:
        par_dict['WIDTH'] = par_dict[key]

    # WAVELENGTH
    speed_of_light = 299792458.0   # meter/second
    key = 'radar_frequency'
    if key in key_list:
        par_dict['WAVELENGTH'] = str(speed_of_light/float(par_dict[key]))

    # HEIGHT & EARTH_RADIUS
    key = 'earth_radius_below_sensor'
    if key in key_list:
        par_dict['EARTH_RADIUS'] = par_dict[key]

        key2 = 'sar_to_earth_center'
        if key2 in key_list:
            par_dict['HEIGHT'] = str(float(par_dict[key2]) - float(par_dict[key]))

    # UTC TIME
    key = 'center_time'
    if key in key_list:
        par_dict['CENTER_LINE_UTC'] = par_dict[key]

    # STARTING_RANGE
    key = 'near_range_slc'
    if key in key_list:
        par_dict['STARTING_RANGE'] = par_dict[key]

    # RANGE_PIXEL_SIZE
    key = 'range_pixel_spacing'
    if key in key_list:
        par_dict['RANGE_PIXEL_SIZE'] = par_dict[key]

    # PLATFORM
    key = 'sensor'
    if key in key_list:
        par_dict['PLATFORM'] = par_dict[key]

    # ORBIT_DIRECTION
    key = 'heading'
    if key in key_list:
        value = float(par_dict[key])
        if 270 < value < 360 or -90 < value < 90:
            par_dict['ORBIT_DIRECTION'] = 'ascending'
        else:
            par_dict['ORBIT_DIRECTION'] = 'descending'

        par_dict['HEADING'] = str(value)

    ##### Optional attributes for PySAR from ROI_PAC
    key = 'azimuth_angle'
    if key in key_list:
        value = float(par_dict[key])
        if 0 < value < 180:
            par_dict['ANTENNA_SIDE'] = '-1'
        else:
            par_dict['ANTENNA_SIDE'] = '1'

    key = 'prf'
    if key in key_list:
        par_dict['PRF'] = par_dict['prf']

    return par_dict


def get_perp_baseline(m_par_file, s_par_file, off_file, atr_dict={}):
    '''Get perpendicular baseline info from master/slave par file and off file.
    Parameters: m_par_file : str, path, master parameter file, i.e. 130118_4rlks.amp.par
                s_par_file : str, path, slave  parameter file, i.e. 130129_4rlks.amp.oar
                off_file   : str, path, interferogram off file, i.e. 130118-130129_4rlks.off
                atr_dict   : dict, optional, attributes dictionary
    Returns:  bperp : str, perpendicular baseline for pixel at [0,0]
    '''
    # Get Path Info
    off_file = os.path.abspath(off_file)
    file_dir = os.path.dirname(off_file)
    file_base = os.path.splitext(os.path.basename(off_file))[0]
    
    base_perp_file = file_dir+'/'+file_base+'.base_perp'
    baseline_file  = file_dir+'/'+file_base+'.baseline'

    # Call Gamma command to calculate Bperp
    if not os.path.isfile(base_perp_file):
        if not os.path.isfile(baseline_file):
            baseCmd = 'base_orbit '+m_par_file+' '+s_par_file+' '+baseline_file
            print baseCmd
            os.system(baseCmd)
        bperpCmd = 'base_perp '+baseline_file+' '+m_par_file+' '+off_file+' > '+base_perp_file
        print bperpCmd
        os.system(bperpCmd)

    # Read bperp txt file
    f = open(base_perp_file,'r')
    line = f.readlines()[12]
    bperp = line.split()[7]
    f.close()

    atr_dict['P_BASELINE_TOP_HDR'] = str(bperp)
    atr_dict['P_BASELINE_BOTTOM_HDR'] = str(bperp)

    return atr_dict


def get_lalo_ref(m_par_file, atr_dict={}):
    '''Get lat/lon of the four corners of the interferograms.
    Parameters: m_par_file : str, path, master date parameter file, i.e. 130118_4rlks.amp.par
                atr_dict   : dict, optional, attributes dictionary
    Returns:  lalo_ref
    '''
    m_par_file = os.path.abspath(m_par_file)
    m_corner_file = os.path.splitext(m_par_file)[0]+'.corner'
    m_corner_full_file = m_corner_file+'_full'

    # Call Gamma command to calculate LAT/LON_REF
    if not os.path.isfile(m_corner_file):
        if not os.path.isfile(m_corner_full_file):
            cornerCmd = 'SLC_corners '+m_par_file+' > '+m_corner_full_file
            print cornerCmd
            os.system(cornerCmd)
        extractCmd = "awk 'NR==3,NR==6 {print $3,$6} ' "+m_corner_full_file+' > '+m_corner_file
        print extractCmd
        os.system(extractCmd)

    # Read corner txt file
    lalo_ref = np.loadtxt(m_corner_file, dtype='str')
    atr_dict['LAT_REF1'] = lalo_ref[0,0]
    atr_dict['LAT_REF2'] = lalo_ref[1,0]
    atr_dict['LAT_REF3'] = lalo_ref[2,0]
    atr_dict['LAT_REF4'] = lalo_ref[3,0]
    atr_dict['LON_REF1'] = lalo_ref[0,1]
    atr_dict['LON_REF2'] = lalo_ref[1,1]
    atr_dict['LON_REF3'] = lalo_ref[2,1]
    atr_dict['LON_REF4'] = lalo_ref[3,1]

    return atr_dict


def extract_attribute(fname):
    '''Read/extract attributes for PySAR from Gamma product
    Parameters: fname : str
                    Gamma interferogram filename or path, i.e. /PopoSLT143TsxD/diff_filt_HDR_130118-130129_4rlks.unw
    Returns:    atr : dict
                    Attributes dictionary
    '''
    atr = {}
    atr['PROCESSOR'] = 'gamma'

    ext = os.path.splitext(fname)[1]
    atr['FILE_TYPE'] = ext

    ## Get info: dir, date12, num of loooks
    file_dir = os.path.dirname(fname)
    file_basename = os.path.basename(fname)
    date12 = str(re.findall('\d{6}[-_]\d{6}', file_basename)[0]).replace('_','-')
    m_date, s_date = date12.split('-')
    lks = os.path.splitext(file_basename.split(date12)[1])[0]
    atr['DATE12'] = date12

    ## Read .off and .par file
    off_file   = file_dir+'/'+date12+lks+'.off'
    m_par_file = file_dir+'/'+m_date+lks+'.amp.par'
    s_par_file = file_dir+'/'+s_date+lks+'.amp.par'

    par_dict = readfile.read_gamma_par(m_par_file)
    off_dict = readfile.read_gamma_par(off_file)

    atr.update(par_dict)
    atr.update(off_dict)
    atr = attribute_gamma2roipac(atr)

    ## Perp Baseline Info
    atr = get_perp_baseline(m_par_file, s_par_file, off_file, atr)

    ## LAT/LON_REF1/2/3/4
    atr = get_lalo_ref(m_par_file, atr)

    return atr



##################################################################################################
EXAMPLE='''example:
  prep_gamma.py  IFGRAM_PopoSLT143TsxD_130118-130129_0011_-0000/diff_filt_HDR_130118-130129_4rlks.unw
  prep_gamma.py  IFGRAM*/diff_filt_*.unw
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Prepare attributes file for Gamma product to be loaded into PySAR.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='Gamma file(s)')
    inps = parser.parse_args()
    return inps


##################################################################################################
def main(argv):
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file, abspath=True)
    print 'number of files: '+str(len(inps.file))

    for File in inps.file:
        atr = extract_attribute(File)
        atr_file = File+'.rsc'
        print 'writing >>> '+atr_file
        writefile.write_roipac_rsc(atr, atr_file)

    print 'Done.'
    return

###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



