############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Yunmeng Cao, 2017                  #
############################################################


import os
import re

import numpy as np

from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, utils as ut, writefile

SPEED_OF_LIGHT = 299792458  # m/s
# list of par file extension for SAR images
PAR_EXT_LIST = ['.amp.par', '.ramp.par', '.mli.par']


######################################## Sub Functions ############################################
def get_perp_baseline(m_par_file, s_par_file, off_file, atr_dict={}):
    """Get perpendicular baseline info from reference/secondary par file and off file.
    Parameters: m_par_file - str, path, reference parameter file, i.e. 130118_4rlks.amp.par
                s_par_file - str, path, secondary parameter file, i.e. 130129_4rlks.amp.oar
                off_file   - str, path, interferogram off file, i.e. 130118-130129_4rlks.off
                atr_dict   - dict, optional, attributes dictionary
    Returns:    bperp      - str, perpendicular baseline for pixel at [0,0]
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

    f = open(base_perp_file)
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
    If it does not exist, call Gamma script - SLC_corners - to generate it from SLC par file
        e.g. 130118_4rlks.amp.par

    Parameters: m_par_file - str, path, reference date parameter file, i.e. 130118_4rlks.amp.par
                atr_dict   - dict, optional, attributes dictionary
    Returns:    atr_dict   - dict, attributes dictionary, contains LAT/LON_REF1/2/3/4
    """
    m_par_file = os.path.abspath(m_par_file)
    m_corner_file = os.path.splitext(m_par_file)[0]+'.corner'
    m_corner_full_file = m_corner_file+'_full'

    # Call Gamma command to calculate LAT/LON_REF
    if not os.path.isfile(m_corner_file):
        # generate *.corner_full file
        if not os.path.isfile(m_corner_full_file):
            cornerCmd = 'SLC_corners {} > {}'.format(m_par_file,
                                                     m_corner_full_file)
            print(cornerCmd)
            os.system(cornerCmd)

        # convert *.corner_full to *.corner
        extractCmd = "awk 'NR==3,NR==6 {print $3,$6} ' "
        extractCmd += f"{m_corner_full_file} > {m_corner_file}"
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


def extract_metadata4interferogram(fname, sensor_name=None):
    """Read/extract attributes from Gamma .unw, .cor and .int file
    Parameters: fname - str, Gamma interferogram filename or path,
                        i.e. /PopoSLT143TsxD/diff_filt_HDR_130118-130129_4rlks.unw
    Returns:    atr   - dict, Attributes dictionary
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
        date12 = str(re.findall(r'\d{8}[-_]\d{8}', file_basename)[0])
    except:
        date12 = str(re.findall(r'\d{6}[-_]\d{6}', file_basename)[0])
    m_date, s_date = date12.replace('-', '_').split('_')
    atr['DATE12'] = ptime.yymmdd(m_date)+'-'+ptime.yymmdd(s_date)
    lks = os.path.splitext(file_basename)[0].split(date12)[1]
    #lks = os.path.splitext(file_basename.split(date12)[1])[0]

    # Read .off and .par file
    off_files = file_dir+'/*'+date12+lks+'.off'
    m_par_files = [file_dir+'/*'+m_date+lks+i for i in PAR_EXT_LIST]
    s_par_files = [file_dir+'/*'+s_date+lks+i for i in PAR_EXT_LIST]

    try:
        m_par_file = ut.get_file_list(m_par_files)[0]
    except FileNotFoundError:
        m_par_file = None
        print('\nERROR: Can not find reference date .par file, it supposed to be like: '+m_par_files)
    try:
        s_par_file = ut.get_file_list(s_par_files)[0]
    except FileNotFoundError:
        s_par_file = None
        print('\nERROR: Can not find secondary date .par file, it supposed to be like: '+s_par_files)

    try:
        off_file = ut.get_file_list(off_files)[0]
    except FileNotFoundError:
        off_file = file_dir+'/'+date12+lks+'.off'
        offCmd = f'create_offset {m_par_file} {s_par_file} {off_file} 1 1 1 0'
        print(offCmd)
        os.system(offCmd)

    par_dict = readfile.read_gamma_par(m_par_file)
    off_dict = readfile.read_gamma_par(off_file)
    atr.update(par_dict)
    atr.update(off_dict)

    # Perp Baseline Info
    atr = get_perp_baseline(m_par_file, s_par_file, off_file, atr)

    # LAT/LON_REF1/2/3/4
    atr = get_lalo_ref(m_par_file, atr)

    # NCORRLOOKS
    if sensor_name:
        rg_bandwidth = float(atr['chirp_bandwidth'])
        rg_resolution = SPEED_OF_LIGHT / (2. * rg_bandwidth)
        rg_pixel_size = float(atr['RANGE_PIXEL_SIZE']) / float(atr['RLOOKS'])
        rg_fact = rg_resolution / rg_pixel_size

        antenna_length = sensor.SENSOR_DICT[sensor_name.lower()]['antenna_length']
        az_resolution = antenna_length / 2
        az_pixel_size = float(atr['AZIMUTH_PIXEL_SIZE']) / float(atr['ALOOKS'])
        az_fact = az_resolution / az_pixel_size

        ncorr_looks = float(atr['RLOOKS']) * float(atr['ALOOKS']) / (rg_fact * az_fact)
        atr['NCORRLOOKS'] = ncorr_looks

    # Write to .rsc file
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = dict()

    if not set(atr.items()).issubset(set(atr_orig.items())):
        print('merge {}, {} and {} into {}'.format(
            os.path.basename(m_par_file),
            os.path.basename(s_par_file),
            os.path.basename(off_file),
            os.path.basename(rsc_file),
        ))
        atr_out = {**atr_orig, **atr}
        writefile.write_roipac_rsc(atr_out, out_file=rsc_file)

    return rsc_file


def extract_metadata4geometry_radar(fname):
    """Read/extract attribute for .hgt_sim file from Gamma to ROI_PAC
    Parameters: fname    - geometry file in radar coordinates, e.g.
                           sim_20070813_20080310.hgt_sim
                           sim_20070813_20080310.rdc.dem
                           which leads to file:
                           sim_20070813_20080310.diff_par
    Returns:    rsc_file - str, metadata file
                           sim_20070813_20080310.hgt_sim.rsc
                           sim_20070813_20080310.rdc.dem.rsc
    """
    # Get/read GAMMA par file
    # for loop to get rid of multiple dot in filename
    fname_base = os.path.splitext(fname)[0]
    for i in range(5):
        fname_base = os.path.splitext(fname_base)[0]
    par_file = fname_base+'.diff_par'
    par_dict = readfile.read_gamma_par(par_file)

    # Get/read LAT/LON_REF1/2/3/4
    msg = 'grab LAT/LON_REF1/2/3/4 from par file: '
    # get date of one acquisition
    try:
        m_date = str(re.findall(r'\d{8}', fname_base)[0])
    except:
        m_date = str(re.findall(r'\d{6}', fname_base)[0])

    ## search existing par file
    # potential directories
    geom_dir = os.path.dirname(fname)
    proj_dir = os.path.dirname(geom_dir)
    par_dirs = [geom_dir,                                              #PROJECT_DIR/geom_reference
                os.path.join(proj_dir, f'*/{m_date[2:]}_*'),   #PROJECT_DIR/interferograms/{m_date}_(20)141225
                os.path.join(proj_dir, f'*/{m_date[2:]}-*')]   #PROJECT_DIR/interferograms/{m_date}-(20)141225
    # potential file patterns
    m_par_files = []
    for par_dir in par_dirs:
        m_par_files += [os.path.join(par_dir, f'*{m_date}*{ext}') for ext in PAR_EXT_LIST]
    # search existing files that meet the file patterns
    m_par_files = ut.get_file_list(m_par_files)

    # read par file
    if len(m_par_files) > 0:
        m_par_file = m_par_files[0]
        msg += m_par_file
        par_dict = get_lalo_ref(m_par_file, par_dict)
    else:
        msg += f' no par file found with date: {m_date}'
    print(msg)

    # initiate ROIPAC dict
    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = os.path.splitext(fname)[1]
    atr.update(par_dict)

    # Write to .rsc file
    rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = dict()

    if not set(atr.items()).issubset(set(atr_orig.items())):
        print('writing >>> '+os.path.basename(rsc_file))
        atr_out = {**atr_orig, **atr}
        writefile.write_roipac_rsc(atr_out, out_file=rsc_file)

    return rsc_file


def extract_metadata4geometry_geo(fname):
    """Read/extract attribute for *.dem / *.UTM_TO_RDC file from Gamma to ROI_PAC
    Parameters: fname    - str, geometry file in geo-coordinates, e.g.
                           sim_20070813_20080310.utm.dem
                           sim_20070813_20080310.UTM_TO_RDC
                           which leads to file:
                           sim_20070813_20080310.utm.dem.par
    Returns:    rsc_file - str, metadata file
                           sim_20070813_20080310.utm.dem.rsc
                           sim_20070813_20080310.UTM_TO_RDC.rsc
    """
    # Get/read GAMMA par file
    ext = os.path.splitext(fname)[1]
    if ext.lower().endswith(('to_rdc', '2_rdc', '2rdc')):
        par_file = os.path.splitext(fname)[0]+'.utm.dem.par'
    elif fname[0].endswith('.utm.dem'):
        par_file = fname+'.par'
    par_dict = readfile.read_gamma_par(par_file)

    # initiate ROIPAC dict
    atr = {}
    atr['PROCESSOR'] = 'gamma'
    atr['FILE_TYPE'] = ext

    if 'post_lat' in par_dict.keys():
        # coordinates in degree
        atr['Y_UNIT'] = 'degrees'
        atr['X_UNIT'] = 'degrees'
    elif 'post_north' in par_dict.keys():
        # coordinates in meter
        atr['Y_UNIT'] = 'm'
        atr['X_UNIT'] = 'm'
        atr['X_FIRST'] = float(par_dict['corner_east']) - int(par_dict['width']) * float(par_dict['post_east'])
    else:
        msg = 'un-recognized coordinates type:'
        for key, value in par_dict.items():
            if key.startswith(('corner','post')):
                msg += f'\n\t{key}: {value}'
        raise ValueError(msg)

    atr.update(par_dict)

    # Write to .rsc file
    rsc_file = fname+'.rsc'
    try:
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    except:
        atr_orig = dict()

    if not set(atr.items()).issubset(set(atr_orig.items())):
        print('writing >>> '+os.path.basename(rsc_file))
        atr_out = {**atr_orig, **atr}
        writefile.write_roipac_rsc(atr_out, out_file=rsc_file)

    return rsc_file


###################################################################################################
def prep_gamma(inps):

    # loop for each file
    for fname in inps.file:
        # interferograms
        if inps.file_ext in ['.unw', '.cor', '.int']:
            extract_metadata4interferogram(fname, sensor_name=inps.sensor)

        # geometry - geo
        elif inps.file_ext.endswith(('to_rdc', '2_rdc', '2rdc')) or fname.endswith('.utm.dem'):
            extract_metadata4geometry_geo(fname)

        # geometry - radar
        elif fname.endswith(('.rdc.dem', '.hgt_sim')):
            extract_metadata4geometry_radar(fname)

    return
