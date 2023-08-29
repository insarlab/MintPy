############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################


import datetime as dt
import os

from mintpy.constants import SPEED_OF_LIGHT
from mintpy.objects import sensor
from mintpy.utils import readfile, utils1 as ut, writefile


#########################################################################
def add_hyp3_metadata(fname, meta, is_ifg=True):
    '''Read/extract attribute data from HyP3 metadata file and add to metadata dictionary
    Inputs:
        *unw_phase.tif, *corr.tif file name, *dem.tif, *inc_map.tif, e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    # read hyp3 metadata file
    # e.g.: burst-wide product using ISCE2: {SAT}_{FRAME}_{SUBSWATH}_{DATE1}_{DATE2}_{POL}_{RES}_{IDS}.txt
    #       scene-wide product using Gamma: {SAT}_{DATE1}_{DATE2}_{POL}_{RES}_{SOFT}_{PROC}_{IDS}.txt
    job_id = '_'.join(os.path.basename(fname).split('_')[:8])
    meta_file = os.path.join(os.path.dirname(fname), f'{job_id}.txt')
    hyp3_meta = {}
    with open(meta_file) as f:
        for line in f:
            key, value = line.strip().replace(' ','').split(':')[:2]
            hyp3_meta[key] = value

    # get date1/2 objects
    if job_id.split('_')[2].startswith('IW'):
        # burst-wide product using ISCE2
        date1_str, date2_str = job_id.split('_')[3:5]
        date1 = dt.datetime.strptime(f'{date1_str}','%Y%m%d')
        date2 = dt.datetime.strptime(f'{date2_str}','%Y%m%d')
    else:
        # scene-wide product using Gamma
        date1_str, date2_str = job_id.split('_')[1:3]
        date1 = dt.datetime.strptime(date1_str,'%Y%m%dT%H%M%S')
        date2 = dt.datetime.strptime(date2_str,'%Y%m%dT%H%M%S')

    # add universal hyp3 metadata
    meta['PROCESSOR'] = 'hyp3'
    meta['CENTER_LINE_UTC'] = hyp3_meta['UTCtime']
    meta['ALOOKS'] = hyp3_meta['Azimuthlooks']
    meta['RLOOKS'] = hyp3_meta['Rangelooks']
    meta['EARTH_RADIUS'] = hyp3_meta['Earthradiusatnadir']
    meta['HEIGHT'] = hyp3_meta['Spacecraftheight']
    meta['STARTING_RANGE'] = hyp3_meta['Slantrangenear']
    # ensure negative value for the heading angle
    meta['HEADING'] = float(hyp3_meta['Heading']) % 360. - 360.

    # add LAT/LON_REF1/2/3/4 based on whether satellite ascending or descending
    meta['ORBIT_DIRECTION'] = 'ASCENDING' if abs(meta['HEADING']) < 90 else 'DESCENDING'
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

    # convert UTM to lat/lon
    N, W = ut.utm2latlon(meta, W, N)
    S, E = ut.utm2latlon(meta, E, S)

    if meta['ORBIT_DIRECTION'] == 'ASCENDING':
        meta['LAT_REF1'] = str(S)
        meta['LAT_REF2'] = str(S)
        meta['LAT_REF3'] = str(N)
        meta['LAT_REF4'] = str(N)
        meta['LON_REF1'] = str(W)
        meta['LON_REF2'] = str(E)
        meta['LON_REF3'] = str(W)
        meta['LON_REF4'] = str(E)
    else:
        meta['LAT_REF1'] = str(N)
        meta['LAT_REF2'] = str(N)
        meta['LAT_REF3'] = str(S)
        meta['LAT_REF4'] = str(S)
        meta['LON_REF1'] = str(E)
        meta['LON_REF2'] = str(W)
        meta['LON_REF3'] = str(E)
        meta['LON_REF4'] = str(W)

    # note: HyP3 currently only supports Sentinel-1 data, so Sentinel-1
    #       configuration is hard-coded.
    if hyp3_meta['ReferenceGranule'].startswith('S1'):
        meta['PLATFORM'] = 'Sen'
        meta['ANTENNA_SIDE'] = -1
        meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
        meta['RANGE_PIXEL_SIZE'] = sensor.SEN['range_pixel_size'] * int(meta['RLOOKS'])
        meta['AZIMUTH_PIXEL_SIZE'] = sensor.SEN['azimuth_pixel_size'] * int(meta['ALOOKS'])

    # note: HyP3 (incidence, azimuth) angle datasets are in the unit of radian
    # which is different from the isce-2 convention of degree
    if any(x in os.path.basename(fname) for x in ['lv_theta', 'lv_phi']):
        meta['UNIT'] = 'radian'

    # add metadata that is only relevant to interferogram files
    if is_ifg:
        meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'
        meta['P_BASELINE_TOP_HDR'] = hyp3_meta['Baseline']
        meta['P_BASELINE_BOTTOM_HDR'] = hyp3_meta['Baseline']

    return(meta)


#########################################################################
def prep_hyp3(inps):
    """Prepare ASF HyP3 metadata files"""

    inps.file = ut.get_file_list(inps.file, abspath=True)

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        is_ifg = any([x in fname for x in ['unw_phase','corr']])
        meta = readfile.read_gdal_vrt(fname)
        meta = add_hyp3_metadata(fname, meta, is_ifg=is_ifg)

        # write
        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return
