############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Zhang Yunjun, Mar 2021         #
############################################################


import datetime as dt
import os

from mintpy.constants import SPEED_OF_LIGHT
from mintpy.objects import sensor
from mintpy.utils import readfile, utils1 as ut, writefile


#########################################################################
def add_hyp3_metadata(fname, meta, is_ifg=True):
    '''Read/extract metadata from HyP3 metadata file and add to metadata dictionary.

    Two types of ASF HyP3 products are supported: isce2_burst, gamma_scene
    1. isce2_burst (burst-wide product using ISCE2) metadata file:
        format: {SAT}_{FRAME}_{SUBSWATH}_{DATE1}_{DATE2}_{POL}_{RES}_{IDS}.txt
        example: S1_213524_IW1_20170411_20170517_VV_INT80_8E81.txt
        content:
            Reference Granule: S1_213524_IW1_20170411T133605_VV_BD30-BURST
            ...

    2. isce2_multi_burst (multi-burst product using ISCE2):
        - Metadata filename format:
            {SAT}_{PATH}_{LON_L}_{LON_U}_{LAT_L}_{LAT_U}_{DATE1}_{DATE2}_{POL}_{RES}_{JOBID}
        - Example:
            S1A_064_E053_1_N27_3_E054_1_N27_8_20200604_20200616_VV_INT80_7EB5.txt
        - Content:
            Reference Granule: S1_290876_IW1_20200604T005150_VV_XXXX-BURST, S1_290875_IW1_20200604T005147_VV_XXXX-BURST
            Secondary Granule: S1_290876_IW1_20200616T005150_VV_YYYY-BURST, S1_290875_IW1_20200616T005147_VV_YYYY-BURST
            ...

    3. gamma_scene (scene-wide product using Gamma) metadata file:
        format: {SAT}_{DATE1}_{DATE2}_{POL}_{RES}_{SOFT}_{PROC}_{IDS}.txt
        example: S1AA_20190610T135156_20190622T135157_VVP012_INT80_G_ueF_F8BF.txt
        content:
            Reference Granule: S1A_IW_SLC__1SDV_20190704T135158_20190704T135225_027968_032877_1C4D
            ...

    Parameters: fname  - str, path to the hyp3 data file, e.g. *unw_phase_clip*.tif, *dem_clip*.tif
                meta   - dict, existing metadata
                is_ifg - bool, is the data file interferogram (unw/corr) or geometry (dem/angles)
    Returns:    meta   - dict, return metadata
    '''

    # job_id -> prod_type and date1/2 objects
    job_id = '_'.join(os.path.splitext(os.path.basename(os.path.dirname(fname)))[0].split('_'))
    parts = job_id.split('_')

    if len(parts) < 8:
        raise ValueError(f"Unexpected format: {job_id}")

    if parts[2].startswith('IW'):
        # burst-wide product using ISCE2
        prod_type = 'isce2_burst'
        date1 = dt.datetime.strptime(parts[3], '%Y%m%d')
        date2 = dt.datetime.strptime(parts[4], '%Y%m%d')    
    
    elif len(parts) >= 12 and parts[10].isdigit() and parts[11].isdigit():
        prod_type = 'isce2_multi_burst'
        date1 = dt.datetime.strptime(parts[10], '%Y%m%d')
        date2 = dt.datetime.strptime(parts[11], '%Y%m%d')

    else:
        # scene-wide product using Gamma
        prod_type = 'gamma_scene'
        date1 = dt.datetime.strptime(parts[1], '%Y%m%dT%H%M%S')
        date2 = dt.datetime.strptime(parts[2], '%Y%m%dT%H%M%S')
    # read hyp3 metadata file
    meta_file = os.path.join(os.path.dirname(fname), f'{job_id}.txt')
    hyp3_meta = {}
    with open(meta_file) as f:
        for line in f:
            key, value = line.strip().replace(' ','').split(':')[:2]
            hyp3_meta[key] = value
    ref_granule = hyp3_meta['ReferenceGranule'].split(',')[0].strip() 

    # add universal hyp3 metadata
    meta['PROCESSOR'] = 'hyp3'
    meta['CENTER_LINE_UTC'] = hyp3_meta['UTCtime']
    meta['ALOOKS'] = hyp3_meta['Azimuthlooks']
    meta['RLOOKS'] = hyp3_meta['Rangelooks']
    meta['EARTH_RADIUS'] = hyp3_meta['Earthradiusatnadir']
    meta['HEIGHT'] = hyp3_meta['Spacecraftheight']
    meta['STARTING_RANGE'] = hyp3_meta['Slantrangenear']
    meta['HEADING'] = float(hyp3_meta['Heading']) % 360. - 360.  # ensure negative value

    # add LAT/LON_REF1/2/3/4 based on whether satellite ascending or descending
    meta['ORBIT_DIRECTION'] = 'ASCENDING' if abs(meta['HEADING']) < 90 else 'DESCENDING'
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

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

    # hard-coded metadata for Sentinel-1
    if ref_granule.startswith('S1'):
        meta['PLATFORM'] = 'Sen'
        meta['ANTENNA_SIDE'] = -1
        meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
        meta['RANGE_PIXEL_SIZE'] = sensor.SEN['range_pixel_size'] * int(meta['RLOOKS'])
        meta['AZIMUTH_PIXEL_SIZE'] = sensor.SEN['azimuth_pixel_size'] * int(meta['ALOOKS'])

    # HyP3 (incidence, azimuth) angle datasets are in the unit of radian,
    # which is different from the isce-2 convention of degree
    if any(x in os.path.basename(fname) for x in ['lv_theta', 'lv_phi']):
        meta['UNIT'] = 'radian'

    # interferogram related metadata
    if is_ifg:
        meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'
        meta['P_BASELINE_TOP_HDR'] = hyp3_meta['Baseline']
        meta['P_BASELINE_BOTTOM_HDR'] = hyp3_meta['Baseline']

    # [optional] HDF-EOS5 metadata, including:
    # beam_mode/swath, relative_orbit, first/last_frame, unwrap_method
    if ref_granule.startswith('S1'):
        # beam_mode
        meta['beam_mode'] = 'IW'

        if prod_type == 'isce2_burst' or prod_type == 'isce2_multi_burst':
            # burst-wide product using ISCE2
            meta['beam_swath'] = job_id.split('_')[2][2:]

            # relative_orbit -> https://forum.step.esa.int/t/sentinel-1-relative-orbit-from-filename/7042
            abs_orbit = int(hyp3_meta['ReferenceOrbitNumber'])
            if job_id.startswith('S1A'):
                meta['relative_orbit'] = ((abs_orbit - 73) % 175) + 1
            elif job_id.startswith('S1B'):
                meta['relative_orbit'] = ((abs_orbit - 202) % 175) + 1

            # first/last_frame [to be added]

        else:
            # scene-wide product using Gamma
            meta['beam_swath'] = '123'

            # relative_orbit
            abs_orbit = int(hyp3_meta['ReferenceOrbitNumber'])
            if ref_granule.startswith('S1A'):
                meta['relative_orbit'] = ((abs_orbit - 73) % 175) + 1
            elif ref_granule.startswith('S1B'):
                meta['relative_orbit'] = ((abs_orbit - 202) % 175) + 1
            else:
                # add equation for Sentinel-C/D in the future
                raise ValueError(f'Un-recognized Sentinel-1 satellite from {ref_granule}!')

            # first/last_frame [to be completed]
            t0, t1 = ref_granule.split('_')[-5:-3]
            meta['startUTC'] = dt.datetime.strptime(t0, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d %H:%M:%S.%f')
            meta['stopUTC']  = dt.datetime.strptime(t1, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d %H:%M:%S.%f')
            # ascendingNodeTime [to be added]

    # unwrap_method
    meta['unwrap_method'] = hyp3_meta['Unwrappingtype']

    return meta


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
