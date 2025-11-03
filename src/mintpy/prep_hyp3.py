############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Zhang Yunjun, Mar 2021         #
############################################################


import datetime as dt
import os
import re

from mintpy.constants import SPEED_OF_LIGHT
from mintpy.objects import sensor
from mintpy.utils import readfile
from mintpy.utils import utils1 as ut
from mintpy.utils import writefile

#########################################################################

def _get_product_name_and_type(filename: str) -> tuple[str, str]:
    if match := re.match(
        r'S1_\d{6}_IW[123](_\d{8}){2}_(VV|HH)_INT\d{2}_[0-9A-F]{4}',
        filename,
    ):
        job_type = 'INSAR_ISCE_BURST'

    elif match := re.match(
        r'S1_\d{3}_\d{6}s1n\d{2}-\d{6}s2n\d{2}-\d{6}s3n\d{2}_IW(_\d{8}){2}_(VV|HH)_INT\d{2}_[0-9A-F]{4}',
        filename,
    ):
        job_type = 'INSAR_ISCE_MULTI_BURST'

    elif match := re.match(
        r'S1[ABC]{2}(_\d{8}T\d{6}){2}_(VV|HH)[PRO]\d{3}_INT\d{2}_G_[uw][ec][123F]_[0-9A-F]{4}',
        filename,
    ):
        job_type = 'INSAR_GAMMA'

    else:
        raise ValueError(f'Failed to parse product name from filename: {filename}')

    return match.group(), job_type


def add_hyp3_metadata(fname, meta, is_ifg=True):
    """Read/extract metadata from HyP3 metadata file and add to metadata dictionary.

    Three types of ASF HyP3 products are supported:

    1. INSAR_ISCE_BURST (legacy single-burst product using ISCE2) metadata file:
        format: https://hyp3-docs.asf.alaska.edu/guides/burst_insar_product_guide/#naming-convention-insar_isce_burst
        example: S1_213524_IW1_20170411_20170517_VV_INT80_8E81.txt

    2. INSAR_ISCE_MULTI_BURST (multi-burst product using ISCE2) metadata file:
        format: https://hyp3-docs.asf.alaska.edu/guides/burst_insar_product_guide/#naming-convention-insar_isce_multi_burst
        example: S1_064_000000s1n00-136231s2n02-000000s3n00_IW_20200604_20200616_VV_INT80_77F1

    3. INSAR_GAMMA (scene-wide product using GAMMA) metadata file:
        format: https://hyp3-docs.asf.alaska.edu/guides/insar_product_guide/#naming-convention
        example: S1AA_20190610T135156_20190622T135157_VVP012_INT80_G_ueF_F8BF.txt

    Parameters: fname  - str, path to the hyp3 data file, e.g. *unw_phase_clip*.tif, *dem_clip*.tif
                meta   - dict, existing metadata
                is_ifg - bool, is the data file interferogram (unw/corr) or geometry (dem/angles)

    Returns:    meta   - dict, return metadata
    """
    product_name, job_type = _get_product_name_and_type(os.path.basename(fname))

    meta_file = os.path.join(os.path.dirname(fname), f'{product_name}.txt')
    hyp3_meta = {}
    with open(meta_file) as f:
        for line in f:
            key, value = line.strip().replace(' ','').split(':')[:2]
            hyp3_meta[key] = value

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
    meta['PLATFORM'] = 'Sen'
    meta['ANTENNA_SIDE'] = -1
    meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
    meta['RANGE_PIXEL_SIZE'] = sensor.SEN['range_pixel_size'] * int(meta['RLOOKS'])
    meta['AZIMUTH_PIXEL_SIZE'] = sensor.SEN['azimuth_pixel_size'] * int(meta['ALOOKS'])

    # HyP3 (incidence, azimuth) angle datasets are in the unit of radian,
    # which is different from the isce-2 convention of degree
    if any(x in os.path.basename(fname) for x in ['lv_theta', 'lv_phi']):
        meta['UNIT'] = 'radian'

    # HDF-EOS5 metadata, including:
    # beam_mode/swath, relative_orbit, first/last_frame, unwrap_method

    meta['beam_mode'] = 'IW'
    meta['unwrap_method'] = hyp3_meta['Unwrappingtype']

    if job_type == 'INSAR_ISCE_BURST':
        date1, date2 = (dt.datetime.strptime(x,'%Y%m%d') for x in product_name.split('_')[3:5])
        meta['beam_swath'] = product_name.split('_')[2][2]

        # relative_orbit [to be added]
        # first/last_frame [to be added]

    elif job_type == 'INSAR_ISCE_MULTI_BURST':
        date1, date2 = (dt.datetime.strptime(x, '%Y%m%d') for x in product_name.split('_')[4:6])
        swath_tokens = product_name.split('_')[2].split('-')
        meta['beam_swath'] = ''.join(s[7] for s in swath_tokens if not s.startswith('000000s'))

    else:
        assert job_type == 'INSAR_GAMMA'

        date1, date2 = (dt.datetime.strptime(x,'%Y%m%dT%H%M%S') for x in product_name.split('_')[1:3])
        meta['beam_swath'] = '123'

        ref_granule = hyp3_meta['ReferenceGranule']
        assert ref_granule.startswith('S1')

        abs_orbit = int(hyp3_meta['ReferenceOrbitNumber'])
        if ref_granule.startswith('S1A'):
            meta['relative_orbit'] = ((abs_orbit - 73) % 175) + 1
        elif ref_granule.startswith('S1B'):
            meta['relative_orbit'] = ((abs_orbit - 202) % 175) + 1
        elif ref_granule.startswith('S1C'):
            meta['relative_orbit'] = ((abs_orbit - 172) % 175) + 1
        else:
            # add equation for Sentinel-D in the future
            raise ValueError(f'Un-recognized Sentinel-1 satellite from {ref_granule}!')

        # first/last_frame [to be completed]
        t0, t1 = ref_granule.split('_')[-5:-3]
        meta['startUTC'] = dt.datetime.strptime(t0, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d %H:%M:%S.%f')
        meta['stopUTC']  = dt.datetime.strptime(t1, '%Y%m%dT%H%M%S').strftime('%Y-%m-%d %H:%M:%S.%f')
        # ascendingNodeTime [to be added]

    # interferogram related metadata
    if is_ifg:
        meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'
        meta['P_BASELINE_TOP_HDR'] = hyp3_meta['Baseline']
        meta['P_BASELINE_BOTTOM_HDR'] = hyp3_meta['Baseline']

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
