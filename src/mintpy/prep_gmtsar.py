############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Xiaohua Xu, Mar 2021               #
############################################################


import glob
import os

import numpy as np

try:
    from osgeo import gdal
except ImportError:
    raise ImportError('Can not import gdal!')

from mintpy.utils import ptime, readfile, utils as ut, writefile


#########################################################################
def get_prm_files(ifg_dir):
    """Get the date1/2.prm files in the interferogram directory."""
    prm_files = sorted(glob.glob(os.path.join(ifg_dir, '*.PRM')))
    prm_files = [i for i in prm_files if os.path.splitext(os.path.basename(i))[0].isdigit()]
    if len(prm_files) != 2:
        raise FileExistsError('NO two *.prm files found!')
    return prm_files


def get_multilook_number(ifg_dir, fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get the number of multilook in range and azimuth direction."""
    # grab an arbitrary file in radar-coordiantes
    rdr_files = [os.path.join(ifg_dir, f'{i}.grd') for i in fbases]
    if len(rdr_files) == 0:
        raise ValueError(f'No radar-coord files found in {ifg_dir} with suffix: {fbases}')

    # read step info into multilook number
    ds = gdal.Open(rdr_files[0], gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()
    rlooks = int(abs(transform[1]))
    alooks = int(abs(transform[5]))
    return alooks, rlooks


def get_lalo_ref(ifg_dir, prm_dict, fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get the LAT/LON_REF1/2/3/4 from *_ll.grd file."""
    # grab an arbitrary file in geo-coordiantes
    geo_files = [os.path.join(ifg_dir, f'{i}_ll.grd') for i in fbases]
    if len(geo_files) == 0:
        print(f'WARNING: No geo-coord files found in {ifg_dir} with suffix: {fbases}')
        return prm_dict

    # read corners lat/lon info
    ds = gdal.Open(geo_files[0], gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()
    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.
    if not (1e-7 < x_step < 1.):
        raise ValueError('File {} is NOT geocoded!')

    W = transform[0]
    N = transform[3]
    E = W + x_step * ds.RasterXSize
    S = N + y_step * ds.RasterYSize

    if prm_dict['ORBIT_DIRECTION'].upper().startswith('ASC'):
        prm_dict['LAT_REF1'] = str(S)
        prm_dict['LAT_REF2'] = str(S)
        prm_dict['LAT_REF3'] = str(N)
        prm_dict['LAT_REF4'] = str(N)
        prm_dict['LON_REF1'] = str(W)
        prm_dict['LON_REF2'] = str(E)
        prm_dict['LON_REF3'] = str(W)
        prm_dict['LON_REF4'] = str(E)
    else:
        prm_dict['LAT_REF1'] = str(N)
        prm_dict['LAT_REF2'] = str(N)
        prm_dict['LAT_REF3'] = str(S)
        prm_dict['LAT_REF4'] = str(S)
        prm_dict['LON_REF1'] = str(E)
        prm_dict['LON_REF2'] = str(W)
        prm_dict['LON_REF3'] = str(E)
        prm_dict['LON_REF4'] = str(W)

    return prm_dict


def get_slant_range_distance(ifg_dir, prm_dict, fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get a constant slant range distance in the image center, for dataset in geo-coord."""
    # grab an arbitrary file in radar-coordiantes
    rdr_files = [os.path.join(ifg_dir, f'{i}.grd') for i in fbases]
    if len(rdr_files) == 0:
        raise ValueError(f'No radar-coord files found in {ifg_dir} with suffix: {fbases}')

    # read width from rdr_file
    ds = gdal.Open(rdr_files[0], gdal.GA_ReadOnly)
    width = ds.RasterXSize

    near_range = float(prm_dict['STARTING_RANGE'])
    range_pixel_size = float(prm_dict['RANGE_PIXEL_SIZE'])
    slant_range_dist = near_range + range_pixel_size * width / 2.
    prm_dict['SLANT_RANGE_DISTANCE'] = slant_range_dist

    return prm_dict


#########################################################################
def extract_gmtsar_metadata(unw_file, template_file, rsc_file=None, update_mode=True):
    """Extract metadata from GMTSAR interferogram stack."""

    # update_mode: check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=unw_file, readable=False) == 'skip':
        return readfile.read_roipac_rsc(rsc_file)

    ifg_dir = os.path.dirname(unw_file)

    # 1. read *.PRM file
    prm_file = get_prm_files(ifg_dir)[0]
    meta = readfile.read_gmtsar_prm(prm_file)
    meta['PROCESSOR'] = 'gmtsar'

    # 2. read template file: HEADING, ORBIT_DIRECTION
    template = readfile.read_template(template_file)
    for key in ['HEADING', 'ORBIT_DIRECTION']:
        if key in template.keys():
            meta[key] = template[key].lower()
        else:
            raise ValueError('Attribute {} is missing! Please manually specify it in the template file.')

    # 3. grab A/RLOOKS from radar-coord data file
    meta['ALOOKS'], meta['RLOOKS'] = get_multilook_number(ifg_dir)
    meta['AZIMUTH_PIXEL_SIZE'] *= meta['ALOOKS']
    meta['RANGE_PIXEL_SIZE'] *= meta['RLOOKS']

    # 4. grab LAT/LON_REF1/2/3/4 from geo-coord data file
    meta = get_lalo_ref(ifg_dir, meta)

    # 5. grab X/Y_FIRST/STEP from unw_file if in geo-coord
    ds = gdal.Open(unw_file, gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()
    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.
    if 1e-7 < x_step < 1.:
        meta['X_STEP'] = x_step
        meta['Y_STEP'] = y_step
        meta['X_FIRST'] = transform[0] - x_step / 2.
        meta['Y_FIRST'] = transform[3] - y_step / 2.
        # constrain longitude within (-180, 180]
        if meta['X_FIRST'] > 180.:
            meta['X_FIRST'] -= 360.

    # 6. extra metadata for the missing geometry dataset: SLANT_RANGE_DISTANCE / INCIDENCE_ANGLE
    # for dataset in geo-coordinates
    if 'Y_FIRST' in meta.keys():
        meta = get_slant_range_distance(ifg_dir, meta)
        Re = float(meta['EARTH_RADIUS'])
        H = float(meta['HEIGHT'])
        Rg = float(meta['SLANT_RANGE_DISTANCE'])
        Inc = (np.pi - np.arccos((Re**2 + Rg**2 - (Re+H)**2) / (2*Re*Rg))) * 180./np.pi
        meta['INCIDENCE_ANGLE'] = Inc

    # convert all value to string format
    for key, value in meta.items():
        meta[key] = str(value)

    # write to .rsc file
    meta = readfile.standardize_metadata(meta)
    if rsc_file:
        print('writing ', rsc_file)
        os.makedirs(os.path.dirname(rsc_file), exist_ok=True)
        writefile.write_roipac_rsc(meta, rsc_file)

    return meta


def prepare_geometry(geom_files, meta, update_mode=True):
    """Prepare .rsc file for all geometry files."""
    num_file = len(geom_files)
    if num_file == 0:
        raise FileNotFoundError('NO geometry file found!')

    # write .rsc file for each geometry file
    for geom_file in geom_files:
        # copy over the common metadata
        geom_meta = {}
        for key, value in meta.items():
            geom_meta[key] = value

        # update from .grd file
        geom_meta.update(readfile.read_gdal_vrt(geom_file))

        # write .rsc file
        rsc_file = geom_file+'.rsc'
        writefile.write_roipac_rsc(
            geom_meta,
            rsc_file,
            update_mode=update_mode,
            print_msg=True,
        )

    return


def prepare_stack(unw_files, meta, update_mode=True):
    """Prepare .rsc file for all unwrapped interferogram files."""
    num_file = len(unw_files)
    if num_file == 0:
        raise FileNotFoundError('NO unwrapped interferogram file found!')

    # write .rsc file for each interferogram file
    prog_bar = ptime.progressBar(maxValue=num_file)
    for i, unw_file in enumerate(unw_files):
        ifg_dir = os.path.dirname(unw_file)
        ifg_meta = {}

        # copy over the common metadata
        for key, value in meta.items():
            ifg_meta[key] = value

        # update from .grd file
        ifg_meta.update(readfile.read_gdal_vrt(unw_file))

        # add DATE12
        prm_files = get_prm_files(ifg_dir)
        date1, date2 = (os.path.splitext(os.path.basename(i))[0] for i in prm_files)
        ifg_meta['DATE12'] = f'{ptime.yymmdd(date1)}-{ptime.yymmdd(date2)}'

        # and P_BASELINE_TOP/BOTTOM_HDR
        baseline_file = os.path.join(ifg_dir, 'baseline.txt')
        if os.path.isfile(baseline_file):
            bDict = readfile.read_template(baseline_file)
            ifg_meta['P_BASELINE_TOP_HDR'] = bDict['B_perpendicular']
            ifg_meta['P_BASELINE_BOTTOM_HDR'] = bDict['B_perpendicular']
        else:
            ifg_meta['P_BASELINE_TOP_HDR'] = '0'
            ifg_meta['P_BASELINE_BOTTOM_HDR'] = '0'
            msg = f'WARNING: No baseline file found in: {baseline_file}. '
            msg += 'Set P_BASELINE* to 0 and continue.'
            print(msg)

        # write .rsc file
        rsc_file = unw_file+'.rsc'
        writefile.write_roipac_rsc(
            ifg_meta,
            rsc_file,
            update_mode=update_mode,
            print_msg=False,
        )

        prog_bar.update(i+1, suffix=f'{date1}_{date2}')
    prog_bar.close()


#########################################################################
def prep_gmtsar(inps):
    # read file path from template file
    template = readfile.read_template(inps.template_file)
    inps.unw_files = sorted(glob.glob(template['mintpy.load.unwFile']))
    inps.cor_files = sorted(glob.glob(template['mintpy.load.corFile']))
    inps.dem_file = glob.glob(template['mintpy.load.demFile'])[0]

    # extract common metadata
    rsc_file = os.path.join(inps.mintpy_dir, 'inputs/data.rsc')
    meta = extract_gmtsar_metadata(
        unw_file=inps.unw_files[0],
        template_file=inps.template_file,
        rsc_file=rsc_file,
        update_mode=inps.update_mode,
    )

    # prepare metadata for geometry files
    prepare_geometry(
        geom_files=[inps.dem_file],
        meta=meta,
        update_mode=inps.update_mode,
    )

    # prepare metadata for interferogram files
    prepare_stack(
        unw_files=inps.unw_files,
        meta=meta,
        update_mode=inps.update_mode,
    )

    print('Done.')
    return
