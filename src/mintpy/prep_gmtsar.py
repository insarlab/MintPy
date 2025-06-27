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
def get_ifg_file(ifg_dir, coord_type='geo', fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get the interferogram data file path for a given interferogram directory.

    Parameters: ifg_dir    - str, path of a given interferogram directory
                coord_type - str, coordinate type in radar or geo
                fbases     - list(str), data file base name to look for
    Returns:    ifg_file   - str, path of the interferogram data file
    """
    ifg_file = None

    # search interferogram data file
    for fbase in fbases:
        if coord_type == 'geo':
            fnames = glob.glob(os.path.join(ifg_dir, f'{fbase}_ll*.grd'))
        else:
            fnames = glob.glob(os.path.join(ifg_dir, f'{fbase}.grd'))

        if len(fnames) > 0:
            ifg_file = fnames[0]
            break

    if not ifg_file:
        print(f'WARNING: No {coord_type}-coord files found in {ifg_dir} with suffix: {fbases}')

    return ifg_file


def get_prm_files(ifg_dir):
    """Get the date1/2.prm files in the interferogram directory.

    Parameters: ifg_dir   - str, path to the interferogram directory
    Returns:    prm_files - list(str), path to the PRM metadata files
    """
    prm_files = sorted(glob.glob(os.path.join(ifg_dir, '*.PRM')))
    prm_files = [i for i in prm_files if os.path.splitext(os.path.basename(i))[0].isdigit()]
    if len(prm_files) != 2:
        raise FileExistsError('NO two *.prm files found!')
    return prm_files


def get_multilook_number(ifg_dir, fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get the number of multilook in range and azimuth direction.

    Parameters: ifg_dir  - str, path to the interferogram directory
    Returns:    a/rlooks - int, number of looks in the azimuth / range direction
    """
    # grab an arbitrary file in radar-coordiantes
    rdr_file = get_ifg_file(ifg_dir, coord_type='radar')

    # read step info into multilook number
    ds = gdal.Open(rdr_file, gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()
    rlooks = int(abs(transform[1]))
    alooks = int(abs(transform[5]))
    return alooks, rlooks


def get_lalo_ref(ifg_dir, prm_dict, fbases=['corr', 'phase', 'phasefilt', 'unwrap']):
    """Get the LAT/LON_REF1/2/3/4 from *_ll.grd file.

    Parameters: ifg_dir  - str, path to the interferogram directory
                prm_dict - dict, metadata from the PRM file
                fbases   - list(str), data file base name to look for
    Returns:    prm_dict - dict, metadata from the PRM file
    """
    # grab an arbitrary file in geo-coordiantes
    geo_file = get_ifg_file(ifg_dir, coord_type='geo')
    if not geo_file:
        return prm_dict

    # read corners lat/lon info
    print(f'grab LAT/LON_REF1/2/3/4 from file: {geo_file}')
    ds = gdal.Open(geo_file, gdal.GA_ReadOnly)
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
    """Get a constant slant range distance in the image center, for dataset in geo-coord.

    Parameters: ifg_dir  - str, path to the interferogram directory
                prm_dict - dict, metadata from the PRM file
                fbases   - list(str), data file base name to look for
    Returns:    prm_dict - dict, metadata from the PRM file
    """
    # grab an arbitrary file in geo/radar-coordiantes
    geo_file = get_ifg_file(ifg_dir, coord_type='geo')
    if not geo_file:
        raise ValueError(f'No radar-coord files found in {ifg_dir} with suffix: {fbases}')

    # read width from rdr_file
    print(f'grab SLANT_RANGE_DISTANCE, INCIDENCE_ANGLE from file: {geo_file}')
    ds = gdal.Open(geo_file, gdal.GA_ReadOnly)
    width = ds.RasterXSize

    near_range = float(prm_dict['STARTING_RANGE'])
    range_pixel_size = float(prm_dict['RANGE_PIXEL_SIZE'])
    slant_range_dist = near_range + range_pixel_size * width / 2.
    prm_dict['SLANT_RANGE_DISTANCE'] = slant_range_dist

    return prm_dict


def read_baseline_table(fname):
    """Read GMTSAR baseline table file.

    Example: baseline_table.dat file
    #  file_ID         yyyyddd.fraction day_cnt    b_para           b_perp
    S1_20141231_ALL_F2 2014364.5885345591 364  -58.814087517826  -16.031404204594
    S1_20150301_ALL_F2 2015059.5885222249 424  -33.273566547687  -17.102628888348
    S1_20150325_ALL_F2 2015083.5885250464 448 -104.664278131966 -142.776284597889

    Parameters: fname - str, path to the baseline table file
    Returns:    bdict - dict, perpendicular baseline dictionary
    """
    print(f'read baseline time-series from file: {fname}')
    fc = np.loadtxt(fname, dtype=str, usecols=(1,4))
    date_list = ptime.yyyyddd2yyyymmdd([x[:7] for x in fc[:,0]])
    bperp_list = fc[:,1].astype(np.float32)
    bdict = {}
    for date_str, bperp in zip(date_list, bperp_list):
        bdict[date_str] = bperp
    return bdict


#########################################################################
def extract_gmtsar_metadata(meta_file, unw_file, template_file, rsc_file=None, update_mode=True):
    """Extract metadata from GMTSAR interferogram stack.

    Parameters: meta_file     - str, path to the reference metadata file in PRM format
                unw_file      - str, path to the unwrapped interferogram file
                template_file - str, path to the template file
                rsc_file      - str, path to the output metadata file in RSC format
    Returns:    meta          - dict, extract common metadata
    """

    # update_mode: check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=unw_file, readable=False) == 'skip':
        return readfile.read_roipac_rsc(rsc_file)

    ifg_dir = os.path.dirname(unw_file)

    # 1. read *.PRM file
    # prm_file = get_prm_files(ifg_dir)[0]
    print(f'extract metadata from metadata file: {meta_file}')
    meta = readfile.read_gmtsar_prm(meta_file)
    meta['PROCESSOR'] = 'gmtsar'

    # 2. grab A/RLOOKS from config or template file
    config_file = os.path.abspath(os.path.join(os.path.dirname(meta_file), 'config.tops.txt'))
    template = readfile.read_template(template_file)

    if os.path.isfile(config_file):
        # read config file
        config = readfile.read_gmtsar_prm(config_file)
        filter_wvl = float(config.get('filter_wavelength', 200))
        # calc a/rlooks
        # Note from Xiaohua Xu: GMTSAR is applying Gaussian filtering, instead of multilooing,
        # thus, the equivalent value is ~0.3 times the mainlobe response of the boxcar filtering.
        print(f'grab A/RLOOKS from config file: {config_file}')
        meta['ALOOKS'] = np.rint(0.3 * filter_wvl / meta['AZIMUTH_PIXEL_SIZE']).astype(int)
        meta['RLOOKS'] = np.rint(0.3 * filter_wvl / meta['RANGE_PIXEL_SIZE']).astype(int)

    else:
        print(f'WARNING: No config file found in {config_file}!')
        # read from template file
        for key in ['ALOOKS', 'RLOOKS']:
            if key in template.keys():
                print(f'grab {key} from template file: {template_file}')
                meta[key] = int(template[key])
            else:
                msg = f'Attribute {key} is missing! '
                msg += f'Please manually specify it in the template file: {template_file}'
                raise ValueError(msg)

    # 2. grab HEADING from template file
    for key in ['HEADING']:
        if key in template.keys():
            print(f'grab {key} from template file: {template_file}')
            meta[key] = template[key].lower()
        else:
            msg = f'Attribute {key} is missing! '
            msg += f'Please manually specify it in the template file: {template_file}'
            raise ValueError(msg)

    # 3. update AZIMUTH/RANGE_PIXEL_SIZE
    meta['AZIMUTH_PIXEL_SIZE'] *= int(meta['ALOOKS'])
    meta['RANGE_PIXEL_SIZE'] *= int(meta['RLOOKS'])

    # 4. grab LAT/LON_REF1/2/3/4 from geo-coord data file
    meta = get_lalo_ref(ifg_dir, meta)

    # 5. grab X/Y_FIRST/STEP from unw_file if in geo-coord
    ds = gdal.Open(unw_file, gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()
    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.
    if 1e-7 < x_step < 1.:
        print(f'grab Y/X_FIRST/STEP from {unw_file}')
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


def prepare_geometry(template, meta, update_mode=True):
    """Prepare .rsc file for all geometry files.

    Parameters: template - dict, input file path/pattern
                meta     - dict, common metadata for the entire stack
    """
    # grab all specified geometry files
    geom_files = []
    key_prefix = 'mintpy.load.'
    for key in ['demFile', 'lookupYFile', 'lookupXFile', 'incAngleFile',
                'azAngleFile', 'shadowMaskFile', 'waterMaskFile']:
        if key_prefix + key in template.keys() and template[key_prefix + key]:
            geom_files.append(glob.glob(template[key_prefix + key])[0])
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


def prepare_stack(unw_files, meta, bperp_dict, update_mode=True):
    """Prepare .rsc file for all unwrapped interferogram files.

    Parameters: unw_files  - list(str), path to the unwrapped interferogram files
                meta       - dict, common metadata for the entire stack
                bperp_dict - dict, perpendicular baseline time series
    """
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
        date1, date2 = os.path.basename(ifg_dir).replace('-', '_').split('_')
        date1 = ptime.yyyyddd2yyyymmdd(date1)
        date2 = ptime.yyyyddd2yyyymmdd(date2)
        ifg_meta['DATE12'] = f'{ptime.yymmdd(date1)}-{ptime.yymmdd(date2)}'

        # and P_BASELINE_TOP/BOTTOM_HDR
        bperp = bperp_dict[date2] - bperp_dict[date1]
        ifg_meta['P_BASELINE_TOP_HDR'] = bperp
        ifg_meta['P_BASELINE_BOTTOM_HDR'] = bperp

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
    template = ut.check_template_auto_value(template)

    inps.meta_file        = glob.glob(template['mintpy.load.metaFile'])[0]
    inps.bperp_file       = glob.glob(template['mintpy.load.baselineDir'])[0]
    inps.unw_files = sorted(glob.glob(template['mintpy.load.unwFile']))
    inps.cor_files = sorted(glob.glob(template['mintpy.load.corFile']))

    # extract common metadata
    rsc_file = os.path.join(os.path.dirname(inps.meta_file), 'data.rsc')
    meta = extract_gmtsar_metadata(
        meta_file=inps.meta_file,
        unw_file=inps.unw_files[0],
        template_file=inps.template_file,
        rsc_file=rsc_file,
        update_mode=inps.update_mode,
    )

    # read baseline
    bperp_dict = read_baseline_table(inps.bperp_file)

    # prepare metadata for geometry files
    prepare_geometry(
        template=template,
        meta=meta,
        update_mode=inps.update_mode,
    )

    # prepare metadata for interferogram files
    prepare_stack(
        unw_files=inps.unw_files,
        meta=meta,
        bperp_dict=bperp_dict,
        update_mode=inps.update_mode,
    )

    print('Done.')
    return
