############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: David Bekaert, March 2026                        #
############################################################


import glob
import netrc
import os
import re
from datetime import datetime, timedelta

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import mintpy.cli.diff
from mintpy.objects import timeseries
from mintpy.utils import ptime, readfile, utils as ut, writefile


OPERA_MODEL_HOURS = (0, 6, 12, 18)
OPERA_MIN_ACQ_DATE = '20160101'
OPERA_TROPO_COLLECTIONS = [
    'C3717139408-ASF',
    'C1273910987-ASF',
]


def _model_time_token(date_str, hour_str):
    """Build OPERA model timestamp token as YYYYMMDDTHHMMSSZ."""
    return f'{date_str}T{int(hour_str):02d}0000Z'


def _is_valid_opera_file(opera_file):
    """Return True if OPERA file is readable and has required datasets."""
    try:
        with h5py.File(opera_file, 'r') as fobj:
            req_dsets = ['latitude', 'longitude', 'height', 'wet_delay', 'hydrostatic_delay']
            return all(k in fobj for k in req_dsets)
    except Exception:
        return False


def get_geom_lat_lon_bounds(geom_file):
    """Get geometry file lat/lon bounds in degrees."""
    from mintpy.tropo_pyaps3 import get_bounding_box

    meta = readfile.read_attribute(geom_file)
    lat0, lat1, lon0, lon1 = get_bounding_box(meta, geom_file=geom_file)
    return min(lat0, lat1), max(lat0, lat1), min(lon0, lon1), max(lon0, lon1)


def get_opera_crop_indices(lat, lon, geom_file, pad_cells=3):
    """Get OPERA cube crop indices from geometry bounds with fixed ±0.2 deg padding.

    The pad_cells argument is kept for backward compatibility and ignored.
    """
    if lat.size < 2 or lon.size < 2:
        raise ValueError('Invalid OPERA latitude/longitude arrays: size < 2')

    lat_min, lat_max, lon_min, lon_max = get_geom_lat_lon_bounds(geom_file)

    buffer_deg = 0.2
    lat_min_pad = lat_min - buffer_deg
    lat_max_pad = lat_max + buffer_deg
    lon_min_pad = lon_min - buffer_deg
    lon_max_pad = lon_max + buffer_deg

    lat_idx = np.where((lat >= lat_min_pad) & (lat <= lat_max_pad))[0]
    lon_idx = np.where((lon >= lon_min_pad) & (lon <= lon_max_pad))[0]
    if lat_idx.size == 0 or lon_idx.size == 0:
        raise ValueError('No overlap between OPERA grid and MintPy geometry bounds')

    i0, i1 = int(lat_idx[0]), int(lat_idx[-1]) + 1
    j0, j1 = int(lon_idx[0]), int(lon_idx[-1]) + 1

    crop_info = {
        'lat_min': lat_min,
        'lat_max': lat_max,
        'lon_min': lon_min,
        'lon_max': lon_max,
        'lat_min_pad': lat_min_pad,
        'lat_max_pad': lat_max_pad,
        'lon_min_pad': lon_min_pad,
        'lon_max_pad': lon_max_pad,
        'buffer_deg': buffer_deg,
    }
    return (i0, i1, j0, j1), crop_info


def read_opera_total_delay_cube(opera_file, geom_file, pad_cells=3):
    """Read and crop OPERA delay cube, returning total zenith delay.

    Returns a dict with keys:
      latitude, longitude, height, wet_delay, hydro_delay, total_delay, crop_info
    """
    with h5py.File(opera_file, 'r') as fobj:
        req_dsets = ['latitude', 'longitude', 'height', 'wet_delay', 'hydrostatic_delay']
        missing = [k for k in req_dsets if k not in fobj]
        if missing:
            raise ValueError(f'OPERA file missing expected datasets: {missing}')

        lat = fobj['latitude'][:]
        lon = fobj['longitude'][:]
        height = fobj['height'][:]

        (i0, i1, j0, j1), crop_info = get_opera_crop_indices(lat, lon, geom_file, pad_cells=pad_cells)

        lat_crop = lat[i0:i1]
        lon_crop = lon[j0:j1]

        wet = fobj['wet_delay'][0, :, i0:i1, j0:j1]
        hydro = fobj['hydrostatic_delay'][0, :, i0:i1, j0:j1]
        total = wet + hydro

    return {
        'latitude': lat_crop,
        'longitude': lon_crop,
        'height': height,
        'wet_delay': wet,
        'hydro_delay': hydro,
        'total_delay': total,
        'crop_info': crop_info,
    }


def get_geom_lat_lon_dem(geom_file):
    """Read DEM and latitude/longitude grids from geometry file.

    Returns: lat2d/lon2d/dem in 2D np.ndarray with same shape
    """
    dem = readfile.read(geom_file, datasetName='height')[0].astype(np.float32)
    atr = readfile.read_attribute(geom_file)

    dset_list = readfile.get_dataset_list(geom_file)
    if 'latitude' in dset_list and 'longitude' in dset_list:
        lat2d = readfile.read(geom_file, datasetName='latitude')[0].astype(np.float64)
        lon2d = readfile.read(geom_file, datasetName='longitude')[0].astype(np.float64)
    elif 'Y_FIRST' in atr and 'X_FIRST' in atr:
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])
        lat0 = float(atr['Y_FIRST'])
        lon0 = float(atr['X_FIRST'])
        dlat = float(atr['Y_STEP'])
        dlon = float(atr['X_STEP'])
        lat1d = lat0 + dlat * np.arange(length, dtype=np.float64)
        lon1d = lon0 + dlon * np.arange(width, dtype=np.float64)
        lon2d, lat2d = np.meshgrid(lon1d, lat1d)
    else:
        raise ValueError('Unable to get latitude/longitude grids from geometry file!')

    return lat2d, lon2d, dem


def calc_zenith_delay_from_opera_file(opera_file, geom_file, pad_cells=3):
    """Calculate 2D zenith tropospheric delay map intersected with DEM."""
    cube = read_opera_total_delay_cube(opera_file, geom_file, pad_cells=pad_cells)
    lat2d, lon2d, dem = get_geom_lat_lon_dem(geom_file)

    z_axis = np.asarray(cube['height'], dtype=np.float64)
    y_axis = np.asarray(cube['latitude'], dtype=np.float64)
    x_axis = np.asarray(cube['longitude'], dtype=np.float64)
    data = np.asarray(cube['total_delay'], dtype=np.float32)

    # RegularGridInterpolator requires strictly ascending axes.
    if z_axis[1] < z_axis[0]:
        z_axis = z_axis[::-1]
        data = data[::-1, :, :]
    if y_axis[1] < y_axis[0]:
        y_axis = y_axis[::-1]
        data = data[:, ::-1, :]
    if x_axis[1] < x_axis[0]:
        x_axis = x_axis[::-1]
        data = data[:, :, ::-1]

    interp = RegularGridInterpolator(
        (z_axis, y_axis, x_axis),
        data,
        method='cubic',
        bounds_error=False,
        fill_value=np.nan,
    )

    points = np.column_stack([
        dem.reshape(-1).astype(np.float64),
        lat2d.reshape(-1).astype(np.float64),
        lon2d.reshape(-1).astype(np.float64),
    ])
    ztd = interp(points).reshape(dem.shape).astype(np.float32)

    # mask invalid pixels from DEM (NaN/Inf); preserve valid zero-elevation pixels
    ztd[~np.isfinite(dem)] = np.nan

    return ztd, cube


############################################################################
def read_inps2date_time(inps):
    """Read acquisition date/time info from input arguments.

    Parameters: inps      - Namespace for input arguments
    Returns:    date_list - list(str), acquisition dates in YYYYMMDD
                utc_sec   - float, acquisition UTC time in seconds
    """
    if not inps.dis_file:
        raise ValueError('input displacement file is required to derive OPERA date/time info.')

    print(f'read date/time info from file: {inps.dis_file}')
    atr = readfile.read_attribute(inps.dis_file)
    ftype = atr['FILE_TYPE']

    if ftype == 'timeseries':
        date_list = timeseries(inps.dis_file).get_date_list()
    elif ftype == '.unw':
        date_list = ptime.yyyymmdd(atr['DATE12'].split('-'))
    else:
        raise ValueError(f'un-supported displacement file type: {ftype}')

    if 'CENTER_LINE_UTC' not in atr:
        raise ValueError(f'CENTER_LINE_UTC is missing in metadata of file: {inps.dis_file}')

    utc_sec = float(atr['CENTER_LINE_UTC'])

    inps.date_list = date_list
    inps.utc_sec = utc_sec

    return date_list, utc_sec


def nearest_opera_product_time(date_str, utc_sec):
    """Round one acquisition date/time to the nearest OPERA model datetime.

    Parameters: date_str   - str, acquisition date in YYYYMMDD
                utc_sec    - float, acquisition UTC time in seconds
    Returns:    out_date   - str, OPERA model date in YYYYMMDD
                out_hour   - str, OPERA model hour in HH format
    """
    acq_dt = datetime.strptime(date_str, '%Y%m%d') + timedelta(seconds=float(utc_sec))

    candidates = []
    for day_offset in (-1, 0, 1):
        base_day = (acq_dt + timedelta(days=day_offset)).replace(hour=0, minute=0, second=0, microsecond=0)
        for hour in OPERA_MODEL_HOURS:
            cand_dt = base_day + timedelta(hours=hour)
            candidates.append(cand_dt)

    opera_dt = min(candidates, key=lambda x: (abs(x - acq_dt), x))
    out_date = opera_dt.strftime('%Y%m%d')
    out_hour = opera_dt.strftime('%H')
    return out_date, out_hour


def get_opera_date_time_list(date_list, utc_sec):
    """Build required OPERA date/hour list for all acquisitions.

    Parameters: date_list       - list(str), acquisition dates in YYYYMMDD
                utc_sec         - float, acquisition UTC time in seconds
    Returns:    opera_date_list - list(str), OPERA model date in YYYYMMDD
                opera_hour_list - list(str), OPERA model hour in HH format
    """
    opera_date_list = []
    opera_hour_list = []
    for date_str in date_list:
        opera_date, opera_hour = nearest_opera_product_time(date_str, utc_sec)
        opera_date_list.append(opera_date)
        opera_hour_list.append(opera_hour)

    return opera_date_list, opera_hour_list


def get_expected_opera_file_patterns(opera_date_list, opera_hour_list):
    """Create expected OPERA ZTD filename patterns for model times.

    The production timestamp is left as wildcard.

    Returns: expected_patterns - list(str), e.g.
             OPERA_L4_TROPO-ZENITH_YYYYMMDDTHHMMSSZ_*_HRES_v*
    """
    expected_patterns = sorted({
        f'OPERA_L4_TROPO-ZENITH_{_model_time_token(date_str, hour_str)}_*_HRES_v*'
        for date_str, hour_str in zip(opera_date_list, opera_hour_list)
    })
    return expected_patterns


def get_opera_file_status(opera_date_list, opera_hour_list, opera_dir):
    """Get matched OPERA files and missing model date/hour list.

    Returns: expected_patterns      - list(str)
             matched_files          - dict(str, list(str))
             missing_date_hour_list - list(tuple(str, str)) in (YYYYMMDD, HH)
    """
    expected_patterns = get_expected_opera_file_patterns(opera_date_list, opera_hour_list)

    matched_files = {}
    for pattern in expected_patterns:
        files = sorted(glob.glob(os.path.join(opera_dir, pattern)))
        valid_files = [f for f in files if _is_valid_opera_file(f)]
        bad_files = [f for f in files if f not in valid_files]
        if len(bad_files) > 0:
            print(f'WARNING: ignore corrupted/invalid OPERA files for pattern {pattern}:')
            for f in bad_files:
                print(f'  {os.path.basename(f)}')
        matched_files[pattern] = valid_files

    missing_list = []
    seen = set()
    for date_str, hour_str in zip(opera_date_list, opera_hour_list):
        pattern = f'OPERA_L4_TROPO-ZENITH_{_model_time_token(date_str, hour_str)}_*_HRES_v*'
        if len(matched_files.get(pattern, [])) == 0:
            key = (date_str, f'{int(hour_str):02d}')
            if key not in seen:
                missing_list.append(key)
                seen.add(key)
    return expected_patterns, matched_files, missing_list


def _asf_product_text(product):
    """Get searchable text from ASF product object."""
    text_items = [str(product)]

    for key in ['sceneName', 'fileName', 'granuleName', 'url', 'downloadUrl', 'fileID']:
        val = getattr(product, key, None)
        if val:
            text_items.append(str(val))

    props = getattr(product, 'properties', None)
    if isinstance(props, dict):
        for key in ['sceneName', 'fileName', 'granuleName', 'url', 'downloadUrl', 'fileID']:
            val = props.get(key, None)
            if val:
                text_items.append(str(val))

        add_urls = props.get('additionalUrls', None)
        if isinstance(add_urls, (list, tuple)):
            for val in add_urls:
                if val:
                    text_items.append(str(val))

    return ' '.join(text_items)


def _parse_opera_time_tokens(text):
    """Parse model/production time tokens from OPERA filename-like text."""
    mobj = re.search(
        r'OPERA_L4_TROPO-ZENITH_(\d{8}T\d{6}Z)_(\d{8}T\d{6}Z)_HRES_v',
        text,
    )
    if not mobj:
        return None, None
    return mobj.group(1), mobj.group(2)


def _create_asf_session(asf):
    """Create authenticated ASF session using ~/.netrc when possible."""
    if not hasattr(asf, 'ASFSession'):
        return None

    session = asf.ASFSession()

    # Newer asf_search API
    if hasattr(session, 'auth_with_netrc'):
        session.auth_with_netrc()
        return session

    # Backward-compatible fallback
    if hasattr(session, 'auth_with_creds'):
        host = 'urs.earthdata.nasa.gov'
        parsed = netrc.netrc().authenticators(host)
        if not parsed:
            raise ValueError(f'No authentication credentials found in ~/.netrc for {host}')

        username = parsed[0]
        password = parsed[2]
        session.auth_with_creds(username, password)
        return session

    return None


def dload_opera_files(missing_date_hour_list, opera_dir):
    """Download missing OPERA TROPO-ZENITH files via ASF Search.

    Uses per-day ASF queries constrained by OPERA collection IDs and
    date windows, then down-filters by required model-time tokens.
    """
    if len(missing_date_hour_list) == 0:
        return

    try:
        import asf_search as asf
    except Exception as exc:
        raise ImportError(
            'Missing dependency: asf_search. Install it to enable OPERA download.'
        ) from exc

    os.makedirs(opera_dir, exist_ok=True)

    print('-' * 50)
    print('search/download missing OPERA TROPO-ZENITH products from ASF ...')

    # Assume user has ~/.netrc configured. Create session once when available.
    session = None
    try:
        session = _create_asf_session(asf)
        if session is not None:
            print('  ASF authentication: enabled (from ~/.netrc)')
    except Exception as exc:
        print(f'  WARNING: ASF authentication setup failed: {exc}')
        session = None

    missing_tokens = sorted({_model_time_token(d, h) for d, h in missing_date_hour_list})
    date_to_tokens = {}
    for token in missing_tokens:
        date_to_tokens.setdefault(token[:8], []).append(token)

    n_selected = 0
    prog_bar = ptime.progressBar(maxValue=len(missing_tokens), prefix='  download ')

    print('  token match summary (token: candidates -> selected):')
    n_done = 0
    for date_str in sorted(date_to_tokens.keys()):
        day_tokens = date_to_tokens[date_str]
        start_str = f'{date_str[0:4]}-{date_str[4:6]}-{date_str[6:8]}T00:00:00Z'
        end_str = f'{date_str[0:4]}-{date_str[4:6]}-{date_str[6:8]}T23:59:59Z'

        opts = asf.ASFSearchOptions(**{
            'maxResults': 250,
            'collections': OPERA_TROPO_COLLECTIONS,
            'collectionAlias': False,
            'start': start_str,
            'end': end_str,
        })
        if session is not None:
            opts.session = session

        results = asf.search(opts=opts)

        token_best_product = {}
        token_best_prod_time = {}
        token_match_count = {token: 0 for token in day_tokens}

        for product in results:
            text = _asf_product_text(product)
            if 'OPERA_L4_TROPO-ZENITH_' not in text:
                continue

            model_token, prod_token = _parse_opera_time_tokens(text)
            if model_token not in token_match_count:
                continue

            token_match_count[model_token] += 1
            new_prod_time = prod_token or ''
            old_prod_time = token_best_prod_time.get(model_token, '')
            if model_token not in token_best_product or new_prod_time > old_prod_time:
                token_best_product[model_token] = product
                token_best_prod_time[model_token] = new_prod_time

        for token in day_tokens:
            if token not in token_best_product:
                n_done += 1
                prog_bar.update(n_done, suffix=f'{token} (no match)')
                print(f'    {token}: {token_match_count[token]} -> no')
                continue

            dload_results = asf.ASFSearchResults([token_best_product[token]])
            if session is not None:
                dload_results.download(path=opera_dir, session=session)
            else:
                dload_results.download(path=opera_dir)

            n_selected += 1
            n_done += 1
            prog_bar.update(n_done, suffix=f'{token} (downloaded)')
            print(f'    {token}: {token_match_count[token]} -> yes')

    prog_bar.close()

    print(f'  missing model-time tokens downloaded     : {n_selected} / {len(missing_tokens)}')
    print(f'  downloaded files to: {opera_dir}')
    print('-' * 50)


############################################################################
def _get_best_local_opera_file(date_str, hour_str, opera_dir):
    token = _model_time_token(date_str, hour_str)
    pattern = os.path.join(opera_dir, f'OPERA_L4_TROPO-ZENITH_{token}_*_HRES_v*.nc')
    files = sorted(f for f in glob.glob(pattern) if _is_valid_opera_file(f))
    if len(files) == 0:
        return None
    return files[-1]


def _project_zenith_to_los(zenith_delay, inc_angle_deg):
    cos_inc = np.cos(np.deg2rad(inc_angle_deg.astype(np.float64)))
    out = zenith_delay.astype(np.float32).copy()
    invalid = (~np.isfinite(cos_inc)) | (np.abs(cos_inc) < 1.0e-8)
    out[invalid] = np.nan
    out[~invalid] = out[~invalid] / cos_inc[~invalid]
    out *= -1
    return out


def _split_supported_dates(date_list, min_date=OPERA_MIN_ACQ_DATE):
    """Split acquisition dates into OPERA-supported and unsupported sets."""
    supported = [d for d in date_list if d >= min_date]
    unsupported = [d for d in date_list if d < min_date]
    return supported, unsupported


def _run_or_skip_opera(opera_files, used_dates, tropo_file, geom_file):
    """Run/skip decision for OPERA delay file generation (MintPy-style update mode)."""

    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        return (atr['LENGTH'], atr['WIDTH'])

    print('update mode: ON')
    print(f'output file: {tropo_file}')
    flag = 'skip'

    if ut.run_or_skip(out_file=tropo_file, in_file=opera_files, print_msg=False) == 'run':
        flag = 'run'
        print('1) output file either does NOT exist or is NOT newer than all OPERA files.')
    else:
        print('1) output file exists and is newer than all OPERA files.')

        if (get_dataset_size(tropo_file) != get_dataset_size(geom_file)
                or any(i not in timeseries(tropo_file).get_date_list() for i in used_dates)):
            flag = 'run'
            print(('2) output file does NOT have the same len/wid as the geometry file {}'
                   ' or does NOT contain all required dates').format(geom_file))
        else:
            print('2) output file has the same len/wid as the geometry file and contains all required dates.')
            with h5py.File(tropo_file, 'r') as fobj:
                if np.all(fobj['timeseries'][-1, :, :] == 0):
                    flag = 'run'
                    print('3) output file is NOT fully written.')
                else:
                    print('3) output file is fully written.')

    print(f'run or skip: {flag}')
    return flag


############################################################################
def calculate_delay_timeseries(tropo_file, dis_file, geom_file, opera_dir):
    """Calculate OPERA delay time-series and write to HDF5 file."""
    atr = readfile.read_attribute(dis_file)
    ftype = atr['FILE_TYPE']

    if ftype == 'timeseries':
        date_list = timeseries(dis_file).get_date_list()
    elif ftype == '.unw':
        date12 = atr['DATE12']
        date_list = ptime.yyyymmdd(date12.split('-'))
    else:
        raise ValueError(f'un-supported displacement file type: {ftype}')

    if 'CENTER_LINE_UTC' not in atr:
        raise ValueError(f'CENTER_LINE_UTC is missing in metadata of file: {dis_file}')
    utc_sec = float(atr['CENTER_LINE_UTC'])

    supported_dates, unsupported_dates = _split_supported_dates(date_list, min_date=OPERA_MIN_ACQ_DATE)
    if len(unsupported_dates) > 0:
        print((f'WARNING: {len(unsupported_dates)} acquisition date(s) are earlier than '
               f'OPERA availability ({OPERA_MIN_ACQ_DATE}) and will be skipped.'))
        for d in unsupported_dates:
            print(f'  {d}')

    if len(supported_dates) == 0:
        raise RuntimeError(
            f'No acquisition dates are supported by OPERA (requires >= {OPERA_MIN_ACQ_DATE}).'
        )

    model_dates, model_hours = get_opera_date_time_list(supported_dates, utc_sec)

    opera_files = []
    used_dates = []
    for acq_date, model_date, model_hour in zip(supported_dates, model_dates, model_hours):
        fname = _get_best_local_opera_file(model_date, model_hour, opera_dir)
        if fname is None:
            print(f'WARNING: NO OPERA file found for acquisition {acq_date} -> {model_date} {model_hour}:00')
            continue
        used_dates.append(acq_date)
        opera_files.append(fname)

    if len(opera_files) == 0:
        raise RuntimeError('No OPERA files found for requested acquisitions. Cannot create delay time-series.')

    if _run_or_skip_opera(opera_files, used_dates, tropo_file, geom_file) == 'skip':
        return tropo_file

    atr_out = atr.copy()
    atr_out['FILE_TYPE'] = 'timeseries'
    atr_out['UNIT'] = 'm'
    for key in ['REF_DATE', 'REF_X', 'REF_Y', 'REF_LAT', 'REF_LON']:
        if key in atr_out:
            atr_out.pop(key)

    length, width = int(atr_out['LENGTH']), int(atr_out['WIDTH'])
    num_date = len(opera_files)
    dates = np.array(used_dates, dtype=np.bytes_)
    ds_name_dict = {
        'date': [dates.dtype, (num_date,), dates],
        'timeseries': [np.float32, (num_date, length, width), None],
    }
    writefile.layout_hdf5(tropo_file, ds_name_dict, metadata=atr_out)

    print(f'read incidenceAngle from file: {geom_file}')
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0].astype(np.float32)

    prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        ztd, _ = calc_zenith_delay_from_opera_file(opera_files[i], geom_file, pad_cells=3)
        delay_los = _project_zenith_to_los(ztd, inc_angle)

        block = [i, i + 1, 0, length, 0, width]
        writefile.write_hdf5_block(
            tropo_file,
            data=delay_los,
            datasetName='timeseries',
            block=block,
            print_msg=False,
        )

        prog_bar.update(i + 1, suffix=os.path.basename(opera_files[i]))
    prog_bar.close()

    return tropo_file


############################################################################
def run_tropo_opera(inps):
    """Run OPERA-based tropospheric correction workflow."""
    print('tropospheric delay correction with OPERA approach')
    print(f'input displacement file : {inps.dis_file}')
    print(f'input geometry file     : {inps.geom_file}')
    print(f'input OPERA directory   : {inps.opera_dir}')
    print(f'output tropo file       : {inps.tropo_file}')
    print(f'output corrected file   : {inps.cor_dis_file}')

    date_list, utc_sec = read_inps2date_time(inps)
    supported_dates, unsupported_dates = _split_supported_dates(date_list, min_date=OPERA_MIN_ACQ_DATE)
    if len(unsupported_dates) > 0:
        print((f'WARNING: {len(unsupported_dates)} acquisition date(s) are earlier than '
               f'OPERA availability ({OPERA_MIN_ACQ_DATE}) and will be skipped.'))
        for d in unsupported_dates:
            print(f'  {d}')

    if len(supported_dates) == 0:
        print(f'No acquisition dates >= {OPERA_MIN_ACQ_DATE}. Skip OPERA correction.')
        return

    opera_date_list, opera_hour_list = get_opera_date_time_list(supported_dates, utc_sec)
    expected_patterns, matched_files, missing_date_hour_list = get_opera_file_status(
        opera_date_list=opera_date_list,
        opera_hour_list=opera_hour_list,
        opera_dir=inps.opera_dir,
    )

    inps.opera_date_list = opera_date_list
    inps.opera_hour_list = opera_hour_list
    inps.expected_opera_patterns = expected_patterns
    inps.matched_opera_files = matched_files
    inps.missing_opera_date_hour_list = missing_date_hour_list

    if len(missing_date_hour_list) > 0:
        dload_opera_files(missing_date_hour_list, inps.opera_dir)
        expected_patterns, matched_files, missing_date_hour_list = get_opera_file_status(
            opera_date_list=opera_date_list,
            opera_hour_list=opera_hour_list,
            opera_dir=inps.opera_dir,
        )
        inps.expected_opera_patterns = expected_patterns
        inps.matched_opera_files = matched_files
        inps.missing_opera_date_hour_list = missing_date_hour_list

    print(f'number of acquisitions (input)    : {len(date_list)}')
    print(f'number of acquisitions (supported): {len(supported_dates)}')
    print(f'acquisition UTC seconds: {utc_sec}')
    if len(set(opera_hour_list)) == 1:
        print(f'nearest OPERA model time: {opera_hour_list[0]}:00 UTC')
    else:
        print(f'nearest OPERA model hours: {sorted(set(opera_hour_list))}')

    n_total = len(expected_patterns)
    n_found = sum(len(i) > 0 for i in matched_files.values())
    n_missing = len(missing_date_hour_list)
    print(f'expected OPERA model files (unique): {n_total}')
    print(f'expected OPERA model files found   : {n_found}')
    print(f'OPERA model files still missing    : {n_missing}')
    if n_missing > 0:
        print('missing model date/hour list (YYYYMMDD, HH):')
        for date_str, hour_str in missing_date_hour_list:
            print(f'  {date_str}, {hour_str}')

    if n_found == 0:
        print('No OPERA files available to calculate delay. Stop.')
        return

    calculate_delay_timeseries(
        tropo_file=inps.tropo_file,
        dis_file=inps.dis_file,
        geom_file=inps.geom_file,
        opera_dir=inps.opera_dir,
    )

    print('correcting delay for using diff.py')
    iargs = [inps.dis_file, inps.tropo_file, '-o', inps.cor_dis_file, '--force']
    print('diff.py', ' '.join(iargs))
    mintpy.cli.diff.main(iargs)

    return
