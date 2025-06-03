"""Utilities for IONEX products."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jun 2022                           #
############################################################
# Recommend import:
#   from mintpy.objects import ionex
# Links:
#   IGS (NASA): https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
#   IMPC (DLR): https://impc.dlr.de/products/total-electron-content/near-real-time-tec/near-real-time-tec-maps-global


import datetime as dt
import netrc
import os
import re
import subprocess

import numpy as np
import requests
from scipy.interpolate import RegularGridInterpolator, interpn

from mintpy.utils import ptime

IGS_SOLUTION_NAMES = {
    'cod' : 'CODE (1-hour)',
    'esa' : 'ESA (2-hour)',
    'igs' : 'IGS (2-hour)',
    'jpl' : 'JPL (2-hour)',
    'upc' : 'UPC (2-hour)',
    'uqr' : 'UPC (15-min; rapid)'
}

# https://cddis.nasa.gov/archive/gnss/products/ionex/2025/009/
# accessed at June 3, 2025
SOL_CODE2SMP = {
    'CAS' : '30M',
    'COD' : '01H',
    'ESA' : '02H',
    'IGS' : '02H',
    'JPL' : '02H',
    'UPC' : '02H',
}

# Constants
EARTH_DATA_AUTH_HOST = "urs.earthdata.nasa.gov"


################################## Download ####################################

def get_ionex_filename(date_str, tec_dir=None, sol_code='jpl', date_fmt='%Y%m%d'):
    """Get the file name of IONEX files.

    Parameters: date_str - str, date in date_fmt format
                tec_dir  - str, path to the local TEC file directory
                           Set to None for the http path.
                sol_code - str, GIM analysis center code in 3 digit
                           https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
                date_fmt - str, date format code
                           https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    Returns:    tec_file - str, path to the local uncompressed TEC file [for tec_dir != None]
                             or url  to the remote  compressed TEC file [for tec_dir == None]
    """
    # basic time info
    dd = dt.datetime.strptime(date_str, date_fmt)
    yy = str(dd.year)[2:4]
    yyyy = str(dd.year)
    doy = f'{dd.timetuple().tm_yday:03d}'

    # file name
    # transition to the IGS long product filename convention sicne GPS week 2238 (Nov 27, 2022)
    # reference:
    #   https://files.igs.org/pub/resource/guidelines/Guidelines_for_Long_Product_Filenames_in_the_IGS_v2.0.pdf
    #   https://igs.org/products/#ionosphere
    #   https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html#iono
    # The following use the IGS PDF file as the official reference, as the NASA website above (accessed on June 3, 2025)
    #   did not describe the new convention correctly.
    GPS_WEEK_2238 = dt.datetime(2022,11,27)
    if dd >= GPS_WEEK_2238:
        # Name convention since GPS week 2238
        #   {YYYY}/{DDD}/{AAA}0OPS{TTT}_{YYYYDDDHHMM}_01D_SMP_CNT.INX.gz
        #   where YYYY - 4-digit year
        #         DDD  - 3-digit day of year
        #         AAA  - analysis center name
        #         TTT  - solution type identifier, FIN or RAP
        #         SMP  - temporal product sampling resolution
        #         CNT  - content type, GIM or ROT
        SMP  = SOL_CODE2SMP[sol_code.upper()]
        remote_name = f"{sol_code.upper()}0OPSFIN_{yyyy}{doy}0000_01D_{SMP}_GIM.INX.gz"
    else:
        # Name format before GPS week 2237
        #   {YYYY}/{DDD}/{AAA}g{DDD}{#}.{YY}i.Z  for vertical TEC maps
        #   {YYYY}/{DDD}/roti{DDD}0.{YY}f.Z      for Daily ROTI (rate of TEC index) product
        #   where YYYY - 4-digit year
        #         DDD  - 3-digit day of year
        #         AAA  - analysis center name
        #         #    - file number for the day, typically 0
        #         YY   - 2-digit year
        remote_name = f"{sol_code.lower()}g{doy}0.{yy}i.Z"

    # full path / url
    if tec_dir:
        # local file path
        tec_file = os.path.join(tec_dir, remote_name)
        # get the uncompressed file name for the local file
        tec_file = os.path.splitext(tec_file)[0]
    else:
        # remote file url
        url_dir = f"https://cddis.nasa.gov/archive/gnss/products/ionex/{yyyy}/{doy}"
        tec_file = os.path.join(url_dir, remote_name)

    return tec_file


class SessionWithHeaderRedirection(requests.Session):
    """Class to maintain headers after EarthData Login redirect.

    This code was adapted from the examples available here:
    https://urs.earthdata.nasa.gov/documentation/for_users/data_access/python
    """

    def __init__(self, username, password, auth_host=EARTH_DATA_AUTH_HOST):
        """Initialize the SessionWithHeaderRedirection class."""
        super().__init__()
        self.auth = (username, password)
        self.auth_host = auth_host


def dload_ionex(date_str: str, tec_dir: str, sol_code='jpl', date_fmt='%Y%m%d'):
    """Download IGS vertical TEC files in IONEX format.

    This code was adapted from the OPERA repo here:
    https://github.com/opera-adt/disp-s1/blob/main/src/disp_s1/ionosphere.py

    Parameters: date_str - str, date in date_fmt format.
                tec_dir  - str, local directory to save the downloaded files.
                sol_code - str, IGS TEC analysis center code.
                date_fmt - str, date format code
    Returns:    out_path - str, path to the local uncompressed IONEX file.
    """

    # get the source (remote) and destination (local) file path/url
    url = get_ionex_filename(date_str, tec_dir=None, sol_code=sol_code, date_fmt=date_fmt)
    fpath = os.path.join(tec_dir, os.path.basename(url))
    out_path = os.path.splitext(fpath)[0]

    # check existing local file
    if os.path.isfile(fpath) and os.path.getsize(fpath) > 100e3:
        return out_path

    # download: grab authentication creditionals
    print(f'downloading {url} to {os.path.dirname(fpath)} ...')
    try:
        parsed = netrc.netrc().authenticators(EARTH_DATA_AUTH_HOST)
        assert parsed is not None
        username, _, password = parsed
    except (FileNotFoundError, TypeError):
        raise ValueError("No authentication credentials found in ~/.netrc!")

    # download: create session with the user credentials
    session = SessionWithHeaderRedirection(username, password)

    try:
        # submit the request using the session
        response = session.get(url, stream=True)
        # raise an exception in case of http errors
        response.raise_for_status()

        with open(fpath, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)

    except requests.exceptions.HTTPError as e:
        # handle any errors here
        print(e)

    # uncompress the downloaded file
    # Note: the .Z format is based on the LZW compression used by the UNIX compress command, and Python
    # does not have built-in support for this format, so we call the uncompress command line tool first,
    # then read the uncompressed file in plain text.
    # Run the command line tools if uncompressed file:
    #   1) does not exist or
    #   2) smaller than 600k in size or
    #   3) older than the compressed file.
    if (not os.path.isfile(out_path)
            or os.path.getsize(out_path) < 600e3
            or os.path.getmtime(out_path) < os.path.getmtime(fpath)):
        result = subprocess.run(
            ["gunzip", "--decompress", "--keep", str(fpath)], check=True, capture_output=True
        )

    return out_path


#################################### Read ######################################

def get_ionex_value(tec_file, utc_sec, lat, lon, interp_method='linear3d', rotate_tec_map=True,
                    print_msg=True):
    """Get the TEC value from input IONEX file for the input lat/lon/datetime.

    Parameters: tec_file       - str, path of local TEC file
                utc_sec        - float or 1D np.ndarray, UTC time of the day in seconds
                lat/lon        - float or 1D np.ndarray, latitude / longitude in degrees
                interp_method  - str, interpolation method
                rotate_tec_map - bool, rotate the TEC map along the SUN direction.
                                 for "interp_method = linear3d" only.
                print_msg      - bool, print out progress bar or not.
    Returns:    tec_val        - float or 1D np.ndarray, vertical TEC value in TECU
    """

    # time info
    utc_min = utc_sec / 60.

    # read TEC file
    mins, lats, lons, tec_maps = read_ionex(tec_file)[:4]

    # resample
    tec_val = interp_3d_maps(
        tec_maps,
        mins, lats, lons,
        utc_min, lat, lon,
        interp_method=interp_method,
        rotate_tec_map=rotate_tec_map,
        print_msg=print_msg,
    )

    return tec_val


def interp_3d_maps(tec_maps, mins, lats, lons, utc_min, lat, lon, interp_method='linear3d',
                   rotate_tec_map=True, print_msg=True):
    """Interpolate the 3D matrix at given time and location.

    Reference:
        Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format
        version 1.1. Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany.

    Parameters: tec_maps       - 3D np.ndarray in size of (num_time, num_lat, num_lon)
                mins           - 1D np.ndarray in size of (num_time), UTC time of the day in minutes
                lats/lons      - 1D np.ndarray in size of (num_lat/lon), latitude/longitude in degrees
                utc_min        - float or 1D np.ndarray, UTC time of the day in minutes
                lat/lon        - float or 1D np.ndarray, latitude / longitude in degrees
                interp_method  - str, interpolation method.
                rotate_tec_map - bool, rotate the TEC map along the SUN direction.
                                 for "interp_method = linear3d" only (Schaer et al., 1998).
                print_msg      - bool, print out progress bar or not.
    Returns:    tec_val        - float or 1D np.ndarray, TEC value at the given time/location(s).
    """

    def interp_3d_rotate(interpfs, mins, lats, lons, utc_min, lat, lon):
        """Linear interpolation in space/time with rotation along longitude direction,
        following Schaer et al. (1998).
        """
        ind0 = np.where((mins - utc_min) <= 0)[0][-1]
        ind1 = ind0 + 1
        lon0 = lon + (utc_min - mins[ind0]) * 360. / (24. * 60.)
        lon1 = lon + (utc_min - mins[ind1]) * 360. / (24. * 60.)
        tec_val0 = interpfs[ind0]((lon0, lat))
        tec_val1 = interpfs[ind1]((lon1, lat))
        tec_val = (  (mins[ind1] - utc_min) / (mins[ind1] - mins[ind0]) * tec_val0
                   + (utc_min - mins[ind0]) / (mins[ind1] - mins[ind0]) * tec_val1 )
        return tec_val

    # resample
    if interp_method == 'nearest':
        lon_ind = np.abs(lons - lon).argmin()
        lat_ind = np.abs(lats - lat).argmin()
        time_ind = np.abs(mins - utc_min).argmin()
        tec_val = tec_maps[time_ind, lat_ind, lon_ind]

    elif interp_method in ['linear', 'linear2d', 'bilinear']:
        time_ind = np.abs(mins.reshape(-1,1) - utc_min).argmin(axis=0)

        if isinstance(utc_min, np.ndarray):
            num_pts = len(utc_min)
            tec_val = np.zeros(num_pts, dtype=np.float32)
            prog_bar = ptime.progressBar(maxValue=num_pts, print_msg=print_msg)
            for i in range(num_pts):
                tec_val[i] = RegularGridInterpolator(
                    points=(lons, lats),
                    values=tec_maps[time_ind[i], :, :].T,
                    method='linear',
                )((lon[i], lat[i]))

                prog_bar.update(i+1, every=200)
            prog_bar.close()

        else:
            tec_val = RegularGridInterpolator(
                points=(lons, lats),
                values=tec_maps[time_ind[0], :, :].T,
                method='linear',
            )((lon, lat))

    elif interp_method in ['linear3d', 'trilinear']:
        if not rotate_tec_map:
            # option 1: interpolate between consecutive TEC maps
            tec_val = interpn(
                points=(mins, np.flip(lats), lons),
                values=np.flip(tec_maps, axis=1),
                xi=(utc_min, lat, lon),
                method='linear',
            )

        else:
            # option 2: interpolate between consecutive rotated TEC maps
            # reference: equation (3) in Schaer et al. (1998)

            # prepare interpolation functions in advance to speed up
            interpfs = []
            for i in range(len(mins)):
                interpfs.append(
                    RegularGridInterpolator(
                        points=(lons, lats),
                        values=tec_maps[i, :, :].T,
                        method='linear',
                    ),
                )

            # interpolate
            if isinstance(utc_min, np.ndarray):
                num_pts = len(utc_min)
                tec_val = np.zeros(num_pts, dtype=np.float32)
                prog_bar = ptime.progressBar(maxValue=num_pts, print_msg=print_msg)
                for i in range(num_pts):
                    tec_val[i] = interp_3d_rotate(
                        interpfs,
                        mins, lats, lons,
                        utc_min[i], lat[i], lon[i],
                    )
                    prog_bar.update(i+1, every=200)
                prog_bar.close()

            else:
                tec_val = interp_3d_rotate(
                    interpfs,
                    mins, lats, lons,
                    utc_min, lat, lon,
                )

    else:
        msg = f'Un-recognized interp_method input: {interp_method}!'
        msg += '\nSupported inputs: nearest, linear2d, linear3d.'
        raise ValueError(msg)

    return tec_val


def datetime_to_gps_week(dd):
    """Calculates the GPS week number from a datetime object.
    Parameters: dd       - datetime, a datetime object.
    Returns:    gps_week - int, the GPS week number.
    """
    time_delta = dd - dt.datetime(1980, 1, 6)
    gps_week = time_delta.days // 7
    return gps_week


def get_ionex_date(tec_file, date_fmt='%Y-%m-%d'):
    """Get the date in the specified format from the input IONEX filename.

    Parameters: tec_file - str, path to the TEC file in IONEX format
                date_fmt - str, date format code
    Returns:    date_str - str, date in date_fmt format
                date_obj - datetime.datetime object
    """
    fbase = os.path.basename(tec_file)
    year = fbase.split('.')[1][:2]
    doy = fbase.split('.')[0].split('g')[1][:3]
    date_obj = dt.datetime.strptime(year, '%y') + dt.timedelta(days=int(doy)-1)
    date_str = dt.datetime.strftime(date_obj, date_fmt)
    return date_str, date_obj


def get_ionex_height(tec_file):
    """Get the height of the thin-shell ionosphere from IONEX file.

    Parameters: tec_file - str, path to the TEC file in IONEX format
    Returns:    ion_hgt  - float, height above the surface in meters
    """

    with open(tec_file) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip().endswith('DHGT'):
                ion_hgt = float(line.split()[0])
                break

    return ion_hgt


def read_ionex(tec_file):
    """Read TEC file in IONEX format.

    Parameters: tec_file - str, path to the TEC file in IONEX format
    Returns:    mins     - 1D np.ndarray in size of (num_map), time of the day in minutes
                lats     - 1D np.ndarray in size of (num_lat), latitude  in degrees
                lons     - 1D np.ndarray in size of (num_lon), longitude in degrees
                tec_maps - 3D np.ndarray in size of (num_map, num_lat, num_lon), vertical TEC in TECU
                rms_maps - 3D np.ndarray in size of (num_map, num_lat, num_lon), vertical TEC RMS in TECU
    Examples:   tec_dir = os.path.expanduser('~/data/aux/IONEX')
                tec_file = get_ionex_filename('20190519', tec_dir=tec_dir, sol_code='jpl')
                mins, lats, lons, tec_maps = read_ionex(tec_file)[:4]
    """

    # functions for parsing ionex file
    # link: https://github.com/daniestevez/jupyter_notebooks/blob/master/IONEX.ipynb
    def parse_map(tec_map, key='TEC', exponent=-1):
        tec_map = re.split(f'.*END OF {key} MAP', tec_map)[0]
        tec_map = [np.fromstring(x, sep=' ') for x in re.split('.*LAT/LON1/LON2/DLON/H\\n', tec_map)[1:]]
        return np.stack(tec_map) * 10**exponent

    # read IONEX file
    with open(tec_file) as f:
        fc = f.read()

        # read header
        header = fc.split('END OF HEADER')[0].split('\n')
        for line in header:
            if line.strip().endswith('# OF MAPS IN FILE'):
                num_map = int(line.split()[0])
            elif line.strip().endswith('DLAT'):
                lat0, lat1, lat_step = (float(x) for x in line.split()[:3])
            elif line.strip().endswith('DLON'):
                lon0, lon1, lon_step = (float(x) for x in line.split()[:3])
            elif line.strip().endswith('EXPONENT'):
                exponent = float(line.split()[0])

        # spatial coordinates
        num_lat = int((lat1 - lat0) / lat_step + 1)
        num_lon = int((lon1 - lon0) / lon_step + 1)
        lats = np.arange(lat0, lat0 + num_lat * lat_step, lat_step)
        lons = np.arange(lon0, lon0 + num_lon * lon_step, lon_step)

        # time stamps
        min_step = 24 * 60 / (num_map - 1)
        mins = np.arange(0, num_map * min_step, min_step)

        # read TEC and its RMS maps
        tec_maps = np.array([parse_map(t, key='TEC', exponent=exponent)
                             for t in fc.split('START OF TEC MAP')[1:]], dtype=np.float32)
        rms_maps = np.array([parse_map(t, key='RMS', exponent=exponent)
                             for t in fc.split('START OF RMS MAP')[1:]], dtype=np.float32)

    return mins, lats, lons, tec_maps, rms_maps


#################################### Plot ######################################

def plot_ionex(tec_file, save_fig=False):
    """Plot the IONEX file as animation.

    Parameters: tec_file - str, path to the local uncompressed IONEX file
                save_fig - bool, save the animation to file
    Returns:    out_fig  - str, path to the output animation file
    """
    from cartopy import crs as ccrs
    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    from mintpy.utils.map import draw_lalo_label, round_to_1

    # read TEC file
    sol_code = os.path.basename(tec_file)[:3]
    sol_name = IGS_SOLUTION_NAMES.get(sol_code, 'Unknown')
    mins, lats, lons, tec_maps = read_ionex(tec_file)[:4]

    # basic info
    num_map = len(mins)
    lat_step = np.median(np.diff(lats))
    lon_step = np.median(np.diff(lons))
    N = np.max(lats) - lat_step / 2;  S = np.min(lats) + lat_step / 2
    W = np.min(lons) - lon_step / 2;  E = np.max(lons) + lon_step / 2
    vmax = round_to_1(np.nanmax(tec_maps) * 0.9)
    date_obj = get_ionex_date(tec_file)[1]

    # init figure
    proj_obj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=[9, 4], subplot_kw=dict(projection=proj_obj))
    im = ax.imshow(
        tec_maps[0,:,:], vmin=0, vmax=vmax,
        extent=(W, E, S, N), origin='upper',
        animated=True, interpolation='nearest',
    )

    ax.coastlines()
    draw_lalo_label(
        ax,
        geo_box=(W, N, E, S),
        projection=proj_obj,
        print_msg=False,
    )

    # colorbar
    cbar = fig.colorbar(im, shrink=0.5)
    cbar.set_label('TECU')
    fig.tight_layout()

    # update image
    global ind
    ind = 0

    def animate(*args):
        global ind
        ind += 1
        if ind >= num_map:
            ind -= num_map

        # update image & title
        im.set_array(tec_maps[ind,:,:])
        dt_obj = date_obj + dt.timedelta(minutes=mins[ind])
        ax.set_title(f'{dt_obj.isoformat()} - {sol_name}')
        return im,

    # play animation
    ani = FuncAnimation(fig, animate, interval=200, blit=True, save_count=num_map)

    # output
    out_fig = f'{os.path.abspath(tec_file)}.gif'
    if save_fig:
        print(f'saving animation to {out_fig}')
        ani.save(out_fig, writer='imagemagick', dpi=300)

    print('showing animation ...')
    plt.show()

    return out_fig
