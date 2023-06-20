"""Ionospheric mapping functions."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Feb 2020                           #
############################################################
# Useful links:
#   IGS (NASA): https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
#   IMPC (DLR): https://impc.dlr.de/products/total-electron-content/near-real-time-tec/nrt-tec-global/
# Recommend usage:
#   from mintpy.simulation import iono


import datetime as dt
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

from mintpy.utils import readfile, utils0 as ut

# constants
SPEED_OF_LIGHT = 299792458 # m/s
EARTH_RADIUS = 6371e3      # Earth radius in meters
K = 40.31                  # m^3/s^2, constant

SAR_BAND = {
    'L' : 1.257e9,  # NISAR-L, ALOS-2
    'S' : 3.2e9,    # NISAR-S
    'C' : 5.405e9,  # Sentinel-1
    'X' : 9.65e9,   # TerraSAR-X
    'Ka': 35.75e9,  # SWOT
}


######################## Ionospheric Mapping Functions #########################
# Reference:
#   Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., & Aoki, Y. (2022). Range
#     Geolocation Accuracy of C-/L-band SAR and its Implications for Operational Stack Coregistration.
#     IEEE Trans. Geosci. Remote Sens., 60, doi:10.1109/TGRS.2022.3168509.

def vtec2range_delay(vtec, inc_angle, freq, obs_type='phase'):
    """Calculate/predict the range delay in SAR from TEC in zenith direction.

    Equation (6-11) from Yunjun et al. (2022).

    Parameters: vtec      - float, zenith TEC in TECU
                inc_angle - float/np.ndarray, incidence angle at the ionospheric shell in deg
                freq      - float, radar carrier frequency in Hz.
                obs_type  - str, given the same iono, the impact on offset (amplitude) and phase is reversed.
    Returns:    rg_delay  - float/np.ndarray, predicted range delay in meters
    """
    # ignore no-data value in inc_angle
    if isinstance(inc_angle, np.ndarray):
        inc_angle[inc_angle == 0] = np.nan

    # convert to TEC in LOS based on equation (3) in Chen and Zebker (2012)
    ref_angle = iono_incidence2refraction_angle(inc_angle, vtec, freq)
    tec = vtec / np.cos(ref_angle * np.pi / 180.0)

    # calculate range delay based on equation (1) in Chen and Zebker (2012)
    range_delay = (tec * 1e16 * K / (freq**2)).astype(np.float32)

    # group delay = phase advance * -1
    if obs_type != 'phase':
        range_delay *= -1.

    return range_delay


def iono_incidence2refraction_angle(inc_angle, vtec, freq):
    """Calculate the refraction angle for the ionospheric shell.

    Reference:
        Equation (8) in Yunjun et al. (2022, TGRS)
        Equation (26) in Bohm & Schuh (2013) Chapter 2.

    Parameters: inc_angle - float / 1/2D np.ndarray, incidence angle in deg
                vtec      - float / 1D np.ndarray, zenith TEC in TECU
                freq      - float, radar carrier frequency in Hz.
    Returns:    ref_angle - float / 1/2/3D np.ndarray, refraction angle in deg
    """
    if isinstance(vtec, np.ndarray):
        # only 1D array is supported
        if vtec.ndim > 1:
            raise ValueError(f'input vtec dimension ({vtec.ndim}) > 1!')

        # tile to ensure same shape between vtec and inc_angle
        if isinstance(inc_angle, np.ndarray) and inc_angle.shape != vtec.shape:
            num_vtec = vtec.size
            if inc_angle.ndim == 1:
                vtec = np.tile(vtec.reshape(-1, 1), (1, inc_angle.size))
                inc_angle = np.tile(inc_angle.reshape(1, -1), (num_vtec, 1))
            elif inc_angle.ndim == 2:
                vtec = np.tile(vtec.reshape(-1, 1, 1), (1, inc_angle.shape[0], inc_angle.shape[1]))
                inc_angle = np.tile(inc_angle[np.newaxis, :, :], (num_vtec, 1, 1))
            else:
                raise ValueError(f'input inc_angle dimension ({inc_angle.ndim}) > 2!')

    # Equation (26) in Bohm & Schuh (2013) Chapter 2.
    n_iono_group = 1 + K * vtec * 1e16 / freq**2

    # Equation (8) in Yunjun et al. (2022, TGRS)
    ref_angle = np.arcsin(1 / n_iono_group * np.sin(inc_angle * np.pi / 180)) * 180 / np.pi

    return ref_angle


def prep_geometry_iono(geom_file, box=None, iono_height=450e3, print_msg=True):
    """Prepare geometry info of LOS vector at thin-shell ionosphere.

    Equation (11-12) in Yunjun et al. (2022, TGRS)

    Parameters: geom_file      - str, path to the geometry file in HDF5/MintPy format
                box            - tuple of 4 int, box of interest in (x0, y0, x1, y1)
                iono_height    - float, height of the assume effective thin-shell ionosphere in m
    Returns:    iono_inc_angle - 2D np.ndarray / float, incidence angle in degree
                iono_lat/lon   - float, latitude/longitude of LOS vector at the thin-shell in deg
                iono_height    - float, height of the assume effective thin-shell ionosphere in m
    """

    def get_center_lat_lon(geom_file, box=None):
        """Get the lat/lon of the scene center"""
        meta = readfile.read_attribute(geom_file)
        if box is None:
            box = (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))

        col_c = int((box[0] + box[2]) / 2)
        row_c = int((box[1] + box[3]) / 2)
        if 'Y_FIRST' in meta.keys():
            lat0 = float(meta['Y_FIRST'])
            lon0 = float(meta['X_FIRST'])
            lat_step = float(meta['Y_STEP'])
            lon_step = float(meta['X_STEP'])
            lat_c = lat0 + lat_step * row_c
            lon_c = lon0 + lon_step * col_c
        else:
            box_c = (col_c, row_c, col_c+1, row_c+1)
            lat_c = float(readfile.read(geom_file, datasetName='latitude',  box=box_c)[0])
            lon_c = float(readfile.read(geom_file, datasetName='longitude', box=box_c)[0])
        return lat_c, lon_c

    # inc_angle on the ground
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle', box=box)[0]
    inc_angle = np.squeeze(inc_angle)
    inc_angle[inc_angle == 0] = np.nan
    inc_angle_center = np.nanmean(inc_angle)
    if print_msg:
        print('incidence angle on the ground     min/max: {:.1f}/{:.1f} deg'.format(np.nanmin(inc_angle),
                                                                                    np.nanmax(inc_angle)))

    # inc_angle on the thin-shell ionosphere - equation (11)
    iono_inc_angle = incidence_angle_ground2iono(inc_angle, iono_height=iono_height)
    if print_msg:
        print('incidence angle on the ionosphere min/max: {:.1f}/{:.1f} deg'.format(np.nanmin(iono_inc_angle),
                                                                                    np.nanmax(iono_inc_angle)))

    # center lat/lon on the ground & thin-shell ionosphere - equation (12)
    lat, lon = get_center_lat_lon(geom_file, box=box)
    az_angle = readfile.read(geom_file, datasetName='azimuthAngle', box=box)[0]
    az_angle[az_angle == 0] = np.nan
    az_angle_center = np.nanmean(az_angle)
    iono_lat, iono_lon = lalo_ground2iono(lat, lon, inc_angle=inc_angle_center, az_angle=az_angle_center)

    if print_msg:
        print(f'center lat/lon  on the ground    : {lat:.4f}/{lon:.4f} deg')
        print(f'center lat/lon  on the ionosphere: {iono_lat:.4f}/{iono_lon:.4f} deg')

    return iono_inc_angle, iono_lat, iono_lon, iono_height


def incidence_angle_ground2iono(inc_angle, iono_height=450e3):
    """Calibrate the incidence angle of LOS vector on the ground surface to the ionosphere shell.

    Equation (11) in Yunjun et al. (2022, TGRS)

    Parameters: inc_angle      - float/np.ndarray, incidence angle on the ground in degrees
                iono_height    - float, effective ionosphere height in meters
                                 under the thin-shell assumption
    Returns:    inc_angle_iono - float/np.ndarray, incidence angle on the iono shell in degrees
    """
    # ignore nodata in inc_angle
    if isinstance(inc_angle, np.ndarray):
        inc_angle[inc_angle == 0] = np.nan

    # deg -> rad & copy to avoid changing input variable
    inc_angle = np.array(inc_angle) * np.pi / 180

    # calculation
    inc_angle_iono = np.arcsin(EARTH_RADIUS * np.sin(inc_angle) / (EARTH_RADIUS + iono_height))
    inc_angle_iono *= 180. / np.pi

    return inc_angle_iono


def lalo_ground2iono(lat, lon, inc_angle, az_angle=None, head_angle=None, iono_height=450e3, method='spherical_distance'):
    """Calculate lat/lon of IPP with the given lat/lon on the ground and LOS geometry

    Equation (12) in Yunjun et al. (2021, TGRS).

    Parameters: lat/lon      - float, latitude/longitude of the point on the ground in degrees
                inc_angle    - float/np.ndarray, incidence angle of the line-of-sight vector on the ground in degrees
                az_angle     - float/np.ndarray, azimuth   angle of the line-of-sight vector
                head_angle   - float, heading angle of the satellite orbit in degrees
                               from the north direction with positive in clockwise direction
                iono_height  - float, height of the ionosphere thin-shell in meters
    Returns:    lat/lon_iono - float, latitude/longitude of the point on the thin-shell ionosphere in degrees
    """
    # LOS incidence angle at ionospheric level
    inc_angle_iono = incidence_angle_ground2iono(inc_angle)

    # geocentric angle distance between the SAR pixel and IPP
    dist_angle = inc_angle - inc_angle_iono

    # LOS azimuth angle (from the ground to the satellite)
    # measured from the north with anti-clockwise as positive (ISCE-2 convention).
    if az_angle is None:
        # heading angle -> azimuth angle
        if head_angle is not None:
            az_angle = ut.heading2azimuth_angle(head_angle)
    if az_angle is None:
        raise ValueError('az_angle can not be None!')

    # In spherical coordinate system, given the starting lat/lon, angular distance and azimuth angle
    #   calculate the ending lat/lon.
    #
    # option 1 - spherical_distance
    # link:
    #   https://gis.stackexchange.com/questions/5821 [there is a typo there]
    #   http://www.movable-type.co.uk/scripts/latlong.html [used]
    if method == 'spherical_distance':
        # deg -> rad & copy to avoid changing input variable
        lat = np.array(lat) * np.pi / 180.
        lon = np.array(lon) * np.pi / 180.
        dist_angle *= np.pi / 180.
        az_angle *= np.pi / 180.

        # calculate
        lat_iono = np.arcsin(  np.sin(lat) * np.cos(dist_angle)
                             + np.cos(lat) * np.sin(dist_angle) * np.cos(az_angle))
        dlon = np.arctan2(np.sin(dist_angle) * np.cos(lat) * np.sin(az_angle) * -1.,
                          np.cos(dist_angle) - np.sin(lat) * np.sin(lat_iono))
        lon_iono = np.mod(lon + dlon + np.pi, 2. * np.pi) - np.pi

        # rad -> deg
        lat_iono *= 180. / np.pi
        lon_iono *= 180. / np.pi

    # option 2 - pyproj
    # link: https://pyproj4.github.io/pyproj/stable/api/geod.html#pyproj.Geod.fwd
    elif method == 'pyproj':
        import pyproj
        geod = pyproj.Geod(ellps='WGS84')
        dist = dist_angle * (np.pi / 180.) * EARTH_RADIUS
        lon_iono, lat_iono = geod.fwd(lon, lat, az=-az_angle, dist=dist)[:2]

    else:
        raise ValueError(f'Un-recognized method: {method}')

    ## obsolete
    ## offset angle from equation (25) in Chen and Zebker (2012)
    ## ~3 degree in magnitude for Sentinel-1 asc track
    #off_iono = inc_angle - np.arcsin(EARTH_RADIUS / (EARTH_RADIUS + iono_height) * np.sin(np.pi - inc_angle))
    #lat += off_iono * np.cos(head_angle) * 180. / np.pi
    #lon += off_iono * np.sin(head_angle) * 180. / np.pi

    return lat_iono, lon_iono



############################ JPL high resolution GIM ###########################
# JPL high resolution GIM
#   Spatial  resolution in latitude / longitude [deg]: 1.0 / 1.0
#   Temporal resolution [min]: 15

def get_gim_tec_list(gim_tec_dir, dt_objs, iono_lat, iono_lon):
    """Get the GIM TEC value at loc (lat/lon) of interest.

    Parameters: gim_tec_dir - str, path to the GIM_TEC directory
                dt_objs     - list of datetime.datetime objects
                iono_lat    - float, latitude  of LOS vector at ionosphere in degree
                iono_lon    - float, longitude of LOS vector at ionosphere in degree
    Returns:    vtec        - 1D np.ndarray of zenith TEC in TECU
    """

    num_date = len(dt_objs)
    vtec = np.zeros(num_date, dtype=np.float32) * np.nan

    for i, dt_obj in enumerate(dt_objs):
        # get gim tec file name
        m_step = 15
        h0 = dt_obj.hour
        m0 = int(np.floor(dt_obj.minute / m_step) * m_step)
        m1 = m0 + m_step
        if m1 < 60:
            h1 = h0
        else:
            m1 = 0
            h1 = h0 + 1
        hm0 = f'{h0:02d}{m0:02d}'
        hm1 = f'{h1:02d}{m1:02d}'
        tec_file = 'gim_{}_{}_{}.tecgrd.txt'.format(dt_obj.strftime('%y%m%d'), hm0, hm1)
        tec_file = os.path.join(gim_tec_dir, tec_file)

        # read tec file
        if os.path.isfile(tec_file):
            lats, lons, vtec_mat = read_gim_tec_grid_file(tec_file)

            # value of the nearest grid
            ind_lat = np.argmin(np.abs(lats - iono_lat))
            ind_lon = np.argmin(np.abs(lons - iono_lon))
            vtec[i] = vtec_mat[ind_lat, ind_lon]
        else:
            print(f'WARNING: NO file found in {tec_file}. Set to NaN.')

    return vtec


def read_gim_tec_grid_file(tec_file, display=False):
    """Read JPL high resolution GIM grid file.

    Parameters: tec_file - str, path of the gim_*.tecgrd.txt file
    Returns:    lat/lon  - 1D np.ndarray in size of 181/361 in degrees
                tec      - 2D np.ndarray in size of (181, 361) in TECU
    """
    lat = np.loadtxt(tec_file, skiprows=2, max_rows=1)
    lon = np.loadtxt(tec_file, skiprows=3, max_rows=1)
    tec = np.loadtxt(tec_file, skiprows=4)

    # plot
    if display:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[10, 4])
        im = ax.imshow(tec, extent=(-180.5, 180.5, -90.5, 90.5))
        # colorbar
        cbar = fig.colorbar(im)
        cbar.set_label('zenith TEC [TECU]')
        # axis format
        ax.set_xlabel('longitude [deg]')
        ax.set_ylabel('latitude [deg]')
        fig.tight_layout()
        # output
        plt.show()

    return lat, lon, tec


def get_sub_tec_list(gim_tec_dir, date_list, iono_lat, iono_lon,
                     interp_method='nearest', fill_value=np.nan):
    """Interpolate / extract SUB/GIM TEC value at loc (lat/lon) of interest.

    Parameters: gim_tec_dir   - str, path of GIM TEC data directorry
                date_list     - list of str, dates in YYYYMMDD format
                iono_lat/lon  - float, latitude / longitude in degree of TOP TEC of interest
                interp_method - str, interpolation method
                fill_value    - float32, the filling value for missing dates
    Returns:    date_list     - list of str, dates in YYYYMMDD format
                tec_ipp       - 1D np.ndarray in float, total       TEC at IPP in TECU
                tec_top_tpp   - 1D np.ndarray in float, topside     TEC at TPP in TECU
                tec_sub_ipp   - 1D np.ndarray in float, sub-orbital TEC at IPP in TECU
                tec_sub_tpp   - 1D np.ndarray in float, sub-orbital TEC at TPP in TECU
                tDicts        - list of dict, each dict contains all info of GIM TEC file
    """
    num_date = len(date_list)
    tec_ipp     = np.zeros(num_date, dtype=np.float32) * fill_value
    tec_top_tpp = np.zeros(num_date, dtype=np.float32) * fill_value
    tec_sub_ipp = np.zeros(num_date, dtype=np.float32) * fill_value
    tec_sub_tpp = np.zeros(num_date, dtype=np.float32) * fill_value

    # read TEC into tDicts (list of dict objects)
    tDicts = []
    for date_str in date_list:
        tec_file = os.path.join(gim_tec_dir, f'subtec_sent1_{date_str}.txt.dt')
        if os.path.isfile(tec_file):
            tDict = read_sub_tec(tec_file, print_msg=False)
        else:
            tDict = None
        tDicts.append(tDict)

    # interpolate TOP TEC at lat/lon of interest
    # tDicts --> gim/top/sub_tec
    for i in range(num_date):
        tDict = tDicts[i]
        if tDict is None:
            continue

        if interp_method == 'nearest':
            ratio_lon2lat = np.cos(iono_lat * np.pi / 180.)
            dist_deg = ( (tDict['lat'] - iono_lat) ** 2 +
                        ((tDict['lon'] - iono_lon) * ratio_lon2lat) ** 2 ) ** 0.5
            ind = np.argmin(dist_deg)
            tec_ipp[i]     = tDict['tec_ipp'][ind]
            tec_top_tpp[i] = tDict['tec_top_tpp'][ind]
            tec_sub_ipp[i] = tDict['tec_sub_ipp'][ind]
            tec_sub_tpp[i] = tDict['tec_sub_tpp'][ind]

        elif interp_method == 'nearestLat':
            ind = np.argmin(np.abs(tDict['lat'] - iono_lat))
            tec_ipp[i]     = tDict['tec_ipp'][ind]
            tec_top_tpp[i] = tDict['tec_top_tpp'][ind]
            tec_sub_ipp[i] = tDict['tec_sub_ipp'][ind]
            tec_sub_tpp[i] = tDict['tec_sub_tpp'][ind]

        elif interp_method == 'linear':
            # NOTE by Yunjun at 2021-08-12: extrapolate give abnormal large value, thus is not good.
            for d_tec, key in zip([ tec_ipp,   tec_top_tpp,   tec_sub_ipp,   tec_sub_tpp ],
                                  ['tec_ipp', 'tec_top_tpp', 'tec_sub_ipp', 'tec_sub_tpp']):
                if tDict['lat'].size <= 2:
                    # use nearest instead
                    ind = np.argmin(np.abs(tDict['lat'] - iono_lat))
                    d_tec[i] = tDict[key][ind]

                else:
                    func = interpolate.interp1d(tDict['lat'], tDict[key],
                                                kind=interp_method,
                                                fill_value='extrapolate',
                                                assume_sorted=False)
                    d_tec[i] = func(iono_lat)

        elif interp_method == 'mean':
            tec_ipp[i]     = np.mean(tDict['tec_ipp'])
            tec_top_tpp[i] = np.mean(tDict['tec_top_tpp'])
            tec_sub_ipp[i] = np.mean(tDict['tec_sub_ipp'])
            tec_sub_tpp[i] = np.mean(tDict['tec_sub_tpp'])

        elif interp_method == 'median':
            tec_ipp[i]     = np.median(tDict['tec_ipp'])
            tec_top_tpp[i] = np.median(tDict['tec_top_tpp'])
            tec_sub_ipp[i] = np.median(tDict['tec_sub_ipp'])
            tec_sub_tpp[i] = np.median(tDict['tec_sub_tpp'])

    return date_list, tec_ipp, tec_top_tpp, tec_sub_ipp, tec_sub_tpp, tDicts


def read_sub_tec(tec_file, version=2.1, print_msg=True):
    """Read the JPL SUB/GIM TEC file.

    Parameters: tec_file - str, path of the sub TEC text file
    Returns:    tDict    - dict, dictionary of GIM/SUB/TOP_TEC data and time/location info

    Acronyms:
        TGT: target  piercing point, where a shell at 450   km is assumed, previously known as IPP.
        TPP: topside piercing point, where a shell at 1,800 km is assumed.
        LEO: low earth orbit, where the SAR satellite is (~700 km).
        LOS: line-of-sight

    Version 2. Format of the data file provided by Xiaoqing Xi in August 6, 2020.
    COL01:  GPSSEC      (GPS seconds)
    COL02:  UThrs       (UTC in hours)
    COL03:  LThrs       (local solar time in hours)
    COL04:  SVN         (space vehicle number)
    COL05:  AZMdeg      (observation azimuth   angle)
    COL06:  ELVdeg      (observation elevation angle)
    COL07:  TPPLATdeg   (TPP latitude)             [*]
    COL08:  TPPLONdeg   (TPP longitude)            [*]
    COL09:  LEOALTkm    (LEO altitude  in km)
    COL10:  LEOLATdeg   (LEO latitude  in deg)
    COL11:  LEOLONdeg   (LEO longitude in deg)
    COL12:  TGTLATdeg   (target latitude)
    COL13:  TGTLONdeg   (target longitude)
    COL14:  STECRMB     (slant TEC, bias removed)
    COL15:  TOPTEC_TPP  (topTEC  at TPP)           [*]
    COL16:  GIMTEC_TPP  (GIM TEC at TPP)
    COL17:  SUBTEC_TPP  (subTEC  at TPP)
    COL18:  GIMTEC_TGT  (GIM TEC at TGT at 450 km) [*]
    COL19:  GIMTEC_LEO  (GIM TEC at the projected LEO coordinates)

    Version 1. Format of the data file provided by Xiaoqing Xi in July 9, 2020.
    COL01: GPSSEC,  GPS seconds past J2000
    COL02: UTC      [sec] within the day
    COL03: LThrs    [hours], local solar time
    COL04: SVN      [integer], GPS space vehicle number
    COL05: AZM      [deg], observation azimuth angle
    COL06: ELV      [deg], observation elevation angle (90 – zenith)
    COL07: SUBLAT   [deg], IPP latitude              [*]
    COL08: SUBLON   [deg], IPP longitude             [*]
    COL09: LEOALT   [km],  LEO altitude
    COL10: LEOLAT   [deg], LEO latitude
    COL11: LEOLON   [deg], LEO longitude
    COL12: STECRMB  [TECU], line-of-sight slant  TEC
    COL13: TOPTEC   [TECU], topside     vertical TEC [*]
    COL14: GIMTRACK [TECU], GIM         vertical TEC [*]
    COL15: SUBTEC   [TECU], sub-orbital vertical TEC [*]
    """

    if not os.path.isfile(tec_file):
        raise FileNotFoundError(f'No file found in {tec_file}')

    if print_msg:
        print(f'read JPL GIM TEC data from file: {tec_file}')

    # read data from text file
    fc = np.loadtxt(tec_file, dtype=bytes).astype(float)
    if fc.ndim == 1:
        fc = fc.reshape(1, -1)

    # read lat/lon of the (topside) ionospheric pierce point (TPP / IPP)
    (lat, lon) = fc[:, 6:8].T
    lon = ut.wrap(lon, wrap_range=(-180, 180))

    # read the topside, sub-orbital and total vertical TEC
    if version == 2.1:
        # version 2.1: GIMTEC_TGT, TOPTEC_TPP and (GIMTEC_TGT - TOPTEC_TPP)
        #tec_tpp     = fc[:, 15].T
        tec_top_tpp = fc[:, 14].T
        tec_sub_tpp = fc[:, 16].T
        tec_ipp     = fc[:, 17].T
        tec_sub_ipp = tec_ipp - tec_top_tpp

    elif version == 2.0:
        # version 2.0: GIMTEC_TPP, TOPTEC_TPP and SUBTEC_TPP
        (TOPTEC,
         GIMTEC,
         SUBTEC) = fc[:, 14:17].T

    elif version == 1.0:
        # version 1.0: GIMTRACK, TOPTEC, SUBTEC
        (TOPTEC,
         GIMTEC,
         SUBTEC) = fc[:, 12:15].T

    # save for external plotting
    date_str = re.findall(r'\d{8}', os.path.basename(tec_file))[0]
    date_obj = dt.datetime.strptime(date_str, '%Y%m%d')
    tDict = {}
    tDict['date'] = date_str
    tDict['time'] = date_obj
    tDict['lat'] = lat
    tDict['lon'] = lon
    tDict['tec_ipp']     = tec_ipp
    tDict['tec_top_tpp'] = tec_top_tpp
    tDict['tec_sub_ipp'] = tec_sub_ipp
    tDict['tec_sub_tpp'] = tec_sub_tpp

    return tDict


def check_date_list_against_reference(date_list, dset_list, date_list_ref, fill_value=np.nan):
    """Check input date/dset_list against the reference date list:
    1. remove dates that are not existing in the reference date list
    2. fill data of the missing dates
    """
    # remove dates not in date_list_ref
    flag = np.array([i in date_list_ref for i in date_list], dtype=np.bool_)
    if np.sum(flag) < len(date_list):
        date_list = np.array(date_list)[flag].tolist()
        dset_list = np.array(dset_list)[flag].tolist()

    # fill missing dates from date_list_ref
    flag = np.array([i in date_list for i in date_list_ref], dtype=np.bool_)
    if np.sum(flag) < len(date_list_ref):
        date_list = list(date_list_ref)
        dset_list_bk = list(dset_list)
        dset_list = np.zeros(len(date_list_ref), dtype=np.float32) * fill_value
        dset_list[flag] = dset_list_bk

    return date_list, dset_list
