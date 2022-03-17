#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, Feb 2020           #
############################################################
# Useful links:
#   IGS (at NASA's archive): https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
#   IMPC (DLR): https://impc.dlr.de/products/total-electron-content/near-real-time-tec/nrt-tec-global/
# Contents
#   Ionospheric Mapping Functions
#   JPL high resolution GIM - I/O
#   IGS low  resolution GIM - download, I/O, plot
#   Test
# Recommend usage:
#   from mintpy.simulation import iono


import os
import sys
import re
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from mintpy.utils import (
    ptime,
    readfile,
    utils0 as ut,
)


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
#   Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., & Aoki, Y. (2022). Range Geolocation Accuracy
#     of C/L-band SAR and its Implications for Operational Stack Coregistration. IEEE Trans. Geosci. Remote Sens. 
#   Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format version 1.1. 
#     Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany, Darmstadt, Germany.

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
    #print('refraction angle on the ionosphere min/max: {:.1f}/{:.1f} deg'.format(np.nanmin(ref_angle), np.nanmax(ref_angle)))
    tec = vtec / np.cos(ref_angle * np.pi / 180.0)

    # calculate range delay based on equation (1) in Chen and Zebker (2012)
    range_delay = (tec * 1e16 * K / (freq**2)).astype(np.float32)

    # group delay = phase advance * -1
    if obs_type != 'phase':
        range_delay *= -1.

    return range_delay


def iono_incidence2refraction_angle(inc_angle, vtec, freq):
    """Calculate the refraction angle for the ionospheric shell.

    Equation (8) in Yunjun et al. (2022, TGRS)

    Parameters: inc_angle - float / np.ndarray, incidence angle in deg
                vtec      - float, zenith TEC in TECU
                freq      - float, radar carrier frequency in Hz.
    Returns:    ref_angle - float / np.ndarray, refraction angle in deg
    """
    Ne = vtec * 1e16
    # equation (26) in Bohm & Schuh (2013) Chapter 2.
    n_iono_group = 1 + K * Ne / freq**2
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
        print('center lat/lon  on the ground    : {:.4f}/{:.4f} deg'.format(lat, lon))
        print('center lat/lon  on the ionosphere: {:.4f}/{:.4f} deg'.format(iono_lat, iono_lon))

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
    if type(inc_angle) is np.ndarray:
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
        hm0 = '{:02d}{:02d}'.format(h0, m0)
        hm1 = '{:02d}{:02d}'.format(h1, m1)
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
            print('WARNING: NO file found in {}. Set to NaN.'.format(tec_file))

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
        tec_file = os.path.join(gim_tec_dir, 'subtec_sent1_{}.txt.dt'.format(date_str))
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
        raise FileNotFoundError('No file found in {}'.format(tec_file))

    if print_msg:
        print('read JPL GIM TEC data from file: {}'.format(tec_file))

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
    date_str = re.findall('\d{8}', os.path.basename(tec_file))[0]
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
    1. remove dates that are not existed in the reference date list
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



################################## IGS TEC #####################################
# Low resolution ionospheric TEC products from IGS
#   Including solutions from various centers, e.g., CODE, JPL, etc.
#   Webiste: https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
#   Spatial  resolution in latitude / longitude [deg]: 2.5 / 5.0
#   Temporal resolution [hour]: 2.0 / 1.0

def calc_igs_iono_ramp(tec_dir, date_str, geom_file, box=None, print_msg=True):
    """Get 2D ionospheric delay from IGS TEC data for one acquisition.
    due to the variation of the incidence angle along LOS.

    Parameters: tec_dir        - str, path of the local TEC directory, i.e. ~/data/aux/IGS_TEC
                date_str       - str, date of interest in YYYYMMDD format
                geom_file      - str, path of the geometry file including incidenceAngle data
                box            - tuple of 4 int, subset in (x0, y0, x1, y1)
    Returns:    range_delay    - 2D np.ndarray for the range delay in meters
                vtec           - float, TEC value in zenith direction in TECU
                iono_lat/lon   - float, latitude/longitude in degrees of LOS vector in ionosphere shell
                iono_height    - float, height             in meter   of LOS vector in ionosphere shell
                iono_inc_angle - float, incidence angle    in degrees of LOS vector in ionosphere shell
    """
    # geometry
    tec_file = dload_igs_tec(date_str, tec_dir, print_msg=print_msg)
    iono_height = grab_ionex_height(tec_file)
    (iono_inc_angle,
     iono_lat,
     iono_lon) = prep_geometry_iono(geom_file, box=box, iono_height=iono_height)[:3]

    # time
    meta = readfile.read_attribute(geom_file)
    utc_sec = float(meta['CENTER_LINE_UTC'])
    if print_msg:
        h, s = divmod(utc_sec, 3600)
        m, s = divmod(s, 60)
        print('UTC time: {:02.0f}:{:02.0f}:{:02.1f}'.format(h, m, s))

    # extract zenith TEC
    freq = SPEED_OF_LIGHT / float(meta['WAVELENGTH'])
    vtec = get_igs_tec_value(tec_file, utc_sec, iono_lon, iono_lat)
    rang_delay = vtec2range_delay(vtec, iono_inc_angle, freq)

    return rang_delay, vtec, iono_lat, iono_lon, iono_height, iono_inc_angle


def get_igs_tec_value(tec_file, utc_sec, lat, lon, interp_method='linear3d', rotate_tec_map=False, print_msg=True):
    """Get the TEC value based on input lat/lon/datetime
    Parameters: tec_file - str, path of local TEC file
                utc_sec  - float or 1D np.ndarray, UTC time of the day in seconds
                lat/lon  - float or 1D np.ndarray, latitude / longitude in degrees
                interp_method  - str, interpolation method
                rotate_tec_map - bool, rotate the TEC map along the SUN direction, for linear3d only.
                print_msg      - bool, print out progress bar or not.
    Returns:    tec_val  - float or 1D np.ndarray, TEC value in TECU
    """

    def interp_3d_rotate(interpfs, lons, lats, mins, lon, lat, utc_min):
        ind0 = np.where((mins - utc_min) <=0)[0][-1]
        ind1 = ind0 + 1
        lon0 = lon + (utc_min - mins[ind0]) * 360. / (24. * 60.)
        lon1 = lon + (utc_min - mins[ind1]) * 360. / (24. * 60.)
        tec_val0 = interpfs[ind0](lon0, lat)
        tec_val1 = interpfs[ind1](lon1, lat)
        tec_val = (  (mins[ind1] - utc_min) / (mins[ind1] - mins[ind0]) * tec_val0
                   + (utc_min - mins[ind0]) / (mins[ind1] - mins[ind0]) * tec_val1 )
        return tec_val

    # read TEC file
    lons, lats, mins, tecs = read_ionex_tec(tec_file)[:4]
    tec_maps = tecs[0]

    # time info
    utc_min = utc_sec / 60.

    # resample
    if interp_method == 'nearest':
        lon_ind = np.abs(lons - lon).argmin()
        lat_ind = np.abs(lats - lat).argmin()
        time_ind = np.abs(mins - utc_min).argmin()
        tec_val = tec_maps[lon_ind, lat_ind, time_ind]

    elif interp_method in ['linear', 'linear2d', 'bilinear']:
        time_ind = np.abs(mins.reshape(-1,1) - utc_min).argmin(axis=0)
        if isinstance(time_ind, np.ndarray):
            num = len(time_ind)
            tec_val = np.zeros(num, dtype=np.float32)
            prog_bar = ptime.progressBar(maxValue=num, print_msg=print_msg)
            for i in range(num):
                tec_val[i] = interpolate.interp2d(lons, lats, tec_maps[:, :, time_ind[i]].T, kind='linear')(lon[i], lat[i])
                prog_bar.update(i+1, every=200)
            prog_bar.close()
        else:
            tec_val = interpolate.interp2d(lons, lats, tec_maps[:, :, time_ind].T, kind='linear')(lon, lat)

    elif interp_method in ['linear3d', 'trilinear']:
        if not rotate_tec_map:
            # option 1: interpolate between consecutive TEC maps
            # testings shows better agreement with SAR obs than option 2.
            tec_val = interpolate.interpn((lons, np.flip(lats), mins),
                                          np.flip(tec_maps, axis=1),
                                          (lon, lat, utc_min),
                                          method='linear')

        else:
            # option 2: interpolate between consecutive rotated TEC maps
            # reference: equation (3) in Schaer and Gurtner (1998)

            # prepare interpolation functions in advance to speed up
            interpfs = []
            for i in range(len(mins)):
                interpfs.append(interpolate.interp2d(lons, lats, tec_maps[:, :, i].T, kind='linear'))

            if isinstance(utc_min, np.ndarray):
                num = len(utc_min)
                tec_val = np.zeros(num, dtype=np.float32)
                prog_bar = ptime.progressBar(maxValue=num, print_msg=print_msg)
                for i in range(num):
                    tec_val[i] = interp_3d_rotate(interpfs, lons, lats, mins, lon[i], lat[i], utc_min[i])
                    prog_bar.update(i+1, every=200)
                prog_bar.close()
            else:
                tec_val = interp_3d_rotate(interpfs, lons, lats, mins, lon, lat, utc_min)

    return tec_val


def get_igs_tec_filename(tec_dir, date_str, sol='jpl', datefmt='%Y%m%d'):
    """Get the local file name of downloaded IGS TEC product."""
    dd = dt.datetime.strptime(date_str, datefmt)
    doy = '{:03d}'.format(dd.timetuple().tm_yday)
    yy = str(dd.year)[2:4]

    fname = "{a}g{d}0.{y}i.Z".format(a=sol.lower(), d=doy, y=yy)
    fbase = fname[:-2]
    tec_file = os.path.join(tec_dir, fbase)
    return tec_file


def dload_igs_tec(d, out_dir, sol='jpl', datefmt='%Y%m%d', print_msg=False):
    """Download IGS vertical TEC files computed by JPL
    Link: https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
    """
    # date info
    dd = dt.datetime.strptime(d, datefmt)
    doy = '{:03d}'.format(dd.timetuple().tm_yday)
    yy = str(dd.year)[2:4]

    fbase = "{a}g{d}0.{y}i.Z".format(a=sol.lower(), d=doy, y=yy)
    src_dir = "https://cddis.nasa.gov/archive/gnss/products/ionex/{0}/{1}".format(dd.year, doy)

    # input/output filename
    fname_src = os.path.join(src_dir, fbase)
    fname_dst = os.path.join(out_dir, fbase)
    fname_dst_uncomp = fname_dst[:-2]

    # download
    cmd = 'wget --continue --auth-no-challenge "{}"'.format(fname_src)
    if os.path.isfile(fname_dst) and os.path.getsize(fname_dst) > 1000:
        cmd += ' --timestamping '
    if not print_msg:
        cmd += ' --quiet '
    else:
        print(cmd)

    # run cmd in output dir
    pwd = os.getcwd()
    os.chdir(out_dir)
    os.system(cmd)
    os.chdir(pwd)

    # uncompress
    # if output file 1) does not exist or 2) smaller than 400k in size or 3) older
    if (not os.path.isfile(fname_dst_uncomp)
            or os.path.getsize(fname_dst_uncomp) < 600e3
            or os.path.getmtime(fname_dst_uncomp) < os.path.getmtime(fname_dst)):
        cmd = "gzip --force --keep --decompress {}".format(fname_dst)
        if print_msg:
            print(cmd)
        os.system(cmd)

    return fname_dst_uncomp


def grab_ionex_height(tec_file):
    """Grab the height of the thin-shell ionosphere from IONEX file"""
    # read ionex file into list of lines
    with open(tec_file, 'r') as f:
        lines = f.readlines()

    # search for height - DHGT
    iono_height = None
    for line in lines:
        c = [i.strip() for i in line.strip().replace('\n', '').split()]
        if c[-1] == 'DHGT':
            iono_height = float(c[0]) * 1e3
            break
    return iono_height


def read_ionex_tec(igs_file):
    """Read IGS TEC file in IONEX format.
    Download from https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html.
    Parameters: igs_file  - str, path of the TEC file
    Returns:    lon       - 1D np.ndarray for the longitude in size of (num_lon,) in degrees
                lat       - 1D np.ndarray for the latitude  in size of (num_lat,) in degrees
                tec_times - 1D np.ndarray for the time of the day in size of (num_map,) in minutes
                tec_array - 4D np.ndarray for the vertical TEC value in size of (2, num_lon, num_lat, num_map) in TECU
                            1 TECU = 10^16 electrons / m^2
                tec_type  - str, TEC solution
    """

    igs_code_dict = {
        'igr' : 'IGS (Rapid)',
        'igs' : 'IGS (Final)',
        'jpr' : 'JPL (Rapid)',
        'jpl' : 'JPL (Final)',
        'cor' : 'CODE (Rapid)',
        'cod' : 'CODE (Final)',
    }
    tec_code = os.path.basename(igs_file)[:3]
    tec_type = igs_code_dict.get(tec_code, None)

    #print(tec_type)
    ## =========================================================================
    ##
    ## The following section reads the lines of the ionex file for 1 day
    ## (13 maps total) into an array a[]. It also retrieves the thin-shell
    ## ionosphere height used by IGS, the lat./long. spacing, etc. for use
    ## later in this script.
    ##
    ## =========================================================================
    #print 'Transfering IONEX data format to a TEC/DTEC array for ',ymd_date

    ## Opening and reading the IONEX file into memory as a list
    linestring = open(igs_file, 'r').read()
    LongList = linestring.split('\n')

    ## Create two lists without the header and only with the TEC and DTEC maps (based on code from ionFR.py)
    AddToList = 0
    TECLongList = []
    DTECLongList = []
    for i in range(len(LongList)-1):
        ## Once LongList[i] gets to the DTEC maps, append DTECLongList
        if LongList[i].split()[-1] == 'MAP':
            if LongList[i].split()[-2] == 'RMS':
                AddToList = 2
        if AddToList == 1:
            TECLongList.append(LongList[i])
        if AddToList == 2:
            DTECLongList.append(LongList[i])
        ## Determine the number of TEC/DTEC maps
        if LongList[i].split()[-1] == 'FILE':
            if LongList[i].split()[-3:-1] == ['MAPS','IN']:
                num_maps = float(LongList[i].split()[0])
        ## Determine the shell ionosphere height (usually 450 km for IGS IONEX files)
        if LongList[i].split()[-1] == 'DHGT':
            ion_H = float(LongList[i].split()[0])
        ## Determine the range in lat. and long. in the ionex file
        if LongList[i].split()[-1] == 'DLAT':
            start_lat = float(LongList[i].split()[0])
            end_lat   = float(LongList[i].split()[1])
            incr_lat  = float(LongList[i].split()[2])
        if LongList[i].split()[-1] == 'DLON':
            start_lon = float(LongList[i].split()[0])
            end_lon   = float(LongList[i].split()[1])
            incr_lon  = float(LongList[i].split()[2])
        ## Find the end of the header so TECLongList can be appended
        if LongList[i].split()[0] == 'END':
            if LongList[i].split()[2] == 'HEADER':
                AddToList = 1

    ## Variables that indicate the number of points in Lat. and Lon.
    points_lon = ((end_lon - start_lon)/incr_lon) + 1
    points_lat = ((end_lat - start_lat)/incr_lat) + 1   ## Note that incr_lat is defined as '-' here
    num_row = int(np.ceil(points_lon/16))               ## Note there are 16 columns of data in IONEX format

    ## 4-D array that will contain TEC & DTEC (a[0] and a[1], respectively) values
    a = np.zeros((2,int(points_lon),int(points_lat),int(num_maps)))

    ## Selecting only the TEC/DTEC values to store in the 4-D array.
    for Titer in range(2):
        counterMaps = 1
        UseList = []
        if Titer == 0:
            UseList = TECLongList
        elif Titer == 1:
            UseList = DTECLongList
        for i in range(len(UseList)):
            ## Pointing to first map (out of 13 maps) then by changing 'counterMaps' the other maps are selected
            if UseList[i].split()[0] == ''+str(counterMaps)+'':
                if UseList[i].split()[-4] == 'START':
                    ## Pointing to the starting Latitude then by changing 'counterLat' we select TEC data
                    ## at other latitudes within the selected map
                    counterLat = 0
                    newstartLat = float(str(start_lat))
                    for iLat in range(int(points_lat)):
                        if UseList[i+2+counterLat].split()[0].split('-')[0] == ''+str(newstartLat)+'':
                            ## Adding to array a[] a line of Latitude TEC data
                            counterLon = 0
                            for row_iter in range(num_row):
                                for item in range(len(UseList[i+3+row_iter+counterLat].split())):
                                    a[Titer,counterLon,iLat,counterMaps-1] = UseList[i+3+row_iter+counterLat].split()[item]
                                    counterLon = counterLon + 1
                        if '-'+UseList[i+2+counterLat].split()[0].split('-')[1] == ''+str(newstartLat)+'':
                            ## Adding to array a[] a line of Latitude TEC data. Same chunk as above but
                            ## in this case we account for the TEC values at negative latitudes
                            counterLon = 0
                            for row_iter in range(num_row):
                                for item in range(len(UseList[i+3+row_iter+counterLat].split())):
                                    a[Titer,counterLon,iLat,counterMaps-1] = UseList[i+3+row_iter+counterLat].split()[item]
                                    counterLon = counterLon + 1
                        counterLat = counterLat + row_iter + 2
                        newstartLat = newstartLat + incr_lat
                    counterMaps = counterMaps + 1

    ## =========================================================================
    ##
    ## The section creates a new array that is a copy of a[], but with the lower
    ## left-hand corner defined as the initial element (whereas a[] has the
    ## upper left-hand corner defined as the initial element).  This also
    ## accounts for the fact that IONEX data is in 0.1*TECU.
    ##
    ## =========================================================================

    ## The native sampling of the IGS maps minutes
    incr_time = 24*60/int(num_maps-1)
    tec_array = np.zeros((2,int(points_lon),int(points_lat),int(num_maps)))

    for Titer in range(2):
        #incr = 0
        for ilat in range(int(points_lat)):
            tec_array[Titer,:,ilat,:] = 0.1*a[Titer,:,int(points_lat-1-ilat),:]

    lat = np.arange(start_lat, start_lat + points_lat*incr_lat, incr_lat)
    lon = np.arange(start_lon, start_lon + points_lon*incr_lon, incr_lon)
    tec_times = np.arange(0, incr_time*num_maps, incr_time)

    return lon, lat, tec_times, tec_array, tec_type


def plot_tec_animation(tec_file, save=False):
    """Plot the input tec file as animation"""
    from cartopy import crs as ccrs
    from matplotlib.animation import FuncAnimation
    from mintpy.utils import plot as pp

    def grab_date(tec_file, datefmt='%Y-%m-%d'):
        """Grab the date in YYYYMMDD format from the TEC filename"""
        tec_file = os.path.basename(tec_file)
        # year
        year = tec_file.split('.')[1][:2]
        if year[0] == '9':
            year = '19' + year
        else:
            year = '20' + year
        year = int(year)

        # month and day
        doy = int(tec_file.split('.')[0].split('g')[1][:3])
        dt_obj = dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)
        date_str = dt_obj.strftime(datefmt)
        return date_str

    # read TEC file
    lon, lat, tec_times, tec_array, tec_type = read_ionex_tec(tec_file)
    vmax = ut.round_to_1(np.nanmax(tec_array[0]) * 0.9)
    # SNWE info
    geo_box = (np.min(lon), np.max(lat), np.max(lon), np.min(lat))  # (W, N, E, S)
    extent = (geo_box[0], geo_box[2], geo_box[3], geo_box[1])       # (W, E, S, N)
    # date/time info
    date_str = grab_date(tec_file)
    num_map = len(tec_times)

    # init figure
    fig, ax = plt.subplots(figsize=[9, 4], subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.coastlines()
    pp.draw_lalo_label(geo_box, ax, projection=ccrs.PlateCarree(), print_msg=False)
    # init image
    im = ax.imshow(tec_array[0,:,:,0].T, vmin=0, vmax=vmax, extent=extent,
                   origin='upper', animated=True, interpolation='nearest')
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

        # update image
        data = tec_array[0,:,:,ind].T
        im.set_array(data)

        # update title
        dt_obj = dt.datetime.strptime(date_str, '%Y-%m-%d') + dt.timedelta(minutes=tec_times[ind])
        title = dt_obj.isoformat()
        if tec_type:
            title += ' - {}'.format(tec_type)
        ax.set_title(title)

        return im,

    # play animation
    ani = FuncAnimation(fig, animate, interval=300, blit=True, save_count=num_map)

    # output
    if save:
        outfig = '{}.gif'.format(os.path.abspath(tec_file))
        print('saving animation to {}'.format(outfig))
        ani.save(outfig, writer='imagemagick', dpi=300)
    print('showing animation ...')
    plt.show()
    return



################################## test ########################################

def main(iargs=None):
    tec_dir = os.path.expanduser('~/data/aux/GIM_IGS')
    tec_file = dload_igs_tec('20190409', tec_dir)
    plot_tec_animation(tec_file, save=True)
    return
if __name__ == '__main__':
    main(sys.argv[1:])
