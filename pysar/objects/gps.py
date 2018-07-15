############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, Jul 2018                          #
############################################################
# Utility scripts for GPS handling
# Recommend import:
#     from pysar.objects import gps

import os
import time
from datetime import datetime as dt
import numpy as np
from pyproj import Geod
from urllib.request import urlretrieve
from pysar.utils import readfile, utils as ut


unr_site_list_file = 'http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt'


def dload_site_list(print_msg=True):
    """download DataHoldings.txt"""
    url = unr_site_list_file
    out_file = os.path.basename(url)
    if print_msg:
        print('downloading site list from UNR Geod Lab: {}'.format(url))
    urlretrieve(url, out_file)
    return out_file


def search_gps(SNWE, site_list_file=None, print_msg=True):
    """Search available GPS sites within the geo bounding box from UNR website
    Parameters: SNWE : tuple of 4 float, indicating (South, North, West, East) in degrees
                site_list_file : string.
    Returns:    site_names : 1D np.array of string for GPS station names
                site_lats  : 2D np.array  for lat/lon per row
    """
    # download site list file if it's not found in current directory
    if site_list_file is None:
        site_list_file = os.path.basename(unr_site_list_file)
        if not os.path.isfile(site_list_file):
            dload_site_list(print_msg=print_msg)

    txt_data = np.loadtxt(site_list_file,
                          dtype=bytes,
                          skiprows=1,
                          usecols=(0,1,2,3)).astype(str)
    site_names = txt_data[:, 0]
    site_lats, site_lons = txt_data[:, 1:3].astype(np.float32).T
    site_lons  -= np.round(site_lons / (360.)) * 360.

    idx = np.multiply(np.multiply(site_lats >= SNWE[0], site_lats <= SNWE[1]), 
                      np.multiply(site_lons >= SNWE[2], site_lons <= SNWE[3]))
    return site_names[idx], site_lats[idx], site_lons[idx]


def read_gps_los_displacement(site, gps_dir, geom_file, ref_site=None):
    """Read GPS displacement in LOS direction
    Parameters: site : string, GPS station ID, e.g. GV05
                gps_dir : string, local path of GPS files
                geom_flie : string, path of geometry files, such as geo_geometryRadar.h5
    Returns:    gps_times : 1D np.array of datetime.datetime object
                gps_dis   : 1D np.array of displacement in meters
                gps_std   : 1D np.array of displacement uncertainty in meters
    """
    # read GPS object
    gps_obj = gps(site=site, data_dir=gps_dir)
    gps_obj.open(print_msg=False)

    # calculate geometry
    atr = readfile.read_attribute(geom_file)
    coord = ut.coordinate(atr, lookup_file=geom_file)
    y, x = coord.geo2radar(gps_obj.site_lat,
                           gps_obj.site_lon,
                           print_msg=False)[0:2]
    box = (x, y, x+1, y+1)
    inc_angle = readfile.read(geom_file,
                              datasetName='incidenceAngle',
                              box=box,
                              print_msg=False)[0].flatten()
    head_angle = readfile.read(geom_file,
                               datasetName='headingAngle',
                               box=box,
                               print_msg=False)[0].flatten()
    head_angle = -1 * (180 + head_angle + 90)

    # convert ENU to LOS direction
    gps_obj.dis_los = ut.enu2los(gps_obj.dis_e,
                                 gps_obj.dis_n,
                                 gps_obj.dis_u,
                                 inc_angle=inc_angle,
                                 head_angle=head_angle)
    gps_obj.std_los = ut.enu2los(gps_obj.std_e,
                                 gps_obj.std_n,
                                 gps_obj.std_u,
                                 inc_angle=inc_angle,
                                 head_angle=head_angle)
    gps_times = gps_obj.times
    gps_dis_los = gps_obj.dis_los
    gps_std_los = gps_obj.std_los
    site_lalo = (gps_obj.site_lat,
                 gps_obj.site_lon)

    # get LOS displacement relative to another GPS site
    if ref_site:
        ref_gps_obj = gps(site=ref_site, data_dir=gps_dir)
        ref_gps_obj.open(print_msg=False)

        try:
            y, x = coord.lalo2yx(ref_gps_obj.site_lat,
                                 ref_gps_obj.site_lon,
                                 print_msg=False)[0:2]
        except:
            y, x = -1, -1
        if 0 <= x < int(atr['WIDTH']) and 0<= y < int(atr['LENGTH']):
            box = (x, y, x+1, y+1)
            inc_angle = readfile.read(geom_file,
                                      datasetName='incidenceAngle',
                                      box=box,
                                      print_msg=False)[0].flatten()
            head_angle = readfile.read(geom_file,
                                       datasetName='headingAngle',
                                       box=box,
                                       print_msg=False)[0].flatten()
            head_angle = -1 * (180 + head_angle + 90)

        ref_gps_obj.dis_los = ut.enu2los(ref_gps_obj.dis_e,
                                         ref_gps_obj.dis_n,
                                         ref_gps_obj.dis_u,
                                         inc_angle=inc_angle,
                                         head_angle=head_angle)
        ref_gps_obj.std_los = ut.enu2los(ref_gps_obj.std_e,
                                         ref_gps_obj.std_n,
                                         ref_gps_obj.std_u,
                                         inc_angle=inc_angle,
                                         head_angle=head_angle)

        # get relative LOS displacement on common dates
        gps_times = np.array(sorted(list(set(gps_obj.times) & set(ref_gps_obj.times))))
        gps_dis_los = np.zeros(gps_times.shape, np.float32)
        gps_std_los = np.zeros(gps_times.shape, np.float32)
        for i in range(len(gps_times)):
            idx1 = np.where(gps_obj.times == gps_times[i])[0][0]
            idx2 = np.where(ref_gps_obj.times == gps_times[i])[0][0]
            gps_dis_los[i] = gps_obj.dis_los[idx1] - ref_gps_obj.dis_los[idx2]
            gps_std_los[i] = np.sqrt(gps_obj.std_los[idx1]**2 + ref_gps_obj.std_los[idx2]**2)

        ref_site_lalo = (ref_gps_obj.site_lat,
                         ref_gps_obj.site_lon)
    else:
        ref_site_lalo = None

    return gps_times, gps_dis_los, gps_std_los, site_lalo, ref_site_lalo


class gps:
    """GPS class for GPS time-series of daily solution

    Example:
      import matplotlib.pyplot as plt
      from pysar.objects import gps
      from pysar.utils import utils as ut
      gps_obj = gps(site='GV05', data_dir='~/insarlab/GPS')
      gps_obj.open()
      dis_los = ut.enu2los(gps_obj.dis_e,
                           gps_obj.dis_n,
                           gps_obj.dis_u)
      times = gps_obj.times
      plt.figure()
      plt.scatter(times, dis_los)
      plt.show()
    """

    def __init__(self, site, data_dir, data_format='env3'):
        self.site = site
        self.data_dir = data_dir
        self.data_format = data_format
        self.source = 'Nevada Geodetic Lab'
        self.file = os.path.join(data_dir, '{}.IGS08.t{}'.format(site, data_format))
        self.img_file = os.path.join(data_dir, 'PIC/{}.png'.format(site))
        self.site_list_file = os.path.join(data_dir, 'DataHoldings.txt')

    def open(self, print_msg=True):
        if not os.path.isfile(self.file):
            self.dload_site()
        self.get_stat_lat_lon(print_msg=print_msg)
        self.read_displacement(print_msg=print_msg)

    def dload_site(self, print_msg=True):
        # download time-series data from Nevada Geodetic Lab
        url = 'http://geodesy.unr.edu/gps_timeseries'
        if self.data_format == 'env3':
            url = os.path.join(url, 'tenv3')
        elif self.data_format == 'xyz2':
            url = os.path.join(url, 'txyz')
        url = os.path.join(url, 'IGS08/{}'.format(os.path.basename(self.file)))

        if print_msg:
            print('downloading {} from {}'.format(self.site, url))
        urlretrieve(url, self.file)

        # download PNG file
        if not os.path.isdir(os.path.dirname(self.img_file)):
            os.makedirs(os.path.dirname(self.img_file))
        url = 'http://geodesy.unr.edu/tsplots/IGS08/TimeSeries/{}.png'.format(self.site)
        if print_msg:
            print('downloading {}.png from {}'.format(self.site, url))
        urlretrieve(url, self.img_file)

        return self.file

    def get_stat_lat_lon(self, print_msg=True):
        """Get station lat/lon"""
        if print_msg:
            print('calculating station lat/lon')
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)
        ref_lon, ref_lat = float(data[0, 6]), 0.
        e0, e_off, n0, n_off = data[0, 7:11].astype(np.float)
        e0 += e_off
        n0 += n_off

        az = np.arctan2(e0, n0) / np.pi * 180.
        dist = np.sqrt(e0**2 + n0**2)
        g = Geod(ellps='WGS84')
        self.site_lon, self.site_lat = g.fwd(ref_lon, ref_lat, az, dist)[0:2]
        return self.site_lat, self.site_lon

    def read_displacement(self, print_msg=True):
        # download file if it's not exists.
        if not os.path.isfile(self.file):
            self.dload_site()

        # read times, dis_e, dis_n, dis_u
        if print_msg:
            print('reading time and displacement in east/north/vertical direction')
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)
        self.times = np.array([dt(*time.strptime(i, "%y%b%d")[0:5])
                               for i in data[:, 1]])
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = data[:, (8,10,12,14,15,16)].astype(np.float32).T

        return self.times, self.dis_e, self.dis_n, self.dis_u, self.std_e, self.std_n, self.std_u





