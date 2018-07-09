############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, Jul 2018                          #
############################################################
# class used for file operation within PySAR
# Recommend import:
#     from pysar.objects import gps

import os
import time
from datetime import datetime as dt
import numpy as np
from pyproj import Geod
from urllib.request import urlretrieve


class gps:
    """GPS class for GPS time-series of daily solution
    
    Check all GPS sites info on:
    http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt

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

    def open(self, print_msg=True):
        if not os.path.isfile(self.file):
            self.download()
        self.get_stat_lat_lon(print_msg=print_msg)
        self.read_displacement(print_msg=print_msg)
    
    def download(self, print_msg=True):
        # get url from Nevada Geodetic Lab
        url = 'http://geodesy.unr.edu/gps_timeseries'
        if self.data_format == 'env3':
            url = os.path.join(url, 'tenv3')
        elif self.data_format == 'xyz2':
            url = os.path.join(url, 'txyz')
        url = os.path.join(url, 'IGS08/{}'.format(os.path.basename(self.file)))

        if print_msg:
            print('downloading {} from {}'.format(self.site, url))
        urlretrieve(url, self.file)
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
            self.download()

        # read times, dis_e, dis_n, dis_u
        if print_msg:
            print('reading time and displacement in east/north/vertical direction')
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)
        self.times = [dt(*time.strptime(i, "%y%b%d")[0:5]) for i in data[:, 1]]
        (self.dis_e,
         self.dis_n,
         self.dis_u) = data[:, (8,10,12)].astype(np.float32).T
        return self.times, self.dis_e, self.dis_n, self.dis_u







