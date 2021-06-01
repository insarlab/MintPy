############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jul 2018                           #
############################################################
# Utility scripts for GPS handling
# Recommend import:
#     from mintpy.objects.gps import GPS


import os
import csv
import codecs
from datetime import datetime as dt
import numpy as np
from pyproj import Geod
from urllib.request import urlretrieve

from mintpy.objects import timeseries
from mintpy.objects.coord import coordinate
from mintpy.utils import ptime, readfile, utils1 as ut


unr_site_list_file = 'http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt'


def dload_site_list(print_msg=True):
    """download DataHoldings.txt"""
    url = unr_site_list_file
    out_file = os.path.basename(url)
    if print_msg:
        print('downloading site list from UNR Geod Lab: {}'.format(url))
    urlretrieve(url, out_file)
    return out_file


def search_gps(SNWE, start_date=None, end_date=None, site_list_file=None, min_num_solution=50, print_msg=True):
    """Search available GPS sites within the geo bounding box from UNR website
    Parameters: SNWE       : tuple of 4 float, indicating (South, North, West, East) in degrees
                start_date : string in YYYYMMDD format
                end_date   : string in YYYYMMDD format
                site_list_file : string.
                min_num_solution : int, minimum number of solutions available
    Returns:    site_names : 1D np.array of string for GPS station names
                site_lats  : 1D np.array for lat
                site_lons  : 1D np.array for lon
    """
    # download site list file if it's not found in current directory
    if site_list_file is None:
        site_list_file = os.path.basename(unr_site_list_file)
        if not os.path.isfile(site_list_file):
            dload_site_list(print_msg=print_msg)

    txt_data = np.loadtxt(site_list_file,
                          dtype=bytes,
                          skiprows=1,
                          usecols=(0,1,2,3,4,5,6,7,8,9,10)).astype(str)
    site_names = txt_data[:, 0]
    site_lats, site_lons = txt_data[:, 1:3].astype(np.float32).T
    site_lons -= np.round(site_lons / (360.)) * 360.
    t_start = np.array([dt.strptime(i, "%Y-%m-%d") for i in txt_data[:, 7].astype(str)])
    t_end   = np.array([dt.strptime(i, "%Y-%m-%d") for i in txt_data[:, 8].astype(str)])
    num_solution = txt_data[:, 10].astype(np.int16)

    # limit on space
    idx = ((site_lats >= SNWE[0]) * (site_lats <= SNWE[1]) *
           (site_lons >= SNWE[2]) * (site_lons <= SNWE[3]))

    # limit on time
    if start_date:
        t0 = ptime.date_list2vector([start_date])[0][0]
        idx *= t_end >= t0
    if end_date:
        t1 = ptime.date_list2vector([end_date])[0][0]
        idx *= t_start <= t1

    # limit on number of solutions
    if min_num_solution is not None:
        idx *= num_solution >= min_num_solution

    return site_names[idx], site_lats[idx], site_lons[idx]


def get_baseline_change(dates1, pos_x1, pos_y1, pos_z1,
                        dates2, pos_x2, pos_y2, pos_z2):
    """Calculate the baseline change between two GPS displacement time-series
    Parameters: dates1/2     : 1D np.array of datetime.datetime object
                pos_x/y/z1/2 : 1D np.ndarray of displacement in meters in np.float32
    Returns:    dates        : 1D np.array of datetime.datetime object for the common dates
                bases        : 1D np.ndarray of displacement in meters in np.float32 for the common dates
    """
    dates = np.array(sorted(list(set(dates1) & set(dates2))))
    bases = np.zeros(dates.shape, dtype=np.float64)
    for i in range(len(dates)):
        idx1 = np.where(dates1 == dates[i])[0][0]
        idx2 = np.where(dates2 == dates[i])[0][0]
        basei = ((pos_x1[idx1] - pos_x2[idx2]) ** 2
               + (pos_y1[idx1] - pos_y2[idx2]) ** 2
               + (pos_z1[idx1] - pos_z2[idx2]) ** 2) ** 0.5
        bases[i] = basei
    bases -= bases[0]
    bases = np.array(bases, dtype=np.float32)
    return dates, bases


def get_gps_los_obs(insar_file, site_names, start_date, end_date,
                    gps_comp='enu2los', print_msg=True, redo=False):
    """Get the GPS LOS observations given the query info.

    Parameters: insar_file - str, InSAR LOS file, e.g. velocity or timeseries
                site_names - list of str, GPS sites, output of search_gps()
                start_date - str, date in YYYYMMDD format
                end_date   - str, date in YYYYMMDD format
                gps_comp   - str, flag of projecting 2/3D GPS into LOS
                             e.g. enu2los, hz2los, up2los
                print_msg  - bool, print verbose info
                redo       - bool, ignore existing CSV file and re-calculate
    Returns:    site_obs   - 1D np.ndarray(), GPS LOS velocity or displacement in m or m/yr
    """

    vprint = print if print_msg else lambda *args, **kwargs: None
    num_site = len(site_names)

    # basic info
    fdir = os.path.dirname(insar_file)
    meta = readfile.read_attribute(insar_file)
    obs_type = meta['FILE_TYPE']

    # GPS CSV file info
    csv_file = os.path.join(fdir, 'gps_{}.csv'.format(gps_comp))
    col_names = ['Site', 'Lon', 'Lat', 'LOS_displacement']
    col_names += ['LOS_velocity'] if obs_type == 'velocity' else []
    num_col = len(col_names)
    col_types = ['U10'] + ['f8'] * (num_col - 1)

    # skip re-calculate GPS if:
    #    redo is False AND
    #    csv_file exists (equivalent to num_row > 0) AND
    #    num_row >= num_site
    num_row = 0
    if os.path.isfile(csv_file):
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        num_row = fc.size

    if not redo and os.path.isfile(csv_file) and num_row >= num_site:
        # read from existing CSV file
        vprint('read GPS LOS observatioin from file: {}'.format(csv_file))
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        site_obs = fc[col_names[-1]]

    else:
        # calculate and save to CSV file
        data_list = []
        vprint('calculating GPS LOS observation ...')

        # get geom_obj (meta / geom_file)
        geom_file = ut.get_geometry_file(['incidenceAngle','azimuthAngle'], work_dir=fdir, coord='geo')
        if geom_file:
            geom_obj = geom_file
            vprint('use incidence / azimuth angle from file: {}'.format(os.path.basename(geom_file)))
        else:
            geom_obj = meta
            vprint('use incidence / azimuth angle from metadata')

        # loop for calculation
        prog_bar = ptime.progressBar(maxValue=num_site, print_msg=print_msg)
        for i, site_name in enumerate(site_names):
            prog_bar.update(i+1, suffix='{}/{} {}'.format(i+1, num_site, site_name))

            # calculate gps data value
            obj = GPS(site_name)

            if obs_type == 'velocity':
                vel, dis_ts = obj.get_gps_los_velocity(geom_obj,
                                                       start_date=start_date,
                                                       end_date=end_date,
                                                       gps_comp=gps_comp)
                data = [dis_ts[-1] - dis_ts[0], vel] if dis_ts.size > 2 else [np.nan, np.nan]

            elif obs_type == 'timeseries':
                dis_ts = obj.read_gps_los_displacement(geom_obj,
                                                       start_date=start_date,
                                                       end_date=end_date,
                                                       gps_comp=gps_comp)[1]
                data = [dis_ts[-1] - dis_ts[0]] if dis_ts.size > 2 else [np.nan]

            # save data to list
            data_list.append([obj.site, obj.site_lon, obj.site_lat] + data)
        prog_bar.close()
        site_obs = np.array([x[-1] for x in data_list])

        # write to CSV file
        vprint('write GPS LOS observations to file: {}'.format(csv_file))
        with open(csv_file, 'w') as fc:
            fcw = csv.writer(fc)
            fcw.writerow(col_names)
            fcw.writerows(data_list)

    return site_obs




#################################### Beginning of GPS-GSI utility functions ########################
def read_pos_file(fname):
    fcp = codecs.open(fname, encoding = 'cp1252')
    fc = np.loadtxt(fcp, skiprows=20, dtype=str, comments=('*','-DATA'))

    dates = [dt(year=y, month=m, day=d) for y,m,d in zip(fc[:,0].astype(int),
                                                         fc[:,1].astype(int),
                                                         fc[:,2].astype(int))]
    X = fc[:,4].astype(np.float64).tolist()
    Y = fc[:,5].astype(np.float64).tolist()
    Z = fc[:,6].astype(np.float64).tolist()
    return dates, X, Y, Z


def get_pos_years(gps_dir, site):
    fnames = glob.glob(os.path.join(gps_dir, '{}.*.pos'.format(site)))
    years = [os.path.basename(i).split('.')[1] for i in fnames]
    years = ptime.yy2yyyy(years)
    return years


def read_GSI_F3(gps_dir, site, start_date=None, end_date=None):
    year0 = int(start_date[0:4])
    year1 = int(end_date[0:4])
    num_year = year1 - year0 + 1

    dates, X, Y, Z = [], [], [], []
    for i in range(num_year):
        yeari = str(year0 + i)
        fname = os.path.join(gps_dir, '{}.{}.pos'.format(site, yeari[2:]))
        datesi, Xi, Yi, Zi = read_pos_file(fname)
        dates += datesi
        X += Xi
        Y += Yi
        Z += Zi
    dates = np.array(dates)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    date0 = dt.strptime(start_date, "%Y%m%d")
    date1 = dt.strptime(end_date, "%Y%m%d")
    flag = np.ones(X.shape, dtype=np.bool_)
    flag[dates < date0] = False
    flag[dates > date1] = False
    return dates[flag], X[flag], Y[flag], Z[flag]

#################################### End of GPS-GSI utility functions ##############################




#################################### Beginning of GPS-UNR class ####################################
class GPS:
    """GPS class for GPS time-series of daily solution

    Example:
      import matplotlib.pyplot as plt
      from mintpy.objects.gps import GPS
      from mintpy.utils import utils as ut
      gps_obj = GPS(site='GV05', data_dir='~/insarlab/GPS')
      gps_obj.open()
      dis_los = ut.enu2los(gps_obj.dis_e,
                           gps_obj.dis_n,
                           gps_obj.dis_u)
      dates = gps_obj.dates
      plt.figure()
      plt.scatter(dates, dis_los)
      plt.show()
    """

    def __init__(self, site, data_dir='./GPS', version='IGS14'):
        self.site = site
        self.data_dir = data_dir
        self.version = version
        self.source = 'Nevada Geodetic Lab'

        # time-series data from Nevada Geodetic Lab
        # example link: http://geodesy.unr.edu/gps_timeseries/tenv3/IGS08/1LSU.IGS08.tenv3
        #               http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/CASU.tenv3
        if version == 'IGS08':
            self.file = os.path.join(data_dir, '{s}.{v}.tenv3'.format(s=site, v=version))
        elif version == 'IGS14':
            self.file = os.path.join(data_dir, '{s}.tenv3'.format(s=site))
        else:
            raise ValueError('un-recognized GPS data version: {}'.format(version))

        url_prefix = 'http://geodesy.unr.edu/gps_timeseries/tenv3'
        self.file_url = os.path.join(url_prefix, version, os.path.basename(self.file))

        # time-series plot from Nevada Geodetic Lab
        # example link: http://geodesy.unr.edu/tsplots/IGS08/TimeSeries/CAMO.png
        #               http://geodesy.unr.edu/tsplots/IGS14/IGS14/TimeSeries/CASU.png
        self.plot_file = os.path.join(data_dir, 'pic/{}.png'.format(site))

        url_prefix = 'http://geodesy.unr.edu/tsplots'
        if version == 'IGS08':
            url_prefix += '/{0}'.format(version)
        elif version == 'IGS14':
            url_prefix += '/{0}/{0}'.format(version)
        self.plot_file_url = os.path.join(url_prefix, 'TimeSeries/{}.png'.format(site))

        # list of stations from Nevada Geodetic Lab
        self.site_list_file = os.path.join(data_dir, 'DataHoldings.txt')

        # directories for data files and plot files
        for fdir in [data_dir, os.path.dirname(self.plot_file)]:
            if not os.path.isdir(fdir):
                os.makedirs(fdir)

    def open(self, print_msg=True):
        if not os.path.isfile(self.file):
            self.dload_site()
        self.get_stat_lat_lon(print_msg=print_msg)
        self.read_displacement(print_msg=print_msg)

    def dload_site(self, print_msg=True):
        if print_msg:
            print('downloading {} from {}'.format(self.site, self.file_url))

        urlretrieve(self.file_url, self.file)
        urlretrieve(self.plot_file_url, self.plot_file)

        return self.file

    def get_stat_lat_lon(self, print_msg=True):
        """Get station lat/lon"""
        if print_msg:
            print('calculating station lat/lon')
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

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

    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """ Read GPS displacement time-series (defined by start/end_date)
        Parameters: start/end_date : str in YYYYMMDD format
        Returns:    dates : 1D np.ndarray of datetime.datetime object
                    dis_e/n/u : 1D np.ndarray of displacement in meters in np.float32
                    std_e/n/u : 1D np.ndarray of displacement STD in meters in np.float32
        """
        # download file if it's not exists.
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg:
            print('reading time and displacement in east/north/vertical direction')
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)

        self.dates = np.array([dt.strptime(i, "%y%b%d") for i in data[:, 1]])
        #self.dates = np.array([ptime.decimal_year2datetime(i) for i in data[:, 2]])

        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = data[:, (8,10,12,14,15,16)].astype(np.float32).T

        # cut out the specified time range
        t_flag = np.ones(len(self.dates), np.bool_)
        if start_date:
            t0 = ptime.date_list2vector([start_date])[0][0]
            t_flag[self.dates < t0] = 0
        if end_date:
            t1 = ptime.date_list2vector([end_date])[0][0]
            t_flag[self.dates > t1] = 0
        self.dates = self.dates[t_flag]
        self.dis_e = self.dis_e[t_flag]
        self.dis_n = self.dis_n[t_flag]
        self.dis_u = self.dis_u[t_flag]
        self.std_e = self.std_e[t_flag]
        self.std_n = self.std_n[t_flag]
        self.std_u = self.std_u[t_flag]

        if display:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)
            ax[0].scatter(self.dates, self.dis_e, s=2**2, label='East')
            ax[1].scatter(self.dates, self.dis_n, s=2**2, label='North')
            ax[2].scatter(self.dates, self.dis_u, s=2**2, label='Up')
            plt.show()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


    #####################################  Utility Functions ###################################
    def displacement_enu2los(self, inc_angle:float, head_angle:float, gps_comp='enu2los'):
        """Convert displacement in ENU to LOS direction
        Parameters: inc_angle  : float, local incidence angle in degree
                    head_angle : float, satellite orbit heading direction in degree
                        from the north, defined as positive in clock-wise direction
                    gps_comp   : string, GPS components used to convert to LOS direction
        Returns:    dis_los : 1D np.array for displacement in LOS direction
                    std_los : 1D np.array for displacement standard deviation in LOS direction
        """
        # get LOS unit vector
        inc_angle *= np.pi/180.
        head_angle *= np.pi/180.
        unit_vec = [np.sin(inc_angle) * np.cos(head_angle) * -1,
                    np.sin(inc_angle) * np.sin(head_angle),
                    np.cos(inc_angle)]

        gps_comp = gps_comp.lower()
        if gps_comp in ['enu2los']:
            pass
        elif gps_comp in ['en2los', 'hz2los']:
            unit_vec[2] = 0.
        elif gps_comp in ['u2los', 'up2los']:
            unit_vec[0] = 0.
            unit_vec[1] = 0.
        else:
            raise ValueError('Un-known input gps components:'+str(gps_comp))

        # convert ENU to LOS direction
        self.dis_los = (self.dis_e * unit_vec[0]
                        + self.dis_n * unit_vec[1]
                        + self.dis_u * unit_vec[2])

        # assuming ENU component are independent with each other
        self.std_los = ((self.std_e * unit_vec[0])**2
                         + (self.std_n * unit_vec[1])**2
                         + (self.std_u * unit_vec[2])**2)**0.5
        return self.dis_los, self.std_los


    def get_los_geometry(self, geom_obj, print_msg=False):
        lat, lon = self.get_stat_lat_lon(print_msg=print_msg)

        # get LOS geometry
        if isinstance(geom_obj, str):
            # geometry file
            atr = readfile.read_attribute(geom_obj)
            coord = coordinate(atr, lookup_file=geom_obj)
            y, x = coord.geo2radar(lat, lon, print_msg=print_msg)[0:2]
            # check against image boundary
            y = max(0, y);  y = min(int(atr['LENGTH'])-1, y)
            x = max(0, x);  x = min(int(atr['WIDTH'])-1, x)
            box = (x, y, x+1, y+1)
            inc_angle = readfile.read(geom_obj, datasetName='incidenceAngle', box=box, print_msg=print_msg)[0][0,0]
            az_angle  = readfile.read(geom_obj, datasetName='azimuthAngle', box=box, print_msg=print_msg)[0][0,0]
            head_angle = ut.azimuth2heading_angle(az_angle)

        elif isinstance(geom_obj, dict):
            # use mean inc/head_angle from metadata
            inc_angle = ut.incidence_angle(geom_obj, dimension=0, print_msg=print_msg)
            head_angle = float(geom_obj['HEADING'])
            # for old reading of los.rdr band2 data into headingAngle directly
            if (head_angle + 180.) > 45.:
                head_angle = ut.azimuth2heading_angle(head_angle)

        else:
            raise ValueError('input geom_obj is neight str nor dict: {}'.format(geom_obj))

        return inc_angle, head_angle


    def read_gps_los_displacement(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                                  gps_comp:str='enu2los', print_msg=False):
        """Read GPS displacement in LOS direction
        Parameters: geom_obj : dict / str, metadata of InSAR file, or geometry file path
                    start_date : string in YYYYMMDD format
                    end_date   : string in YYYYMMDD format
                    ref_site   : string, reference GPS site
                    gps_comp   : string, GPS components used to convert to LOS direction
        Returns:    dates : 1D np.array of datetime.datetime object
                    dis   : 1D np.array of displacement in meters
                    std   : 1D np.array of displacement uncertainty in meters
                    site_lalo : tuple of 2 float, lat/lon of GPS site
                    ref_site_lalo : tuple of 2 float, lat/lon of reference GPS site
        """
        # read GPS object
        inc_angle, head_angle = self.get_los_geometry(geom_obj)
        dates = self.read_displacement(start_date, end_date, print_msg=print_msg)[0]
        dis, std = self.displacement_enu2los(inc_angle, head_angle, gps_comp=gps_comp)
        site_lalo = self.get_stat_lat_lon(print_msg=print_msg)

        # get LOS displacement relative to another GPS site
        if ref_site:
            ref_obj = GPS(site=ref_site, data_dir=self.data_dir)
            ref_obj.read_displacement(start_date, end_date, print_msg=print_msg)
            inc_angle, head_angle = ref_obj.get_los_geometry(geom_obj)
            ref_obj.displacement_enu2los(inc_angle, head_angle, gps_comp=gps_comp)
            ref_site_lalo = ref_obj.get_stat_lat_lon(print_msg=print_msg)

            # get relative LOS displacement on common dates
            dates = np.array(sorted(list(set(self.dates) & set(ref_obj.dates))))
            dis = np.zeros(dates.shape, np.float32)
            std = np.zeros(dates.shape, np.float32)
            for i in range(len(dates)):
                idx1 = np.where(self.dates == dates[i])[0][0]
                idx2 = np.where(ref_obj.dates == dates[i])[0][0]
                dis[i] = self.dis_los[idx1] - ref_obj.dis_los[idx2]
                std[i] = (self.std_los[idx1]**2 + ref_obj.std_los[idx2]**2)**0.5
        else:
            ref_site_lalo = None

        return dates, dis, std, site_lalo, ref_site_lalo


    def get_gps_los_velocity(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                             gps_comp='enu2los'):
        dates, dis = self.read_gps_los_displacement(geom_obj,
                                                    start_date=start_date,
                                                    end_date=end_date,
                                                    ref_site=ref_site,
                                                    gps_comp=gps_comp)[:2]

        # displacement -> velocity
        date_list = [dt.strftime(i, '%Y%m%d') for i in dates]
        if len(date_list) > 2:
            A = timeseries.get_design_matrix4time_func(date_list)
            self.velocity = np.dot(np.linalg.pinv(A), dis)[1]
        else:
            self.velocity = np.nan

        return self.velocity, dis

#################################### End of GPS-UNR class ####################################
