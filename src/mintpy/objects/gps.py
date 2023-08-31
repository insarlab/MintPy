"""Class / utilities for GPS download / operations."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jul 2018                           #
############################################################
# Recommend import:
#     from mintpy.objects.gps import GPS


import csv
import datetime as dt
import glob
import os
from urllib.request import urlretrieve

import numpy as np
from pyproj import Geod

from mintpy.objects.coord import coordinate
from mintpy.utils import ptime, readfile, time_func, utils1 as ut

UNR_SITE_LIST_FILE = 'http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt'


def dload_site_list(out_file=None, url=UNR_SITE_LIST_FILE, print_msg=True):
    """download DataHoldings.txt.

    Parameters: out_file - str, path to the local file.
    Returns:    out_file - str, path to the local file.
                           default: './GPS/DataHoldings.txt'
    """
    # output file path
    if not out_file:
        out_dir = os.path.abspath('./GPS')
        out_file = os.path.join(out_dir, os.path.basename(url))
    else:
        out_dir = os.path.abspath(os.path.dirname(out_file))

    # output file directory
    if not os.path.isdir(out_dir):
        if print_msg:
            print('create directory:', out_dir)
        os.makedirs(out_dir)

    if print_msg:
        print(f'downloading site list from UNR Geod Lab: {url} to {out_file}')
    urlretrieve(url, out_file)
    return out_file


def search_gps(SNWE, start_date=None, end_date=None, site_list_file=None, min_num_solution=50, print_msg=True):
    """Search available GPS sites within the geo bounding box from UNR website.
    Parameters: SNWE             - tuple of 4 float, indicating (South, North, West, East) in degrees
                start_date       - str in YYYYMMDD format
                end_date         - str in YYYYMMDD format
                site_list_file   - str.
                min_num_solution - int, minimum number of solutions available
    Returns:    site_names       - 1D np.array of string for GPS station names
                site_lats        - 1D np.array for lat
                site_lons        - 1D np.array for lon
    """
    # download site list file if it's not found in current directory
    if site_list_file is None:
        site_list_file = os.path.basename(UNR_SITE_LIST_FILE)

    if not os.path.isfile(site_list_file):
        dload_site_list(site_list_file, print_msg=print_msg)

    txt_data = np.loadtxt(site_list_file,
                          dtype=bytes,
                          skiprows=1,
                          usecols=(0,1,2,3,4,5,6,7,8,9,10)).astype(str)
    site_names = txt_data[:, 0]
    site_lats, site_lons = txt_data[:, 1:3].astype(np.float32).T
    site_lons -= np.round(site_lons / (360.)) * 360.
    t0s = np.array([dt.datetime.strptime(i, "%Y-%m-%d") for i in txt_data[:, 7].astype(str)])
    t1s = np.array([dt.datetime.strptime(i, "%Y-%m-%d") for i in txt_data[:, 8].astype(str)])
    num_solution = txt_data[:, 10].astype(np.int16)

    # limit in space
    idx = ((site_lats >= SNWE[0]) * (site_lats <= SNWE[1]) *
           (site_lons >= SNWE[2]) * (site_lons <= SNWE[3]))

    # limit in time
    t0 = ptime.date_list2vector([start_date])[0][0] if start_date else None
    t1 = ptime.date_list2vector([end_date])[0][0] if end_date else None
    if start_date:
        idx *= t1s >= t0
    if end_date:
        idx *= t0s <= t1

    # limit on number of solutions
    if min_num_solution is not None:
        idx *= num_solution >= min_num_solution

    return site_names[idx], site_lats[idx], site_lons[idx]


def get_baseline_change(dates1, pos_x1, pos_y1, pos_z1,
                        dates2, pos_x2, pos_y2, pos_z2):
    """Calculate the baseline change between two GPS displacement time-series.
    Parameters: dates1/2     - 1D np.array of datetime.datetime object
                pos_x/y/z1/2 - 1D np.ndarray of displacement in meters in np.float32
    Returns:    dates        - 1D np.array of datetime.datetime object for the common dates
                bases        - 1D np.ndarray of displacement in meters in np.float32 for the common dates
    """
    dates = np.array(sorted(list(set(dates1) & set(dates2))))
    bases = np.zeros(dates.shape, dtype=np.float64)
    for i, date_str in enumerate(dates):
        idx1 = np.where(dates1 == date_str)[0][0]
        idx2 = np.where(dates2 == date_str)[0][0]
        basei = ((pos_x1[idx1] - pos_x2[idx2]) ** 2
               + (pos_y1[idx1] - pos_y2[idx2]) ** 2
               + (pos_z1[idx1] - pos_z2[idx2]) ** 2) ** 0.5
        bases[i] = basei
    bases -= bases[0]
    bases = np.array(bases, dtype=np.float32)
    return dates, bases


def get_gps_los_obs(meta, obs_type, site_names, start_date, end_date, gps_comp='enu2los',
                    horz_az_angle=-90., model=None, print_msg=True, redo=False):
    """Get the GPS LOS observations given the query info.

    Parameters: meta          - dict, dictionary of metadata of the InSAR file
                obs_type      - str, GPS observation data type, displacement or velocity.
                site_names    - list of str, GPS sites, output of search_gps()
                start_date    - str, date in YYYYMMDD format
                end_date      - str, date in YYYYMMDD format
                gps_comp      - str, flag of projecting 2/3D GPS into LOS
                                e.g. enu2los, hz2los, up2los
                horz_az_angle - float, azimuth angle of the horizontal motion in degree
                                measured from the north with anti-clockwise as positive
                model         - dict, time function model, e.g. {'polynomial': 1, 'periodic': [1.0, 0.5]}
                print_msg     - bool, print verbose info
                redo          - bool, ignore existing CSV file and re-calculate
    Returns:    site_obs      - 1D np.ndarray(), GPS LOS velocity or displacement in m or m/yr
    Examples:   from mintpy.objects import gps
                from mintpy.utils import readfile, utils as ut
                meta = readfile.read_attribute('geo/geo_velocity.h5')
                SNWE = ut.four_corners(meta)
                site_names = gps.search_gps(SNWE, start_date='20150101', end_date='20190619')
                vel = gps.get_gps_los_obs(meta, 'velocity',     site_names, start_date='20150101', end_date='20190619')
                dis = gps.get_gps_los_obs(meta, 'displacement', site_names, start_date='20150101', end_date='20190619')
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    num_site = len(site_names)

    # obs_type --> obs_ind
    obs_types = ['displacement', 'velocity']
    if obs_type not in obs_types:
        raise ValueError(f'un-supported obs_type: {obs_type}')
    obs_ind = 3 if obs_type.lower() == 'displacement' else 4

    # GPS CSV file info
    file_dir = os.path.dirname(meta['FILE_PATH'])
    csv_file = os.path.join(file_dir, f'gps_{gps_comp}')
    csv_file += f'{horz_az_angle:.0f}' if gps_comp == 'horz' else ''
    csv_file += '.csv'
    col_names = ['Site', 'Lon', 'Lat', 'Displacement', 'Velocity']
    col_types = ['U10'] + ['f8'] * (len(col_names) - 1)
    vprint(f'default GPS observation file name: {csv_file}')

    # skip re-calculate GPS if:
    # 1. redo is False AND
    # 2. csv_file exists (equivalent to num_row > 0) AND
    # 3. num_row >= num_site
    num_row = 0
    if os.path.isfile(csv_file):
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        num_row = fc.size

    if not redo and os.path.isfile(csv_file) and num_row >= num_site:
        # read from existing CSV file
        vprint(f'read GPS observations from file: {csv_file}')
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        site_obs = fc[col_names[obs_ind]]

        # get obs for the input site names only
        # in case the site_names are not consistent with the CSV file.
        if num_row != num_site:
            temp_names = fc[col_names[0]]
            temp_obs = np.array(site_obs, dtype=np.float32)
            site_obs = np.zeros(num_site, dtype=np.float32) * np.nan
            for i, site_name in enumerate(site_names):
                if site_name in temp_names:
                    site_obs[i] = temp_obs[temp_names == site_name][0]

    else:
        # calculate and save to CSV file
        data_list = []
        vprint('calculating GPS observation ...')

        # get geom_obj (meta / geom_file)
        geom_file = ut.get_geometry_file(['incidenceAngle','azimuthAngle'], work_dir=file_dir, coord='geo')
        if geom_file:
            geom_obj = geom_file
            vprint(f'use incidence / azimuth angle from file: {os.path.basename(geom_file)}')
        else:
            geom_obj = meta
            vprint('use incidence / azimuth angle from metadata')

        # loop for calculation
        prog_bar = ptime.progressBar(maxValue=num_site, print_msg=print_msg)
        for i, site_name in enumerate(site_names):
            prog_bar.update(i+1, suffix=f'{i+1}/{num_site} {site_name}')

            # calculate gps data value
            obj = GPS(site_name)
            vel, dis_ts = obj.get_gps_los_velocity(
                geom_obj,
                start_date=start_date,
                end_date=end_date,
                gps_comp=gps_comp,
                horz_az_angle=horz_az_angle,
                model=model)

            # ignore time-series if the estimated velocity is nan
            dis = np.nan if np.isnan(vel) else dis_ts[-1] - dis_ts[0]

            # save data to list
            data_list.append([obj.site, obj.site_lon, obj.site_lat, dis, vel])
        prog_bar.close()

        # # discard invalid sites
        # flag = np.isnan([x[-1] for x in data_list])
        # vprint('discard extra {} stations due to limited overlap/observations in time:'.format(np.sum(flag)))
        # vprint('  {}'.format(np.array(data_list)[flag][:,0].tolist()))
        # data_list = [x for x in data_list if not np.isnan(x[-1])]

        # write to CSV file
        vprint(f'write GPS observations to file: {csv_file}')
        with open(csv_file, 'w') as fc:
            fcw = csv.writer(fc)
            fcw.writerow(col_names)
            fcw.writerows(data_list)

        # prepare API output
        site_obs = np.array([x[obs_ind] for x in data_list])

    return site_obs




#################################### Beginning of GPS-GSI utility functions ########################
def read_pos_file(fname):
    import codecs
    fcp = codecs.open(fname, encoding = 'cp1252')
    fc = np.loadtxt(fcp, skiprows=20, dtype=str, comments=('*','-DATA'))

    ys = fc[:,0].astype(int)
    ms = fc[:,1].astype(int)
    ds = fc[:,2].astype(int)
    dates = [dt.datetime(year=y, month=m, day=d) for y,m,d in zip(ys, ms, ds)]

    X = fc[:,4].astype(np.float64).tolist()
    Y = fc[:,5].astype(np.float64).tolist()
    Z = fc[:,6].astype(np.float64).tolist()
    return dates, X, Y, Z


def get_pos_years(gps_dir, site):
    fnames = glob.glob(os.path.join(gps_dir, f'{site}.*.pos'))
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
        fname = os.path.join(gps_dir, f'{site}.{yeari[2:]}.pos')
        datesi, Xi, Yi, Zi = read_pos_file(fname)
        dates += datesi
        X += Xi
        Y += Yi
        Z += Zi
    dates = np.array(dates)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    date0 = dt.datetime.strptime(start_date, "%Y%m%d")
    date1 = dt.datetime.strptime(end_date, "%Y%m%d")
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
        self.data_dir = os.path.abspath(data_dir)
        self.version = version
        self.source = 'Nevada Geodetic Lab'

        # time-series data from Nevada Geodetic Lab
        # example link: http://geodesy.unr.edu/gps_timeseries/tenv3/IGS08/1LSU.IGS08.tenv3
        #               http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/CASU.tenv3
        if version == 'IGS08':
            self.file = os.path.join(data_dir, f'{site}.{version}.tenv3')
        elif version == 'IGS14':
            self.file = os.path.join(data_dir, f'{site}.tenv3')
        else:
            raise ValueError(f'un-recognized GPS data version: {version}')

        url_prefix = 'http://geodesy.unr.edu/gps_timeseries/tenv3'
        self.file_url = os.path.join(url_prefix, version, os.path.basename(self.file))

        # time-series plot from Nevada Geodetic Lab
        # example link: http://geodesy.unr.edu/tsplots/IGS08/TimeSeries/CAMO.png
        #               http://geodesy.unr.edu/tsplots/IGS14/IGS14/TimeSeries/CASU.png
        self.plot_file = os.path.join(data_dir, f'pic/{site}.png')

        url_prefix = 'http://geodesy.unr.edu/tsplots'
        if version == 'IGS08':
            url_prefix += f'/{version}'
        elif version == 'IGS14':
            url_prefix += '/{0}/{0}'.format(version)
        self.plot_file_url = os.path.join(url_prefix, f'TimeSeries/{site}.png')

        # list of stations from Nevada Geodetic Lab
        self.site_list_file = os.path.join(data_dir, 'DataHoldings.txt')
        if not os.path.isfile(self.site_list_file):
            dload_site_list(self.site_list_file)
        site_names = np.loadtxt(self.site_list_file, dtype=bytes, skiprows=1, usecols=(0)).astype(str)
        if site not in site_names:
            raise ValueError(f'Site {site} NOT found in file: {UNR_SITE_LIST_FILE}')

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
            print(f'downloading {self.site} from {self.file_url}')

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
        e0, e_off, n0, n_off = data[0, 7:11].astype(np.float32)
        e0 += e_off
        n0 += n_off

        az = np.arctan2(e0, n0) / np.pi * 180.
        dist = np.sqrt(e0**2 + n0**2)
        g = Geod(ellps='WGS84')
        self.site_lon, self.site_lat = g.fwd(ref_lon, ref_lat, az, dist)[0:2]
        return self.site_lat, self.site_lon

    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """ Read GPS displacement time-series (defined by start/end_date).
        Parameters: start/end_date - str in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in np.float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in np.float32
        """
        # download file if it's not exists.
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg:
            print('reading time and displacement in east/north/vertical direction')
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)

        self.dates = np.array([dt.datetime.strptime(i, "%y%b%d") for i in data[:, 1]])
        #self.dates = np.array([ptime.decimal_year2datetime(i) for i in data[:, 2]])
        self.date_list = [x.strftime('%Y%m%d') for x in self.dates]

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
            _, ax = plt.subplots(nrows=3, ncols=1, sharex=True)
            ax[0].scatter(self.dates, self.dis_e, s=2**2, label='East')
            ax[1].scatter(self.dates, self.dis_n, s=2**2, label='North')
            ax[2].scatter(self.dates, self.dis_u, s=2**2, label='Up')
            plt.show()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


    #####################################  Utility Functions ###################################
    def displacement_enu2los(self, inc_angle:float, az_angle:float, gps_comp='enu2los', horz_az_angle=-90.):
        """Convert displacement in ENU to LOS direction.

        Parameters: inc_angle     - float, LOS incidence angle in degree
                    az_angle      - float, LOS aziuth angle in degree
                                    from the north, defined as positive in clock-wise direction
                    gps_comp      - str, GPS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal to fault-parallel
                                    measured from the north with anti-clockwise as positive
        Returns:    dis_los       - 1D np.array for displacement in LOS direction
                    std_los       - 1D np.array for displacement standard deviation in LOS direction
        """
        # get unit vector for the component of interest
        unit_vec = ut.get_unit_vector4component_of_interest(
            los_inc_angle=inc_angle,
            los_az_angle=az_angle,
            comp=gps_comp.lower(),
            horz_az_angle=horz_az_angle,
        )

        # convert ENU to LOS direction
        self.dis_los = (  self.dis_e * unit_vec[0]
                        + self.dis_n * unit_vec[1]
                        + self.dis_u * unit_vec[2])
        # assuming ENU component are independent with each other
        self.std_los = (   (self.std_e * unit_vec[0])**2
                         + (self.std_n * unit_vec[1])**2
                         + (self.std_u * unit_vec[2])**2 ) ** 0.5

        return self.dis_los, self.std_los


    def get_los_geometry(self, geom_obj, print_msg=False):
        """Get the Line-of-Sight geometry info in incidence and azimuth angle in degrees."""
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
            az_angle  = readfile.read(geom_obj, datasetName='azimuthAngle',   box=box, print_msg=print_msg)[0][0,0]

        elif isinstance(geom_obj, dict):
            # use mean inc/az_angle from metadata
            inc_angle = ut.incidence_angle(geom_obj, dimension=0, print_msg=print_msg)
            az_angle  = ut.heading2azimuth_angle(float(geom_obj['HEADING']))

        else:
            raise ValueError(f'input geom_obj is neither str nor dict: {geom_obj}')

        return inc_angle, az_angle


    def read_gps_los_displacement(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                                  gps_comp='enu2los', horz_az_angle=-90., print_msg=False):
        """Read GPS displacement in LOS direction.

        Parameters: geom_obj      - dict / str, metadata of InSAR file, or geometry file path
                    start_date    - str in YYYYMMDD format
                    end_date      - str in YYYYMMDD format
                    ref_site      - str, reference GPS site
                    gps_comp      - str, GPS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal to fault-parallel
        Returns:    dates         - 1D np.array of datetime.datetime object
                    dis/std       - 1D np.array of displacement / uncertainty in meters
                    site_lalo     - tuple of 2 float, lat/lon of GPS site
                    ref_site_lalo - tuple of 2 float, lat/lon of reference GPS site
        """
        # read GPS object
        inc_angle, az_angle = self.get_los_geometry(geom_obj)
        dates = self.read_displacement(start_date, end_date, print_msg=print_msg)[0]
        dis, std = self.displacement_enu2los(inc_angle, az_angle, gps_comp=gps_comp, horz_az_angle=horz_az_angle)
        site_lalo = self.get_stat_lat_lon(print_msg=print_msg)

        # get LOS displacement relative to another GPS site
        if ref_site:
            ref_obj = GPS(site=ref_site, data_dir=self.data_dir)
            ref_obj.read_displacement(start_date, end_date, print_msg=print_msg)
            inc_angle, az_angle = ref_obj.get_los_geometry(geom_obj)
            ref_obj.displacement_enu2los(inc_angle, az_angle, gps_comp=gps_comp, horz_az_angle=horz_az_angle)
            ref_site_lalo = ref_obj.get_stat_lat_lon(print_msg=print_msg)

            # get relative LOS displacement on common dates
            dates = np.array(sorted(list(set(self.dates) & set(ref_obj.dates))))
            dis = np.zeros(dates.shape, np.float32)
            std = np.zeros(dates.shape, np.float32)
            for i, date_i in enumerate(dates):
                idx1 = np.where(self.dates == date_i)[0][0]
                idx2 = np.where(ref_obj.dates == date_i)[0][0]
                dis[i] = self.dis_los[idx1] - ref_obj.dis_los[idx2]
                std[i] = (self.std_los[idx1]**2 + ref_obj.std_los[idx2]**2)**0.5
        else:
            ref_site_lalo = None

        return dates, dis, std, site_lalo, ref_site_lalo


    def get_gps_los_velocity(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                             gps_comp='enu2los', horz_az_angle=-90., model=None):

        dates, dis = self.read_gps_los_displacement(
            geom_obj,
            start_date=start_date,
            end_date=end_date,
            ref_site=ref_site,
            gps_comp=gps_comp,
            horz_az_angle=horz_az_angle)[:2]

        # displacement -> velocity
        # if:
        # 1. num of observations > 2 AND
        # 2. time overlap > 1/4
        dis2vel = True
        if len(dates) <= 2:
            dis2vel = False
        elif start_date and end_date:
            t0 = ptime.date_list2vector([start_date])[0][0]
            t1 = ptime.date_list2vector([end_date])[0][0]
            if dates[-1] - dates[0] <= (t1 - t0) / 4:
                dis2vel = False

        if dis2vel:
            # specific time_func model
            date_list = [dt.datetime.strftime(i, '%Y%m%d') for i in dates]
            A = time_func.get_design_matrix4time_func(date_list, model=model)
            self.velocity = np.dot(np.linalg.pinv(A), dis)[1]
        else:
            self.velocity = np.nan

        return self.velocity, dis

#################################### End of GPS-UNR class ####################################
