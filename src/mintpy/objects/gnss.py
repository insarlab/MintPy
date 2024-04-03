############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jul 2018                           #
############################################################
# Utility scripts for GNSS handling
# Recommend import:
#     from mintpy.objects import gnss


import csv
import datetime as dt
import glob
import os
import zipfile
from urllib.request import urlopen, urlretrieve

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Geod

from mintpy.objects.coord import coordinate
from mintpy.utils import ptime, readfile, time_func, utils1 as ut

supported_sources = ['UNR', 'ESESES', 'JPL-SIDESHOW', 'Generic']


def dload_site_list(out_file=None, source='UNR', print_msg=True) -> str:
    """Download single file with list of GNSS site locations.
    """
    # check source is supported
    assert source in supported_sources, \
        f'Source {source:s} not supported. Use one of {supported_sources}'

    # determine URL
    if source == 'UNR':
        UNR_site_list_file_url = 'http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt'
        site_list_file_url = UNR_site_list_file_url
    elif source == 'ESESES':
        ESESES_site_list_file_url = 'http://garner.ucsd.edu/pub/measuresESESES_products/Velocities/ESESES_Velocities.txt'
        site_list_file_url = ESESES_site_list_file_url
    elif source == 'JPL-SIDESHOW':
        JPL_SIDESHOW_site_list_file_url = 'https://sideshow.jpl.nasa.gov/post/tables/table2.html'
        site_list_file_url = JPL_SIDESHOW_site_list_file_url

    # handle output file
    if out_file is None:
        out_file = os.path.basename(site_list_file_url)

    # download file
    if not os.path.exists(out_file):
        if print_msg:
            print(f'Downloading site list from {source:s}: {site_list_file_url:s} to {out_file:s}')
        urlretrieve(site_list_file_url, out_file)  #nosec

    return out_file


def search_gnss(SNWE, source='UNR', start_date=None, end_date=None,
               site_list_file=None, min_num_solution=None, print_msg=True):
    """Search available GNSS sites within the geo bounding box from UNR website
    Parameters: SNWE             - tuple of 4 float, indicating (South, North, West, East) in degrees
                source           - str, program or institution that processed the GNSS data
                start_date       - str, date in YYYYMMDD format
                end_date         - str, date in YYYYMMDD format
                site_list_file   - str
                min_num_solution - int, minimum number of solutions available
    Returns:    site_names       - 1D np.array of string, GNSS station names
                site_lats        - 1D np.array, lat
                site_lons        - 1D np.array, lon
    """
    # check start and end dates if provided
    if start_date is not None:
        start_date = dt.datetime.strptime(start_date, '%Y%m%d')
    if end_date is not None:
        end_date = dt.datetime.strptime(end_date, '%Y%m%d')
    if start_date is not None and end_date is not None:
        assert(start_date < end_date), 'Start date must be before end date'

    # check file name
    if site_list_file is None:
        if source == 'Generic':
            raise ValueError('Site list file must be specified for generic inputs')
        else:
            site_list_file = dload_site_list(source=source, print_msg=print_msg)

    # check whether site list file is in current directory
    if not os.path.isfile(site_list_file):
        # Download file
        dload_site_list(site_list_file, source=source, print_msg=print_msg)

    # parse data from file
    if source == 'UNR':
        site_data = read_UNR_station_list(site_list_file)
    elif source == 'ESESES':
        site_data = read_ESESES_station_list(site_list_file)
    elif source == 'JPL-SIDESHOW':
        site_data = read_JPL_SIDESHOW_station_list(site_list_file)
    elif source == 'Generic':
        site_data = read_Generic_station_list(site_list_file)

    if print_msg == True:
        print('Loaded data for fields: {:s}'.\
            format(' '.join(list(site_data.columns))))

    # ensure that station name is consistent
    site_data['site'] = [site_data.iloc[i,:].site.upper() for i in range(site_data.shape[0])]

    # parse bounding box
    lat_min, lat_max, lon_min, lon_max = SNWE
    assert (lon_min < lon_max) and (lat_min < lat_max), 'Check bounding box'

    if print_msg == True:
        print('Cropping to')
        print(f'lon range: {lon_min:.5f} to {lon_max:.5f}')
        print(f'lat range: {lat_min:.5f} to {lat_max:.5f}')

    # Ensure lon values in (-180, 180]
    site_data['lon'] = [lon - 360 if lon > 180 else lon for lon in site_data['lon']]

    # limit in space
    drop_ndx = (site_data.lat < lat_min) \
                | (site_data.lat > lat_max) \
                | (site_data.lon < lon_min) \
                | (site_data.lon > lon_max)
    site_data.drop(site_data[drop_ndx].index, inplace=True)

    # limit in time
    if start_date is not None:
        if hasattr(site_data, 'start_date'):
            drop_ndx = site_data.start_date > start_date
            site_data.drop(site_data[drop_ndx].index, inplace=True)
        else:
            print('No date information available--date range not applied to GNSS site selection')

    if end_date is not None:
        if hasattr(site_data, 'end_date'):
            drop_ndx = site_data.end_date < end_date
            site_data.drop(site_data[drop_ndx].index, inplace=True)
        else:
            print('No date information available--date range not applied to GNSS site selection')

    # limit based on number of solutions
    if hasattr(site_data, 'num_solution'):
        drop_ndx = site_data.num_solution < min_num_solution
        site_data.drop(site_data[drop_ndx].index, inplace=True)

    # final reporting
    if print_msg == True:
        print(f'{site_data.shape[0]:d} stations available')

    return (site_data.site.to_numpy(),
            site_data.lat.to_numpy(),
            site_data.lon.to_numpy())

def read_UNR_station_list(site_list_file:str, print_msg=True) -> pd.DataFrame:
    """Return names and lon/lat values for UNR GNSS stations.
    """
    if print_msg == True:
        print('Parsing UNR site list file')

    # Read file contents
    site_data = pd.read_fwf(site_list_file,
                            widths=(4, 9, 12, 9, 14, 14, 14, 11, 11, 11, 7))

    # Rename columns for uniformity
    site_data.rename(columns={'Sta': 'site',
        'Lat(deg)': 'lat', 'Long(deg)': 'lon',
        'Dtbeg': 'start_date', 'Dtend': 'end_date',
        'NumSol': 'num_solution'}, inplace=True)

    # Format dates
    site_data['start_date'] = [dt.datetime.strptime(date, '%Y-%m-%d') \
                                for date in site_data.start_date]
    site_data['end_date'] = [dt.datetime.strptime(date, '%Y-%m-%d') \
                                for date in site_data.end_date]

    return site_data

def read_ESESES_station_list(site_list_file:str, print_msg=True) -> pd.DataFrame:
    """Return names and lon/lat values for JPL GNSS stations.
    """
    if print_msg == True:
        print('Parsing ESESES site list file')

    # Read file contents
    site_data = pd.read_csv(site_list_file, header = 14, sep=r'\s+')

    # Rename columns for uniformity
    site_data.rename(columns={'Site': 'site',
        'Latitude': 'lat', 'Longitude': 'lon'}, inplace=True)

    return site_data

def read_JPL_SIDESHOW_station_list(site_list_file:str, print_msg=True) -> pd.DataFrame:
    """Return names and lon/lat values for JPL-SIDESHOW GNSS stations.
    """
    if print_msg == True:
        print('Parsing JPL-SIDESHOW site list file')

    # read file contents
    with open(site_list_file, 'r') as site_list:
        lines = site_list.readlines()

    # find lines containing position and velocity data
    line_len = 8  # number of entries in a data line
    name_len = 4  # number of letters in a station name
    data_lines = [line.strip('\n') for line in lines if (len(line.split()) == line_len) \
                    and (len(line.split()[0]) == name_len)]

    # transform format from (POS \n VEL) to (POS VEL)
    pos_lines = data_lines[0::2]
    vel_lines = data_lines[1::2]
    combo_lines = list(zip(pos_lines, vel_lines))

    # empty lists
    sites = []
    lats = []
    lons = []
    elevs = []
    Nvels = []
    Evels = []
    Uvels = []
    Nerrs = []
    Eerrs = []
    Uerrs = []

    # parse data
    for line in combo_lines:
        pos_info, vel_info = line

        # parse line values
        site,    _,  lat,  lon, elev,    _,    _,    _ = pos_info.split()
        _   ,    _, Nvel, Evel, Uvel, Nerr, Eerr, Uerr = vel_info.split()

        # format data
        sites.append(site)
        lats.append(float(lat))
        lons.append(float(lon))
        elevs.append(float(elev))
        Nvels.append(float(Nvel))
        Evels.append(float(Evel))
        Uvels.append(float(Uvel))
        Nerrs.append(float(Nerr))
        Eerrs.append(float(Eerr))
        Uerrs.append(float(Uerr))

    # format data frame
    data = {'site': sites,
            'lat': lats,  'lon': lons, 'elev': elevs,
            'vel_n': Nvels, 'vel_e': Evels, 'vel_u': Uvels,
            'err_n': Nerrs, 'err_e': Eerrs, 'err_u': Uerrs}
    site_data = pd.DataFrame(data)

    return site_data

def read_Generic_station_list(site_list_file:str, print_msg=True) -> pd.DataFrame:
    """Return names and lon/lat values for GNSS stations processed by an
    otherwise-unsupported source.

    The user must format the station position data in a file named
    GenericList.txt The file should have three, nine, or eleven space-
    separated columns:

     site lat lon [vel_e vel_n vel_u err_e err_n err_u] [start_date end_date]

    where site is the four-digit, alphanumeric (uppercase) site code; and
    lat/lon are in decimal degrees. If included, vel should be in units of
    m/yr; and dates should be in format YYYYMMDD.
    """
    if print_msg == True:
        print('Parsing JPL-SIDESHOW site list file')

    # read file contents
    site_data = pd.read_csv(site_list_file, delimiter=' ', names=('site', 'lat', 'lon'))

    return site_data


def get_baseline_change(dates1, pos_x1, pos_y1, pos_z1,
                        dates2, pos_x2, pos_y2, pos_z2):
    """Calculate the baseline change between two GNSS displacement time-series
    Parameters: dates1/2     - 1D np.array, datetime.datetime object
                pos_x/y/z1/2 - 1D np.ndarray, displacement in meters in float32
    Returns:    dates        - 1D np.array, datetime.datetime object for the
                               common dates
                bases        - 1D np.ndarray, displacement in meters in float32
                               for the common dates
    """
    dates = np.array(sorted(list(set(dates1) & set(dates2))))
    bases = np.zeros(dates.shape, dtype=float)
    for i, date in enumerate(dates):
        idx1 = np.where(dates1 == date)[0][0]
        idx2 = np.where(dates2 == date)[0][0]
        basei = ((pos_x1[idx1] - pos_x2[idx2]) ** 2
               + (pos_y1[idx1] - pos_y2[idx2]) ** 2
               + (pos_z1[idx1] - pos_z2[idx2]) ** 2) ** 0.5
        bases[i] = basei
    bases -= bases[0]
    bases = np.array(bases, dtype=float)

    return dates, bases


def get_gnss_los_obs(meta, obs_type, site_names, start_date, end_date, source='UNR',
                    gnss_comp='enu2los', horz_az_angle=-90., model=None,
                    print_msg=True, redo=False):
    """Get the GNSS LOS observations given the query info.

    Parameters: meta       - dict, dictionary of metadata of the InSAR file
                obs_type   - str, GNSS observation data type, displacement or velocity.
                site_names - list of str, GNSS sites, output of search_gnss()
                start_date - str, date in YYYYMMDD format
                end_date   - str, date in YYYYMMDD format
                source     - str, program or institution that processed the GNSS data
                gnss_comp  - str, flag of projecting 2/3D GNSS into LOS
                             e.g. enu2los, hz2los, up2los
                horz_az_angle - float, azimuth angle of the horizontal motion in degree
                             measured from the north with anti-clockwise as positive
                model         - dict, time function model, e.g. {'polynomial': 1, 'periodic': [1.0, 0.5]}
                print_msg  - bool, print verbose info
                redo       - bool, ignore existing CSV file and re-calculate
    Returns:    site_obs   - 1D np.ndarray(), GNSS LOS velocity or displacement in m or m/yr
    Examples:   from mintpy.objects import gnss
                from mintpy.utils import readfile, utils as ut
                meta = readfile.read_attribute('geo/geo_velocity.h5')
                SNWE = ut.four_corners(meta)
                site_names = gnss.search_gnss(SNWE, start_date='20150101', end_date='20190619')
                vel = gnss.get_gnss_los_obs(meta, 'velocity',     site_names, start_date='20150101', end_date='20190619')
                dis = gnss.get_gnss_los_obs(meta, 'displacement', site_names, start_date='20150101', end_date='20190619')
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    num_site = len(site_names)

    # obs_type --> obs_ind
    obs_types = ['displacement', 'velocity']
    if obs_type not in obs_types:
        raise ValueError(f'un-supported obs_type: {obs_type}')
    obs_ind = 3 if obs_type.lower() == 'displacement' else 4

    # GNSS CSV file info
    file_dir = os.path.dirname(meta['FILE_PATH'])
    csv_file = os.path.join(file_dir, f'gnss_{gnss_comp:s}')
    csv_file += f'{horz_az_angle:.0f}' if gnss_comp == 'horz' else ''
    csv_file += '.csv'
    col_names = ['Site', 'Lon', 'Lat', 'Displacement', 'Velocity']
    col_types = ['U10'] + ['f8'] * (len(col_names) - 1)
    vprint(f'default GNSS observation file name: {csv_file:s}')

    # skip re-calculate GNSS if:
    # 1. redo is False AND
    # 2. csv_file exists (equivalent to num_row > 0) AND
    # 3. num_row >= num_site
    num_row = 0
    if os.path.isfile(csv_file):
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        num_row = fc.size

    if not redo and os.path.isfile(csv_file) and num_row >= num_site:
        # read from existing CSV file
        vprint(f'read GNSS observations from file: {csv_file:s}')
        fc = np.genfromtxt(csv_file, dtype=col_types, delimiter=',', names=True)
        site_obs = fc[col_names[obs_ind]]

        # get obs for the input site names only
        # in case the site_names are not consistent with the CSV file.
        if num_row != num_site:
            temp_names = fc[col_names[0]]
            temp_obs = np.array(site_obs, dtype=float)
            site_obs = np.zeros(num_site, dtype=float) * np.nan
            for i, site_name in enumerate(site_names):
                if site_name in temp_names:
                    site_obs[i] = temp_obs[temp_names == site_name][0]

    else:
        # calculate and save to CSV file
        data_list = []
        vprint('calculating GNSS observation ...')

        # get geom_obj (meta / geom_file)
        geom_file = ut.get_geometry_file(['incidenceAngle','azimuthAngle'],
                                         work_dir=file_dir, coord='geo')
        if geom_file:
            geom_obj = geom_file
            vprint('use incidence / azimuth angle from file: {}'.\
                   format(os.path.basename(geom_file)))
        else:
            geom_obj = meta
            vprint('use incidence / azimuth angle from metadata')

        # loop for calculation
        prog_bar = ptime.progressBar(maxValue=num_site, print_msg=print_msg)
        for i, site_name in enumerate(site_names):
            prog_bar.update(i+1, suffix=f'{i+1}/{num_site} {site_name:s}')

            # calculate GNSS data value
            GNSSclass = GNSS.get_gnss_obj_by_source(source)
            obj = GNSSclass(site_name)
            obj.open(print_msg=print_msg)
            vel, dis_ts = obj.get_gnss_los_velocity(
                geom_obj,
                start_date=start_date,
                end_date=end_date,
                gnss_comp=gnss_comp,
                horz_az_angle=horz_az_angle,
                model=model)

            # ignore time-series if the estimated velocity is nan
            dis = np.nan if np.isnan(vel) else dis_ts[-1] - dis_ts[0]

            # save data to list
            data_list.append([obj.site, obj.site_lon, obj.site_lat, dis, vel])
        prog_bar.close()

        # write to CSV file
        vprint(f'write GNSS observations to file: {csv_file:s}')
        with open(csv_file, 'w') as fc:
            fcw = csv.writer(fc)
            fcw.writerow(col_names)
            fcw.writerows(data_list)

        # prepare API output
        site_obs = np.array([x[obs_ind] for x in data_list])

    return site_obs



#################################### Beginning of GNSS-GSI utility functions ########################
def read_pos_file(fname):
    import codecs
    fcp = codecs.open(fname, encoding = 'cp1252')
    fc = np.loadtxt(fcp, skiprows=20, dtype=str, comments=('*','-DATA'))

    ys = fc[:,0].astype(int)
    ms = fc[:,1].astype(int)
    ds = fc[:,2].astype(int)
    dates = [dt.datetime(year=y, month=m, day=d) for y,m,d in zip(ys, ms, ds)]

    X = fc[:,4].astype(float).tolist()
    Y = fc[:,5].astype(float).tolist()
    Z = fc[:,6].astype(float).tolist()

    return dates, X, Y, Z


def get_pos_years(gnss_dir, site):
    fnames = glob.glob(os.path.join(gnss_dir, f'{site:s}.*.pos'))
    years = [os.path.basename(i).split('.')[1] for i in fnames]
    years = ptime.yy2yyyy(years)
    return years


def read_GSI_F3(gnss_dir, site, start_date=None, end_date=None):
    year0 = int(start_date[0:4])
    year1 = int(end_date[0:4])
    num_year = year1 - year0 + 1

    dates, X, Y, Z = [], [], [], []
    for i in range(num_year):
        yeari = str(year0 + i)
        fname = os.path.join(gnss_dir, f'{site:d}.{yeari[2:]:s}.pos')
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


#################################### End of GNSS-GSI utility functions ##############################



#################################### Beginning of GNSS class ########################################
class GNSS:
    """GNSS class for time-series of daily solution.

    The GNSS class is solely meant to be a parent class. Child classes, defined
    below, support functions for downloading and parsing GNSS position based on
    the processing source (e.g., UNR, etc.). Use the `get_gnss_obj_by_source`
    method to determine appropriate child class.
    """

    def __init__(self, site: str, data_dir='./GNSS', version='IGS14'):
        # Record properties
        self.site = site
        self.version = version

        # create data directory if not exist
        self.data_dir = os.path.abspath(data_dir)
        if not os.path.exists(self.data_dir):
            os.mkdir(self.data_dir)

        # variables to be filled by child classes
        self.dates = None
        self.date_list = None
        self.dis_e = None
        self.dis_n = None
        self.dis_u = None
        self.std_e = None
        self.std_n = None
        self.std_u = None

        return None

    def open(self, file=None, print_msg=True):
        """Read the lat/lon and displacement data of the station.
        Download if necessary.
        """
        # download file if not present
        if not hasattr(self, 'file'):
            self.dload_site(print_msg=print_msg)

        # retrieve data from file
        self.get_stat_lat_lon(print_msg=print_msg)
        self.read_displacement(print_msg=print_msg)

        return None

    @staticmethod
    def get_gnss_obj_by_source(source:str):
        """Return the appropriate GNSS child class based on processing source.
        """
        if source == 'UNR':
            return UNR_GNSS
        elif source == 'ESESES':
            return ESESES_GNSS
        elif source == 'JPL-SIDESHOW':
            return JPL_SIDESHOW_GNSS
        elif source == 'Generic':
            return Generic_GNSS
        else:
            raise ValueError(f'{source:s} source not supported.')

    def dload_site(self, print_msg=True):
        raise NotImplementedError('Func. dload_site not implemented. Override with child class.')

    def get_stat_lat_lon(self, print_msg=True):
        raise NotImplementedError('Func. get_stat_lat_lon not implemented. Override with child class.')

    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        raise NotImplementedError('Func. read_displacement not implemented. Override with child class.')

    @staticmethod
    def lon_360to180(lon: float) -> float:
        """Convert longitude in the range [0, 360) to
        range (-180, 180].
        """
        if lon > 180:
            lon -= 360
        return lon

    def __crop_to_date_range__(self, start_date: str, end_date: str):
        """Cut out the specified time range.
        start/end_date in format YYYYMMDD
        """
        t_flag = np.ones(len(self.dates), bool)
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

        return None


    def display_data(self, marker_size=2, marker_color='k', plot_errors=True):
        """Display displacement data.
        """
        # instantiate figure and axes
        fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)

        # plot displacement data
        ax[0].scatter(self.dates, self.dis_e, s=marker_size**2, c=marker_color)
        ax[1].scatter(self.dates, self.dis_n, s=marker_size**2, c=marker_color)
        ax[2].scatter(self.dates, self.dis_u, s=marker_size**2, c=marker_color)

        # plot displacement errors
        if plot_errors == True:
            ax[0].errorbar(self.dates, self.dis_e, yerr=self.std_e,
                           linestyle='none', color=marker_color)
            ax[1].errorbar(self.dates, self.dis_n, yerr=self.std_n,
                           linestyle='none', color=marker_color)
            ax[2].errorbar(self.dates, self.dis_u, yerr=self.std_u,
                           linestyle='none', color=marker_color)

        # format plot
        ax[0].set_ylabel('East (m)')
        ax[1].set_ylabel('North (m)')
        ax[2].set_ylabel('Up (m)')
        fig.suptitle(f'{self.site:s} ({self.source:s})')

        plt.show()

        return fig, ax


    #####################################  Utility Functions ###################################
    def displacement_enu2los(self, inc_angle:float, az_angle:float, gnss_comp='enu2los',
                             horz_az_angle=-90., display=False, model=None):
        """Convert displacement in ENU to LOS direction.

        Parameters: inc_angle     - float, LOS incidence angle in degree
                    az_angle      - float, LOS aziuth angle in degree
                                    from the north, defined as positive in clock-wise direction
                    gnss_comp     - str, GNSS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal to fault-parallel
                                    measured from the north with anti-clockwise as positive
        Returns:    dis_los       - 1D np.array for displacement in LOS direction
                    std_los       - 1D np.array for displacement standard deviation in LOS direction
        """
        # get unit vector for the component of interest
        unit_vec = ut.get_unit_vector4component_of_interest(
            los_inc_angle=inc_angle,
            los_az_angle=az_angle,
            comp=gnss_comp.lower(),
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

        # display if requested
        if display == True:
            # instantiate figure and axes
            _, ax = plt.subplots(sharex=True)

            # plot LOS displacement
            ax.scatter(self.dates, self.dis_los, s=2**2,
                       c='k', label='LOS')

            # plot fit if model specified
            if model is not None:
                # specific time_func model
                A = time_func.get_design_matrix4time_func(self.date_list, model=model)
                estm_dis = np.dot(np.linalg.pinv(A), self.dis_los)
                ax.plot(self.dates, estm_dis, 'b', label='model')
            ax.legend()

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
            inc_angle = readfile.read(geom_obj, datasetName='incidenceAngle', box=box,
                                      print_msg=print_msg)[0][0,0]
            az_angle  = readfile.read(geom_obj, datasetName='azimuthAngle', box=box,
                                      print_msg=print_msg)[0][0,0]

        elif isinstance(geom_obj, dict):
            # use mean inc/az_angle from metadata
            inc_angle = ut.incidence_angle(geom_obj, dimension=0, print_msg=print_msg)
            az_angle  = ut.heading2azimuth_angle(float(geom_obj['HEADING']))

        else:
            raise ValueError(f'input geom_obj is neither str nor dict: {geom_obj}')

        return inc_angle, az_angle


    def read_gnss_los_displacement(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                                  gnss_comp='enu2los', horz_az_angle=-90., print_msg=False):
        """Read GNSS displacement in LOS direction.

        Parameters: geom_obj      - dict / str, metadata of InSAR file, or geometry file path
                    start_date    - str, dates in YYYYMMDD format
                    end_date      - str, dates in YYYYMMDD format
                    ref_site      - str, reference GNSS site
                    gnss_comp     - str, GNSS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal
                                    to fault-parallel
        Returns:    dates         - 1D np.array of datetime.datetime object
                    dis/std       - 1D np.array of displacement / uncertainty in meters
                    site_lalo     - tuple of 2 float, lat/lon of GNSS site
                    ref_site_lalo - tuple of 2 float, lat/lon of reference GNSS site
        """
        # read GNSS object
        inc_angle, az_angle = self.get_los_geometry(geom_obj)
        dates = self.read_displacement(start_date, end_date, print_msg=print_msg)[0]
        dis, std = self.displacement_enu2los(inc_angle, az_angle, gnss_comp=gnss_comp,
                                             horz_az_angle=horz_az_angle)
        site_lalo = self.get_stat_lat_lon(print_msg=print_msg)

        # define GNSS station object based on processing source
        GNSS = self.get_gnss_obj_by_source(self.source)

        # get LOS displacement relative to another GNSS site
        if ref_site:
            ref_obj = GNSS(site=ref_site, data_dir=self.data_dir)
            ref_obj.open()
            ref_obj.read_displacement(start_date, end_date, print_msg=print_msg)
            inc_angle, az_angle = ref_obj.get_los_geometry(geom_obj)
            ref_obj.displacement_enu2los(inc_angle, az_angle, gnss_comp=gnss_comp,
                                         horz_az_angle=horz_az_angle)
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


    def get_gnss_los_velocity(self, geom_obj, start_date=None, end_date=None,
                             ref_site=None, gnss_comp='enu2los',
                             horz_az_angle=-90., model=None,
                             print_msg=True):
        """Convert the three-component displacement data into LOS
        velocity.

        Parameters: geom_obj        - dict / str, metadata of InSAR file, or
                                      geometry file path
                    start_date      - str, YYYYMMDD format
                    end_date        - str, YYYYMMDD format
                    ref_site        - str, reference GNSS site
                    gnss_comp       - str, GNSS components used to convert to
                                      LOS direction
                    horz_az_angle   - float, fault azimuth angle used to convert
                                 horizontal to fault-parallel
                    model           - dict, time function model, e.g.
                                      {'polynomial': 1, 'periodic': [1.0, 0.5]}
        Returns:    dates         - 1D np.array, datetime.datetime object
                    dis           - 1D np.array, displacement in meters
                    std           - 1D np.array, displacement uncertainty in meters
                    site_lalo     - tuple of 2 float, lat/lon of GNSS site
                    ref_site_lalo - tuple of 2 float, lat/lon of reference GNSS site
        """
        # retrieve displacement data
        dates, dis = self.read_gnss_los_displacement(geom_obj,
                                                    start_date=start_date,
                                                    end_date=end_date,
                                                    ref_site=ref_site,
                                                    gnss_comp=gnss_comp,
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
            if print_msg == True:
                print(f'Velocity calculation failed for site {self.site:s}')

        return self.velocity, dis



class UNR_GNSS(GNSS):
    """GNSS class for daily solutions processed by UNR NGL.

    This object will assign the attributes:
            site          - str, four-digit site code
            site_lat/lon  - float
            dates         - 1D np.ndarray
            date_list     - list
            dis_e/n/u     - 1D np.ndarray
            std_e,n,u     - 1D np.ndarray

    based on the specific formats of the data source, using the functions:
            dload_site
            get_stat_lat_lon
            read_displacement
    """
    source = 'UNR'

    def dload_site(self, print_msg=True) -> str:
        """Download the station displacement data from the specified source.

        Modifies:   self.file     - str, local file path/name
                    self.file_url - str, file URL
        Returns:    self.file     - str, local file path/name
        """
        if print_msg == True:
            print(f"Downloading data for site {self.site:s} from UNR NGL source")

        # URL and file name specs
        url_prefix = 'http://geodesy.unr.edu/gps_timeseries/tenv3'
        if self.version == 'IGS08':
            self.file = os.path.join(self.data_dir,
                                     '{site:s}.{version:s}.tenv3'.\
                                     format(site=self.site, version=self.version))
        elif self.version == 'IGS14':
            self.file = os.path.join(self.data_dir,
                                    '{site:s}.tenv3'.\
                                    format(site=self.site))
        self.file_url = os.path.join(url_prefix, self.version,
                                     os.path.basename(self.file))

        # download file if not present
        if os.path.exists(self.file):
            if print_msg == True:
                print(f'file {self.file:s} exists--reading')
        else:
            if print_msg == True:
                print(f'... downloading {self.file_url:s} to {self.file:s}')
            urlretrieve(self.file_url, self.file)  #nosec

        return self.file

    def get_stat_lat_lon(self, print_msg=True) -> (float, float):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        if print_msg == True:
            print('calculating station lat/lon')

        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)
        ref_lon, ref_lat = float(data[0, 6]), 0.
        e0, e_off, n0, n_off = data[0, 7:11].astype(float)
        e0 += e_off
        n0 += n_off

        az = np.arctan2(e0, n0) / np.pi * 180.
        dist = np.sqrt(e0**2 + n0**2)
        g = Geod(ellps='WGS84')
        self.site_lon, self.site_lat = g.fwd(ref_lon, ref_lat, az, dist)[0:2]

        if print_msg == True:
            print(f'\t{self.site_lat:f}, {self.site_lon:f}')

        return self.site_lat, self.site_lon

    def read_displacement(self, start_date=None, end_date=None, print_msg=True,
                          display=False):
        """Read GNSS displacement time-series (defined by start/end_date)
        Parameters: start/end_date - str, date in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in float32
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg == True:
            print('reading time and displacement in east/north/vertical direction')

        # read data from file
        data = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)

        # Parse dates
        self.dates = np.array([dt.datetime.strptime(i, "%y%b%d") \
                               for i in data[:,1]])

        # parse displacement data
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = data[:, (8,10,12,14,15,16)].astype(np.float32).T

        # cut out the specified time range
        self.__crop_to_date_range__(start_date, end_date)

        # formulate date list
        self.date_list = [date.strftime('%Y%m%d') for date in self.dates]

        # display if requested
        if display == True:
            self.display_data()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class ESESES_GNSS(GNSS):
    """GNSS class for daily solutions processed by ESESES.

    This object will assign the attributes:
            site          - str, four-digit site code
            site_lat/lon  - float
            dates         - 1D np.ndarray
            date_list     - list
            dis_e/n/u     - 1D np.ndarray
            std_e,n,u     - 1D np.ndarray

    based on the specific formats of the data source, using the functions:
            dload_site
            get_stat_lat_lon
            read_displacement
    """
    source = 'ESESES'

    def dload_site(self, print_msg=True) -> str:
        """Download the station displacement data from the
        specified source.

        Modifies:   self.file     - str, local file path/name
                    self.file_url - str, file URL
        Returns:    self.file     - str, local file path/name
        """
        if print_msg == True:
            print(f'downloading data for site {self.site:s} from the ESESES source')

        # determine proper URL
        url_fmt = 'http://garner.ucsd.edu/pub/measuresESESES_products/Timeseries/CurrentUntarred/Clean_TrendNeuTimeSeries_comb_{:s}'

        # start with today and check back in time
        today = dt.date.today()
        day_lim = 21
        for days in range(day_lim):
            # formulate "days ago"
            days_ago = dt.timedelta(days=days)

            # formulate URL based on date
            url_prefix = url_fmt.format((today - days_ago).strftime('%Y%m%d'))

            # check if page exists
            try:
                urlopen(url_prefix)  #nosec
                break
            except Exception:
                if days_ago.days == (day_lim - 1):
                    raise FileNotFoundError('The ESESES source repository cannot be found.')
                else:
                    pass

        # file name and full url
        self.file = os.path.join(self.data_dir,
                                 '{site:s}CleanTrend.neu.Z'.\
                                 format(site=self.site.lower()))
        self.file_url = os.path.join(url_prefix, os.path.basename(self.file))

        # download file if not present
        if os.path.exists(self.file):
            if print_msg == True:
                print(f'file {self.file:s} exists--reading')
        else:
            if print_msg == True:
                print(f'... downloading {self.file_url:s} to {self.file:s}')
            urlretrieve(self.file_url, self.file)  #nosec

        # unzip file
        with zipfile.ZipFile(self.file, 'r') as Zfile:
            Zfile.extractall(self.data_dir)

        # update file name
        self.file = self.file.strip('.Z')
        if print_msg == True:
            print(f'... extracted to {self.file:s}')

        return self.file

    def get_stat_lat_lon(self, print_msg=True) -> (float, float):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        if print_msg == True:
            print('calculating station lat/lon')

        with open(self.file) as data_file:
            # Read raw file contents
            lines = data_file.readlines()

            # Determine reference latitude
            lat_line = [line for line in lines \
                        if line.find('# Latitude') != -1]
            lat_line = lat_line[0].strip('\n')
            self.site_lat = float(lat_line.split()[-1])

            # Determine reference longitude
            lon_line = [line for line in lines \
                        if line.find('# East Longitude') != -1]
            lon_line = lon_line[0].strip('\n')
            site_lon = float(lon_line.split()[-1])
            self.site_lon = self.lon_360to180(site_lon)

        if print_msg == True:
            print(f'\t{self.site_lat:f}, {self.site_lon:f}')

        return self.site_lat, self.site_lon

    def read_displacement(self, start_date=None, end_date=None, print_msg=True,
                          display=False):
        """Read GNSS displacement time-series (defined by start/end_date)
        Parameters: start/end_date - str, date in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in float32
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg == True:
            print('reading time and displacement in east/north/vertical direction')

        # read data from file
        data = np.loadtxt(self.file, usecols=tuple(range(0,12)))
        n_data = data.shape[0]

        # parse dates
        dates = [dt.datetime(int(data[i,1]), 1, 1) \
                 + dt.timedelta(days=int(data[i,2])) \
                 for i in range(n_data)]
        self.dates = np.array(dates)

        # parse displacement data
        (self.dis_n,
         self.dis_e,
         self.dis_u,
         self.std_n,
         self.std_e,
         self.std_u) = data[:, 3:9].astype(np.float32).T / 1000

        # cut out the specified time range
        self.__crop_to_date_range__(start_date, end_date)

        # formulate date list
        self.date_list = [date.strftime('%Y%m%d') for date in self.dates]

        # display if requested
        if display == True:
            self.display_data()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class JPL_SIDESHOW_GNSS(GNSS):
    """GNSS class for daily solutions processed by JPL-SIDESHOW.

    This object will assign the attributes:
            site          - str, four-digit site code
            site_lat/lon  - float
            dates         - 1D np.ndarray
            date_list     - list
            dis_e/n/u     - 1D np.ndarray
            std_e,n,u     - 1D np.ndarray

    based on the specific formats of the data source, using the functions:
            dload_site
            get_stat_lat_lon
            read_displacement
    """
    source = 'JPL-SIDESHOW'

    def dload_site(self, print_msg=True) -> str:
        """Download the station displacement data from the
        specified source.

        Modifies:   self.file     - str, local file path/name
                    self.file_url - str, file URL
        Returns:    self.file     - str, local file path/name
        """
        if print_msg == True:
            print(f'downloading data for site {self.site:s} from the JPL-SIDESHOW source')

        # URL and file name specs
        url_prefix = 'https://sideshow.jpl.nasa.gov/pub/JPL_GPS_Timeseries/repro2018a/post/point/'
        self.file = os.path.join(self.data_dir, f'{self.site:s}.series')
        self.file_url = os.path.join(url_prefix, os.path.basename(self.file))

        # download file if not present
        if os.path.exists(self.file):
            if print_msg == True:
                print(f'file {self.file:s} exists--reading')
        else:
            if print_msg == True:
                print(f'... downloading {self.file_url:s} to {self.file:s}')
            urlretrieve(self.file_url, self.file)  #nosec

        return self.file

    def get_stat_lat_lon(self, print_msg=True) -> (float, float):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        if print_msg == True:
            print('calculating station lat/lon')

        # need to refer to the site list
        site_list_file = dload_site_list(source='JPL-SIDESHOW')

        # find site in site list file
        with open(site_list_file, 'r') as site_list:
            for line in site_list:
                if (line[:4] == self.site) and (line[5:8] == 'POS'):
                    site_lat, site_lon = line.split()[2:4]

        # format
        self.site_lat = float(site_lat)
        self.site_lon = float(site_lon)

        if print_msg == True:
            print(f'\t{self.site_lat:f}, {self.site_lon:f}')

        return self.site_lat, self.site_lon

    def read_displacement(self, start_date=None, end_date=None, print_msg=True,
                          display=False):
        """Read GNSS displacement time-series (defined by start/end_date)
        Parameters: start/end_date - str, date in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in float32
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg == True:
            print('reading time and displacement in east/north/vertical direction')

        # read data from file
        data = np.loadtxt(self.file)
        n_data = data.shape[0]

        # parse dates
        self.dates = np.array([dt.datetime(*data[i,-6:].astype(int)) for i in range(n_data)])

        # parse displacement data
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = data[:, 1:7].astype(float).T

        # cut out the specified time range
        self.__crop_to_date_range__(start_date, end_date)

        # formulate date list
        self.date_list = [date.strftime('%Y%m%d') for date in self.dates]

        # display if requested
        if display == True:
            self.display_data()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class Generic_GNSS(GNSS):
    """GNSS class for daily solutions of an otherwise-unsupported source.
    The user should format the station position data in a file called
    <sitename>.dat The file should have seven space-separated columns:

        date dis_e dis_n dis_u std_e std_n std_u

    where date is in the format <YYYYMMDD> or <YYYYMMDD>T<HHMMSS>; and
    displacement values are in meters.

    For the generic type, it is necessary to have an accompanying file with
    the site reference coordinates in the current folder. The specifications
    for the GenericList.txt file are given above.

    This object will assign the attributes:
            site          - str, four-digit site code
            site_lat/lon  - float
            dates         - 1D np.ndarray
            date_list     - list
            dis_e/n/u     - 1D np.ndarray
            std_e,n,u     - 1D np.ndarray

    based on the specific formats of the data source, using the functions:
            dload_site
            get_stat_lat_lon
            read_displacement
    """
    source = 'Generic'

    def dload_site(self, print_msg=True) -> str:
        """Download the station displacement data from the
        specified source.

        Modifies:   self.file     - str, local file path/name
                    self.file_url - str, file URL
        Returns:    self.file     - str, local file path/name
        """
        if print_msg == True:
            print(f'reading data for site {self.site:s}')

        # URL and file name specs
        self.file = os.path.join(self.data_dir, f'{self.site:s}.txt')
        self.file_url = ''

        # download file if not present
        if print_msg == True:
            print(f'reading file {self.file:s}')

        return self.file

    def get_stat_lat_lon(self, print_msg=True) -> (str, str):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - str
        Returns:    self.lat/lon - str
        """
        if print_msg == True:
            print('calculating station lat/lon')

        # need to refer to the site list
        site_list_file = 'GenericList.txt'

        # find site in site list file
        with open(site_list_file, 'r') as site_list:
            for line in site_list:
                if line[:4] == self.site:
                    site_lat, site_lon = line.split()[1:3]

        # format
        self.site_lat = float(site_lat)
        self.site_lon = float(site_lon)

        if print_msg == True:
            print(f'\t{self.site_lat:f}, {self.site_lon:f}')

        return self.site_lat, self.site_lon

    def read_displacement(self, start_date=None, end_date=None, print_msg=True,
                          display=False):
        """Read GNSS displacement time-series (defined by start/end_date)
        Parameters: start/end_date - str, date in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in float32
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        if print_msg == True:
            print('reading time and displacement in east/north/vertical direction')

        # parse dates
        with open(self.file, 'r') as data_file:
            lines = data_file.readlines()
        self.dates = []
        for line in lines:
            date = line.split()[0]
            date_len = len(date)

            # format
            if date_len == 8:
                datetime = dt.datetime.strptime(date, '%Y%m%d')
            elif date_len == 15:
                datetime = dt.datetime.strptime(date, '%Y%m%dT%H%M%S')
            else:
                raise ValueError('Date/time format not recognized')

            self.dates.append(datetime)
        self.dates = np.array(self.dates)

        # parse displacement data
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = np.loadtxt(self.file, usecols=tuple(range(1,7))).T

        # cut out the specified time range
        self.__crop_to_date_range__(start_date, end_date)

        # formulate date list
        self.date_list = [date.strftime('%Y%m%d') for date in self.dates]

        # display if requested
        if display == True:
            self.display_data()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


#################################### End of GNSS class ####################################
