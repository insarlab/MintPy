"""Class / utilities for GNSS download / operations."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Robert Zinke, Jul 2018             #
############################################################
# Utility scripts for GNSS handling
# Recommend import:
#     from mintpy.objects import gnss


import csv
import datetime as dt
import glob
import os
from urllib.request import urlopen, urlretrieve

import numpy as np

from mintpy.objects.coord import coordinate
from mintpy.utils import ptime, readfile, time_func, utils1 as ut

GNSS_SITE_LIST_URLS = {
    'UNR'      : 'https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt',
    'ESESES'   : 'http://garner.ucsd.edu/pub/measuresESESES_products/Velocities/ESESES_Velocities.txt',
    'SIDESHOW' : 'https://sideshow.jpl.nasa.gov/post/tables/table2.html',
    'GENERIC'  : None,
}
GNSS_SOURCES = list(GNSS_SITE_LIST_URLS.keys())



######################################### Search GNSS ###############################################

def search_gnss(SNWE, start_date=None, end_date=None, source='UNR', site_list_file=None,
                min_num_solution=50, print_msg=True):
    """Search available GNSS sites within the geo bounding box for a given GNSS source.

    Parameters: SNWE             - tuple of 4 float, indicating (South, North, West, East) in degrees
                source           - str, program or institution that processed the GNSS data
                start_date       - str, date in YYYYMMDD format
                end_date         - str, date in YYYYMMDD format
                site_list_file   - str, site list file name
                min_num_solution - int, minimum number of solutions available
    Returns:    site_names       - 1D np.ndarray in string, GNSS station names
                site_lats        - 1D np.ndarray in float32, latitude
                site_lons        - 1D np.ndarray in float32, longitude
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

    # check: site_list_file name
    if site_list_file is None:
        if source == 'GENERIC':
            raise ValueError('Site list file must be specified for GENERIC GNSS source!')
        else:
            site_list_file = os.path.basename(GNSS_SITE_LIST_URLS[source])

    # download site_list_file (if it does not exist in current directory)
    if not os.path.isfile(site_list_file):
        dload_site_list(site_list_file, source=source, print_msg=print_msg)

    # read site_list_file
    if source == 'UNR':
        sites = read_UNR_site_list(site_list_file)
    elif source == 'ESESES':
        sites = read_ESESES_site_list(site_list_file)
    elif source == 'SIDESHOW':
        sites = read_SIDESHOW_site_list(site_list_file)
    elif source == 'GENERIC':
        sites = read_GENERIC_site_list(site_list_file)

    # ensure that site data formatting is consistent
    sites['site'] = np.array([site.upper() for site in sites['site']])
    # ensure longitude values in (-180, 180]
    sites['lon'] = ut.standardize_longitude(sites['lon'], limit='-180to180')
    vprint(f'load {len(sites["site"]):d} GNSS sites with fields: {" ".join(sites.keys())}')

    # limit in space
    idx = ((sites['lat'] >= SNWE[0]) * (sites['lat'] <= SNWE[1]) *
           (sites['lon'] >= SNWE[2]) * (sites['lon'] <= SNWE[3]))
    vprint(f'keep sites within SNWE of {SNWE}: [{np.sum(idx)}]')

    # limit in time
    if start_date and 'end_date' in sites.keys():
        start_dt = ptime.date_list2vector([start_date])[0][0]
        idx *= sites['end_date'] >= start_dt
        vprint(f'keep sites with end_date >= {start_date}: [{np.sum(idx)}]')
    if end_date and 'start_date' in sites.keys():
        end_dt = ptime.date_list2vector([end_date])[0][0]
        idx *= sites['start_date'] <= end_dt
        vprint(f'keep sites with start_date <= {end_date}: [{np.sum(idx)}]')

    # limit based on number of solutions
    if min_num_solution is not None and 'num_solution' in sites.keys():
        idx *= sites['num_solution'] >= min_num_solution
        vprint(f'keep sites with # of solutions >= {min_num_solution}: [{np.sum(idx)}]')

    # print remaining site names
    vprint(sites['site'][idx])

    return sites['site'][idx], sites['lat'][idx], sites['lon'][idx]


def dload_site_list(out_file=None, source='UNR', print_msg=True) -> str:
    """Download single file with list of GNSS site locations.
    """
    # check source is supported
    assert source in GNSS_SOURCES, f'{source:s} GNSS is NOT supported! Use one of {GNSS_SOURCES}.'

    # determine URL
    site_list_url = GNSS_SITE_LIST_URLS[source]

    # handle output file
    if out_file is None:
        out_file = os.path.basename(site_list_url)

    # download file
    if not os.path.isfile(out_file):
        if print_msg:
            print(f'downloading site list from {source:s}: {site_list_url:s} to {out_file:s}')
        urlretrieve(site_list_url, out_file)  #nosec

    return out_file


def read_UNR_site_list(site_list_file:str):
    """Return names and lon/lat values for UNR GNSS stations.
    """
    fc = np.loadtxt(site_list_file, dtype=str, skiprows=1, usecols=(0,1,2,3,4,5,6,7,8,9,10))
    sites = {
        'site'         : fc[:,0],
        'lat'          : fc[:,1].astype(np.float32),
        'lon'          : fc[:,2].astype(np.float32),
        'start_date'   : fc[:,7],
        'end_date'     : fc[:,8],
        'num_solution' : fc[:,10].astype(np.int16),
    }

    # re-format dates
    sites['start_date'] = np.array([dt.datetime.strptime(x, '%Y-%m-%d') for x in sites['start_date']])
    sites['end_date']   = np.array([dt.datetime.strptime(x, '%Y-%m-%d') for x in sites['end_date']])

    return sites


def read_ESESES_site_list(site_list_file:str):
    """Return names and lon/lat values for JPL/SOPAC ESESES GNSS stations.
    """
    fc = np.loadtxt(site_list_file, skiprows=17, dtype=str)
    sites = {
        'site' : fc[:,0],
        'lon'  : fc[:,1].astype(np.float32),
        'lat'  : fc[:,2].astype(np.float32),
    }
    return sites


def read_SIDESHOW_site_list(site_list_file:str):
    """Return names and lon/lat values for JPL SIDESHOW GNSS stations.
    """
    fc = np.loadtxt(site_list_file, comments='<', skiprows=9, dtype=str)
    sites = {
        'site' : fc[::2, 0],
        'lat'  : fc[::2, 2].astype(np.float32),
        'lon'  : fc[::2, 3].astype(np.float32),
    }
    return sites


def read_GENERIC_site_list(site_list_file:str):
    """Return names and lon/lat values for GNSS stations processed by an
    otherwise-unsupported source.

    The user must format the station position data in a file named
    GenericList.txt The file should have three, nine, or eleven space-
    separated columns:

    SITE lat lon [vel_e vel_n vel_u err_e err_n err_u] [start_date end_date]

    where site is the four-digit, alphanumeric (uppercase) site code; and
    lat/lon are in decimal degrees. If included, vel should be in units of
    m/yr; and dates should be in format YYYYMMDD.
    """
    fc = np.loadtxt(site_list_file, dtype=str)
    sites = {
        'site' : fc[:,0],
        'lon'  : fc[:,1].astype(np.float32),
        'lat'  : fc[:,2].astype(np.float32),
    }

    return sites



######################################### Utils Functions ###########################################

def get_los_obs(meta, obs_type, site_names, start_date, end_date, source='UNR', gnss_comp='enu2los',
                horz_az_angle=-90., model=None, print_msg=True, redo=False):
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
                model      - dict, time function model, e.g. {'polynomial': 1, 'periodic': [1.0, 0.5]}
                print_msg  - bool, print verbose info
                redo       - bool, ignore existing CSV file and re-calculate
    Returns:    site_obs   - 1D np.ndarray(), GNSS LOS velocity or displacement in m or m/yr
    Examples:   from mintpy.objects import gnss
                from mintpy.utils import readfile, utils as ut
                meta = readfile.read_attribute('geo/geo_velocity.h5')
                SNWE = ut.four_corners(meta)
                site_names = gnss.search_gnss(SNWE, start_date='20150101', end_date='20190619')[0]
                vel = gnss.get_los_obs(meta, 'velocity',     site_names, start_date='20150101', end_date='20190619')
                dis = gnss.get_los_obs(meta, 'displacement', site_names, start_date='20150101', end_date='20190619')
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    num_site = len(site_names)

    # obs_type --> obs_ind
    assert obs_type in ['displacement', 'velocity'], f'un-supported obs_type: {obs_type}'
    obs_ind = 3 if obs_type.lower() == 'displacement' else 4

    # GNSS CSV file info
    file_dir = os.path.dirname(meta['FILE_PATH'])
    csv_file = os.path.join(file_dir, f'gnss_{gnss_comp:s}')
    csv_file += f'{horz_az_angle:.0f}' if gnss_comp == 'horz' else ''
    csv_file += f'_{source.upper()}.csv'
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
            vprint(f'use incidence / azimuth angle from file: {os.path.basename(geom_file)}')
        else:
            geom_obj = meta
            vprint('use incidence / azimuth angle from metadata')

        # get url_prefix [to speed up downloading for ESESES]
        url_prefix = get_ESESES_url_prefix() if source == 'ESESES' else None

        # loop for calculation
        prog_bar = ptime.progressBar(maxValue=num_site, print_msg=print_msg)
        for i, site_name in enumerate(site_names):
            prog_bar.update(i+1, suffix=f'{i+1}/{num_site} {site_name:s}')

            # calculate GNSS data value
            gnss_obj = get_gnss_class(source)(site_name, url_prefix=url_prefix)
            vel, dis_ts = gnss_obj.get_los_velocity(
                geom_obj,
                start_date=start_date,
                end_date=end_date,
                gnss_comp=gnss_comp,
                horz_az_angle=horz_az_angle,
                model=model,
            )

            # ignore time-series if the estimated velocity is nan
            dis = np.nan if np.isnan(vel) else dis_ts[-1] - dis_ts[0]

            # save data to list
            data_list.append([gnss_obj.site, gnss_obj.site_lon, gnss_obj.site_lat, dis, vel])
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


def get_baseline_change(dates1, pos_x1, pos_y1, pos_z1,
                        dates2, pos_x2, pos_y2, pos_z2):
    """Calculate the baseline change between two GNSS displacement time-series.

    Parameters: dates1/2   - 1D np.ndarray, datetime.datetime object
                pos_x/y/z1 - 1D np.ndarray in float32, displacement in meters of the 1st site
                pos_x/y/z2 - 1D np.ndarray in float32, displacement in meters of the 2nd site
    Returns:    dates      - 1D np.ndarray in dt.datetime object for the common dates
                bases      - 1D np.ndarray in float32, baseline in meters for the common dates
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


def get_gnss_class(source:str):
    """Return the appropriate GNSS child class based on processing source.
    """
    if source == 'UNR':
        return GNSS_UNR
    elif source == 'ESESES':
        return GNSS_ESESES
    elif source == 'SIDESHOW':
        return GNSS_SIDESHOW
    elif source == 'GENERIC':
        return GNSS_GENERIC
    else:
        raise ValueError(f'GNSS source {source:s} is NOT supported!')


def get_ESESES_url_prefix():
    """Get the url prefix for the ESESES source, which updates regularly.
    [Poor design of ESESES website].
    """
    print('searching for ESESES url_prefix ...')
    # url prefix format
    url_fmt = 'http://garner.ucsd.edu/pub/measuresESESES_products/Timeseries'
    url_fmt += '/CurrentUntarred/Clean_TrendNeuTimeSeries_comb_{:s}'

    # start with today and check back in time
    today = dt.date.today()
    max_day = 21
    num_day = 0
    while num_day < max_day:
        # formulate URL based on date
        day_str = (today - dt.timedelta(days=num_day)).strftime('%Y%m%d')
        url_prefix = url_fmt.format(day_str)

        # check if page exists
        try:
            urlopen(url_prefix)
            print(f'{url_prefix} [YES!]')
        except:
            if num_day == max_day - 1:
                raise FileNotFoundError(f'The ESESES repository {url_fmt} CANNOT be found!')
            else:
                num_day += 1
                print(f'{url_prefix} [no]')
                continue
        else:
            break

    return url_prefix



#################################### GNSS-GSI utility functions #####################################

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



#################################### GNSS parent/child classes ######################################

class GNSS:
    """GNSS parent class for time-series of daily solution.

    The GNSS class is solely meant to be a parent class. Child classes, defined
    below, support functions for downloading and parsing GNSS position based on
    the processing source (e.g., UNR, etc.). Use the `get_gnss_class`
    method to determine appropriate child class.

    The parent class/object will assign the following attributes:
        source        - str, GNSS solution source
        version       - str, GNSS solution version
        url_prefix    - str, GNSS data file url prefix

    The chile class/object will assign the following attributes:
        file          - str, path of the local data file
        url           - str, path of the remote data file
        site          - str, four-digit site code
        site_lat/lon  - float, site latitude/longitude in degree
        dates         - 1D np.ndarray, dt.datetime object
        date_list     - list(str), dates in YYYYMMDD format
        dis_e/n/u     - 1D np.ndarray, displacement in meters
        std_e,n,u     - 1D np.ndarray, displacement STD in meters
    """

    def __init__(self, site: str, data_dir=None, version='IGS20', source='UNR', url_prefix=None):
        # site info
        self.site = site
        self.source = source
        self.version = version
        self.url_prefix = url_prefix
        self.url = None

        # local file/dir info
        self.file = None
        self.data_dir = data_dir if data_dir else os.path.abspath(f'GNSS-{source.upper()}')

        # ensure local dir exists
        if not os.path.exists(self.data_dir):
            print('create directory:', self.data_dir)
            os.mkdir(self.data_dir)

        # displacement data
        self.dates = None
        self.date_list = None
        self.dis_e = None
        self.dis_n = None
        self.dis_u = None
        self.std_e = None
        self.std_n = None
        self.std_u = None


    def open(self, file=None, print_msg=True):
        """Read the lat/lon and displacement data of the station.
        Download if necessary.
        """
        # download file if not present
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # retrieve data from file
        self.get_site_lat_lon()
        self.read_displacement(print_msg=print_msg)


    def dload_site(self, overwrite=False, total_tries=5, print_msg=True):
        """Download GNSS site data file.

        Parameters: overwrite   - bool, overwrite existing data file
                    total_tries - int, number of tries to download if failed
                    print_msg   - bool, verbose print out msg
        Returns:    self.file   - str, path to the local data file
        """
        vprint = print if print_msg else lambda *args, **kwargs: None
        if self.url and overwrite or not os.path.isfile(self.file):
            vprint(f"downloading site {self.site:s} from {self.source} to {self.file:s}")
            # retry on download fail
            # https://stackoverflow.com/questions/31529151
            remain_tries = total_tries
            while remain_tries > 0 :
                try:
                    urlretrieve(self.url, self.file)
                    vprint(f'successfully downloaded: {self.url}')
                except:
                    vprint(f'error downloading {self.url} on trial no. {total_tries-remain_tries}')
                    remain_tries -= 1
                    continue
                else:
                    break
        return self.file


    def get_site_lat_lon(self, print_msg=False):
        """Get the GNSS site latitude & longitude into:
        Returns: site_lat/lon - float, site latitude/longitude in degree
        """
        raise NotImplementedError('get_site_lat_lon() is NOT implemented. Override with child class.')

    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """Get the GNSS time/displacement(Std) into:
        Returns: dates      - 1D np.ndarray in datetime.datetime object
                 date_list  - list(str), date in YYYYMMDD format
                 dis_e/n/u  - 1D np.ndarray in float32, displacement in meters
                 std_e/n/u  - 1D np.ndarray in float32, displacement STD in meters
        """
        raise NotImplementedError('read_displacement() is NOT implemented. Override with child class.')


    def _crop_to_date_range(self, start_date: str, end_date: str):
        """Crop the time-series given the start/end_date in format YYYYMMDD,
        and create date_list from dates.
        """
        flag = np.ones(len(self.dates), dtype=bool)
        if start_date:
            t0 = ptime.date_list2vector([start_date])[0][0]
            flag[self.dates < t0] = 0
        if end_date:
            t1 = ptime.date_list2vector([end_date])[0][0]
            flag[self.dates > t1] = 0

        self.dates = self.dates[flag]
        self.dis_e = self.dis_e[flag]
        self.dis_n = self.dis_n[flag]
        self.dis_u = self.dis_u[flag]
        self.std_e = self.std_e[flag]
        self.std_n = self.std_n[flag]
        self.std_u = self.std_u[flag]

        # create member var: date_list
        self.date_list = [x.strftime('%Y%m%d') for x in self.dates]


    #####################################  Utility Functions ###################################
    def plot(self, marker_size=2, marker_color='k', plot_error_bar=True):
        """Plot the displacement time-series.
        """
        import matplotlib.pyplot as plt

        if self.dis_e is None:
            self.open()

        # instantiate figure and axes
        fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)

        # plot displacement data
        kwargs = dict(s=marker_size**2, c=marker_color)
        ax[0].scatter(self.dates, self.dis_e, **kwargs)
        ax[1].scatter(self.dates, self.dis_n, **kwargs)
        ax[2].scatter(self.dates, self.dis_u, **kwargs)

        # plot displacement errors
        if plot_error_bar:
            kwargs = dict(linestyle='none', color=marker_color)
            ax[0].errorbar(self.dates, self.dis_e, yerr=self.std_e, **kwargs)
            ax[1].errorbar(self.dates, self.dis_n, yerr=self.std_n, **kwargs)
            ax[2].errorbar(self.dates, self.dis_u, yerr=self.std_u, **kwargs)

        # format plot
        for i, label in enumerate(['East', 'North', 'Up']):
            ax[i].set_ylabel(f'{label} [m]')
        fig.suptitle(f'{self.site:s} ({self.source:s})')
        fig.tight_layout()

        plt.show()

        return fig, ax


    def displacement_enu2los(self, inc_angle:float, az_angle:float, gnss_comp='enu2los',
                             horz_az_angle=-90.):
        """Convert displacement in ENU to LOS direction.

        Parameters: inc_angle     - float, LOS incidence angle in degree
                    az_angle      - float, LOS aziuth    angle in degree from the north,
                                    defined as positive in clock-wise direction
                    gnss_comp     - str, GNSS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal to fault-parallel
                                    measured from the north with anti-clockwise as positive
        Returns:    dis_los       - 1D np.ndarray for displacement in LOS direction
                    std_los       - 1D np.ndarray for displacement standard deviation in LOS direction
        """
        if self.dis_e is None:
            self.open()

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

        return self.dis_los, self.std_los


    def get_los_geometry(self, geom_obj, print_msg=False):
        """Get the Line-of-Sight geometry info in incidence and azimuth angle in degrees."""
        lat, lon = self.get_site_lat_lon()

        # get LOS geometry
        if isinstance(geom_obj, str):
            # geometry file
            atr = readfile.read_attribute(geom_obj)
            coord = coordinate(atr, lookup_file=geom_obj)
            y, x = coord.geo2radar(lat, lon, print_msg=print_msg)[0:2]
            # check against image boundary
            y = max(0, y);  y = min(int(atr['LENGTH'])-1, y)
            x = max(0, x);  x = min(int(atr['WIDTH'])-1, x)
            kwargs = dict(box=(x,y,x+1,y+1), print_msg=print_msg)
            inc_angle = readfile.read(geom_obj, datasetName='incidenceAngle', **kwargs)[0][0,0]
            az_angle  = readfile.read(geom_obj, datasetName='azimuthAngle',   **kwargs)[0][0,0]

        elif isinstance(geom_obj, dict):
            # use mean inc/az_angle from metadata
            inc_angle = ut.incidence_angle(geom_obj, dimension=0, print_msg=print_msg)
            az_angle  = ut.heading2azimuth_angle(float(geom_obj['HEADING']))

        else:
            raise ValueError(f'input geom_obj is neither str nor dict: {geom_obj}')

        return inc_angle, az_angle


    def get_los_displacement(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                             gnss_comp='enu2los', horz_az_angle=-90., print_msg=False):
        """Get GNSS displacement in LOS direction.

        Parameters: geom_obj      - dict / str, metadata of InSAR file, or geometry file path
                    start_date    - str, dates in YYYYMMDD format
                    end_date      - str, dates in YYYYMMDD format
                    ref_site      - str, reference GNSS site
                    gnss_comp     - str, GNSS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal
                                    to fault-parallel
        Returns:    dates         - 1D np.ndarray of datetime.datetime object
                    dis/std       - 1D np.ndarray of displacement / uncertainty in meters
                    site_lalo     - tuple of 2 float, lat/lon of GNSS site
                    ref_site_lalo - tuple of 2 float, lat/lon of reference GNSS site
        """
        # read GNSS object
        site_lalo = self.get_site_lat_lon()
        dates = self.read_displacement(start_date, end_date, print_msg=print_msg)[0]
        inc_angle, az_angle = self.get_los_geometry(geom_obj)
        dis, std = self.displacement_enu2los(
            inc_angle, az_angle,
            gnss_comp=gnss_comp,
            horz_az_angle=horz_az_angle,
        )

        # get LOS displacement relative to another GNSS site
        if ref_site:
            ref_obj = get_gnss_class(self.source)(site=ref_site, data_dir=self.data_dir)
            ref_site_lalo = ref_obj.get_site_lat_lon()
            ref_obj.read_displacement(start_date, end_date, print_msg=print_msg)
            inc_angle, az_angle = ref_obj.get_los_geometry(geom_obj)
            ref_obj.displacement_enu2los(
                inc_angle, az_angle,
                gnss_comp=gnss_comp,
                horz_az_angle=horz_az_angle,
            )

            # get relative LOS displacement on common dates
            dates = np.array(sorted(list(set(self.dates) & set(ref_obj.dates))))
            dis = np.zeros(dates.shape, dtype=np.float32)
            std = np.zeros(dates.shape, dtype=np.float32)
            for i, date_i in enumerate(dates):
                idx1 = np.where(self.dates == date_i)[0][0]
                idx2 = np.where(ref_obj.dates == date_i)[0][0]
                dis[i] = self.dis_los[idx1] - ref_obj.dis_los[idx2]
                std[i] = (self.std_los[idx1]**2 + ref_obj.std_los[idx2]**2)**0.5
        else:
            ref_site_lalo = None

        return dates, dis, std, site_lalo, ref_site_lalo


    def get_los_velocity(self, geom_obj, start_date=None, end_date=None, ref_site=None,
                         gnss_comp='enu2los', horz_az_angle=-90., model=None, print_msg=True):
        """Convert the three-component displacement data into LOS velocity.

        Parameters: geom_obj      - dict / str, metadata of InSAR file, or
                                    geometry file path
                    start_date    - str, YYYYMMDD format
                    end_date      - str, YYYYMMDD format
                    ref_site      - str, reference GNSS site
                    gnss_comp     - str, GNSS components used to convert to LOS direction
                    horz_az_angle - float, fault azimuth angle used to convert horizontal
                                    to fault-parallel
                    model         - dict, time function model, e.g.
                                      {'polynomial': 1, 'periodic': [1.0, 0.5]}
        Returns:    dates         - 1D np.ndarray, datetime.datetime object
                    dis           - 1D np.ndarray, displacement in meters
        """
        # retrieve displacement data
        dates, dis = self.get_los_displacement(
            geom_obj,
            start_date=start_date,
            end_date=end_date,
            ref_site=ref_site,
            gnss_comp=gnss_comp,
            horz_az_angle=horz_az_angle,
        )[:2]

        # displacement -> velocity
        # if 1. num of observations > 2 AND
        #    2. time overlap > 1/4
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
            if print_msg:
                print(f'\nVelocity calculation failed for site {self.site:s}')

        return self.velocity, dis


class GNSS_UNR(GNSS):
    """GNSS child class for daily solutions processed by Nevada Geodetic Lab
    at University of Nevada, Reno (UNR).

    Website: http://geodesy.unr.edu/NGLStationPages/GlobalStationList

    Reference:
      Blewitt, G., Hammond, W., & Kreemer, C. (2018). Harnessing the GPS data
        explosion for interdisciplinary science. Eos, 99. doi:10.1029/2018EO104623

    """
    def __init__(self, site: str, data_dir=None, version='IGS20', url_prefix=None):
        super().__init__(
            site=site,
            data_dir=data_dir,
            version=version,
            source='UNR',
            url_prefix=url_prefix,
        )

        # get file
        if version == 'IGS08':
            self.file = os.path.join(self.data_dir, f'{self.site:s}.{version:s}.tenv3')
        elif version == 'IGS14' or version == 'IGS20':
            self.file = os.path.join(self.data_dir, f'{self.site:s}.tenv3')
        else:
            raise ValueError(f'Un-supported GNSS version: {version}!')

        # get url
        # examples: http://geodesy.unr.edu/gps_timeseries/tenv3/IGS08/1LSU.IGS08.tenv3
        #           http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/CASU.tenv3
        #           https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/CAKG.tenv3
        if not self.url_prefix:
            if version == 'IGS08' or version == 'IGS14':
                self.url_prefix = f'https://geodesy.unr.edu/gps_timeseries/tenv3/{self.version}'
            if version == 'IGS20':
                self.url_prefix = f'https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/{self.version}'

        self.url = os.path.join(self.url_prefix, os.path.basename(self.file))


    def dload_site(self, overwrite=False, total_tries=5, print_msg=True):
        """Download GNSS site data file.

        Parameters: overwrite   - bool, overwrite existing data file
                    total_tries - int, number of tries to download if failed
                    print_msg   - bool, verbose print out msg
        Returns:    self.file   - str, path to the local data file
        """
        # download data file via the parent class member function
        super().dload_site(overwrite=overwrite, print_msg=print_msg)

        # download time-series plot file
        # example link: http://geodesy.unr.edu/tsplots/IGS08/TimeSeries/CAMO.png
        #               http://geodesy.unr.edu/tsplots/IGS14/IGS14/TimeSeries/CASU.png
        #               https://geodesy.unr.edu/gps_timeseries/IGS20/tsplots/IGS20/TimeSeries/HAND.png
        plot_file = os.path.join(self.data_dir, f'pic/{self.site}.png')

        # ensure local plot directory exists
        if not os.path.exists(os.path.dirname(plot_file)):
            os.makedirs(os.path.dirname(plot_file), exist_ok=True)

        # get plot file url
        url_prefix = {
            'IGS08' : 'https://geodesy.unr.edu/tsplots/IGS08/TimeSeries',
            'IGS14' : 'https://geodesy.unr.edu/tsplots/IGS14/IGS14/TimeSeries',
            'IGS20' : 'https://geodesy.unr.edu/gps_timeseries/IGS20/tsplots/IGS20/TimeSeries',
        }[self.version]
        plot_file_url = os.path.join(url_prefix, f'{self.site}.png')

        # download
        urlretrieve(plot_file_url, plot_file)

        return self.file


    def get_site_lat_lon(self, print_msg=False) -> (float, float):
        """Get station lat/lon from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        data = np.loadtxt(self.file, dtype=bytes, skiprows=1, max_rows=10)
        self.site_lat, self.site_lon = data[0, 20:22].astype(float)
        # ensure longitude in the range of (-180, 180]
        self.site_lon = ut.standardize_longitude(self.site_lon, limit='-180to180')

        return self.site_lat, self.site_lon


    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """Read GNSS displacement time-series (defined by start/end_date)
        Parameters: start_date - str, start date in YYYYMMDD format
                    end_date   - str, end_date   in YYYYMMDD format
        Returns:    dates      - 1D np.ndarray in datetime.datetime object
                    dis_e/n/u  - 1D np.ndarray in float32, displacement in meters
                    std_e/n/u  - 1D np.ndarray in float32, displacement STD in meters
        """
        vprint = print if print_msg else lambda *args, **kwargs: None

        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # read dates, dis_e, dis_n, dis_u
        vprint('reading time and displacement in east/north/vertical direction')
        try:
            fc = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)
        except:
            msg = 'Error occurred while reading, probably due to interuptions during previous downloading. '
            msg += 'Remove the file and re-download.'
            print(msg)
            self.dload_site(overwrite=True, print_msg=print_msg)
            fc = np.loadtxt(self.file, dtype=bytes, skiprows=1).astype(str)

        self.dates = np.array([dt.datetime.strptime(x, "%y%b%d") for x in fc[:, 1]])
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = fc[:, (8,10,12,14,15,16)].astype(np.float32).T

        # cut out the specified time range
        self._crop_to_date_range(start_date, end_date)

        # display if requested
        if display:
            self.plot()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class GNSS_ESESES(GNSS):
    """GNSS child class for daily solutions processed for the Enhanced Solid
    Earth Science ESDR System (ESESES) project by JPL and SOPAC.

    Website: https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/ESESES_products.html
             http://garner.ucsd.edu/pub/measuresESESES_products/
    """
    def __init__(self, site: str, data_dir=None, version='IGS14', url_prefix=None):
        super().__init__(
            site=site,
            data_dir=data_dir,
            version=version,
            source='ESESES',
            url_prefix=url_prefix,
        )

        # get file
        self.file = os.path.join(self.data_dir, f'{self.site.lower():s}CleanTrend.neu.Z')

        # get url
        # moved to GNSS_ESESES.dload_site() to avoid searching url_prefix
        # when downloading is not needed.


    def dload_site(self, overwrite=False, total_tries=5, print_msg=True):
        """Download GNSS data file.
        """
        from zipfile import ZipFile

        # get url
        if not self.url_prefix:
            self.url_prefix = get_ESESES_url_prefix()
        self.url = os.path.join(self.url_prefix, os.path.basename(self.file))

        # call parent class to download
        super().dload_site(overwrite=overwrite, print_msg=print_msg)

        # uncompress the downloaded *.z file [for ESESES only]
        with ZipFile(self.file, 'r') as fz:
            fz.extractall(self.data_dir)
        self.file = self.file.strip('.Z')    # update file name
        if print_msg:
            print(f'... extracted to {self.file:s}')

        return self.file


    def get_site_lat_lon(self, print_msg=False) -> (float, float):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # use the uncompressed data file
        if self.file.endswith('.Z'):
            self.file = self.file[:-2]

        with open(self.file) as f:
            lines = f.readlines()

            # latitude
            lat_line = [x for x in lines if x.startswith('# Latitude')][0].strip('\n')
            self.site_lat = float(lat_line.split()[-1])

            # longitude
            lon_line = [x for x in lines if x.startswith('# East Longitude')][0].strip('\n')
            self.site_lon = float(lon_line.split()[-1])

        # ensure longitude in the range of (-180, 180]
        self.site_lon = ut.standardize_longitude(self.site_lon, limit='-180to180')

        return self.site_lat, self.site_lon


    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """Read GNSS displacement time-series (defined by start/end_date).

        Parameters: start/end_date - str, date in YYYYMMDD format
        Returns:    dates          - 1D np.ndarray of datetime.datetime object
                    dis_e/n/u      - 1D np.ndarray of displacement in meters in float32
                    std_e/n/u      - 1D np.ndarray of displacement STD in meters in float32
        """
        vprint = print if print_msg else lambda *args, **kwargs: None

        # download file if it does not exist
        if not os.path.isfile(self.file):
            self.dload_site(print_msg=print_msg)

        # use the uncompressed data file
        if self.file.endswith('.Z'):
            self.file = self.file[:-2]

        # read data file
        # use the first 9 cols only, as some epoches miss 10-13 cols: CorrNE/NU/EU, Chi-Squared
        vprint('reading time and displacement in east/north/vertical direction')
        fc = np.loadtxt(self.file, usecols=tuple(range(0,9)))
        num_solution = fc.shape[0]

        # parse dates
        dates = [dt.datetime(int(fc[i, 1]), 1, 1) + dt.timedelta(days=int(fc[i, 2]))
                 for i in range(num_solution)]
        self.dates = np.array(dates)

        # parse displacement data
        (self.dis_n,
         self.dis_e,
         self.dis_u,
         self.std_n,
         self.std_e,
         self.std_u) = fc[:, 3:9].astype(np.float32).T / 1000

        # cut out the specified time range
        self._crop_to_date_range(start_date, end_date)

        # display if requested
        if display:
            self.plot()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class GNSS_SIDESHOW(GNSS):
    """GNSS class for daily solutions processed by JPL SIDESHOW,
    funded by NASA's Space Geodesy Task.

    Website: https://sideshow.jpl.nasa.gov/pub/
             https://sideshow.jpl.nasa.gov/post/series.html

    Reference:
      Heflin, M., Donnellan, A., Parker, J., Lyzenga, G., Moore, A., Ludwig, L. G., et al.
        (2020). Automated Estimation and Tools to Extract Positions, Velocities, Breaks, and
        Seasonal Terms From Daily GNSS Measurements: Illuminating Nonlinear Salton Trough
        Deformation. Earth and Space Science, 7(7), e2019EA000644, doi:10.1029/2019EA000644
    """
    def __init__(self, site: str, data_dir=None, version='IGS14', url_prefix=None):
        super().__init__(
            site=site,
            data_dir=data_dir,
            version=version,
            source='SIDESHOW',
            url_prefix=url_prefix,
        )

        # get file
        self.file = os.path.join(self.data_dir, f'{self.site:s}.series')

        # get url
        if not self.url_prefix:
            self.url_prefix = 'https://sideshow.jpl.nasa.gov/pub/JPL_GPS_Timeseries/repro2018a/post/point'
        self.url = os.path.join(self.url_prefix, os.path.basename(self.file))


    def get_site_lat_lon(self, print_msg=False) -> (float, float):
        """Get station lat/lon based on processing source.
        Retrieve data from the displacement file.

        Modifies:   self.lat/lon - float
        Returns:    self.lat/lon - float
        """
        # need to refer to the site list
        site_list_file = os.path.basename(GNSS_SITE_LIST_URLS['SIDESHOW'])
        if not os.path.exists(site_list_file):
            dload_site_list(out_file=site_list_file, source=self.source, print_msg=True)

        # find site in site list file
        with open(site_list_file) as site_list:
            for line in site_list:
                if (line[:4] == self.site) and (line[5:8] == 'POS'):
                    site_lat, site_lon = line.split()[2:4]

        # format
        self.site_lat = float(site_lat)
        self.site_lon = float(site_lon)
        # ensure longitude in the range of (-180, 180]
        self.site_lon = ut.standardize_longitude(self.site_lon, limit='-180to180')

        if print_msg == True:
            print(f'\t{self.site_lat:f}, {self.site_lon:f}')

        return self.site_lat, self.site_lon


    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
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
         self.std_u) = data[:, 1:7].astype(np.float32).T

        # cut out the specified time range
        self._crop_to_date_range(start_date, end_date)

        # display if requested
        if display == True:
            self.plot()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)


class GNSS_GENERIC(GNSS):
    """GNSS class for daily solutions of an otherwise-unsupported source.

    Required local files:
    1. GenericList.txt: site list file, with the following 3, 9 or 11
        space-separated columns:

    SITE lat lon [vel_e vel_n vel_u err_e err_n err_u] [start_date end_date]

    where site is the four-digit, alphanumeric (uppercase) site code; and
    lat/lon are in decimal degrees. If included, vel should be in units of
    m/yr; and dates should be in format YYYYMMDD.

    2. <sitename>.txt: data file for each site, with 7 space-separated columns:

    date dis_e dis_n dis_u std_e std_n std_u

    where date is in the format <YYYYMMDD> or <YYYYMMDD>T<HH:MM:SS> or <YYYYMMDD>:<HHMMSS>,
    displacement values are in meters.
    """
    def __init__(self, site: str, data_dir=None, version='IGS14', url_prefix=None):
        super().__init__(
            site=site,
            data_dir=data_dir,
            version=version,
            source='GENERIC',
            url_prefix=url_prefix,
        )

        # get file
        self.file = os.path.join(self.data_dir, f'{self.site:s}.txt')

        # get url
        self.url_prefix = ''
        self.url = ''


    def get_site_lat_lon(self, print_msg=False) -> (str, str):
        """Get station lat/lon based on processing source.
        """
        sites = read_GENERIC_site_list('GenericList.txt')
        ind = sites['site'].tolist().index(self.site)
        site_lat, site_lon = sites['lat'][ind], sites['lon'][ind]
        # ensure longitude in the range of (-180, 180]
        site_lon = ut.standardize_longitude(site_lon, limit='-180to180')

        return site_lat, site_lon


    def read_displacement(self, start_date=None, end_date=None, print_msg=True, display=False):
        """Read GNSS displacement time-series (defined by start/end_date).

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
        fc = np.loadtxt(self.file)

        # parse dates
        self.dates = np.array(ptime.date_list2vector(fc[:, 0]))

        # parse displacement data
        (self.dis_e,
         self.dis_n,
         self.dis_u,
         self.std_e,
         self.std_n,
         self.std_u) = fc[:, tuple(range(1,7))].astype(np.float32).T

        # cut out the specified time range
        self._crop_to_date_range(start_date, end_date)

        # display if requested
        if display:
            self.plot()

        return (self.dates,
                self.dis_e, self.dis_n, self.dis_u,
                self.std_e, self.std_n, self.std_u)
