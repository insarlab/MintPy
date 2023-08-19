"""Class for comparing InSAR with GPS."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#   from mintpy.objects.insar_vs_gps import insar_vs_gps


import datetime as dt
import sys

import numpy as np
from scipy import stats
from scipy.interpolate import griddata

from mintpy.defaults.plot import *
from mintpy.objects import giantTimeseries, timeseries
from mintpy.objects.gps import GPS
from mintpy.utils import readfile, utils as ut


############################## beginning of insar_vs_gps class ##############################
class insar_vs_gps:
    """ Comparing InSAR time-series with GPS time-series in LOS direction
    Parameters: ts_file        : str, time-series HDF5 file
                geom_file      : str, geometry HDF5 file
                temp_coh_file  : str, temporal coherence HDF5 file
                site_names     : list of str, GPS site names
                gps_dir        : str, directory of the local GPS data files
                ref_site       : str, common reference site in space for InSAR and GPS
                start/end_date : str, date in YYYYMMDD format for the start/end date
                min_ref_date   : str, date in YYYYMMDD format for the earliest common
                    reference date between InSAR and GPS
    Returns:    ds : dict, each element has the following components:
                    'GV03': {
                      'name': 'GV03',
                      'lat': -0.7977926892712729,
                      'lon': -91.13294444114553,
                      'gps_datetime': array([datetime.datetime(2014, 11, 1, 0, 0),
                             datetime.datetime(2014, 11, 2, 0, 0),
                             ...,
                             datetime.datetime(2018, 6, 25, 0, 0)], dtype=object),
                      'gps_dis': array([-2.63673663e-02, ..., 6.43612206e-01], dtype=float32),
                      'gps_std': array([0.00496152, ..., 0.00477411], dtype=float32),
                      'reference_site': 'GV01',
                      'insar_datetime': array([datetime.datetime(2014, 12, 13, 0, 0),
                             datetime.datetime(2014, 12, 25, 0, 0),
                             ...,
                             datetime.datetime(2018, 6, 19, 0, 0)], dtype=object),
                      'insar_dis_linear': array([-0.01476493, ...,  0.62273948]),
                      'temp_coh': 0.9961861392598478,
                      'gps_std_mean': 0.004515478,
                      'comm_dis_gps': array([-0.02635017, ..., 0.61315614], dtype=float32),
                      'comm_dis_insar': array([-0.01476493, ..., 0.60640174], dtype=float32),
                      'r_square': 0.9993494518609801,
                      'dis_rmse': 0.008023425326946351
                    }
    """

    def __init__(self, ts_file, geom_file, temp_coh_file,
                 site_names, gps_dir='./GPS', ref_site='GV01',
                 start_date=None, end_date=None, min_ref_date=None):
        self.insar_file = ts_file
        self.geom_file = geom_file
        self.temp_coh_file = temp_coh_file
        self.site_names = site_names
        self.gps_dir = gps_dir
        self.ref_site = ref_site
        self.num_site = len(site_names)
        self.ds = {}
        self.start_date = start_date
        self.end_date = end_date
        self.min_ref_date = min_ref_date

    def open(self):
        atr = readfile.read_attribute(self.insar_file)
        k = atr['FILE_TYPE']
        if k == 'timeseries':
            ts_obj = timeseries(self.insar_file)
        elif k == 'giantTimeseries':
            ts_obj = giantTimeseries(self.insar_file)
        else:
            raise ValueError(f'Un-supported time-series file: {k}')
        ts_obj.open(print_msg=False)
        self.metadata = dict(ts_obj.metadata)
        self.num_date = ts_obj.numDate
        # remove time info from insar_datetime to be consistent with gps_datetime
        self.insar_datetime = np.array([i.replace(hour=0, minute=0, second=0, microsecond=0)
                                        for i in ts_obj.times])

        # default start/end_date & min_ref_date
        dt_buffer = dt.timedelta(days=30)
        self.start_date = self.start_date if self.start_date else (ts_obj.times[0] - dt_buffer).strftime('%Y%m%d')
        self.end_date = self.end_date if self.end_date else (ts_obj.times[-1] + dt_buffer).strftime('%Y%m%d')
        self.min_ref_date = self.min_ref_date if self.min_ref_date else ts_obj.times[5].strftime('%Y%m%d')
        if self.min_ref_date not in ts_obj.dateList:
            msg = f'min_ref_date {self.min_ref_date} does NOT exist in InSAR file: {self.insar_file}'
            raise ValueError(msg)

        self.read_gps()
        self.read_insar()
        self.calculate_rmse()
        return

    def read_gps(self):
        for sname in self.site_names:
            site = {}
            site['name'] = sname
            gps_obj = GPS(sname, data_dir=self.gps_dir)
            gps_obj.open(print_msg=False)
            site['lat'] = gps_obj.site_lat
            site['lon'] = gps_obj.site_lon

            dates, dis, dis_std = gps_obj.read_gps_los_displacement(
                self.geom_file,
                start_date=self.start_date,
                end_date=self.end_date,
                ref_site=self.ref_site,
                gps_comp='enu2los',
            )[0:3]
            site['gps_datetime'] = dates
            site['gps_dis'] = dis
            site['gps_std'] = dis_std

            site['reference_site'] = self.ref_site
            self.ds[sname] = site

            sys.stdout.write(f'\rreading GPS {sname}')
            sys.stdout.flush()
        print()
        return

    def read_insar(self):
        # 2.1 prepare interpolation
        coord = ut.coordinate(self.metadata, lookup_file=self.geom_file)
        lats = [self.ds[k]['lat'] for k in self.ds.keys()]
        lons = [self.ds[k]['lon'] for k in self.ds.keys()]
        geo_box = (min(lons), max(lats), max(lons), min(lats))     #(W, N, E, S)
        pix_box = coord.bbox_geo2radar(geo_box)     #(400, 1450, 550, 1600)
        src_lat = readfile.read(self.geom_file, datasetName='latitude', box=pix_box)[0].reshape(-1,1)
        src_lon = readfile.read(self.geom_file, datasetName='longitude', box=pix_box)[0].reshape(-1,1)
        src_pts = np.hstack((src_lat, src_lon))

        dest_pts = np.zeros((self.num_site, 2))
        for i in range(self.num_site):
            site = self.ds[self.site_names[i]]
            dest_pts[i,:] = site['lat'], site['lon']

        # 2.2 interpolation - displacement / temporal coherence
        interp_method = 'linear'   #nearest, linear, cubic
        src_value, atr = readfile.read(self.insar_file, box=pix_box)
        src_value = src_value.reshape(self.num_date, -1)
        if atr['FILE_TYPE'] == 'giantTimeseries':
            src_value *= 0.001
        insar_dis = np.zeros((self.num_site, self.num_date))
        for i in range(self.num_date):
            insar_dis[:,i] = griddata(src_pts, src_value[i,:], dest_pts, method=interp_method)
            sys.stdout.write(('\rreading InSAR acquisition {}/{}'
                              ' with {} interpolation').format(i+1, self.num_date, interp_method))
            sys.stdout.flush()
        print()

        print('reading temporal coherence')
        src_value = readfile.read(self.temp_coh_file, box=pix_box)[0].flatten()
        temp_coh = griddata(src_pts, src_value, dest_pts, method=interp_method)

        # 2.3 write interpolation result
        self.insar_dis_name = f'insar_dis_{interp_method}'
        insar_dis_ref = insar_dis[self.site_names.index(self.ref_site),:]
        for i in range(self.num_site):
            site = self.ds[self.site_names[i]]
            site['insar_datetime'] = self.insar_datetime
            # reference insar to the precise location in space
            site[self.insar_dis_name] = insar_dis[i,:] - insar_dis_ref
            site['temp_coh'] = temp_coh[i]

        # 2.4 reference insar and gps to a common date
        print('reference insar and gps to a common date')
        for i in range(self.num_site):
            site = self.ds[self.site_names[i]]
            gps_date = site['gps_datetime']
            insar_date = site['insar_datetime']

            # find common reference date
            ref_date = dt.datetime.strptime(self.min_ref_date, "%Y%m%d")
            ref_idx = insar_date.tolist().index(ref_date)
            while ref_idx < self.num_date:
                if insar_date[ref_idx] not in gps_date:
                    ref_idx += 1
                else:
                    break
            if ref_idx == self.num_date:
                msg = f"InSAR and GPS do not share ANY date for site: {site['name']}"
                raise RuntimeError(msg)
            comm_date = insar_date[ref_idx]

            # reference insar in time
            site[self.insar_dis_name] -= site[self.insar_dis_name][ref_idx]
            # reference gps dis/std in time
            ref_idx_gps = np.where(gps_date == comm_date)[0][0]
            site['gps_dis'] -= site['gps_dis'][ref_idx_gps]
            site['gps_std'] = np.sqrt(site['gps_std']**2 + site['gps_std'][ref_idx_gps]**2)
            site['gps_std_mean'] = np.mean(site['gps_std'])
        return


    def calculate_rmse(self):
        ## 3. calculate RMSE
        for i in range(self.num_site):
            site = self.ds[self.site_names[i]]
            gps_date = site['gps_datetime']
            insar_date = site['insar_datetime']
            comm_dates = np.array(sorted(list(set(gps_date) & set(insar_date))))
            num_comm_date = len(comm_dates)

            # get displacement at common dates
            comm_dis_insar = np.zeros(num_comm_date, np.float32)
            comm_dis_gps   = np.zeros(num_comm_date, np.float32)
            for j in range(num_comm_date):
                idx1 = np.where(gps_date   == comm_dates[j])[0][0]
                idx2 = np.where(insar_date == comm_dates[j])[0][0]
                comm_dis_gps[j]   = site['gps_dis'][idx1]
                comm_dis_insar[j] = site[self.insar_dis_name][idx2]
            site['comm_dis_gps'] = comm_dis_gps
            site['comm_dis_insar'] = comm_dis_insar
            site['r_square'] = stats.linregress(comm_dis_gps, comm_dis_insar)[2]
            site['dis_rmse'] = np.sqrt(np.sum(np.square(comm_dis_gps - comm_dis_insar)) / (num_comm_date - 1))
            #print('site: {}, RMSE: {:.1f} cm'.format(self.site_names[i], dis_rmse*100.))


    def sort_by_velocity(ds):
        ## 4. calculate velocity to sort plotting order
        site_vel = {}
        site_names = sorted(list(ds.keys()))
        for sname in site_names:
            site = ds[sname]
            # design matrix
            yr_diff = np.array([i.year + (i.timetuple().tm_yday - 1) / 365.25 for i in site['gps_datetime']])
            yr_diff -= yr_diff[0]
            A = np.ones([len(site['gps_datetime']), 2], dtype=np.float32)
            A[:, 0] = yr_diff
            # LS estimation
            ts = np.array(site['gps_dis'])
            ts -= ts[0]
            X = np.dot(np.linalg.pinv(A), ts)[0]
            site_vel[sname] = X

        site_names2plot = [i[0] for i in sorted(site_vel.items(), key=lambda kv: kv[1], reverse=True)]
        site_names2plot = [i for i in site_names2plot if site_vel[i] != 0]
        return site_names2plot

    def print_stats(ds):
        site_names = sorted(list(ds.keys()))
        for sname in site_names:
            site = ds[sname]
            print('{}, rmse: {:.1f} cm, r_square: {:.2f}, temp_coh: {:.2f}'.format(
                sname,
                site['dis_rmse']*100.,
                site['r_square'],
                site['temp_coh'],
            ))
        return

    def plot_one_site(ax, site, offset=0.):
        # GPS
        ax.errorbar(site['gps_datetime'],
                    site['gps_dis']-offset,
                    yerr=site['gps_std']*3.,
                    ms=marker_size*0.2, lw=0, alpha=1., fmt='-o',
                    elinewidth=edge_width*0.5, ecolor='C0',
                    capsize=marker_size*0.25, markeredgewidth=edge_width*0.5,
                    label='GPS', zorder=1)
        # InSAR
        ecolor = 'gray' if site['temp_coh'] < 0.7 else 'C1'
        insar_dis_name = [i for i in site.keys() if i.startswith('insar_dis')][0]
        ax.scatter(site['insar_datetime'],
                   site[insar_dis_name]-offset,
                   s=5**2, label='InSAR',
                   facecolors='none', edgecolors=ecolor, linewidth=1., alpha=0.7, zorder=2)
        # Label
        ax.annotate('{:.1f} / {:.2f} / {:.2f}'.format(site['dis_rmse']*100., site['r_square'], site['temp_coh']),
                    xy=(1.03, site[insar_dis_name][-1] - offset - 0.02),
                    xycoords=ax.get_yaxis_transform(),  # y in data untis, x in axes fraction
                    color='k', fontsize=font_size)
        ax.annotate('{}'.format(site['name']),
                    xy=(0.05, site[insar_dis_name][0] - offset + 0.1),
                    xycoords=ax.get_yaxis_transform(),  # y in data untis, x in axes fraction
                    color='k', fontsize=font_size)
        return ax

############################## end of insar_vs_gps class ####################################
