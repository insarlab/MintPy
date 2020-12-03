############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#     from mintpy.objects.resample import resample


try:
    import pyresample as pr
except ImportError:
    raise ImportError('Can not import pyresample!')

import numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator as RGI
from mintpy.utils import readfile, ptime, utils0 as ut


class resample:
    """
    Geometry Definition objects for geocoding using:
    1) pyresample (http://pyresample.readthedocs.org)
    2) scipy.interpolate.RegularGridInterpolator:
       (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html)

    Example:
        # prepare resample object
        res_obj = resample(lut_file='./inputs/geometryRadar.h5', src_file='velocity.h5')
        res_obj.open()

        # run geocoding
        rdr_data = readfile.read('velocity.h5')[0]
        geo_data = res_obj.run_resample(src_data=rdr_data, interp_method='nearest', fill_value=np.nan)
    """

    def __init__(self, lut_file, src_file, SNWE=None, lalo_step=None, software='pyresample'):
        """
        Parameters: lut_file  - str, path of lookup table file, containing datasets:
                                latitude / longitude      for lut_file in radar-coord
                                azimuthCoord / rangeCoord for lut_file in geo-coord
                    src_file  - str, path of data file to be resampled.
                    SNWE      - tuple of 4 float, indicating south/north/west/east
                    lalo_step - list of 2 float, output step in lat/lon direction in degree / meter
                    software  - str, interpolation software
        """
        self.lut_file = lut_file
        self.src_file = src_file
        self.SNWE = SNWE
        self.lalo_step = lalo_step
        self.software = software
        self.valid_index = None


    def open(self):
        """Prepare aux data before interpolation operation"""
        self.lut_meta = readfile.read_attribute(self.lut_file)
        self.src_meta = readfile.read_attribute(self.src_file)

        if not self.software:
            if 'Y_FIRST' in self.lut_meta.keys():
                self.software = 'scipy'
            else:
                self.software = 'pyresample'

        # get src_def and dest_def, update SNWE and lalo_step
        if self.software == 'pyresample':
            if 'Y_FIRST' in self.lut_meta.keys():
                # gamma / roipac
                self.prepare_geometry_definition_geo()
            else:
                # isce / doris
                self.prepare_geometry_definition_radar()
            self.length, self.width = self.dest_def.lats.shape

        elif self.software == 'scipy' and 'Y_FIRST' in self.lut_meta.keys():
            # gamma / roipac
            self.prepare_regular_grid_interpolator()


    def run_resample(self, src_data, interp_method='nearest', fill_value=np.nan, nprocs=1, print_msg=True):
        """Run interpolation operation for input 2D/3D data
        Parameters: src_data      : 2D/3D np.array, source data to be geocoded
                    interp_method : string, nearest | linear
                    fill_value    : NaN or number
                    nprocs        : int, number of processes to be used
                    print_msg     : bool
        Returns:    geo_data      : 2D/3D np.array
        """
        # use pyresample
        if self.software == 'pyresample':
            if len(src_data.shape) == 3:
                src_data = np.moveaxis(src_data, 0, -1)

            if src_data.dtype == np.bool_:
                fill_value = False
                print('restrict fill value to False for bool type source data')

            # resample source data into target data
            geo_data = self.run_pyresample(src_data=src_data,
                                           interp_method=interp_method,
                                           fill_value=fill_value,
                                           nprocs=nprocs,
                                           radius=None,
                                           print_msg=True)

            if len(geo_data.shape) == 3:
                geo_data = np.moveaxis(geo_data, -1, 0)

        # use scipy.interpolate.RegularGridInterpolator
        else:
            if print_msg:
                print('{} resampling using scipy.interpolate.RegularGridInterpolator ...'.format(interp_method))
            if len(src_data.shape) == 3:
                geo_data = np.empty((src_data.shape[0], self.length, self.width), src_data.dtype)
                prog_bar = ptime.progressBar(maxValue=src_data.shape[0], print_msg=print_msg)
                for i in range(src_data.shape[0]):
                    geo_data[i, :, :] = self.run_regular_grid_interpolator(src_data=src_data[i, :, :],
                                                                           interp_method=interp_method,
                                                                           fill_value=fill_value,
                                                                           print_msg=True)
                    prog_bar.update(i+1)
                prog_bar.close()
            else:
                geo_data = self.run_regular_grid_interpolator(src_data=src_data,
                                                              interp_method=interp_method,
                                                              fill_value=fill_value,
                                                              print_msg=True)
        return geo_data


    ##### resample using scipy.interpolate
    def prepare_regular_grid_interpolator(self):
        """Prepare aux data for RGI module"""
        # source points in regular grid
        src_length = int(self.src_meta['LENGTH'])
        src_width = int(self.src_meta['WIDTH'])
        self.src_pts = (np.arange(src_length) + 0.5,
                        np.arange(src_width) + 0.5)

        # destination points
        dest_y = readfile.read(self.lut_file, datasetName='azimuthCoord')[0]
        dest_x = readfile.read(self.lut_file, datasetName='rangeCoord')[0]
        if 'SUBSET_XMIN' in self.src_meta.keys():
            print('input data file was cropped before.')
            dest_y[dest_y != 0.] -= float(self.src_meta['SUBSET_YMIN'])
            dest_x[dest_x != 0.] -= float(self.src_meta['SUBSET_XMIN'])
        self.interp_mask = np.multiply(np.multiply(dest_y > 0, dest_y < src_length),
                                       np.multiply(dest_x > 0, dest_x < src_width))
        self.dest_pts = np.hstack((dest_y[self.interp_mask].reshape(-1, 1),
                                   dest_x[self.interp_mask].reshape(-1, 1)))

        # destimation data size
        self.length = int(self.lut_meta['LENGTH'])
        self.width = int(self.lut_meta['WIDTH'])
        lat0 = float(self.lut_meta['Y_FIRST'])
        lon0 = float(self.lut_meta['X_FIRST'])
        lat_step = float(self.lut_meta['Y_STEP'])
        lon_step = float(self.lut_meta['X_STEP'])
        self.lalo_step = (lat_step, lon_step)
        self.SNWE = (lat0 + lat_step * self.length,
                     lat0,
                     lon0,
                     lon0 + lon_step * self.width)


    def run_regular_grid_interpolator(self, src_data, interp_method='nearest', fill_value=np.nan, print_msg=True):
        """Interpolate 2D matrix"""
        # prepare interpolation function
        rgi_func = RGI(self.src_pts,
                       src_data,
                       method=interp_method,
                       bounds_error=False,
                       fill_value=fill_value)

        # prepare output matrix
        geo_data = np.empty((self.length, self.width), src_data.dtype)
        geo_data.fill(fill_value)

        # interpolate output matrix
        geo_data[self.interp_mask] = rgi_func(self.dest_pts)
        return geo_data


    ##### resample using pyresample.kd_tree
    def prepare_geometry_definition_radar(self):
        """Get src_def and dest_def for lookup table from ISCE, DORIS"""

        def mark_lalo_anomoly(lat, lon):
            """mask pixels with abnormal values (0, etc.)
            This is found on sentinelStack multiple swath lookup table file.
            """
            # ignore pixels with zero value
            zero_mask = np.multiply(lat != 0., lon != 0.)

            # ignore anomaly non-zero values 
            # by get the most common data range (d_min, d_max) based on histogram
            mask = np.array(zero_mask, np.bool_)
            for data in [lat, lon]:
                bin_value, bin_edge = np.histogram(data[mask], bins=10)                
                # if there is anomaly, histogram won't be evenly distributed
                while np.max(bin_value) > np.sum(zero_mask) * 0.3:
                    # find the continous bins where the largest bin is --> normal data range
                    bin_value_thres = ut.median_abs_deviation_threshold(bin_value, cutoff=3)
                    bin_label = ndimage.label(bin_value > bin_value_thres)[0]
                    idx = np.where(bin_label == bin_label[np.argmax(bin_value)])[0]
                    # convert to min/max data value
                    bin_step = bin_edge[1] - bin_edge[0]
                    d_min = bin_edge[idx[0]] - bin_step / 2.
                    d_max = bin_edge[idx[-1]+1] + bin_step / 2.
                    mask *= np.multiply(data >= d_min, data <= d_max)
                    bin_value, bin_edge = np.histogram(data[mask], bins=10)
            lat[mask == 0] = 90.
            lon[mask == 0] = 0.
            return lat, lon, mask

        # radar2geo
        if 'Y_FIRST' not in self.src_meta.keys():

            ## source
            # src_lat/lon at pixel center
            src_lat = readfile.read(self.lut_file, datasetName='latitude')[0]
            src_lon = readfile.read(self.lut_file, datasetName='longitude')[0]
            src_lat, src_lon, mask = mark_lalo_anomoly(src_lat, src_lon)

            # SNWE at pixel center
            SNWE = (np.nanmin(src_lat[mask]), np.nanmax(src_lat[mask]),
                    np.nanmin(src_lon[mask]), np.nanmax(src_lon[mask]))

            ## destination
            # 1. lalo_step
            if self.lalo_step is None:
                self.lalo_step = ((SNWE[0] - SNWE[1]) / (src_lat.shape[0] - 1),
                                  (SNWE[3] - SNWE[2]) / (src_lat.shape[1] - 1))
            else:
                # ensure lat/lon step sign
                self.lalo_step = (abs(self.lalo_step[0]) * -1.,
                                  abs(self.lalo_step[1]) * 1.)
            print('output pixel size in (lat, lon) in degree: {}'.format(self.lalo_step))

            # SNWE at pixel center --> SNWE at pixel outer boundary
            # to be consistent with input bounding box
            SNWE = (SNWE[0] + self.lalo_step[0] / 2.0,
                    SNWE[1] - self.lalo_step[0] / 2.0,
                    SNWE[2] - self.lalo_step[1] / 2.0,
                    SNWE[3] + self.lalo_step[1] / 2.0)

            # 2. SNWE at pixel outer boundary
            # consistent with Y/X_FIRST
            if self.SNWE is None:
                self.SNWE = SNWE
            # adjust SNWE ending coordinate (S, E) for precise alignment
            lat_num = int(np.rint((self.SNWE[0] - self.SNWE[1]) / self.lalo_step[0]))
            lon_num = int(np.rint((self.SNWE[3] - self.SNWE[2]) / self.lalo_step[1]))
            self.SNWE = (self.SNWE[1] + self.lalo_step[0] * lat_num,
                         self.SNWE[1],
                         self.SNWE[2],
                         self.SNWE[2] + self.lalo_step[1] * lon_num)
            print('output area extent in (S, N, W, E) in degree: {}'.format(self.SNWE))

            # 3. dest_lat/lon at pixel center
            lat0 = self.SNWE[1] + self.lalo_step[0]/2
            lat1 = self.SNWE[0] - self.lalo_step[0]/2
            lon0 = self.SNWE[2] + self.lalo_step[1]/2
            lon1 = self.SNWE[3] - self.lalo_step[1]/2
            dest_lat, dest_lon = np.mgrid[lat0:lat1:lat_num*1j,
                                          lon0:lon1:lon_num*1j]

            ## reduction of swath data
            # https://pyresample.readthedocs.io/en/latest/data_reduce.html
            src_size_deg = (SNWE[1] - SNWE[0]) * (SNWE[3] - SNWE[2])
            dest_size_deg = (self.SNWE[1] - self.SNWE[0]) * (self.SNWE[3] - self.SNWE[2])
            if dest_size_deg < src_size_deg * 0.5:
                self.valid_index = pr.data_reduce.get_valid_index_from_lonlat_grid(dest_lon,
                                                                                   dest_lat,
                                                                                   src_lon,
                                                                                   src_lat,
                                                                                   radius_of_influence=3000)
                src_lon = src_lon[self.valid_index]
                src_lat = src_lat[self.valid_index]

                # bounding box [can be used to read src data; not used afterwards yet]
                idx_row, idx_col = np.where(self.valid_index)
                self.valid_box = (np.min(idx_col), np.min(idx_row),
                                  np.max(idx_col), np.max(idx_row))

            ## geometry definition
            self.src_def = pr.geometry.SwathDefinition(lons=src_lon, lats=src_lat)
            self.dest_def = pr.geometry.GridDefinition(lons=dest_lon, lats=dest_lat)

        # geo2radar
        else:
            # dest_lat/lon at pixel center
            dest_lat = readfile.read(self.lut_file, datasetName='latitude')[0]
            dest_lon = readfile.read(self.lut_file, datasetName='longitude')[0]
            dest_lat, dest_lon, mask = mark_lalo_anomoly(dest_lat, dest_lon)

            # src_lat/lon at pixel center
            lat_step = float(self.src_meta['Y_STEP'])
            lon_step = float(self.src_meta['X_STEP'])
            lat_num = int(self.src_meta['LENGTH'])
            lon_num = int(self.src_meta['WIDTH'])
            lat0 = float(self.src_meta['Y_FIRST']) + lat_step / 2.0
            lon0 = float(self.src_meta['X_FIRST']) + lon_step / 2.0
            lat1 = lat0 + lat_step * (lat_num - 1)
            lon1 = lon0 + lon_step * (lon_num - 1)

            src_lat, src_lon = np.mgrid[lat0:lat1:lat_num * 1j,
                                        lon0:lon1:lon_num * 1j]

            # geometry definition
            self.src_def = pr.geometry.GridDefinition(lons=src_lon, lats=src_lat)
            self.dest_def = pr.geometry.SwathDefinition(lons=dest_lon, lats=dest_lat)

        return src_def, dest_def, src_valid_index


    ## Currently NOT used by default, as GAMMA/ROI_PAC geocoding is using RGI from scipy
    def prepare_geometry_definition_geo(self):
        """Get src_def and dest_def for lookup table from Gamma and ROI_PAC."""
        def project_yx2lalo(yy, xx, SNWE, laloScale):
            """scale/project coordinates in pixel number into lat/lon
            based on bbox and step
            """
            lats = SNWE[1] + yy * laloScale[0]
            lons = SNWE[2] + xx * laloScale[1]
            lats[np.isnan(lats)] = 90.
            lats[np.isnan(lons)] = 90.
            lons[lats == 90.] = 0
            return lats, lons

        # radar2geo
        if 'Y_FIRST' not in self.src_meta.keys():
            lat_step = float(self.lut_meta['Y_STEP'])
            lon_step = float(self.lut_meta['X_STEP'])
            lat_num = int(self.lut_meta['LENGTH'])
            lon_num = int(self.lut_meta['WIDTH'])
            lat0 = float(self.lut_meta['Y_FIRST']) + lat_step * 0.5
            lon0 = float(self.lut_meta['X_FIRST']) + lon_step * 0.5

            ## lalo_step
            if self.lalo_step is None:
                self.lalo_step = (lat_step, lon_step)
            else:
                self.lalo_step = (abs(self.lalo_step[0]) * -1.,
                                  abs(self.lalo_step[1]) * 1.)
            print('output pixel size in (lat, lon) in degree: {}'.format(self.lalo_step))

            ## SNWE at pixel outer boundary
            if self.SNWE is None:
                self.SNWE = (lat0 + lat_step * (lat_num - 0.5),
                             lat0 - lat_step * 0.5,
                             lon0 - lon_step * 0.5,
                             lon0 + lon_step * (lon_num - 0.5))
            # SNWE --> dest_box
            dest_box = (int(np.rint((self.SNWE[2] + lat_step * 0.5 - lon0) / lon_step)),  # x0 - W
                        int(np.rint((self.SNWE[1] + lat_step * 0.5 - lat0) / lat_step)),  # y0 - N
                        int(np.rint((self.SNWE[3] + lon_step * 0.5 - lon0) / lon_step)),  # x1 - E
                        int(np.rint((self.SNWE[0] + lat_step * 0.5 - lat0) / lat_step)))  # y1 - S
            dest_box = (max(dest_box[0], 0),
                        max(dest_box[1], 0),
                        min(dest_box[2], lon_num),
                        min(dest_box[3], lat_num))
            # dest_box --> SNWE
            # adjust SNWE to precisely align with lookup table in geo-coord, based on dest_box,
            # including the starting coordinates (N, W).
            # This is different/disadvantage from the procedure for lookup table in radar-coord
            self.SNWE = (lat0 + lat_step * (dest_box[3] - 0.5),
                         lat0 + lat_step * (dest_box[1] - 0.5),
                         lon0 + lon_step * (dest_box[0] - 0.5),
                         lon0 + lon_step * (dest_box[2] - 0.5))
            print('output area extent in (S, N, W, E) in degree: {}'.format(self.SNWE))

            # src_y/x at pixel center
            length = int(self.src_meta['LENGTH'])
            width = int(self.src_meta['WIDTH'])
            src_y, src_x = np.mgrid[0.5:length-0.5:length*1j,
                                    0.5:width-0.5:width*1j]

            # dest_y/x at pixel center
            dest_y = readfile.read(self.lut_file, datasetName='azimuthCoord', box=dest_box)[0]
            dest_x = readfile.read(self.lut_file, datasetName='rangeCoord', box=dest_box)[0]
            if 'SUBSET_XMIN' in self.src_meta.keys():
                print('input data file was cropped before.')
                dest_y[dest_y != 0.] -= float(self.src_meta['SUBSET_YMIN'])
                dest_x[dest_x != 0.] -= float(self.src_meta['SUBSET_XMIN'])

            # Convert y/x to lat/lon to be able to use pyresample.geometryDefinition
            commMask = np.multiply(np.multiply(dest_y > 0, dest_y < length),
                                   np.multiply(dest_x > 0, dest_x < width))
            idx_row, idx_col = np.where(commMask)
            commSNWE = (self.SNWE[1] + lat_step * np.max(idx_row),
                        self.SNWE[1] + lat_step * np.min(idx_row),
                        self.SNWE[2] + lon_step * np.min(idx_col),
                        self.SNWE[2] + lon_step * np.max(idx_col))
            laloScale = ((commSNWE[0] - commSNWE[1]) / length,
                         (commSNWE[3] - commSNWE[2]) / width)

            src_lat, src_lon = project_yx2lalo(src_y, src_x, commSNWE, laloScale)
            dest_y[dest_y == 0.] = np.nan
            dest_x[dest_x == 0.] = np.nan
            dest_lat, dest_lon = project_yx2lalo(dest_y, dest_x, commSNWE, laloScale)

            # src_def and dest_def
            self.src_def = pr.geometry.GridDefinition(lons=src_lon, lats=src_lat)
            self.dest_def = pr.geometry.SwathDefinition(lons=dest_lon, lats=dest_lat)

        # geo2radar [not finished]
        else:
            raise ValueError('geo2radar with lookup table in geo-coord it NOT supported yet!')

            # provide dest_lat/dest_lon to get dest_x, dest_y
            dest_y_lt = readfile.read(self.lut_file, datasetName='azimuthCoord', box=dest_box)[0]
            dest_x_lt = readfile.read(self.lut_file, datasetName='rangeCoord', box=dest_box)[0]

            # please check this is right, I am not sure what self is. Just copy the above code
            lat0 = float(self.lut_meta['Y_FIRST'])
            lon0 = float(self.lut_meta['X_FIRST'])
            lat_step = float(self.lut_meta['Y_STEP'])
            lon_step = float(self.lut_meta['X_STEP'])

            yy = int((dest_lat - lat0)/lat_step)
            xx = int((dest_lon - lon0)/lon_step)

            row, col = yy.shape
            yy = yy.flatten()
            xx = xx.flatten()

            dest_y = dest_y_lt[yy,xx] # rows in radar coord
            dest_x = dest_x_lt[yy,xx] # column in radar coord
            
            dest_y = dest_y.reshape(row,col)
            dest_x = dest_x.reshape(row,col)
            # please Yunjun to finish the left part to update self

            # src_y/x
            #raise NotImplementedError('Not implemented yet for GAMMA and ROIPAC products')
        return


    def get_radius_of_influence(self, ratio=3):
        """Get radius of influence based on the lookup table resolution in lat/lon direction"""
        earth_radius = 6371.0e3
        # Get lat/lon resolution/step in meter
        lat_c = np.nanmean(self.dest_def.lats)
        lat_step = self.lalo_step[1] * np.pi/180.0 * earth_radius
        lon_step = self.lalo_step[0] * np.pi/180.0 * earth_radius * np.cos(lat_c * np.pi/180)
        radius = np.max(np.abs([lat_step, lon_step])) * ratio
        return radius


    def get_segment_number(self, unit_size=1e6):
        num_segment = int(self.dest_def.size / unit_size + 0.5)
        return num_segment


    def run_pyresample(self, src_data, interp_method='nearest', fill_value=np.nan, nprocs=1,
                       radius=None, print_msg=True):
        """
        Resample input src_data into dest_data
        Parameters: src_data : 2D np.array
                    interp_method : str,
                    fill_value : number
                    nprocs : int
                    print_msg : bool
        Returns:    dest_data
        Example:    dest_data = reObj.run_pyresample(src_data, src_def, dest_def,
                                                     interp_method=inps.interpMethod,
                                                     fill_value=np.fillValue, nprocs=4)
        """
        if not radius:
            # geo2radar
            if 'Y_FIRST' in self.src_meta.keys():
                radius = 100e3
            # radar2geo
            else:
                radius = self.get_radius_of_influence()

        if np.isnan(fill_value) and src_data.dtype not in [np.float32, np.float64, np.float128,
                                                           np.complex64, np.complex128]:
            fill_value = 0
            if print_msg:
                print('input source data is not float, change fill_value from NaN to 0.')

        # reduction of swath data
        if self.valid_index is not None:
            src_data = src_data[self.valid_index]

        # get number of segments
        num_segment = self.get_segment_number()

        if interp_method.startswith('near'):
            if print_msg:
                msg = '{} resampling with pyresample.kd_tree '.format(interp_method)
                msg += 'using {} CPU cores in {} segments ...'.format(nprocs, num_segment)
                print(msg)
            dest_data = pr.kd_tree.resample_nearest(self.src_def,
                                                    src_data,
                                                    self.dest_def,
                                                    nprocs=nprocs,
                                                    fill_value=fill_value,
                                                    radius_of_influence=radius,
                                                    segments=num_segment,
                                                    epsilon=0.5)

        elif interp_method.endswith('linear'):
            if print_msg:
                msg = '{} resampling with pyresample.bilinear '.format(interp_method)
                msg += 'using {} CPU cores ...'.format(nprocs)
                print(msg)
            dest_data = pr.bilinear.resample_bilinear(src_data,
                                                      self.src_def,
                                                      self.dest_def,
                                                      nprocs=nprocs,
                                                      fill_value=fill_value,
                                                      radius=radius,
                                                      neighbours=32,
                                                      segments=num_segment,
                                                      epsilon=0)

        # for debug
        debug_mode = False
        if debug_mode:
            import matplotlib.pyplot as plt
            fig, ((ax11, ax12, ax13), (ax21, ax22, ax23)) = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
            dest_lats = np.array(self.dest_def.lats); dest_lats[dest_lats == 90.] = np.nan
            dest_lons = np.array(self.dest_def.lons); dest_lons[dest_lons == 0.] = np.nan
            im = ax11.imshow(self.src_def.lats);    fig.colorbar(im, ax=ax11);   ax11.set_title('src_lats')
            im = ax12.imshow(self.src_def.lons);    fig.colorbar(im, ax=ax12);   ax12.set_title('src_lons')
            im = ax13.imshow(src_data);             fig.colorbar(im, ax=ax13);   ax13.set_title('src_data')
            im = ax21.imshow(dest_lats);            fig.colorbar(im, ax=ax21);   ax21.set_title('dest_lats')
            im = ax22.imshow(dest_lons);            fig.colorbar(im, ax=ax22);   ax22.set_title('dest_lons')
            im = ax23.imshow(dest_data);            fig.colorbar(im, ax=ax23);   ax23.set_title('dest_data')
            plt.show()

        return dest_data
