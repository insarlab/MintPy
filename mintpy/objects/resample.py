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

import os
import h5py
import numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator as RGI
from mintpy.objects.cluster import split_box2sub_boxes
from mintpy.utils import readfile, ptime, utils0 as ut

EARTH_RADIUS = 6378122.65   # m


class resample:
    """
    Geometry Definition objects for geocoding using:
    1) pyresample (http://pyresample.readthedocs.org)
    2) scipy.interpolate.RegularGridInterpolator:
       (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html)

    Example:
        from mintpy.objects.resample import resample
        from mintpy.utils import readfile, attribute as attr

        ##### opt 1 - entire matrix (by not changing max_memory=0)
        src_file = 'velocity.h5'
        res_obj = resample(lut_file='./inputs/geometryRadar.h5', src_file=src_file)
        res_obj.open()
        res_obj.prepare()

        # read data / attribute
        rdr_data, atr = readfile.read(src_file)
        # resample data
        box = res_obj.src_box_list[0]
        geo_data = res_obj.run_resample(src_data=rdr_data[box[1]:box[3], box[0]:box[2]])
        # update attribute
        atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)

        ##### opt 2 - block-by-block IO (by setting max_memory=4)
        src_file = 'timeseries.h5'
        res_obj = resample(lut_file='./inputs/geometryRadar.h5', src_file=src_file, max_memory=4)
        res_obj.open()
        res_obj.prepare()

        # prepare output: metadata and initial file
        atr = readfile.read_attribute(src_file)
        atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)
        writefile.layout_hdf5(outfile, metadata=atr, ref_file=src_file)

        # block-by-block IO
        for i in range(res_obj.num_bix):
            src_box = res_obj.src_box_list[i]
            dest_box = res_obj.dest_box_list[i]
            # read
            data = readfile.read(src_file, box=src_box)[0]
            # resample
            data = res_obj.run_resample(src_data=rdr_data, box_ind=i)
            # write
            block = [0, data.shape[0],
                     dest_box[1], dest_box[3],
                     dest_box[0], dest_box[2]]
            writefile.write_hdf5_block(outfile,
                                       data=data,
                                       datasetName='timeseries',
                                       block=block)
    """

    def __init__(self, lut_file, src_file=None, SNWE=None, lalo_step=None, interp_method='nearest', fill_value=np.nan,
                 nprocs=1, max_memory=0, software='pyresample', print_msg=True):
        """
        Parameters: lut_file      - str, path of lookup table file, containing datasets:
                                    latitude / longitude      for lut_file in radar-coord
                                    azimuthCoord / rangeCoord for lut_file in geo-coord
                    src_file      - str, path of data file to be resampled.
                    SNWE          - tuple of 4 float, indicating south/north/west/east
                                    coordinates at pixel outer boundary (consistent with Y/X_FIRST; not pixel center)
                    lalo_step     - list of 2 float, output step in lat/lon direction in degree / meter
                                    input lalo_step is used ONLY IF 1) radar2geo = True AND 2) lut_file is in radar-coord
                    interp_method - str, interpolation / resampling method, nearest / linear
                    fill_value    - number, fill value for extrapolation pixels
                    nprocs        - int, number of processors used in parallel
                    max_memory    - float, maximum memory to use
                                    set to 0 or negative value to disable block-by-block IO (default)
                    software      - str, interpolation software, pyresample / scipy
        """
        # input variables
        self.lut_file = lut_file
        self.src_file = src_file
        self.SNWE = SNWE
        self.lalo_step = lalo_step
        self.interp_method = interp_method
        self.fill_value = fill_value

        self.nprocs = nprocs
        self.max_memory = max_memory
        self.software = software
        self.print_msg = print_msg

        # initial variables
        self.num_box = 1
        self.radius = None


    def open(self):
        """Read metadata

        Note: separate open() from prepare() to handle numpy matrix directly w/o src_file
        by assigning src_meta after open() and before prepare()
        """
        # read metadata: lookup table and source data
        if self.lut_file:
            self.lut_meta = readfile.read_attribute(self.lut_file)
        else:
            self.lut_meta = None

        if self.src_file:
            self.src_meta = readfile.read_attribute(self.src_file)
        else:
            self.src_meta = None

        # get num_box
        if self.software == 'pyresample':
            self.num_box = self.get_num_box(self.src_file, self.max_memory)


    def prepare(self):
        """Prepare aux data before interpolation operation"""

        # check metadata: lookup table and source data
        if self.lut_meta is None or self.src_meta is None:
            raise ValueError('lookup table or source data metadata is None!')

        # prepare geometry for resampling
        print('resampling software: {}'.format(self.software))
        if self.software == 'scipy':
            if 'Y_FIRST' in self.lut_meta.keys():
                # gamma / roipac
                self.prepare_regular_grid_interpolator()
            else:
                raise ValueError('resampling using scipy with lookup table in radar-coord (ISCE / DORIS) is NOT supported!')

        elif self.software == 'pyresample':
            if 'Y_FIRST' in self.lut_meta.keys():
                # gamma / roipac
                self.prepare_geometry_definition_geo()

            else:
                # isce / doris
                self.prepare_geometry_definition_radar()

            # aux info
            self.radius = self.get_radius_of_influence(self.lalo_step, self.src_meta)


    def run_resample(self, src_data, box_ind=0, print_msg=True):
        """Run interpolation operation for input 2D/3D data
        Parameters: src_data   - 2D/3D np.array, source data to be resampled
                    box_ind    - int, index of the current box of interest
                                 for multiple boxes with pyresample only
                    print_msg  - bool
        Returns:    dest_data  - 2D/3D np.array, resampled data
        """
        # adjust fill_value for each source data / block
        fill_value = self.fill_value
        float_types = [np.float32, np.float64, np.float128, np.complex64, np.complex128]
        if src_data.dtype == np.bool_:
            fill_value = False
            if print_msg:
                print('input source data is bool type, restrict fill_value to False.')

        elif src_data.dtype not in float_types and np.isnan(fill_value):
            fill_value = 0
            if print_msg:
                print('input source data is NOT float, change fill_value from NaN to 0.')

        ## pyresample
        if self.software == 'pyresample':
            # move 1st/time dimension to the last
            # so that rows/cols axis are the frist, as required by pyresample
            if len(src_data.shape) == 3:
                src_data = np.moveaxis(src_data, 0, -1)

            # resample source data into target data
            dest_data = self.run_pyresample(src_data=src_data,
                                            src_def=self.src_def_list[box_ind],
                                            dest_def=self.dest_def_list[box_ind],
                                            radius=self.radius,
                                            interp_method=self.interp_method,
                                            fill_value=fill_value,
                                            nprocs=self.nprocs,
                                            print_msg=self.print_msg)

            # move 1st/time dimension back
            if len(dest_data.shape) == 3:
                dest_data = np.moveaxis(dest_data, -1, 0)

        ## scipy
        else:
            if print_msg:
                print('{} resampling using scipy.interpolate.RegularGridInterpolator ...'.format(self.interp_method))
            if len(src_data.shape) == 3:
                dest_data = np.empty((src_data.shape[0], self.length, self.width), src_data.dtype)
                prog_bar = ptime.progressBar(maxValue=src_data.shape[0], print_msg=print_msg)
                for i in range(src_data.shape[0]):
                    dest_data[i, :, :] = self.run_regular_grid_interpolator(src_data=src_data[i, :, :],
                                                                            interp_method=self.interp_method,
                                                                            fill_value=fill_value,
                                                                            print_msg=True)
                    prog_bar.update(i+1)
                prog_bar.close()
            else:
                dest_data = self.run_regular_grid_interpolator(src_data=src_data,
                                                               interp_method=self.interp_method,
                                                               fill_value=fill_value,
                                                               print_msg=True)
        return dest_data


    @staticmethod
    def get_num_box(src_file=None, max_memory=0, scale_fac=3.0):
        """Get the number of boxes to split to not exceed the max memory
        Parameters: src_file   - str, path of source data file (largest one if multiple)
                    max_memory - float, memory size in GB
                    scale_fac  - float, scale factor from data size to memory used
                                 empirically estimated.
        Returns:    num_box    - int, number of boxes to be splitted
        """
        num_box = 1

        # auto split into list of boxes ONLY IF:
        # 1. source file is in HDF5 format AND
        # 2. max_memory > 0
        if (src_file and os.path.isfile(src_file)
                and os.path.splitext(src_file)[1] in ['.h5', '.he5'] 
                and max_memory > 0):
            # get max dataset shape
            with h5py.File(src_file, 'r') as f:
                ds_shapes = [f[i].shape for i in f.keys()
                             if isinstance(f[i], h5py.Dataset)]
                max_ds_size = max([np.prod(i) for i in ds_shapes])

            # calc num_box
            num_box = int(np.ceil((max_ds_size * 4 * 4) / (max_memory * 1024**3)))

        return num_box


    @staticmethod
    def get_radius_of_influence(lalo_step, src_meta, ratio=3):
        """Get radius of influence based on the lookup table resolution in lat/lon direction"""
        if 'Y_FIRST' in src_meta.keys():
            # geo2radar
            radius = 100e3
        else:
            # radar2geo
            step_deg = max(np.abs(lalo_step))
            step_m = step_deg * np.pi / 180.0 * EARTH_RADIUS
            radius = step_m * ratio
        return radius


    ##------------------------------ resample using pyresample -------------------------------------##

    def prepare_geometry_definition_radar(self):
        """Get src_def and dest_def for lookup table in radar-coord (from ISCE, DORIS)"""

        def mark_lat_lon_anomoly(lat, lon):
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

            # set invalid pixels to fixed values
            lat[mask == 0] = 90.
            lon[mask == 0] = 0.
            return lat, lon, mask


        # read lookup table: lat/lon at pixel center
        # src  for radar2geo
        # dest for geo2radar
        print('read latitude / longitude from lookup table file: {}'.format(self.lut_file))
        lut_lat = readfile.read(self.lut_file, datasetName='latitude')[0]
        lut_lon = readfile.read(self.lut_file, datasetName='longitude')[0]
        lut_lat, lut_lon, mask = mark_lat_lon_anomoly(lut_lat, lut_lon)

        # radar2geo (with block-by-block support)
        if 'Y_FIRST' not in self.src_meta.keys():

            # src_lat/lon0/1
            src_lat0 = np.nanmax(lut_lat[mask])
            src_lat1 = np.nanmin(lut_lat[mask])
            src_lon0 = np.nanmin(lut_lon[mask])
            src_lon1 = np.nanmax(lut_lon[mask])

            # parameter 1 - lalo_step (output grid)
            if self.lalo_step is None:
                self.lalo_step = ((src_lat1 - src_lat0) / (lut_lat.shape[0] - 1),
                                  (src_lon1 - src_lon0) / (lut_lat.shape[1] - 1))
            else:
                # ensure lat/lon step sign
                self.lalo_step = (abs(self.lalo_step[0]) * -1.,
                                  abs(self.lalo_step[1]) * 1.)
            print('output pixel size in (lat, lon) in degree: {}'.format(self.lalo_step))

            # parameter 2 / 3 - SNWE (at pixel outer boundary; output grid) / length & width
            if self.SNWE is None:
                self.SNWE = (src_lat1 + self.lalo_step[0] / 2.0,
                             src_lat0 - self.lalo_step[0] / 2.0,
                             src_lon0 - self.lalo_step[1] / 2.0,
                             src_lon1 + self.lalo_step[1] / 2.0)
            self.length = int(np.rint((self.SNWE[0] - self.SNWE[1]) / self.lalo_step[0]))
            self.width  = int(np.rint((self.SNWE[3] - self.SNWE[2]) / self.lalo_step[1]))
            # adjust SNWE ending coordinate (S, E) for precise alignment
            self.SNWE = (self.SNWE[1] + self.lalo_step[0] * self.length,
                         self.SNWE[1],
                         self.SNWE[2],
                         self.SNWE[2] + self.lalo_step[1] * self.width)
            print('output area extent in (S, N, W, E) in degree: {}'.format(self.SNWE))
            print('output file row / column number: ({}, {})'.format(self.length, self.width))

            # parameter 4 - list of boxes & geometry definitions
            self.src_box_list = []
            self.src_def_list = []
            self.dest_box_list = []
            self.dest_def_list = []

            # split dest_box (in grid)
            self.dest_box_list = split_box2sub_boxes(box=(0, 0, self.width, self.length),
                                                     num_split=self.num_box,
                                                     dimension='y',
                                                     print_msg=True)

            # dest_box --> src_box / src_def / dest_def
            for i, dest_box in enumerate(self.dest_box_list):
                if self.num_box > 1:
                    print('preparing geometry for dest_box {}/{}: {}'.format(i+1, self.num_box, dest_box))

                # dest_lat/lon at pixel center
                lat_num = dest_box[3] - dest_box[1]
                lon_num = dest_box[2] - dest_box[0]
                lat0 = self.SNWE[1] + self.lalo_step[0] * (dest_box[1] + 0.5)
                lat1 = self.SNWE[1] + self.lalo_step[0] * (dest_box[3] - 0.5)
                lon0 = self.SNWE[2] + self.lalo_step[1] * (dest_box[0] + 0.5)
                lon1 = self.SNWE[2] + self.lalo_step[1] * (dest_box[2] - 0.5)
                dest_lat, dest_lon = np.mgrid[lat0:lat1:lat_num*1j,
                                              lon0:lon1:lon_num*1j]

                # src_box
                src_area = (src_lat1 - src_lat0) * (src_lon1 - src_lon0)
                dest_area = (lat1 - lat0) * (lon1 - lon0)
                if dest_area < src_area * 0.5:
                    # reduction of swath data
                    # https://pyresample.readthedocs.io/en/latest/data_reduce.html
                    # get src_box (in swath) from lat/lon (from dest_box in grid)
                    print('searching relevant box covering the current SNWE')
                    flag = pr.data_reduce.get_valid_index_from_lonlat_grid(dest_lon,
                                                                           dest_lat,
                                                                           lut_lon,
                                                                           lut_lat,
                                                                           radius_of_influence=3000)
                    idx_row, idx_col = np.where(flag)
                    src_box = (np.min(idx_col), np.min(idx_row),
                               np.max(idx_col), np.max(idx_row))
                else:
                    src_box = (0, 0, lut_lat.shape[1], lut_lat.shape[0])

                # geometry definition
                src_def = pr.geometry.SwathDefinition(lons=lut_lon[src_box[1]:src_box[3],
                                                                   src_box[0]:src_box[2]],
                                                      lats=lut_lat[src_box[1]:src_box[3],
                                                                   src_box[0]:src_box[2]])
                dest_def = pr.geometry.GridDefinition(lons=dest_lon, lats=dest_lat)

                self.src_box_list.append(src_box)
                self.src_def_list.append(src_def)
                self.dest_def_list.append(dest_def)


        # geo2radar (WITHOUT block-by-block support)
        else:
            # parameter 1 - lalo_step (input grid)
            self.lalo_step = [float(self.src_meta['Y_STEP']),
                              float(self.src_meta['X_STEP'])]
            print('input pixel size in (lat, lon) in degree: {}'.format(self.lalo_step))

            # parameter 2 - SNWE (input grid)
            lat0 = float(self.src_meta['Y_FIRST'])
            lon0 = float(self.src_meta['X_FIRST'])
            if not self.SNWE:
                # default SNWE --> src_box
                src_box = (0, 0, int(self.src_meta['WIDTH']), int(self.src_meta['LENGTH']))
            else:
                # custom input SNWE --> src_box
                # to align SNWE to precisely to source file in geo-coord
                src_box = (int(np.rint((self.SNWE[2] - lon0) / self.lalo_step[1])),  # x0 - W
                           int(np.rint((self.SNWE[1] - lat0) / self.lalo_step[0])),  # y0 - N
                           int(np.rint((self.SNWE[3] - lon0) / self.lalo_step[1])),  # x1 - E
                           int(np.rint((self.SNWE[0] - lat0) / self.lalo_step[0])))  # y1 - S
            # src_box --> SNWE
            self.SNWE = (lat0 + self.lalo_step[0] * src_box[3],  # S - y1
                         lat0 + self.lalo_step[0] * src_box[1],  # N - y0
                         lon0 + self.lalo_step[1] * src_box[0],  # W - x0
                         lon0 + self.lalo_step[1] * src_box[2])  # E - x1
            print('input area extent in (S, N, W, E) in degree: {}'.format(self.SNWE))

            # parameter 3 - length / width (output grid)
            self.length, self.width = lut_lat.shape

            # src_lat/lon (at pixel center)
            src_len = src_box[3] - src_box[1]
            src_wid = src_box[2] - src_box[0]
            src_lat0 = self.SNWE[1] + self.lalo_step[0] * (src_box[1] + 0.5)
            src_lat1 = self.SNWE[1] + self.lalo_step[0] * (src_box[3] - 0.5)
            src_lon0 = self.SNWE[2] + self.lalo_step[1] * (src_box[0] + 0.5)
            src_lon1 = self.SNWE[2] + self.lalo_step[1] * (src_box[2] - 0.5)
            src_lat, src_lon = np.mgrid[src_lat0:src_lat1:src_len*1j,
                                        src_lon0:src_lon1:src_wid*1j]

            # parameter 4 - list of boxes & geometry definitions
            self.src_box_list = [src_box]
            self.src_def_list = [pr.geometry.GridDefinition(lons=src_lon, lats=src_lat)]
            self.dest_box_list = [(0, 0, self.width, self.length)]
            self.dest_def_list = [pr.geometry.SwathDefinition(lons=lut_lon, lats=lut_lat)]
            self.num_box = 1

        return


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
        # block-by-block support is NOT implemented because:
        # writefile.write_hdf5_block requires split in dest_box for writing
        # while data_reduce.get_valid_index_from_lonlat_grid() requires split in grid (src_box)
        #   resulting in uneven / overlapped flag matrix in dest_box
        if 'Y_FIRST' not in self.src_meta.keys():

            # parameter 1 - lalo_step (output grid)
            # custom input lalo_step is not supported
            # however, it could be added by rounding the input step to the nearest
            # multiple of Y/X_STEP
            self.lalo_step = (float(self.lut_meta['Y_STEP']),
                              float(self.lut_meta['X_STEP']))
            print('output pixel size in (lat, lon) in degree: {}'.format(self.lalo_step))

            # parameter 2 - SNWE (output grid)
            lat0 = float(self.lut_meta['Y_FIRST'])
            lon0 = float(self.lut_meta['X_FIRST'])
            if not self.SNWE:
                # default SNWE --> dest_box
                dest_box = (0, 0, int(self.lut_meta['WIDTH']), int(self.lut_meta['LENGTH']))
            else:
                # custom input SNWE --> dest_box
                # to align SNWE to precisely to lookup table file in geo-coord
                dest_box = (int(np.rint((self.SNWE[2] - lon0) / self.lalo_step[1])),  # x0 - W
                            int(np.rint((self.SNWE[1] - lat0) / self.lalo_step[0])),  # y0 - N
                            int(np.rint((self.SNWE[3] - lon0) / self.lalo_step[1])),  # x1 - E
                            int(np.rint((self.SNWE[0] - lat0) / self.lalo_step[0])))  # y1 - S
                # check coverage box
                dest_box = (max(dest_box[0], 0),
                            max(dest_box[1], 0),
                            min(dest_box[2], int(self.lut_meta['WIDTH'])),
                            min(dest_box[3], int(self.lut_meta['LENGTH'])))
            # dest_box --> SNWE
            self.SNWE = (lat0 + self.lalo_step[0] * dest_box[3],  # S - y1
                         lat0 + self.lalo_step[0] * dest_box[1],  # N - y0
                         lon0 + self.lalo_step[1] * dest_box[0],  # W - x0
                         lon0 + self.lalo_step[1] * dest_box[2])  # E - x1
            print('output area extent in (S, N, W, E) in degree: {}'.format(self.SNWE))

            # parameter 3 - length / width (output grid)
            self.length = dest_box[3] - dest_box[1]
            self.width = dest_box[2] - dest_box[0]

            # parameter 4 - list of boxes & geometry definitions

            # src_y/x at pixel center
            src_len = int(self.src_meta['LENGTH'])
            src_wid = int(self.src_meta['WIDTH'])
            src_y, src_x = np.mgrid[0.5:src_len-0.5:src_len*1j,
                                    0.5:src_wid-0.5:src_wid*1j]

            # dest_y/x at pixel center
            dest_y = readfile.read(self.lut_file, datasetName='azimuthCoord', box=dest_box)[0]
            dest_x = readfile.read(self.lut_file, datasetName='rangeCoord', box=dest_box)[0]
            if 'SUBSET_XMIN' in self.src_meta.keys():
                print('input data file was cropped before.')
                dest_y[dest_y != 0.] -= float(self.src_meta['SUBSET_YMIN'])
                dest_x[dest_x != 0.] -= float(self.src_meta['SUBSET_XMIN'])

            # src/desc_y/x --> src/dest_lat/lon
            # scale y/x to the range of lat/lon to be able to use pyresample.geometryDefinition
            commMask = np.multiply(np.multiply(dest_y > 0, dest_y < src_len),
                                   np.multiply(dest_x > 0, dest_x < src_wid))
            idx_row, idx_col = np.where(commMask)
            commSNWE = (self.SNWE[1] + self.lalo_step[0] * np.max(idx_row),
                        self.SNWE[1] + self.lalo_step[0] * np.min(idx_row),
                        self.SNWE[2] + self.lalo_step[1] * np.min(idx_col),
                        self.SNWE[2] + self.lalo_step[1] * np.max(idx_col))
            laloScale = ((commSNWE[0] - commSNWE[1]) / src_len,
                         (commSNWE[3] - commSNWE[2]) / src_wid)

            src_lat, src_lon = project_yx2lalo(src_y, src_x, commSNWE, laloScale)
            dest_y[dest_y == 0.] = np.nan
            dest_x[dest_x == 0.] = np.nan
            dest_lat, dest_lon = project_yx2lalo(dest_y, dest_x, commSNWE, laloScale)

            # src_def and dest_def
            self.src_box_list = [(0, 0, src_wid, src_len)]
            self.src_def_list = [pr.geometry.GridDefinition(lons=src_lon, lats=src_lat)]
            self.dest_box_list = [dest_box]
            self.dest_def_list = [pr.geometry.SwathDefinition(lons=dest_lon, lats=dest_lat)]
            self.num_box = 1


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


    @staticmethod
    def run_pyresample(src_data, src_def, dest_def, radius, interp_method='nearest', fill_value=np.nan,
                       nprocs=1, print_msg=True):
        """Resample input src_data into dest_data

        Parameters: src_data      - 2/3D np.ndarray, source data to be resampled
                    src_def       - pyresample geometry definition for source data
                    dest_def      - pyresample geometry definition for destination data
                    radius        - float, radius of influence
                    interp_method - str, interpolation method, nearest / linear
                    fill_value    - number, fill value for extrapolation
                    nprocs        - int, number of processors
        Returns:    dest_data     - 2/3D np.ndarray, source data to be resampled
        Example:    dest_data = reObj.run_pyresample(src_data, src_def, dest_def, radius,
                                                     interp_method=inps.interpMethod,
                                                     fill_value=inps.fillValue)
        """
        # get number of segments
        num_segment = int(np.ceil(dest_def.size / 1e6 + 0.5))

        # resample
        if interp_method.startswith('near'):
            if print_msg:
                msg = '{} resampling with pyresample.kd_tree '.format(interp_method)
                msg += 'using {} CPU cores in {} segments ...'.format(nprocs, num_segment)
                print(msg)
            dest_data = pr.kd_tree.resample_nearest(src_def,
                                                    src_data,
                                                    dest_def,
                                                    radius_of_influence=radius,
                                                    fill_value=fill_value,
                                                    nprocs=nprocs,
                                                    segments=num_segment,
                                                    epsilon=0.5)

        elif interp_method.endswith('linear'):
            if print_msg:
                msg = '{} resampling with pyresample.bilinear '.format(interp_method)
                msg += 'using {} CPU cores ...'.format(nprocs)
                print(msg)
            dest_data = pr.bilinear.resample_bilinear(src_data,
                                                      src_def,
                                                      dest_def,
                                                      radius=radius,
                                                      fill_value=fill_value,
                                                      neighbours=32,
                                                      nprocs=nprocs,
                                                      segments=num_segment,
                                                      epsilon=0)

        # for debug
        debug_mode = False
        if debug_mode:
            import matplotlib.pyplot as plt
            fig, ((ax11, ax12, ax13), (ax21, ax22, ax23)) = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
            dest_lats = np.array(dest_def.lats);    dest_lats[dest_lats == 90.] = np.nan
            dest_lons = np.array(dest_def.lons);    dest_lons[dest_lons == 0.] = np.nan
            im = ax11.imshow(src_def.lats);         fig.colorbar(im, ax=ax11);   ax11.set_title('src_lats')
            im = ax12.imshow(src_def.lons);         fig.colorbar(im, ax=ax12);   ax12.set_title('src_lons')
            im = ax13.imshow(src_data);             fig.colorbar(im, ax=ax13);   ax13.set_title('src_data')
            im = ax21.imshow(dest_lats);            fig.colorbar(im, ax=ax21);   ax21.set_title('dest_lats')
            im = ax22.imshow(dest_lons);            fig.colorbar(im, ax=ax22);   ax22.set_title('dest_lons')
            im = ax23.imshow(dest_data);            fig.colorbar(im, ax=ax23);   ax23.set_title('dest_data')
            plt.show()

        return dest_data


    ##------------------------------ resample using scipy.interpolate ------------------------------##

    def prepare_regular_grid_interpolator(self):
        """Prepare aux data for RGI module"""

        # radar2geo
        if 'Y_FIRST' not in self.src_meta.keys():
            # source points in regular grid
            src_len = int(self.src_meta['LENGTH'])
            src_wid = int(self.src_meta['WIDTH'])
            self.src_pts = (np.arange(src_len) + 0.5,
                            np.arange(src_wid) + 0.5)

            # destination points
            dest_y = readfile.read(self.lut_file, datasetName='azimuthCoord')[0]
            dest_x = readfile.read(self.lut_file, datasetName='rangeCoord')[0]
            if 'SUBSET_XMIN' in self.src_meta.keys():
                print('input data file was cropped before.')
                dest_y[dest_y != 0.] -= float(self.src_meta['SUBSET_YMIN'])
                dest_x[dest_x != 0.] -= float(self.src_meta['SUBSET_XMIN'])
            self.interp_mask = np.multiply(np.multiply(dest_y > 0, dest_y < src_len),
                                           np.multiply(dest_x > 0, dest_x < src_wid))
            self.dest_pts = np.hstack((dest_y[self.interp_mask].reshape(-1, 1),
                                       dest_x[self.interp_mask].reshape(-1, 1)))

            # parameter 1/2/3 - lalo_step / SNWE / length & width
            lat_num = int(self.lut_meta['LENGTH'])
            lon_num = int(self.lut_meta['WIDTH'])
            lat_step = float(self.lut_meta['Y_STEP'])
            lon_step = float(self.lut_meta['X_STEP'])
            lat0 = float(self.lut_meta['Y_FIRST'])
            lon0 = float(self.lut_meta['X_FIRST'])
            self.lalo_step = (lat_step, lon_step)
            self.SNWE = (lat0 + lat_step * lat_num,
                         lat0,
                         lon0,
                         lon0 + lon_step * lon_num)
            self.length = lat_num
            self.width = lon_num

            # parameter 4 - list of boxes
            self.src_box_list = [(0, 0, src_wid, src_len)]
            self.dest_box_list = [(0, 0, lon_num, lat_num)]

            # WITHOUT block-by-block support
            # however, date-by-date resampling is already used
            # one could add date-by-date IO in geocoded.py to memory efficiency


        # geo2radar
        else:
            raise ValueError('geo2radar with lookup table in geo-coord it NOT supported yet!')


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
