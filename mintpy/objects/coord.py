############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################
# Recommend import:
#   from mintpy.utils import utils as ut


from argparse import Namespace
import numpy as np
from mintpy.utils import readfile
from mintpy.utils.utils0 import *
from mintpy.utils.utils1 import *


#####################################  coordinate class begin ##############################################
class coordinate:
    """
    Coordinates conversion based lookup file in pixel-level accuracy

    Convention:
        radar coordinates or row and column numbers represent the whole pixel, starting from 0 as a python convention
        geo coordinates represent an infinite precise point without size:
            (list) of points indicate the lat/lon of the pixel center.
            (bounding) box of points indicate the lat/lon of the pixel UL corner.

    Example:
        atr = readfile.read('velocity.h5')
        coord = ut.coordinate(atr, lookup_file='inputs/geometryRadar.h5')  # for radar coord file
        coord = ut.coordinate(atr)                                         # for geo   coord file
        y, x = coord.geo2radar(lat, lon)[0:2]
        lat, lon = coord.radar2geo(y, x)[0:2]
    """

    def __init__(self, metadata, lookup_file=None):
        """Define a coordinate object
        Parameters: metadata    - dict, source metadata
                    lookup_file - list of 2 strings, or string, lookup table file(s)
        Example:    from mintpy.utils import readfile, utils as ut
                    atr = readfile.read_attribute('./velocity.h5')
                    coord = ut.coordinate(atr, './inputs/geometryRadar.h5')
                    coord.geo2radar(33.450, -90.22)
                    coord.radar2geo(50, 200)
        """
        self.src_metadata = metadata
        if lookup_file is None:
            lookup_file = get_lookup_file(lookup_file, print_msg=False)
        if isinstance(lookup_file, str):
            lookup_file = [lookup_file, lookup_file]
        self.lookup_file = lookup_file
        self.lut_y = None
        self.lut_x = None

    def open(self):
        try:
            self.earth_radius = float(self.src_metadata['EARTH_RADIUS'])
        except:
            self.earth_radius = 6371.0e3

        if 'Y_FIRST' in self.src_metadata.keys():
            self.geocoded = True
            self.lat0 = float(self.src_metadata['Y_FIRST'])
            self.lon0 = float(self.src_metadata['X_FIRST'])
            self.lat_step = float(self.src_metadata['Y_STEP'])
            self.lon_step = float(self.src_metadata['X_STEP'])
        else:
            self.geocoded = False
            if self.lookup_file:
                self.lut_metadata = readfile.read_attribute(self.lookup_file[0])

    def lalo2yx(self, coord_in, coord_type):
        """convert geo coordinates into radar coordinates for Geocoded file only
        Parameters: geoCoord   - list / tuple / 1D np.ndarray / float, coordinate(s) in latitude or longitude
                    metadata   - dict, dictionary of file attributes
                    coord_type - str, coordinate type: latitude or longitude
        Example:    300 = coordinate.lalo2yx(32.104990,    metadata,'lat')
                    [1000,1500] = coordinate.lalo2yx([130.5,131.4],metadata,'lon')
        """
        self.open()
        if not self.geocoded:
            raise ValueError('Input file is not geocoded.')

        # input format
        if isinstance(coord_in, np.ndarray):
            coord_in = coord_in.tolist()
        if isinstance(coord_in, float):
            coord_in = [coord_in]
        coord_in = list(coord_in)

        # convert coordinates
        coord_type = coord_type.lower()
        coord_out = []
        for i in range(len(coord_in)):
            if coord_type.startswith('lat'):
                coord = int(np.floor((coord_in[i] - self.lat0) / self.lat_step + 0.01)) #plus 0.01 to be more robust in practice
            elif coord_type.startswith('lon'):
                coord = int(np.floor((coord_in[i] - self.lon0) / self.lon_step + 0.01))
            else:
                raise ValueError('Unrecognized coordinate type: '+coord_type)
            coord_out.append(coord)

        # output format
        if len(coord_out) == 1:
            coord_out = coord_out[0]
        elif isinstance(coord_in, tuple):
            coord_out = tuple(coord_out)
        return coord_out


    def yx2lalo(self, coord_in, coord_type):
        """convert radar coordinates into geo coordinates (pixel center)
            for Geocoded file only
        Parameters: coord_in   _ list / tuple / 1D np.ndarray / int, coordinate(s) in row or col in int
                    metadata   _ dict, dictionary of file attributes
                    coord_type _ str, coordinate type: row / col / y / x
        Example:    32.104990 = coord_yx2lalo(300, metadata, 'y')
                    [130.5,131.4] = coord_yx2lalo([1000,1500], metadata, 'x')
        """
        self.open()
        if not self.geocoded:
            raise ValueError('Input file is not geocoded.')

        # Convert to List if input is String
        if isinstance(coord_in, np.ndarray):
            coord_in = coord_in.tolist()
        if isinstance(coord_in, int):
            coord_in = [coord_in]
        coord_in = list(coord_in)

        coord_type = coord_type.lower()
        coord_out = []
        for i in range(len(coord_in)):
            if coord_type.startswith(('row', 'y', 'az', 'azimuth')):
                coord = (coord_in[i] + 0.5) * self.lat_step + self.lat0
            elif coord_type.startswith(('col', 'x', 'rg', 'range')):
                coord = (coord_in[i] + 0.5) * self.lon_step + self.lon0
            else:
                raise ValueError('Unrecognized coordinate type: '+coord_type)
            coord_out.append(coord)

        if len(coord_out) == 1:
            coord_out = coord_out[0]
        elif isinstance(coord_in, tuple):
            coord_out = tuple(coord_out)
        return coord_out


    def _get_lookup_row_col(self, y, x, y_factor=10, x_factor=10, geo_coord=False, debug_mode=False):
        """Get row/col number in y/x value matrix from input y/x
        Use overlap mean value between y and x buffer;
        To support point outside of value pool/matrix, could use np.polyfit to fit a line
        for y and x value buffer and return the intersection point row/col
        """
        ymin = y - y_factor;  ymax = y + y_factor
        xmin = x - x_factor;  xmax = x + x_factor
        if not geo_coord:
            ymin = max(ymin, 0.5)
            xmin = max(xmin, 0.5)
        mask_y = np.multiply(self.lut_y >= ymin, self.lut_y <= ymax)
        mask_x = np.multiply(self.lut_x >= xmin, self.lut_x <= xmax)
        mask_yx = np.multiply(mask_y, mask_x)
        row, col = np.nanmean(np.where(mask_yx), axis=1)

        # for debugging only
        if debug_mode:
            print('Debug mode is ON.\nShow the row/col number searching result.')
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(nrows=1, ncols=3, figsize=[12, 5])
            axs[0].imshow(mask_y);  axs[0].set_title('Buffer in Y direction')
            axs[1].imshow(mask_x);  axs[1].set_title('Buffer in X direction')
            axs[2].imshow(mask_yx); axs[2].set_title('Y & X overlap (zoom in)')

            idx = np.where(np.sum(mask_yx, axis=0))[0]
            idy = np.where(np.sum(mask_yx, axis=1))[0]
            axs[2].set_xlim(idx[0], idx[-1])
            axs[2].set_ylim(idy[0], idy[-1])
            axs[1].set_yticklabels([])
            plt.show()

        # Error message
        if any(np.isnan(i) for i in [row, col]):
            raise RuntimeError('No coresponding coordinate found for y/x: {}/{}'.format(y, x))

        return row, col


    def read_lookup_table(self, print_msg=True):
        if 'Y_FIRST' in self.lut_metadata.keys():
            self.lut_y = readfile.read(self.lookup_file[0],
                                       datasetName='azimuthCoord',
                                       print_msg=print_msg)[0]
            self.lut_x = readfile.read(self.lookup_file[1],
                                       datasetName='rangeCoord',
                                       print_msg=print_msg)[0]
        else:
            self.lut_y = readfile.read(self.lookup_file[0],
                                       datasetName='latitude',
                                       print_msg=print_msg)[0]
            self.lut_x = readfile.read(self.lookup_file[1],
                                       datasetName='longitude',
                                       print_msg=print_msg)[0]
        return self.lut_y, self.lut_x

    def _read_geo_lut_metadata(self):
        """Read lat/lon0, lat/lon_step_deg, lat/lon_step into a Namespace - lut"""
        lut = Namespace()
        lut.lat0 = float(self.lut_metadata['Y_FIRST'])
        lut.lon0 = float(self.lut_metadata['X_FIRST'])
        lut.lat_step_deg = float(self.lut_metadata['Y_STEP'])
        lut.lon_step_deg = float(self.lut_metadata['X_STEP'])

        # Get lat/lon resolution/step in meter
        length = int(self.lut_metadata['LENGTH'])
        lat_c = lut.lat0 + lut.lat_step_deg * length / 2.
        lut.lat_step = lut.lat_step_deg * np.pi/180.0 * self.earth_radius
        lut.lon_step = lut.lon_step_deg * np.pi/180.0 * self.earth_radius * np.cos(lat_c * np.pi/180)
        return lut

    def geo2radar(self, lat, lon, print_msg=True, debug_mode=False):
        """Convert geo coordinates into radar coordinates.
        Parameters: lat/lon   - np.array / float, latitude/longitude
        Returns:    az/rg     - np.array / int, range/azimuth pixel number
                    az/rg_res - float, residul/uncertainty of coordinate conversion
        """
        self.open()
        if self.geocoded:
            az = self.lalo2yx(lat, coord_type='lat')
            rg = self.lalo2yx(lon, coord_type='lon')
            return az, rg, 0, 0

        if not isinstance(lat, np.ndarray):
            lat = np.array(lat)
            lon = np.array(lon)

        # read lookup table
        if self.lookup_file is None:
            raise FileNotFoundError('No lookup table file found!')
        if self.lut_y is None or self.lut_x is None:
            self.read_lookup_table(print_msg=print_msg)

        # For lookup table in geo-coord, read value directly (GAMMA and ROI_PAC)
        if 'Y_FIRST' in self.lut_metadata.keys():
            lut = self._read_geo_lut_metadata()

            # if source data file is subsetted before
            az0 = 0
            rg0 = 0
            if 'SUBSET_YMIN' in self.src_metadata.keys():
                az0 = int(self.src_metadata['SUBSET_YMIN'])
            if 'SUBSET_XMIN' in self.src_metadata.keys():
                rg0 = int(self.src_metadata['SUBSET_XMIN'])

            # uncertainty due to different resolution between src and lut file
            try:
                az_step = azimuth_ground_resolution(self.src_metadata)
                rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
                x_factor = np.ceil(abs(lut.lon_step) / rg_step).astype(int)
                y_factor = np.ceil(abs(lut.lat_step) / az_step).astype(int)
            except:
                x_factor = 2
                y_factor = 2

            # read y/x value from lookup table
            row = np.rint((lat - lut.lat0) / lut.lat_step_deg).astype(int)
            col = np.rint((lon - lut.lon0) / lut.lon_step_deg).astype(int)
            rg = np.rint(self.lut_x[row, col]).astype(int) - rg0
            az = np.rint(self.lut_y[row, col]).astype(int) - az0

        # For lookup table in radar-coord, search the buffer and use center pixel (ISCE)
        else:
            # get resolution in degree in range/azimuth direction
            az_step = azimuth_ground_resolution(self.src_metadata)
            rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
            lat_c = (np.nanmax(lat) + np.nanmin(lat)) / 2.
            az_step_deg = 180./np.pi * az_step / (self.earth_radius)
            rg_step_deg = 180./np.pi * rg_step / (self.earth_radius * np.cos(lat_c * np.pi/180.))

            az, rg = np.zeros(lat.shape), np.zeros(lat.shape)
            x_factor = 10
            y_factor = 10

            # search the overlap area of buffer in x/y direction and use the cross center
            if lat.size == 1:
                az, rg = self._get_lookup_row_col(lat, lon,
                                                  y_factor*az_step_deg,
                                                  x_factor*rg_step_deg,
                                                  geo_coord=True,
                                                  debug_mode=debug_mode)
            else:
                for i in range(rg.size):
                    az[i], rg[i] = self._get_lookup_row_col(lat[i], lon[i],
                                                            y_factor*az_step_deg,
                                                            x_factor*rg_step_deg,
                                                            geo_coord=True,
                                                            debug_mode=debug_mode)
            az = np.floor(az).astype(int)
            rg = np.floor(rg).astype(int)

        rg_resid = x_factor
        az_resid = y_factor
        return az, rg, az_resid, rg_resid


    def radar2geo(self, az, rg, print_msg=True, debug_mode=False):
        """Convert radar coordinates into geo coordinates (pixel center)
        Parameters: rg/az      - np.array / int, range/azimuth pixel number
        Returns:    lon/lat    - np.array / float, longitude/latitude of input point (rg,az);
                                 nan if not found.
                    latlon_res - float, residul/uncertainty of coordinate conversion
        """
        self.open()
        if self.geocoded:
            lat = self.yx2lalo(az, coord_type='az')
            lon = self.yx2lalo(rg, coord_type='rg')
            return lat, lon, 0, 0

        if not isinstance(az, np.ndarray):
            az = np.array(az)
            rg = np.array(rg)

        # read lookup table file
        if self.lookup_file is None:
            raise FileNotFoundError('No lookup table file found!')
        if self.lut_y is None or self.lut_x is None:
            self.read_lookup_table(print_msg=print_msg)

        # For lookup table in geo-coord, search the buffer and use center pixel
        if 'Y_FIRST' in self.lut_metadata.keys():
            lut = self._read_geo_lut_metadata()

            # Get buffer ratio from range/azimuth ground resolution/step
            try:
                az_step = azimuth_ground_resolution(self.src_metadata)
                rg_step = range_ground_resolution(self.src_metadata, print_msg=False)
                x_factor = 2 * np.ceil(abs(lut.lon_step) / rg_step)
                y_factor = 2 * np.ceil(abs(lut.lat_step) / az_step)
            except:
                x_factor = 10
                y_factor = 10

            if 'SUBSET_XMIN' in self.src_metadata.keys():
                rg += int(self.src_metadata['SUBSET_XMIN'])
                az += int(self.src_metadata['SUBSET_YMIN'])

            lut_row = np.zeros(rg.shape)
            lut_col = np.zeros(rg.shape)
            if rg.size == 1:
                lut_row, lut_col = self._get_lookup_row_col(az, rg, y_factor, x_factor,
                                                            debug_mode=debug_mode)
            else:
                for i in range(rg.size):
                    (lut_row[i],
                     lut_col[i]) = self._get_lookup_row_col(az[i], rg[i],
                                                            y_factor, x_factor,
                                                            debug_mode=debug_mode)
            lat = (lut_row + 0.5) * lut.lat_step_deg + lut.lat0
            lon = (lut_col + 0.5) * lut.lon_step_deg + lut.lon0
            lat_resid = abs(y_factor * lut.lat_step_deg)
            lon_resid = abs(x_factor * lut.lon_step_deg)

        # For lookup table in radar-coord, read the value directly.
        else:
            lat = self.lut_y[az, rg]
            lon = self.lut_x[az, rg]
            lat_resid = 0.
            lon_resid = 0.
        return lat, lon, lat_resid, lon_resid

    def box_pixel2geo(self, pixel_box):
        """Convert pixel_box to geo_box in UL corner
        Parameters: pixel_box - list/tuple of 4 int   in (x0, y0, x1, y1)
        Returns:    geo_box   - tuple      of 4 float in (W, N, E, S)
        """
        try:
            lat = self.yx2lalo([pixel_box[1], pixel_box[3]], coord_type='y')
            lon = self.yx2lalo([pixel_box[0], pixel_box[2]], coord_type='x')
            # shift lat from pixel center to the UL corner
            lat = [i - self.lat_step / 2.0 for i in lat]
            lon = [i - self.lon_step / 2.0 for i in lon]
            geo_box = (lon[0], lat[0], lon[1], lat[1])
        except:
            geo_box = None
        return geo_box

    def box_geo2pixel(self, geo_box):
        """Convert geo_box to pixel_box
        Parameters: geo_box   - tuple      of 4 float in (W, N, E, S)
        Returns:    pixel_box - list/tuple of 4 int   in (x0, y0, x1, y1)
        """
        try:
            y = self.lalo2yx([geo_box[1], geo_box[3]], coord_type='latitude')
            x = self.lalo2yx([geo_box[0], geo_box[2]], coord_type='longitude')
            pixel_box = (x[0], y[0], x[1], y[1])
        except:
            pixel_box = None
        return pixel_box

    def bbox_radar2geo(self, pix_box, print_msg=False):
        """Calculate bounding box in lat/lon for file in geo coord, based on input radar/pixel box
        Parameters: pix_box - tuple of 4 int, in (x0, y0, x1, y1)
        Returns:    geo_box - tuple of 4 float, in (W, N, E, S)
        """
        x = np.array([pix_box[0], pix_box[2], pix_box[0], pix_box[2]])
        y = np.array([pix_box[1], pix_box[1], pix_box[3], pix_box[3]])
        lat, lon, lat_res, lon_res = self.radar2geo(y, x, print_msg=print_msg)
        buf = 2*(np.max(np.abs([lat_res, lon_res])))
        geo_box = (np.min(lon) - buf, np.max(lat) + buf,
                   np.max(lon) + buf, np.min(lat) - buf)
        return geo_box

    def bbox_geo2radar(self, geo_box, print_msg=False):
        """Calculate bounding box in x/y for file in radar coord, based on input geo box.
        Parameters: geo_box - tuple of 4 float, indicating the UL/LR lon/lat 
        Returns:    pix_box - tuple of 4 int, indicating the UL/LR x/y of the bounding box in radar coord
                              for the corresponding lat/lon coverage.
        """
        lat = np.array([geo_box[3], geo_box[3], geo_box[1], geo_box[1]])
        lon = np.array([geo_box[0], geo_box[2], geo_box[0], geo_box[2]])
        y, x, y_res, x_res = self.geo2radar(lat, lon, print_msg=print_msg)
        buf = 2 * np.max(np.abs([x_res, y_res]))
        pix_box = (np.min(x) - buf, np.min(y) - buf,
                   np.max(x) + buf, np.max(y) + buf)
        return pix_box

    def check_box_within_data_coverage(self, pixel_box, print_msg=True):
        """Check the subset box's conflict with data coverage
        Parameters:  pixel_box - 4-tuple of int, indicating y/x coordinates of subset
        Returns:     out_box   - 4-tuple of int
        """
        self.open()
        length, width = int(self.src_metadata['LENGTH']), int(self.src_metadata['WIDTH'])
        if pixel_box is None:
            pixel_box = (0, 0, width, length)
        sub_x = [pixel_box[0], pixel_box[2]]
        sub_y = [pixel_box[1], pixel_box[3]]

        if sub_y[0] >= length or sub_y[1] <= 0 or sub_x[0] >= width or sub_x[1] <= 0:
            data_box = (0, 0, width, length)
            msg = 'ERROR: input index is out of data range!\n'
            msg += '\tdata   range in x/y: {}\n'.format(data_box)
            msg += '\tsubset range in x/y: {}\n'.format(pixel_box)
            msg += '\tdata   range in lat/lon: {}\n'.format(self.box_pixel2geo(data_box))
            msg += '\tsubset range in lat/lon: {}\n'.format(self.box_pixel2geo(pixel_box))
            raise ValueError(msg)

        # Check Y/Azimuth/Latitude subset range
        if sub_y[0] < 0:
            sub_y[0] = 0
            if print_msg:
                print('WARNING: input y < min (0)! Set it to min.')
        if sub_y[1] > length:
            sub_y[1] = length
            if print_msg:
                print('WARNING: input y > max ({})! Set it to max.'.format(length))

        # Check X/Range/Longitude subset range
        if sub_x[0] < 0:
            sub_x[0] = 0
            if print_msg:
                print('WARNING: input x < min (0)! Set it to min.')
        if sub_x[1] > width:
            sub_x[1] = width
            if print_msg:
                print('WARNING: input x > max ({})! Set it to max.'.format(width))

        out_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])
        return out_box

#####################################  coordinate class end ##############################################
