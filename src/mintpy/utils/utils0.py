"""Miscellaneous utilities - independent."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Contents
#   InSAR
#   File Operation
#   Coordinate
#   Orbit
#   Geometry
#   Image Processing
#   User Interaction
#   Math / Statistics
# Recommend import:
#   from mintpy.utils import utils as ut


import math
import os

import h5py
import numpy as np

# global variables
SPEED_OF_LIGHT = 299792458  # m/s
EARTH_RADIUS = 6371.0088e3  # Earth radius in meters
K = 40.31                   # m^3/s^2, constant


#################################### InSAR ##########################################
def misregistration2coherence(mu):
    """Calculate the resulting coherence due to mis-registration (coregistration error).

    Reference:
      Equation (30) in Just and Bamler (1994); Equation (4.4.27) in Hanssen (2001).

    Parameters: mu  - float / np.ndarray, mis-registration in the unit of resolution cell
                      NOTE: the unit is resolution cell, NOT pixel size/spacing.
                      Applicable to both SAR range and azimuth directions.
    Returns:    coh - float / np.ndarray, coherence.
    """

    # https://numpy.org/doc/stable/reference/generated/numpy.sinc.html
    coh = np.sinc(mu)

    # for coregistration errors >1, set coherence to zero
    if isinstance(mu, np.ndarray):
        coh[mu > 1] = 0
    else:
        if mu > 1:
            coh = 0

    return coh


def range_distance(atr, dimension=2, print_msg=True):
    """Calculate slant range distance from input attribute dict
    Parameters: atr : dict, including the following ROI_PAC attributes:
                    STARTING_RANGE
                    RANGE_PIXEL_SIZE
                    LENGTH
                    WIDTH
                dimension : int, choices = [0,1,2]
                    2 for 2d matrix, vary in range direction, constant in az direction,
                        for radar coord only
                    1 for 1d matrix, in range direction, for radar coord file
                    0 for center value
    Returns:    np.array (0, 1 or 2 D) : range distance between antenna and ground target in meters
    """
    # return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        if print_msg:
            print('input file is geocoded, return center range distance for the whole area')

    range_n, dR = float(atr['STARTING_RANGE']), float(atr['RANGE_PIXEL_SIZE'])
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])

    range_f = range_n + dR*(width-1)
    range_c = (range_f + range_n)/2.0
    if print_msg:
        print('near   range : %.2f m' % (range_n))
        print('center range : %.2f m' % (range_c))
        print('far    range : %.2f m' % (range_f))

    if dimension == 0:
        return np.array(range_c, np.float32)

    range_x = np.linspace(range_n, range_f, num=width)
    if dimension == 1:
        return np.array(range_x, np.float32)
    else:
        range_xy = np.tile(range_x, (length, 1))
        return np.array(range_xy, np.float32)


def incidence_angle(atr, dem=None, dimension=2, print_msg=True):
    """Calculate 2D matrix of incidence angle from ROI_PAC attributes, very accurate.
    Parameters: atr : dict - ROI_PAC attributes including the following items:
                     STARTING_RANGE
                     RANGE_PIXEL_SIZE
                     EARTH_RADIUS
                     HEIGHT
                     WIDTH
                     LENGTH     #for dimension=2
                dem : 2D array for height to calculate local incidence angle
                dimension : int,
                            2 for 2d matrix
                            1 for 1d array
                            0 for one center value
                print_msg : bool
    Returns:    inc_angle : 2D np.array, incidence angle in degree for each pixel
    Example:    dem = readfile.read('hgt.rdr')[0]
                atr = readfile.read_attribute('filt_fine.unw')
                inc_angle = ut.incidence_angle(atr, dem=dem)
    """
    # Return center value for geocoded input file
    if 'Y_FIRST' in atr.keys() and dimension > 0:
        dimension = 0
        if print_msg:
            print('input file is geocoded, return center incident angle only')

    # Read Attributes
    range_n = float(atr['STARTING_RANGE'])
    dR = float(atr['RANGE_PIXEL_SIZE'])
    r = float(atr.get('EARTH_RADIUS', EARTH_RADIUS))
    H = float(atr['HEIGHT'])
    width = int(atr['WIDTH'])

    # Calculation
    range_f = range_n+dR*width
    inc_angle_n = (np.pi - np.arccos((r**2 + range_n**2 - (r+H)**2)/(2*r*range_n))) * 180.0/np.pi
    inc_angle_f = (np.pi - np.arccos((r**2 + range_f**2 - (r+H)**2)/(2*r*range_f))) * 180.0/np.pi
    inc_angle_c = (inc_angle_n + inc_angle_f) / 2.0
    if print_msg:
        print(f'near   incidence angle : {inc_angle_n:.4f} degree')
        print(f'center incidence angle : {inc_angle_c:.4f} degree')
        print(f'far    incidence angle : {inc_angle_f:.4f} degree')

    if dimension == 0:
        inc_angle = inc_angle_c

    elif dimension == 1:
        inc_angle = np.linspace(inc_angle_n, inc_angle_f, num=width,
                                endpoint='FALSE', dtype=np.float32)

    elif dimension == 2:
        length = int(atr['LENGTH'])

        # consider the local variable due to topography
        if dem is not None:
            range_dist = range_distance(atr, dimension=2, print_msg=False)
            inc_angle = (np.pi - np.arccos(((r+dem)**2 + range_dist**2 - (r+H)**2) /
                                           (2*(r+dem)*range_dist))) * 180.0/np.pi
        else:
            inc_angle = np.tile(np.linspace(inc_angle_n, inc_angle_f, num=width,
                                            endpoint='FALSE', dtype=np.float32), (length, 1))
    else:
        raise Exception(f'un-supported dimension input: {dimension}')
    return inc_angle


def incidence_angle2slant_range_distance(atr, inc_angle):
    """Calculate the corresponding slant range distance given an incidence angle

    Law of sines:
               r + H                   r               range_dist
       --------------------- = ----------------- = ------------------ = 2R
        sin(pi - inc_angle)     sin(look_angle)     sin(range_angle)

    where range_angle = inc_angle - look_angle
          R is the radius of the circumcircle.

    link: http://www.ambrsoft.com/TrigoCalc/Triangle/BasicLaw/BasicTriangle.htm

    Parameters: atr         - dict, metadata including the following items:
                                  EARTH_RADIUS
                                  HEIGHT
                inc_angle   - float / np.ndarray, incidence angle in degree
    Returns:    slant_range - float, slant range distance
    """
    if isinstance(inc_angle, str):
        inc_angle = float(inc_angle)
    inc_angle = np.array(inc_angle, dtype=np.float32) / 180 * np.pi
    r = float(atr['EARTH_RADIUS'])
    H = float(atr['HEIGHT'])

    # calculate 2R based on the law of sines
    R2 = (r + H) / np.sin(np.pi - inc_angle)

    look_angle = np.arcsin( r / R2 )
    range_angle = inc_angle - look_angle
    range_dist = R2 * np.sin(range_angle)

    return range_dist


def range_ground_resolution(atr, print_msg=False):
    """Get range resolution on the ground in meters,
        from ROI_PAC attributes, for file in radar coord
    """
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no range resolution info.')
        return
    inc_angle = incidence_angle(atr, dimension=0, print_msg=print_msg)
    rg_step = float(atr['RANGE_PIXEL_SIZE'])/np.sin(inc_angle/180.0*np.pi)
    return rg_step


def azimuth_ground_resolution(atr):
    """Get azimuth resolution on the ground in meters,
        from ROI_PAC attributes, for file in radar coord
    """
    if 'X_FIRST' in atr.keys():
        print('Input file is in geo coord, no azimuth resolution info.')
        return

    proc = atr.get('PROCESSOR', 'isce')
    if proc in ['roipac', 'isce']:
        Re = float(atr['EARTH_RADIUS'])
        height = float(atr['HEIGHT'])
        az_step = float(atr['AZIMUTH_PIXEL_SIZE']) * Re / (Re + height)
    elif proc == 'gamma':
        az_step = float(atr['AZIMUTH_PIXEL_SIZE'])

    return az_step



#################################### File Operation ##########################################
def touch(fname_list, times=None):
    """python equivalent function to Unix utility - touch
    It sets the modification and access times of files to the current time of day.
    If the file doesn't exist, it is created with default permissions.
    Inputs/Output:
        fname_list - string / list of string
    """
    if not fname_list:
        return None

    if isinstance(fname_list, str):
        fname_list = [fname_list]

    fname_list = [x for x in fname_list if x is not None]
    for fname in fname_list:
        if os.path.isfile(fname):
            with open(fname, 'a'):
                os.utime(fname, times)
                print('touch '+fname)

    if len(fname_list) == 1:
        fname_list = fname_list[0]
    return fname_list



################################## Coordinate ##########################################
def utm_zone2epsg_code(utm_zone):
    """Convert UTM Zone string to EPSG code.

    Parameters: utm_zone  - str, atr['UTM_ZONE'], comprises a zone number
                            and a hemisphere, e.g. 11N, 36S, etc.
    Returns:    epsg_code - str, EPSG code
    Examples:   epsg_code = utm_zone2epsg_code('11N')
    """
    from pyproj import CRS
    crs = CRS.from_dict({
        'proj': 'utm',
        'zone': int(utm_zone[:-1]),
        'south': utm_zone[-1].upper() == 'S',
    })
    epsg_code = crs.to_authority()[1]
    return epsg_code


def to_latlon(infile, x, y):
    """Convert x, y in the projection coordinates of the file to lat/lon in degree.

    Similar functionality also exists in utm.to_latlon() at:
        https://github.com/Turbo87/utm#utm-to-latitudelongitude

    Parameters: infile - str, GDAL supported file path
                x/y    - scalar or 1/2D np.ndarray, coordinates in x and y direction
    Returns:    y/x    - scalar or 1/2D np.ndarray, coordinates in latitutde and longitude
    """
    from osgeo import gdal
    from pyproj import Proj, Transformer

    # read projection info using gdal
    ds = gdal.Open(infile)
    srs = ds.GetSpatialRef()

    # if input file is already in lat/lon, do nothing and return
    if (not srs.IsProjected()) and (srs.GetAttrValue('unit') == 'degree'):
        return y, x

    # convert coordinates using pyproj
    # note that Transform.from_proj(x, y, always_xy=True) convert the x, y to lon, lat
    p_in = Proj(ds.GetProjection())
    p_out = Proj('epsg:4326')
    transformer = Transformer.from_proj(p_in, p_out)
    y, x = transformer.transform(x, y)
    return y, x


def utm2latlon(meta, easting, northing):
    """Convert UTM easting/northing in meters to lat/lon in degrees.

    Parameters: meta     - dict, mintpy attributes that includes:
                           UTM_ZONE
                easting  - scalar or 1/2D np.ndarray, UTM    coordinates in x direction
                northing - scalar or 1/2D np.ndarray, UTM    coordinates in y direction
    Returns:    lat      - scalar or 1/2D np.ndarray, WGS 84 coordinates in y direction
                lon      - scalar or 1/2D np.ndarray, WGS 84 coordinates in x direction
    """
    import utm
    zone_num = int(meta['UTM_ZONE'][:-1])
    northern = meta['UTM_ZONE'][-1].upper() == 'N'
    # set 'strict=False' to allow coordinates outside the range of a typical single UTM zone,
    # which can be common for large area analysis, e.g. the Norwegian mapping authority
    # publishes a height data in UTM zone 33 coordinates for the whole country, even though
    # most of it is technically outside zone 33.
    lat, lon = utm.to_latlon(easting, northing, zone_num, northern=northern, strict=False)
    return lat, lon


def latlon2utm(lat, lon):
    """Convert latitude/longitude in degrees to UTM easting/northing in meters.

    Parameters: lat      - scalar or 1/2D np.ndarray, WGS 84 coordinates in y direction
                lon      - scalar or 1/2D np.ndarray, WGS 84 coordinates in x direction
    Returns:    easting  - scalar or 1/2D np.ndarray, UTM    coordinates in x direction
                northing - scalar or 1/2D np.ndarray, UTM    coordinates in y direction
    """
    import utm
    return utm.from_latlon(lat, lon)[:2]


def snwe_to_wkt_polygon(snwe):
    """Convert the input bounding box in SNWE into WKT format POLYGON.

    Parameters: snwe    - list of 4 float, south, north, west and east in degrees/meters
    Returns:    polygon - str, WKT format POLYGON
    """
    S, N, W, E = snwe
    lats = [N, N, S, S, N]
    lons = [W, E, E, W, W]
    polygon = "POLYGON((" + ",".join([f"{lon} {lat}" for lon, lat in zip(lons, lats)])  + "))"
    return polygon


def get_lat_lon(meta, geom_file=None, box=None, dimension=2, ystep=1, xstep=1):
    """Extract precise pixel-wise lat/lon.

    For meta dict in geo-coordinates OR geom_file with latitude/longitude dataset

    Returned lat/lon are corresponds to the pixel center

    Parameters: meta      - dict, including LENGTH, WIDTH and Y/X_FIRST/STEP
                box       - 4-tuple of int for (x0, y0, x1, y1)
                dimension - int, output lat/lon matrix dimension, 1 or 2
                y/xstep   - int, number of pixels to skip for each output pixel
    Returns:    lats      - 1/2D np.array for latitude  in size of (length, _width_)
                lons      - 1/2D np.array for longitude in size of (_length_, width)
    """
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    if box is None:
        box = (0, 0, width, length)

    ds_list = []
    if geom_file is not None:
        with h5py.File(geom_file, 'r') as f:
            ds_list = list(f.keys())

    if 'latitude' in ds_list:
        # read 2D matrices from geometry file
        with h5py.File(geom_file, 'r') as f:
            lats = f['latitude'][box[1]:box[3], box[0]:box[2]]
            lons = f['longitude'][box[1]:box[3], box[0]:box[2]]

    elif 'Y_FIRST' in meta.keys():
        # get lat/lon0/1
        lat_step = float(meta['Y_STEP'])
        lon_step = float(meta['X_STEP'])
        lat0 = float(meta['Y_FIRST']) + lat_step * (box[1] + 0.5)
        lon0 = float(meta['X_FIRST']) + lon_step * (box[0] + 0.5)

        lat_num = box[3] - box[1]
        lon_num = box[2] - box[0]

        lat1 = lat0 + lat_step * (lat_num - 1)
        lon1 = lon0 + lon_step * (lon_num - 1)

        # get matrix of lat/lon
        if dimension == 2:
            lats, lons = np.mgrid[lat0:lat1:lat_num*1j,
                                  lon0:lon1:lon_num*1j]

        elif dimension == 1:
            lats = np.linspace(lat0, lat1, num=lat_num, endpoint=True)
            lons = np.linspace(lon0, lon1, num=lon_num, endpoint=True)

        else:
            raise ValueError(f'un-supported dimension = {dimension}')

        # UTM to lat/lon
        if not meta['Y_UNIT'].startswith('deg') and 'UTM_ZONE' in meta.keys():
            print('UTM coordinates detected, convert UTM into lat/lon')
            lats, lons = utm2latlon(meta, easting=lons, northing=lats)

    else:
        msg = 'Can not get pixel-wise lat/lon!'
        msg += '\nmeta dict is not geocoded and/or geometry file does not contains latitude/longitude dataset.'
        raise ValueError(msg)

    lats = np.array(lats, dtype=np.float32)
    lons = np.array(lons, dtype=np.float32)

    # y/xstep
    if ystep * xstep > 1:
        if lats.ndim == 1:
            lats = lats[::ystep]
            lons = lons[::xstep]
        elif lats.ndim == 2:
            lats = lats[::ystep, ::xstep]
            lons = lons[::ystep, ::xstep]

    return lats, lons


def get_lat_lon_rdc(meta):
    """Get 2D array of lat and lon for metadata dict in radar-coord.

    WARNING: This is a rough lat/lon value, NOT accurate!

    Parameters: meta : dict, including LENGTH, WIDTH and LAT/LON_REF1/2/3/4
    Returns:    lats : 2D np.array for latitude  in size of (length, width)
                lons : 2D np.array for longitude in size of (length, width)
    """
    if 'Y_FIRST' in meta.keys():
        raise Exception('Input file is in geo-coordinates, use more accurate get_lat_lon() instead.')

    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    lats = [float(meta[f'LAT_REF{i}']) for i in [1,2,3,4]]
    lons = [float(meta[f'LON_REF{i}']) for i in [1,2,3,4]]

    lat = np.zeros((length,width),dtype = np.float32)
    lon = np.zeros((length,width),dtype = np.float32)

    for i in range(length):
        for j in range(width):
            lat[i,j] = lats[0] + j*(lats[1] - lats[0])/width + i*(lats[2] - lats[0])/length
            lon[i,j] = lons[0] + j*(lons[1] - lons[0])/width + i*(lons[2] - lons[0])/length
    return lat, lon


def four_corners(atr):
    """Get the 4 corners coordinates from metadata dict in geo-coordinates.
    Parameters: atr - dict
    Returns:    south, north, west, east - float, in degrees or meters
    Examples:   S, N, W, E = ut.four_corners(atr)
                SNWE = ut.four_corners(atr)
    """
    width  = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    west  = float(atr['X_FIRST'])
    north = float(atr['Y_FIRST'])
    south = north + lat_step * length
    east  = west  + lon_step * width
    return south, north, west, east


def get_lalo_digit4display(meta, coord_unit='degree'):
    """Get the digit of the decimal place for the lat/lon info for display (e.g., at the status bar).

    Parameters: meta       - dict, metadata for the following attributes:
                             X_STEP
                             Y_STEP
                             RANGE_PIXEL_SIZE
                             AZIMUTH_PIXEL_SIZE
                coord_unit - str, coordinate unit, degree or meter
    Returns:    digit      - int, the digit for the decimal places of lat/lon
    """

    if coord_unit.startswith('meter'):
        digit = 2

    else:
        geo_step_keys = ['X_STEP', 'Y_STEP']
        rdr_step_keys = ['RANGE_PIXEL_SIZE', 'AZIMUTH_PIXEL_SIZE']

        # get step size in degree
        if all(x in meta.keys() for x in geo_step_keys):
            min_step = min([abs(float(meta[x])) for x in geo_step_keys])

        elif all(x in meta.keys() for x in rdr_step_keys):
            # default scaling for spaceborne system to ground range / azimuth
            rg_pix_size = float(meta['RANGE_PIXEL_SIZE']) / np.cos(np.deg2rad(30))
            az_pix_size = float(meta['AZIMUTH_PIXEL_SIZE']) * 0.9
            min_step = min([rg_pix_size, az_pix_size]) / 108e3
        else:
            # default pixel size 30 m
            min_step = 30 / 108e3

        # set the decimal place one order   (precision for step-range)
        if   min_step >= 2e-2:  digit = 3  # 110  m for >2.2 km
        elif min_step >= 2e-3:  digit = 4  # 11   m for >220 m
        elif min_step >= 2e-4:  digit = 5  # 1    m for >22  m
        elif min_step >= 2e-5:  digit = 6  # 0.1  m for >2.2 m
        else:                   digit = 7  # 0.01 m for <2.2 m

    return digit


###################################### Orbit ###########################################
def xyz_to_local_radius(xyz):
    """Calculate satellite height and ellipsoid local radius from orbital state vector.

    This is a simplified version of the following functions from ISCE-2:
    + isce.isceobj.Planet.xyz_to_llh()
    + isce.isceobj.Ellipsoid.localRadius()

    Parameters: xyz    - tuple of 3 float, orbital state vector
    Reference:  height - float, satellite altitude in m
                radius - float, Earth radius in m
    """

    # parameters from isce.isceobj.Planet.AstronomicalHandbook
    a = 6378137.000       # WGS84 semimajor
    e2 = 0.0066943799901  # WGS84 eccentricity squared

    # xyz --> llh
    a2 = a**2
    e4 = e2**2
    r_llh =  [0]*3
    d_llh = [[0]*3 for i in range(len(xyz))]
    xy2 = xyz[0]**2 + xyz[1]**2
    p = (xy2) / a2
    q = (1. - e2) * xyz[2]**2 / a2
    r = (p + q - e4) / 6.
    s = e4 * p * q / (4. * r**3)
    t = (1. + s + math.sqrt(s * (2. + s))) **(1. / 3.)
    u = r * (1. + t + 1./t)
    v = math.sqrt(u**2 + e4 * q)
    w = e2 * (u + v - q) / (2. * v)
    k = math.sqrt(u + v + w**2) - w
    D = k * math.sqrt(xy2) / (k + e2)
    r_llh[0] = math.atan2(xyz[2], D)
    r_llh[1] = math.atan2(xyz[1], xyz[0])
    r_llh[2] = (k + e2 - 1.) * math.sqrt(D**2 + xyz[2]**2) / k
    d_llh[0] = math.degrees(r_llh[0])
    d_llh[1] = math.degrees(r_llh[1])
    d_llh[2] = r_llh[2]
    height = r_llh[2]

    # calculate local radius
    b = a * (1.0 - e2)**0.5
    latg = math.atan( math.tan(math.radians(d_llh[0])) * a**2 / b**2 )
    arg = (math.cos(latg)**2 / a**2) + (math.sin(latg)**2 / b**2)
    radius = 1.0 / math.sqrt(arg)

    return height, radius



#################################### Geometry ##########################################
# Definition of angles:
# (los_)inc_angle - the incidence angle of the LOS vector (from the ground to the SAR platform)
#                   measured from vertical. Used in isce2.
# (los_)az_angle  - the azimuth   angle of the LOS vector (from the ground to the SAR platform)
#                   measured from the north, with anti-clockwise as positive. Used in isce2.
# orb_az_angle    - the azimuth   angle of the SAR platform's orbit (along-track direction)
#                   measured from the north, with anti-clockwise as positive
# head_angle      - the azimuth   angle of the SAR platform's orbit (along-track direction)
#                   measured from the north, with      clockwise as positive. Used in roipac.
#
# Typical values in deg for satellites with near-polar orbit:
#     AlosA: los_inc_angle = 34, los_az_angle =  102, orb_az_angle =  12, head_angle =  -12
#     AlosD: los_inc_angle = 34, los_az_angle = -102, orb_az_angle = 168, head_angle = -168
#     SenA : los_inc_angle = 40, los_az_angle =  102, orb_az_angle =  12, head_angle =  -12
#     SenD : los_inc_angle = 40, los_az_angle = -102, orb_az_angle = 168, head_angle = -168
#     NiA  : los_inc_angle = 42, los_az_angle =  -78, orb_az_angle =  12, head_angle =  -12
#     NiD  : los_inc_angle = 42, los_az_angle =   78, orb_az_angle = 168, head_angle = -168

def los2orbit_azimuth_angle(los_az_angle, look_direction='right'):
    """Convert the azimuth angle of the LOS vector to the one of the orbit flight vector.
    Parameters: los_az_angle - np.ndarray or float, azimuth angle of the LOS vector from the ground to the SAR platform
                               measured from the north with anti-clockwise direction as positive, in the unit of degrees
    Returns:    orb_az_angle - np.ndarray or float, azimuth angle of the SAR platform along track/orbit direction
                               measured from the north with anti-clockwise direction as positive, in the unit of degrees
    """
    if look_direction == 'right':
        orb_az_angle = los_az_angle - 90
    else:
        orb_az_angle = los_az_angle + 90
    orb_az_angle -= np.round(orb_az_angle / 360.) * 360.
    return orb_az_angle


def azimuth2heading_angle(az_angle, look_direction='right'):
    """Convert azimuth angle from ISCE los.rdr band2 into satellite orbit heading angle

    ISCE-2 los.* file band2 is azimuth angle of LOS vector from ground target to the satellite
        measured from the north in anti-clockwise as positive

    Below are typical values in deg for satellites with near-polar orbit:
        ascending  orbit: heading angle of -12  and azimuth angle of 102
        descending orbit: heading angle of -168 and azimuth angle of -102
    """
    if look_direction == 'right':
        head_angle = (az_angle - 90) * -1
    else:
        head_angle = (az_angle + 90) * -1
    head_angle -= np.round(head_angle / 360.) * 360.
    return head_angle


def heading2azimuth_angle(head_angle, look_direction='right'):
    """Convert satellite orbit heading angle into azimuth angle as defined in ISCE-2."""
    if look_direction == 'right':
        az_angle = (head_angle - 90) * -1
    else:
        az_angle = (head_angle + 90) * -1
    az_angle -= np.round(az_angle / 360.) * 360.
    return az_angle


def enu2los(v_e, v_n, v_u, inc_angle, head_angle=None, az_angle=None):
    """Project east/north/up motion into the line-of-sight (LOS) direction defined by incidence/azimuth angle.
    Parameters: v_e        - np.ndarray or float, displacement in east-west   direction, east  as positive
                v_n        - np.ndarray or float, displacement in north-south direction, north as positive
                v_u        - np.ndarray or float, displacement in vertical    direction, up    as positive
                inc_angle  - np.ndarray or float, incidence angle from vertical, in the unit of degrees
                head_angle - np.ndarray or float, azimuth angle of the SAR platform along track direction
                             measured from the north with clockwise direction as positive, in the unit of degrees
                az_angle   - np.ndarray or float, azimuth angle of the LOS vector from the ground to the SAR platform
                             measured from the north with anti-clockwise direction as positive, in the unit of degrees
                             head_angle = 90 - az_angle
    Returns:    v_los      - np.ndarray or float, displacement in LOS direction, motion toward satellite as positive
    """
    # unite (los_)head/az_angle into (los_)az_angle
    if az_angle is None:
        if head_angle is not None:
            az_angle = heading2azimuth_angle(head_angle)
        else:
            raise ValueError(f'invalid az_angle: {az_angle}!')

    # project ENU onto LOS
    v_los = (  v_e * np.sin(np.deg2rad(inc_angle)) * np.sin(np.deg2rad(az_angle)) * -1
             + v_n * np.sin(np.deg2rad(inc_angle)) * np.cos(np.deg2rad(az_angle))
             + v_u * np.cos(np.deg2rad(inc_angle)))

    return v_los


def en2az(v_e, v_n, orb_az_angle):
    """Project east/north motion into the radar azimuth direction.
    Parameters: v_e          - np.ndarray or float, displacement in east-west   direction, east  as positive
                v_n          - np.ndarray or float, displacement in north-south direction, north as positive
                orb_az_angle - np.ndarray or float, azimuth angle of the SAR platform along track/orbit direction
                               measured from the north with anti-clockwise direction as positive, in the unit of degrees
                               orb_az_angle = los_az_angle - 90 for right-looking radar.
    Returns:    v_az         - np.ndarray or float, displacement in azimuth direction,
                               motion along flight direction as positive
    """
    # project EN onto azimuth
    v_az = (  v_e * np.sin(np.deg2rad(orb_az_angle)) * -1
            + v_n * np.cos(np.deg2rad(orb_az_angle)))
    return v_az


def calc_azimuth_from_east_north_obs(east, north):
    """Calculate the azimuth angle of the given horizontal observation (in east and north)

    Parameters: east     - float,  eastward motion
                north    - float, northward motion
    Returns:    az_angle - float, azimuth angle in degree
                           measured from the north with anti-clockwise as positive
    """
    az_angle = -1 * np.rad2deg(np.arctan2(east, north)) % 360
    return az_angle


def get_unit_vector4component_of_interest(los_inc_angle, los_az_angle, comp='enu2los', horz_az_angle=None):
    """Get the unit vector for the component of interest.
    Parameters: los_inc_angle - np.ndarray or float, incidence angle from vertical, in the unit of degrees
                los_az_angle  - np.ndarray or float, azimuth angle of the LOS vector from the ground to the SAR platform
                                measured from the north with anti-clockwise direction as positive, in the unit of degrees
                comp          - str, component of interest, choose among the following values:
                                enu2los, en2los, hz2los, u2los, up2los, orb(it)_az, vert, horz
                horz_az_angle - np.ndarray or float, azimuth angle of the horizontal direction of interest
                                measured from the north with anti-clockwise direction as positive, in the unit of degrees
    Returns:    unit_vec      - list(np.ndarray/float), unit vector of the ENU component for the component of interest
    """
    # check input arguments
    comps = [
        'enu2los', 'en2los', 'hz2los', 'horz2los', 'u2los', 'vert2los',   # radar LOS / cross-track
        'en2az', 'hz2az', 'orb_az', 'orbit_az',                           # radar azimuth / along-track
        'vert', 'vertical', 'horz', 'horizontal',                         # vertical / arbitrary horizontal
    ]

    if comp not in comps:
        raise ValueError(f'un-recognized comp input: {comp}.\nchoose from: {comps}')

    if comp == 'horz' and horz_az_angle is None:
        raise ValueError('comp=horz requires horz_az_angle input!')

    # initiate output
    unit_vec = None

    if comp in ['enu2los']:
        unit_vec = [
            np.sin(np.deg2rad(los_inc_angle)) * np.sin(np.deg2rad(los_az_angle)) * -1,
            np.sin(np.deg2rad(los_inc_angle)) * np.cos(np.deg2rad(los_az_angle)),
            np.cos(np.deg2rad(los_inc_angle)),
        ]

    elif comp in ['en2los', 'hz2los', 'horz2los']:
        unit_vec = [
            np.sin(np.deg2rad(los_inc_angle)) * np.sin(np.deg2rad(los_az_angle)) * -1,
            np.sin(np.deg2rad(los_inc_angle)) * np.cos(np.deg2rad(los_az_angle)),
            np.zeros_like(los_inc_angle),
        ]

    elif comp in ['u2los', 'vert2los']:
        unit_vec = [
            np.zeros_like(los_inc_angle),
            np.zeros_like(los_inc_angle),
            np.cos(np.deg2rad(los_inc_angle)),
        ]

    elif comp in ['en2az', 'hz2az', 'orb_az', 'orbit_az']:
        orb_az_angle = los2orbit_azimuth_angle(los_az_angle)
        unit_vec = [
            np.sin(np.deg2rad(orb_az_angle)) * -1,
            np.cos(np.deg2rad(orb_az_angle)),
            np.zeros_like(orb_az_angle),
        ]

    elif comp in ['vert', 'vertical']:
        unit_vec = [0, 0, 1]

    elif comp in ['horz', 'horizontal']:
        unit_vec = [
            np.sin(np.deg2rad(horz_az_angle)) * -1,
            np.cos(np.deg2rad(horz_az_angle)),
            np.zeros_like(horz_az_angle),
        ]

    return unit_vec



#################################### Image Processing ##########################################
def wrap(data_in, wrap_range=[-1.*np.pi, np.pi]):
    """Wrap data into a range.
    Parameters: data_in    : np.array, array to be wrapped
                wrap_range : list of 2 float, range to be wrapped into
    Returns:    data       : np.array, data after wrapping
    """
    w0, w1 = wrap_range
    data = np.array(data_in)
    data = w0 + np.mod(data - w0, w1 - w0)
    return data


def diff_wrapped_phase(pha1, pha2):
    """Calculate the difference between two input wrapped phase.

    Parameters: pha1     - np.ndarray, (un)wrapped phase
                pha2     - np.ndarray / float, (un)wrapped phase
    Returns:    pha_diff - np.ndarray, wrapped phase difference for (pha1 - pha2)
    """
    pha_diff = np.angle(np.exp(-1j * pha1) * np.conj(np.exp(-1j * pha2)))
    return pha_diff


def get_all_conn_components(mask_in, min_num_pixel=1e4):
    """Get all connected component with number of pixels larger than threshold
    Parameters: mask_in  : 2D np.array with zero as background and non-zero as foreground
                min_num_pixel : int/float, minimum number of pixels to be identified and output
    Returns:    mask_out : list of 2D np.array in np.bool_ format
    """
    mask_in = np.array(mask_in)
    mask_out = []  # list of 2D np.array in bool
    mask_cc = get_largest_conn_component(mask_in, min_num_pixel=1e4)
    while not np.all(~mask_cc):
        mask_out.append(mask_cc)
        mask_in ^= mask_cc
        mask_cc = get_largest_conn_component(mask_in, min_num_pixel=1e4)
    return mask_out


def get_largest_conn_component(mask_in, min_num_pixel=1e4, display=False):
    """Extract the largest connected component from an 2D array
       with zero as background value
    Parameters: mask_in  : 2D np.array with zero as background and non-zero as foreground
                min_num_pixel : int/float, minimum number of pixels to be identified and output
                display : bool, display the result or not.
    Returns:    mask_out : 2D np.array in np.bool_ format
    """
    import matplotlib.pyplot as plt
    from scipy import ndimage

    mask_out = np.zeros(mask_in.shape, np.bool_)
    labels = ndimage.label(mask_in)[0]
    num_pixel = np.max(np.bincount(labels.flatten())[1:])
    if num_pixel < min_num_pixel:
        return mask_out

    max_label = np.argmax(np.bincount(labels.flatten())[1:]) + 1
    mask_out = labels == max_label
    if display:
        _, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 5])
        ax[0].imshow(mask_in)
        ax[1].imshow(mask_out)
        ax[2].imshow(mask_in ^ mask_out)
        plt.show()
    return mask_out


def min_region_distance(mask1, mask2, display=False):
    """Calculate the min distance between two regions of pixels marked by mask1 and mask2
    Parameters: mask1/2 : 2D np.array in size of (length, width) in np.bool_ format
    Returns:    pts1 : tuple of 2 int, bridge point in mask1, in (x, y)
                pts2 : tuple of 2 int, bridge point in mask2, in (x, y)
                min_dist : float, min euclidean distance
    """
    import matplotlib.pyplot as plt
    from scipy.spatial import cKDTree

    y, x = np.where(mask1 != 0)
    pts1 = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    tree = cKDTree(pts1)

    y, x = np.where(mask2 != 0)
    pts2 = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    dist, idx = tree.query(pts2)

    idx_min = np.argmin(dist)
    xy2 = pts2[idx_min]
    xy1 = pts1[idx[idx_min]]
    min_dist = dist[idx_min]

    if display:
        plt.figure()
        plt.imshow(mask1 * 1 + mask2 * 2)
        plt.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], '-o')
        plt.show()

    return xy1, xy2, min_dist


def interpolate_data(in_data, out_shape, interp_method='linear'):
    """Interpolate input 2D matrix into different shape.

    Used to get full resolution perp baseline from ISCE coarse grid baseline file.

    Parameters: in_data       : 2D np.ndarray
                out_shape     : tuple of 2 int in (length, width)
                interp_method : string, choose in [nearest, linear, cubic]
    Returns:    out_data      : 2D np.ndarray in out_shape
    """
    from scipy.interpolate import RegularGridInterpolator

    # prepare interpolation function
    in_shape = in_data.shape
    in_pts = (np.arange(in_shape[0]), np.arange(in_shape[1]))
    interp_func = RegularGridInterpolator(
        in_pts,
        in_data,
        method=interp_method,
        bounds_error=False,
    )

    # prepare output coordinates
    xx, yy = np.meshgrid(np.linspace(0, in_shape[1]-1, out_shape[1], endpoint=False),
                         np.linspace(0, in_shape[0]-1, out_shape[0], endpoint=False))
    out_pts = np.hstack((yy.reshape(-1, 1), xx.reshape(-1, 1)))

    # run interpolation
    out_data = interp_func(out_pts).reshape(out_shape)

    return out_data


def polygon2mask(polygon, shape):
    """Create a 2D mask (numpy array in binary) from a polygon.

    Link: https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
    Parameters: polygon - list of tuples of 2 int, e.g. [(x1, y1), (x2, y2), ...]
                shape   - list/tuple of 2 int, for length and width
    Returns:    mask    - 2D np.ndarray in bool in size of (length, width)
    """
    from PIL import Image, ImageDraw
    length, width = shape

    img = Image.new('L', (width, length), 0)
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
    mask = np.array(img, dtype=np.bool_)

    return mask


def get_circular_mask(x, y, radius, shape):
    """Get mask of pixels within circle defined by (x, y, r)"""
    length, width = shape
    yy, xx = np.ogrid[-y:length-y,
                      -x:width-x]
    cmask = (xx**2 + yy**2 <= radius**2)
    return cmask


def circle_index(atr, circle_par):
    """Return Index of Elements within a Circle centered at input pixel
    Parameters: atr : dictionary
                    containing the following attributes:
                    WIDT
                    LENGTH
                circle_par : string in the format of 'y,x,radius'
                    i.e. '200,300,20'          for radar coord
                         '31.0214,130.5699,20' for geo   coord
    Returns:    idx : 2D np.array in bool type
                    mask matrix for those pixel falling into the circle
                    defined by circle_par
    Examples:   idx_mat = ut.circle_index(atr, '200,300,20')
                idx_mat = ut.circle_index(atr, '31.0214,130.5699,20')
    """

    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])

    if isinstance(circle_par, tuple):
        cir_par = circle_par
    elif isinstance(circle_par, list):
        cir_par = circle_par
    else:
        cir_par = circle_par.replace(',', ' ').split()
    cir_par = [str(i) for i in cir_par]

    try:
        c_y = int(cir_par[0])
        c_x = int(cir_par[1])
        radius = int(float(cir_par[2]))
        print('Input circle index in y/x coord: %d, %d, %d' % (c_y, c_x, radius))
    except:
        try:
            c_lat = float(cir_par[0])
            c_lon = float(cir_par[1])
            radius = int(float(cir_par[2]))
            c_y = np.rint((c_lat-float(atr['Y_FIRST'])) / float(atr['Y_STEP']))
            c_x = np.rint((c_lon-float(atr['X_FIRST'])) / float(atr['X_STEP']))
            print('Input circle index in lat/lon coord: '
                  '{:.4f}, {:.4f}, {}'.format(c_lat, c_lon, radius))
        except:
            print('\nERROR: Unrecognized circle index format: '+circle_par)
            print('Supported format:')
            print('--circle 200,300,20            for radar coord input')
            print('--circle 31.0214,130.5699,20   for geo   coord input\n')
            return 0

    y, x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx



#################################### User Interaction #####################################
def yes_or_no(question):
    """garrettdreyfus on Github: https://gist.github.com/garrettdreyfus/8153571"""
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    elif reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhhh... please enter ")


def update_attribute_or_not(atr_new, atr_orig):
    """Compare new attributes with existing ones"""
    for key in atr_new.keys():
        value = str(atr_new[key])
        if not ((key in atr_orig.keys() and value == str(atr_orig[key]) and value != 'None')
                or (key not in atr_orig.keys() and value == 'None')):
            return True
    return False


def which(program):
    """Test if executable exists"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if os.path.split(program)[0]:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def check_parallel(file_num=1, print_msg=True, maxParallelNum=8):
    """Check parallel option based file num and installed module
    Link: https://joblib.readthedocs.io/en/latest/parallel.html
    Examples:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(file_list))
        Parallel(n_jobs=num_cores)(delayed(subset_file)(file, vars(inps)) for file in file_list)
    """
    enable_parallel = True

    # Disable parallel option for one input file
    if file_num <= 1:
        enable_parallel = False
        if print_msg:
            print('parallel processing is disabled for one input file')
        return 1, enable_parallel, None, None

    # Check required python module
    try:
        from joblib import Parallel, delayed
    except:
        print('Can not import joblib')
        print('parallel is disabled.')
        enable_parallel = False
        return 1, enable_parallel, None, None

    # Find proper number of cores for parallel processing
    num_cores = min(os.cpu_count(), file_num, maxParallelNum)
    if num_cores <= 1:
        enable_parallel = False
        print('parallel processing is disabled because min of the following two numbers <= 1:')
        print(f'available cpu number of the computer: {os.cpu_count()}')
    elif print_msg:
        print('parallel processing using %d cores ...' % (num_cores))

    try:
        return num_cores, enable_parallel, Parallel, delayed
    except:
        return num_cores, enable_parallel, None, None


def print_command_line(script_name, args):
    """print the command line with "" for arguments containing *.

    Parameters: script_name - str, e.g. prep_isce.py
                args        - list(str), list of input arguments
    """
    cmd = script_name
    for arg in args:
        # for option values containing *, add "" in the print out msg
        # so that it can be copied and pasted to run directly.
        if not arg.startswith('-') and '*' in arg:
            cmd += f' "{arg}"'
        else:
            cmd += f' {arg}'

    print(cmd)
    return


#################################### Math / Statistics ###################################
def median_abs_deviation(data, center=None, scale=0.67449):
    """Compute the median absolute deviation (MAD) of the data along the LAST axis.

    This function is also available in the following packages:
    + statsmodels.robust.mad()
      https://www.statsmodels.org/dev/generated/statsmodels.robust.scale.mad.html
    + scipy.stats.median_abs_deviation() since v1.5.0
      https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html

    The implementation here is preferred because we would like to:
    1. omit the NaN value in the data for the median and MAD calculation
    2. scale the returned value to be comparable with standard deviation (STD)
       for easy interpretation with 1/2/3-sigma rule.
    While statsmodels could not handle the NaN value in the data; and scipy does not hebave
    as desired in the default setting (propagate nan and do not scale).

    The following two are equivalent for a 2D np.ndarray X:
    + mintpy.utils.utils0.median_abs_deviation(X)
    + scipy.stats.median_abs_deviation(X, axis=-1, center=np.nanmedian, scale=0.67449, nan_policy='omit')

    Parameters: data   - 1/2D np.ndarray, input array
                center - 0/1D np.ndarray or None
                scale  - float, the normalization constant
    Returns:    mad    - 0/1D np.ndarray
    """
    # check input: compute MAD over 1/2D matrix only
    data = np.array(data)
    if data.ndim > 2:
        ntime = data.shape[0]
        data = data.reshape(ntime, -1)

    # calculate the center value: median along the last axis
    if center is None:
        center = np.nanmedian(data, axis=-1)

        # replicate center matrix for 2D matrix calculation later on
        if data.ndim == 2:
            center = np.tile(center.reshape(-1, 1), (1, data.shape[1]))

    # calculate median absolute deviation
    mad = np.nanmedian(np.abs(data - center), axis=-1) / scale

    return mad


def median_abs_deviation_threshold(data, center=None, cutoff=3.):
    """calculate rms_threshold based on the standardized residual

    Outlier detection with median absolute deviation.
    """
    if center is None:
        center = np.nanmedian(data)
    mad = median_abs_deviation(data, center=center)
    threshold = center + cutoff * mad
    return threshold


def root_mean_sq_error(x, y=None):
    """Calculate the root-mean-square error between x and y."""
    # make a copy & ensure 1D numpy array format
    x = np.array(x).flatten()
    if y is not None:
        y = np.array(y).flatten()
        if x.size != y.size:
            raise ValueError(f'Input x & y have different size: {x.size} vs. {y.size}!')
        x -= y

    # omit NaN values
    x = x[~np.isnan(x)]
    rmse = np.sqrt(np.sum(x**2) / (x.size - 1))
    return rmse


def ceil_to_1(x):
    """Return the most significant digit of input number and ceiling it"""
    digit = int(np.floor(np.log10(abs(x))))
    x_round = round(x, -digit)
    # round to ceil
    if x_round >= x:
        x_ceil = x_round
    else:
        x_ceil = x_round + 10**digit
    return x_ceil


def round_to_1(x):
    """Return the most significant digit of input number"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)


def round_up_to_odd(x):
    """Round a float up to the next odd integer .
    Link: https://stackoverflow.com/questions/31648729
    """
    y = np.ceil(x) // 2 * 2 + 1
    return y.astype(np.int16)


def highest_power_of_2(x):
    """Given a number x, find the highest power of 2 that <= x"""
    res = np.power(2, np.floor(np.log2(x)))
    res = np.int16(res)
    return res


def most_common(L, k=1):
    """Return the k most common item in the list L.
    Examples:
        5, 8 = most_common([4,5,5,5,5,8,8,8,9], k=2)
        'duck' = most_common(['goose','duck','duck','dog'])
        'goose' = most_common(['goose','duck','duck','goose'])
    """
    from collections import Counter
    cnt = Counter(L)
    item_mm = [i[0] for i in cnt.most_common(k)]
    if k == 1:
        item_mm = item_mm[0]
    return item_mm


def is_number(string):
    """Check string is a number.
    Not using str.isnumeric() because it can not handle floating point nor negative sign.
    Link: https://elearning.wsldp.com/python3/python-check-string-is-a-number/
    Parameters: string - str, a string
    Returns:    True/False
    """
    try:
        float(string)
        return True
    except ValueError:
        return False
