############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Low level utilities script (independent)
# Contents
#   InSAR
#   File Operation
#   Geometry
#   Image Processing
#   User Interaction
#   Math / Statistics
# Recommend import:
#   from mintpy.utils import utils as ut


import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from scipy import ndimage


# global variables
SPEED_OF_LIGHT = 299792458 # m/s
EARTH_RADIUS = 6371e3      # Earth radius in meters
K = 40.31                  # m^3/s^2, constant


#################################### InSAR ##########################################
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
    r = float(atr['EARTH_RADIUS'])
    H = float(atr['HEIGHT'])
    width = int(atr['WIDTH'])

    # Calculation
    range_f = range_n+dR*width
    inc_angle_n = (np.pi - np.arccos((r**2 + range_n**2 - (r+H)**2)/(2*r*range_n))) * 180.0/np.pi
    inc_angle_f = (np.pi - np.arccos((r**2 + range_f**2 - (r+H)**2)/(2*r*range_f))) * 180.0/np.pi
    inc_angle_c = (inc_angle_n + inc_angle_f) / 2.0
    if print_msg:
        print('near   incidence angle : {:.4f} degree'.format(inc_angle_n))
        print('center incidence angle : {:.4f} degree'.format(inc_angle_c))
        print('far    incidence angle : {:.4f} degree'.format(inc_angle_f))

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
        raise Exception('un-supported dimension input: {}'.format(dimension))
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
    try:
        proc = atr['PROCESSOR']
    except:
        proc = 'isce'
    if proc in ['roipac', 'isce']:
        Re = float(atr['EARTH_RADIUS'])
        Height = float(atr['HEIGHT'])
        az_step = float(atr['AZIMUTH_PIXEL_SIZE']) * Re/(Re+Height)
    elif proc == 'gamma':
        az_step = float(atr['AZIMUTH_PIXEL_SIZE'])
    return az_step



#################################### File Operation ##########################################
def touch(fname_list, times=None):
    """python equivalent function to Unix utily - touch
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



#################################### Geometry ##########################################

def lalo_ground2iono_shell_along_los(lat, lon, inc_angle=30, head_angle=-168, iono_height=450e3):
    """Convert the lat/lon of a point on the ground to the ionosphere thin-shell 
    along the line-of-sight (LOS) direction.

    Reference: Jingyi, C., and H. A. Zebker (2012), Ionospheric Artifacts in Simultaneous
        L-Band InSAR and GPS Observations, Geoscience and Remote Sensing, IEEE Transactions on,
        50(4), 1227-1239, doi:10.1109/TGRS.2011.2164805.

    Parameters: lat/lon     - float, latitude/longitude of the point on the ground in degrees
                inc_angle   - float, incidence angle of the line-of-sight on the ground in degrees
                head_angle  - float, heading angle of the satellite orbit in degrees
                              from the north direction with positive in clockwise direction
                iono_height - float, height of the ionosphere thin-shell in meters
    """
    # degrees to radians
    inc_angle /= 180 / np.pi
    head_angle /= 180 / np.pi

    # offset angle from equation (25) in Chen and Zebker (2012)
    off_iono = inc_angle - np.arcsin(EARTH_RADIUS / (EARTH_RADIUS + iono_height) * np.sin(np.pi - inc_angle))

    # update lat/lon
    lat += off_iono * np.cos(head_angle) * 180 / np.pi
    lon += off_iono * np.sin(head_angle) * 180 / np.pi
    return lat, lon


def incidence_angle_ground2iono_shell_along_los(inc_angle, iono_height=450e3):
    """Calibrate the incidence angle of LOS vector on the ground surface to the ionosphere shell
    based on equation (6) in Chen ang Zebker (2012, TGRS)

    Reference: Jingyi, C., and H. A. Zebker (2012), Ionospheric Artifacts in Simultaneous
        L-Band InSAR and GPS Observations, Geoscience and Remote Sensing, IEEE Transactions on,
        50(4), 1227-1239, doi:10.1109/TGRS.2011.2164805.

    Parameters: inc_angle      - float/np.ndarray, incidence angle on the ground in degrees
                iono_height    - float, effective ionosphere height in meters
                                 under the thin-shell assumption
    Returns:    inc_angle_iono - float/np.ndarray, incidence angle on the iono shell in degrees
    """
    # ignore nodata in inc_angle
    if type(inc_angle) is np.ndarray:
        inc_angle[inc_angle == 0] = np.nan
    inc_angle *= np.pi / 180.

    # calculation
    cos_inc_angle_iono = np.sqrt(1 - (EARTH_RADIUS * np.sin(inc_angle) / (EARTH_RADIUS + iono_height))**2)
    inc_angle_iono = np.arccos(cos_inc_angle_iono) / np.pi * 180.0
    return inc_angle_iono


def get_lat_lon(meta, geom_file=None, box=None):
    """Extract precise pixel-wise lat/lon.

    For meta dict in geo-coordinates OR geom_file with latitude/longitude dataset

    Returned lat/lon are corresponds to the pixel center

    Parameters: meta : dict, including LENGTH, WIDTH and Y/X_FIRST/STEP
                box  : 4-tuple of int for (x0, y0, x1, y1)
    Returns:    lats : 2D np.array for latitude  in size of (length, width)
                lons : 2D np.array for longitude in size of (length, width)
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
        # generate 2D matrices for lat/lon
        lat_step = float(meta['Y_STEP'])
        lon_step = float(meta['X_STEP'])
        lat0 = float(meta['Y_FIRST']) + lat_step * (box[1] + 0.5)
        lon0 = float(meta['X_FIRST']) + lon_step * (box[0] + 0.5)

        lat_num = box[3] - box[1]
        lon_num = box[2] - box[0]

        lat1 = lat0 + lat_step * lat_num
        lon1 = lon0 + lon_step * lon_num
        lats, lons = np.mgrid[lat0:lat1:lat_num*1j,
                              lon0:lon1:lon_num*1j]

    else:
        msg = 'Can not get pixel-wise lat/lon!'
        msg += '\nmeta dict is not geocoded and/or geometry file does not contains latitude/longitude dataset.'
        raise ValueError(msg)

    lats = np.array(lats, dtype=np.float32)
    lons = np.array(lons, dtype=np.float32)
    return lats, lons


def get_lat_lon_rdc(meta):
    """Get 2D array of lat and lon.
    For metadata dict in radar-coord
    Parameters: meta : dict, including LENGTH, WIDTH and LAT/LON_REF1/2/3/4
    Returns:    lats : 2D np.array for latitude  in size of (length, width)
                lons : 2D np.array for longitude in size of (length, width)
    """
    if 'Y_FIRST' in meta.keys():
        raise Exception('Input file is in geo-coordinates, use more accurate get_lat_lon() instead.')

    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    lats = [float(meta['LAT_REF{}'.format(i)]) for i in [1,2,3,4]]
    lons = [float(meta['LON_REF{}'.format(i)]) for i in [1,2,3,4]]
    
    lat = np.zeros((length,width),dtype = np.float32)
    lon = np.zeros((length,width),dtype = np.float32)
    
    for i in range(length):
        for j in range(width):
            lat[i,j] = lats[0] + j*(lats[1] - lats[0])/width + i*(lats[2] - lats[0])/length
            lon[i,j] = lons[0] + j*(lons[1] - lons[0])/width + i*(lons[2] - lons[0])/length
    return lat, lon


def azimuth2heading_angle(az_angle):
    """Convert azimuth angle from ISCE los.rdr band2 into satellite orbit heading angle

    ISCE-2 los.* file band2 is azimuth angle of LOS vector from ground target to the satellite 
        measured from the north in anti-clockwise as positive

    Below are typical values in deg for satellites with near-polar orbit:
        ascending  orbit: heading angle of -12  and azimuth angle of 102
        descending orbit: heading angle of -168 and azimuth angle of -102
    """

    head_angle = -1 * (180 + az_angle + 90)
    head_angle -= np.round(head_angle / 360.) * 360.
    return head_angle


def enu2los(e, n, u, inc_angle=34., head_angle=-168.):
    """
    Parameters: e          - np.array or float, displacement in east-west direction, east as positive
                n          - np.array or float, displacement in north-south direction, north as positive
                u          - np.array or float, displacement in vertical direction, up as positive
                inc_angle  - np.array or float, local incidence angle from vertical
                head_angle - np.array or float, satellite orbit from the north in clock-wise direction as positive
    Returns:    v_los      - np.array or float, displacement in line-of-sight direction, moving toward satellite as positive

    Typical values in deg for satellites with near-polar orbit:
        For AlosA: inc_angle = 34, head_angle = -12.9,  az_angle = 102.9
        For AlosD: inc_angle = 34, head_angle = -167.2, az_angle = -12.8
        For  SenD: inc_angle = 34, head_angle = -168.0, az_angle = -102
    """
    # if input angle is azimuth angle
    if np.abs(np.abs(head_angle) - 90) < 30:
        head_angle = azimuth2heading_angle(head_angle)

    inc_angle *= np.pi/180.
    head_angle *= np.pi/180.
    v_los = (  e * np.sin(inc_angle) * np.cos(head_angle) * -1
             + n * np.sin(inc_angle) * np.sin(head_angle) 
             + u * np.cos(inc_angle))
    return v_los

def four_corners(atr):
    """Return 4 corners lat/lon"""
    width  = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    west  = float(atr['X_FIRST'])
    north = float(atr['Y_FIRST'])
    south = north + lat_step * length
    east  = west  + lon_step * width
    return west, east, south, north


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
                    containging the following attributes:
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
            print(('Input circle index in lat/lon coord: '
                   '{:.4f}, {:.4f}, {}'.format(c_lat, c_lon, radius)))
        except:
            print('\nERROR: Unrecognized circle index format: '+circle_par)
            print('Supported format:')
            print('--circle 200,300,20            for radar coord input')
            print('--circle 31.0214,130.5699,20   for geo   coord input\n')
            return 0

    y, x = np.ogrid[-c_y:length-c_y, -c_x:width-c_x]
    idx = x**2 + y**2 <= radius**2

    return idx


def subset_attribute(atr_dict, subset_box, print_msg=True):
    """Update attributes dictionary due to subset
    Parameters: atr_dict : dict, data attributes to update
                subset_box : 4-tuple of int, subset box defined in (x0, y0, x1, y1)
    Returns:    atr : dict, updated data attributes
    """
    if subset_box is None:
        return atr_dict

    sub_x = [subset_box[0], subset_box[2]]
    sub_y = [subset_box[1], subset_box[3]]

    # Update attribute variable
    atr = dict(atr_dict)
    atr['LENGTH'] = str(sub_y[1]-sub_y[0])
    atr['WIDTH'] = str(sub_x[1]-sub_x[0])
    atr['YMAX'] = str(sub_y[1]-sub_y[0] - 1)
    atr['XMAX'] = str(sub_x[1]-sub_x[0] - 1)
    if print_msg:
        print('update LENGTH, WIDTH, Y/XMAX')

    # Subset atribute
    atr['SUBSET_YMAX'] = str(sub_y[1] + int(atr_dict.get('SUBSET_YMIN', '0')))
    atr['SUBSET_YMIN'] = str(sub_y[0] + int(atr_dict.get('SUBSET_YMIN', '0')))
    atr['SUBSET_XMAX'] = str(sub_x[1] + int(atr_dict.get('SUBSET_XMIN', '0')))
    atr['SUBSET_XMIN'] = str(sub_x[0] + int(atr_dict.get('SUBSET_XMIN', '0')))
    if print_msg:
        print(('update/add SUBSET_XMIN/YMIN/XMAX/YMAX: '
               '{x0}/{y0}/{x1}/{y1}').format(x0=atr['SUBSET_XMIN'],
                                             y0=atr['SUBSET_YMIN'],
                                             x1=atr['SUBSET_XMAX'],
                                             y1=atr['SUBSET_YMAX']))

    # Geo coord
    if 'Y_FIRST' in atr.keys():
        atr['Y_FIRST'] = str(float(atr['Y_FIRST']) + sub_y[0]*float(atr['Y_STEP']))
        atr['X_FIRST'] = str(float(atr['X_FIRST']) + sub_x[0]*float(atr['X_STEP']))
        if print_msg:
            print('update Y/X_FIRST')

    # Reference in space
    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str(int(atr['REF_Y']) - sub_y[0])
        atr['REF_X'] = str(int(atr['REF_X']) - sub_x[0])
        if print_msg:
            print('update REF_Y/X')

    # Starting Range for file in radar coord
    if 'Y_FIRST' not in atr_dict.keys():
        try:
            atr['STARTING_RANGE'] = float(atr['STARTING_RANGE'])
            atr['STARTING_RANGE'] += float(atr['RANGE_PIXEL_SIZE'])*sub_x[0]
            if print_msg:
                print('update STARTING_RANGE')
        except:
            pass

    return atr


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
    mask_out = np.zeros(mask_in.shape, np.bool_)
    labels, n_features = ndimage.label(mask_in)
    num_pixel = np.max(np.bincount(labels.flatten())[1:])
    if num_pixel < min_num_pixel:
        return mask_out

    max_label = np.argmax(np.bincount(labels.flatten())[1:]) + 1
    mask_out = labels == max_label
    if display:
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 5])
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


def interpolate_data(inData, outShape, interpMethod='linear'):
    """Interpolate input 2D matrix into different shape.
    Used to get full resolution perp baseline from ISCE coarse grid baseline file.
    Parameters: inData : 2D array
                outShape : tuple of 2 int in (length, width)
                interpMethod : string, choose in [nearest, linear, cubic]
    Returns:    outData : 2D array in outShape
    """
    from scipy.interpolate import RegularGridInterpolator as RGI
    inShape = inData.shape
    inPts = (np.arange(inShape[0]), np.arange(inShape[1]))
    xx, yy = np.meshgrid(np.linspace(0, inShape[1]-1, outShape[1], endpoint=False),
                         np.linspace(0, inShape[0]-1, outShape[0], endpoint=False))
    outPts = np.hstack((yy.reshape(-1, 1), xx.reshape(-1, 1)))
    outData = RGI(inPts, inData, method=interpMethod,
                  bounds_error=False)(outPts).reshape(outShape)
    return outData




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
    """Compare new attributes with exsiting ones"""
    update = False
    for key in atr_new.keys():
        value = str(atr_new[key])
        if ((key in atr_orig.keys() and value == str(atr_orig[key]) and value != 'None')
                or (key not in atr_orig.keys() and value == 'None')):
            next
        else:
            update = True
    return update


def which(program):
    """Test if executable exists"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
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
    Examples:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(file_list))
        Parallel(n_jobs=num_cores)(delayed(subset_file)(file, vars(inps)) for file in file_list)
    """
    enable_parallel = True

    # Disable parallel option for one input file
    if file_num <= 1:
        enable_parallel = False
        if print_msg:
            print('parallel processing is diabled for one input file')
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
    num_cores = min(multiprocessing.cpu_count(), file_num, maxParallelNum)
    if num_cores <= 1:
        enable_parallel = False
        print('parallel processing is disabled because min of the following two numbers <= 1:')
        print('available cpu number of the computer: {}'.format(multiprocessing.cpu_count()))
    elif print_msg:
        print('parallel processing using %d cores ...' % (num_cores))

    try:
        return num_cores, enable_parallel, Parallel, delayed
    except:
        return num_cores, enable_parallel, None, None



#################################### Math / Statistics ###################################
def median_abs_deviation_threshold(data, center=None, cutoff=3.):
    """calculate rms_threshold based on the standardised residual
    outlier detection with median absolute deviation.
    https://www.statsmodels.org/dev/generated/statsmodels.robust.scale.mad.html
    """
    data = np.array(data)
    if center is None:
        center = np.median(data)
    #from statsmodels.robust import mad
    #data_mad = mad(data, center=center)
    data_mad = np.median(np.abs(data - center)) / 0.67448975019608171
    threshold = center + cutoff * data_mad
    return threshold


def ceil_to_1(x):
    """Return the most significant digit of input number and ceiling it"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)+10**digit


def round_to_1(x):
    """Return the most significant digit of input number"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -1*digit)


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



