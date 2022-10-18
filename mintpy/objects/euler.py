############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Oct 2022                           #
############################################################
# Utility scripts for Euler pole handling
# Recommend import:
#     from mintpy.objects.euler import EulerPole


import sys
import pyproj
import numpy as np
from mintpy.objects.constants import EARTH_RADIUS

################################################################################################
## global variables
MAS2RAD      = np.pi/3600000/180    #  1 mas (milliarcsecond) = yy radian
MASY2DMY     = 1e6 / 3600000        #  1 mas per year         = xx degree per million year


## Euler pole class object

class EulerPole(object):

    def __init__(self, pole, coord='cartesian', unit='mas/yr'):
        """Initialize the Euler pole
        The unit of angular velocity must be 'mas per year'.
        Parameters: pole - list or 1D array
                    (positive: rotational vector pointing outward of a sphere;
                    counterclockwise rotation when looking from outside of a sphere)
                        (1) [wx, wy, wz] when coord='cartesian'
                            wx        float, angular velocity in x-axis      unit
                            wy        float, angular velocity in y-axis      unit
                            wz        float, angular velocity in z-axis      unit
                        (2) [plat, plon, omega] when coord='spherical'
                            plat      float, Euler pole latitude             [degree]
                            plon      float, Euler pole longitude            [degree]
                            rate      float, magnitude of angular velocity   unit
                    coord -  str,   ['cartesian', 'spherical']
                    unit  -  str,   ['mas/yr', 'deg/Ma']
        """
        ## Input the Euler pole
        if coord.lower() == 'cartesian':
            wx, wy, wz = np.array(pole, dtype=float)
            if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
                wx /= MASY2DMY
                wy /= MASY2DMY
                wz /= MASY2DMY
            plat, plon, rate = cart2sph(wx, wy, wz)

        elif coord.lower() == 'spherical':
            plat, plon, rate = np.array(pole).astype('float')
            if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
                rate /= MASY2DMY
            wx, wy, wz = sph2cart(plat, plon, r=rate)

        else:
            print('Some user input is wrong here...')
            sys.exit(1)

        self.coord = coord.lower()
        self.plat  = plat
        self.plon  = plon
        self.rate  = rate
        self.wx    = wx
        self.wy    = wy
        self.wz    = wz


    def __sub__(self, other):
        pole = np.array(self.get_omega()) - np.array(other.get_omega())
        return EulerPole(pole)


    def __add__(self, other):
        pole = np.array(self.get_omega()) + np.array(other.get_omega())
        return EulerPole(pole)


    def __neg__(self):
        pole = -np.array(self.get_omega())
        return EulerPole(pole)


    def print_msg(self):
        print('\n------------------ Euler Pole description ------------------')
        print('Spherical expression:')
        print(f'   Pole Latitude  : {self.plat:9.4f} DEG')
        print(f'   Pole Longitude : {self.plon:9.4f} DEG')
        print(f'   Rotation rate  : {self.rate*MASY2DMY:9.4f} DEG/MA = {self.rate:9.4f} MAS/Y')
        print('Cartesian expression (angular velocity vector):')
        print(f'   wx             : {self.wx*MASY2DMY:9.4f} DEG/MA = {self.wx:9.4f} MAS/Y')
        print(f'   wy             : {self.wy*MASY2DMY:9.4f} DEG/MA = {self.wy:9.4f} MAS/Y')
        print(f'   wz             : {self.wz*MASY2DMY:9.4f} DEG/MA = {self.wz:9.4f} MAS/Y')
        print('------------------------------------------------------------\n')


    def get_omega(self):
        """ Return the angular velocity vector
        Returns:    wx in MAS/Y
                    wy in MAS/Y
                    wz in MAS/Y
        """
        return self.wx, self.wy, self.wz


    def get_pole(self):
        """ Return the Euler pole
        Returns:    plat in DEG
                    plon in DEG
                    rate in MAS/Y
        """
        return self.plat, self.plon, self.rate


    def velocity(self, lat, lon):
        """
        Calculates the azimuth (in degrees) and rate of plate motion
        (in millimeters per year) at a given point.
        Parameters:  lat     - Latitude of the point                [DEG]
                     lon     - Longitude of the point               [DEG]
        Returns:     azi     - Azimuth angle clockwise from north   [DEG]
                     v_rate  - Velocity                             [mm/yr]
        """
        ve, vn, vu = self.velocity_enu(lat, lon)
        azi        = azimuth(ve, vn)
        v_rate     = np.sqrt(ve**2 + vn**2 + vu**2)
        return azi, v_rate


    def velocity_enu(self, lat, lon, alt=0.0, ellps=True):
        ## Convert to local ENU linear velocity
        """Given Euler pole object, compute V_xyz for given pixel(s) at (lat,lon) of interest.
        Only supports a constant altitude now.
        Parameters: lat   -      points of interest (latitude)      [degree]
                    lon   -      points of interest (longitude)     [degree]
                    alt   -      points of interest (altitude)      [meter]
                    ellps -      True/False; consider ellipsoidal Earth projection of point(s) of interest
        Returns:    ve    -      east  linear velocity              [meter/year]
                    vn    -      north linear velocity              [meter/year]
                    vu    -      up    linear velocity              [meter/year]

        V = T * Vxyz, where V is the local ENU velocity
        """
        vx, vy, vz = self.velocity_xyz(lat, lon, alt=alt, ellps=ellps)
        ve, vn, vu = rotate_coord_cart2enu(lat, lon, vx, vy, vz)
        return ve, vn, vu


    def velocity_xyz(self, lat, lon, alt=0.0, ellps=True):
        """Given Euler pole object, compute V_xyz for given pixel(s) at (lat,lon) of interest.
        Only supports a constant altitude now.
        Parameters: lat   -      points of interest (latitude)      [degree]
                    lon   -      points of interest (longitude)     [degree]
                    alt   -      points of interest (altitude)      [meter]
                    ellps -      True/False; consider ellipsoidal Earth projection of point(s) of interest
        Returns:    vx    -      ECEF x linear velocity             [meter/year]
                    vy    -      ECEF y linear velocity             [meter/year]
                    vz    -      ECEF z linear velocity             [meter/year]
        """
        ## report how many points of interest
        if isinstance(lat, (list, tuple, np.ndarray)):
            npts = len(lat)
        elif isinstance(lat, (int, float)):
            npts = 1
        print(f'number of points to compute: {npts}')

        ## Local coordinates handling for location(s) of interest
        if not ellps:
            # a perfect sphere
            print(f'Assume a perfect spherical Earth with radius={EARTH_RADIUS} m')
            x, y, z = sph2cart(lat, lon, EARTH_RADIUS)   # meter
        else:
            # WGS84 ellips; only supports uniform altitude now, but can change later
            print('Assume WGS84 ellipse from pyproj')
            x, y, z = coord_lalo2xyz(lat, lon, alt)     # meter
        xyz = np.array([x, y, z], dtype=float)        # meter
        shp = xyz.shape

        ## Compute the cartesian linear velocity (i.e., ECEF) in km/year
        # Vxyz = omega x R_i, where R_i is location vector at pixel i
        omega = np.array([self.wx, self.wy, self.wz]) * MAS2RAD
        V     = np.cross(omega, xyz.T).T.reshape(shp)
        vx, vy, vz = V[0], V[1], V[2]          # meter/year
        return vx, vy, vz



################################################################################################
## Sub-functions for the math of Euler Pole and linear velocity

def sph2cart(lat, lon, r=1):
    """Convert spherical coordinates to cartesian.

    Parameters: lat    - float / np.ndarray, latitude  [degree]
                lon    - float / np.ndarray, longitude [degree]
                r      - float / np.ndarray, radius [any units of angular distance]
    Returns:    rx/y/z - float / np.ndarray, angular distance in X/Y/Z direction [same unit as r]
    Examples:
        # spherical coord to xyz coord
        x, y, z = sph2cart(lat, lon, r=radius)
        # spherical Euler pole to cartesian rotational vector
        wx, wy, wz = sph2cart(plat, plon, r=rate)
    """
    rx = r * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    ry = r * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    rz = r * np.sin(np.deg2rad(lat))
    return rx, ry, rz


def cart2sph(rx, ry, rz):
    """Convert cartesian coordinates to spherical.

    Parameters: rx/y/z - float / np.ndarray, angular distance in X/Y/Z direction [any units of distance]
    Returns:    lat    - float / np.ndarray, latitude  [degree]
                lon    - float / np.ndarray, longitude [degree]
                r      - float / np.ndarray, radius [same unit as rx/y/z]
    Examples:
        # xyz coord to spherical coord
        lat, lon, r = cart2sph(x, y, z)
        # cartesian rotational vector to spherical Euler pole
        plat, plon, rate = cart2sph(wx, wy, wz)
    """
    r = np.sqrt(rx**2 + ry**2 + rz**2)
    lat = np.rad2deg(np.arcsin(rz / r))
    lon = np.rad2deg(np.arctan2(ry, rx))
    return lat, lon, r


def rotate_coord_cart2enu(lat, lon, x, y, z, reverse=False):
    """Rotation matrix to convert cartesian coordinates (global ECEF) to ENU components (local cartesian)
    at given (lat, lon)

    Reference:
        Navipedia, https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications,
          Palo Alto, doi: 10.4236/ojapps.2015.54016. Page 145-156

    Parameters: lat     - 1D np.ndarray, latitude [degree]
                lon     - 1D np.ndarray, longitude [degree]
                reverse - bool, revert to the enu2cart conversion
    Returns:    e       - 1D np.ndarray, east  component
                n       - 1D np.ndarray, north component
                u       - 1D np.ndarray, up    component
    """
    lat, lon = check_lalo(lat, lon)
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # construct rotation matrix
    if not reverse:
        # cart2enu
        e = - np.sin(lon) * x + np.cos(lon) * y
        n = - np.sin(lat) * np.cos(lon) * x \
            - np.sin(lat) * np.sin(lon) * y \
            + np.cos(lat) * z
        u =   np.cos(lat) * np.cos(lon) * x \
            + np.cos(lat) * np.sin(lon) * y \
            + np.sin(lat) * z

    else:
        # enu2cart
        e = - np.sin(lon) * x \
            - np.cos(lon) * np.sin(lat) * y \
            + np.cos(lon) * np.cos(lat) * z
        n =   np.cos(lon) * x \
            - np.sin(lon) * np.sin(lat) * y \
            + np.sin(lon) * np.cos(lat) * z
        u =   np.cos(lat) * y \
            + np.sin(lat) * z

    return e, n, u


def coord_lalo2xyz(lat, lon, alt):
    """Convert coordinates from WGS84 lat/long to ECEF x/y/z.

    Parameters: lat   - float / np.ndarray, latitude  [degree]
                lon   - float / np.ndarray, longitude [degree]
                alt   - float / np.ndarray, altitude  [meter]
    Returns:    x/y/z - float / np.ndarray, ECEF coordinate [meter]
    """
    lat, lon = check_lalo(lat, lon)
    if isinstance(lat, np.ndarray) and not isinstance(alt, np.ndarray):
        alt *= np.ones_like(lat)

    # construct pyproj transform object
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    )

    # apply transformation
    x, y, z = transformer.transform(lon, lat, alt, radians=False)

    return x, y, z


def azimuth(east, north):
    """Azimuth in degrees clockwise from North
    Parameter:  east  - float, eastward motion
                north - float, northward motion
    """
    azi = np.rad2deg(np.arctan2(east, north)) % 360
    return azi


def check_lalo(lat, lon):
    """ ensure numpy array format
    """
    if isinstance(lat, (list, tuple, np.ndarray)):
        lat = np.array(lat).flatten()
        lon = np.array(lon).flatten()

        if not lat.size == lon.size:
            raise ValueError(f'Inconsistent size between input lat ({lat.size}) and lon ({lon.size})!')
    return lat, lon
