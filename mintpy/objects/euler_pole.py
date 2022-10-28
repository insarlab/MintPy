############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Oct 2022                           #
############################################################
# Utility scripts for Euler pole handling
# Recommend import:
#     from mintpy.objects.euler_pole import EulerPole


import sys

import numpy as np
import pyproj

from mintpy.objects.constants import EARTH_RADIUS

################################################################################################
## global variables
MAS2RAD      = np.pi/3600000/180    #  1 mas (milliarcsecond) = yy radian
MASY2DMY     = 1e6 / 3600000        #  1 mas per year         = xx degree per million year


## Define the Euler pole class

class EulerPole:
    """ Class to compute velocity from given Euler pole vector
    EXAMPLE:
        1. Build an Euler pole object given pole and rotation rate
           (e.g., Eurasian; check ITRF2014-PMM defined in Altamimi et al. (2017)):
            Pole = EulerPole(lat=55.070, lon=-99.095, rotationRate=0.939, unit='mas/yr')
        2. Print the pole information:
            Pole.printMsg()
        3. Compute velocity on the ground given location(s) of interest. Assume flat topography and WGS84 ellipse.
            vx, vy, vz = Pole.getVelocityXYZ(lats, lons, alt=0.0, ellps=True)   # in ECEF xyz coordinate
            ve, vn, vu = Pole.getVelocityENU(lats, lons, alt=0.0, ellps=True)   # in local ENU coordinate
            azi, v_abs = Pole.getVelocityAzi(lats, lons, alt=0.0, ellps=True)   # in local ENU, get the azimuth and the absolute velocity
    """

    def __init__(self, wx=None, wy=None, wz=None, lat=None, lon=None, rotationRate=None, unit='mas/yr'):
        EXAMPLE = """Initialize the Euler pole
        Specify either [wx, wy,wz] or [lat, lon, rotationRate]
        Parameters: (1) Euler vector: [wx, wy, wz]
                        wx           -     float, angular velocity in x-axis      [unit]
                        wy           -     float, angular velocity in y-axis      [unit]
                        wz           -     float, angular velocity in z-axis      [unit]
                    (2) Pole and rotationRate: [lat, lon, omega]
                        lat          -     float, Euler pole latitude             [degree]
                        lon          -     float, Euler pole longitude            [degree]
                        rotationRate -     float, magnitude of angular velocity   [unit]
                      CONVENTION: positive = rotation vector pointing outward from the center of sphere;
                                           = counterclockwise rotation when looking from outside of sphere
                    unit             -  str,   ['mas/yr', 'deg/Ma']
        Example:    equivalent ways to describe the Eurasian plate in the ITRF2014 plate motion model
                        EulerPole(wx=-0.085,  wy=-0.531,   wz=0.770)
                        EulerPole(wx=-0.024,  wy=-0.148,   wz=0.214,           unit='deg/Ma')
                        EulerPole(lat=55.070, lon=-99.095, rotationRate=0.939)
                        EulerPole(lat=55.070, lon=-99.095, rotationRate=0.261, unit='deg/Ma')
                        EulerPole(lat=-55.070, lon=80.905, rotationRate=-0.939)     # flip the sign!
        """
        # check unit
        if unit.lower().startswith('mas'):
            unit = 'mas/yr'
        elif unit.lower().startswith('deg'):
            unit = 'deg/Ma'
        else:
            unit = 'mas/yr'
            print(f'Unrecognized unit, use defualt: {unit}')

        # input Euler vector or pole
        if all([wx, wy, wz]):
            if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
                wx /= MASY2DMY
                wy /= MASY2DMY
                wz /= MASY2DMY
            lat, lon, rotationRate = cart2sph(wx, wy, wz)
        elif all([lat, lon, rotationRate]):
            if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
                rotationRate /= MASY2DMY
            wx, wy, wz = sph2cart(lat, lon, r=rotationRate)
        else:
            print(EXAMPLE)
            sys.exit('User input Euler vector is incomplete!')

        # save to self variable
        self.plat          = lat            # Euler pole lat      [degree]
        self.plon          = lon            # Euler pole lon      [degree]
        self.rotationRate  = rotationRate   # Rotation rate       [mas/yr]
        self.wx            = wx             # Angular velocity x  [mas/yr]
        self.wy            = wy             # Angular velocity y  [mas/yr]
        self.wz            = wz             # Angular velocity z  [mas/yr]


    def __add__(self, other):
        """Add two Euler pole objects
            pole1 = EulerPole(...)
            pole2 = EulerPole(...)
            pole3 = pol2 + pol1
        """
        new_wx = self.wx + other.wx
        new_wy = self.wy + other.wy
        new_wz = self.wz + other.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def __sub__(self, other):
        """Subtract two Euler pole objects
            pole1 = EulerPole(...)
            pole2 = EulerPole(...)
            pole3 = pol2 - pol1
        """
        new_wx = self.wx - other.wx
        new_wy = self.wy - other.wy
        new_wz = self.wz - other.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def __neg__(self):
        """Negative of an Euler pole object
            pole1 = EulerPole(...)
            pole2 = -pol1
        """
        new_wx = -self.wx
        new_wy = -self.wy
        new_wz = -self.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def printMsg(self):
        """Print the Euler pole information
            pole = EulerPole(...)
            pole.printMsg()
        """
        print('\n------------------ Euler Pole description ------------------')
        print('Spherical expression:')
        print(f'   Pole Latitude  : {self.plat:9.4f} DEG')
        print(f'   Pole Longitude : {self.plon:9.4f} DEG')
        print(f'   Rotation rate  : {self.rotationRate * MASY2DMY:9.4f} DEG/MA = {self.rotationRate:9.4f} MAS/Y')
        print('Cartesian expression (angular velocity vector):')
        print(f'   wx             : {self.wx * MASY2DMY:9.4f} DEG/MA = {self.wx:9.4f} MAS/Y')
        print(f'   wy             : {self.wy * MASY2DMY:9.4f} DEG/MA = {self.wy:9.4f} MAS/Y')
        print(f'   wz             : {self.wz * MASY2DMY:9.4f} DEG/MA = {self.wz:9.4f} MAS/Y')
        print('------------------------------------------------------------\n')


    def getVelocityAzi(self, lat, lon, alt=0.0, ellps=True):
        """Calculates the azimuth (in degrees) and rate of plate motion
        (in millimeters per year) at a given point.
        Parameters:  lat     - Latitude of the point                [DEG]
                     lon     - Longitude of the point               [DEG]
        Returns:     azi     - Azimuth angle clockwise from north   [DEG]
                     v_abs   - Velocity                             [mm/yr]
        """
        ve, vn, vu = self.getVelocityENU(lat, lon, alt=alt, ellps=ellps)
        azi        = azimuth_from_EN(ve, vn)
        v_abs      = np.sqrt(ve**2 + vn**2 + vu**2)
        return azi, v_abs


    def getVelocityENU(self, lat, lon, alt=0.0, ellps=True):
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
        vx, vy, vz = self.getVelocityXYZ(lat, lon, alt=alt, ellps=ellps)
        ve, vn, vu = rotate_xyz_enu(lat, lon, x=vx, y=vy, z=vz)
        return ve, vn, vu


    def getVelocityXYZ(self, lat, lon, alt=0.0, ellps=True):
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
        ## Local coordinates handling for location(s) of interest
        if not ellps:
            # a perfect sphere
            print(f'Assume a perfect spherical Earth with radius={EARTH_RADIUS} m')
            x, y, z = sph2cart(lat, lon, EARTH_RADIUS)   # meter
        else:
            # WGS84 ellips; only supports uniform altitude now, but can change later
            print('Assume WGS84 ellipse from pyproj')
            x, y, z = wgslalo2xyz(lat, lon, alt)        # meter
        xyz = np.array([x, y, z], dtype=float)          # meter
        shp = xyz.shape

        ## Compute the cartesian linear velocity (i.e., ECEF) in km/year
        # Vxyz = omega x R_i, where R_i is location vector at pixel i
        omega = np.array([self.wx, self.wy, self.wz]) * MAS2RAD
        V     = np.cross(omega, xyz.T).T.reshape(shp)
        vx, vy, vz = V[0], V[1], V[2]          # meter/year
        return vx, vy, vz



################################################################################################
## Sub-functions for the math of Euler Pole and linear velocity
#  https://yuankailiu.github.io/assets/docs/Euler_pole_doc.pdf

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


def rotate_xyz_enu(lat, lon, x=None, y=None, z=None, e=None, n=None, u=None):
    """Matrix rotation to convert between global xyz coord (ECEF) and local ENU coord at given (lat, lon)
    Reference:
        Navipedia, https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications,
          Palo Alto, doi: 10.4236/ojapps.2015.54016. Page 145-156

    Parameters: lat     - 1D np.ndarray, latitude at location(s)        [degree]
                lon     - 1D np.ndarray, longitude at location(s)       [degree]
                x       - 1D np.ndarray, x component at location(s)     [e.g., length, velocity]
                y       - 1D np.ndarray, y component at location(s)     [e.g., length, velocity]
                z       - 1D np.ndarray, z component at location(s)     [e.g., length, velocity]
                e       - 1D np.ndarray, east component at location(s)  [e.g., length, velocity]
                n       - 1D np.ndarray, north component at location(s) [e.g., length, velocity]
                u       - 1D np.ndarray, up  component at location(s)   [e.g., length, velocity]
    Returns:
                Given x,y,z       => returns e,n,u
                Given e,n,u       => returns x,y,z
    """
    lat, lon = check_lalo(lat, lon, print_msg=True)
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # matrix rotation
    if any(x) and any(y) and any(z):
        # cart2enu
        e = - np.sin(lon) * x \
            + np.cos(lon) * y
        n = - np.sin(lat) * np.cos(lon) * x \
            - np.sin(lat) * np.sin(lon) * y \
            + np.cos(lat) * z
        u =   np.cos(lat) * np.cos(lon) * x \
            + np.cos(lat) * np.sin(lon) * y \
            + np.sin(lat) * z
        return e, n, u
    elif any(e) and any(n) and any(u):
        # enu2cart
        x = - np.sin(lon) * e \
            - np.cos(lon) * np.sin(lat) * n \
            + np.cos(lon) * np.cos(lat) * u
        y =   np.cos(lon) * e \
            - np.sin(lon) * np.sin(lat) * n \
            + np.sin(lon) * np.cos(lat) * u
        z =   np.cos(lat) * n \
            + np.sin(lat) * u
        return x, y, z
    else:
        sys.exit('User input (x,y,z) or (e,n,u) is incomplete!')



def wgslalo2xyz(lat, lon, alt):
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


def azimuth_from_EN(east, north):
    """Azimuth in degrees clockwise from North
    Parameter:  east  - float, eastward motion
                north - float, northward motion
    """
    azi = np.rad2deg(np.arctan2(east, north)) % 360
    return azi


def check_lalo(lat, lon, print_msg=False):
    """ ensure numpy array format
    """
    if isinstance(lat, (list, tuple, np.ndarray)):
        lat = np.array(lat).flatten()
        lon = np.array(lon).flatten()
        nla, nlo = lat.size, lon.size

        if not nla == nlo:
            raise ValueError(f'Inconsistent size between input lat ({nla}) and lon ({nlo})!')

    elif isinstance(lat, (int, float)):
        nla, nlo = 1, 1

    if print_msg:
        ## report how many points of interest
        print(f'number of points to compute: {nla}')

    return lat, lon
