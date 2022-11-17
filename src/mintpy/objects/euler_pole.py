"""Euler Pole class definition and utility functions."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Oct 2022                           #
############################################################
# Recommend import:
#     from mintpy.objects.euler_pole import EulerPole
#
# Reference:
#  Pichon, X. L., Francheteau, J. & Bonnin, J. Plate Tectonics; Developments in Geotectonics 6;
#    Hardcover – January 1, 1973. Page 28-29
#  Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications,
#    Palo Alto. DOI: 10.4236/ojapps.2015.54016. Page 145-156.
#  Navipedia, Transformations between ECEF and ENU coordinates. [Online].
#    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
#  Goudarzi, M. A., Cocard, M. & Santerre, R. (2014), EPC: Matlab software to estimate Euler
#    pole parameters, GPS Solutions, 18, 153-162, doi: 10.1007/s10291-013-0354-4


import numpy as np
import pyproj

from mintpy.constants import EARTH_RADIUS

# global variables
MAS2RAD = np.pi / 3600000 / 180    # 1 mas (milli arc second) = x radian
MASY2DMY = 1e6 / 3600000           # 1 mas per year = x degree per million year


####################################  EulerPole class begin  #############################################
# Define the Euler pole class
EXAMPLE = """Define an Euler pole:
  Method 1 - Use an Euler vector [wx, wy, wz]
             wx/y/z   - float, angular velocity in x/y/z-axis [mas/yr or deg/Ma]
  Method 2 - Use an Euler Pole lat/lon and rotation rate [lat, lon, rot_rate]
             lat/lon  - float, Euler pole latitude/longitude [degree]
             rot_rate - float, magnitude of the angular velocity [deg/Ma or mas/yr]
             1) define rotation vection as from the center of sphere to the pole lat/lon and outward;
             2) positive for conterclockwise rotation when looking from outside along the rotation vector.

  Example:
    # equivalent ways to describe the Eurasian plate in the ITRF2014 plate motion model
    EulerPole(wx=-0.085, wy=-0.531, wz=0.770, unit='mas/yr')
    EulerPole(wx=-0.024, wy=-0.148, wz=0.214, unit='deg/Ma')
    EulerPole(pole_lat=55.070, pole_lon=-99.095, rot_rate=0.939, unit='mas/yr')
    EulerPole(pole_lat=55.070, pole_lon=-99.095, rot_rate=0.261, unit='deg/Ma')
    EulerPole(pole_lat=-55.070, pole_lon=80.905, rot_rate=-0.939, unit='mas/yr')
"""

class EulerPole:
    """EulerPole object to compute velocity for a given tectonic plate.

    Example:
        # compute velocity of the Eurasia plate in ITRF2014-PMM from Altamimi et al. (2017)
        pole_obj = EulerPole(pole_lat=55.070, pole_lon=-99.095, rot_rate=0.939, unit='mas/yr')
        pole_obj.print_info()
        vx, vy, vz = pole_obj.get_velocity_xyz(lats, lons, alt=0.0) # in  ECEF xyz coordinate
        ve, vn, vu = pole_obj.get_velocity_enu(lats, lons, alt=0.0) # in local ENU coordinate
    """

    def __init__(self, wx=None, wy=None, wz=None, pole_lat=None, pole_lon=None, rot_rate=None,
                 unit='mas/yr', name=None):
        # check - unit
        if unit.lower().startswith('mas'):
            unit = 'mas/yr'

        elif unit.lower().startswith('deg'):
            unit = 'deg/Ma'
            # convert input deg/Ma to mas/yr for internal calculation
            wx = wx / MASY2DMY if wx else None
            wy = wy / MASY2DMY if wy else None
            wz = wz / MASY2DMY if wz else None
            rot_rate = rot_rate / MASY2DMY if rot_rate else None

        else:
            raise ValueError(f'Unrecognized rotation rate unit: {unit}! Use mas/yr or deg/Ma')

        # calculate Euler vector and pole
        if all([wx, wy, wz]):
            # calc Euler pole from vector
            pole_lat, pole_lon, rot_rate = cart2sph(wx, wy, wz)

        elif all([pole_lat, pole_lon, rot_rate]):
            # calc Euler vector from pole
            wx, wy, wz = sph2cart(pole_lat, pole_lon, r=rot_rate)

        else:
            raise ValueError(f'Incomplete Euler Pole input!\n{EXAMPLE}')

        # save member variables
        self.name = name
        self.poleLat = pole_lat   # Euler pole latitude   [degree]
        self.poleLon = pole_lon   # Euler pole longitude  [degree]
        self.rotRate = rot_rate   # angular rotation rate [mas/yr]
        self.wx = wx              # angular velocity x    [mas/yr]
        self.wy = wy              # angular velocity y    [mas/yr]
        self.wz = wz              # angular velocity z    [mas/yr]


    def __repr__(self):
        msg = f'{self.__class__.__name__}(name={self.name}, poleLat={self.poleLat}, poleLon={self.poleLon}, '
        msg += f'rotRate={self.rotRate}, wx={self.wx}, wy={self.wy}, wz={self.wz}, unit=mas/yr)'
        return msg


    def __add__(self, other):
        """Add two Euler pole objects.

        Example:
            pole1 = EulerPole(...)
            pole2 = EulerPole(...)
            pole3 = pol2 + pol1
        """
        new_wx = self.wx + other.wx
        new_wy = self.wy + other.wy
        new_wz = self.wz + other.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def __sub__(self, other):
        """Subtract two Euler pole objects.

        Example:
            pole1 = EulerPole(...)
            pole2 = EulerPole(...)
            pole3 = pol2 - pol1
        """
        new_wx = self.wx - other.wx
        new_wy = self.wy - other.wy
        new_wz = self.wz - other.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def __neg__(self):
        """Negative of an Euler pole object.

        Example:
            pole1 = EulerPole(...)
            pole2 = -pol1
        """
        new_wx = -self.wx
        new_wy = -self.wy
        new_wz = -self.wz
        return EulerPole(wx=new_wx, wy=new_wy, wz=new_wz)


    def print_info(self):
        """Print the Euler pole information.
        """
        # maximum digit
        vals = [self.poleLat, self.poleLon, self.rotRate, self.wx, self.wy, self.wz]
        md = len(str(int(np.max(np.abs(vals))))) + 5
        md += 1 if any(x < 0 for x in vals) else 0

        print('\n------------------ Euler Pole description ------------------')
        print('Spherical expression:')
        print(f'   Pole Latitude  : {self.poleLat:{md}.4f} deg')
        print(f'   Pole Longitude : {self.poleLon:{md}.4f} deg')
        print(f'   Rotation rate  : {self.rotRate * MASY2DMY:{md}.4f} deg/Ma   = {self.rotRate:{md}.4f} mas/yr')
        print('Cartesian expression (angular velocity vector):')
        print(f'   wx             : {self.wx * MASY2DMY:{md}.4f} deg/Ma   = {self.wx:{md}.4f} mas/yr')
        print(f'   wy             : {self.wy * MASY2DMY:{md}.4f} deg/Ma   = {self.wy:{md}.4f} mas/yr')
        print(f'   wz             : {self.wz * MASY2DMY:{md}.4f} deg/Ma   = {self.wz:{md}.4f} mas/yr')
        print('------------------------------------------------------------\n')


    def get_velocity_xyz(self, lat, lon, alt=0.0, ellps=True, print_msg=True):
        """Compute cartesian velocity (vx, vy, vz) of the Euler Pole at point(s) of interest.

        Parameters: lat   - float / 1D/2D np.ndarray, points of interest (latitude)  [degree]
                    lon   - float / 1D/2D np.ndarray, points of interest (longitude) [degree]
                    alt   - float / 1D/2D np.ndarray, points of interest (altitude)      [meter]
                    ellps - bool, consider ellipsoidal Earth projection
        Returns:    vx    - float / 1D/2D np.ndarray, ECEF x linear velocity [meter/year]
                    vy    - float / 1D/2D np.ndarray, ECEF y linear velocity [meter/year]
                    vz    - float / 1D/2D np.ndarray, ECEF z linear velocity [meter/year]
        """
        # check input lat/lon data type (scalar / array) and shape
        poi_shape = lat.shape if isinstance(lat, np.ndarray) else None

        # convert lat/lon into x/y/z
        # Note: the conversion assumes either a spherical or spheroidal Earth, tests show that
        # using a ellipsoid as defined in WGS84 produce results closer to the UNAVCO website
        # calculator, which also uses the WGS84 ellipsoid.
        if ellps:
            if print_msg:
                print('assume a spheroidal Earth as defined in WGS84')
            x, y, z = coord_llh2xyz(lat, lon, alt)
        else:
            if print_msg:
                print(f'assume a spherical Earth with radius={EARTH_RADIUS} m')
            x, y, z = sph2cart(lat, lon, alt+EARTH_RADIUS)

        # ensure matrix is flattened
        if poi_shape is not None:
            x = x.flatten()
            y = y.flatten()
            z = z.flatten()

        # compute the cartesian linear velocity (i.e., ECEF) in meter/year as:
        #
        #     V_xyz = Omega x R_i
        #
        # where R_i is location vector at point i

        xyz = np.array([x, y, z], dtype=np.float32)
        omega = np.array([self.wx, self.wy, self.wz]) * MAS2RAD
        vx, vy, vz = np.cross(omega, xyz.T).T.reshape(xyz.shape)

        # reshape to the original shape of lat/lon
        if poi_shape is not None:
            vx = vx.reshape(poi_shape)
            vy = vy.reshape(poi_shape)
            vz = vz.reshape(poi_shape)

        return vx, vy, vz


    def get_velocity_enu(self, lat, lon, alt=0.0, ellps=True, print_msg=True):
        """Compute the spherical velocity (ve, vn, vu) of the Euler Pole at point(s) of interest.

        Parameters: lat   - float / 1D/2D np.ndarray, points of interest (latitude)  [degree]
                    lon   - float / 1D/2D np.ndarray, points of interest (longitude) [degree]
                    alt   - float / 1D/2D np.ndarray, points of interest (altitude) [meter]
                    ellps - bool, consider ellipsoidal Earth projection
        Returns:    ve    - float / 1D/2D np.ndarray, east  linear velocity [meter/year]
                    vn    - float / 1D/2D np.ndarray, north linear velocity [meter/year]
                    vu    - float / 1D/2D np.ndarray, up    linear velocity [meter/year]
        """
        # calculate ECEF velocity
        vx, vy, vz = self.get_velocity_xyz(lat, lon, alt=alt, ellps=ellps, print_msg=print_msg)

        # convert ECEF to ENU velocity via matrix rotation: V_enu = T * V_xyz
        ve, vn, vu = transform_xyz_enu(lat, lon, x=vx, y=vy, z=vz)

        # enforce zero vertical velocitpy when ellps=False
        # to avoid artefacts due to numerical precision
        if not ellps:
            if isinstance(lat, np.ndarray):
                vu[:] = 0
            else:
                vu = 0

        return ve, vn, vu

####################################  EulerPole class end  ###############################################


####################################  Utility functions  #################################################
# Utility functions for the math/geometry operations of Euler Pole and linear velocity
# reference: https://yuankailiu.github.io/assets/docs/Euler_pole_doc.pdf
def cart2sph(rx, ry, rz):
    """Convert cartesian coordinates to spherical.

    Parameters: rx/y/z  - float / np.ndarray, angular distance in X/Y/Z direction [any units of distance]
    Returns:    lat/lon - float / np.ndarray, latitude / longitude  [degree]
                r       - float / np.ndarray, radius [same unit as rx/y/z]
    Examples:
        # convert xyz coord to spherical coord
        lat, lon, r = cart2sph(x, y, z)
        # convert Euler vector (in cartesian) to Euler pole (in spherical)
        pole_lat, pole_lon, rot_rate = cart2sph(wx, wy, wz)
    """
    r = np.sqrt(rx**2 + ry**2 + rz**2)
    lat = np.rad2deg(np.arcsin(rz / r))
    lon = np.rad2deg(np.arctan2(ry, rx))
    return lat, lon, r


def sph2cart(lat, lon, r=1):
    """Convert spherical coordinates to cartesian.

    Parameters: lat/lon - float / np.ndarray, latitude / longitude [degree]
                r       - float / np.ndarray, radius [any units of angular distance]
    Returns:    rx/y/z  - float / np.ndarray, angular distance in X/Y/Z direction [same unit as r]
    Examples:
        # convert spherical coord to xyz coord
        x, y, z = sph2cart(lat, lon, r=radius)
        # convert Euler pole (in spherical) to Euler vector (in cartesian)
        wx, wy, wz = sph2cart(pole_lat, pole_lon, r=rot_rate)
    """
    rx = r * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))
    ry = r * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon))
    rz = r * np.sin(np.deg2rad(lat))
    return rx, ry, rz


def coord_llh2xyz(lat, lon, alt):
    """Convert coordinates from WGS84 lat/long/hgt to ECEF x/y/z.

    Parameters: lat   - float / list(float) / np.ndarray, latitude  [degree]
                lon   - float / list(float) / np.ndarray, longitude [degree]
                alt   - float / list(float) / np.ndarray, altitude  [meter]
    Returns:    x/y/z - float / list(float) / np.ndarray, ECEF coordinate [meter]
    """
    # ensure same type between alt and lat/lon
    if isinstance(lat, np.ndarray) and not isinstance(alt, np.ndarray):
        alt *= np.ones_like(lat)
    elif isinstance(lat, list) and not isinstance(alt, list):
        alt = [alt] * len(lat)

    # construct pyproj transform object
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    )

    # apply coordinate transformation
    x, y, z = transformer.transform(lon, lat, alt, radians=False)

    return x, y, z


def transform_xyz_enu(lat, lon, x=None, y=None, z=None, e=None, n=None, u=None):
    """Transform between ECEF (global xyz) and ENU at given locations (lat, lon) via matrix rotation.

    Reference:
        Navipedia, https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications,
          Palo Alto, doi: 10.4236/ojapps.2015.54016. Page 145-156

    Parameters: lat/lon - float / np.ndarray, latitude/longitude      at location(s) [degree]
                x/y/z   - float / np.ndarray, x/y/z         component at location(s) [e.g., length, velocity]
                e/n/u   - float / np.ndarray, east/north/up component at location(s) [e.g., length, velocity]
    Returns:    e/n/u if given x/y/z
                x/y/z if given e/n/u
    """
    # convert the unit from degree to radian
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # transformation via matrix rotation:
    #     V_enu = T * V_xyz
    #     V_xyz = T^-1 * V_enu
    #
    # Equilevant 3D matrix code is as below:
    #     V_enu = np.diagonal(
    #         np.matmul(
    #             T.reshape([-1,3]),
    #             V_xyz.T,
    #         ).reshape([3, npts, npts], order='F'),
    #         axis1=1,
    #         axis2=2,
    #     ).T
    # To avoid this complex matrix operation above, we calculate for each element as below:

    if all(i is not None for i in [x, y, z]):
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

    elif all(i is not None for i in [e, n, u]):
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
        raise ValueError('Input (x,y,z) or (e,n,u) is NOT complete!')
