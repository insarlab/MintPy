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

import collections
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pyproj
from shapely import geometry

import mintpy
from mintpy.cli.plate_motion import GSRM_PMM, MORVEL56_PMM
from mintpy.constants import EARTH_RADIUS

# global variables
MAS2RAD = np.pi / 3600000 / 180    # 1 mas (milli arc second) = x radian
MASY2DMY = 1e6 / 3600000           # 1 mas per year = x degree per million year

# path for plate boundary files
PLATE_BOUNDARY_FILE = {
    'GSRM'     : os.path.join(mintpy.__path__[0], 'data/plate_boundary/GSRM/GSRM_plate_outlines.gmt'),
    'MORVEL56' : os.path.join(mintpy.__path__[0], 'data/plate_boundary/MORVEL56/All_boundaries'),
}

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



####################################  Utility for plotting  ##############################################
# Utility for plotting the plate motion on a globe
# check usage: https://github.com/yuankailiu/utils/blob/main/notebooks/PMM_plot.ipynb

def read_plate_attributes(filename):
    """read the plate names and abbreviated names (obsolete)
    Parameters:
        filename - filename of the csv file
    Returns:
        pDict    - a dictionary that contains abbreviations and Euler pole attributes (not used)
    """
    pDict = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('Plate,'):
                continue
            key  = line.split('\n')[0]
            name = key.split(',')[0]
            abbv = key.split(',')[1]
            lat  = key.split(',')[2]
            lon  = key.split(',')[3]
            rate = key.split(',')[4]
            pDict[abbv] = (name, lat, lon, rate)
    return pDict


def read_plate_outlines(filename, order='lola'):
    """Load the plate boundary text files
    Paramters:
        filename - filename of the boundary text file
        order    - the order of columns, 'lalo' or 'lola', default: (lola)
    Returns:
        Bnds     - a dictionary that contains a list of vertices of the plate polygon (lat, lon)
    """
    datatype = filename.split('data/plate_boundary/')[-1].split('/')[0]
    pDict = {}
    if 'GSRM' in datatype:
        for key, val in GSRM_PMM.items():
            pDict[val.Abbrev] = key
    elif 'MORVEL56' in datatype:
        for key, val in MORVEL56_PMM.items():
            pDict[val.Abbrev] = key
    Bnds = {}
    with open(filename) as f:
        lines = f.readlines()
        key, vertices = None, None
        for line in lines:
            if line.startswith('> ') or line.startswith('# ') or len(line.split())==1:
                if key and vertices:
                    name = pDict[key]
                    Bnds[name] = np.array(vertices)
                if   line.startswith('> '):
                    key = line.split('> ')[1]
                elif line.startswith('# '):
                    key = line.split('# ')[1]
                else:
                    key = str(line)
                if key.endswith('\n'):
                    key = key.split('\n')[0]
                vertices = []
            else:
                if order == 'lalo':
                    vert = np.array(line.split()).astype(float)
                elif order == 'lola':
                    vert = np.flip(np.array(line.split()).astype(float))
                vertices.append(vert)
    return Bnds


def generate_sample_in_polygon(polygon, nx=10, ny=10):
    """Make a set of points inside the defined sphericalpolygon object
    Parameters:
        polygon - the shapely polygon object
        nx      - the number of points in the x (lon) direction
        ny      - the number of points in the y (lat) direction
    Returns:
        sample  - tuple that contains the coordinate (lats, lons) of the generated points
    """
    vlats = np.array(polygon.exterior.coords)[:,0]
    vlons = np.array(polygon.exterior.coords)[:,1]
    lats = np.linspace(np.min(vlats), np.max(vlats), ny)
    lons = np.linspace(np.min(vlons), np.max(vlons), nx)
    lats, lons = np.meshgrid(lats, lons)
    sample_candidates = np.vstack([lats.flatten(), lons.flatten()]).T

    # check which points are inside the polygon
    in_idx = np.full(len(sample_candidates), False)
    for i, candi in enumerate(sample_candidates):
        if polygon.contains(geometry.Point(candi[0], candi[1])):
            in_idx[i] = True

    in_idx = np.array(in_idx).reshape(lats.shape)
    lons   = lons[in_idx]
    lats   = lats[in_idx]

    Point = collections.namedtuple('Points', 'lats lons')
    sample = Point(lats=lats, lons=lons)
    return sample


def plot_map(polygon, plate_obj, center_lat=None, center_lon=None, qscale=100, unit_q=5, zoom=100, dpi=100, outfile=None, **kwargs):
    """Plot the globe map wityh plate boundary, quivers on some points
    Parameters:
        polygon   - a sphericalpolygon object created from an instance of class Sphericalpolygon
                    (https://pypi.org/project/sphericalpolygon/)
        plate_obj - an EulerPole object created from an instance of class EulerPole
                    (mintpy.objects.euler_pole)
        center_lat - center the map at this latitute
        center_lon - center the map at this longitude
        qscale     - scaling factor of the quiver
        unit_q     - what value of quiver legend
        zoom       - zoom the globe to see a narrower region
        dpi        - output figure dpi
        outfile    - output figure name
        kwargs     - dictionary for plotting

    Returns:
        Draw a map and save it
    """
    kwargs['c_ocean']       = kwargs.get('c_ocean',     'w')
    kwargs['c_land']        = kwargs.get('c_land',      'lightgray')
    kwargs['c_plate']       = kwargs.get('c_plate',     'mistyrose')
    kwargs['lw_coast']      = kwargs.get('lw_coast',    0.5)
    kwargs['lw_pbond']      = kwargs.get('lw_pbond',    1)
    kwargs['lc_pbond']      = kwargs.get('lc_pbond',    'coral')
    kwargs['alpha_plate']   = kwargs.get('alpha_plate', 0.4)
    kwargs['grid_ls']       = kwargs.get('grid_ls',     '--')
    kwargs['grid_lw']       = kwargs.get('grid_lw',     0.3)
    kwargs['grid_lc']       = kwargs.get('grid_lc',     'gray')
    kwargs['Nq']            = kwargs.get('Nq', 10)
    kwargs['pts']           = kwargs.get('pts',         None)
    kwargs['pts_marker']    = kwargs.get('pts_marker',  '^')
    kwargs['pts_ms']        = kwargs.get('pts_ms',      20)
    kwargs['pts_mfc']       = kwargs.get('pts_mfc',     'r')
    kwargs['pts_mec']       = kwargs.get('pts_mec',     'k')
    kwargs['pts_mew']       = kwargs.get('pts_mew',     1)


    # center the map
    if (center_lat is None) or (center_lon is None):
        center_lat, center_lon = np.array(polygon.centroid)

    # make a base map from cartopy
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(projection=ccrs.NearsidePerspective(center_lon, center_lat, satellite_height=6.5e7/zoom))
    ax.set_global()
    ax.gridlines(color=kwargs['grid_lc'], linestyle=kwargs['grid_ls'], xlocs=np.arange(-180,180,30), ylocs=np.linspace(-80,80,10), linewidth=kwargs['grid_lw'])
    ax.add_feature(cfeature.OCEAN, color=kwargs['c_ocean'])
    ax.add_feature(cfeature.LAND,  color=kwargs['c_land'])
    ax.add_feature(cfeature.COASTLINE, linewidth=kwargs['lw_coast'])

    # add the plate polygon
    if polygon:
        vlats = np.array(polygon.exterior.coords)[:,0]
        vlons = np.array(polygon.exterior.coords)[:,1]
        ax.plot(vlons, vlats, color=kwargs['lc_pbond'], transform=ccrs.Geodetic(), linewidth=kwargs['lw_pbond'])
        ax.fill(vlons, vlats, color=kwargs['c_plate'],  transform=ccrs.Geodetic(), alpha=kwargs['alpha_plate'])

        # compute the plate motion from Euler rotation
        if plate_obj:
            # select some points inside the polygon and calc PMM, ENU
            sample = generate_sample_in_polygon(polygon, nx=kwargs['Nq'], ny=kwargs['Nq'])

            ve, vn, vu = plate_obj.get_velocity_enu(lat=sample.lats, lon=sample.lons, alt=0.0, ellps=True, print_msg=True)
            ve *= 1e3
            vn *= 1e3
            vu *= 1e3
            norm = np.sqrt(ve**2 + vn**2)

            # correcting for "East" further toward polar region; re-normalize ve, vn
            ve /= np.cos(np.deg2rad(sample.lats))
            renorm = np.sqrt(ve**2 + vn**2)/norm
            ve /= renorm
            vn /= renorm

            # ---------- plot inplate dots, vectors --------------
            ax.scatter(sample.lons, sample.lats, marker='.', s=20, color='k', transform=ccrs.PlateCarree())
            q = ax.quiver(sample.lons, sample.lats, ve,  vn, transform=ccrs.PlateCarree(), scale=qscale, width=.0075, color='coral', angles="xy")
            ax.set_title('  ', pad=10) # put an empty title for extra whitepace at the top
            ax.quiverkey(q, X=0.3, Y=0.9, U=unit_q, label=f'{unit_q} mm/yr', labelpos='E', coordinates='figure')

    # add custom points (e.g., show some points of interest)
    if not kwargs['pts'] is None:
        ax.scatter(kwargs['pts'][1], kwargs['pts'][0], marker=kwargs['pts_marker'], s=kwargs['pts_ms'],
                   fc=kwargs['pts_mfc'], ec=kwargs['pts_mec'], lw=kwargs['pts_mew'], transform=ccrs.PlateCarree())
    # save figure
    if outfile:
        print('save figure to file:', outfile)
        fig.savefig(outfile, bbox_inches='tight', transparent=True, dpi=300)

    plt.show()
