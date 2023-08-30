"""Class / utilities for Euler Pole / plate motion (models)."""
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

import matplotlib.pyplot as plt
import numpy as np
import pyproj
from cartopy import crs as ccrs, feature as cfeature
from shapely import geometry

import mintpy
from mintpy.constants import EARTH_RADIUS

# global variables
MAS2RAD = np.pi / 3600000 / 180    # 1 mas (milli arc second) = x radian
MASY2DMY = 1e6 / 3600000           # 1 mas per year = x degree per million year

#################################  Plate boundary files  #######################################
PLATE_BOUNDARY_FILE = {
    'GSRM'   : os.path.join(mintpy.__path__[0], 'data/plate_boundary/GSRM/plate_outlines.lola'),
    'MORVEL' : os.path.join(mintpy.__path__[0], 'data/plate_boundary/MORVEL/plate_outlines.lalo'),
}

#################################  Plate Motion Models  ########################################
# Later will be moved to a separate script `pmm.py` in plate motion package

# 1). ITRF2014-PMM defined in Altamimi et al. (2017)
# Reference frame: ITRF2014
Tag = collections.namedtuple('Tag', 'Abbrev num_site omega_x omega_y omega_z omega wrms_e wrms_n')
ITRF2014_PMM = {
    'Antartica'     : Tag('ANTA'  ,   7,  -0.248,  -0.324,   0.675,  0.219,  0.20,  0.16),
    'Arabia'        : Tag('ARAB'  ,   5,   1.154,  -0.136,   1.444,  0.515,  0.36,  0.43),
    'Australia'     : Tag('AUST'  ,  36,   1.510,   1.182,   1.215,  0.631,  0.24,  0.20),
    'Eurasia'       : Tag('EURA'  ,  97,  -0.085,  -0.531,   0.770,  0.261,  0.23,  0.19),
    'India'         : Tag('INDI'  ,   3,   1.154,  -0.005,   1.454,  0.516,  0.21,  0.21),
    'Nazca'         : Tag('NAZC'  ,   2,  -0.333,  -1.544,   1.623,  0.629,  0.13,  0.19),
    'NorthAmerica'  : Tag('NOAM'  ,  72,   0.024,  -0.694,  -0.063,  0.194,  0.23,  0.28),
    'Nubia'         : Tag('NUBI'  ,  24,   0.099,  -0.614,   0.733,  0.267,  0.28,  0.36),
    'Pacific'       : Tag('PCFC'  ,  18,  -0.409,   1.047,  -2.169,  0.679,  0.36,  0.31),
    'SouthAmerica'  : Tag('SOAM'  ,  30,  -0.270,  -0.301,  -0.140,  0.119,  0.34,  0.35),
    'Somalia'       : Tag('SOMA'  ,   3,  -0.121,  -0.794,   0.884,  0.332,  0.32,  0.30),
}
PMM_UNIT = {
    'omega'   : 'deg/Ma',  # degree per megayear or one-million-year
    'omega_x' : 'mas/yr',  # milli-arcsecond per year
    'omega_y' : 'mas/yr',  # milli-arcsecond per year
    'omega_z' : 'mas/yr',  # milli-arcsecond per year
    'wrms_e'  : 'mm/yr',   # milli-meter per year, weighted root mean scatter
    'wrms_n'  : 'mm/yr',   # milli-meter per year, weighted root mean scatter
}


# 2). GSRMv2.1 defined in Kreemer et al. (2014)
# Reference frame: IGS08
# (unit: Lat: °N; Lon: °E; omega: °/Ma)
Tag = collections.namedtuple('Tag', 'Abbrev Lat Lon omega')
GSRM_V21_PMM = {
    'Africa'          : Tag('AF'  , 49.66   ,  -78.08   , 0.285),
    'Amur'            : Tag('AM'  , 61.64   ,  -101.29  , 0.287),
    'Antarctica'      : Tag('AN'  , 60.08   ,  -120.14  , 0.234),
    'Arabia'          : Tag('AR'  , 51.12   ,  -19.87   , 0.484),
    'AegeanSea'       : Tag('AS'  , 47.78   ,  59.86    , 0.253),
    'Australia'       : Tag('AU'  , 33.31   ,  36.38    , 0.639),
    'BajaCalifornia'  : Tag('BC'  , -63.04  ,  104.02   , 0.640),
    'Bering'          : Tag('BG'  , -40.62  ,  -53.84   , 0.333),
    'Burma'           : Tag('BU'  , -4.38   ,  -76.17   , 2.343),
    'Caribbean'       : Tag('CA'  , 37.84   ,  -96.49   , 0.290),
    'Caroline'        : Tag('CL'  , -76.41  ,  30.22    , 0.552),
    'Cocos'           : Tag('CO'  , 27.21   ,  -124.02  , 1.169),
    'Capricorn'       : Tag('CP'  , 42.13   ,  24.28    , 0.622),
    'Danakil'         : Tag('DA'  , 21.80   ,  36.05    , 2.497),
    'Easter'          : Tag('EA'  , 25.14   ,  67.55    , 11.331),
    'Eurasia'         : Tag('EU'  , 55.38   ,  -95.41   , 0.271),
    'Galapagos'       : Tag('GP'  , 2.83    ,  81.26    , 5.473),
    'Gonave'          : Tag('GV'  , 23.89   ,  -84.86   , 0.476),
    'India'           : Tag('IN'  , 50.95   ,  -8.00    , 0.524),
    'JuandeFuca'      : Tag('JF'  , -37.71  ,  59.44    , 0.977),
    'JuanFernandez'   : Tag('JZ'  , 34.33   ,  70.76    , 22.370),
    'Lwandle'         : Tag('LW'  , 52.20   ,  -60.68   , 0.273),
    'Mariana'         : Tag('MA'  , 11.20   ,  142.82   , 2.165),
    'NorthAmerica'    : Tag('NA'  , 2.19    ,  -83.75   , 0.219),
    'NorthBismarck'   : Tag('NB'  , -30.20  ,  135.30   , 1.201),
    'Niuafo`ou'       : Tag('NI'  , -3.51   ,  -174.04  , 3.296),
    'Nazca'           : Tag('NZ'  , 49.05   ,  -102.13  , 0.611),
    'Okhotsk'         : Tag('OK'  , 28.80   ,  -90.91   , 0.209),
    'Okinawa'         : Tag('ON'  , 39.11   ,  145.94   , 1.361),
    'Pacific'         : Tag('PA'  , -63.09  ,  109.63   , 0.663),
    'Panama'          : Tag('PM'  , 16.55   ,  -84.30   , 1.392),
    'PuertoRico'      : Tag('PR'  , 27.81   ,  -81.51   , 0.502),
    'PhilippineSea'   : Tag('PS'  , -46.62  ,  -28.39   , 0.895),
    'Rivera'          : Tag('RI'  , 20.27   ,  -107.10  , 4.510),
    'Rovuma'          : Tag('RO'  , 51.72   ,  -69.88   , 0.270),
    'SouthAmerica'    : Tag('SA'  , -14.10  ,  -117.86  , 0.123),
    'SouthBismarck'   : Tag('SB'  , 6.91    ,  -32.41   , 6.665),
    'Scotia'          : Tag('SC'  , 23.02   ,  -98.78   , 0.122),
    'Sinai'           : Tag('SI'  , 53.34   ,  -7.27    , 0.476),
    'Sakishima'       : Tag('SK'  , 27.31   ,  128.68   , 7.145),
    'Shetland'        : Tag('SL'  , 66.05   ,  134.03   , 1.710),
    'Somalia'         : Tag('SO'  , 47.59   ,  -94.36   , 0.346),
    'SolomonSea'      : Tag('SS'  , -3.33   ,  130.60   , 1.672),
    'Satunam'         : Tag('ST'  , 36.68   ,  135.30   , 2.846),
    'Sunda'           : Tag('SU'  , 51.11   ,  -91.75   , 0.350),
    'Sandwich'        : Tag('SW'  , -30.11  ,  -35.58   , 1.369),
    'Tonga'           : Tag('TO'  , 26.38   ,  4.27     , 8.853),
    'Victoria'        : Tag('VI'  , 44.96   ,  -102.19  , 0.330),
    'Woodlark'        : Tag('WL'  , -1.62   ,  130.63   , 1.957),
    'Yangtze'         : Tag('YA'  , 64.76   ,  -109.19  , 0.335),
}


# 3). NNR-MORVEL56 defined in Argus et al. (2011)
# (unit: Lat: °N; Lon: °E; omega: °/Ma)
# Note: we use "NU", instead of "nb" from Argus et al. (2011), for Nubia plate
#   to distinguish from "NB" for North Bismarck plate.
Tag = collections.namedtuple('Tag', 'Abbrev Lat Lon omega')
NNR_MORVEL56_PMM = {
    'Amur'            : Tag('AM'  , 63.17   , -122.82   , 0.297),
    'Antarctica'      : Tag('AN'  , 65.42   , -118.11   , 0.250),
    'Arabia'          : Tag('AR'  , 48.88   , -8.49     , 0.559),
    'Australia'       : Tag('AU'  , 33.86   , 37.94     , 0.632),
    'Capricorn'       : Tag('CP'  , 44.44   , 23.09     , 0.608),
    'Caribbean'       : Tag('CA'  , 35.20   , -92.62    , 0.286),
    'Cocos'           : Tag('CO'  , 26.93   , -124.31   , 1.198),
    'Eurasia'         : Tag('EU'  , 48.85   , -106.50   , 0.223),
    'India'           : Tag('IN'  , 50.37   , -3.29     , 0.544),
    'JuandeFuca'      : Tag('JF'  , -38.31  , 60.04     , 0.951),
    'Lwandle'         : Tag('LW'  , 51.89   , -69.52    , 0.286),
    'Macquarie'       : Tag('MQ'  , 49.19   , 11.05     , 1.144),
    'Nazca'           : Tag('NZ'  , 46.23   , -101.06   , 0.696),
    'NorthAmerica'    : Tag('NA'  , -4.85   , -80.64    , 0.209),
    'Nubia'           : Tag('NU'  , 47.68   , -68.44    , 0.292),
    'Pacific'         : Tag('PA'  , -63.58  , 114.70    , 0.651),
    'PhilippineSea'   : Tag('PS'  , -46.02  , -31.36    , 0.910),
    'Rivera'          : Tag('RI'  , 20.25   , -107.29   , 4.536),
    'Sandwich'        : Tag('SW'  , -29.94  , -36.87    , 1.362),
    'Scotia'          : Tag('SC'  , 22.52   , -106.15   , 0.146),
    'Somalia'         : Tag('SM'  , 49.95   , -84.52    , 0.339),
    'SouthAmerica'    : Tag('SA'  , -22.62  , -112.83   , 0.109),
    'Sunda'           : Tag('SU'  , 50.06   , -95.02    , 0.337),
    'Sur'             : Tag('SR'  , -32.50  , -111.32   , 0.107),
    'Yangtze'         : Tag('YZ'  , 63.03   , -116.62   , 0.334),
    'AegeanSea'       : Tag('AS'  , 19.43   , 122.87    , 0.124),
    'Altiplano'       : Tag('AP'  , -6.58   , -83.98    , 0.488),
    'Anatolia'        : Tag('AT'  , 40.11   , 26.66     , 1.210),
    'BalmoralReef'    : Tag('BR'  , -63.74  , 142.06    , 0.490),
    'BandaSea'        : Tag('BS'  , -1.49   , 121.64    , 2.475),
    'BirdsHead'       : Tag('BH'  , -40.00  , 100.50    , 0.799),
    'Burma'           : Tag('BU'  , -6.13   , -78.10    , 2.229),
    'Caroline'        : Tag('CL'  , -72.78  , 72.05     , 0.607),
    'ConwayReef'      : Tag('CR'  , -20.40  , 170.53    , 3.923),
    'Easter'          : Tag('EA'  , 24.97   , 67.53     , 11.334),
    'Futuna'          : Tag('FT'  , -16.33  , 178.07    , 5.101),
    'Galapagos'       : Tag('GP'  , 2.53    , 81.18     , 5.487),
    'JuanFernandez'   : Tag('JZ'  , 34.25   , 70.74     , 22.368),
    'Kermadec'        : Tag('KE'  , 39.99   , 6.46      , 2.347),
    'Manus'           : Tag('MN'  , -3.67   , 150.27    , 51.569),
    'Maoke'           : Tag('MO'  , 14.25   , 92.67     , 0.774),
    'Mariana'         : Tag('MA'  , 11.05   , 137.84    , 1.306),
    'MoluccaSea'      : Tag('MS'  , 2.15    , -56.09    , 3.566),
    'NewHebrides'     : Tag('NH'  , 0.57    , -6.60     , 2.469),
    'Niuafo`ou'       : Tag('NI'  , -3.29   , -174.49   , 3.314),
    'NorthAndes'      : Tag('ND'  , 17.73   , -122.68   , 0.116),
    'NorthBismarck'   : Tag('NB'  , -45.04  , 127.64    , 0.856),
    'Okhotsk'         : Tag('OK'  , 30.30   , -92.28    , 0.229),
    'Okinawa'         : Tag('ON'  , 36.12   , 137.92    , 2.539),
    'Panama'          : Tag('PM'  , 31.35   , -113.90   , 0.317),
    'Shetland'        : Tag('SL'  , 50.71   , -143.47   , 0.268),
    'SolomonSea'      : Tag('SS'  , -2.87   , 130.62    , 1.703),
    'SouthBismarck'   : Tag('SB'  , 6.88    , -31.89    , 8.111),
    'Timor'           : Tag('TI'  , -4.44   , 113.50    , 1.864),
    'Tonga'           : Tag('TO'  , 25.87   , 4.48      , 8.942),
    'Woodlark'        : Tag('WL'  , 0.10    , 128.52    , 1.744),
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
        # to avoid artifacts due to numerical precision
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
# Later will be moved to a separate script `plot_utils.py` in plate motion package

def read_plate_outline(pmm_name='GSRM', plate_name=None):
    """Read the plate boundaries for the given plate motion model.

    Parameters: pmm_name   - str, plate motion (model) name
                plate_name - str, plate name of interest, return all plates if None
    Returns:    outline    - dict, a dictionary that contains lists of vertices in lat/lon for all plates
                             OR shapely.geometry.polygon.Polygon object, boundary of the given "plate".
    """

    # check input
    if 'GSRM' in pmm_name:
        pmm_name = 'GSRM'
        pmm_dict = GSRM_V21_PMM

    elif 'MORVEL' in pmm_name:
        pmm_name = 'MORVEL'
        pmm_dict = NNR_MORVEL56_PMM

    else:
        msg = f'Un-recognized plate motion model: {pmm_name}!'
        msg += '\nAvailable models: GSRM, MORVEL.'
        raise ValueError(msg)

    # plate boundary file
    plate_boundary_file = PLATE_BOUNDARY_FILE[pmm_name]
    coord_order = os.path.basename(plate_boundary_file).split('.')[-1]
    if coord_order not in ['lalo', 'lola']:
        raise ValueError(f'Can NOT recognize the lat/lon order from the file extension: .{coord_order}!')

    # dict to convert plate abbreviation to name
    plate_abbrev2name = {}
    for key, val in pmm_dict.items():
        plate_abbrev2name[val.Abbrev.upper()] = key

    # read the plate outlines file, save them to a dictionary {plate_A: [vertices], ..., ..., ...}
    outlines = {}
    with open(plate_boundary_file) as f:
        lines = f.readlines()
        key, vertices = None, None
        # loop over lines to read
        for line in lines:
            # whether we meet a new plate name abbreviation
            if line.startswith('> ') or line.startswith('# ') or len(line.split()) == 1:
                # whether to add the previous plate to the dictionary
                if key and vertices:
                    pname = plate_abbrev2name[key]
                    outlines[pname] = np.array(vertices)
                # identify the new plate name abbreviation
                if line.startswith('> '):
                    key = line.split('> ')[1]
                elif line.startswith('# '):
                    key = line.split('# ')[1]
                else:
                    key = str(line)
                # remove the line change string
                if key.endswith('\n'):
                    key = key.split('\n')[0]
                # new vertices for the new plate
                vertices = []

            # get plate outline vertices
            else:
                vert = np.array(line.split()).astype(float)
                if coord_order == 'lola':
                    vert = np.flip(vert)
                vertices.append(vert)

    # outline of a specific plate
    if plate_name:
        if plate_name not in plate_abbrev2name.values():
            plate_abbrev = pmm_dict[plate_name].Abbrev
            raise ValueError(f'Can NOT found plate {plate_name} ({plate_abbrev}) in file: {plate_boundary_file}!')

        # convert list into shapely polygon object
        # for easy use
        outline = geometry.Polygon(outlines[plate_name])

    else:
        outline = outlines

    return outline


def plot_plate_motion(plate_boundary, epole_obj, center_lalo=None, qscale=200, qunit=50,
                      satellite_height=1e6, figsize=[5, 5], **kwargs):
    """Plot the globe map wityh plate boundary, quivers on some points.

    Parameters: plate_boundary   - shapely.geometry.Polygon object
                epole_obj        - mintpy.objects.euler_pole.EulerPole object
                center_lalo      - list of 2 float, center the map at this latitude, longitude
                qscale           - float, scaling factor of the quiver
                qunit            - float, length of the quiver legend in mm/yr
                satellite_height - height of the perspective view looking in meters
                kwargs           - dict, dictionary for plotting
    Returns:    fig, ax          - matplotlib figure and axes objects
    Examples:
        from matplotlib import pyplot as plt
        from mintpy.objects import euler_pole
        from shapely import geometry

        # build EulerPole object
        plate_pmm = euler_pole.ITRF2014_PMM['Arabia']
        epole_obj = euler_pole.EulerPole(wx=plate_pmm.omega_x, wy=plate_pmm.omega_y, wz=plate_pmm.omega_z)

        # read plate boundary
        plate_boundary = euler_pole.read_plate_outline('GSRM', 'Arabia')

        # plot plate motion
        fig, ax = euler_pole.plot_plate_motion(plate_boundary, epole_obj)
        plt.show()
    """

    def _sample_coords_within_polygon(polygon_obj, ny=10, nx=10):
        """Make a set of points inside the defined sphericalpolygon object.

        Parameters: polygon_obj - shapely.geometry.Polygon, a polygon object in lat/lon.
                    ny          - int, number of initial sample points in the y (lat) direction.
                    nx          - int, number of initial sample points in the x (lon) direction.
        Returns:    sample_lats - 1D np.ndarray, sample coordinates   in the y (lat) direction.
                    sample_lons - 1D np.ndarray, sample coordinates   in the x (lon) direction.
        """
        # generate sample point grid
        poly_lats = np.array(polygon_obj.exterior.coords)[:,0]
        poly_lons = np.array(polygon_obj.exterior.coords)[:,1]
        cand_lats, cand_lons = np.meshgrid(
            np.linspace(np.min(poly_lats), np.max(poly_lats), ny),
            np.linspace(np.min(poly_lons), np.max(poly_lons), nx),
        )
        cand_lats = cand_lats.flatten()
        cand_lons = cand_lons.flatten()

        # select points inside the polygon
        flag = np.zeros(cand_lats.size, dtype=np.bool_)
        for i, (cand_lat, cand_lon) in enumerate(zip(cand_lats, cand_lons)):
            if polygon_obj.contains(geometry.Point(cand_lat, cand_lon)):
                flag[i] = True
        sample_lats = cand_lats[flag]
        sample_lons = cand_lons[flag]

        return sample_lats, sample_lons

    # default plot settings
    kwargs['c_ocean']     = kwargs.get('c_ocean', 'w')
    kwargs['c_land']      = kwargs.get('c_land', 'lightgray')
    kwargs['c_plate']     = kwargs.get('c_plate', 'mistyrose')
    kwargs['lw_coast']    = kwargs.get('lw_coast', 0.5)
    kwargs['lw_pbond']    = kwargs.get('lw_pbond', 1)
    kwargs['lc_pbond']    = kwargs.get('lc_pbond', 'coral')
    kwargs['alpha_plate'] = kwargs.get('alpha_plate', 0.4)
    kwargs['grid_ls']     = kwargs.get('grid_ls', '--')
    kwargs['grid_lw']     = kwargs.get('grid_lw', 0.3)
    kwargs['grid_lc']     = kwargs.get('grid_lc', 'gray')
    kwargs['qnum']        = kwargs.get('qnum', 6)
    # point of interest
    kwargs['pts_lalo']    = kwargs.get('pts_lalo', None)
    kwargs['pts_marker']  = kwargs.get('pts_marker', '^')
    kwargs['pts_ms']      = kwargs.get('pts_ms', 20)
    kwargs['pts_mfc']     = kwargs.get('pts_mfc', 'r')
    kwargs['pts_mec']     = kwargs.get('pts_mec', 'k')
    kwargs['pts_mew']     = kwargs.get('pts_mew', 1)

    # map projection
    # based on: 1) map center and 2) satellite_height
    if not center_lalo:
        if kwargs['pts_lalo']:
            center_lalo = kwargs['pts_lalo']
        else:
            center_lalo = np.array(plate_boundary.centroid.coords)[0]
    map_proj = ccrs.NearsidePerspective(center_lalo[1], center_lalo[0], satellite_height=satellite_height)

    # make a base map from cartopy
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection=map_proj))
    ax.set_global()
    ax.gridlines(color=kwargs['grid_lc'],
                 linestyle=kwargs['grid_ls'],
                 linewidth=kwargs['grid_lw'],
                 xlocs=np.arange(-180,180,30),
                 ylocs=np.linspace(-80,80,10))
    ax.add_feature(cfeature.OCEAN, color=kwargs['c_ocean'])
    ax.add_feature(cfeature.LAND,  color=kwargs['c_land'])
    ax.add_feature(cfeature.COASTLINE, linewidth=kwargs['lw_coast'])

    # add the plate polygon
    if plate_boundary:
        poly_lats = np.array(plate_boundary.exterior.coords)[:, 0]
        poly_lons = np.array(plate_boundary.exterior.coords)[:, 1]
        ax.plot(poly_lons, poly_lats, color=kwargs['lc_pbond'], transform=ccrs.Geodetic(), linewidth=kwargs['lw_pbond'])
        ax.fill(poly_lons, poly_lats, color=kwargs['c_plate'],  transform=ccrs.Geodetic(), alpha=kwargs['alpha_plate'])

        # compute the plate motion from Euler rotation
        if epole_obj:
            # select sample points inside the polygon
            sample_lats, sample_lons = _sample_coords_within_polygon(plate_boundary, ny=kwargs['qnum'], nx=kwargs['qnum'])

            # calculate plate motion on sample points
            ve, vn = epole_obj.get_velocity_enu(lat=sample_lats, lon=sample_lons)[:2]

            # scale from m/yr to mm/yr
            ve *= 1e3
            vn *= 1e3
            norm = np.sqrt(ve**2 + vn**2)

            # correcting for "East" further toward polar region; re-normalize ve, vn
            ve /= np.cos(np.deg2rad(sample_lats))
            renorm = np.sqrt(ve**2 + vn**2)/norm
            ve /= renorm
            vn /= renorm

            # ---------- plot inplate vectors --------------
            q = ax.quiver(sample_lons, sample_lats, ve, vn,
                          transform=ccrs.PlateCarree(), scale=qscale,
                          width=.0075, color='coral', angles="xy")
            # legend
            # put an empty title for extra whitepace at the top
            ax.set_title('  ', pad=10)
            ax.quiverkey(q, X=0.3, Y=0.9, U=qunit, label=f'{qunit} mm/yr', labelpos='E', coordinates='figure')

    # add custom points (e.g., show some points of interest)
    if kwargs['pts_lalo']:
        ax.scatter(kwargs['pts_lalo'][1], kwargs['pts_lalo'][0],
                   marker=kwargs['pts_marker'], s=kwargs['pts_ms'],
                   fc=kwargs['pts_mfc'], ec=kwargs['pts_mec'],
                   lw=kwargs['pts_mew'], transform=ccrs.PlateCarree())

    return fig, ax
