#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Yuan-Kai Liu, Aug 2022        #
############################################################


import collections
import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

# ITRF2014-PMM defined in Altamimi et al. (2017)
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

# GSRMv2.1 defined in Kreemer et al. (2014)
# (unit: Lat: °N; Lon: °E; omega: °/Ma)
Tag = collections.namedtuple('Tag', 'Abbrev Lat Lon omega')
GSRM_PMM = {
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


# MORVEL56 defined in DeMets et al. (2010)
# (unit: Lat: °N; Lon: °E; omega: °/Ma)
Tag = collections.namedtuple('Tag', 'Abbrev Lat Lon omega')
MORVEL56_PMM = {
    'Amur'            : Tag('am'  , 63.17   , -122.82   , 0.297),
    'Antarctica'      : Tag('an'  , 65.42   , -118.11   , 0.250),
    'Arabia'          : Tag('ar'  , 48.88   , -8.49     , 0.559),
    'Australia'       : Tag('au'  , 33.86   , 37.94     , 0.632),
    'Capricorn'       : Tag('cp'  , 44.44   , 23.09     , 0.608),
    'Caribbean'       : Tag('ca'  , 35.20   , -92.62    , 0.286),
    'Cocos'           : Tag('co'  , 26.93   , -124.31   , 1.198),
    'Eurasia'         : Tag('eu'  , 48.85   , -106.50   , 0.223),
    'India'           : Tag('IN'  , 50.37   , -3.29     , 0.544),
    'JuandeFuca'      : Tag('jf'  , -38.31  , 60.04     , 0.951),
    'Lwandle'         : Tag('lw'  , 51.89   , -69.52    , 0.286),
    'Macquarie'       : Tag('mq'  , 49.19   , 11.05     , 1.144),
    'Nazca'           : Tag('nz'  , 46.23   , -101.06   , 0.696),
    'NorthAmerica'    : Tag('na'  , -4.85   , -80.64    , 0.209),
    'Nubia'           : Tag('nb'  , 47.68   , -68.44    , 0.292),
    'Pacific'         : Tag('pa'  , -63.58  , 114.70    , 0.651),
    'PhilippineSea'   : Tag('ps'  , -46.02  , -31.36    , 0.910),
    'Rivera'          : Tag('ri'  , 20.25   , -107.29   , 4.536),
    'Sandwich'        : Tag('SW'  , -29.94  , -36.87    , 1.362),
    'Scotia'          : Tag('SC'  , 22.52   , -106.15   , 0.146),
    'Somalia'         : Tag('sm'  , 49.95   , -84.52    , 0.339),
    'SouthAmerica'    : Tag('sa'  , -22.62  , -112.83   , 0.109),
    'Sunda'           : Tag('su'  , 50.06   , -95.02    , 0.337),
    'Sur'             : Tag('sr'  , -32.50  , -111.32   , 0.107),
    'Yangtze'         : Tag('yz'  , 63.03   , -116.62   , 0.334),
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

#########################################  Usage  ##############################################
REFERENCE = """reference:
  Stephenson, O. L., Liu, Y. K., Yunjun, Z., Simons, M., Rosen, P. and Xu, X., (2022),
    The Impact of Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields,
    Geophys. Res. Lett. 49, e2022GL099835, doi:10.1029/2022GL099835.

  # list of no-net-rotation (NNR) plate motion models (PMMs):
  # ONLY ITRF14 should be used, as Sentinel-1's orbit is in ITRF2014 reference frame.
  # Other values, e.g. MORVEL56, should be converted into ITR2014 before use.
  ITRF14 - Table 1 of Altamimi et al. (2017) - 11 plates
    Altamimi, Z., Métivier, L., Rebischung, P., Rouby, H., & Collilieux, X. (2017).
    ITRF2014 plate motion model. Geophysical Journal International, 209(3), 1906-1912.
    doi:10.1093/gji/ggx136
  MORVEL - Table 1 of Argus et al. (2011) - 56 plates
    Argus, D. F., Gordon, R. G., & DeMets, C. (2011). Geologically current motion of 56
    plates relative to the no-net-rotation reference frame. Geochemistry, Geophysics,
    Geosystems, 12(11). doi:10.1029/2011GC003751
"""

EXAMPLE = """example:
  # Use build-in plate motion model of Table 1 from Altamimi et al. (2017)
  plate_motion.py -g inputs/geometryGeo.h5   --plate Arabia
  plate_motion.py -g inputs/geometryRadar.h5 --plate Eurasia

  # Cartesian form of Euler pole rotation in [wx, wy, wz] in unit of mas/year [milli arc second per year]
  # e.g., Arabia plate in ITRF14-PMM (Table 1 in Altamimi et al., 2017)
  plate_motion.py -g inputs/geometryRadar.h5 --om-cart 1.154 -0.136  1.444
  plate_motion.py -g inputs/geometryGeo.h5   --om-cart 1.154 -0.136  1.444 -v velocity.h5
  plate_motion.py -g inputs/geometryRadar.h5 --om-cart 1.154 -0.136  1.444 --comp en2az

  # Simple constant local ENU translation (based on one GNSS vector) in [ve, vn, vu] in unit of m/year
  #   E.g., https://www.unavco.org/software/visualization/GPS-Velocity-Viewer/GPS-Velocity-Viewer.html
  #   -> select 'GNSS Data source' as 'World, IGS08/NNR, GEM GSRM' (referenced to ITRF2008, NNR PMM)
  #   -> check box `Station labels and data download` and click `Draw Map`
  #   -> navigate to the region of interest,
  #   -> click on a representative station,
  #   -> get the "Speed components" in mm/yr.
  plate_motion.py -g inputs/geometryGeo.h5 --enu 25.0 30.5 0.0 -v velocity.h5
"""

NOTE = """
  Removing the effect of rigid plate motion (translation and rotation) using plate motion model (PMM).
  For Sentinel-1, its orbit is measured with respect to ITRF2014 (Table 3-2 of Peter et al., 2021,
  Copernicus POD Product Handbook), which is an Earth-centered, Earth-fixed (ECEF) reference frame
  in which there is no net rotation (NNR) of the Earth surface.
"""

def create_parser(subparsers=None):
    synopsis = 'Plate Motion Correction.'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    # input files
    parser.add_argument('-g', '--geom', dest='geom_file', type=str, required=True,
                        help='Input geometry file in geo-coordinates, e.g., geometryGeo.h5')
    parser.add_argument('-v', '--velo', dest='vel_file', type=str,
                        help='Input velocity file to be corrected.')
    parser.add_argument('-o', '--output', dest='cor_vel_file', type=str,
                        help='Output velocity file after the correction, default: add "_ITRF14" suffix.')

    # plate motion calculation
    parser.add_argument('--comp', dest='pmm_comp', choices={'enu2los', 'en2az'}, default='enu2los',
                        help='Convert the ENU components into the component of interest for radar (default: %(default)s).')
    parser.add_argument('--step','--pmm-step', dest='pmm_step', type=float, default=10.,
                        help='Ground step/resolution in km for computing PMM to ENU velocity (default: %(default)s).')

    # plate motion configurations
    pmm = parser.add_argument_group('plate motion model')
    pmg = pmm.add_mutually_exclusive_group(required=True)
    pmg.add_argument('--plate', dest='plate_name', type=str, choices=ITRF2014_PMM.keys(), default=None,
                     help='Tectonic plate in ITRF14 (Table 1 in Altamimi et al., 2017).')
    pmg.add_argument('--om-cart', dest='omega_cart', type=float, nargs=3, metavar=('WX', 'WY', 'WZ'), default=None,
                     help='Cartesian form of Euler Pole rotation (unit: mas/yr) (default: %(default)s).')
    pmg.add_argument('--om-sph', dest='omega_sph', type=float, nargs=3, metavar=('LAT', 'LON', 'W'), default=None,
                     help='Spherical form of Euler Pole rotation (unit: deg, deg, deg/Ma) (default: %(default)s).')
    pmg.add_argument('--enu', dest='const_vel_enu', type=float, nargs=3, metavar=('VE', 'VN', 'VU'), default=None,
                     help='Constant local ground translation (unit: m/year) (default: %(default)s).')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default: output PMM filenames
    geom_dir = os.path.dirname(inps.geom_file)
    inps.pmm_enu_file = os.path.join(geom_dir, 'ITRF14enu.h5')
    inps.pmm_file = os.path.join(geom_dir, 'ITRF14.h5')
    if inps.pmm_comp.endswith('2az'):
        inps.pmm_file = os.path.join(geom_dir, 'ITRF14az.h5')

    # default: --output option
    if inps.vel_file and not inps.cor_vel_file:
        vbase = os.path.splitext(inps.vel_file)[0]
        inps.cor_vel_file = os.path.abspath(f'{vbase}_ITRF14.h5')

    # check: --plate option (and convert to --om-cart)
    if inps.plate_name:
        plate = ITRF2014_PMM[inps.plate_name]
        inps.omega_cart = [plate.omega_x, plate.omega_y, plate.omega_z]
        msg = f'get rotation parameters for {inps.plate_name} plate from Table 1 in Altamimi et al. (2017): '
        msg += f'wx, wy, wz = {plate.omega_x}, {plate.omega_y}, {plate.omega_z} mas/yr'
        print(msg)

    return inps


#######################################  Main Function  ########################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.plate_motion import run_plate_motion

    # run
    run_plate_motion(inps)


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
