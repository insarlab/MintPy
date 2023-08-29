#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Yuan-Kai Liu, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

ITRF2014_PMM_PLATES = [
    'Antartica',  'Arabia',       'Australia',    'Eurasia',
    'India',      'Nazca',        'NorthAmerica', 'Nubia',
    'Pacific',    'SouthAmerica', 'Somalia'
]

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
  # Use built-in plate motion model of Table 1 from Altamimi et al. (2017)
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
    pmg.add_argument('--plate', dest='plate_name', type=str, choices=ITRF2014_PMM_PLATES, default=None,
                     help='Tectonic plate in ITRF2014-PMM (Table 1 in Altamimi et al., 2017) (default: %(default)s).')
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
