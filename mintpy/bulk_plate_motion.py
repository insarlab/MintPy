#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Zhang Yunjun, May 2022             #
############################################################
#
# Extra dependency:
#   + platemotion: https://github.com/lcx366/PlateTectonic
#   + astropy
#   How to install both:
#      option (1) pip install platemotion
#      option (2) git clone git@github.com:lcx366/PlateTectonic.git $TOOL_DIR/PlateTectonic
#                 echo 'export PYTHONPATH=$PYTHONPATH:$TOOL_DIR/PlateTectonic' >> ~/.bashrc
#                 somehow install other dependencies in setup.py using your conda
#
# To-Do List (updated 2022.5.30 Yuan-Kai Liu):
#   + Potentially, we can make built-in PMM tables/dictionaries for easier user string input of the plate name
#   + Calculate Euler rotation to multi-points ENU motion is slow (called by pmm2enu_at() here)
#       In `platemotion` package, use array operation rather than for loops
#           https://github.com/lcx366/PlateTectonic/blob/main/platemotion/classes/plate.py#L153
#   + Replace platemotion package by equations of Euler trasnformation to relax this dependency at all?


import os
import sys
import collections
import numpy as np
from skimage.transform import resize

from mintpy.objects.resample import resample
from mintpy.utils import readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser
from mintpy.diff import diff_file

# https://docs.astropy.org/en/stable/units/index.html
try:
    from platemotion import Plate
    from astropy import units as u
except ImportError:
    msg = 'Can NOT import platemotion!'
    msg += '\nCheck more details at https://github.com/lcx366/PlateTectonic.'
    raise ImportError(msg)


# ITRF2014-PMM defined in Altamimi et al. (2017)
# units:
#   omega_x/y/z in mas/yr (milli-second of arc per year)
#   omega       in deg/Ma (degree per megayear or one-million-year)
#   wrms_e/n    in mm/yr  (milli-meter per year), WRMS: weighted root mean scatter
Tag = collections.namedtuple('Tag', 'name num_site omega_x omega_y omega_z omega wrms_e wrms_n')
ITRF2014_PMM = {
    'ANTA' : Tag('Antartica'  ,   7,  -0.248,  -0.324,   0.675,  0.219,  0.20,  0.16),
    'ARAB' : Tag('Arabia'     ,   5,   1.154,  -0.136,   1.444,  0.515,  0.36,  0.43),
    'AUST' : Tag('Australia'  ,  36,   1.510,   1.182,   1.215,  0.631,  0.24,  0.20),
    'EURA' : Tag('Eurasia'    ,  97,  -0.085,  -0.531,   0.770,  0.261,  0.23,  0.19),
    'INDI' : Tag('India'      ,   3,   1.154,  -0.005,   1.454,  0.516,  0.21,  0.21),
    'NAZC' : Tag('Nazca'      ,   2,  -0.333,  -1.544,   1.623,  0.629,  0.13,  0.19),
    'NOAM' : Tag('N. America' ,  72,   0.024,  -0.694,  -0.063,  0.194,  0.23,  0.28),
    'NUBI' : Tag('Nubia'      ,  24,   0.099,  -0.614,   0.733,  0.267,  0.28,  0.36),
    'PCFC' : Tag('Pacific'    ,  18,  -0.409,   1.047,  -2.169,  0.679,  0.36,  0.31),
    'SOAM' : Tag('S. America' ,  30,  -0.270,  -0.301,  -0.140,  0.119,  0.34,  0.35),
    'SOMA' : Tag('Somalia'    ,   3,  -0.121,  -0.794,   0.884,  0.332,  0.32,  0.30),
}



#########################################  Usage  ##############################################
REFERENCE = """reference:
  Stephenson, O. L., Liu, Y. K., Yunjun, Z., Simons, M., Rosen, P. and Xu, X., (2022), The Impact of
    Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields, Geophys. Res. Lett. (under review)
    doi:10.1002/essoar.10511538.2
  Peter, H., Fernández, M., Aguilar, J., & Fernández, J. (2021). Copernicus POD Product Handbook:
    Copernicus Sentinel-1, -2 and -3 Precise orbit Determination Serivice (CPOD) (GMV-CPOD-TN-0009).
    https://sentinels.copernicus.eu/documents/247904/3372484/Sentinels-POD-Product-Handbook-1.19.pdf

  # list of no-net-rotation (NNR) plate motion models (PMMs):
  # ONLY ITRF14 should be used, as Sentinel-1's orbit is in ITRF2014 reference frame.
  # Other values, e.g. MORVEL56, should be converted into ITR2014 before use.
  ITRF14 - Table 1 of Altamimi et al. (2017) - 11 plates
    Altamimi, Z., Métivier, L., Rebischung, P., Rouby, H., & Collilieux, X. (2017).
    ITRF2014 plate motion model. Geophysical Journal International, 209(3), 1906-1912.
    doi:10.1093/gji/ggx136
  MORVEL56 - Table 1 of Argus et al. (2011) - 56 plates
    Argus, D. F., Gordon, R. G., & DeMets, C. (2011). Geologically current motion of 56
    plates relative to the no-net-rotation reference frame. Geochemistry, Geophysics, Geosystems, 12(11).
    doi:10.1029/2011GC003751
"""

EXAMPLE = """example:
  # Cartesian form of Euler pole rotation in [wx, wy, wz] in unit of mas/year [milli arc second per year]
  # e.g., Arabia plate in ITRF14-PMM (Table 1 in Altamimi et al., 2017)
  bulk_plate_motion.py -g inputs/geometryGeo.h5   --om-cart 1.154 -0.136  1.444 -v velocity.h5
  bulk_plate_motion.py -g inputs/geometryRadar.h5 --om-cart 1.154 -0.136  1.444

  # Simple constant local ENU translation (based on one GNSS vector) in [ve, vn, vu] in unit of m/year
  #   E.g., https://www.unavco.org/software/visualization/GPS-Velocity-Viewer/GPS-Velocity-Viewer.html
  #   -> select 'GNSS Data source' as 'World, IGS08/NNR, GEM GSRM' (referenced to ITRF2008, NNR PMM)
  #   -> check box `Station labels and data download` and click `Draw Map`
  #   -> navigate to the region of interest,
  #   -> click on a representative station,
  #   -> get the "Speed components" in mm/yr.
  bulk_plate_motion.py -g inputs/geometryGeo.h5 --enu 25.0 30.5 0.0 -v velocity.h5
"""

NOTE = """
  Removing the effect of bulk traslation and rotation based on a given plate motion model (PMM).
  For Sentinel-1, its orbit is measured with respect to ITRF2014 (Table 3-2 of Peter et al., 2021), which is an
  Earth-centered, Earth-fixed (ECEF) reference frame in which there is no net rotation (NNR) of the Earth surface.
"""

def create_parser(subparsers=None):
    synopsis = 'Bulk Plate Motion Correction.'
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

    # plate motion configurations
    pmm = parser.add_mutually_exclusive_group(required=True)
    pmm.add_argument('--om-cart', dest='omega_cart', type=float, nargs=3, metavar=('WX', 'WY', 'WZ'), default=None,
                     help='Cartesian form of Euler Pole rotation (unit: mas/yr) (default: %(default)s).')
    pmm.add_argument('--om-sph', dest='omega_sph', type=float, nargs=3, metavar=('LAT', 'LON', 'W'), default=None,
                     help='Spherical form of Euler Pole rotation (unit: deg, deg, deg/Ma) (default: %(default)s).')
    pmm.add_argument('--enu', dest='const_vel_enu', type=float, nargs=3, metavar=('VE', 'VN', 'VU'), default=None,
                     help='Constant local ground translation (unit: m/year) (default: %(default)s).')

    parser.add_argument('--step','--pmm-step', dest='pmm_step', type=float, default=10.,
                        help='Ground step/resolution in km for computing PMM to ENU velocity (default: %(default)s).')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default output filenames
    geom_dir = os.path.dirname(inps.geom_file)
    inps.pmm_enu_file = os.path.join(geom_dir, 'ITRF14ENU.h5')
    inps.pmm_los_file = os.path.join(geom_dir, 'ITRF14.h5')

    if inps.vel_file and not inps.cor_vel_file:
        vbase = os.path.splitext(inps.vel_file)[0]
        inps.cor_vel_file = os.path.abspath(f'{vbase}_ITRF14.h5')

    return inps


########################################## Sub Functions #############################################

def build_plate_motion_model(omega_cart=None, omega_sph=None):
    """Build a plate motion model based on the given Euler roation vector
    Parameters: omega_sph  - list or np.array, Spherical representation [lat, lon, w] (deg, deg, deg/Ma)
                omega_cart - list or np.array, Cartesian representation [wx, wy, wz] (mas/yr)
    Returns:    plate      - platemotion.Plate object
    """
    # Check input variables
    if (omega_cart is None) and (omega_sph is None):
        raise ValueError('Neither omega_cart (wxyz) nor omega_sph (Euler Pole) are given! At least one is required.')

    # Set a NaN moment of inertia (to get class Plate running)
    iner_null = np.zeros([3,3]) * np.nan * (u.kg * u.km**2)

    # Create an instrance of class Plate from `platemotion` pkg
    plate = Plate(info={'inertia_tensor': iner_null})
    if omega_cart is not None:
        print('input omega_cartesian in [wx, wy, wz] (mas/yr)')
        omega = np.array(omega_cart) * u.mas/u.yr
        plate.set_omega(omega, 'cartesian')
    else:
        print('input omega_spherical in [lat, lon, w] (deg, deg, deg/Ma)')
        omega = [omega_sph[0]*u.deg,
                 omega_sph[1]*u.deg,
                 omega_sph[2]*u.deg/u.Ma]
        plate.set_omega(omega, 'spherical')

    print('\n--------------------------------------')
    print('Euler Pole and Rotation Vector')
    print('in spherical coordinates:')
    print('  Pole Latitude : {:10.4f} deg'.format(plate.omega_spherical[0].degree))
    print('  Pole Longitude: {:10.4f} deg'.format(plate.omega_spherical[1].degree))
    print('  Rotation rate : {:10.4f}  \n'.format(plate.omega_spherical[2].to(u.deg/u.Ma)))
    print('in Cartesian coordinates:')
    print('  wx: {:10.4f}'.format(plate.omega_cartesian[0]))
    print('  wy: {:10.4f}'.format(plate.omega_cartesian[1]))
    print('  wz: {:10.4f}'.format(plate.omega_cartesian[2]))
    print('--------------------------------------\n')
    return plate


def pmm2enu_at(pmm_obj, lats, lons):
    """Evaluate the PMM at given lats/lons for the motion in ENU.

    Parameters: pmm_obj - plate motion model instance
                lats    - 0/1/2D array in float32, latitudes
                lons    - 0/1/2D array in float32, longitudes
    Returns:    ve/n/u  - 0/1/2D array in float32, plate motion in east / north / up
                          in meter/year.
    """
    if isinstance(lats, float) or isinstance(lats, int):
        loc = np.array([lats, lons, 0])
        v = pmm_obj.velocity_at(loc,'geodetic')
        ve = v.en[0]
        vn = v.en[1]
        vu = 0

    elif lats.ndim in [1, 2]:
        # prepare locations as array in size of (3, num_pts)
        elev = np.zeros_like(lats)
        locs = np.vstack((
            lats.flatten(),
            lons.flatten(),
            elev.flatten(),
        ))
        # run PMM
        v = pmm_obj.velocity_at(locs, 'geodetic')
        ve = v.en[:, 0].reshape(lats.shape)
        vn = v.en[:, 1].reshape(lats.shape)
        vu = np.zeros(lats.shape, dtype=np.float32)

    else:
        raise ValueError(f'Un-recognized lat/lon grid dimension: {lats.ndim}!')

    # convert from mm/year to meter/year
    #     and from astropy.units.quantity.Quantity to np.ndarray
    ve = np.array(ve, dtype=np.float32) * 1e-3
    vn = np.array(vn, dtype=np.float32) * 1e-3
    vu = np.array(vu, dtype=np.float32) * 1e-3

    return ve, vn, vu


####################################### Higher-level Sub Functions ##########################################

def calc_bulk_plate_motion(geom_file, omega_cart=None, omega_sph=None, const_vel_enu=None,
                           pmm_enu_file=None, pmm_los_file=None, pmm_step=10.):
    """Estimate LOS motion due to pure bulk tranlation or due to plate rotation
    Parameters: geom_file     - str, path to the input geometry file
                omega_cart    - list or 1D array, Cartesian representation of plate rotation
                                in [wx, wy, wz]  (mas/yr)
                omega_sph     - list or 1D array, Spherical representation of plate rotation
                                in [lat, lon, w] (deg, deg, deg/Ma)
                const_vel_enu - list or 1D array, a single-vector [ve, vn, vu] (meter/year)
                                simulating the bulk translation of the ground (e.g., from GNSS)
                pmm_enu_file  - str, path to the output bulk plate motion in east, north, up direction
                pmm_los_file  - str, path to the output bulk plate motion in LOS direction
                pmm_reso      - float, ground resolution for computing Plate rotation to ENU velocity (km)
    Returns:    ve/vn/vu/vlos - 2D np.ndarray, bulk plate motion in east / north / up / LOS direction
    """

    # Get LOS geometry
    atr_geo = ut.prepare_geo_los_geometry(geom_file, unit='deg')[2]
    shape_geo = [int(atr_geo['LENGTH']), int(atr_geo['WIDTH'])]

    ## Bulk motion model in the region
    print('-'*50)
    if omega_cart or omega_sph:
        print('compute the bulk plate motion using a plate motion model (PMM; translation & rotation)')
        if omega_cart is not None:
            pmm_obj = build_plate_motion_model(omega_cart=omega_cart)
        else:
            pmm_obj = build_plate_motion_model(omega_sph=omega_sph)

        # prepare the coarse grid
        latc = float(atr_geo['Y_FIRST']) + float(atr_geo['Y_STEP']) * shape_geo[0] / 2
        ystep = abs(int(pmm_step * 1000 / (float(atr_geo['Y_STEP']) * 108e3)))
        xstep = abs(int(pmm_step * 1000 / (float(atr_geo['X_STEP']) * 108e3 * np.cos(np.deg2rad(latc)))))
        ystep, xstep = max(ystep, 5), max(xstep, 5)
        lats, lons = ut.get_lat_lon(atr_geo, dimension=2, ystep=ystep, xstep=xstep)

        # transform PMM to ENU velocity on a coarse grid
        # to-do: multi-pixel rotation is slow; change the `platemotion` code to Big Array operation rather than For Loops
        print(f'compute PMM via platemotion.Plate: grid_size = {pmm_step} km, grid_shape = {lats.shape} ...')
        ve_low, vn_low, vu_low = pmm2enu_at(pmm_obj, lats, lons)

        # interpolate back to the initial grid
        print(f'interpolate corase PMM to the full resolution: {lats.shape} -> {shape_geo}'
              ' via skimage.transform.resize ...')
        kwargs = dict(order=1, mode='edge', anti_aliasing=True, preserve_range=True)
        ve = resize(ve_low, shape_geo, **kwargs)
        vn = resize(vn_low, shape_geo, **kwargs)
        vu = resize(vu_low, shape_geo, **kwargs)

    elif const_vel_enu:
        print(f'compute the bulk plate motion using a single vector (translation): {const_vel_enu}')
        ve = const_vel_enu[0] * np.ones(shape_geo, dtype=np.float32)
        vn = const_vel_enu[1] * np.ones(shape_geo, dtype=np.float32)
        vu = const_vel_enu[2] * np.ones(shape_geo, dtype=np.float32)


    # radar-code PMM if input geometry is in radar coordinates
    atr = readfile.read_attribute(geom_file)
    if 'Y_FIRST' not in atr.keys():
        print('radar-coding the bulk plate motion in ENU ...')
        res_obj = resample(lut_file=geom_file)
        res_obj.open()
        res_obj.src_meta = atr_geo
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        ve = res_obj.run_resample(src_data=ve[box[1]:box[3], box[0]:box[2]])
        vn = res_obj.run_resample(src_data=vn[box[1]:box[3], box[0]:box[2]])
        vu = res_obj.run_resample(src_data=vu[box[1]:box[3], box[0]:box[2]])


    # Project model to LOS velocity
    print('project the bulk plate motion from ENU onto LOS direction')
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    az_angle = readfile.read(geom_file, datasetName='azimuthAngle')[0]
    vlos = ut.enu2los(ve, vn, vu, inc_angle=inc_angle, az_angle=az_angle)

    # Save the bulk motion model velocity into HDF5 files
    atr['FILE_TYPE'] = 'velocity'
    atr['DATA_TYPE'] = 'float32'
    atr['UNIT'] = 'm/year'
    dsDict = {
        'east'  : ve,
        'north' : vn,
        'up'    : vu,
    }
    writefile.write(dsDict, out_file=pmm_enu_file, metadata=atr)
    writefile.write(vlos,   out_file=pmm_los_file, metadata=atr)

    return ve, vn, vu, vlos


def correct_bulk_plate_motion(vel_file, mfile, ofile):
    """Apply the bulk motion correction from files.
    """
    file1 = vel_file       # input uncorrected LOS velocity file
    file2 = [mfile]        # PMM LOS velocity file
    diff_file(file1, file2, ofile)
    return


#######################################  Main Function  ########################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    calc_bulk_plate_motion(
        geom_file=inps.geom_file,
        omega_cart=inps.omega_cart,
        omega_sph=inps.omega_sph,
        const_vel_enu=inps.const_vel_enu,
        pmm_enu_file=inps.pmm_enu_file,
        pmm_los_file=inps.pmm_los_file,
        pmm_step=inps.pmm_step,
    )

    if inps.vel_file and inps.pmm_los_file and inps.cor_vel_file:
        print('-'*50)
        print('Correct input velocity for the bulk plate motion')
        correct_bulk_plate_motion(inps.vel_file, inps.pmm_los_file, inps.cor_vel_file)

    return

################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
