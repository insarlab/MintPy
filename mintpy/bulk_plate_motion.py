#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, May 2022                           #
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
#   + Replace scipy.interpolate with alternatives for efficiency. E.g.:
#       skimage.resize https://scikit-image.org/docs/stable/auto_examples/transform/plot_rescale.html
#       check mintpy usage https://github.com/insarlab/MintPy/blob/7adb3a11f875b832488a0c8e44c174d98b1df254/mintpy/tropo_gacos.py#L129
#   + Calculate Euler rotation to multi-points ENU motion is slow (called by pmm2enu_at() here)
#       In `platemotion` package, use array operation rather than for loops (at https://github.com/lcx366/PlateTectonic/blob/main/platemotion/classes/plate.py#L153)
#   + Replace platemotion package by equations of Euler trasnformation to relax this dependency at all?

import os
import sys
import scipy
import pyproj
import argparse
import numpy as np

from mintpy.utils import readfile, writefile, utils as ut
from mintpy.diff import diff_file
from mintpy.save_gmt import get_geo_lat_lon
from mintpy.solid_earth_tides import prepare_los_geometry

# https://docs.astropy.org/en/stable/units/index.html
try:
    from platemotion import Plate
    from astropy import units as u
except ImportError:
    msg = 'Can NOT import platemotion!'
    msg += '\nCheck more details at https://github.com/lcx366/PlateTectonic.'
    raise ImportError(msg)



#########################################  Usage  ##############################################
REFERENCE = """reference:
  Stephenson, O. L., Liu, Y. K., Yunjun, Z., Simons, M., Rosen, P. and Xu, X., (2022), The Impact of
    Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields, Geophys. Res. Lett. (under review)
    doi:10.1002/essoar.10511538.2
  Peter, H., Fernández, M., Aguilar, J., & Fernández, J. (2021). Copernicus POD Product Handbook: 
    Copernicus Sentinel-1, -2 and -3 Precise orbit Determination Serivice (CPOD) (GMV-CPOD-TN-0009). 
    https://sentinels.copernicus.eu/documents/247904/3372484/Sentinels-POD-Product-Handbook-1.19.pdf
"""

EXAMPLE = """example:
  # Spherical form of Euler Pole rotation in [lat, lon, w] in unit of deg, deg, deg/Ma
  #   Africa  plate (NNR-NUVEL1A)  - Table 2 in Argus & Gordon (1991, GRL), doi:10.1029/91GL01532
  #   Eurasia plate (NNR-MORVEL56) - Table 1 in Argus, Gordon & DeMets (2011, G3), doi:10.1029/2011GC003751
  bulk_plate_motion.py -g inputs/geometryGeo.h5 --om-sph  50.6  -74.0   0.30
  bulk_plate_motion.py -g inputs/geometryGeo.h5 --om-sph  48.85 -106.50 0.223 -v velocity.h5 

  # Cartesian form of Euler Pole rotation in [wx, wy, wz] in unit of mas/year [milli arc second per year]
  #   Arabia plate (NNR-ITRF14) - Table 1 in Altamimi et al. (2017, GJI), doi:10.1093/gji/ggx136
  bulk_plate_motion.py -g inputs/geometryGeo.h5 --om-cart 1.154 -0.136  1.444 -v velocity.h5

  # Simple constant local ENU translation (based on one GNSS vector) in [ve, vn, vu] in unit of meter/year
  #   E.g., https://www.unavco.org/software/visualization/GPS-Velocity-Viewer/GPS-Velocity-Viewer.html
  #   Step 1: select `GNSS Data source` as `World, IGS08/NNR, GEM GSRM` (referenced to ITRF2008, NNR PMM)
  #   Step 2: check box `Station labels and data download` and click `Draw Map`
  #   Step 3: Navigate to the region of interest, click on a representative station, get the "Speed components" in mm/yr.
  bulk_plate_motion.py -g inputs/geometryGeo.h5 --enu 25.0 30.5 0.0 -v velocity.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(
        description='Bulk Plate Motion Correction.\n'
                    '  Removing the effect of bulk traslation and rotation in velocity field based on a given plate motion model (PMM).\n'
                    '  E.g., Sentinel-1 orbit is measured with respect to ITRF2014 (Table 3-2 of Peter et al., 2021), which is an\n'
                    '  Earth-centered, Earth-fixed reference frame in which there is no net rotation (NNR) of the Earth surface.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=REFERENCE+'\n'+EXAMPLE,
    )

    # input files
    parser.add_argument('-g', '--geom', dest='geom_file', type=str, required=True,
                        help='Input geometry file in geo-coordinates, e.g., geometryGeo.h5')
    parser.add_argument('-v', '--velo', dest='vel_file', type=str,
                        help='Input velocity file to be corrected.')
    parser.add_argument('-o', '--output', dest='cor_vel_file', type=str,
                        help='Output velocity file after the correction, default: add "_BPM" suffix.')

    # plate motion configurations
    pmms = parser.add_mutually_exclusive_group(required=True)
    pmms.add_argument('--om-sph', dest='omega_sph', type=float, nargs=3, metavar=('LAT', 'LON', 'W'), default=None,
                      help='Spherical form of Euler Pole rotation; [lat, lon, w] (unit: deg, deg, deg/Ma) (default: %(default)s).')
    pmms.add_argument('--om-cart', dest='omega_cart', type=float, nargs=3, metavar=('WX', 'WY', 'WZ'), default=None,
                      help='Cartesian form of Euler Pole rotation; [wx, wy, wz] (unit: mas/yr) (default: %(default)s).')
    pmms.add_argument('--enu', dest='const_vel_enu', type=float, nargs=3, metavar=('VE', 'VN', 'VU'), default=None,
                      help='Simple constant local ENU translation of ground [ve, vn, vu] unit: meter/year (default: %(default)s).')

    parser.add_argument('--reso','--pmm-reso', dest='pmm_reso', type=float, default=10.,
                        help='Ground resolution for computing Plate rotation to ENU velocity (unit: km) (default: %(default)s).')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default output filenames
    inps.bpm_enu_file = os.path.join(os.path.dirname(inps.geom_file), 'BulkPlateMotion3D.h5')
    inps.bpm_los_file = os.path.join(os.path.dirname(inps.geom_file), 'BulkPlateMotion.h5')

    if inps.vel_file and not inps.cor_vel_file:
        vbase = os.path.splitext(vel_file)[0]
        inps.cor_vel_file  = os.path.abspath(f'{vbase}_BPM.h5')

    return inps


########################################## Sub Functions #############################################

def build_plate_motion_model(omega_sph=None, omega_cart=None):
    """Build a plate motion model based on the given Euler roation vector
    Parameters: omega_sph  - list or np.array, Spherical representation [lat, lon, w] (deg, deg, deg/Ma)
                omega_cart - list or np.array, Cartesian representation [wx, wy, wz] (mas/yr)
    Returns:    plate      - platemotion.Plate object
    """

    # Set a NaN moment of inertia (to get class Plate running)
    iner_null = np.zeros([3,3]) * np.nan * (u.kg * u.km**2)

    # Create an instrance of class Plate from `platemotion` pkg
    plate = Plate(info={'inertia_tensor': iner_null})

    # Check input variables
    if (omega_cart is None) and (omega_sph is None):
        raise ValueError('Neither omega_cart (wxyz) nor omega_sph (Euler Pole) are given! At least one is required.')

    elif omega_cart is not None:
        print('your input: omega_cartesian; [wx, wy, wz] (mas/yr)')
        omega = np.array(omega_cart) * u.mas/u.yr
        plate.set_omega(omega, 'cartesian')

    else:
        print('your input: omega_spherical; Euler pole; [lat, lon, w] (deg, deg, deg/Ma)')
        omega = [omega_sph[0]*u.deg,
                 omega_sph[1]*u.deg,
                 omega_sph[2]*u.deg/u.Ma]
        plate.set_omega(omega,'spherical')

    print('--- Euler Pole and rotation vector ---')
    print('\nSpherical representation:')
    print(' Latitude:      {:.4f} deg'.format(plate.omega_spherical[0].degree))
    print(' Longitude:     {:.4f} deg'.format(plate.omega_spherical[1].degree))
    print(' Rotation rate: {:.4f}  \n'.format(plate.omega_spherical[2].to(u.deg/u.Ma)))
    print('\nCartesian representation:')
    print(' wx:            {:.4f}'.format(plate.omega_cartesian[0]))
    print(' wy:            {:.4f}'.format(plate.omega_cartesian[1]))
    print(' wz:            {:.4f}'.format(plate.omega_cartesian[2]))
    print('--------------------------------------')
    return plate


def interp_2d3l(data, X, Y, nx, ny, kind):
    """Interpolate 3-layer 2D arrays individually
    can be very slow
    skimage.resize: https://scikit-image.org/docs/stable/auto_examples/transform/plot_rescale.html
    check mintpy: https://github.com/insarlab/MintPy/blob/7adb3a11f875b832488a0c8e44c174d98b1df254/mintpy/tropo_gacos.py#L129
    """
    # Check and flip Y array if needed
    if Y[0,0] > Y[-1,0]:
        Y = np.flipud(Y)

    f0 = scipy.interpolate.interp2d(X, Y, data[:,:,0], kind=kind)
    f1 = scipy.interpolate.interp2d(X, Y, data[:,:,1], kind=kind)
    f2 = scipy.interpolate.interp2d(X, Y, data[:,:,2], kind=kind)

    y_new = np.linspace(np.min(Y), np.max(Y), ny)
    x_new = np.linspace(np.min(X), np.max(X), nx)

    data_new = np.empty([len(y_new), len(x_new), 3])
    data_new[:,:,0] = f0(x_new, y_new)
    data_new[:,:,1] = f1(x_new, y_new)
    data_new[:,:,2] = f2(x_new, y_new)

    return data_new


def pmm2enu_at(pmm, Lats, Lons):
    """
    Input:
        pmm     plate motion model instance
        Lats    2D array of latitudes;                 dim = (length, width)
        Lons    2D array of longitudes;                dim = (length, width)

    Output:
        enu     3D array of {east, north, up} motions; dim = (length, width, 3)
    """
    if isinstance(Lats, float) or isinstance(Lats, int):
        print('Single location')
        loc = np.array([Lats, Lons, 0])
        v   = pmm.velocity_at(loc,'geodetic')   # default: mm/yr
        en  = np.array(v.en)
        enu = np.concatenate((en, [0]))

    elif len(Lats.shape) == 1:
        print('1D array locations: {}'.format(len(Lats)))
        ele  = np.zeros_like(Lats)
        locs = np.dstack((Lats, Lons, ele))[0].T
        v    = pmm.velocity_at(locs,'geodetic')   # default: mm/yr
        en   = np.array(v.en).reshape([-1,2])
        enu  = np.concatenate((en, np.zeros([en.shape[0],1])),1)

    elif len(Lats.shape) > 1:
        print('2D array locations: {}'.format(Lats.size))
        length, width = Lats.shape
        ele  = np.zeros_like(Lats)
        locs = np.dstack((Lats, Lons, ele))
        locs = locs.reshape([-1, 3]).T
        v    = pmm.velocity_at(locs,'geodetic')   # default: mm/yr
        en   = np.array(v.en).reshape([-1,2])
        enu  = np.concatenate((en, np.zeros([en.shape[0],1])),1)
        enu  = enu.reshape([length, width, -1])

    else:
        print('Weird lat lon grid input')

    enu *= 1e-3     # convert to meter/year
    return enu


def get_geobox_width_length(geo_box):
    """Get width and length of the geo box in km
    geo_box     tuple of 4 floats (lon_min, lat_max, lon_max, lat_min) in decimal degrees
    """
    geod = pyproj.Geod(ellps='WGS84')
    width_km  = geod.inv(geo_box[0], geo_box[3], geo_box[2], geo_box[3])[2] * 1e-3
    length_km = geod.inv(geo_box[0], geo_box[3], geo_box[0], geo_box[1])[2] * 1e-3
    return width_km, length_km


####################################### Higher-level Sub Functions ##########################################

def estimate_bulk_motion(geom_file, omega_sph=None, omega_cart=None, const_vel_enu=None,
                         bpm_enu_file=None, bpm_los_file=None, pmm_reso=5.):
    """Estimate LOS motion due to pure bulk tranlation or due to plate rotation
    Parameters: geom_file     - str, path to the input geometry file
                omega_sph     - list or 1D array, Spherical representation of plate rotation [lat, lon, w] (deg, deg, deg/Ma)
                omega_cart    - list or 1D array, Cartesian representation of plate rotation [wx, wy, wz]  (mas/yr)
                const_vel_enu - list or 1D array, a single-vector [ve, vn, vu] (meter/year) 
                                simulating the bulk translation of the ground (e.g., from GNSS)
                bpm_enu_file  - str, path to the output BPM (bulk plate motion) east, north, up velocity field
                bpm_los_file  - str, path to the output BPM (bulk plate motion) LOS velocity field
                pmm_reso      - float, ground resolution for computing Plate rotation to ENU velocity (km)
    Returns:    ve/vn/vu/vlos - 2D np.ndarray, bulk plate motion in east / north / up / line-of-sight direction
    """
    # Read attributes
    atr = readfile.read_attribute(geom_file)
    width, length = int(atr['WIDTH']), int(atr['LENGTH'])
    if 'Y_FIRST' not in atr.keys():
        raise ValueError('Input geometry file is NOT in geo-coordinates!')

    # Get LOS geometry
    inc_rad, head_rad, atr_geo = prepare_los_geometry(geom_file)
    lats, lons                 = get_geo_lat_lon(atr_geo)
    width_km, length_km        = get_geobox_width_length((np.min(lons), np.max(lats), np.max(lons), np.min(lats)))
    inc_deg, head_deg          = np.rad2deg(inc_rad), np.rad2deg(head_rad)

    # Turn default null value to nan
    inc_deg[inc_deg==0]    = np.nan
    inc_rad[inc_deg==0]    = np.nan
    head_deg[head_deg==90] = np.nan
    head_rad[head_deg==90] = np.nan

    ## Bulk motion model in the region
    print('\nForm a bulk motion model velocity field')
    if any([omega for omega in [omega_cart, omega_sph]]):
        print('estimate bulk motion by building a plate motion model (e.g., ITRF2014 MORVEL)')
        if omega_cart is not None:
            plate = build_plate_motion_model(omega_cart=omega_cart)
        elif omega_sph is None:
            plate = build_plate_motion_model(omega_sph=omega_sph)

        # Use the Euler rotation to compute the ENU for each pixel in the geometry
        skx, sky = int(width*pmm_reso/width_km), int(length*pmm_reso/length_km)
        Lats = np.meshgrid(lons, lats)[1][::sky, ::skx]
        Lons = np.meshgrid(lons, lats)[0][::sky, ::skx]
        print('Compute ENU from PMM with a ~{} km ground resolution'.format(pmm_reso))
        print('low resolution / original grid dimension: ({}, {}) / ({}, {})'.format(Lats.shape[0], Lats.shape[1], length, width))

        # transform PMM to ENU velocity on a coarse grid
        # to-do: multi-pixel rotation is slow; change the `platemotion` code to Big Array operation rather than For Loops
        v_enu_low = pmm2enu_at(plate, Lats, Lons)

        # interpolate back to the initial grid
        print('interpolate back to the original grid')
        v_enu = interp_2d3l(v_enu_low, Lons, Lats, width, length, 'linear')

    elif const_vel_enu:
        print('estimate bulk motion by assuming a single pure translation vector: {}'.format(const_vel_enu))
        v_enu = const_vel_enu * np.ones([length, width, 3])

    else:
        print('Error: Please specify either the Plate Rotation or a Single Translation Vector; See -h option')
        print('Need to give at least ONE of [--enu, --om_cart, --om_sph]')
        sys.exit(1)

    # Project model to LOS velocity
    print('Project the bulk motion ENU onto radar LOS velocity')
    ve = np.array(v_enu[:,:,0], dtype=np.float32)
    vn = np.array(v_enu[:,:,1], dtype=np.float32)
    vu = np.array(v_enu[:,:,2], dtype=np.float32)
    vlos = ut.enu2los(ve, vn, vu, inc_angle=inc_deg, head_angle=head_deg)

    # Save the bulk motion model velocity into HDF5 files
    atr['FILE_TYPE'] = 'velocity'
    atr['DATA_TYPE'] = 'float32'
    atr['UNIT'] = 'm/year'
    dsDict = {
        'east'  : ve,
        'north' : vn,
        'up'    : vu,
    }
    writefile.write(dsDict, out_file=bpm_enu_file, metadata=atr)
    writefile.write(vlos,   out_file=bpm_los_file, metadata=atr)

    return ve, vn, vu, vlos


def remove_bulk_motion(vel_file, mfile, ofile):
    """
    Apply the bulk motion correction from files
    """
    file1 = vel_file       # input uncorrected LOS velocity file
    file2 = [mfile]        # BPM model LOS velocity file
    diff_file(file1, file2, ofile)
    return


#######################################  Main Function  ########################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n-----------------------------------')
    print('Evaluate bulk motion')
    print('-----------------------------------\n')
    estimate_bulk_motion(
        geom_file=inps.geom_file,
        omega_sph=inps.omega_sph,
        omega_cart=inps.omega_cart,
        const_vel_enu=inps.const_vel_enu,
        bpm_enu_file=inps.bpm_enu_file,
        bpm_los_file=inps.bpm_los_file,
        pmm_reso=inps.pmm_reso,
    )

    if inps.vel_file and inps.bpm_los_file and inps.cor_vel_file:
        print('\n-----------------------------------')
        print('Remove bulk motion')
        print('-----------------------------------\n')
        remove_bulk_motion(inps.vel_file, inps.bpm_los_file, inps.cor_vel_file)
        print('The BPM corrected velocity : {}'.format(inps.cor_vel_file))

    return

################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
