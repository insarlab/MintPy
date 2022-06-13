#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, May 2022                           #
############################################################

# 1. Citation:
# Oliver L. Stephenson, Yuan-Kai Liu, Zhang Yunjun, Mark Simons, and Paul Rosen. (2022)
# The Impact of Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields. Manuscript submitted to Geophysical Research Letters.
# preprint available: https://www.essoar.org/doi/abs/10.1002/essoar.10511538.1

# 2. Extra dependency:
#   + platemotion: https://github.com/lcx366/PlateTectonic
#   + astropy
#   How to install both:
#      option (1) pip install platemotion
#      option (2) git clone git@github.com:lcx366/PlateTectonic.git $TOOL_DIR/PlateTectonic
#                 echo 'export PYTHONPATH=$PYTHONPATH:$TOOL_DIR/PlateTectonic' >> ~/.bashrc
#                 somehow install other dependencies in setup.py using your conda

# 3. To-Do List (updated 2022.5.30 Yuan-Kai Liu):
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
from mintpy.mask import mask_matrix
from mintpy.diff import diff_file
from mintpy import reference_point
from mintpy.save_gmt import get_geo_lat_lon
from mintpy.solid_earth_tides import prepare_los_geometry

try:
    from platemotion import Plate
    from astropy import units as u
except ImportError:
    msg = 'Can NOT import platemotion!'
    msg += '\nCheck more details at https://github.com/lcx366/PlateTectonic.'
    raise ImportError(msg)



#########################################  Usage  ##############################################
REFERENCE = """reference:
  Stephenson, O. L., Liu, Y. K., Yunjun, Z., Simons, M., and Rosen, P., (2022),
  The Impact of Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields, submitted to GRL.
    preprint: https://www.essoar.org/doi/abs/10.1002/essoar.10511538.1
"""
EXAMPLE = """example:
  bulk_plate_motion.py -g inputs/geometryGeo.h5                --om_sph  50.6  -74.0   0.30   -m None
  \testimating for NNR-NUVEL1A PMM of the Africa plate\n\t(Table 2 of Argus & Gordon, GRL 1991; doi: 10.1029/91GL01532)\n
  bulk_plate_motion.py -g inputs/geometryGeo.h5 -v velocity.h5 --om_sph  48.85 -106.50 0.223  -m waterMask.h5
  \testimating/correcting for NNR-MORVEL56 PMM of the Eurasia plate\n\t(Table 1 of Argus, Gordon, and DeMets, GGG 2011; doi: 10.1029/2011GC003751)\n
  bulk_plate_motion.py -g inputs/geometryGeo.h5 -v velocity.h5 --om_cart 1.154 -0.136  1.444  -m maskTempCoh.h5
  \testimating/correcting for NNR-ITRF14 PMM of the Arabia plate\n\t(Table 1 of Altamimi et al., GJI 2017; doi: 10.1093/gji/ggx136)\n
  bulk_plate_motion.py -g inputs/geometryGeo.h5 -v velocity.h5 --enu     25.0 30.5 0.0        -m zero
  \testimating/correcting for a simple constant local ENU translating field based on one GNSS vector
  \t(e.g., https://www.unavco.org/software/visualization/GPS-Velocity-Viewer/GPS-Velocity-Viewer.html)
  \tyou can: 1. select `GNSS Data source` as `World, IGS08/NNR, GEM GSRM` (referenced to ITRF2008, NNR PMM)
  \t         2. check box `Station labels and data download` and click `Draw Map`
  \t         3. Navigate to the region of interest and click on a representative station.
  \t            You get the `Speed components` in mm/yr.
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Bulk plate motion correction. Removing the effect of bulk traslation and rotation in velocity field based on a given plate motion model (PMM).'+
                                                 ' Sentinel-1 orbit is measured with respect to ITRF2014 (page 25 of Peter, 2021: https://sentinel.esa.int/documents/247904/4599719/Copernicus-POD-Product-Handbook.pdf),'
                                                 ' an Earth-centered, Earth-fixed reference frame in which there is no net rotation (NNR) of the Earth surface.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)
    # input files
    parser.add_argument('-g', '--geom', dest='geomfile', type=str, required=True,
                        help = 'Input geometry file (required), e.g., geometryGeo.h5')
    parser.add_argument('-v', '--velo', dest='vfile', type=str, default=None,
                        help='Input velocity file (optional), e.g., velocity.h5 (default: %(default)s).')

    # plate motion configurations (use ONE AND ONLY ONE of the below; --om_sph overwrites --om_cart overwrites --enu)
    pmms = parser.add_argument_group('pmms', 'Plate motion models (required); Use ONE AND ONLY ONE of the below; '+
                                     '\nIf more than one are given, --om_sph overwrites --om_cart overwrites --enu')
    pmms.add_argument('--enu', dest='venu', type=float, nargs=3, metavar=('VE', 'VN', 'VU'), default=None,
                        help = 'Simple constant local ENU translation of ground [ve, vn, vu] unit: meter/year (default: %(default)s).')
    pmms.add_argument('--om_cart', dest='omega_cart', type=float, nargs=3, metavar=('WX', 'WY', 'WZ'), default=None,
                        help = 'Cartesian form of Euler Pole rotation; [wx, wy, wz] (unit: mas/yr) (default: %(default)s).')
    pmms.add_argument('--om_sph', dest='omega_sph', type=float, nargs=3, metavar=('LAT', 'LON', 'W'), default=None,
                        help = 'Spherical form of Euler Pole rotation; [lat, lon, w] (unit: deg, deg, deg/Ma) (default: %(default)s).')

    # others
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None,
                        help = 'Mask file to apply to the inout and predicted velocity fields. \n' +
                                '   zero: masking 0.0 values of the input velocity \n' +
                                '   None: no mask is applied (default: %(default)s).')
    parser.add_argument('--resol', dest='resol', type=float, default=10.,
            help = 'Ground resolution for computing Plate rotation to ENU velocity (unit: km) (default: %(default)s km grid).')

    return parser


def get_filenames(geomfile, vfile):
    inputsDir = os.path.dirname(geomfile)
    BPM3d     = os.path.join(inputsDir, 'BulkPlateMotion3D.h5')
    BPMlos    = os.path.join(inputsDir, 'BulkPlateMotion.h5')
    if vfile:
        vbase = os.path.abspath(vfile).split('.')[0]
        vout  = os.path.abspath('{}_BPM.h5'.format(vbase))
    else:
        vout  = None
    return BPM3d, BPMlos, vout


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if (inps.venu is None) and (inps.omega_cart is None) and (inps.omega_sph is None):
        print('Error: Please specify either the Plate Rotation or a Single Translation Vector; See -h option')
        print('Need to give at least ONE of [--enu, --om_cart, --om_sph]')
        sys.exit(1)
    inps.BPM3d, inps.BPMlos, inps.vout = get_filenames(inps.geomfile, inps.vfile)
    return inps

########################################## Sub Functions #############################################

def build_PMM(omega_cart=None, omega_sph=None):
    """
    Build a plate motion model by giving a Euler roation vector
    omega_cart      list or np.array;   Cartesian representation [wx, wy, wz]  (mas/yr)
    omega_sph       list or np.array;   Spherical representation [lat, lon, w] (deg, deg, deg/Ma)
    """

    # Set a NaN moment of inertia (to get class Plate running)
    iner_null = np.zeros([3,3]) * np.nan * (u.kg * u.km**2)

    # Create an instrance of class Plate from `platemotion` pkg
    plate = Plate(info={'inertia_tensor': iner_null})

    # Check input variables
    if (omega_cart is None) and (omega_sph is None):
        print('Need to give either omega_cartesian (wxyz) or omega_spherical (Euler Pole)!')
        sys.exit(1)

    elif omega_cart is not None:
        print('your input: omega_cartesian; [wx, wy, wz] (mas/yr)')
        omega = np.array(omega_cart) * u.mas/u.yr
        plate.set_omega(omega,'cartesian')

    else:
        print('your input: omega_spherical; Euler pole; [lat, lon, w] (deg, deg, deg/Ma)')
        omega = [omega_sph[0]*u.deg, omega_sph[1]*u.deg, omega_sph[2]*u.deg/u.Ma]
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
    """
    Interpolate 3-layer 2D arrays individually
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
    """
    Get width and length of the geo box in km
    geo_box     tuple of 4 floats (lon_min, lat_max, lon_max, lat_min) in decimal degrees
    """
    geod = pyproj.Geod(ellps='WGS84')
    width_km  = geod.inv(geo_box[0], geo_box[3], geo_box[2], geo_box[3])[2] * 1e-3
    length_km = geod.inv(geo_box[0], geo_box[3], geo_box[0], geo_box[1])[2] * 1e-3
    return width_km, length_km


####################################### Higher-level Sub Functions ##########################################

def estimate_bulkMotion(geomfile, vfile, venu=None, omega_cart=None, omega_sph=None, mask=None, pmmResol=5., BPM3d=None, BPMlos=None):
    """
    Estimate LOS motion due to pure bulk tranlation or due to plate rotation
    geomfile        str:                path to the input geometry file
    vfile           str;                path to the input velocity file
    venu            list or 1D array;   a single-vector [ve, vn, vu] (meter/year) simulating the bulk translation of the ground (e.g., from GNSS)
    omega_cart      list or 1D array;   Cartesian representation of plate rotation [wx, wy, wz]  (mas/yr)
    omega_sph       list or 1D array;   Spherical representation of plate rotation [lat, lon, w] (deg, deg, deg/Ma)
    mask            str;                a mask file for the bulk motion. Can use 'zero' for masking 0.0; or None for no masking
    pmmResol        float;              ground resolution for computing Plate rotation to ENU velocity (km); default is at 5.0 km grid
    BPM3d           str;                path to the output BPM (bulk plate motion) east, north, up velocity field
    BPMlos          str;                path to the output BPM (bulk plate motion) LOS velocity field
    """
    # Read attributes / reference file
    if vfile:
        vel_in = readfile.read(vfile, datasetName='velocity')[0]
        atr    = readfile.read(vfile)[1]
    else:
        atr = readfile.read(geomfile)[1]
        atr['FILE_TYPE'] = 'velocity'

    width, length = int(atr['WIDTH']), int(atr['LENGTH'])

    # Get LOS geometry
    inc_rad, head_rad, atr_geo = prepare_los_geometry(geomfile)
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
            plate = build_PMM(omega_cart=omega_cart)
        elif omega_sph is None:
            plate = build_PMM(omega_sph=omega_sph)

        # Use the Euler rotation to compute the ENU for each pixel in the geometry
        skx, sky = int(width*pmmResol/width_km), int(length*pmmResol/length_km)
        Lats = np.meshgrid(lons, lats)[1][::sky, ::skx]
        Lons = np.meshgrid(lons, lats)[0][::sky, ::skx]
        print('Compute ENU from PMM with a ~{} km ground resolution'.format(pmmResol))
        print('low resolution / original grid dimension: ({}, {}) / ({}, {})'.format(Lats.shape[0], Lats.shape[1], length, width))

        # transform PMM to ENU velocity on a coarse grid
        # to-do: multi-pixel rotation is slow; change the `platemotion` code to Big Array operation rather than For Loops
        V_enu_low = pmm2enu_at(plate, Lats, Lons)

        # interpolate back to the initial grid
        print('interpolate back to the original grid')
        V_enu = interp_2d3l(V_enu_low, Lons, Lats, width, length, 'linear')

    elif venu:
        print('estimate bulk motion by assuming a single pure translation vector: {}'.format(venu))
        V_enu = venu * np.ones([length, width, 3])

    else:
        print('Error: Please specify either the Plate Rotation or a Single Translation Vector; See -h option')
        print('Need to give at least ONE of [--enu, --om_cart, --om_sph]')
        sys.exit(1)

    # Project model to LOS velocity
    print('Project the bulk motion ENU onto radar LOS velocity')
    ve, vn, vu = np.array((V_enu[:,:,0])), np.array((V_enu[:,:,1])), np.array((V_enu[:,:,2]))
    vlos = ut.enu2los(ve, vn, vu, inc_angle=inc_deg, head_angle=head_deg)

    ## Masking the bulk motion?
    if mask:
        if mask == 'zero':
            if vel_in:
                print('Maksing zeros based on the input reference velocity file')
                msk_bool = 1*(vel_in!=0.0)
                msk_bool[int(atr['REF_Y']), int(atr['REF_X'])] = 1
            else:
                msk_bool = 1*(vlos!=0.0)
        else:
            print('Masking with the given mask file: {}'.format(mask))
            msk_bool = readfile.read(mask)[0]
        vlos = mask_matrix(vlos, msk_bool, fill_value=np.nan)

    # Save the bulk motion model velocity into two files
    dsDict = {'east': ve, 'north': vn, 'up': vu}        # 3D velocity
    writefile.write(dsDict, out_file=BPM3d, metadata=atr)
    dsDict = {'velocity': vlos}                         # LOS velocity
    writefile.write(dsDict, out_file=BPMlos, metadata=atr)

    # Reset the reference point for both files; they are absolute velocity models
    for file in [BPM3d, BPMlos]:
        iargs = [file, '--reset']
        reference_point.main(iargs)
    return


def remove_bulkMotion(vfile, mfile, ofile):
    """
    Apply the bulk motion correction from files
    """
    file1 = vfile       # input uncorrected LOS velocity file
    file2 = [mfile]     # BPM model LOS velocity file
    diff_file(file1, file2, ofile)
    return


#######################################  Main Function  ########################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n-----------------------------------')
    print('Evaluate bulk motion')
    print('-----------------------------------\n')
    estimate_bulkMotion(inps.geomfile, inps.vfile,
                        inps.venu, inps.omega_cart, inps.omega_sph,
                        inps.mask, inps.resol, inps.BPM3d, inps.BPMlos)

    if inps.vfile and inps.BPMlos and inps.vout:
        print('\n-----------------------------------')
        print('Remove bulk motion')
        print('-----------------------------------\n')
        remove_bulkMotion(inps.vfile, inps.BPMlos, inps.vout)
        print('The BPM corrected velocity           : {}'.format(inps.vout))

    print('The absolute BPM velocity model (3D) : {}'.format(inps.BPM3d))
    print('The absolute BPM velocity model (LOS): {}'.format(inps.BPMlos))
    print('Normal complete')
    return

################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])