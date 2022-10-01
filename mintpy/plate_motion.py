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


import collections

import numpy as np
from skimage.transform import resize

from mintpy.diff import diff_file
from mintpy.objects.resample import resample
from mintpy.utils import readfile, utils as ut, writefile

# https://docs.astropy.org/en/stable/units/index.html
try:
    from astropy import units as u
    from platemotion import Plate
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


####################################### Sub Functions ##########################################
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
    print(f'  Pole Latitude : {plate.omega_spherical[0].degree:10.4f} deg')
    print(f'  Pole Longitude: {plate.omega_spherical[1].degree:10.4f} deg')
    print(f'  Rotation rate : {plate.omega_spherical[2].to(u.deg/u.Ma):10.4f}  \n')
    print('in Cartesian coordinates:')
    print(f'  wx: {plate.omega_cartesian[0]:10.4f}')
    print(f'  wy: {plate.omega_cartesian[1]:10.4f}')
    print(f'  wz: {plate.omega_cartesian[2]:10.4f}')
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


####################################### Sub Functions ##########################################
def calc_plate_motion(geom_file, omega_cart=None, omega_sph=None, const_vel_enu=None,
                      pmm_enu_file=None, pmm_file=None, pmm_comp='enu2los', pmm_step=10.):
    """Estimate LOS motion due to the rigid plate motion (translation and/or rotation).

    Parameters: geom_file     - str, path to the input geometry file
                omega_cart    - list or 1D array, Cartesian representation of plate rotation
                                in [wx, wy, wz]  (mas/yr)
                omega_sph     - list or 1D array, Spherical representation of plate rotation
                                in [lat, lon, w] (deg, deg, deg/Ma)
                const_vel_enu - list or 1D array, a single-vector [ve, vn, vu] (meter/year)
                                simulating the rigid translation of the ground (e.g., from GNSS)
                pmm_enu_file  - str, path to the output plate motion in east, north, up direction
                pmm_file      - str, path to the output plate motion in LOS direction
                set_comp      - str, output PMM in the given component of interest
                pmm_reso      - float, ground resolution for computing Plate rotation to ENU velocity (km)
    Returns:    ve/vn/vu/vlos - 2D np.ndarray, ridig plate motion in east / north / up / LOS direction
    """

    # Get LOS geometry
    atr_geo = ut.prepare_geo_los_geometry(geom_file, unit='deg')[2]
    shape_geo = [int(atr_geo['LENGTH']), int(atr_geo['WIDTH'])]

    ## plate motion model in the region
    print('-'*50)
    if omega_cart or omega_sph:
        print('compute the rigid plate motion using a plate motion model (PMM; translation & rotation)')
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
        print(f'compute the rigid plate motion using a single vector (translation): {const_vel_enu}')
        ve = const_vel_enu[0] * np.ones(shape_geo, dtype=np.float32)
        vn = const_vel_enu[1] * np.ones(shape_geo, dtype=np.float32)
        vu = const_vel_enu[2] * np.ones(shape_geo, dtype=np.float32)


    # radar-code PMM if input geometry is in radar coordinates
    atr = readfile.read_attribute(geom_file)
    if 'Y_FIRST' not in atr.keys():
        print('radar-coding the rigid plate motion in ENU ...')
        res_obj = resample(lut_file=geom_file)
        res_obj.open()
        res_obj.src_meta = atr_geo
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        ve = res_obj.run_resample(src_data=ve[box[1]:box[3], box[0]:box[2]])
        vn = res_obj.run_resample(src_data=vn[box[1]:box[3], box[0]:box[2]])
        vu = res_obj.run_resample(src_data=vu[box[1]:box[3], box[0]:box[2]])


    # Project model from ENU to direction of interest
    c0, c1 = pmm_comp.split('2')
    print(f'project the ridig plate motion from {c0.upper()} onto {c1.upper()} direction')
    los_inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    los_az_angle = readfile.read(geom_file, datasetName='azimuthAngle')[0]
    unit_vec = ut.get_unit_vector4component_of_interest(los_inc_angle, los_az_angle, comp=pmm_comp)
    vlos = (  ve * unit_vec[0]
            + vn * unit_vec[1]
            + vu * unit_vec[2])

    # Save the plate motion model velocity into HDF5 files
    # metadata
    atr['FILE_TYPE'] = 'velocity'
    atr['DATA_TYPE'] = 'float32'
    atr['UNIT'] = 'm/year'
    for key in ['REF_Y', 'REF_X', 'REF_DATE']:
        if key in atr.keys():
            atr.pop(key)

    if pmm_enu_file:
        # dataset
        dsDict = {
            'east'  : ve,
            'north' : vn,
            'up'    : vu,
        }

        # write
        writefile.write(dsDict, out_file=pmm_enu_file, metadata=atr)

    if pmm_file:
        writefile.write(vlos, out_file=pmm_file, metadata=atr)

    return ve, vn, vu, vlos


def correct_plate_motion(vel_file, mfile, ofile):
    """Apply the plate motion correction from files.
    """
    file1 = vel_file       # input uncorrected LOS velocity file
    file2 = [mfile]        # PMM LOS velocity file
    diff_file(file1, file2, ofile)
    return


################################################################################################
def run_plate_motion(inps):
    """Calculate and/or correct for the rigid motion from tectonic plates."""

    calc_plate_motion(
        geom_file=inps.geom_file,
        omega_cart=inps.omega_cart,
        omega_sph=inps.omega_sph,
        const_vel_enu=inps.const_vel_enu,
        pmm_enu_file=inps.pmm_enu_file,
        pmm_file=inps.pmm_file,
        pmm_comp=inps.pmm_comp,
        pmm_step=inps.pmm_step,
    )

    if inps.vel_file and inps.pmm_file and inps.cor_vel_file:
        print('-'*50)
        print('Correct input velocity for the rigid plate motion')
        correct_plate_motion(inps.vel_file, inps.pmm_file, inps.cor_vel_file)

    return
