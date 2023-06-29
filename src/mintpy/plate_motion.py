############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Zhang Yunjun, May 2022             #
############################################################
# Recommend usage:
#   from mintpy import plate_motion as pmm
#
# Reference:
#   Stephenson, O. L., Liu, Y. K., Yunjun, Z., Simons, M., Rosen, P. and Xu, X., (2022),
#     The Impact of Plate Motions on Long-Wavelength InSAR-Derived Velocity Fields,
#     Geophys. Res. Lett. 49, e2022GL099835, doi:10.1029/2022GL099835.


import numpy as np
from skimage.transform import resize

from mintpy.diff import diff_file
from mintpy.objects.euler_pole import ITRF2014_PMM, EulerPole
from mintpy.objects.resample import resample
from mintpy.utils import readfile, utils as ut, writefile

####################################### Major Function ###########################################

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
                pmm_comp      - str, output PMM in the given component of interest
                pmm_step      - float, ground resolution for computing Plate rotation to ENU velocity (km)
    Returns:    ve/vn/vu/vlos - 2D np.ndarray, ridig plate motion in east / north / up / LOS direction
    """

    # Get LOS geometry
    atr_geo = ut.prepare_geo_los_geometry(geom_file, unit='deg')[2]
    shape_geo = [int(atr_geo['LENGTH']), int(atr_geo['WIDTH'])]

    ## calc plate motion in the region
    print('-'*50)
    if omega_cart or omega_sph:
        print('compute the rigid plate motion defined as an Euler Pole')

        # construct Euler Pole object
        if omega_cart is not None:
            print(f'input omega_cartesian in [wx, wy, wz]: {omega_cart} [mas/yr]')
            pole_obj = EulerPole(
                wx=omega_cart[0],
                wy=omega_cart[1],
                wz=omega_cart[2],
                unit='mas/yr',
            )

        else:
            print(f'input omega_spherical in [lat, lon, w]: {omega_sph} [deg, deg, deg/Ma]')
            pole_obj = EulerPole(
                pole_lat=omega_sph[0],
                pole_lon=omega_sph[1],
                rot_rate=omega_sph[2],
                unit='deg/Ma',
            )
        pole_obj.print_info()

        # prepare the coarse grid (for the points of interest)
        latc = float(atr_geo['Y_FIRST']) + float(atr_geo['Y_STEP']) * shape_geo[0] / 2
        ystep = abs(int(pmm_step * 1000 / (float(atr_geo['Y_STEP']) * 108e3)))
        xstep = abs(int(pmm_step * 1000 / (float(atr_geo['X_STEP']) * 108e3 * np.cos(np.deg2rad(latc)))))
        ystep, xstep = max(ystep, 5), max(xstep, 5)
        lats, lons = ut.get_lat_lon(atr_geo, dimension=2, ystep=ystep, xstep=xstep)
        print(f'calculate plate motion on the coarse grid: size = ~{pmm_step} km, shape = {lats.shape}')

        # calculate plate motion in ENU at the coarse grid
        ve_low, vn_low, vu_low = pole_obj.get_velocity_enu(lats, lons, alt=0.0, ellps=True)

        # for debugging purpose
        debug_mode = False
        if debug_mode:
            from matplotlib import pyplot as plt

            # calculate plate motion in ECEF (XYZ) coordinates
            vx, vy, vz = pole_obj.get_velocity_xyz(lats, lons, alt=0.0, ellps=True)

            # plot
            fig, axs = plt.subplots(nrows=2, ncols=3, figsize=[12, 6])
            vlist = [vx, vy, vz, ve_low, vn_low, vu_low]
            titles = ['X', 'Y', 'Z', 'E', 'N', 'U']
            for ax, data, title in zip(axs.flatten(), vlist, titles):
                im = ax.imshow(data, interpolation='nearest')
                fig.colorbar(im, ax=ax)
                ax.set_title(title)
            fig.tight_layout()
            plt.show()

        # resample coarse grid back to the initial fine grid
        print(f'resample plate motion from corase back to original grid: {lats.shape} -> {shape_geo}'
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

    else:
        raise ValueError('No plate motion configuration (--om-cart/sph or --enu) found!')


    # radar-code the plate motion if input geometry is in radar coordinates
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


    ## project Plate motion from ENU to direction of interest, e.g. LOS or az
    c0, c1 = pmm_comp.split('2')
    print(f'project the ridig plate motion from {c0.upper()} onto {c1.upper()} direction')
    los_inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    los_az_angle = readfile.read(geom_file, datasetName='azimuthAngle')[0]
    unit_vec = ut.get_unit_vector4component_of_interest(los_inc_angle, los_az_angle, comp=pmm_comp)
    vlos = (  ve * unit_vec[0]
            + vn * unit_vec[1]
            + vu * unit_vec[2])

    # save the plate motion model velocity into HDF5 files
    # metadata
    atr['FILE_TYPE'] = 'velocity'
    atr['DATA_TYPE'] = 'float32'
    atr['UNIT'] = 'm/year'
    for key in ['REF_Y', 'REF_X', 'REF_DATE']:
        if key in atr.keys():
            atr.pop(key)

    if pmm_enu_file:
        # dataset
        dsDict = {'east'  : ve,
                  'north' : vn,
                  'up'    : vu}
        # write
        writefile.write(dsDict, out_file=pmm_enu_file, metadata=atr)

    if pmm_file:
        writefile.write(vlos, out_file=pmm_file, metadata=atr)

    return ve, vn, vu, vlos


################################################################################################
def run_plate_motion(inps):
    """Calculate and/or correct for the rigid motion from tectonic plates."""

    # check: --plate option (convert plate_name to omega_cart)
    if inps.plate_name:
        plate = ITRF2014_PMM[inps.plate_name]
        inps.omega_cart = [plate.omega_x, plate.omega_y, plate.omega_z]
        msg = f'get rotation parameters for {inps.plate_name} plate from Table 1 in Altamimi et al. (2017): '
        msg += f'wx, wy, wz = {plate.omega_x}, {plate.omega_y}, {plate.omega_z} mas/yr'
        print(msg)

    # calculate plate motion
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

    # correct plate motion from input velocity
    if inps.vel_file and inps.pmm_file and inps.cor_vel_file:
        print('-'*50)
        print('Correct input velocity for the rigid plate motion')
        diff_file(inps.vel_file, [inps.pmm_file], inps.cor_vel_file)

    return
