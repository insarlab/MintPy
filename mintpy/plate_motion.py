############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Zhang Yunjun, May 2022             #
############################################################
# Recommend usage:
#   from mintpy import plate_motion as pmm
#
# To-Do List (updated 2022.10.12 Yuan-Kai Liu):
#   + Potentially, we can make built-in PMM tables/dictionaries for easier user string input of the plate name
#
# Reference
#   + Pichon, X. L., Francheteau, J. & Bonnin, J. Plate Tectonics; Developments in Geotectonics 6; Hardcover – January 1, 1973. Page 28-29
#   + Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications, Palo Alto. DOI: 10.4236/ojapps.2015.54016. Page 145-156.
#   + [Transformations between ECEF and ENU coordinates](https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates)



import collections
import sys

import numpy as np
import pyproj
from skimage.transform import resize

from mintpy.diff import diff_file
from mintpy.objects.resample import resample
from mintpy.utils import readfile, utils as ut, writefile

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

        # prepare the coarse grid
        latc = float(atr_geo['Y_FIRST']) + float(atr_geo['Y_STEP']) * shape_geo[0] / 2
        ystep = abs(int(pmm_step * 1000 / (float(atr_geo['Y_STEP']) * 108e3)))
        xstep = abs(int(pmm_step * 1000 / (float(atr_geo['X_STEP']) * 108e3 * np.cos(np.deg2rad(latc)))))
        ystep, xstep = max(ystep, 5), max(xstep, 5)
        lats, lons = ut.get_lat_lon(atr_geo, dimension=2, ystep=ystep, xstep=xstep)

        # transform PMM to ENU velocity on a coarse grid
        print(f'compute PMM via matrix rotation: grid_size = {pmm_step} km, grid_shape = {lats.shape} ...')
        lats_flt = lats.flatten()
        lons_flt = lons.flatten()
        if omega_cart is not None:
            au = 'mas/yr' # angular velo unit
            print(f'input omega_cartesian in [wx, wy, wz]: {omega_cart} ({au})')
            Omega = set_Omega(omega_cart, coord='cartesian', unit=au)
        else:
            au = 'deg/Ma' # angular velo unit
            print(f'input omega_spherical in [lat, lon, w]: {omega_sph} (deg, deg, {au})')
            Omega = set_Omega(omega_sph, coord='spherical', unit=au)

        V_enu = Omega2Venu(Omega, lats_flt, lons_flt, alts=0.0, ellps=False)[0]
        print('Done with PMM to ENU calculation')
        ve_low = V_enu[:,0].reshape(lats.shape)
        vn_low = V_enu[:,1].reshape(lats.shape)
        vu_low = V_enu[:,2].reshape(lats.shape)


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
## Sub-functions for the math of Euler Pole and linear velocity

EARTH_RADIUS = 6371.009             #  km (Mean radius from IUGG)
MAS2RAD      = np.pi/3600000/180    #  1 mas (milliarcsecond) = yy radian
MASY2DMY     = 1e6 / 3600000        #  1 mas per year         = xx degree per million year


def sph2cart(lat, lon, r=1):
    """Convert spherical coordinates to cartesian. Default raduis is 1 (unit length).
    INPUT
        lat     latitude   [degree]
        lon     longitude  [degree]
        r       radius     [any units of distance]
    OUTPUT
        u = ndarray([u1, u2, u3])      cartesian vector    [same unit as r]
    EXAMPLE
     (1) spherical coord to xyz coord
        x, y, z = sph2cart(lat, lon, r=radius)
     (2) spherical Euler pole to cartesian rotational vector
        Omega = sph2cart(plat, plon, r=omega)
        where:
            plat    latitude of the pole    [degree]
            plon    longitude of the pole   [degree]
            omega   scalar angular velocity [any units, e.g., milliarcsec (mas) per year]
            Omega   angular velocity vector [same unit as omega]
    """
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    u1 = r * np.cos(lat) * np.cos(lon)       # in x-axis
    u2 = r * np.cos(lat) * np.sin(lon)       # in y-axis
    u3 = r * np.sin(lat)                     # in z-axis
    u  = np.array([u1, u2, u3]).T
    return u


def cart2sph(u):
    """Convert cartesian coordinates to spherical.
    INPUT
        u = ndarray([u1, u2, u3])     cartesian vector     [any units of distance]
    OUTPUT
        lat     latitude   [degree]
        lon     longitude  [degree]
        r       radius     [same unit as u]
    EXAMPLE
     (1) xyz coord to spherical coord
        lat, lon, r = cart2sph(u=[x,y,z])
     (2) cartesian rotational vector to spherical Euler pole
        plat, plon, omega = cart2sph(Omega=[wx, wy, wz])
        where:
            Omega   angular velocity vector [any units, e.g., milliarcsec (mas) per year]
            plat    latitude of the pole    [degree]
            plon    longitude of the pole   [degree]
            omega   scalar angular velocity [same unit as Omega]
    """
    u1, u2, u3 = u
    r   = np.sqrt(u1**2+u2**2+u3**2)
    lat = np.arcsin(u3/r)
    lon = np.arctan2(u2, u1)
    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon)
    return lat, lon, r



def T_cart2enu(lats, lons):
    """Rotation matrix to convert cartesian coordinates (global ECEF) at given (lats,lons) to ENU components (local cartesian).
       Input lats lons are in degree. If input N points of (lats,lons), T.shape = (N, 3, 3)
    INPUT
        lats    1d-array of latitude      [degree]
        lons    1d-array of longitude     [degree]
    OUTPUT
        T       rotation matrix (N, 3, 3) [-]
    REFERENCE
        + https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        + Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications, Palo Alto.
          (DOI: 10.4236/ojapps.2015.54016. Page 145-156)
    """
    lats = np.deg2rad(lats)
    lons = np.deg2rad(lons)
    if not isinstance(lats,(list,np.ndarray)):
        lats = [lats]
        lons = [lons]
    if not len(lats) == len(lons):
        print('Dimension is not the same')
        sys.exit(1)
    else:
        npts = len(lats)
        Te = np.vstack((-np.sin(lons), np.cos(lons), np.zeros(npts))).T
        Tn = np.vstack((-np.sin(lats)*np.cos(lons), -np.sin(lats)*np.sin(lons), np.cos(lats))).T
        Tu = np.vstack((np.cos(lats)*np.cos(lons), np.cos(lats)*np.sin(lons), np.sin(lats))).T
        T  = np.dstack((Te, Tn, Tu))
        T  = np.transpose(T, (0, 2, 1))
    return T


def T_enu2cart(lats, lons):
    """Rotation matrix to convert ENU components (local cartesian) at given (lats,lons) to cartesian coordinates (global ECEF).
       Input lats lons are in degree. If input N points of (lats,lons), T.shape = (N, 3, 3)
    INPUT
        lats    1d-array of latitude      [degree]
        lons    1d-array of longitude     [degree]
    OUTPUT
        T       rotation matrix (N, 3, 3) [-]
    REFERENCE
        + https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        + Cox, A., and Hart, R.B. (1986) Plate tectonics: How it works. Blackwell Scientific Publications, Palo Alto.
          (DOI: 10.4236/ojapps.2015.54016. Page 145-156)
    """
    lats = np.deg2rad(lats)
    lons = np.deg2rad(lons)
    if not isinstance(lats,(list,np.ndarray)):
        lats = [lats]
        lons = [lons]
    if not len(lats) == len(lons):
        print('Dimension is not the same')
        sys.exit(1)
    else:
        npts = len(lats)
        print(f'{npts} points')
        Te = np.vstack((-np.sin(lons), -np.cos(lons)*np.sin(lats), np.cos(lons)*np.cos(lats))).T
        Tn = np.vstack((np.cos(lons), -np.sin(lons)*np.sin(lats), np.sin(lons)*np.cos(lats))).T
        Tu = np.vstack((np.zeros(npts), np.cos(lats), np.sin(lats))).T
        T  = np.dstack((Te, Tn, Tu))
        T  = np.transpose(T, (0, 2, 1))
    return T


def azimuth(east, north):
    """ (Not being used anywhere)
    Returns azimuth in degrees counterclockwise from North given north and
    east components
    """
    azi = np.rad2deg(np.arctan2(north, east))
    azi = 90 - azi
    if azi <= 0:
        azi +=360
    return azi


def geo2ecef_pyproj(lat, lon, alt):
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    x ,y, z = transformer.transform(lon, lat, alt, radians=False)
    return x, y, z


def set_Omega(pole, coord='cartesian', unit='mas/yr'):
    """Given Euler pole, convert to Omega
       The unit of angular velocity must be 'mas per year'.
    INPUT
        pole:
         (1) coord = 'cartesian' description
            Omega     angular velocity vector         unit
         (2) coord = 'spherical' description
            plat      Euler pole latitude             [degree]
            plon      Euler pole longitude            [degree]
            omega     scalar angular velocity         unit

        coord:     ['cartesian', 'spherical']
        unit:      ['mas/yr', 'deg/Ma']
    OUTPUT
        Omega      Cartesian angular velocity vector  [mas per year]
    """
    ## Input the Euler pole
    if coord.lower() == 'cartesian':
        Omega = np.array(pole).astype('float')
        if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
            Omega /= MASY2DMY
        plat, plon, omega = cart2sph(Omega)
    elif coord.lower() == 'spherical':
        plat, plon, omega = np.array(pole).astype('float')
        if unit == 'deg/Ma': #  (consistently using mas/yr for calculation)
            omega /= MASY2DMY
        Omega = sph2cart(plat, plon, r=omega)
    else:
        print('Some user input is wrong here...')
        sys.exit(1)

    print('\n---------------- Full Euler Pole description ---------------')
    print('Euler Pole')
    print(' (1) In spherical expression:')
    print(f'   Pole Latitude : {plat:10.4f} DEG')
    print(f'   Pole Longitude: {plon:10.4f} DEG')
    print(f'   Rotation rate : {omega*MASY2DMY:10.4f} DEG/MA    {omega:10.4f} MAS/Y')
    print(' (2) In Cartesian expression (angular velocity vector):')
    print(f'   wx:             {Omega[0]*MASY2DMY:10.4f} DEG/MA    {Omega[0]:10.4f} MAS/Y')
    print(f'   wy:             {Omega[1]*MASY2DMY:10.4f} DEG/MA    {Omega[1]:10.4f} MAS/Y')
    print(f'   wz:             {Omega[2]*MASY2DMY:10.4f} DEG/MA    {Omega[2]:10.4f} MAS/Y')
    print('------------------------------------------------------------\n')
    return Omega


def Omega2Venu(Omega, lats, lons, alts=0.0, ellps=False):
    """Given Euler pole Omega, compute V_enu for given pixel(s) at (lat,lon) of interest.
       Only supports a constant altitude now.
    INPUT
        Omega      Cartesian angular velocity vector  [mas per year]
        lats:      points of interest (latitude)      [degree]
        lons:      points of interest (longitude)     [degree]
        alts:      points of interest (altitude)      [meter]
        ellps:     True/False; consider ellipsoidal Earth projection of point(s) of interest
    OUTPUT
        V_enu     east, north, up linear velocity     [meter/year]
    """

    ## report how many points of interest
    if isinstance(lats, (list, tuple, np.ndarray)):
        npts = len(lats)
    elif isinstance(lats, (int, float)):
        npts = 1
    print(f'number of points to compute: {npts}')


    ## Local coordinates handling for location(s) of interest
    if not ellps:
        # a perfect sphere
        print(f'Assume a perfect spherical Earth, radius: {EARTH_RADIUS} km')
        locs_xyz = sph2cart(lats, lons, EARTH_RADIUS)                    # unit is km
    else:
        # WGS84 ellips; only supports uniform altitude now, but can change later
        if npts == 1:
            alts = float(alts)
        else:
            alts = alts * np.ones_like(lats)
        print('Assume WGS84 ellipse from pyproj')
        locs_xyz = 1e-3 * np.array(geo2ecef_pyproj(lats, lons, alts)).T  # set unit as km

    ## Compute the cartesian linear velocity (i.e., ECEF)
    # V_ecef = Omega x Ri    where Ri is location vector at pixel i
    V_ecef = np.cross(Omega*MAS2RAD, locs_xyz)                           # watch out unit here!

    ## Convert to local ENU linear velocity
    """
    To-do:
        (1) speed of this big matrix multiplication?
        (2) memory of this big diagonal matrix V_tmp?
    """
    T = T_cart2enu(lats, lons)
    # T is the rotation matrix (can save it to avoid computing again)
    # V = T V_ecef    where V is the local ENU velocity
    V_tmp = np.matmul(T.reshape([-1,3]) , V_ecef.T).reshape([3,npts,npts], order='F')
    V_enu = 1e3 * np.diagonal(V_tmp, axis1=1, axis2=2).T                 # convert km/year to m/year
    del V_tmp
    return V_enu, T


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
