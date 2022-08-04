#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import sys
import numpy as np
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser


################################################################################
REFERENCE = """reference:
  Fialko, Y., Simons, M., & Agnew, D. (2001). The complete (3-D) surface displacement
    field in the epicentral area of the 1999 MW7.1 Hector Mine Earthquake, California,
    from space geodetic observations. Geophysical Research Letters, 28(16), 3063-3066.
    doi:10.1029/2001GL013174
  Wright, T. J., B. E. Parsons, and Z. Lu (2004), Toward mapping surface deformation
    in three dimensions using InSAR, Geophysical Research Letters, 31(1), L01607,
    doi:10.1029/2003GL018827.
"""

EXAMPLE = """example:
  # for data with different spatial resolution and coverage
  # use geocode.py -x/y --bbox option to make them consistent
  cd AlosAT424/mintpy
  mask.py velocity.h5 -m maskTempCoh.h5
  geocode.py velocity_msk.h5 -l inputs/geometryRadar.h5 -x 0.00027778 -y -0.00027778 --bbox 32.0 32.5 130.1 130.5
  cd AlosDT73/mintpy
  mask.py velocity.h5 -m maskTempCoh.h5
  geocode.py velocity_msk.h5 -l inputs/geometryRadar.h5 -x 0.00027778 -y -0.00027778 --bbox 32.0 32.5 130.1 130.5
  asc_desc2horz_vert.py AlosAT424/mintpy/geo_velocity_msk.h5 AlosDT73/mintpy/geo_velocity_msk.h5

  # write horz/vert to two files
  asc_desc2horz_vert.py AlosAT424/mintpy/velocity_msk.h5 AlosDT73/mintpy/velocity_msk.h5
  asc_desc2horz_vert.py AlosAT424/mintpy/velocity_msk.h5 AlosDT73/mintpy/velocity_msk.h5  --azimuth 16
  asc_desc2horz_vert.py AlosAT424/mintpy/velocity_msk.h5 AlosDT73/mintpy/velocity_msk.h5  --dset step20200107

  # write all asc/desc/horz/vert datasets into one file
  asc_desc2horz_vert.py Alos2AT131/mintpy/20171219_20190702.unw Alos2DT23/mintpy/20171211_20190819.unw --oo Kirishima2017post.h5
  view.py Kirishima2017post.h5 -u cm --wrap --wrap-range -5 5  #check deformation signal with multiple viewing geometries.

  # pixel-wise decomposition [for large area analysis]
  asc_desc2horz_vert.py asc_velocity.h5 desc_velocity.h5 -g asc_geometry.h5 desc_geometry.h5
"""


def create_parser(subparsers=None):
    synopsis = 'Project Asc and Desc LOS displacement to Horizontal and Vertical direction'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)
 
    # input files
    parser.add_argument('file', nargs=2,
                        help='Ascending and descending files\n'
                             'Both files need to be geocoded in the same spatial resolution.')
    parser.add_argument('-d', '--dset', dest='ds_name', type=str, help='dataset to use, default: 1st dataset')
    parser.add_argument('-g','--geom-file', dest='geom_file', nargs=2, help='Geometry files for the input data files.')

    # inputs - checking
    parser.add_argument('--max-ref-yx-diff', dest='max_ref_yx_diff', type=int, default=3,
                        help='Maximum difference between REF_Y/X (derived from REF_LAT/LON) of input files '+
                             '(default: %(default)s).')

    # outputs - horizontal direction of interest
    parser.add_argument('--az','--horz-az-angle', dest='horz_az_angle', type=float, default=-90.0,
                        help='Azimuth angle in degrees of the interested horizontal direction (default: %(default)s).\n'
                             'Measured from the north with positive for anti-clockwise direction.\n'
                             'E.g.: -90 for East direction\n'
                             '      0   for North direction\n'
                             'Set to the azimuth angle of the strike-slip fault to measure the fault-parallel displacement.\n'
                             'Note:\n'
                             'a. This assumes no deformation in its perpendicular direction\n'
                             'b. Near north direction can not be well resolved due to the lack of\n'
                             '   diversity in viewing geometry. Check exact dilution of precision for \n'
                             '   each component in Wright et al. (2004, GRL)')

    # output - data files
    parser.add_argument('-o', '--output', dest='outfile', nargs=2, metavar=('HZ_FILE','UP_FILE'), default=['hz.h5', 'up.h5'],
                        help='output file name for vertical and horizontal components')
    parser.add_argument('--oo','--one-output', dest='one_outfile',
                        help='Stack the input/output files into one HDF5 file.\n' +
                             'This will disable the HZ/UP_FILE output option.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr1 = readfile.read_attribute(inps.file[0])
    atr2 = readfile.read_attribute(inps.file[1])

    # check 1 - geo-coordinates
    if any('X_FIRST' not in i for i in [atr1, atr2]):
        raise Exception('Not all input files are geocoded.')

    # check 2 - spatial resolution
    if any(atr1[i] != atr2[i] for i in ['X_STEP','Y_STEP']):
        msg  = '\tfile1: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[0], atr1['Y_STEP'], atr1['X_STEP'], atr1.get('X_UNIT', 'degrees'))
        msg += '\tfile2: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[1], atr2['Y_STEP'], atr2['X_STEP'], atr2.get('X_UNIT', 'degrees'))
        msg += '\tRe-run geocode.py --lat-step --lon-step to make them consistent.'
        raise ValueError('input files do NOT have the same spatial resolution\n{}'.format(msg))

    # check 3 - reference point
    ref_lat1, ref_lon1 = [float(atr1[i]) for i in ['REF_LAT', 'REF_LON']]
    ref_lat2, ref_lon2 = [float(atr2[i]) for i in ['REF_LAT', 'REF_LON']]
    ref_y_diff = abs((ref_lat1 - ref_lat2) / float(atr1['Y_STEP']))
    ref_x_diff = abs((ref_lon1 - ref_lon2) / float(atr1['X_STEP']))
    if any(ref_diff > inps.max_ref_yx_diff for ref_diff in [ref_y_diff, ref_x_diff]):
        msg = 'REF_LAT/LON difference between input files > {} pixels!\n'.format(inps.max_ref_yx_diff)
        for fname, ref_lat, ref_lon in zip(inps.file, [ref_lat1, ref_lat2], [ref_lon1, ref_lon2]):
            msg += 'file1: {}\n'.format(fname)
            msg += '\tREF_LAT/LON: [{:.8f}, {:.8f}]\n'.format(ref_lat, ref_lon)
        raise ValueError(msg)

    return inps


################################################################################
def get_overlap_lalo(atr_list):
    """Find overlap area in lat/lon of geocoded files based on their metadata.
    Parameters: atr_list - list of dict, attribute dictionary of two input files in geo coord
    Returns:    S/N/W/E  - float, West/East/South/North in deg
    """
    S, N, W, E = None, None, None, None
    for i, atr in enumerate(atr_list):
        Si, Ni, Wi, Ei = ut.four_corners(atr)
        if i == 0:
            S, N, W, E = Si, Ni, Wi, Ei
        else:
            S = max(Si, S)
            N = min(Ni, N)
            W = max(Wi, W)
            E = min(Ei, E)

    return S, N, W, E


def get_design_matrix4east_north_up(los_inc_angle, los_az_angle, obs_direction):
    """Design matrix G to convert multi-track range/azimuth displacement into east/north/up direction.
    Parameters: los_inc_angle - 1D np.ndarray in size of (num_obs,) in float32, LOS incidence angle in degree
                los_az_angle  - 1D np.ndarray in size of (num_obs,) in float32, LOS azimuth   angle in degree
                obs_direction - 1D np.ndarray in size of (num_obs,) in str, observation direction: range or azimuth
    Returns:    G             - 2D np.ndarray in size of (num_obs, 3) in float32, design matrix
    """
    num_obs = los_inc_angle.shape[0]
    G = np.zeros((num_obs, 3), dtype=np.float32)

    for i, (inc_angle, az_angle, obs_dir) in enumerate(zip(los_inc_angle, los_az_angle, obs_direction)):
        # calculate the unit vector
        if obs_dir == 'range':
            # for range offset / InSAR phase [with positive value for motion toward the satellite]
            ve = np.sin(np.deg2rad(inc_angle)) * np.sin(np.deg2rad(az_angle)) * -1
            vn = np.sin(np.deg2rad(inc_angle)) * np.cos(np.deg2rad(az_angle))
            vu = np.cos(np.deg2rad(inc_angle))

        elif obs_dir == 'azimuth':
            # for azimuth offset [with positive value for motion same as flight]
            ve = np.sin(np.deg2rad(az_angle - 90)) * -1
            vn = np.cos(np.deg2rad(az_angle - 90)) * 1
            vu = 0.

        else:
            raise ValueError(f'un-recognized observation direction: {obs_dir}')

        # fill the design matrix
        G[i, :] = [ve, vn, vu]

    return G


def get_design_matrix4horz_vert(los_inc_angle, los_az_angle, horz_az_angle=-90):
    """Design matrix G to convert asc/desc range displacement into horz/vert direction.
    Only asc + desc -> hz + up is implemented for now.

    Project displacement from LOS to Horizontal and Vertical components:
    Math for 3D:
        dLOS =   dE * sin(inc_angle) * sin(az_angle) * -1
               + dN * sin(inc_angle) * cos(az_angle)
               + dU * cos(inc_angle)
    Math for 2D:
        dLOS =   dH * sin(inc_angle) * cos(az_angle - az)
               + dV * cos(inc_angle)
        with dH_perp = 0.0
    This could be easily modified to support multiple view geometry
        (e.g. two adjcent tracks from asc & desc) to resolve 3D

    Parameters: los_inc_angle - 1D np.ndarray in size of (num_file), LOS incidence angle in degree.
                los_az_angle  - 1D np.ndarray in size of (num_file), LOS azimuth   angle in degree.
                horz_az_angle - float, azimuth angle for the horizontal direction of interest in degree.
                                Measured from the north with anti-clockwise direction as positive.
    Returns:    G             - 2D matrix in size of (num_file, 2)
    """
    num_file = los_inc_angle.shape[0]
    G = np.zeros((num_file, 2), dtype=np.float32)
    for i in range(num_file):
        G[i, 0] = np.sin(np.deg2rad(los_inc_angle[i])) * np.cos(np.deg2rad(los_az_angle[i] - horz_az_angle))
        G[i, 1] = np.cos(np.deg2rad(los_inc_angle[i]))

    return G


def asc_desc2horz_vert(dlos, los_inc_angle, los_az_angle, horz_az_angle=-90, step=20):
    """Decompose asc / desc LOS data into horz / vert data.
    Parameters: dlos          - 3D np.ndarray in size of (num_file, length, width), LOS displacement in meters.
                los_inc_angle - 1/3D np.ndarray in size of (num_file), length, width), LOS incidence angle in degree.
                los_az_angle  - 1/3D np.ndarray in size of (num_file), length, width), LOS azimuth   angle in degree.
                horz_az_angle - float, horizontal azimuth angle of interest in degree.
                step          - int, geometry step size
    Returns:    dhorz         - 2D np.ndarray in size of (length, width), horizontal displacement in meters.
                dvert         - 2D np.ndarray in size of (length, width), vertical   displacement in meters.
    """
    # initiate output
    (num_file, length, width) = dlos.shape
    dhorz = np.zeros((length, width), dtype=np.float32) * np.nan
    dvert = np.zeros((length, width), dtype=np.float32) * np.nan

    # 0D (constant) incidence / azimuth angle --> invert once for all pixels
    if los_inc_angle.ndim == 1:
        G = get_design_matrix4horz_vert(los_inc_angle, los_az_angle, horz_az_angle)
        print('decomposing asc/desc into horz/vert direction ...')
        dhv = np.dot(np.linalg.pinv(G), dlos.reshape(num_file, -1)).astype(np.float32)
        dhorz = dhv[0, :].reshape(length, width)
        dvert = dhv[1, :].reshape(length, width)

    # 2D incidence / azimuth angle --> invert window-by-window
    elif los_inc_angle.ndim == 3:
        num_row = np.ceil(length / step).astype(int)
        num_col = np.ceil(width / step).astype(int)

        print(f'decomposing asc/desc into horz/vert direction in windows of {step}x{step} ...')
        prog_bar = ptime.progressBar(maxValue=num_row)
        for i in range(num_row):
            y0, y1 = step * i, min(step * (i + 1), length)
            for j in range(num_col):
                x0, x1 = step * j, min(step * (j + 1), width)

                # calculate the median geometry for the local window
                med_los_inc_angle = np.nanmedian(los_inc_angle[:, y0:y1, x0:x1], axis=(1,2))
                med_los_az_angle  = np.nanmedian( los_az_angle[:, y0:y1, x0:x1], axis=(1,2))
                if np.all(~np.isnan(med_los_inc_angle)):

                    G = get_design_matrix4horz_vert(med_los_inc_angle, med_los_az_angle, horz_az_angle)
                    dhv = np.dot(np.linalg.pinv(G), dlos[:, y0:y1, x0:x1].reshape(num_file, -1))
                    dhorz[y0:y1, x0:x1] = dhv[0].reshape(y1-y0, x1-x0)
                    dvert[y0:y1, x0:x1] = dhv[1].reshape(y1-y0, x1-x0)

            prog_bar.update(i+1, suffix=f'{i+1}/{num_row}')
        prog_bar.close()

    else:
        raise ValueError(f'un-supported incidence angle matrix dimension ({los_inc_angle.ndim})!')

    return dhorz, dvert


def run_asc_desc2horz_vert(inps):
    """Decompose asc / desc LOS files into horz / vert file(s).
    Parameters: inps         - namespace, input parameters
    Returns:    inps.outfile - str(s) output file(s)
    """

    ## 1. calculate the overlaping area in lat/lon
    atr_list = [readfile.read_attribute(fname, datasetName=inps.ds_name) for fname in inps.file]
    S, N, W, E = get_overlap_lalo(atr_list)
    lat_step = float(atr_list[0]['Y_STEP'])
    lon_step = float(atr_list[0]['X_STEP'])
    length = int(round((S - N) / lat_step))
    width  = int(round((E - W) / lon_step))
    print('overlaping area in SNWE: {}'.format((S, N, W, E)))


    ## 2. read LOS data and geometry
    num_file = len(inps.file)
    num_pixel = length * width
    dlos = np.zeros((num_file, length, width), dtype=np.float32)
    if inps.geom_file:
        los_inc_angle = np.zeros((num_file, length, width), dtype=np.float32)
        los_az_angle  = np.zeros((num_file, length, width), dtype=np.float32)
    else:
        los_inc_angle = np.zeros(num_file, dtype=np.float32)
        los_az_angle  = np.zeros(num_file, dtype=np.float32)

    for i, (atr, fname) in enumerate(zip(atr_list, inps.file)):
        # overlap SNWE --> box to read for each specific file
        coord = ut.coordinate(atr)
        x0 = coord.lalo2yx(W, coord_type='lon')
        y0 = coord.lalo2yx(N, coord_type='lat')
        box = (x0, y0, x0 + width, y0 + length)

        # read data
        dlos[i, :] = readfile.read(fname, box=box, datasetName=inps.ds_name)[0]
        msg = f'{inps.ds_name} ' if inps.ds_name else ''
        print(f'read {msg} from file: {fname}')

        # read geometry
        if inps.geom_file:
            los_inc_angle[i, :] = readfile.read(inps.geom_file[i], box=box, datasetName='incidenceAngle')[0]
            los_az_angle[i, :]  = readfile.read(inps.geom_file[i], box=box, datasetName='azimuthAngle')[0]
            print(f'read 2D LOS incidence / azimuth angles from file: {inps.geom_file[i]}')
        else:
            los_inc_angle[i] = ut.incidence_angle(atr, dimension=0, print_msg=False)
            los_az_angle[i] = ut.heading2azimuth_angle(float(atr['HEADING']))
            print('calculate the constant LOS incidence / azimuth angles from metadata as:')
            print(f'LOS incidence angle: {los_inc_angle[i]:.1f} deg')
            print(f'LOS azimuth   angle: {los_az_angle[i]:.1f} deg')


    ## 3. decompose LOS displacements into horizontal / Vertical displacements
    print('---------------------')
    dhorz, dvert = asc_desc2horz_vert(dlos, los_inc_angle, los_az_angle, inps.horz_az_angle)


    ## 4. write outputs
    print('---------------------')
    # Update attributes
    atr = atr_list[0].copy()
    if inps.ds_name and atr['FILE_TYPE'] in ['ifgramStack', 'timeseries', 'HDFEOS']:
        atr['FILE_TYPE'] = 'displacement'

    atr['WIDTH']  = str(width)
    atr['LENGTH'] = str(length)
    atr['X_STEP'] = str(lon_step)
    atr['Y_STEP'] = str(lat_step)
    atr['X_FIRST'] = str(W)
    atr['Y_FIRST'] = str(N)

    # update REF_X/Y
    ref_lat, ref_lon = float(atr['REF_LAT']), float(atr['REF_LON'])
    [ref_y, ref_x] = ut.coordinate(atr).geo2radar(ref_lat, ref_lon)[0:2]
    atr['REF_Y'] = int(ref_y)
    atr['REF_X'] = int(ref_x)

    # use ref_file for time-series file writing
    ref_file = inps.file[0] if atr_list[0]['FILE_TYPE'] == 'timeseries' else None

    if inps.one_outfile:
        print('write asc/desc/horz/vert datasets into {}'.format(inps.one_outfile))
        dsDict = {}
        for i, atr_i in enumerate(atr_list):
            # dataset name for LOS data
            track_num = atr_i.get('trackNumber', None)
            proj_name = atr_i.get('PROJECT_NAME', None)
            if proj_name in ['none', 'None', None]:
                proj_name = atr_i.get('FILE_PATH', None)
            proj_name = sensor.project_name2sensor_name(proj_name)[0]

            ds_name = proj_name if proj_name else ''
            ds_name += 'A' if atr_i['ORBIT_DIRECTION'].lower().startswith('asc') else 'D'
            ds_name += f'T{track_num}' if track_num else ''
            ds_name += '_{}'.format(atr_i['DATE12'])

            # assign dataset value
            dsDict[ds_name] = dlos[i]
        dsDict['horizontal'] = dhorz
        dsDict['vertical'] = dvert
        writefile.write(dsDict, out_file=inps.one_outfile, metadata=atr, ref_file=ref_file)

    else:
        print('writing horizontal component to file: '+inps.outfile[0])
        writefile.write(dhorz, out_file=inps.outfile[0], metadata=atr, ref_file=ref_file)
        print('writing vertical   component to file: '+inps.outfile[1])
        writefile.write(dvert, out_file=inps.outfile[1], metadata=atr, ref_file=ref_file)

    return inps.outfile


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    run_asc_desc2horz_vert(inps)

    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
