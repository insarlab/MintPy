#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import sys
import argparse
import numpy as np
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile, utils as ut


################################################################################
REFERENCE = """reference:
  Wright, T. J., B. E. Parsons, and Z. Lu (2004), Toward mapping surface deformation
    in three dimensions using InSAR, Geophysical Research Letters, 31(1), n/a-n/a,
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
  asc_desc2horz_vert.py Alos2AT131/mintpy/20171219_20190702.unw Alos2DT23/mintpy/20171211_20190819.unw --output-one Kirishima2017post.h5
  view.py Kirishima2017post.h5 -u cm --wrap --wrap-range -5 5  #check deformation signal with multiple viewing geometries.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Project Asc and Desc LOS displacement to Horizontal and Vertical direction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('file', nargs=2,
                        help='Ascending and descending files\n'
                             'Both files need to be geocoded in the same spatial resolution.')
    parser.add_argument('-d', '--dset', dest='dsname', type=str, help='dataset to use, default: 1st dataset')
    parser.add_argument('--azimuth', '--az', dest='azimuth', type=float, default=-90.0,
                        help='Azimuth angle in degrees (anti-clockwise) of the direction of the horizontal movement (default: %(default)s).\n'
                             'E.g.: -90 for East direction\n'
                             '      0   for North direction\n'
                             'Set to the azimuth angle of the strike-slip fault to measure the fault-parallel displacement.\n'
                             'Note:\n'
                             'a. This assumes no deformation in its perpendicular direction\n'
                             'b. Near north direction can not be well resolved due to the lack of\n'
                             '   diversity in viewing geometry. Check exact dilution of precision for \n'
                             '   each component in Wright et al. (2004, GRL)')
    parser.add_argument('--max-ref-yx-diff', dest='max_ref_yx_diff', type=int, default=3,
                        help='Maximum difference between REF_Y/X (derived from REF_LAT/LON) of input files '+
                             '(default: %(default)s).')

    # output
    parser.add_argument('-o', '--output', dest='outfile', nargs=2, metavar=('HZ_FILE','UP_FILE'), default=['hz.h5', 'up.h5'],
                        help='output file name for vertical and horizontal components')
    parser.add_argument('--one-output','--oo', dest='one_outfile',
                        help='Stack the input/output files into one HDF5 file.\n' +
                             'This will disable the HZ/UP_FILE output option.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr1 = readfile.read_attribute(inps.file[0])
    atr2 = readfile.read_attribute(inps.file[1])

    # check coordinates
    if any('X_FIRST' not in i for i in [atr1, atr2]):
        raise Exception('Not all input files are geocoded.')

    # check spatial resolution
    if any(atr1[i] != atr2[i] for i in ['X_STEP','Y_STEP']):
        msg  = '\tfile1: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[0], atr1['Y_STEP'], atr1['X_STEP'], atr1.get('X_UNIT', 'degrees'))
        msg += '\tfile2: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[1], atr2['Y_STEP'], atr2['X_STEP'], atr2.get('X_UNIT', 'degrees'))
        msg += '\tRe-run geocode.py --lat-step --lon-step to make them consistent.'
        raise ValueError('input files do NOT have the same spatial resolution\n{}'.format(msg))

    # check reference point
    coord1 = ut.coordinate(atr1)
    coord2 = ut.coordinate(atr2)
    ref_lalo1 = [float(atr1[i]) for i in ['REF_LAT', 'REF_LON']]
    ref_lalo2 = [float(atr2[i]) for i in ['REF_LAT', 'REF_LON']]
    ref_yx1 = coord1.geo2radar(ref_lalo1[0], ref_lalo1[1])[:2]
    ref_yx2 = coord2.geo2radar(ref_lalo2[0], ref_lalo2[1])[:2]

    if any(abs(i1 - i2) > inps.max_ref_yx_diff for i1,i2 in zip(ref_yx1, ref_yx2)):
        msg = 'REF_Y/X difference between input files > {}!\n'.format(inps.max_ref_yx_diff)
        for fname, ref_lalo, ref_yx in zip(inps.file, [ref_lalo1, ref_lalo2], [ref_yx1, ref_yx2]):
            msg += 'file1: {}\n'.format(fname)
            msg += '\tREF_LAT/LON: {}\n'.format(ref_lalo)
            msg += '\tREF_Y/X: {}\n'.format(ref_yx)
        raise ValueError(msg)

    # use ref_file for time-series file writing
    if atr1['FILE_TYPE'] == 'timeseries':
        inps.ref_file = inps.file[0]
    else:
        inps.ref_file = None
    return inps


################################################################################
def get_overlap_lalo(atr1, atr2):
    """Find overlap area in lat/lon of two geocoded files
    Parameters: atr1/2  - dict, attribute dictionary of two input files in geo coord
    Returns:    W/E/S/N - float, West/East/South/North in deg
    """
    W1, E1, S1, N1 = ut.four_corners(atr1)
    W2, E2, S2, N2 = ut.four_corners(atr2)

    west = max(W1, W2)
    east = min(E1, E2)
    north = min(N1, N2)
    south = max(S1, S2)

    return west, east, south, north


def get_design_matrix(atr1, atr2, azimuth=-90):
    """Get the design matrix A to convert asc/desc to hz/up.
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

    Parameters: atr1/2   - dict, metadata of input LOS files
                azimuth  - float, azimuth angle for the horizontal direction of interest in degrees.
                           Measured from the north with anti-clockwise direction as positive.
    Returns:    A        - 2D matrix in size of (2, 2)

    """
    # degree to radian
    azimuth *= np.pi / 180.

    atr_list = [atr1, atr2]
    A = np.zeros((2, 2))
    for i in range(len(atr_list)):
        atr = atr_list[i]

        # LOS incidence angle
        los_inc_angle = float(ut.incidence_angle(atr, dimension=0, print_msg=False))
        print('LOS incidence angle: {} deg'.format(los_inc_angle))
        los_inc_angle *= np.pi / 180.

        # LOS azimuth angle
        los_az_angle = ut.heading2azimuth_angle(float(atr['HEADING']))
        print('LOS azimuth angle: {} deg'.format(los_az_angle))
        los_az_angle *= np.pi / 180.

        # construct design matrix
        A[i, 0] = np.cos(los_inc_angle)
        A[i, 1] = np.sin(los_inc_angle) * np.cos(los_az_angle - azimuth)

    return A


def asc_desc2horz_vert(data_asc, data_desc, atr_asc, atr_desc, azimuth=-90):
    """Decompose asc / desc LOS data into horz / vert data.
    Parameters: data_asc/desc - 2D np.ndarray, displacement in LOS
                atr_asc/desc  - dict, metadata
                azimuth       - float, azimuth angle for the horizontal direction of interest in degrees.
                                Measured from the north with anti-clockwise direction as positive.
    Returns:    data_h/v      - 2D np.ndarray, displacement in horizontal / vertical
    """
    length, width = data_asc.shape
    # prepare LOS data
    data_los = np.vstack((data_asc.flatten(), data_desc.flatten()))

    # get design matrix
    print('get design matrix')
    A = get_design_matrix(atr_asc, atr_desc, azimuth=azimuth)

    # decompose
    print('project asc/desc into horz/vert direction')
    data_vh = np.dot(np.linalg.pinv(A), data_los).astype(np.float32)
    data_v = np.reshape(data_vh[0, :], (length, width))
    data_h = np.reshape(data_vh[1, :], (length, width))

    return data_h, data_v


def asc_desc_files2horz_vert(fname1, fname2, dsname=None, azimuth=-90):
    """Decompose asc / desc LOS files into horz / vert data.
    Parameters: fname1/2  - str, LOS data
                dsname    - str, dataset name
                azimuth   - float, azimuth angle for the horizontal direction of interest in degrees.
                            Measured from the north with anti-clockwise direction as positive.
    Returns:    dH/dV     - 2D matrix
                atr       - dict, metadata with updated size and resolution.
                dLOS_list - list of 2D matrices
                atr_list  - list of dict
    """
    print('---------------------')
    fnames = [fname1, fname2]
    # 1. Extract the common area of two input files
    # Basic info
    atr_list = []
    for fname in fnames:
        atr_list.append(readfile.read_attribute(fname, datasetName=dsname))

    # Common AOI in lalo
    west, east, south, north = get_overlap_lalo(atr_list[0], atr_list[1])
    lon_step = float(atr_list[0]['X_STEP'])
    lat_step = float(atr_list[0]['Y_STEP'])
    width = int(round((east - west) / lon_step))
    length = int(round((south - north) / lat_step))
    print('common area in SNWE: {}'.format((south, north, west, east)))

    # 2. Read LOS data in common AOI
    dLOS_list = []
    for i in range(len(fnames)):
        fname = fnames[i]
        atr = readfile.read_attribute(fname, datasetName=dsname)

        # get box2read for the current file
        coord = ut.coordinate(atr)
        x0 = coord.lalo2yx(west, coord_type='lon')
        y0 = coord.lalo2yx(north, coord_type='lat')
        box = (x0, y0, x0 + width, y0 + length)

        # read
        dLOS_list.append(readfile.read(fname, box=box, datasetName=dsname)[0])

        # msg
        msg = 'read '
        if dsname:
            msg += '{} '.format(dsname)
        msg += 'from file: {}'.format(fname)
        print(msg)

    # 3. Project displacement from LOS to Horizontal and Vertical components
    print('---------------------')
    dH, dV = asc_desc2horz_vert(dLOS_list[0],
                                dLOS_list[1],
                                atr_list[0],
                                atr_list[1],
                                azimuth=azimuth)

    # 4. Update Attributes
    atr = atr_list[0].copy()
    if dsname and atr['FILE_TYPE'] in ['ifgramStack', 'timeseries', 'HDFEOS']:
        atr['FILE_TYPE'] = 'displacement'

    atr['WIDTH']  = str(width)
    atr['LENGTH'] = str(length)
    atr['X_STEP'] = str(lon_step)
    atr['Y_STEP'] = str(lat_step)
    atr['X_FIRST'] = str(west)
    atr['Y_FIRST'] = str(north)

    # update REF_X/Y
    ref_lat, ref_lon = float(atr['REF_LAT']), float(atr['REF_LON'])
    coord = ut.coordinate(atr)
    [ref_y, ref_x] = coord.geo2radar(ref_lat, ref_lon)[0:2]
    atr['REF_Y'] = int(ref_y)
    atr['REF_X'] = int(ref_x)

    return dH, dV, atr, dLOS_list, atr_list


def write_to_one_file(outfile, dH, dV, atr, dLOS_list, atr_list, ref_file=None):
    """Write all datasets into one HDF5 file"""

    print('write all datasets into {}'.format(outfile))
    dsDict = {}
    for i in range(len(atr_list)):
        # auto dataset name
        atr_i = atr_list[i]
        dsName = sensor.project_name2sensor_name(atr_i['PROJECT_NAME'])[0]
        if atr['ORBIT_DIRECTION'].lower().startswith('asc'):
            dsName += 'A'
        else:
            dsName += 'D'
        if 'trackNumber' in atr_i.keys():
            dsName += 'T{}'.format(atr_i['trackNumber'])
        dsName += '_{}'.format(atr_i['DATE12'])

        dsDict[dsName] = dLOS_list[i]
    dsDict['vertical'] = dV
    dsDict['horizontal'] = dH

    writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=ref_file)
    return outfile


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    (dH, dV, atr,
     dLOS_list, atr_list) = asc_desc_files2horz_vert(inps.file[0],
                                                     inps.file[1],
                                                     dsname=inps.dsname,
                                                     azimuth=inps.azimuth)

    print('---------------------')
    if inps.one_outfile:
        write_to_one_file(inps.one_outfile,
                          dH, dV, atr,
                          dLOS_list, atr_list,
                          ref_file=inps.ref_file)

    else:
        print('writing horizontal component to file: '+inps.outfile[0])
        writefile.write(dH, out_file=inps.outfile[0], metadata=atr, ref_file=inps.ref_file)

        print('writing   vertical component to file: '+inps.outfile[1])
        writefile.write(dV, out_file=inps.outfile[1], metadata=atr, ref_file=inps.ref_file)

    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
