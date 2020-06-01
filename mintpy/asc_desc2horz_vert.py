#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import sys
import argparse
import numpy as np
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

  # write all asc/desc/horz/vert datasets into one file
  asc_desc2horz_vert.py Alos2AT131/mintpy/20171219_20190702.unw Alos2DT23/mintpy/20171211_20190819.unw --output-one Kirishima2017post.h5
  view.py Kirishima2017post.h5 -u cm --wrap --wrap-range -5 5  #check deformation signal with multiple viewing geometries.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Project Asc and Desc LOS displacement to Horizontal and Vertical direction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('file', nargs=2,
                        help='ascending and descending files\n' +
                             'Both files need to be geocoded in the same spatial resolution.')
    parser.add_argument('--azimuth', '--az', dest='azimuth', type=float, default=90.0,
                        help='azimuth angle in degree (clockwise) of the direction of the horizontal movement\n' +
                             'default is 90.0 for E-W component, assuming no N-S displacement.\n' +
                             'i.e. azimuth angle of strike-slip fault\n\n' +
                             'Note:\n' +
                             'a. This assumes no deformation in its perpendicular direction\n' +
                             'b. Near north direction can not be well resolved due to the lack of\n' +
                             '   diversity in viewing geometry. Check exact dilution of precision for \n' +
                             '   each component in Wright et al., 2004, GRL')
    parser.add_argument('-o', '--output', dest='outfile', nargs=2, metavar=('HZ_FILE','UP_FILE'), default=['hz.h5', 'up.h5'],
                        help='output file name for vertical and horizontal components')
    parser.add_argument('--one-output','--oo', dest='one_outfile',
                        help='Stack the input/output files into one HDF5 file.\n' +
                             'This will disable the HZ/UP_FILE output option.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input azimuth angle
    if inps.azimuth < 0.:
        inps.azimuth += 360.
    inps.azimuth *= np.pi/180.

    atr1 = readfile.read_attribute(inps.file[0])
    atr2 = readfile.read_attribute(inps.file[1])

    # check coordinates
    if any('X_FIRST' not in i for i in [atr1, atr2]):
        raise Exception('Not all input files are geocoded.')

    # check spatial resolution
    if any(atr1[i] != atr2[i] for i in ['X_STEP','Y_STEP']):
        msg  = '\tfile1: {}, Y/X_STEP: {} / {} m\n'.format(inps.file[0], atr1['Y_STEP'], atr1['X_STEP'])
        msg += '\tfile2: {}, Y/X_STEP: {} / {} m\n'.format(inps.file[1], atr2['Y_STEP'], atr2['X_STEP'])
        msg += '\tRe-run geocode.py --lat-step --lon-step to make them consistent.'
        raise ValueError('input files do not have the same spatial resolution\n{}'.format(msg))

    # check reference point
    # round to 3 decimal digits (~100 m)
    ref_lalo_1 = ['{:.3f}'.format(float(atr1[i])) for i in ['REF_LAT','REF_LON']]
    ref_lalo_2 = ['{:.3f}'.format(float(atr2[i])) for i in ['REF_LAT','REF_LON']]
    if ref_lalo_1 != ref_lalo_2:
        msg  = '\tfile1: {}, REF_LAT/LON: {} {}\n'.format(inps.file[0], ref_lalo_1, atr1.get('X_UNIT', 'deg'))
        msg += '\tfile2: {}, REF_LAT/LON: {} {}\n'.format(inps.file[1], ref_lalo_2, atr2.get('X_UNIT', 'deg'))
        msg += '\tRe-run reference_point.py --lat --lon to make them consistent.'
        raise ValueError('input files do not have the same reference point from REF_LAT/LON values!\n{}'.format(msg))

    # use ref_file for time-series file writing
    if atr1['FILE_TYPE'] == 'timeseries':
        inps.ref_file = inps.file[0]
    else:
        inps.ref_file = None
    return inps


################################################################################
def get_overlap_lalo(atr1, atr2):
    """Find overlap area in lat/lon of two geocoded files
    Inputs:
        atr1/2 - dict, attribute dictionary of two input files in geo coord
    Outputs:
        W/E/S/N - float, West/East/South/North in deg
    """
    W1, E1, S1, N1 = ut.four_corners(atr1)
    W2, E2, S2, N2 = ut.four_corners(atr2)

    west = max(W1, W2)
    east = min(E1, E2)
    north = min(N1, N2)
    south = max(S1, S2)

    return west, east, south, north


def get_design_matrix(atr1, atr2, az_angle=90):
    """Get the design matrix A to convert asc/desc to hz/up.
    Only asc + desc -> hz + up is implemented for now.

    Project displacement from LOS to Horizontal and Vertical components
        math for 3D: cos(theta)*Uz - cos(alpha)*sin(theta)*Ux + sin(alpha)*sin(theta)*Uy = Ulos
        math for 2D: cos(theta)*Uv - sin(alpha-az)*sin(theta)*Uh = Ulos   #Uh_perp = 0.0
    This could be easily modified to support multiple view geometry
        (e.g. two adjcent tracks from asc & desc) to resolve 3D

    Parameters: atr1/2 : dict, metadata of input LOS files
    Returns:    A : 2D matrix in size of (2,2)

    """
    atr_list = [atr1, atr2]
    A = np.zeros((2, 2))
    for i in range(len(atr_list)):
        atr = atr_list[i]

        # incidence angle
        inc_angle = float(ut.incidence_angle(atr, dimension=0, print_msg=False))
        print('incidence angle: '+str(inc_angle))
        inc_angle *= np.pi/180.

        # heading angle
        head_angle = float(atr['HEADING'])
        if head_angle < 0.:
            head_angle += 360.
        print('heading angle: '+str(head_angle))
        head_angle *= np.pi/180.

        # construct design matrix
        A[i, 0] = np.cos(inc_angle)
        A[i, 1] = np.sin(inc_angle) * np.sin(head_angle - az_angle)
    return A


def asc_desc2horz_vert(fname1, fname2):
    """Decompose asc / desc LOS data into horz / vert data.
    Parameters: fname1/2 : str, LOS data
    Returns:    dH/dV    : 2D matrix
                atr      : dict, metadata with updated size and resolution.
    """
    fnames = [fname1, fname2]
    # 1. Extract the common area of two input files
    # Basic info
    atr_list = []
    for fname in fnames:
        atr_list.append(readfile.read_attribute(fname))

    # Common AOI in lalo
    west, east, south, north = get_overlap_lalo(atr_list[0], atr_list[1])
    lon_step = float(atr_list[0]['X_STEP'])
    lat_step = float(atr_list[0]['Y_STEP'])
    width = int(round((east - west) / lon_step))
    length = int(round((south - north) / lat_step))

    # 2. Read data in common AOI: LOS displacement, heading angle, incident angle
    print('---------------------')
    dLOS = np.zeros((2, width*length), dtype=np.float32)
    for i in range(len(fnames)):
        fname = fnames[i]
        print('reading '+fname)
        atr = readfile.read_attribute(fname)

        # get box2read for the current file
        coord = ut.coordinate(atr)
        [x0, x1] = coord.lalo2yx([west, east], coord_type='lon')
        [y0, y1] = coord.lalo2yx([north, south], coord_type='lat')
        dLOS[i, :] = readfile.read(fname, box=(x0, y0, x0 + width, y0 + length))[0].flatten()

    # 3. Project displacement from LOS to Horizontal and Vertical components
    print('---------------------')
    print('get design matrix')
    A = get_design_matrix(atr_list[0], atr_list[1])
    print('project asc/desc into horz/vert direction')
    dVH = np.dot(np.linalg.pinv(A), dLOS).astype(np.float32)
    dV = np.reshape(dVH[0, :], (length, width))
    dH = np.reshape(dVH[1, :], (length, width))

    # 4. Update Attributes
    atr = atr_list[0].copy()
    atr['WIDTH'] = str(width)
    atr['LENGTH'] = str(length)
    atr['X_FIRST'] = str(west)
    atr['Y_FIRST'] = str(north)
    atr['X_STEP'] = str(lon_step)
    atr['Y_STEP'] = str(lat_step)

    # update REF_X/Y
    ref_lat, ref_lon = float(atr['REF_LAT']), float(atr['REF_LON'])
    coord = ut.coordinate(atr)
    [ref_y, ref_x] = coord.geo2radar(ref_lat, ref_lon)[0:2]
    atr['REF_Y'] = int(ref_y)
    atr['REF_X'] = int(ref_x)

    return dH, dV, atr, dLOS, atr_list


def write_to_one_file(outfile, dH, dV, atr, dLOS, atr_list, ref_file=None):
    """Write all datasets into one HDF5 file"""
    from mintpy.objects import sensor

    print('write all datasets into {}'.format(outfile))
    length, width = dH.shape
    dsDict = {}
    for i in range(len(atr_list)):
        # auto dataset name
        atr = atr_list[i]
        dsName = sensor.project_name2sensor_name(atr['FILE_PATH'])[0]
        if atr['ORBIT_DIRECTION'].lower().startswith('asc'):
            dsName += 'A'
        else:
            dsName += 'D'
        if 'trackNumber' in atr.keys():
            dsName += 'T{}'.format(atr['trackNumber'])
        dsName += '_{}'.format(atr['DATE12'])

        dsDict[dsName] = dLOS[i,:].reshape(length, width)
    dsDict['vertical'] = dV
    dsDict['horizontal'] = dH

    writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=ref_file)
    return outfile


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    dH, dV, atr, dLOS, atr_list = asc_desc2horz_vert(inps.file[0], inps.file[1])

    print('---------------------')
    if inps.one_outfile:
        write_to_one_file(inps.one_outfile, dH, dV, atr, dLOS, atr_list, ref_file=inps.ref_file)
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
