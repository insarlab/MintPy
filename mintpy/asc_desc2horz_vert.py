#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2013-2019, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import argparse
import numpy as np
from mintpy.utils import readfile, writefile, utils as ut


################################################################################
REFERENCE = """reference:
  Wright, T. J., B. E. Parsons, and Z. Lu (2004), Toward mapping 
  surface deformation in three dimensions using InSAR, GRL, 31(1),
"""

EXAMPLE = """example:
  #For asc / desc data with different spatial resolution and coverage
  cd AlosAT424/mintpy
  mask.py velocity.h5 -m maskTempCoh.h5
  geocode.py velocity_msk.h5 -l inputs/geometryRadar.h5 -x 0.00027778 -y -0.00027778 --bbox 32.0 32.5 130.1 130.5

  cd AlosDT73/mintpy
  mask.py velocity.h5 -m maskTempCoh.h5
  geocode.py velocity_msk.h5 -l inputs/geometryRadar.h5 -x 0.00027778 -y -0.00027778 --bbox 32.0 32.5 130.1 130.5

  asc_desc2horz_vert.py  AlosAT424/mintpy/geo_velocity_msk.py  AlosDT73/mintpy/geo_velocity_msk.py

  asc_desc2horz_vert.py  vel_AlosAT424_msk.h5  vel_AlosDT73_msk.h5
  asc_desc2horz_vert.py  vel_EnvAT134_msk.h5   vel_EnvAT256_msk.h5  --azimuth 16
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
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input azimuth angle
    if inps.azimuth < 0.:
        inps.azimuth += 360.
    inps.azimuth *= np.pi/180.
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


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # 1. Extract the common area of two input files
    # Basic info
    atr1 = readfile.read_attribute(inps.file[0])
    atr2 = readfile.read_attribute(inps.file[1])

    # check coordinates
    if any('X_FIRST' not in i for i in [atr1, atr2]):
        raise Exception('Not all input files are geocoded.')
    # check spatial resolution
    if any(atr1[i] != atr2[i] for i in ['X_STEP','Y_STEP']):
        msg = 'file1: {}\n'.format(inps.file[0])
        msg += 'Y_STEP: {} m, X_STEP: {} m\n'.format(atr1['Y_STEP'], atr1['X_STEP'])
        msg += 'file2: {}\n'.format(inps.file[1])
        msg += 'Y_STEP: {} m, X_STEP: {} m'.format(atr2['Y_STEP'], atr2['X_STEP'])
        raise ValueError('input files do not have the same spatial resolution\n{}'.format(msg))

    k1 = atr1['FILE_TYPE']
    print('Input 1st file is '+k1)

    # Common AOI in lalo
    west, east, south, north = get_overlap_lalo(atr1, atr2)
    lon_step = float(atr1['X_STEP'])
    lat_step = float(atr1['Y_STEP'])
    width = int(round((east - west) / lon_step))
    length = int(round((south - north) / lat_step))

    # Read data in common AOI: LOS displacement, heading angle, incident angle
    u_los = np.zeros((2, width*length))
    heading = []
    incidence = []
    for i in range(len(inps.file)):
        fname = inps.file[i]
        print('---------------------')
        print('reading '+fname)
        atr = readfile.read_attribute(fname)

        coord = ut.coordinate(atr)
        [x0, x1] = coord.lalo2yx([west, east], coord_type='lon')
        [y0, y1] = coord.lalo2yx([north, south], coord_type='lat')
        V = readfile.read(fname, box=(x0, y0, x1, y1))[0]
        u_los[i, :] = V.flatten(0)

        heading_angle = float(atr['HEADING'])
        if heading_angle < 0.:
            heading_angle += 360.
        print('heading angle: '+str(heading_angle))
        heading_angle *= np.pi/180.
        heading.append(heading_angle)

        inc_angle = float(ut.incidence_angle(atr, dimension=0))
        inc_angle *= np.pi/180.
        incidence.append(inc_angle)

    # 2. Project displacement from LOS to Horizontal and Vertical components
    # math for 3D: cos(theta)*Uz - cos(alpha)*sin(theta)*Ux + sin(alpha)*sin(theta)*Uy = Ulos
    # math for 2D: cos(theta)*Uv - sin(alpha-az)*sin(theta)*Uh = Ulos   #Uh_perp = 0.0
    # This could be easily modified to support multiple view geometry
    # (e.g. two adjcent tracks from asc & desc) to resolve 3D

    # Design matrix
    A = np.zeros((2, 2))
    for i in range(len(inps.file)):
        A[i, 0] = np.cos(incidence[i])
        A[i, 1] = np.sin(incidence[i]) * np.sin(heading[i]-inps.azimuth)

    A_inv = np.linalg.pinv(A)
    u_vh = np.dot(A_inv, u_los)

    u_v = np.reshape(u_vh[0, :], (length, width))
    u_h = np.reshape(u_vh[1, :], (length, width))

    # 3. Output
    # Attributes
    atr = atr1.copy()
    atr['WIDTH'] = str(width)
    atr['LENGTH'] = str(length)
    atr['X_FIRST'] = str(west)
    atr['Y_FIRST'] = str(north)
    atr['X_STEP'] = str(lon_step)
    atr['Y_STEP'] = str(lat_step)

    print('---------------------')
    horz_file = inps.outfile[0]
    print('writing horizontal component to file: '+horz_file)
    writefile.write(u_h, out_file=horz_file, metadata=atr)

    vert_file = inps.outfile[1]
    print('writing   vertical component to file: '+vert_file)
    writefile.write(u_v, out_file=vert_file, metadata=atr)

    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main()
