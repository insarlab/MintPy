#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check
    atr1 = readfile.read_attribute(inps.file[0])
    atr2 = readfile.read_attribute(inps.file[1])

    # check: if input file type is supported
    ts_ftypes = ['timeseries', 'HDFEOS']
    ftype1, ftype2 = atr1['FILE_TYPE'], atr2['FILE_TYPE']
    if any(x in ts_ftypes for x in [ftype1, ftype2]):
        raise Exception(f'input file types ({ftype1}, {ftype2}) contains UN-supported file types: {ts_ftypes}!')

    # check: if input is in geo-coordinates
    if any('X_FIRST' not in i for i in [atr1, atr2]):
        raise Exception('Not all input files are geocoded.')

    # check: if input spatial resolutions are consistent
    if any(atr1[i] != atr2[i] for i in ['X_STEP','Y_STEP']):
        msg  = '\tfile1: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[0], atr1['Y_STEP'], atr1['X_STEP'], atr1.get('X_UNIT', 'degrees'))
        msg += '\tfile2: {}, Y/X_STEP: {} / {} {}\n'.format(inps.file[1], atr2['Y_STEP'], atr2['X_STEP'], atr2.get('X_UNIT', 'degrees'))
        msg += '\tRe-run geocode.py --lat-step --lon-step to make them consistent.'
        raise ValueError(f'input files do NOT have the same spatial resolution\n{msg}')

    # check: if input reference points are consistent
    ref_lat1, ref_lon1 = (float(atr1[i]) for i in ['REF_LAT', 'REF_LON'])
    ref_lat2, ref_lon2 = (float(atr2[i]) for i in ['REF_LAT', 'REF_LON'])
    ref_y_diff = abs((ref_lat1 - ref_lat2) / float(atr1['Y_STEP']))
    ref_x_diff = abs((ref_lon1 - ref_lon2) / float(atr1['X_STEP']))
    if any(ref_diff > inps.max_ref_yx_diff for ref_diff in [ref_y_diff, ref_x_diff]):
        msg = f'REF_LAT/LON difference between input files > {inps.max_ref_yx_diff} pixels!\n'
        for fname, ref_lat, ref_lon in zip(inps.file, [ref_lat1, ref_lat2], [ref_lon1, ref_lon2]):
            msg += f'file: {fname}\n'
            msg += f'\tREF_LAT/LON: [{ref_lat:.8f}, {ref_lon:.8f}]\n'
        raise ValueError(msg)

    return inps


################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.asc_desc2horz_vert import run_asc_desc2horz_vert

    # run
    run_asc_desc2horz_vert(inps)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
