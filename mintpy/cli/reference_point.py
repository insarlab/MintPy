#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

#########################################  Usage  ##############################################
TEMPLATE = get_template_content('reference_point')

NOTE = """note: Reference value cannot be nan, thus, all selected reference point must be:
  a. non zero in mask, if mask is given
  b. non nan  in data (stack)

  Priority:
      input reference_lat/lon
      input reference_y/x
      input selection_method
      existing REF_Y/X attributes (can be ignored by --force option)
      default selection methods:
          maxCoherence
          random

  The recommended reference pixel should meets the following criteria:
  1) not in deforming areas
  2) not in areas affected by strong atmospheric turbulence, such as ionospheric streaks
  3) close but outside of deforming area of interest with similar elevation, to minimize
     the spatial correlation effect of atmosspheric delay, especially for shot-wavelength
     deformation (Chaussard et al., 2013; Morales-Rivera et al., 2016)
  4) in high coherent area to minimize the decorrelation effect
"""

EXAMPLE = """example:
  # for ifgramStack file, update metadata only
  # add --write-data to update data matrix value
  reference_point.py  inputs/ifgramStack.h5  -t smallbaselineApp.cfg  -c avgSpatialCoh.h5
  reference_point.py  inputs/ifgramStack.h5 --method manual
  reference_point.py  inputs/ifgramStack.h5 --method random

  # for all the other files, update both metadata and data matrix value
  reference_point.py  091120_100407.unw -y 257    -x 151      -m Mask.h5
  reference_point.py  geo_velocity.h5   -l 34.45  -L -116.23  -m Mask.h5
"""


def create_parser(subparsers=None):
    synopsis = 'Reference to the same pixel in space.'
    epilog = NOTE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='file to be referenced.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template with reference info')
    parser.add_argument('-m', '--mask', dest='maskFile', help='mask file')

    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='output file name (default: %(default)s). This option is disabled for ifgramStack file.\n'
                             'None (default) for update data value directly without writing to a new file.\n')

    parser.add_argument('--write-data', dest='write_data', action='store_true',
                        help='(option for ifgramStack file only) update data value, in addition to update metadata.')

    parser.add_argument('--reset', action='store_true',
                        help='remove reference pixel information from attributes in the file')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the re-selection of reference point.')

    # coordinates
    coord = parser.add_argument_group('input coordinates')
    coord.add_argument('-y', '--row', dest='ref_y', type=int,
                       help='row/azimuth  number of reference pixel')
    coord.add_argument('-x', '--col', dest='ref_x', type=int,
                       help='column/range number of reference pixel')
    coord.add_argument('-l', '--lat', dest='ref_lat',
                       type=float, help='latitude  of reference pixel')
    coord.add_argument('-L', '--lon', dest='ref_lon',
                       type=float, help='longitude of reference pixel')

    coord.add_argument('-r', '--reference', dest='reference_file',
                       help='use reference/seed info of this file')
    coord.add_argument('--lookup', '--lookup-file', dest='lookup_file',
                       help='Lookup table file from SAR to DEM, i.e. geomap_4rlks.trans\n' +
                            'Needed for radar coord input file with --lat/lon seeding option.')

    # selection method
    parser.add_argument('-c', '--coherence', dest='coherenceFile', default='averageSpatialCoherence.h5',
                        help='use input coherence file to find the pixel with max coherence for reference pixel.')
    parser.add_argument('--min-coherence', dest='minCoherence', type=float, default=0.85,
                        help='minimum coherence of reference pixel for max-coherence method.')
    parser.add_argument('--method', type=str, choices=['maxCoherence', 'manual', 'random'],
                        help='methods to select reference pixel if not given in specific y/x or lat/lon:\n' +
                             'maxCoherence : select pixel with highest coherence value as reference point\n' +
                             '               enabled when there is --coherence option input\n' +
                             'manual       : display stack of input file and manually select reference point\n' +
                             'random       : random select pixel as reference point\n')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']

    # check: turn OFF --output option for ifgramStack file
    if ftype == 'ifgramStack' and inps.outfile:
        inps.outfile = None
        print('WARNING: --outfile is NOT supported for "ifgramStack" file! Ignore it and continue.')

    # check: turn ON --write-data option for non-ifgramStack file
    if ftype != 'ifgramStack' and not inps.write_data:
        inps.write_data = True
        print(f'WARNING: auto turn ON --write-data for file tpye {ftype}')

    return inps


#######################################  Main Function  ########################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.reference_point import read_reference_input, reference_file

    # run
    inps = read_reference_input(inps)

    if inps.go_reference:
        reference_file(inps)


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
