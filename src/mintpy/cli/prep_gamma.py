#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

# from mintpy.objects.sensor import SENSOR_NAMES
# copy as below to avoid importing the non-empty mintpy.objects.__init__.py
SENSOR_NAMES = [
    'alos',
    'alos2',
    'csk',
    'env',
    'ers',
    'gf3',
    'jers',
    'ksat5',
    'ni',
    'rs1',
    'rs2',
    'rcm',
    'sen',
    'tsx',
    'uav',
]


##################################################################################################
NOTE = """
  For each interferogram, including unwrapped/wrapped interferograms and coherence, 3 metadata files are required:
  1) reference .par file, e.g. 130118_4rlks.amp.par
  2) secondary .par file, e.g. 130129_4rlks.amp.par
  3) interferogram .off file, e.g. 130118-130129_4rlks.off

  Other metadata files are recommended and can be generated from the above 3 if not existing, more specifically:
  4) baseline files, e.g. 130118-130129_4rlks.baseline and 130118-130129_4rlks.base_perp,
      which can be generated from file 1-3 with Gamma command base_orbit and base_perp.
  5) corner files, e.g. 130118_4rlks.amp.corner_full and 130118_4rlks.amp.corner,
      which can be generated from file 1 with Gamma command SLC_corners.

  This script will read all these files (generate 4 and 5 if not existing), merge them into one, convert their name from
  Gamma style to ROI_PAC style, and write to an metadata file, same name as input binary data file with suffix .rsc,
  e.g. diff_filt_HDR_130118-130129_4rlks.unw.rsc


  For DEM file in radar/geo coordinates (.hgt_sim or .rdc.dem / .utm.dem) and
      lookup table file for geocoding (.UTM_TO_RDC), 2 metadata files are required:
  1) .par      file, for DEM in geo   coordinates and lookup table, e.g.: sim_150911_4rlks.utm.dem.par
  2) .diff_par file, for DEM in radar coordinates, e.g. sim_150911_4rlks.diff_par


  Here is an example of how your Gamma files should look like:
  Before loading:
      For each interferogram, 5 files are needed:
          130118-130129_4rlks.off
          130118_4rlks.amp.par
          130129_4rlks.amp.par
          filt_130118-130129_4rlks.cor
          diff_130118-130129_4rlks.unw
      For each dataset, only one sim* folder with 5 files are needed,
          sim_150911_4rlks.UTM_TO_RDC
          sim_150911_4rlks.diff_par
          sim_150911_4rlks.hgt_sim or sim_150911.rdc.dem
          sim_150911_4rlks.utm.dem
          sim_150911_4rlks.utm.dem.par
  After running prep_gamma.py:
      For each interferogram:
          130118-130129_4rlks.base_perp
          130118-130129_4rlks.baseline
          130118-130129_4rlks.off
          130118_4rlks.ramp.corner
          130118_4rlks.ramp.corner_full
          130118_4rlks.ramp.par
          130129_4rlks.ramp.par
          filt_130118-130129_4rlks.cor
          filt_130118-130129_4rlks.cor.rsc
          diff_130118-130129_4rlks.unw
          diff_130118-130129_4rlks.unw.rsc
      For the geometry files in each dataset:
          sim_150911_4rlks.UTM_TO_RDC
          sim_150911_4rlks.UTM_TO_RDC.rsc
          sim_150911_4rlks.diff_par
          sim_150911_4rlks.rdc.dem      or sim_150911_4rlks.hgt_sim
          sim_150911_4rlks.rdc.dem.rsc  or sim_150911_4rlks.hgt_sim.rsc
          sim_150911_4rlks.utm.dem
          sim_150911_4rlks.utm.dem.par
          sim_150911_4rlks.utm.dem.rsc

  Notes: both - and _ are supported;
         both YYMMDD and YYYYMMDD naming are also supported;
         if no multilooking applied, do not add "_4rlks" in your file names.
"""

EXAMPLE = """example:
  prep_gamma.py  diff_filt_HDR_20130118_20130129_4rlks.unw
  prep_gamma.py  interferograms/*/diff_*rlks.unw --sensor sen
  prep_gamma.py  interferograms/*/filt_*rlks.cor
  prep_gamma.py  interferograms/*/diff_*rlks.int
  prep_gamma.py  sim_20150911_20150922.hgt_sim
  prep_gamma.py  sim_20150911_20150922.utm.dem
  prep_gamma.py  sim_20150911_20150922.UTM_TO_RDC
"""


def create_parser(subparsers=None):
    synopsis = 'Prepare attributes file for Gamma product.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='Gamma file(s)')
    parser.add_argument('--sensor', dest='sensor', type=str, choices=SENSOR_NAMES,
                        help='SAR sensor')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import utils1 as ut

    # check
    inps.file = ut.get_file_list(inps.file, abspath=True)
    inps.file_ext = os.path.splitext(inps.file[0])[1].lower()

    # check: input file extension
    ext_list = ['.unw', '.cor', '.int', '.dem', '.hgt_sim']
    ext_ends = ['to_rdc', '2_rdc', '2rdc']
    if inps.file_ext not in ext_list and not inps.file_ext.endswith(tuple(ext_ends)):
        msg = f'unsupported input file extension: {inps.file_ext}'
        msg += f'\nsupported file extensions: {ext_list + ["*"+x for x in ext_ends]}'
        raise ValueError(msg)

    return inps


##################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_gamma import prep_gamma

    # run
    prep_gamma(inps)


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
