#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2024                               #
############################################################


import glob
import sys
import os
from mintpy.utils.arg_utils import create_argument_parser

#########################################################################
# Default geometry files to extract from static_layers HDF5
GEOMETRY_FILENAMES = [
    'height.tif',
    'los_east.tif',
    'los_north.tif',
    'layover_shadow_mask.tif',
    'local_incidence_angle.tif',
]

EXAMPLE = """example:
  ## Dolphin/ISCE-3 topsStack (auto‑generate metadata)
  prep_isce3.py -f "../../dolphin/unwrapped/*.unw.tif" -b ../baselines -g ../merged/geom/

  ## with existing metadata file
  prep_isce3.py -f "../../dolphin/unwrapped/*.unw.tif" -m ../reference/IW1.xml -g ../merged/geom/

  ## force overwrite existing .rsc files
  prep_isce3.py -f "../../dolphin/unwrapped/*.unw.tif" -g ../merged/geom/ --force
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare ISCE-3 / Dolphin metadata files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # observations
    parser.add_argument('-f', dest='obs_files', type=str, nargs='+',
                        default=['../../dolphin/unwrapped/*.unw.tif'],
                        help='Wildcard path pattern(s) for observation files.\n'
                             'E.g.: unwrapped phase:   ../../dolphin/unwrapped/*.unw.tif\n'
                             '      coherence:        ../../dolphin/interferograms/*.int.cor.tif\n'
                             '      wrapped phase:    ../../dolphin/interferograms/*.int.tif\n'
                             '      connected comp:   ../../dolphin/unwrapped/*.unw.conncomp.tif')

    # metadata
    parser.add_argument('-m', '--meta-file', dest='meta_file', type=str, default=None,
                        help='Metadata file to extract common metadata for the stack.\n'
                             'E.g.: reference burst XML (e.g., reference/IW1.xml). '
                             'If not provided, one will be generated from static_layers.h5.')

    # geometry and baseline
    parser.add_argument('-b', '--baseline-dir', dest='baseline_dir', type=str, default=None,
                        help='Directory with baseline files (e.g., ../baselines). '
                             'If omitted, baseline info will not be added to metadata.')
    parser.add_argument('-g', '--geometry-dir', dest='geom_dir', type=str, required=True,
                        help='Directory containing burst subdirectories with static_layers*.h5 files.\n'
                             'E.g.: ../merged/geom/')
    parser.add_argument('--geom-files', dest='geom_files', type=str, nargs='*',
                        default=GEOMETRY_FILENAMES,
                        help='List of geometry file basenames to extract/merge. Default: %(default)s.')

    # processing flag
    parser.add_argument('--force', dest='update_mode', action='store_false',
                        help='Force to overwrite all .rsc metadata files (disable update mode).')

    parser.add_argument('--out-dir', dest='out_dir', type=str, default=None,
                        help='Output directory for merged geometry files. '
                             'If not provided, defaults to (geometry_dir)/../merged_geom')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.meta_file and '*' in inps.meta_file:
        fnames = glob.glob(inps.meta_file)
        if fnames:
            inps.meta_file = fnames[0]
        else:
            raise FileNotFoundError(inps.meta_file)

    # Expand glob patterns in geometry directory (e.g. "../../t124*/20210104/")
    inps.geom_dirs = [inps.geom_dir]
    if inps.geom_dir and ('*' in inps.geom_dir or '?' in inps.geom_dir):
        matches = sorted(glob.glob(inps.geom_dir))
        if matches:
            inps.geom_dir = matches[0]
            inps.geom_dirs = matches

    # Set default output directory if not provided
    if inps.out_dir is None:
        inps.out_dir = os.path.join(os.path.dirname(inps.geom_dir), 'merged_geom')

    inps.processor = 'tops'
    return inps


#########################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import core function
    from mintpy.prep_isce3 import prep_isce3

    # run
    prep_isce3(inps)


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])