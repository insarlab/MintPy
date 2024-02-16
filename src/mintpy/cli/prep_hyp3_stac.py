#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Forrest Williams, Aug 2022    #
############################################################


import re
import sys
from pathlib import Path

from mintpy.utils import readfile
from mintpy.utils import utils1 as ut
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################
NOTE = """
  prep_hyp3_stac is designed to work with STAC collections of HyP3 products generated using the HyP3 SDK.
  These collections can be created using the following syntax:

  from hyp3_sdk import HyP3, stac
  hyp3 = HyP3(prompt=True)
  batch = hyp3.find_jobs(name='PROJECT_NAME_HERE', status_code='SUCCEEDED')
  stac.create_stac_collection(batch=batch, out_path='stac')


  Notes:
    HyP3 currently only supports generation of Sentinel-1 interferograms, so
    some Sentinel-1 metadata is hard-coded. If HyP3 adds processing of interferograms
    from other satellites, changes will be needed.
"""

EXAMPLE = """example:
  prep_hyp3_stac.py -t smallbaselineApp.cfg    # recommended
  prep_hyp3_stac.py stac/collection.json
"""


def create_parser(subparsers=None):
    synopsis = 'Create input ifgramStack.h5 and geometryGeo.h5 files using HyP3 STAC datasets.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis + NOTE, epilog=epilog, subparsers=subparsers
    )

    parser.add_argument('-t', '--template', dest='template_file', type=str, help='template file with the options')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        nargs=2,
        dest='outfile',
        default=['./inputs/ifgramStack.h5', './inputs/geometryGeo.h5'],
        help='output ifgramStack and geometryGeo HDF5 files',
    )
    parser.add_argument(
        '--yx',
        type=str,
        dest='yx',
        default=None,
        help='Subset data to bounding box in format "y0 y1 x0 x1" using array coordinates. Default is the full frame.',
    )
    parser.add_argument(
        '--lalo',
        type=str,
        dest='lalo',
        default=None,
        help='Subset data to bounding box in format "S N W E" using geographic coordinates. Default is the full frame.',
    )
    parser.add_argument(
        '--update',
        dest='updateMode',
        action='store_true',
        help='Enable the update mode: checking dataset already loaded.',
    )
    parser.add_argument(
        '--compression', choices={'gzip', 'lzf', None}, default=None, help='HDF5 file compression, default: %(default)s'
    )
    return parser


def template_subset_to_flat_string(subset_str):
    y_range, x_range = subset_str.split(',')
    y_range1, y_range2 = [int(i) for i in y_range.split(':')]
    x_range1, x_range2 = [int(i) for i in x_range.split(':')]
    return f'{y_range1} {y_range2} {x_range1} {x_range2}'


def read_template2inps(template_file, inps):
    """Read input template file into inps"""
    print(f'read options from template file: {Path(template_file).name}')
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)
    for key in list(template.keys()):
        if template[key] == 'auto':
            template.pop(key)

    to_load = {}
    prefix = 'mintpy.load.'
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        to_load[key] = template[prefix + key]

    prefix = 'mintpy.subset.'
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        to_load[key] = template_subset_to_flat_string(template[prefix + key])

    iDict = vars(inps)
    for key, value in to_load.items():
        if key not in iDict:
            iDict[key] = value
        elif iDict[key] is None:
            iDict[key] = value

    return inps


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    if inps.yx:
        inps.yx = [int(i) for i in inps.yx.split(' ')]

    if inps.lalo:
        inps.lalo = [float(i) for i in inps.lalo.split(' ')]

    return inps


#########################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_hyp3_stac import prep_hyp3_stac

    # run
    prep_hyp3_stac(inps)


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
