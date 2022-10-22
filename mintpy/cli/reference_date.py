#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_utils

##################################################################
TEMPLATE = get_template_content('reference_date')

EXAMPLE = """example:
  reference_date.py timeseries.h5 timeseries_ERA5.h5 timeseries_ERA5_demErr.h5 --template smallbaselineApp.cfg
  reference_date.py timeseries_ERA5_demErr.h5 --ref-date 20050107
"""


def create_parser(subparsers=None):
    synopsis = 'Change reference date of time-series HDF5 file.'
    epilog = TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file', nargs='+', help='timeseries file(s)')
    parser.add_argument('-r', '--ref-date', dest='refDate', default='minRMS',
                        help='reference date or method, default: auto. e.g.\n' +
                             '20101120\n' +
                             'time-series HDF5 file with REF_DATE in its attributes\n' +
                             'reference_date.txt - text file with date in YYYYMMDD format in it\n' +
                             'minRMS             - choose date with min residual standard deviation')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('--force', action='store_true',
                        help='Force updating the data matrix.')

    # computing
    parser = arg_utils.add_memory_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: input file type
    ftype = readfile.read_attribute(inps.timeseries_file[0])['FILE_TYPE']
    if ftype != 'timeseries':
        raise ValueError(f'input file type ({ftype}) is NOT timeseries.')

    # check: -t / --template option
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    return inps


def read_template2inps(template_file, inps):
    """Update inps with options from template_file"""
    from mintpy.utils import readfile, utils1 as ut

    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    key = 'mintpy.reference.date'
    if key in template.keys() and template[key]:
        inps.refDate = template[key]

    key = 'mintpy.compute.maxMemory'
    if key in template.keys() and template[key]:
        inps.maxMemory = float(template[key])

    return inps


##################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.reference_date import run_reference_date

    # run
    run_reference_date(inps)


##################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
