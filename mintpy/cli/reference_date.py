############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
import time

from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_utils


##################################################################
TEMPLATE = get_template_content('reference_date')

EXAMPLE = """example:
  reference_date.py timeseries.h5 timeseries_ERA5.h5 timeseries_ERA5_demErr.h5 --template smallbaselineApp.cfg
  reference_date.py timeseries_ERA5_demErr.h5 --ref-date 20050107
"""


def create_parser(subparsers=None):
    synopsis = 'Change reference date of timeseries.'
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
    from mintpy.utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.timeseries_file[0])
    if 'timeseries' not in atr['FILE_TYPE'].lower():
        raise ValueError('input file type: {} is not timeseries.'.format(atr['FILE_TYPE']))
    return inps


##################################################################
def main(iargs=None):
    from mintpy.reference_date import change_timeseries_ref_date, read_ref_date, read_template2inps

    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # read reference date
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    inps.refDate = read_ref_date(inps)

    # run referencing in time
    if inps.refDate:
        for ts_file in inps.timeseries_file:
            change_timeseries_ref_date(ts_file,
                                       ref_date=inps.refDate,
                                       outfile=inps.outfile,
                                       max_memory=inps.maxMemory,
                                       force=inps.force)

            #to distinguish the modification time of input files
            time.sleep(1)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))


##################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
