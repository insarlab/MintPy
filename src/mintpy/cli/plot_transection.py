#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import argparse
import sys

from mintpy.utils import arg_utils

#####################################################################
# Only one line is supported right now.
GMT_FILE = """GMT xy file, i.e. transect_lonlat.xy:
>
131.1663    33.1157
131.2621    33.0860
"""

EXAMPLE = """example:
  plot_transection.py velocity.h5 --start-yx 5290 5579 --end-yx 12177 482
  plot_transection.py velocity.h5 --start-lalo 30.125 129.988 --end-lalo 30.250 130.116
  plot_transection.py velocity.h5 --line-file  transect_lonlat.xy --dem gsi10m.dem

  # multiple files
  plot_transection.py AlosA*/velocity.h5 AlosD*/velocity.h5 --off 2
  plot_transection.py Kirishima2017*.h5 Kirishima2008*.h5 --off 0 0 10 10
  plot_transection.py Kirishima2017*.h5 Kirishima2008*.h5 --off 0 0 10 10 --lalo0 31.947 130.843 --lalo1 31.947 130.860

  # interactive plot: click two points to draw a profile
"""


def create_parser(subparsers=None):
    synopsis = 'Generate transect/profile along a line'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+',
                        help='input file to show transection')
    parser.add_argument('--dset', dest='dset', help='Dataset name to read')
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                        help='Display limits for matrix plotting.')
    parser.add_argument('--offset','--off', dest='offset', type=float, nargs='+', default=[0.05],
                        help='offset between transects [for multiple files only; default: %(default)s m].\n'
                             'number of input offsets should be:\n'
                             '    1 - same (sequential) offset between adjacent transects OR\n'
                             '    num_file - different (cumulative) offset for each file, starting from 0.')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')

    lines = parser.add_argument_group('Profile location', 'Start/end points of profile')
    lines.add_argument('--start-yx','--yx0', dest='start_yx', metavar=('Y0', 'X0'), type=int, nargs=2,
                       help='start point of the profile in pixel number [y, x]')
    lines.add_argument('--end-yx','--yx1', dest='end_yx', metavar=('Y1', 'X1'), type=int, nargs=2,
                       help='end   point of the profile in pixel number [y, x]')
    lines.add_argument('--start-lalo','--lalo0', dest='start_lalo', metavar=('LAT0', 'LON0'), type=float, nargs=2,
                       help='start point of the profile in [lat, lon]')
    lines.add_argument('--end-lalo','--lalo1', dest='end_lalo', metavar=('LAT1', 'LON1'), type=float, nargs=2,
                       help='end   point of the profile in [lat, lon]')
    lines.add_argument('--line-file', dest='lola_file',
                       help='file with start and end point info in lon lat, same as GMT format.\n'+GMT_FILE)

    lines.add_argument('--interpolation', default='nearest', choices=['nearest', 'bilinear', 'cubic'],
                       help='interpolation method while extacting profile along the line. Default: nearest.')
    lines.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=2.0,
                       help='Point marker size. Default: 2.0')

    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_save_argument(parser)
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    import numpy as np

    from mintpy.utils import readfile, utils as ut

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check: coupled options
    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True

    # check: --offset option
    inps.file = ut.get_file_list(inps.file)
    inps.atr = readfile.read_attribute(inps.file[0])
    inps.coord = ut.coordinate(inps.atr)
    inps.num_file = len(inps.file)

    if inps.num_file > 1:
        # a) one input: it's interval between adjacent files
        if len(inps.offset) == 1:
            inps.offset = np.ones(inps.num_file, dtype=np.float32) * inps.offset
            inps.offset[0] = 0.
            inps.offset = np.cumsum(inps.offset)

        # b) multiple input: it's exact offset of all files
        elif len(inps.offset) == inps.num_file:
            inps.offset = np.array(inps.offset, dtype=np.float32)

        # c) do not support any other numbers of inputs
        else:
            msg = f'input number of offsets: {len(inps.offset)}.'
            msg += f'\nIt should be 1 or number of files: {inps.num_file}'
            raise ValueError(msg)
    else:
        # disable offset for single input file
        inps.offset = np.array([0], dtype=np.float32)

    # check: --line-file option (lola_file --> start/end_lalo)
    if inps.lola_file:
        inps.start_lalo, inps.end_lalo = read_lonlat_file(inps.lola_file)

    # check: --start/end-lalo (start/end_lalo --> start/end_yx)
    if inps.start_lalo and inps.end_lalo:
        [y0, y1] = inps.coord.lalo2yx([inps.start_lalo[0], inps.end_lalo[0]], coord_type='lat')
        [x0, x1] = inps.coord.lalo2yx([inps.start_lalo[1], inps.end_lalo[1]], coord_type='lon')
        inps.start_yx = [y0, x0]
        inps.end_yx = [y1, x1]

    # default: --dset option
    if not inps.dset:
        inps.dset = readfile.get_slice_list(inps.file[0])[0]

    return inps


#####################################################################
def read_lonlat_file(lonlat_file):
    """Read Start/End lat/lon from lonlat text file in gmt format.

    Parameters: lonlat_file    - str, text file in gmt lonlat point file
    Returns:    start/end_lalo - list of 2 float
    """
    fll = open(lonlat_file)
    lines = fll.read().splitlines()
    [lon0, lat0] = [float(i) for i in lines[1].split()]
    [lon1, lat1] = [float(i) for i in lines[2].split()]
    fll.close()

    start_lalo = [lat0, lon0]
    end_lalo = [lat1, lon1]
    return start_lalo, end_lalo


def get_view_cmd(iargs):
    """Assemble view.py command line from input arguments"""

    # define ALL parsing options from create_parser() that are common to view.py
    parser = argparse.ArgumentParser(description='view.py parser')
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                        help='Display limits for matrix plotting.')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')
    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_save_argument(parser)

    # get args that are applicable to view.py
    unique_args = parser.parse_known_args(iargs)[1]
    view_args = [i for i in iargs if i not in unique_args] if iargs else []

    # assemble view.py command line
    inps = cmd_line_parse(iargs)
    view_cmd = f'view.py {inps.file[0]} '
    view_cmd += f' {inps.dset} ' if inps.dset else ''
    view_cmd += ' '.join(view_args)

    return view_cmd


############################ Main ###################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.plot_transection import transectionViewer

    # run
    view_cmd = get_view_cmd(iargs)
    obj = transectionViewer(inps, view_cmd)
    obj.open()
    obj.plot()
    obj.fig.canvas.mpl_disconnect(obj.cid)


#####################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
