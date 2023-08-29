#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils import arg_utils

###########################################################################################
EXAMPLE = """example:
  tsview.py timeseries.h5
  tsview.py timeseries.h5  --wrap
  tsview.py timeseries.h5  --yx 300 400 --zero-first  --nodisplay
  tsview.py geo_timeseries.h5  --lalo 33.250 131.665  --nodisplay
  tsview.py slcStack.h5 -u dB -v 20 60 -c gray

  # press left / right key to slide images

  # multiple time-series files
  tsview.py timeseries_ERA5_ramp_demErr.h5 timeseries_ERA5_ramp.h5 timeseries_ERA5.h5 timeseries.h5 --off 5
  tsview.py timeseries_ERA5_ramp_demErr.h5 ../GIANT/Stack/LS-PARAMS.h5 --off 5 --label mintpy giant
"""


def create_parser(subparsers=None):
    synopsis = 'Interactive time-series viewer'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+',
                        help='time-series file to display, e.g.:\n'
                             'timeseries_ERA5_ramp_demErr.h5 (MintPy)\n'
                             'S1_IW12_128_0593_0597_20141213_20180619.he5 (HDF-EOS5)')
    parser.add_argument('--label', dest='file_label', nargs='*',
                        help='labels to display for multiple input files')
    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                        help='Y limits for point plotting.')
    parser.add_argument('--tick-right', dest='tick_right', action='store_true',
                        help='set tick and tick label to the right')
    parser.add_argument('-l','--lookup', dest='lookup_file', type=str,
                        help='lookup table file')
    parser.add_argument('--no-show-img','--not-show-image', dest='disp_fig_img', action='store_false',
                        help='do NOT show the map figure.\n'
                             'Useful for plotting a point time series only.\n'
                             'This option requires --yx/lalo input.')

    parser.add_argument('-n', dest='idx', metavar='NUM', type=int,
                        help='Epoch/slice number for initial display.')
    parser.add_argument('--error', dest='error_file',
                        help='txt file with error for each date.')

    # time info
    parser.add_argument('--start-date', dest='start_date', type=str,
                        help='start date of displacement to display')
    parser.add_argument('--end-date', dest='end_date', type=str,
                        help='end date of displacement to display')
    parser.add_argument('--exclude', '--ex', dest='ex_date_list', nargs='*', default=['exclude_date.txt'],
                        help='Exclude date shown as gray.')
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                        help='Set displacement at first acquisition to zero.')
    parser.add_argument('--off','--offset', dest='offset', type=float,
                        help='Offset for each timeseries file.')

    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')

    # temporal model fitting
    parser.add_argument('--nomodel', '--nofit', dest='plot_model', action='store_false',
                        help='Do not plot the prediction of the time function (deformation model) fitting.')
    parser.add_argument('--plot-model-conf-int', '--plot-fit-conf-int',
                        dest='plot_model_conf_int', action='store_true',
                        help='Plot the time function prediction confidence intervals.\n'
                             '[!-- Preliminary feature alert! --!]\n'
                             '[!-- This feature is NOT thoroughly checked. '
                             'Read the code before use. Interpret at your own risk! --!]')

    parser = arg_utils.add_timefunc_argument(parser)

    # pixel of interest
    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                       help='initial pixel to plot in lat/lon coord')

    pixel.add_argument('--marker', type=str, default='o',
                       help='marker style (default: %(default)s).')
    pixel.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=6.0,
                       help='marker size (default: %(default)s).')
    pixel.add_argument('--lw', '--linewidth', dest='linewidth', type=float, default=0,
                       help='line width (default: %(default)s).')
    pixel.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0,
                       help='Edge width for the error bar (default: %(default)s)')

    # other groups
    parser = arg_utils.add_data_disp_argument(parser)
    parser = arg_utils.add_dem_argument(parser)
    parser = arg_utils.add_figure_argument(parser, figsize_img=True)
    parser = arg_utils.add_gps_argument(parser)
    parser = arg_utils.add_mask_argument(parser)
    parser = arg_utils.add_map_argument(parser)
    parser = arg_utils.add_memory_argument(parser)
    parser = arg_utils.add_reference_argument(parser)
    parser = arg_utils.add_save_argument(parser)
    parser = arg_utils.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check: --gps-comp option (not implemented for tsview yet)
    if inps.gps_component:
        msg = f'--gps-comp is not supported for {os.path.basename(__file__)}'
        raise NotImplementedError(msg)

    # check: --label option (same number as input files)
    if inps.file_label:
        if len(inps.file_label) != len(inps.file):
            raise Exception('input number of labels != number of files.')

    # check: coupled options
    if not inps.save_fig and (inps.outfile or not inps.disp_fig):
        inps.save_fig = True

    if inps.flip_lr or inps.flip_ud:
        inps.auto_flip = False

    if inps.ylim:
        inps.ylim = sorted(inps.ylim)

    if inps.zero_mask:
        inps.mask_file = 'no'

    if not inps.disp_fig_img:
        if not inps.yx and not inps.lalo:
            inps.disp_fig_img = True
            msg = 'WARNING: --yx/lalo is required for --no-show-img but NOT found! '
            msg += 'Ignore it and continue'
            print(msg)

    # default: -u / -c / --fig-size options
    inps.disp_unit = inps.disp_unit if inps.disp_unit else 'cm'
    inps.colormap = inps.colormap if inps.colormap else 'jet'
    inps.fig_size = inps.fig_size if inps.fig_size else [8.0, 4.5]

    return inps


###########################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.tsview import timeseriesViewer

    # run
    obj = timeseriesViewer(inps)
    obj.open()
    obj.plot()
    #obj.fig_img.canvas.mpl_disconnect(obj.cid)


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
