#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################
# Recommend import:
#   from mintpy.cli import view


import os
import sys

from mintpy.utils import arg_utils

##################################################################################################
EXAMPLE = """example:
  view.py velocity.h5
  view.py velocity.h5 velocity --wrap --wrap-range -2 2 -c cmy --lalo-label
  view.py velocity.h5 --ref-yx  210 566                              #change reference pixel for display
  view.py velocity.h5 --sub-lat 31.05 31.10 --sub-lon 130.05 130.10  #subset in lalo / yx
  view.py velocity.h5 velocity --mask waterBody.h5 --mask-vmax 1

  view.py timeseries.h5
  view.py timeseries.h5 --ref-date 20101120     #change reference date
  view.py timeseries.h5 --ex drop_date.txt      #exclude dates to plot
  view.py timeseries.h5 '*2017*' '*2018*'       #all acquisitions in 2017 and 2018
  view.py timeseries.h5 20200616_20200908       #reconstruct interferogram on the fly

  view.py ifgramStack.h5 coherence
  view.py ifgramStack.h5 unwrapPhase-           #unwrapPhase only in the presence of unwrapPhase_bridging
  view.py ifgramStack.h5 -n 6                   #the 6th slice
  view.py ifgramStack.h5 20171010_20171115      #all data      related with 20171010_20171115
  view.py ifgramStack.h5 'coherence*20171010*'  #all coherence related with 20171010

  # GPS (for one subplot in geo-coordinates only)
  view.py geo_velocity_msk.h5 velocity --show-gps --gps-label   #show locations of available GPS
  view.py geo_velocity_msk.h5 velocity --show-gps --gps-comp enu2los --ref-gps GV01
  view.py geo_timeseries_ERA5_ramp_demErr.h5 20180619 --ref-date 20141213 --show-gps --gps-comp enu2los --ref-gps GV01

  # Faults
  view.py filt_dense_offsets.bil range --faultline simple_fault_confident.lonlat

  # Save and Output
  view.py velocity.h5 --save
  view.py velocity.h5 --nodisplay
  view.py geo_velocity.h5 velocity --nowhitespace
"""


def create_parser(subparsers=None):
    synopsis = 'Plot InSAR Product in 2D'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    infile = parser.add_argument_group('Input File', 'File/Dataset to display')
    infile.add_argument('file', type=str, help='file for display')
    infile.add_argument('dset', type=str, nargs='*', default=[],
                        help='optional - dataset(s) to display (default: %(default)s).')
    infile.add_argument('-n', '--dset-num', dest='dsetNumList', metavar='NUM', type=int, nargs='*', default=[],
                        help='optional - order number of date/dataset(s) to display (default: %(default)s).')
    infile.add_argument('--nosearch', dest='search_dset', action='store_false',
                        help='Disable glob search for input dset.')
    infile.add_argument('--ex', '--exclude', dest='exDsetList', metavar='Dset', nargs='*', default=[],
                        help='dates will not be displayed (default: %(default)s).')
    parser.add_argument('--show-kept','--show-kept-ifgram', dest='plot_drop_ifgram', action='store_false',
                        help='display kept interferograms only, without dropped interferograms')

    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing (default: %(default)s).')

    parser.add_argument('--math', dest='math_operation', choices={'square','sqrt','reverse','inverse','rad2deg','deg2rad'},
                        help='Apply the math operation before displaying [for single subplot ONLY].\n'
                             'E.g. plot the std. dev. of the variance file.\n'
                             '  square  = x^2\n'
                             '  sqrt    = x^1/2\n'
                             '  reverse = x * -1\n'
                             '  inverse = 1 / x')

    parser = arg_utils.add_data_disp_argument(parser)
    parser = arg_utils.add_dem_argument(parser)
    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_gps_argument(parser)
    parser = arg_utils.add_mask_argument(parser)
    parser = arg_utils.add_map_argument(parser)
    parser = arg_utils.add_memory_argument(parser)
    parser = arg_utils.add_point_argument(parser)
    parser = arg_utils.add_reference_argument(parser)
    parser = arg_utils.add_save_argument(parser)
    parser = arg_utils.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import ptime, readfile

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check: invalid file inputs
    for key in ['file','dem_file','mask_file','pts_file']:
        fname = vars(inps)[key]
        if fname not in [None, 'no'] and not os.path.isfile(fname):
            raise FileNotFoundError(f'input {key} file {fname} NOT exist!')

    # check: --exclude option
    if inps.exDsetList:
        inps.exDsetList = ptime.read_date_list(inps.exDsetList)

    # check: coupled option behaviors
    if not inps.save_fig and (inps.outfile or not inps.disp_fig):
        # turn ON save_fig if a) --output OR b) figure not shown
        inps.save_fig = True

    if inps.lalo_step:
        inps.lalo_label = True

    if inps.zero_mask:
        # turn OFF default mask file detection for --zero-mask
        # extra manual mask file is still supported
        if not inps.mask_file:
            inps.mask_file = 'no'

    if not inps.disp_whitespace:
        inps.disp_axis = False
        inps.disp_title = False
        inps.disp_cbar = False

    if not inps.disp_axis:
        inps.disp_tick = False

    if inps.flip_lr or inps.flip_ud:
        inps.auto_flip = False

    if inps.disp_dem_blend:
        inps.disp_dem_shade = False
        # --dem-blend option requires --dem option
        if inps.dem_file is None:
            parser.error("--dem-blend requires -d/-dem.")
        # --cbar-ext option is ignored
        if '--cbar-ext' in inps.argv:
            print('WARNING: --cbar-ext is NOT compatible with --dem-blend, ignore --cbar-ext and continue.')

    # check: conflicted options (geo-only options if inpput file is in radar-coordinates)
    geo_opt_names = ['--coord', '--show-gps', '--coastline', '--lalo-label', '--lalo-step', '--scalebar', '--faultline']
    geo_opt_names = list(set(geo_opt_names) & set(inps.argv))
    if geo_opt_names and 'Y_FIRST' not in readfile.read_attribute(inps.file).keys():
        for opt_name in geo_opt_names:
            print(f'WARNING: {opt_name} is NOT supported for files in radar-coordinate, ignore it and continue.')

    # check: --noverbose option
    # print view.py command line if --noverbose (used in smallbaselineApp.py)
    if not inps.print_msg:
        print('view.py', ' '.join(inps.argv))

    return inps


#########################################  Main Function  ########################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.view import viewer

    # run
    obj = viewer(iargs=iargs)
    obj.configure(inps)
    if obj.flag == 'run':
        obj.plot()


##################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
