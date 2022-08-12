############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
from mintpy.utils import arg_utils


##################################################################################################
EXAMPLE = """example:
  view.py velocity.h5
  view.py velocity.h5 velocity --wrap --wrap-range -2 2 -c cmy --lalo-label
  view.py velocity.h5 --ref-yx  210 566                              #change reference pixel for display
  view.py velocity.h5 --sub-lat 31.05 31.10 --sub-lon 130.05 130.10  #subset in lalo / yx

  view.py timeseries.h5
  view.py timeseries.h5 -m no                   #do not use auto mask
  view.py timeseries.h5 --ref-date 20101120     #change reference date
  view.py timeseries.h5 --ex drop_date.txt      #exclude dates to plot
  view.py timeseries.h5 '*2017*' '*2018*'       #all acquisitions in 2017 and 2018
  view.py timeseries.h5 20200616_20200908       #reconstruct interferogram on the fly

  view.py ifgramStack.h5 coherence
  view.py ifgramStack.h5 unwrapPhase-           #unwrapPhase only in the presence of unwrapPhase_bridging
  view.py ifgramStack.h5 -n 6                   #the 6th slice
  view.py ifgramStack.h5 20171010_20171115      #all data      related with 20171010_20171115
  view.py ifgramStack.h5 'coherence*20171010*'  #all coherence related with 20171010
  view.py ifgramStack.h5 unwrapPhase-20070927_20100217 --zero-mask --wrap     #wrapped phase
  view.py ifgramStack.h5 unwrapPhase-20070927_20100217 --mask ifgramStack.h5  #mask using connected components

  # GPS (for one subplot in geo-coordinates only)
  view.py geo_velocity_msk.h5 velocity --show-gps --gps-label   #show locations of available GPS
  view.py geo_velocity_msk.h5 velocity --show-gps --gps-comp enu2los --ref-gps GV01
  view.py geo_timeseries_ERA5_ramp_demErr.h5 20180619 --ref-date 20141213 --show-gps --gps-comp enu2los --ref-gps GV01

  # Save and Output
  view.py velocity.h5 --save
  view.py velocity.h5 --nodisplay
  view.py geo_velocity.h5 velocity --nowhitespace
"""

PLOT_TEMPLATE = """Plot Setting:
  plot.name          = 'Yunjun et al., 2016, AGU, Fig 4f'
  plot.type          = LOS_VELOCITY
  plot.startDate     =
  plot.endDate       =
  plot.displayUnit   = cm/yr
  plot.displayMin    = -2
  plot.displayMax    = 2
  plot.colormap      = jet
  plot.subset.lalo   = 33.05:33.15, 131.15:131.27
  plot.seed.lalo = 33.0651, 131.2076
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

    parser.add_argument('--plot-setting', dest='disp_setting_file',
                        help='Template file with plot setting.\n'+PLOT_TEMPLATE)

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
    from ..utils import ptime, readfile
    from ..view import update_inps_with_display_setting_file

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check invalid file inputs
    for key in ['file','dem_file','mask_file','pts_file']:
        fname = vars(inps)[key]
        if fname not in [None, 'no'] and not os.path.isfile(fname):
            raise FileNotFoundError('input {} file {} NOT exist!'.format(key, fname))

    # --exclude
    if inps.exDsetList:
        inps.exDsetList = ptime.read_date_list(inps.exDsetList)

    # If output flie name assigned or figure shown is turned off, turn on the figure save
    if inps.outfile or not inps.disp_fig:
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

    # check geo-only options for files in radar-coordinates
    geo_opt_names = ['--coord', '--show-gps', '--coastline', '--lalo-label', '--lalo-step', '--scalebar']
    geo_opt_names = list(set(geo_opt_names) & set(inps.argv))
    if geo_opt_names and 'Y_FIRST' not in readfile.read_attribute(inps.file).keys():
        for opt_name in geo_opt_names:
            print(f'WARNING: {opt_name} is NOT supported for files in radar-coordinate, ignore it and continue.')

    # verbose print using --noverbose option
    from ..import view
    view.vprint = print if inps.print_msg else lambda *args, **kwargs: None
    # print view.py command line if --noverbose (used in smallbaselineApp.py)
    if not inps.print_msg:
        print('view.py', ' '.join(inps.argv))

    if inps.disp_setting_file:
        inps = update_inps_with_display_setting_file(inps, inps.disp_setting_file)

    # Backend setting
    if not inps.disp_fig:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')

    return inps


def prep_slice(cmd, auto_fig=False):
    """Prepare data from command line as input, for easy call plot_slice() externally
    Parameters: cmd  - string, command to be run in terminal
    Returns:    data - 2D np.ndarray, data to be plotted
                atr  - dict, metadata
                inps - namespace, input argument for plot setup
    Example:
        subplot_kw = dict(projection=ccrs.PlateCarree())
        fig, ax = plt.subplots(figsize=[4, 3], subplot_kw=subplot_kw)
        W, N, E, S = (-91.670, -0.255, -91.370, -0.515)    # geo_box
        cmd = 'view.py geo_velocity.h5 velocity --mask geo_maskTempCoh.h5 --dem srtm1.dem --dem-nocontour '
        cmd += f'--sub-lon {W} {E} --sub-lat {S} {N} -c jet -v -3 10 '
        cmd += '--cbar-loc bottom --cbar-nbins 3 --cbar-ext both --cbar-size 5% '
        cmd += '--lalo-step 0.2 --lalo-loc 1 0 1 0 --scalebar 0.3 0.80 0.05 --notitle'
        data, atr, inps = view.prep_slice(cmd)
        ax, inps, im, cbar = view.plot_slice(ax, data, atr, inps)
        plt.show()
    """
    import numpy as np
    from .. import view
    from mintpy.utils import readfile, plot as pp
    from ..view import read_input_file_info, update_inps_with_file_metadata, update_data_with_plot_inps

    inps = cmd_line_parse(cmd.split()[1:])
    view.vprint(cmd)
    inps, atr = read_input_file_info(inps)
    inps = update_inps_with_file_metadata(inps, atr)

    inps.msk, inps.mask_file = pp.read_mask(inps.file,
                                            mask_file=inps.mask_file,
                                            datasetName=inps.dset[0],
                                            box=inps.pix_box,
                                            vmin=inps.mask_vmin,
                                            vmax=inps.mask_vmax,
                                            print_msg=inps.print_msg)

    # read data
    data, atr = readfile.read(inps.file,
                              datasetName=inps.dset[0],
                              box=inps.pix_box,
                              print_msg=inps.print_msg)
    # reference in time
    if inps.ref_date:
        data -= readfile.read(inps.file,
                              datasetName=inps.ref_date,
                              box=inps.pix_box,
                              print_msg=False)[0]

    # reference in space for unwrapPhase
    if (inps.key in ['ifgramStack']
            and inps.dset[0].split('-')[0].startswith('unwrapPhase')
            and 'REF_Y' in atr.keys()):
        ref_y, ref_x = int(atr['REF_Y']), int(atr['REF_X'])
        ref_data = readfile.read(inps.file,
                                 datasetName=inps.dset[0],
                                 box=(ref_x, ref_y, ref_x+1, ref_y+1),
                                 print_msg=False)[0]
        data[data != 0.] -= ref_data

    # masking
    if inps.zero_mask:
        data = np.ma.masked_where(data == 0., data)
    if inps.msk is not None:
        data = np.ma.masked_where(inps.msk == 0., data)
    else:
        inps.msk = np.ones(data.shape, dtype=np.bool_)
    # update/save mask info
    if np.ma.is_masked(data):
        inps.msk *= ~data.mask
        inps.msk *= ~np.isnan(data.data)
    else:
        inps.msk *= ~np.isnan(data)

    data, inps = update_data_with_plot_inps(data, atr, inps)

    # matplotlib.Axes
    if auto_fig == True:
        import matplotlib.pyplot as plt
        figsize = [i/2.0 for i in inps.fig_size]
        subplot_kw = dict(projection=inps.map_proj_obj) if inps.map_proj_obj is not None else {}
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=subplot_kw)
        return data, atr, inps, ax
    else:
        return data, atr, inps


#########################################  Main Function  ########################################
def main(iargs=None):
    from ..view import viewer
    inps = cmd_line_parse(iargs)
    obj = viewer(iargs=iargs)
    obj.configure(inps)
    if obj.flag == 'run':
        obj.plot()


##################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
