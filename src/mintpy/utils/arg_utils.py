"""Argument parsers."""
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Zhang Yunjun, Nov 2020                            #
#############################################################
# Recommend import:
#   from mintpy.utils import arg_utils
#   from mintpy.utils.arg_utils import create_argument_parser


import argparse
import math


##################################  generic parser  ####################################
def create_argument_parser(name=None, synopsis=None, description=None, epilog=None,
                           subparsers=None, formatter_class=argparse.RawTextHelpFormatter):
    """Create an argument parser.

    Parameters: name            - str, sub-command name, for sub-parser
                synopsis        - str, a brief summary of the script, for sub-parser
                description     - str, same as synopsis, plus optional note
                epilog          - str, reference, template options and example usage
                subparsers      - argparse._SubParsersAction
                                  https://docs.python.org/3/library/argparse.html#sub-commands
                formatter_class - argparse formatting class object
                                  https://docs.python.org/3/library/argparse.html#formatter-class
    Returns:    parser          - argparse.ArgumentParser object
    Examples:
        def create_parser(subparsers=None):
            synopsis = 'Resample radar-coded files into geo-coordinates or vice versa.'
            epilog = REFERENCE + '\n' + EXAMPLE
            name = __name__.split('.')[-1]
            parser = create_argument_parser(
                name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)
    """
    if subparsers:
        # for mintpy sub-command [used in linux with apt install]
        parser = subparsers.add_parser(
            name, description=description, formatter_class=formatter_class, epilog=epilog, help=synopsis)

    else:
        # for regular command usage
        parser = argparse.ArgumentParser(
            description=description, formatter_class=formatter_class, epilog=epilog)

    return parser



##################################  argument group  ####################################
def add_data_disp_argument(parser):
    """Argument group parser for data display options"""
    data = parser.add_argument_group('Data Display Options', 'Options to adjust the dataset display')
    data.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')
    data.add_argument('-u', '--unit', dest='disp_unit', metavar='UNIT',
                      help='unit for display.  Its priority > wrap')
    data.add_argument('--nd','--no-data-val','--no-data-value', dest='no_data_value', type=float,
                      help='Specify the no-data-value to be ignored and masked.')

    data.add_argument('--interp','--interpolation', dest='interpolation', default='nearest',
                      help='matplotlib interpolation method for imshow, e.g.:\n'
                           'none, antialiased, nearest, bilinear, bicubic, spline16, sinc, etc. Check more at:\n'
                           'https://matplotlib.org/stable/gallery/images_contours_and_fields/'
                           'interpolation_methods.html')
    data.add_argument('--wrap', action='store_true',
                      help='re-wrap data to display data in fringes.')
    data.add_argument('--wrap-range', dest='wrap_range', type=float, nargs=2,
                      default=[-1.*math.pi, math.pi], metavar=('MIN', 'MAX'),
                      help='range of one cycle after wrapping (default: %(default)s).')

    data.add_argument('--flip-lr', dest='flip_lr',
                      action='store_true', help='flip left-right')
    data.add_argument('--flip-ud', dest='flip_ud',
                      action='store_true', help='flip up-down')
    data.add_argument('--noflip', dest='auto_flip', action='store_false',
                      help='turn off auto flip for radar coordinate file')

    data.add_argument('--nmli','--num-multilook','--multilook-num', dest='multilook_num',
                      type=int, default=1, metavar='NUM',
                      help='multilook data in X and Y direction with a factor for display '
                           '(default: %(default)s).')
    data.add_argument('--nomultilook', '--no-multilook', dest='multilook', action='store_false',
                      help='do not multilook, for high quality display. \n'
                           'If multilook is True and multilook_num=1, '
                           'multilook_num will be estimated automatically.\n'
                           'Useful when displaying big datasets.')
    data.add_argument('--alpha', dest='transparency', type=float,
                      help='Data transparency. \n'
                           '0.0 - fully transparent, 1.0 - no transparency.')
    return parser


def add_dem_argument(parser):
    """Argument group parser for DEM display options"""
    dem = parser.add_argument_group('DEM', 'display topography in the background')
    dem.add_argument('-d', '--dem', dest='dem_file', metavar='DEM_FILE',
                     help='DEM file to show topography as background')
    dem.add_argument('--mask-dem', dest='mask_dem', action='store_true',
                     help='Mask out DEM pixels not coincident with valid data pixels')
    dem.add_argument('--dem-noshade', dest='disp_dem_shade', action='store_false',
                     help='do not show DEM shaded relief')
    dem.add_argument('--dem-nocontour', dest='disp_dem_contour', action='store_false',
                     help='do not show DEM contour lines')
    dem.add_argument('--dem-blend', dest='disp_dem_blend', action='store_true',
                     help='blend the DEM shade with input image to have a GMT-like impression.')

    # DEM contours
    dem.add_argument('--contour-smooth', dest='dem_contour_smooth', type=float, default=3.0, metavar='NUM',
                     help='[Contour] Topography contour smooth factor - sigma of Gaussian filter. \n'
                          'Set to 0.0 for no smoothing; (default: %(default)s).')
    dem.add_argument('--contour-step', dest='dem_contour_step', metavar='NUM', type=float, default=200.0,
                     help='[Contour] Topography contour step in meters (default: %(default)s).')
    dem.add_argument('--contour-lw','--contour-linewidth', dest='dem_contour_linewidth',
                     metavar='NUM', type=float, default=0.5,
                     help='[Contour] Topography contour linewidth (default: %(default)s).')

    # DEM shade
    dem.add_argument('--shade-az', dest='shade_azdeg', type=float, default=315., metavar='DEG',
                     help='[Shade] Azimuth angle (0-360, degrees clockwise from North) of the light source '
                          '(default: %(default)s).')
    dem.add_argument('--shade-alt', dest='shade_altdeg', type=float, default=45., metavar='DEG',
                     help='[Shade] Altitude (0-90, degrees up from horizontal) of the light source '
                          '(default: %(default)s).')
    dem.add_argument('--shade-min', dest='shade_min', type=float, default=-4000., metavar='MIN',
                     help='[Shade] Minimum height of shaded relief topography (default: %(default)s m).')
    dem.add_argument('--shade-max', dest='shade_max', type=float, default=999., metavar='MAX',
                     help='[Shade] Maximum height of shaded relief topography (default: max(DEM)+2000 m).')
    dem.add_argument('--shade-exag', dest='shade_exag', type=float, default=0.5,  metavar='NUM',
                     help='[Shade] Vertical exaggeration ratio (default: %(default)s).')

    # DEM-blended image
    dem.add_argument('--shade-frac', dest='shade_frac', type=float, default=0.5, metavar='NUM',
                     help='[Blend] Increases/decreases the contrast of the hillshade (default: %(default)s).')
    dem.add_argument('--base-color', dest='base_color', type=float, default=0.7, metavar='NUM',
                     help='[Blend] Topograhpy basemap greyish color ranges in [0,1] (default: %(default)s).')
    dem.add_argument('--blend-mode', dest='blend_mode', type=str, default='overlay',
                     choices={'hsv','overlay','soft'}, metavar='STR',
                     help='[Blend] Type of blending used to combine the colormapped data with illumated '
                          'topography.\n(choices: %(choices)s; default: %(default)s).\n'
                          'https://matplotlib.org/stable/gallery/specialty_plots/topographic_hillshading.html')
    return parser


def add_figure_argument(parser, figsize_img=False):
    """Argument group parser for figure options"""
    fig = parser.add_argument_group('Figure', 'Figure settings for display')
    fig.add_argument('--fontsize', dest='font_size',
                     type=int, help='font size')
    fig.add_argument('--fontcolor', dest='font_color',
                     default='k', help='font color (default: %(default)s).')

    # axis format
    fig.add_argument('--nowhitespace', dest='disp_whitespace',
                     action='store_false', help='do not display white space')
    fig.add_argument('--noaxis', dest='disp_axis',
                     action='store_false', help='do not display axis')
    fig.add_argument('--notick', dest='disp_tick',
                     action='store_false', help='do not display tick in x/y axis')
    fig.add_argument('--ylabel-rot', dest='ylabel_rot', type=float,
                     help='Y-axis tick label rotation in degree anti-clockwisely (default: %(default)s).\n'
                          'Set to 90 for a vertical y-axis tick labels')

    # colormap
    fig.add_argument('-c', '--colormap', dest='colormap',
                     help='colormap used for display, i.e. jet, cmy, RdBu, hsv, jet_r, viridis, etc.\n'
                          'More at https://mintpy.readthedocs.io/en/latest/api/colormaps/')
    fig.add_argument('--cm-lut','--cmap-lut', dest='cmap_lut', type=int, default=256, metavar='NUM',
                     help='number of increment of colormap lookup table (default: %(default)s).')
    fig.add_argument('--cm-vlist','--cmap-vlist', dest='cmap_vlist', type=float, nargs=3, default=[0.0, 0.7, 1.0],
                     help='list of 3 float numbers, for truncated colormap only (default: %(default)s).')

    # colorbar
    fig.add_argument('--nocbar', '--nocolorbar', dest='disp_cbar',
                     action='store_false', help='do not display colorbar')
    tik = fig.add_mutually_exclusive_group(required=False)
    tik.add_argument('--cbar-nbins', dest='cbar_nbins', metavar='NUM', type=int,
                     help='number of bins for colorbar.')
    tik.add_argument('--cbar-ticks', dest='cbar_ticks', nargs='+', metavar='NUM', type=float,
                     help='colorbar ticks.')
    fig.add_argument('--cbar-ext', dest='cbar_ext', default=None,
                     choices={'neither', 'min', 'max', 'both', None},
                     help='Extend setting of colorbar; based on data stat by default.')
    fig.add_argument('--cbar-label', dest='cbar_label', default=None, help='colorbar label')
    fig.add_argument('--cbar-loc', dest='cbar_loc', type=str, default='right',
                     help='colorbar location for single plot (default: %(default)s).')
    fig.add_argument('--cbar-size', dest='cbar_size', type=str, default="2%",
                     help='colorbar size and pad (default: %(default)s).')

    # title
    fig.add_argument('--notitle', dest='disp_title',
                     action='store_false', help='do not display title')
    fig.add_argument('--title-in', dest='fig_title_in',
                     action='store_true', help='draw title in/out of axes')
    fig.add_argument('--title','--fig-title','--figtitle', dest='fig_title',
                     help='Title shown in the figure.')
    fig.add_argument('--title4sen','--title4sentinel1', dest='disp_title4sentinel1', action='store_true',
                     help='display Sentinel-1 A/B and IPF info in title.')

    # size, subplots number and space
    if figsize_img:
        fig.add_argument('--figsize-img', dest='fig_size_img', metavar=('WID', 'LEN'), type=float, nargs=2,
                         help='figure size in inches for the image/map (for tsview.py) figure')
    fig.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='figure size in inches - width and length')
    fig.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=300,
                     help='DPI - dot per inch - for display/write (default: %(default)s).')
    fig.add_argument('--figext', dest='fig_ext', default='.png',
                     choices=['.emf', '.eps', '.pdf', '.png', '.ps', '.raw', '.rgba', '.svg', '.svgz'],
                     help='File extension for figure output file (default: %(default)s).')

    fig.add_argument('--fignum', dest='fig_num', type=int, metavar='NUM',
                     help='number of figure windows')
    fig.add_argument('--nrows', dest='fig_row_num', type=int, default=1, metavar='NUM',
                     help='subplot number in row')
    fig.add_argument('--ncols', dest='fig_col_num', type=int, default=1, metavar='NUM',
                     help='subplot number in column')

    fig.add_argument('--wspace', dest='fig_wid_space', type=float,
                     help='width space between subplots in inches')
    fig.add_argument('--hspace', dest='fig_hei_space', type=float,
                     help='height space between subplots in inches')
    fig.add_argument('--no-tight-layout', dest='fig_tight_layout', action='store_false',
                     help='disable automatic tight layout for multiple subplots')

    fig.add_argument('--coord', dest='fig_coord', choices=['geo','radar','yx'], default='geo',
                     help='Display axes in geo or yx coordinates '
                          '(for geocoded file only; default: %(default)s).')
    fig.add_argument('--animation', action='store_true',
                     help='enable animation mode')

    return parser


def add_gps_argument(parser):
    """Argument group parser for GPS options"""
    gps = parser.add_argument_group('GPS', 'GPS data to display')
    gps.add_argument('--show-gps', dest='disp_gps', action='store_true',
                     help='Show UNR GPS location within the coverage.')
    gps.add_argument('--mask-gps', dest='mask_gps', action='store_true',
                     help='Mask out GPS stations not coincident with valid data pixels')
    gps.add_argument('--gps-label', dest='disp_gps_label', action='store_true',
                     help='Show GPS site name')
    gps.add_argument('--gps-ms', dest='gps_marker_size', type=float, default=6,
                     help='Plot GPS value as scatter in size of ms**2 (default: %(default)s).')
    gps.add_argument('--gps-comp', dest='gps_component',
                     choices={'enu2los', 'hz2los', 'up2los', 'horz', 'vert'},
                     help='Plot GPS in color indicating deformation velocity direction')
    gps.add_argument('--gps-redo', dest='gps_redo', action='store_true',
                     help='Re-calculate GPS observations in LOS direction, '
                          'instead of read from existing CSV file.')
    gps.add_argument('--ref-gps', dest='ref_gps_site', type=str, help='Reference GPS site')
    gps.add_argument('--ex-gps', dest='ex_gps_sites', type=str, nargs='*',
                     help='Exclude GPS sites, require --gps-comp.')

    gps.add_argument('--gps-start-date', dest='gps_start_date', type=str, metavar='YYYYMMDD',
                     help='start date of GPS data, default is date of the 1st SAR acquisition')
    gps.add_argument('--gps-end-date', dest='gps_end_date', type=str, metavar='YYYYMMDD',
                     help='start date of GPS data, default is date of the last SAR acquisition')
    gps.add_argument('--horz-az','--hz-az', dest='horz_az_angle', type=float, default=-90.,
                     help='Azimuth angle (anti-clockwise from the north) of the horizontal movement in degrees\n'
                          'E.g.: -90. for east  direction [default]\n'
                          '       0.  for north direction\n'
                          'Set to the azimuth angle of the strike-slip fault to '
                          'show the fault-parallel displacement.')
    return parser


def add_mask_argument(parser):
    """Argument group parser for mask options"""
    mask = parser.add_argument_group('Mask', 'Mask file/options')
    mask.add_argument('-m','--mask', dest='mask_file', metavar='FILE',
                      help='mask file for display. "no" to turn OFF masking.')
    mask.add_argument('--mask-vmin', dest='mask_vmin', type=float,
                      help='hide pixels with mask value < vmin (default: %(default)s).')
    mask.add_argument('--mask-vmax', dest='mask_vmax', type=float,
                      help='hide pixels with mask value > vmax (default: %(default)s).')

    mask.add_argument('--zm','--zero-mask', dest='zero_mask', action='store_true',
                      help='mask pixels with zero value.')
    return parser


def add_map_argument(parser):
    """Argument group parser for map options"""
    mapg = parser.add_argument_group('Map', 'for one subplot in geo-coordinates only')

    # coastline
    mapg.add_argument('--coastline', dest='coastline', type=str, choices={'10m', '50m', '110m'},
                      help="Draw coastline with specified resolution (default: %(default)s).\n"
                           "This will enable --lalo-label option.\n"
                           "https://scitools.org.uk/cartopy/docs/latest/reference/generated/"
                           "cartopy.mpl.geoaxes.GeoAxes.html")
    mapg.add_argument('--coastline-lw', '--coastline-linewidth', dest='coastline_linewidth',
                      metavar='NUM', type=float, default=1,
                      help='Coastline linewidth (default: %(default)s).')

    # faultline
    mapg.add_argument('--faultline', dest='faultline_file', type=str,
                      help='Draw fault line using specified GMT lonlat file.')
    mapg.add_argument('--faultline-lw', '--faultline-linewidth', dest='faultline_linewidth',
                      metavar='NUM', type=float, default=0.5,
                      help='Faultline linewidth (default: %(default)s).')
    mapg.add_argument('--faultline-min-dist','--faultline-min-len', dest='faultline_min_dist',
                      metavar='NUM', type=float, default=0.1,
                      help='Show fault segments with length >= X km (default: %(default)s).')

    # lalo label
    mapg.add_argument('--lalo-label', dest='lalo_label', action='store_true',
                      help='Show N, S, E, W tick label for plot in geo-coordinate.\n'
                           'Useful for final figure output.')
    mapg.add_argument('--lalo-step', dest='lalo_step', metavar='DEG',
                      type=float, help='Lat/lon step for lalo-label option.')
    mapg.add_argument('--lalo-max-num', dest='lalo_max_num', type=int, default=3, metavar='NUM',
                      help='Maximum number of lalo tick label (default: %(default)s).')
    mapg.add_argument('--lalo-loc', dest='lalo_loc', type=int, nargs=4, default=[1, 0, 0, 1],
                      metavar=('left', 'right', 'top', 'bottom'),
                      help='Draw lalo label in [left, right, top, bottom] (default: %(default)s).')
    mapg.add_argument('--lalo-off','--lalo-offset', dest='lalo_offset', type=float, nargs=2,
                      help='Distance between tick and label in points (default: %(default)s).\n'
                           'Set to negative value, e.g. -36 -18, to move the ticklabel inside the plot.')
    mapg.add_argument('--lalo-fs','--lalo-fontsize', dest='lalo_font_size', type=float,
                      help='Lalo label font size in points (default: %(default)s).')

    #mapg.add_argument('--proj', '--projection', '--map-proj', dest='map_projection', metavar='NAME',
    #                  help='map projection when plotting in geo-coordinate.\n'
    #                       'Default: PlateCarree / UTM for units in degrees / meters.\n'
    #                       'Check the link below for the full list of supported projections:\n'
    #                       'https://scitools.org.uk/cartopy/docs/latest/crs/projections.html\n\n')

    # scale bar
    mapg.add_argument('--scalebar', nargs=3, metavar=('LEN', 'X', 'Y'), type=float,
                      default=[0.2, 0.2, 0.1],
                      help='scale bar distance and location in ratio (default: %(default)s).\n' +
                           '\tdistance in ratio of total width\n' +
                           '\tlocation in X/Y in ratio with respect to the lower left corner\n' +
                           '--scalebar 0.2 0.2 0.1  #for lower left  corner\n' +
                           '--scalebar 0.2 0.2 0.8  #for upper left  corner\n' +
                           '--scalebar 0.2 0.8 0.1  #for lower right corner\n' +
                           '--scalebar 0.2 0.8 0.8  #for upper right corner\n')
    mapg.add_argument('--noscalebar', '--nosbar', dest='disp_scalebar',
                      action='store_false', help='do not display scale bar.')
    mapg.add_argument('--scalebar-pad','--sbar-pad', dest='scalebar_pad', type=float, default=0.05,
                      help='scale bar label pad in ratio of scalebar width (default: %(default)s).')
    mapg.add_argument('--scalebar-lw','--scalebar-linewidth', dest='scalebar_linewidth', type=float,
                      default=2.0, help='scale bar symbol line width (default: %(default)s).')
    return parser


def add_memory_argument(parser):
    """Argument parser for memory usage options"""
    parser.add_argument('--ram', '--memory', dest='maxMemory', type=float, default=4.0,
                        help='Max amount of memory in GB to use (default: %(default)s).\n' +
                             'Adjust according to your computer memory.')
    return parser


def add_parallel_argument(parser):
    """Argument group parser for parallel computing options"""
    # from mintpy.objects.cluster import CLUSTER_LIST
    CLUSTER_LIST = ['lsf', 'pbs', 'slurm', 'local']

    par = parser.add_argument_group('parallel', 'parallel processing using dask')
    par.add_argument('-c', '--cluster', '--cluster-type', dest='cluster', type=str,
                     choices=CLUSTER_LIST,
                     help='Cluster to use for parallel computing (default: %(default)s to turn OFF).')
    par.add_argument('--num-worker', dest='numWorker', type=str, default='4',
                     help='Number of workers to use (default: %(default)s).')
    par.add_argument('--config', '--config-name', dest='config', type=str, default=None,
                     help='Configuration name to use in dask.yaml (default: %(default)s).')
    return parser


def add_point_argument(parser):
    """Argument group parser for point display options"""
    ppt = parser.add_argument_group('Point', 'Plot points defined by y/x or lat/lon')
    ppt.add_argument('--pts-marker', dest='pts_marker', type=str, default='k^',
                     help='Marker of points of interest (default: %(default)s).')
    ppt.add_argument('--pts-ms', dest='pts_marker_size', type=float, default=6.,
                     help='Marker size for points of interest (default: %(default)s).')

    pts = ppt.add_mutually_exclusive_group(required=False)
    pts.add_argument('--pts-yx', dest='pts_yx', type=int, nargs=2, metavar=('Y', 'X'),
                     help='Point in Y/X')
    pts.add_argument('--pts-lalo', dest='pts_lalo', type=float, nargs=2, metavar=('LAT', 'LON'),
                     help='Point in Lat/Lon')
    pts.add_argument('--pts-file', dest='pts_file', type=str,
                     help='Text file for point(s) in lat/lon column')
    return parser


def add_reference_argument(parser, plot=True):
    """Argument group parser for (spatial / temporal) referencing options"""

    goal = 'display' if plot else 'estimation'
    ref = parser.add_argument_group('Reference date / point',
                                    f'Modify reference in time / space for {goal}')

    # reference date / pixel
    ref.add_argument('--ref-date', dest='ref_date', metavar='DATE',
                     help=f'Change reference date for {goal}')
    ref.add_argument('--ref-lalo', dest='ref_lalo', metavar=('LAT', 'LON'), type=float, nargs=2,
                     help=f'Change reference point in LAT/LON for {goal}')
    ref.add_argument('--ref-yx', dest='ref_yx', metavar=('Y', 'X'), type=int, nargs=2,
                     help=f'Change reference point in Y/X for {goal}')

    # reference pixel - plotting style
    if plot:
        ref.add_argument('--noreference', dest='disp_ref_pixel',
                         action='store_false', help='do not show reference point')
        ref.add_argument('--ref-marker', dest='ref_marker', default='ks',
                         help='marker of reference pixel (default: %(default)s).')
        ref.add_argument('--ref-size', dest='ref_marker_size', metavar='NUM', type=int, default=6,
                         help='marker size of reference point (default: %(default)s).')
    return parser


def add_save_argument(parser):
    """Argument group parser for figure save options"""
    save = parser.add_argument_group('Save/Output', 'Save figure and write to file(s)')
    save.add_argument('-o', '--outfile', type=str, nargs='*',
                      help="save the figure with assigned filename.\n"
                           "By default, it's calculated based on the input file name.")
    save.add_argument('--save', dest='save_fig', action='store_true',
                      help='save the figure')
    save.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                      help='save and do not display the figure')
    save.add_argument('--update', dest='update_mode', action='store_true',
                      help='enable update mode for save figure: skip running if\n'+
                           '\t1) output file already exists AND\n'+
                           '\t2) output file is newer than input file.')
    return parser


def add_subset_argument(parser, geo=True):
    """Argument group parser for subset options"""
    sub = parser.add_argument_group('Subset', 'Display dataset in subset range')
    sub.add_argument('--sub-x','--subx','--subset-x', dest='subset_x', type=int, nargs=2,
                     metavar=('XMIN', 'XMAX'), help='subset display in x/cross-track/range direction')
    sub.add_argument('--sub-y','--suby','--subset-y', dest='subset_y', type=int, nargs=2,
                     metavar=('YMIN', 'YMAX'), help='subset display in y/along-track/azimuth direction')
    if geo:
        sub.add_argument('--sub-lat','--sublat','--subset-lat', dest='subset_lat', type=float, nargs=2,
                         metavar=('LATMIN', 'LATMAX'), help='subset display in latitude')
        sub.add_argument('--sub-lon','--sublon','--subset-lon', dest='subset_lon', type=float, nargs=2,
                         metavar=('LONMIN', 'LONMAX'), help='subset display in longitude')
    return parser


def add_timefunc_argument(parser):
    """Argument group parser for time functions"""
    model = parser.add_argument_group('Deformation Model', 'A suite of time functions')

    model.add_argument('--poly', '--polynomial', '--poly-order', dest='polynomial', type=int, default=1,
                       help='a polynomial function with the input degree (default: %(default)s). E.g.:\n'
                            '--poly 1                               # linear\n'
                            '--poly 2                               # quadratic\n'
                            '--poly 3                               # cubic\n')

    model.add_argument('--periodic', '--period', '--peri', dest='periodic', type=float, nargs='+', default=[],
                       help='periodic function(s) with period in decimal years (default: %(default)s). E.g.:\n'
                            '--periodic 1.0                         # an annual cycle\n'
                            '--periodic 1.0 0.5                     # an annual cycle plus a semi-annual cycle\n')

    model.add_argument('--step','--step-date', dest='stepDate', type=str, nargs='+', default=[],
                       help='step function(s) at YYYYMMDD (default: %(default)s). E.g.:\n'
                            '--step 20061014                        # coseismic step  at 2006-10-14T00:00\n'
                            '--step 20110311 20120928T1733          # coseismic steps at 2011-03-11T00:00 and 2012-09-28T17:33\n')

    model.add_argument('--polyline', dest='polyline', type=str, nargs='+', default=[],
                       help='polyline segment(s) starting at YYYYMMDD (default: %(default)s). E.g.:\n'
                            '--polyline 20190101                    # extra velocity   since 2019-01-01T00:00\n'
                            '--polyline 20190101 20200501T1725      # extra velocities since 2019-01-01T00:00 and 2020-05-01T17:25\n')

    model.add_argument('--exp', '--exponential', dest='exp', type=str, nargs='+', action='append', default=[],
                       help='exponential function(s) defined by onset time(s) and characteristic time(s) tau in days (default: %(default)s). E.g.:\n'
                            '--exp  20181026 60                     # one exp w/ onset at 2018-10-26       w/ tau=60  days\n'
                            '--exp  20181026T1355 60 120            # 1st exp w/ onset at 2018-10-26T13:55 w/ tau=60  days\n'
                            '                                       # 2nd exp w/ onset at 2018-10-26T13:55 w/ tau=120 days\n'
                            '--exp  20161231 80 --exp 20190125 100  # 1st exp w/ onset at 2016-12-31       w/ tau=80  days\n'
                            '                                       # 2nd exp w/ onset at 2019-01-25       w/ tau=100 days')

    model.add_argument('--log', '--logarithmic', dest='log', type=str, nargs='+', action='append', default=[],
                       help='logarithmic function(s) defined by onset time(s) and characteristic time(s) tau in days (default: %(default)s). E.g.:\n'
                            '--log  20181016 90                     # one log w/ onset at 2018-10-16       w/ tau=90  days\n'
                            '--log  20181016T1733 90 240            # 1st log w/ onset at 2018-10-16T17:33 w/ tau=90  days\n'
                            '                                       # 2nd log w/ onset at 2018-10-16T17:33 w/ tau=240 days\n'
                            '--log  20161231 60 --log 20190125 180  # 1st log w/ onset at 2016-12-31       w/ tau=60  days\n'
                            '                                       # 2nd log w/ onset at 2019-01-25       w/ tau=180 days\n')

    return parser
