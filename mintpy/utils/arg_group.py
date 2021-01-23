#!/usr/bin/env python3
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Zhang Yunjun, Nov 2020                            #
#############################################################
# Recommend import:
#     from mintpy.utils import arg_group


import argparse
import numpy as np



def add_data_disp_argument(parser):
    """Argument group parser for data display options"""
    data = parser.add_argument_group('Data Display Options', 'Options to adjust the dataset display')
    data.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')
    data.add_argument('-u', '--unit', dest='disp_unit', metavar='UNIT',
                      help='unit for display.  Its priority > wrap')

    data.add_argument('--wrap', action='store_true',
                      help='re-wrap data to display data in fringes.')
    data.add_argument('--wrap-range', dest='wrap_range', type=float, nargs=2,
                      default=[-1.*np.pi, np.pi], metavar=('MIN', 'MAX'),
                      help='range of one cycle after wrapping (default: %(default)s).')

    data.add_argument('--flip-lr', dest='flip_lr',
                      action='store_true', help='flip left-right')
    data.add_argument('--flip-ud', dest='flip_ud',
                      action='store_true', help='flip up-down')
    data.add_argument('--noflip', dest='auto_flip', action='store_false',
                      help='turn off auto flip for radar coordinate file')

    data.add_argument('--multilook-num', dest='multilook_num', type=int, default=1, metavar='NUM',
                      help='multilook data in X and Y direction with a factor for display (default: %(default)s).')
    data.add_argument('--nomultilook', '--no-multilook', dest='multilook', action='store_false',
                      help='do not multilook, for high quality display. \n'
                           'If multilook and multilook_num=1, multilook_num will be estimated automatically.\n'
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
    dem.add_argument('--dem-noshade', dest='disp_dem_shade', action='store_false',
                     help='do not show DEM shaded relief')
    dem.add_argument('--dem-nocontour', dest='disp_dem_contour', action='store_false',
                     help='do not show DEM contour lines')

    dem.add_argument('--contour-smooth', dest='dem_contour_smooth', type=float, default=3.0,
                     help='Background topography contour smooth factor - sigma of Gaussian filter. \n'
                          'Set to 0.0 for no smoothing; (default: %(default)s).')
    dem.add_argument('--contour-step', dest='dem_contour_step', metavar='NUM', type=float, default=200.0,
                     help='Background topography contour step in meters (default: %(default)s).')
    dem.add_argument('--contour-linewidth', dest='dem_contour_linewidth', metavar='NUM', type=float, default=0.5,
                     help='Background topography contour linewidth (default: %(default)s).')

    dem.add_argument('--shade-az', dest='shade_azdeg', type=float, default=315., metavar='DEG',
                     help='The azimuth (0-360, degrees clockwise from North) of the light source (default: %(default)s).')
    dem.add_argument('--shade-alt', dest='shade_altdeg', type=float, default=45., metavar='DEG',
                     help='The altitude (0-90, degrees up from horizontal) of the light source (default: %(default)s).')

    dem.add_argument('--shade-min', dest='shade_min', type=float, default=-4000., metavar='MIN',
                     help='The min height in m of colormap of shaded relief topography (default: %(default)s).')
    dem.add_argument('--shade-max', dest='shade_max', type=float, default=999., metavar='MAX',
                     help='The max height of colormap of shaded relief topography (default: max(DEM)+2000).')
    dem.add_argument('--shade-exag', dest='shade_exag', type=float, default=0.5,
                     help='Vertical exaggeration ratio (default: %(default)s).')
    return parser


def add_figure_argument(parser):
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

    # colormap
    fig.add_argument('-c', '--colormap', dest='colormap',
                     help='colormap used for display, i.e. jet, cmy, RdBu, hsv, jet_r, temperature, viridis, etc.\n'
                          'More at https://mintpy.readthedocs.io/en/latest/api/colormaps/')
    fig.add_argument('--cm-lut','--cmap-lut', dest='cmap_lut', type=int, default=256, metavar='NUM',
                     help='number of increment of colormap lookup table (default: %(default)s).')
    fig.add_argument('--cm-vlist','--cmap-vlist', dest='cmap_vlist', type=float, nargs=3, default=[0.0, 0.7, 1.0],
                     help='list of 3 float numbers, for truncated colormap only (default: %(default)s).')

    # colorbar
    fig.add_argument('--nocbar', '--nocolorbar', dest='disp_cbar',
                     action='store_false', help='do not display colorbar')
    fig.add_argument('--cbar-nbins', dest='cbar_nbins', metavar='NUM',
                     type=int, help='number of bins for colorbar.')
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
    fig.add_argument('--figtitle', dest='fig_title',
                     help='Title shown in the figure.')
    fig.add_argument('--title4sen','--title4sentinel1', dest='disp_title4sentinel1', action='store_true',
                     help='display Sentinel-1 A/B and IPF info in title.')

    # size, subplots number and space
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

    fig.add_argument('--coord', dest='fig_coord', choices=['radar', 'geo'], default='geo',
                     help='Display in radar/geo coordination system, for geocoded file only (default: %(default)s).')
    fig.add_argument('--animation', action='store_true',
                     help='enable animation mode')

    return parser


def add_gps_argument(parser):
    """Argument group parser for GPS options"""
    gps = parser.add_argument_group('GPS', 'GPS data to display')
    gps.add_argument('--show-gps', dest='disp_gps', action='store_true',
                     help='Show UNR GPS location within the coverage.')
    gps.add_argument('--gps-label', dest='disp_gps_label', action='store_true',
                     help='Show GPS site name')
    gps.add_argument('--gps-comp', dest='gps_component', choices={'enu2los', 'hz2los', 'up2los'},
                     help='Plot GPS in color indicating deformation velocity direction')
    gps.add_argument('--ref-gps', dest='ref_gps_site', type=str, help='Reference GPS site')

    gps.add_argument('--gps-start-date', dest='gps_start_date', type=str, metavar='YYYYMMDD',
                     help='start date of GPS data, default is date of the 1st SAR acquisiton')
    gps.add_argument('--gps-end-date', dest='gps_end_date', type=str, metavar='YYYYMMDD',
                     help='start date of GPS data, default is date of the last SAR acquisiton')
    return parser


def add_mask_argument(parser):
    """Argument group parser for mask options"""
    mask = parser.add_argument_group('Mask', 'Mask file/options')
    mask.add_argument('-m','--mask', dest='mask_file', metavar='FILE',
                      help='mask file for display. "no" to turn OFF masking.')
    mask.add_argument('--zm','--zero-mask', dest='zero_mask', action='store_true',
                      help='mask pixels with zero value.')
    return parser


def add_map_argument(parser):
    """Argument group parser for map options"""
    mapg = parser.add_argument_group('Map', 'for one subplot in geo-coordinates only')
    mapg.add_argument('--coastline', dest='coastline', type=str, choices={'10m', '50m', '110m'},
                      help="Draw coastline with specified resolution (default: %(default)s).\n"
                           "This will enable --lalo-label option.\n"
                           "Link: https://scitools.org.uk/cartopy/docs/latest/matplotlib/geoaxes.html"
                           "#cartopy.mpl.geoaxes.GeoAxes.coastlines")

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

    mapg.add_argument('--projection', dest='map_projection', metavar='NAME', default='PlateCarree',
                      choices={'PlateCarree', 'LambertConformal'},
                      help='map projection when plotting in geo-coordinate (default: %(default)s).\n'
                           'https://scitools.org.uk/cartopy/docs/latest/crs/projections.html\n\n')

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
    mapg.add_argument('--scalebar-pad','--sbar-pad', dest='scalebar_pad', type=float,
                      default=0.05, help='scale bar label pad in ratio of scalebar width (default: %(default)s).')
    return parser


def add_memory_argument(parser):
    """Argument parser for memory usage options"""
    parser.add_argument('--ram', '--memory', dest='maxMemory', type=float, default=4.0,
                        help='Max amount of memory in GB to use (default: %(default)s).\n' +
                             'Adjust according to your computer memory.')
    return parser


def add_parallel_argument(parser):
    """Argument group parser for parallel computing options"""
    from mintpy.objects.cluster import CLUSTER_LIST

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


def add_reference_argument(parser):
    """Argument group parser for (spatial / temporal) referencing options"""
    ref = parser.add_argument_group('Reference', 'Show / Modify reference in time and space for display')
    # reference date
    ref.add_argument('--ref-date', dest='ref_date', metavar='DATE',
                     help='Change reference date for display')

    # reference pixel
    ref.add_argument('--ref-lalo', dest='ref_lalo', metavar=('LAT', 'LON'), type=float, nargs=2,
                     help='Change referene point LAT LON for display')
    ref.add_argument('--ref-yx', dest='ref_yx', metavar=('Y', 'X'), type=int, nargs=2,
                     help='Change referene point Y X for display')

    # reference pixel style
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
                           '\t1) output file already exists\n'+
                           '\t2) output file is newer than input file.')
    return parser


def add_subset_argument(parser):
    """Argument group parser for subset options"""
    sub = parser.add_argument_group('Subset', 'Display dataset in subset range')
    sub.add_argument('--sub-x','--subx', dest='subset_x', type=int, nargs=2, metavar=('XMIN', 'XMAX'),
                     help='subset display in x/cross-track/range direction')
    sub.add_argument('--sub-y','--suby', dest='subset_y', type=int, nargs=2, metavar=('YMIN', 'YMAX'),
                     help='subset display in y/along-track/azimuth direction')
    sub.add_argument('--sub-lat','--sublat', dest='subset_lat', type=float, nargs=2, metavar=('LATMIN', 'LATMAX'),
                     help='subset display in latitude')
    sub.add_argument('--sub-lon','--sublon', dest='subset_lon', type=float, nargs=2, metavar=('LONMIN', 'LONMAX'),
                     help='subset display in longitude')
    return parser

