#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend import:
#   import pysar.view as view


import os
import sys
import argparse
from datetime import datetime as dt

import h5py
import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import cm, pyproj

from pysar.objects import (geometryDatasetNames,
                           geometry,
                           ifgramDatasetNames,
                           ifgramStack,
                           timeseriesDatasetNames,
                           timeseriesKeyNames,
                           timeseries,
                           HDFEOS)
from pysar.utils import (ptime,
                         readfile,
                         utils as ut,
                         plot as pp)
from pysar.multilook import multilook_data
from pysar import subset

fig = None


##################################################################################################
EXAMPLE = """example:
  view.py velocity.h5
  view.py velocity.h5  velocity  --vlim -2 2  -c RdBu
  view.py velocity.h5  --ref-yx  210 566                               #Change reference pixel for display
  view.py velocity.h5  --sub-x 100 600  --sub-y 200 800                #plot subset in yx
  view.py velocity.h5  --sub-lat 31.05 31.10  --sub-lon 130.05 130.10  #plot subset in lalo

  view.py timeseries.h5 
  view.py timeseries.h5 -m no                          #Do not use auto mask
  view.py timeseries.h5 --ref-date 20101120            #Change reference date
  view.py timeseries.h5 --ex drop_date.txt             #Exclude dates to plot

  view.py INPUTS/ifgramStack.h5 coherence
  view.py INPUTS/ifgramStack.h5 unwrapPhase-20070927_20100217 --zero-mask --wrap
  view.py INPUTS/ifgramStack.h5 -n 6
  view.py INPUTS/ifgramStack.h5 20171010_20171115      #Display all data related with one interferometric pair

  view.py GEOCODE/geo_velocity.h5  --pts-file pts.yx   #Plot points with lat/lon in file

  # Save and Output:
  view.py velocity.h5 --save
  view.py velocity.h5 --nodisplay
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


def create_parser():
    parser = argparse.ArgumentParser(description='Plot InSAR Product in 2D',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    infile = parser.add_argument_group('Input File', 'File/Dataset to display')
    infile.add_argument('file', type=str, help='file for display')
    infile.add_argument('dset', type=str, nargs='*', default=[],
                        help='optional - dataset(s) to display')
    infile.add_argument('--exact', '--no-glob', dest='globSearch', action='store_false',
                        help='Disable glob search for input dset')
    infile.add_argument('-n', '--dset-num', dest='dsetNumList', metavar='NUM', type=int, nargs='*', default=[],
                        help='optional - order number of date/dataset(s) to display')
    infile.add_argument('--ex', '--exclude', dest='exDsetList', metavar='Dset', nargs='*', default=[],
                        help='dates will not be displayed')
    infile.add_argument('--plot-setting', dest='disp_setting_file',
                        help='Template file with plot setting.\n'+PLOT_TEMPLATE)

    parser = pp.add_data_disp_argument(parser)
    parser = pp.add_dem_argument(parser)
    parser = pp.add_figure_argument(parser)
    parser = pp.add_gps_argument(parser)
    parser = pp.add_mask_argument(parser)
    parser = pp.add_map_argument(parser)
    parser = pp.add_point_argument(parser)
    parser = pp.add_reference_argument(parser)
    parser = pp.add_save_argument(parser)
    parser = pp.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # If output flie name assigned or figure shown is turned off, turn on the figure save
    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.coastline and inps.resolution in ['c', 'l']:
        inps.resolution = 'i'
    if inps.lalo_step:
        inps.lalo_label = True
    return inps


##################################################################################################
def check_multilook_input(pixel_box, row_num, col_num):
    # Estimate multilook_num
    box_size = (pixel_box[2] - pixel_box[0]) * (pixel_box[3] - pixel_box[1])
    pixel_num_per_figure = box_size * row_num * col_num
    if   pixel_num_per_figure > (8e6*160):   multilook_num=16;      ## 2k * 2k image with 120 subplots
    elif pixel_num_per_figure > (4e6*80) :   multilook_num=8;       ## 2k * 2k image with 80  subplots
    elif pixel_num_per_figure > (4e6*20) :   multilook_num=4;       ## 2k * 2k image with 40  subplots
    elif pixel_num_per_figure > (1e6*20) :   multilook_num=2;       ## 2k * 2k image with 40  subplots
    else: multilook_num=1
    # Update multilook based on multilook_num
    if multilook_num > 1:
        multilook = True
        print('number of data points per figure: '+'%.1E' %
              (pixel_num_per_figure))
        print('multilook with a factor of '+str(multilook_num)+' for display')
    else:
        multilook = False
    return multilook, multilook_num


##################################################################################################
def update_inps_with_display_setting_file(inps, disp_set_file):
    """Update inps using values from display setting file"""
    disp_set_dict = readfile.read_template(disp_set_file)
    if not inps.disp_unit and 'plot.displayUnit' in disp_set_dict.keys():
        inps.disp_unit = disp_set_dict['plot.displayUnit']
    if not inps.disp_min and 'plot.displayMin' in disp_set_dict.keys():
        inps.disp_min = float(disp_set_dict['plot.displayMin'])
    if not inps.disp_max and 'plot.displayMax' in disp_set_dict.keys():
        inps.disp_max = float(disp_set_dict['plot.displayMax'])

    if not inps.colormap and 'plot.colormap' in disp_set_dict.keys():
        inps.colormap = disp_set_dict['plot.colormap']

    if not inps.subset_lat and 'plot.subset.lalo' in disp_set_dict.keys():
        inps.subset_lat = [float(n) for n in disp_set_dict['plot.subset.lalo'].replace(',', ' ').split()[0:2]]
    if not inps.subset_lon and 'plot.subset.lalo' in disp_set_dict.keys():
        inps.subset_lon = [float(n) for n in disp_set_dict['plot.subset.lalo'].replace(',', ' ').split()[2:4]]
    if not inps.ref_lalo and 'plot.seed.lalo' in disp_set_dict.keys():
        inps.ref_lalo = [float(n) for n in disp_set_dict['plot.referenceLalo'].replace(',', ' ').split()]
    return inps


def update_inps_with_file_metadata(inps, metadata, print_msg=True):
    # Subset
    # Convert subset input into bounding box in radar / geo coordinate
    # geo_box = None if atr is not geocoded.
    coord = ut.coordinate(metadata)
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), metadata)
    inps.pix_box = coord.check_box_within_data_coverage(inps.pix_box)
    inps.geo_box = coord.box_pixel2geo(inps.pix_box)
    # Out message
    data_box = (0, 0, inps.width, inps.length)
    if print_msg:
        print('data   coverage in y/x: '+str(data_box))
        print('subset coverage in y/x: '+str(inps.pix_box))
        print('data   coverage in lat/lon: '+str(coord.box_pixel2geo(data_box)))
        print('subset coverage in lat/lon: '+str(inps.geo_box))
        print('------------------------------------------------------------------------')

    # DEM contour display
    if max(inps.pix_box[3] - inps.pix_box[1],
           inps.pix_box[2] - inps.pix_box[0]) > 2e3:
        inps.disp_dem_contour = False
        if print_msg:
            print('area exceed 2000 pixels, turn off default DEM contour display')

    # Multilook, if too many subplots in one figure for less memory and faster speed
    if inps.multilook_num > 1:
        inps.multilook = True

    # Colormap
    inps.colormap = pp.check_colormap_input(metadata,
                                            inps.colormap,
                                            datasetName=inps.dset[0],
                                            print_msg=print_msg)

    # Reference Point
    # Convert ref_lalo if existed, to ref_yx, and use ref_yx for the following
    # ref_yx is referenced to input data coverage, not subseted area for display
    if inps.ref_lalo and inps.geo_box:
        inps.ref_yx = [coord.lalo2yx(inps.ref_lalo[0], coord_type='lat'),
                       coord.lalo2yx(inps.ref_lalo[1], coord_type='lon')]
        if print_msg:
            print('input reference point in lat/lon: {}'.format(inps.ref_lalo))
            print('input reference point in y  /x  : {}'.format(inps.ref_yx))

    # ref_lalo
    if inps.ref_yx and inps.geo_box:
        inps.ref_lalo = [coord.yx2lalo(inps.ref_yx[0], coord_type='y'),
                         coord.yx2lalo(inps.ref_yx[1], coord_type='x')]
    elif 'REF_LAT' in metadata.keys():
        inps.ref_lalo = [float(metadata['REF_LAT']),
                         float(metadata['REF_LON'])]
    else:
        inps.ref_lalo = None

    # Points of interest
    inps = pp.read_point2inps(inps, coord)

    # Unit and Wrap
    inps.disp_unit, inps.wrap = pp.check_disp_unit_and_wrap(metadata,
                                                            disp_unit=inps.disp_unit,
                                                            wrap=inps.wrap)

    # Min / Max - Display
    if not inps.vlim:
        if (any(i in inps.key.lower() for i in ['coherence', '.cor'])
                or (inps.key == 'ifgramStack'
                        and inps.dset[0].split('-')[0] in ['coherence', 'connectComponent'])):
            inps.vlim = [0.0, 1.0]
        elif inps.key in ['.int'] or inps.wrap:
            inps.vlim = [-np.pi, np.pi]

    # Transparency - Alpha
    if not inps.transparency:
        # Auto adjust transparency value when showing shaded relief DEM
        if inps.dem_file and inps.disp_dem_shade:
            inps.transparency = 0.8
        else:
            inps.transparency = 1.0

    # Flip Left-Right / Up-Down
    if not inps.flip_lr and not inps.flip_ud:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(metadata)

    # Figure Title
    if not inps.fig_title:
        inps.fig_title = pp.auto_figure_title(metadata['FILE_PATH'],
                                              datasetNames=inps.dset,
                                              inps_dict=vars(inps))
    if print_msg:
        print('figure title: '+inps.fig_title)

    # Figure output file name
    if not inps.outfile:
        inps.outfile = inps.fig_title+inps.fig_ext

    inps = update_figure_setting(inps, print_msg=print_msg)
    return inps


##################################################################################################
def update_data_with_plot_inps(data, metadata, inps, print_msg=True):
    # Seed Point
    if inps.ref_yx:   # and inps.ref_yx != [int(metadata['REF_Y']), int(metadata['REF_X'])]:
        try:
            ref_y = inps.ref_yx[0] - inps.pix_box[1]
            ref_x = inps.ref_yx[1] - inps.pix_box[0]
        except:
            pass
        if len(data.shape) == 2:
            data -= data[ref_y, ref_x]
        elif len(data.shape) == 3:
            data -= np.tile(data[:, ref_y, ref_x].reshape(-1, 1, 1),
                            (1, data.shape[1], data.shape[2]))
        if print_msg:
            print('set reference pixel to: {}'.format(inps.ref_yx))
    else:
        inps.ref_yx = None

    # Convert data to display unit and wrap
    (data,
     inps.disp_unit,
     inps.disp_scale,
     inps.wrap) = pp.scale_data4disp_unit_and_rewrap(data,
                                                     metadata=metadata,
                                                     disp_unit=inps.disp_unit,
                                                     wrap=inps.wrap,
                                                     print_msg=print_msg)
    if inps.wrap:
        inps.vlim = [-np.pi, np.pi]

    # 1.6 Min / Max - Data/Display
    inps.dlim = [np.nanmin(data), np.nanmax(data)]
    if not inps.vlim:
        inps.vlim = [np.nanmin(data), np.nanmax(data)]
    if print_msg:
        print('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
        print('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    return data, inps


##################################################################################################
def plot_2d_matrix(ax, data, metadata, inps=None, print_msg=True):
    """Plot 2D matrix 
    Parameters: ax : matplot.pyplot axes object
                data : 2D np.array, 
                metadata : dictionary, attributes of data
                inps : Namespace, optional, input options for display
    Returns:    ax  : matplot.pyplot axes object
    Example:    import matplotlib.pyplot as plt
                import pysar.utils.readfile as readfile
                import pysar.view as pv
                data, atr = readfile.read('velocity.h5')
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                ax = pv.plot_2d_matrix(ax, data, atr)[0]
                plt.show()
    """

    #----------------------- 0. Initial a inps Namespace if no inps input --------------------#
    if not inps:
        inps = cmd_line_parse([''])
        inps = update_inps_with_file_metadata(inps, metadata)

    #----------------------- 1. Update plot inps/data with data matrix -----------------------#
    #data, inps = update_data_with_plot_inps(data, metadata, inps)

    # 1.7 DEM 
    if inps.dem_file:
        dem = pp.read_dem(inps.dem_file,
                          pix_box=inps.pix_box,
                          geo_box=inps.geo_box,
                          print_msg=print_msg)

    if print_msg:
        print('display data in transparency: '+str(inps.transparency))

    #-------------------- 2.1 Plot in Geo-coordinate using Basemap --------------------------------#
    num_row, num_col = data.shape
    if inps.geo_box and inps.fig_coord == 'geo':
        # Map Setup
        if print_msg:
            print('plot in Lat/Lon coordinate')
            print('map projection: '+inps.map_projection)
            print('boundary database resolution: '+inps.resolution)
        if inps.map_projection in ['cyl', 'merc', 'mill', 'cea', 'gall']:
            m = pp.BasemapExt(llcrnrlon=inps.geo_box[0], llcrnrlat=inps.geo_box[3],
                              urcrnrlon=inps.geo_box[2], urcrnrlat=inps.geo_box[1],
                              projection=inps.map_projection,
                              resolution=inps.resolution, area_thresh=1.,
                              suppress_ticks=False, ax=ax)
        elif inps.map_projection in ['ortho']:
            m = pp.BasemapExt(lon_0=(inps.geo_box[0]+inps.geo_box[2])/2.0,
                              lat_0=(inps.geo_box[3]+inps.geo_box[1])/2.0,
                              projection=inps.map_projection,
                              resolution=inps.resolution, area_thresh=1.,
                              suppress_ticks=False, ax=ax)
        else:
            m = pp.BasemapExt(lon_0=(inps.geo_box[0]+inps.geo_box[2])/2.0,
                              lat_0=(inps.geo_box[3]+inps.geo_box[1])/2.0,
                              llcrnrlon=inps.geo_box[0], llcrnrlat=inps.geo_box[3],
                              urcrnrlon=inps.geo_box[2], urcrnrlat=inps.geo_box[1],
                              projection=inps.map_projection,
                              resolution=inps.resolution, area_thresh=1.,
                              suppress_ticks=False, ax=ax)

        # Draw coastline
        if inps.coastline:
            if print_msg:
                print('draw coast line')
            m.drawcoastlines()

        # Plot DEM
        if inps.dem_file:
            if print_msg:
                print('plotting DEM background ...')
            m = pp.plot_dem_background(ax=m, geo_box=inps.geo_box,
                                       dem=dem, inps_dict=vars(inps),
                                       print_msg=print_msg)

        # Plot Data
        if print_msg:
            print('plotting Data ...')
        im = m.imshow(data, cmap=inps.colormap, origin='upper',
                      vmin=inps.vlim[0], vmax=inps.vlim[1],
                      alpha=inps.transparency, interpolation='nearest',
                      animated=inps.animation)

        # Scale Bar
        if inps.disp_scalebar:
            if print_msg:
                print('plot scale bar')
            if not inps.scalebar:
                inps.scalebar = [999, 999, 999]

            # Default Distance - 20% of data width
            if inps.scalebar[0] == 999.0:
                gc = pyproj.Geod(a=m.rmajor, b=m.rminor)
                wid_dist = gc.inv(inps.geo_box[0], inps.geo_box[3],
                                  inps.geo_box[2], inps.geo_box[3])[2]
                inps.scalebar[0] = ut.round_to_1(wid_dist * 0.2)

            # Default center - Lower Left Corner
            if inps.scalebar[1] == 999.0:
                inps.scalebar[1] = inps.geo_box[3] + 0.1 * (inps.geo_box[1] - inps.geo_box[3])
            if inps.scalebar[2] == 999.0:
                inps.scalebar[2] = inps.geo_box[0] + 0.2 * (inps.geo_box[2] - inps.geo_box[0])

            # Draw scale bar
            m.draw_scale_bar(inps.scalebar[1], inps.scalebar[2], inps.scalebar[0], ax=ax,
                             font_size=inps.font_size, color=inps.font_color)

        # Lat Lon labels
        if inps.lalo_label:
            if print_msg:
                print('plot lat/lon labels')
            m.draw_lalo_label(inps.geo_box, ax=ax, lalo_step=inps.lalo_step, labels=inps.lalo_label_loc,
                              font_size=inps.font_size, color=inps.font_color, print_msg=print_msg)
        else:
            ax.tick_params(labelsize=inps.font_size, colors=inps.font_color)

        # Plot Reference Point
        if inps.disp_ref_pixel and inps.ref_lalo:
            ax.plot(inps.ref_lalo[1], inps.ref_lalo[0],
                    inps.ref_marker, ms=inps.ref_size)
            if print_msg:
                print('plot reference point')

        # Plot points of interest
        if inps.pts_lalo is not None:
            ax.plot(inps.pts_lalo[:, 1], inps.pts_lalo[:, 0],
                    inps.pts_marker, ms=inps.ref_size,
                    markeredgecolor='black')
            if print_msg:
                print('plot points of interest')

        # Show UNR GPS stations
        if inps.disp_gps:
            from pysar.objects import gps
            SNWE = (inps.geo_box[3], inps.geo_box[1], inps.geo_box[0], inps.geo_box[2])
            site_names, site_lats, site_lons = gps.search_gps(SNWE)
            ax.scatter(site_lons, site_lats, s=7**2, color='w', edgecolors='k')
            for i in range(len(site_names)):
                ax.annotate(site_names[i], xy=(site_lons[i], site_lats[i]), fontsize=inps.font_size)
            if print_msg:
                print('displaying GPS stations')

        # Status bar
        coord = ut.coordinate(metadata)
        def format_coord(x, y):
            col = coord.lalo2yx(x, coord_type='lon') - inps.pix_box[0]
            row = coord.lalo2yx(y, coord_type='lat') - inps.pix_box[1]
            msg = 'Lon={:.4f}, Lat={:.4f}'.format(x, y)
            if 0 <= col < num_col and 0 <= row < num_row:
                try:
                    z = data[row, col]
                    msg += ', value={:.4f}'.format(z)
                except:
                    msg += ', value=[]'
                if inps.dem_file:
                    dem_col = coord.lalo2yx(x, coord_type='lon') - inps.dem_pix_box[0]
                    dem_row = coord.lalo2yx(y, coord_type='lat') - inps.dem_pix_box[1]
                    h = dem[dem_row, dem_col]
                    msg += ', elev={:.1f}'.format(h)
                msg += ', x={:.1f}, y={:.1f}'.format(col+inps.pix_box[0],
                                                     row+inps.pix_box[1])
            return msg
        ax.format_coord = format_coord

    #-------------------- 2.2 Plot in Y/X-coordinate ------------------------------------------------#
    else:
        inps.fig_coord = 'radar'
        if print_msg:
            print('plotting in Y/X coordinate ...')
        # Plot DEM
        if inps.dem_file:
            if print_msg:
                print('plotting DEM background ...')
            ax = pp.plot_dem_background(ax=ax,
                                        geo_box=None,
                                        dem=dem,
                                        inps_dict=vars(inps),
                                        print_msg=print_msg)

        # Plot Data
        if print_msg:
            print('plotting Data ...')
        im = ax.imshow(data, cmap=inps.colormap,
                       vmin=inps.vlim[0], vmax=inps.vlim[1],
                       alpha=inps.transparency, interpolation='nearest')
        ax.tick_params(labelsize=inps.font_size)

        # Plot Seed Point
        if inps.disp_ref_pixel:
            ref_y, ref_x = None, None
            if inps.ref_yx:
                ref_y, ref_x = inps.ref_yx[0], inps.ref_yx[1]
            elif 'REF_Y' in metadata.keys():
                ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])

            if ref_y and ref_x:            
                ax.plot(ref_x - inps.pix_box[0],
                        ref_y - inps.pix_box[1],
                        inps.ref_marker, ms=inps.ref_size)
                if print_msg:
                    print('plot reference point')

        # Plot points of interest
        if inps.pts_yx is not None:
            ax.plot(inps.pts_yx[0, :], inps.pts_yx[1, :],
                    inps.pts_marker, ms=inps.ref_size,
                    markeredgecolor='black')
            if print_msg:
                print('plot points of interest')

        ax.set_xlim(-0.5, num_col-0.5)
        ax.set_ylim(num_row-0.5, -0.5)

        # Status bar
        def format_coord(x, y):
            col = int(x+0.5)
            row = int(y+0.5)
            msg = 'x={:.1f}, y={:.1f}'.format(x, y)
            if 0 <= col < num_col and 0 <= row < num_row:
                try:
                    z = data[row, col]
                    msg += ', value={:.4f}'.format(z)
                except:
                    msg += ', value=[]'
                if inps.dem_file:
                    h = dem[row, col]
                    msg += ', elev={:.1f} m'.format(h)
            return msg
        ax.format_coord = format_coord

    #-------------------- 3 Figure Setting --------------------------------------------------------#
    # 3.1 Colorbar
    cbar = None
    if inps.disp_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "2%", pad="2%")
        inps, cbar = pp.plot_colorbar(inps, im, cax)

    # 3.2 Title
    if inps.disp_title:
        ax.set_title(inps.fig_title, fontsize=inps.font_size,
                     color=inps.font_color)

    # 3.3 Flip Left-Right / Up-Down
    if inps.flip_lr:
        if print_msg:
            print('flip figure left and right')
        ax.invert_xaxis()

    if inps.flip_ud:
        if print_msg:
            print('flip figure up and down')
        ax.invert_yaxis()

    # 3.4 Turn off axis
    if not inps.disp_axis:
        ax.axis('off')
        if print_msg:
            print('turn off axis display')

    # 3.5 Turn off tick label
    if not inps.disp_tick:
        # ax.set_xticklabels([])
        # ax.set_yticklabels([])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    # Figure Output
    if inps.save_fig:
        if print_msg:
            print('save figure to '+inps.outfile)
        plt.savefig(inps.outfile, bbox_inches='tight',
                    transparent=True, dpi=inps.fig_dpi)

    return ax, inps, im, cbar


def check_input_file_info(inps, print_msg=True):
    # File Baic Info
    atr = readfile.read_attribute(inps.file)
    msg = 'input file is '
    if not inps.file.endswith(('.h5', '.he5')):
        msg += '{} '.format(atr['PROCESSOR'])
    msg += '{} file: {}'.format(atr['FILE_TYPE'], inps.file)
    if 'DATA_TYPE' in atr.keys():
        msg += ' in {} format'.format(atr['DATA_TYPE'])
    if print_msg:
        print('\n******************** Display ********************')
        print(msg)

    ## size and name
    inps.length = int(atr['LENGTH'])
    inps.width = int(atr['WIDTH'])
    inps.key = atr['FILE_TYPE']
    inps.fileBase = os.path.splitext(os.path.basename(inps.file))[0]
    inps.fileExt = os.path.splitext(inps.file)[1]
    if print_msg:
        print('file size in y/x: {}'.format((inps.length, inps.width)))

    # File dataset List
    inps.fileDatasetList = readfile.get_2d_dataset_list(inps.file)

    # Read input list of dataset to display
    inps, atr = read_dataset_input(inps, print_msg=print_msg)

    return inps, atr


def check_dataset_input(allList, inList=[], inNumList=[], globSearch=True):
    """Get dataset(es) from input dataset / dataset_num"""
    # inList + inNumList --> outNumList --> outList
    if inList:
        if isinstance(inList, str):
            inList = [inList]
        tempList = []
        if globSearch:
            for i in inList:
                tempList += [e for e in allList if i in e]
        else:
            tempList += [i for i in inList if i in allList]
        tempList = sorted(list(set(tempList)))
        inNumList += [allList.index(e) for e in tempList]
    outNumList = sorted(list(set(inNumList)))
    outList = [allList[i] for i in outNumList]
    return outList, outNumList


def read_dataset_input(inps, print_msg=True):
    """Check input / exclude / reference dataset input with file dataset list"""
    # read inps.dset + inps.dsetNumList --> inps.dsetNumList
    if len(inps.dset) > 0 or len(inps.dsetNumList) > 0:
        if inps.key == 'velocity':
            inps.globSearch = False
            if print_msg:
                print('turning glob search OFF for {} file'.format(inps.key))
        inps.dsetNumList = check_dataset_input(inps.fileDatasetList,
                                               inps.dset,
                                               inps.dsetNumList,
                                               inps.globSearch)[1]
    elif inps.key in ['geometry', 'ifgramStack', 'HDFEOS']:
        # default dataset to display for certain type of files
        if inps.key == 'geometry':
            inps.dset = geometryDatasetNames
            inps.dset.remove('bperp')
        elif inps.key == 'ifgramStack':
            inps.dset = ['unwrapPhase']
        elif inps.key == 'HDFEOS':
            inps.dset = ['displacement']
        inps.dsetNumList = check_dataset_input(inps.fileDatasetList,
                                               inps.dset,
                                               inps.dsetNumList,
                                               inps.globSearch)[1]
    else:
        inps.dsetNumList = range(len(inps.fileDatasetList))

    # read inps.exDsetList
    inps.exDsetList, inps.exDsetNumList = check_dataset_input(inps.fileDatasetList,
                                                              inps.exDsetList,
                                                              [],
                                                              inps.globSearch)

    # get inps.dset
    inps.dsetNumList = sorted(list(set(inps.dsetNumList) - set(inps.exDsetNumList)))
    inps.dset = [inps.fileDatasetList[i] for i in inps.dsetNumList]
    inps.dsetNum = len(inps.dset)

    if inps.ref_date:
        if inps.key not in timeseriesKeyNames:
            inps.ref_date = None
        ref_date = check_dataset_input(inps.fileDatasetList,
                                       [inps.ref_date],
                                       [],
                                       inps.globSearch)[0][0]
        if not ref_date:
            if print_msg:
                print('WARNING: input reference date is not included in input file!')
                print('input reference date: '+inps.ref_date)
            inps.ref_date = None
        else:
            inps.ref_date = ref_date

    if print_msg:
        if inps.key in ['ifgramStack']:
            print('num of datasets in file {}: {}'.format(os.path.basename(inps.file),
                                                          len(inps.fileDatasetList)))
            print('num of datasets to exclude: {}'.format(len(inps.exDsetList)))
            print('num of datasets to display: {}'.format(len(inps.dset)))
        else:
            print('num of datasets in file {}: {}'.format(os.path.basename(inps.file),
                                                          len(inps.fileDatasetList)))
            print('datasets to exclude ({}):\n{}'.format(len(inps.exDsetList),
                                                         inps.exDsetList))
            print('datasets to display ({}):\n{}'.format(len(inps.dset),
                                                         inps.dset))
        if inps.ref_date and inps.key in timeseriesKeyNames:
            print('input reference date: {}'.format(inps.ref_date))

    if inps.dsetNum == 0:
        print('ERROR: no input dataset found!')
        print('available datasets:\n{}'.format(inps.fileDatasetList))
        raise Exception()

    atr = readfile.read_attribute(inps.file, datasetName=inps.dset[0].split('-')[0])
    return inps, atr


def update_figure_setting(inps, print_msg=True):
    """Update figure setting based on number of subplots/datasets
    1) fig_size and font_size
    2) for multi: figure/row/column number
    3) for multi: output file name
    """
    length = float(inps.pix_box[3]-inps.pix_box[1])
    width = float(inps.pix_box[2]-inps.pix_box[0])

    # One Plot
    if inps.dsetNum == 1:
        if not inps.font_size:
            inps.font_size = 16
        if not inps.fig_size:
            if inps.geo_box and inps.fig_coord == 'geo':
                length = abs(inps.geo_box[3] - inps.geo_box[1])
                width = abs(inps.geo_box[2] - inps.geo_box[0])
                plot_shape = []
            plot_shape = [width*1.25, length]
            fig_scale = min(pp.min_figsize_single/min(plot_shape),
                            pp.max_figsize_single/max(plot_shape),
                            pp.max_figsize_height/plot_shape[1])
            inps.fig_size = [np.floor(i*fig_scale*2)/2 for i in plot_shape]
            if print_msg:
                print('create figure in size: '+str(inps.fig_size))

    # Multiple Plots
    else:
        if not inps.fig_size:
            inps.fig_size = pp.default_figsize_multi
            if print_msg:
                print('create figure in size: '+str(inps.fig_size))

        # Figure number (<= 200 subplots per figure)
        if not inps.fig_num:
            inps.fig_num = 1
            while inps.dsetNum/float(inps.fig_num) > 200.0:
                inps.fig_num += 1

        # Row/Column number
        if inps.fig_row_num == 1 and inps.fig_col_num == 1:
            # calculate row and col number based on input info
            data_shape = [length*1.1, width]
            fig_size4plot = [inps.fig_size[0]*0.95, inps.fig_size[1]]
            (inps.fig_row_num,
             inps.fig_col_num) = pp.auto_row_col_num(inps.dsetNum,
                                                     data_shape,
                                                     fig_size4plot,
                                                     inps.fig_num)
        inps.fig_num = np.ceil(float(inps.dsetNum) / float(inps.fig_row_num * 
                                                           inps.fig_col_num)).astype(int)
        if print_msg:
            print('dataset number: '+str(inps.dsetNum))
            print('row     number: '+str(inps.fig_row_num))
            print('column  number: '+str(inps.fig_col_num))
            print('figure  number: '+str(inps.fig_num))

        if not inps.font_size:
            inps.font_size = 12
            if inps.fig_row_num * inps.fig_col_num > 50:
                inps.font_size = 8

        # Output File Name
        if inps.outfile:
            inps.fig_ext = os.path.splitext(inps.outfile)[1].lower()
            inps.outfile_base = os.path.basename(inps.outfile).split(inps.fig_ext)[0]
        else:
            inps.outfile_base = os.path.splitext(inps.file)[0]

            if (inps.pix_box[2]-inps.pix_box[0])*(inps.pix_box[3]-inps.pix_box[1]) < width*length:
                inps.outfile_base += '_sub'
            if inps.wrap:
                inps.outfile_base += '_wrap'
            if inps.ref_date:
                inps.outfile_base += '_ref'+inps.ref_date
            if inps.exDsetList:
                inps.outfile_base += '_ex'
    return inps


def read_data4figure(i_start, i_end, inps, metadata):
    """Read multiple datasets for one figure into 3D matrix based on i_start/end"""
    data = np.zeros((i_end - i_start,
                     inps.pix_box[3] - inps.pix_box[1],
                     inps.pix_box[2] - inps.pix_box[0]))

    # fast reading for single dataset type
    if (len(inps.dsetFamilyList) == 1
            and inps.key in ['timeseries', 'giantTimeseries', 'ifgramStack', 'HDFEOS', 'geometry']):
        dset_list = [inps.dset[i] for i in range(i_start, i_end)]
        data = readfile.read(inps.file, datasetName=dset_list, box=inps.pix_box)[0]

        if inps.key == 'ifgramStack':
            # reference pixel info in unwrapPhase
            if inps.dsetFamilyList[0] == 'unwrapPhase' and inps.file_ref_yx:
                for i in range(data.shape[0]):
                    mask = data[i, :, :] != 0.
                    data[i, mask] -= data[i, inps.file_ref_yx[0], inps.file_ref_yx[1]]

    # slow reading with one 2D matrix at a time
    else:
        print('reading data ...')
        prog_bar = ptime.progressBar(maxValue=i_end-i_start)
        for i in range(i_start, i_end):
            d = readfile.read(inps.file,
                              datasetName=inps.dset[i],
                              box=inps.pix_box,
                              print_msg=False)[0]
            data[i - i_start, :, :] = d
            prog_bar.update(i - i_start + 1, suffix=inps.dset[i])
        prog_bar.close()

    # ref_date for timeseries
    if inps.ref_date:
        print('consider input reference date: '+inps.ref_date)
        ref_data = readfile.read(inps.file,
                                 datasetName=inps.ref_date,
                                 box=inps.pix_box,
                                 print_msg=False)[0]
        data -= ref_data

    # v/dlim, adjust data if all subplots are the same type
    #metadata = readfile.read_attribute(inps.file)
    if len(inps.dsetFamilyList) == 1 or inps.key in ['velocity']:
        data, inps = update_data_with_plot_inps(data, metadata, inps, print_msg=False)
        if (not inps.vlim 
                and not (inps.dsetFamilyList[0].startswith('unwrap') and not inps.file_ref_yx)
                and inps.dsetFamilyList[0] not in ['bperp']):
            data_mli = multilook_data(data, 10, 10)
            inps.vlim = [np.nanmin(data_mli), np.nanmax(data_mli)]
            del data_mli
    inps.dlim = [np.nanmin(data), np.nanmax(data)]

    # multilook
    if inps.multilook:
        data = multilook_data(data, inps.multilook_num, inps.multilook_num)

    # mask
    if inps.msk is not None:
        print('masking data')
        msk = np.tile(inps.msk, (data.shape[0], 1, 1))
        data = np.ma.masked_where(msk == 0., data)
    if inps.zero_mask:
        print('masking pixels with zero value')
        data = np.ma.masked_where(data == 0., data)
    return data


def plot_subplot4figure(i, inps, ax, data, metadata):
    """Plot one subplot for one 3D array
    1) Plot DEM, data and reference pixel
    2) axes setting: tick, ticklabel, title, axis etc.
    """
    # Plot DEM
    if inps.dem_file:
        ax = pp.plot_dem_background(ax=ax, geo_box=None, dem_shade=inps.dem_shade,
                                    dem_contour=inps.dem_contour,
                                    dem_contour_seq=inps.dem_contour_seq)
    # Plot Data
    try:
        im = ax.imshow(data, cmap=inps.colormap, interpolation='nearest',
                       alpha=inps.transparency,
                       vmin=inps.vlim[0], vmax=inps.vlim[1])
    except:
        im = ax.imshow(data, cmap=inps.colormap, interpolation='nearest',
                       alpha=inps.transparency)
    # Plot Seed Point
    if inps.disp_ref_pixel:
        ref_y, ref_x = None, None
        if inps.ref_yx:
            ref_y, ref_x = inps.ref_yx[0], inps.ref_yx[1]
        elif 'REF_Y' in metadata.keys():
            ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])

        if ref_y and ref_x:
            ax.plot(ref_x - inps.pix_box[0],
                    ref_y - inps.pix_box[1],
                    inps.ref_marker, ms=inps.ref_size)

    ax.set_xlim(-0.5, np.shape(data)[1]-0.5)
    ax.set_ylim(np.shape(data)[0]-0.5, -0.5)

    # Subplot Setting
    ## Tick and Label
    #ax.set_yticklabels([])
    #ax.set_xticklabels([])
    #ax.set_xticks([])
    #ax.set_yticks([])
    if not inps.disp_tick or inps.fig_row_num * inps.fig_col_num > 10:
        # ax.set_xticklabels([])
        # ax.set_yticklabels([])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    # status bar
    def format_coord(x, y):
        row = int(y + 0.5)
        col = int(x + 0.5)
        return 'x={:.1f}, y={:.1f}, value={:.4f}'.format(x, y, data[row, col])
    ax.format_coord = format_coord

    # Title
    if inps.disp_title:
        # get title
        if inps.key in timeseriesKeyNames or inps.dset[0].startswith('bperp'):
            try:
                subplot_title = dt.strptime(inps.dset[i].split('-')[1], '%Y%m%d').isoformat()[0:10]
            except:
                subplot_title = str(inps.dset[i])
        elif inps.key in ['ifgramStack', 'interferograms', 'coherence', 'wrapped']:
            subplot_title = str(i)
            if inps.fig_row_num * inps.fig_col_num < 50:
                subplot_title += '\n{}'.format(inps.dset[i])
            elif inps.fig_row_num * inps.fig_col_num > 200:
                subplot_title = ''
        else:
            subplot_title = str(inps.dset[i])

        # plot title
        if not inps.fig_title_in:
            if inps.dset[i] in inps.dropDatasetList:
                ax.set_title(subplot_title, fontsize=inps.font_size,
                             color='crimson', fontweight='bold')
            else:
                ax.set_title(subplot_title, fontsize=inps.font_size)
        else:
            pp.add_inner_title(ax, subplot_title, loc=1)

    # Flip Left-Right / Up-Down
    if inps.flip_lr:
        ax.invert_xaxis()
    if inps.flip_ud:
        ax.invert_yaxis()

    # Turn off axis
    if not inps.disp_axis:
        ax.axis('off')
    return im


def plot_figure(j, inps, metadata):
    global fig
    """Plot one figure with multiple subplots
    1) create figure
    2) read all data into 3D array
    3) loop to plot each subplot using plot_subplot4figure()
    4) common colorbar and save
    """
    # Output file name for current figure
    if inps.fig_num > 1:
        inps.outfile = '{}_{}{}'.format(inps.outfile_base, str(j), inps.fig_ext)
    else:
        inps.outfile = '{}{}'.format(inps.outfile_base, inps.fig_ext)
    fig_title = 'Figure {} - {}'.format(str(j), inps.outfile)
    print('----------------------------------------')
    print(fig_title)

    # Open a new figure object
    fig = plt.figure(j, figsize=inps.fig_size)
    fig.canvas.set_window_title(fig_title)

    # Read all data for the current figure into 3D np.array
    i_start = (j - 1) * inps.fig_row_num * inps.fig_col_num
    i_end = min([inps.dsetNum, i_start + inps.fig_row_num * inps.fig_col_num])
    data = read_data4figure(i_start, i_end, inps, metadata)

    # Loop - Subplots
    print('plotting ...')
    prog_bar = ptime.progressBar(maxValue=i_end - i_start)
    for i in range(i_start, i_end):
        idx = i - i_start
        ax = fig.add_subplot(inps.fig_row_num, inps.fig_col_num, idx + 1)
        im = plot_subplot4figure(i, inps, ax=ax,
                                 data=data[idx, :, :],
                                 metadata=metadata)
        prog_bar.update(idx+1, suffix=inps.dset[i])
    prog_bar.close()
    del data
    fig.tight_layout()

    # Min and Max for this figure
    inps.dlim_all = [np.nanmin([inps.dlim_all[0], inps.dlim[0]]), 
                     np.nanmax([inps.dlim_all[1], inps.dlim[1]])]
    print('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
    if inps.vlim:
        print('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    # Colorbar
    if not inps.vlim:
        print('Note: different color scale for EACH subplot!')
    else:
        print('show colorbar')
        fig.subplots_adjust(right=0.93)
        cax = fig.add_axes([0.94, 0.3, 0.005, 0.4])
        inps, cax = pp.plot_colorbar(inps, im, cax)

    # Save Figure
    if inps.save_fig:
        print('save figure to '+inps.outfile)
        fig.savefig(inps.outfile, bbox_inches='tight',
                    transparent=True, dpi=inps.fig_dpi)
        if not inps.disp_fig:
            fig.clf()
    return


def prepare4multi_subplots(inps, metadata):
    """Prepare for multiple subplots:
    1) check multilook to save memory
    2) read existed reference pixel info for unwrapPhase
    3) read dropIfgram info
    4) read and prepare DEM for background
    """
    inps.dsetFamilyList = list(set(i.split('-')[0] for i in inps.dset))

    # Update multilook parameters with new num and col number
    if inps.multilook and inps.multilook_num == 1:
        # Do not auto multilook mask and lookup table file
        auto_multilook = True
        for dsFamily in inps.dsetFamilyList:
            if any(i in dsFamily.lower() for i in ['mask', 'coord']):
                auto_multilook = False
        if auto_multilook:
            inps.multilook, inps.multilook_num = check_multilook_input(inps.pix_box,
                                                                       inps.fig_row_num,
                                                                       inps.fig_col_num)
        if inps.msk is not None:
            inps.msk = multilook_data(inps.msk, inps.multilook_num, inps.multilook_num)

    # Reference pixel for timeseries and ifgramStack
    #metadata = readfile.read_attribute(inps.file)
    inps.file_ref_yx = None
    if inps.key in ['ifgramStack'] and 'REF_Y' in metadata.keys():
        ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])
        length, width = int(metadata['LENGTH']), int(metadata['WIDTH'])
        if 0 <= ref_y < length and 0 <= ref_x < width:
            inps.file_ref_yx = [ref_y, ref_x]
            print('consider reference pixel in y/x: {}'.format(inps.file_ref_yx))

    if inps.dsetNum > 10:
        inps.disp_ref_pixel = False
        print('turn off reference pixel plot for more than 10 datasets to display')

    # Check dropped interferograms
    inps.dropDatasetList = []
    if inps.key == 'ifgramStack' and inps.disp_title:
        obj = ifgramStack(inps.file)
        obj.open(print_msg=False)
        dropDate12List = obj.get_drop_date12_list()
        for i in inps.dsetFamilyList:
            inps.dropDatasetList += ['{}-{}'.format(i, j) for j in dropDate12List]
        print("mark interferograms with 'dropIfgram=False' in red colored title")

    # Read DEM
    if inps.dem_file:
        dem_metadata = readfile.read_attribute(inps.dem_file)
        if all(dem_metadata[i] == metadata[i] for i in ['LENGTH', 'WIDTH']):
            print('reading DEM: {} ... '.format(os.path.basename(inps.dem_file)))
            dem = readfile.read(inps.dem_file,
                                datasetName='height',
                                box=inps.pix_box,
                                print_msg=False)[0]
            if inps.multilook:
                dem = multilook_data(dem, inps.multilook_num, inps.multilook_num)
            (inps.dem_shade,
             inps.dem_contour,
             inps.dem_contour_seq) = pp.prepare_dem_background(dem=dem, inps_dict=vars(inps))
        else:
            inps.dem_file = None
            inps.transparency = 1.0
            print('Input DEM file has different size than data file, ignore it.')
    return inps


#########################################  Main Function  ########################################
def main(iargs=None):
    global fig

    inps = cmd_line_parse(iargs)
    if not inps.disp_fig:
        plt.switch_backend('Agg')  # Backend setting

    inps, atr = check_input_file_info(inps)

    if inps.disp_setting_file:
        inps = update_inps_with_display_setting_file(inps, inps.disp_setting_file)

    inps = update_inps_with_file_metadata(inps, atr)

    inps.msk, inps.mask_file = pp.read_mask(inps.file,
                                            mask_file=inps.mask_file,
                                            datasetName=inps.dset[0],
                                            box=inps.pix_box)

    ############################### One Subplot ###############################
    if inps.dsetNum == 1:
        print('reading data ...')
        data, atr = readfile.read(inps.file,
                                  datasetName=inps.dset[0],
                                  box=inps.pix_box,
                                  print_msg=False)
        if inps.ref_date:
            data -= readfile.read(inps.file,
                                  datasetName=inps.ref_date,
                                  box=inps.pix_box,
                                  print_msg=False)[0]
        if inps.zero_mask:
            print('masking pixels with zero value')
            data = np.ma.masked_where(data == 0., data)
        if inps.msk is not None:
            print('masking data')
            data = np.ma.masked_where(inps.msk == 0., data)

        data, inps = update_data_with_plot_inps(data, atr, inps)

        fig, ax = plt.subplots(figsize=inps.fig_size, num='Figure')

        ax, inps = plot_2d_matrix(ax, data, atr, inps)[0:2]

    ############################### Multiple Subplots #########################
    else:
        inps = prepare4multi_subplots(inps, metadata=atr)

        inps.dlim_all = [0., 0.]
        for j in range(1, inps.fig_num + 1):
            plot_figure(j, inps, metadata=atr)

        if inps.fig_num > 1:
            print('----------------------------------------')
            print('all data range: {} {}'.format(inps.dlim_all, inps.disp_unit))
            if inps.vlim:
                print('display  range: {} {}'.format(inps.vlim, inps.disp_unit))

    # Display Figure
    if inps.disp_fig:
        print('showing ...')
        plt.show()
    return


##################################################################################################
if __name__ == '__main__':
    main()
