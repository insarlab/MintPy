############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.view import prep_slice, plot_slice, viewer


import datetime as dt
import os
import re
import warnings  # suppress UserWarning from matplotlib

warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mintpy import subset, version
from mintpy.multilook import multilook_data
from mintpy.objects import TIMESERIES_KEY_NAMES, giantIfgramStack, ifgramStack
from mintpy.objects.gps import GPS
from mintpy.utils import plot as pp, ptime, readfile, utils as ut


##################################################################################################
def run_or_skip(inps):
    vprint('update mode: ON')
    flag = 'skip'

    # run if showing figures
    if inps.disp_fig:
        flag = 'run'
        return flag

    # get existed output file names
    outfiles = []
    for fname in inps.outfile:
        fnames = [fname, os.path.join(os.path.dirname(fname), 'pic', os.path.basename(fname))]
        fnames = [i for i in fnames if os.path.isfile(i)]
        if len(fnames) > 0:
            outfiles.append(fnames[0])
        else:
            flag = 'run'

    if flag == 'skip':
        ti = os.path.getmtime(inps.file)
        to = min(os.path.getmtime(i) for i in outfiles)
        if ti > to:
            flag = 'run'
        else:
            vprint(f'{outfiles} exist and are newer than input file: {inps.file} --> skip.')

    return flag


##################################################################################################
def update_inps_with_file_metadata(inps, metadata):
    # Subset
    # Convert subset input into bounding box in radar / geo coordinate
    # geo_box = None if atr is not geocoded.
    coord = ut.coordinate(metadata)
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), metadata)
    inps.pix_box = coord.check_box_within_data_coverage(inps.pix_box)
    inps.geo_box = coord.box_pixel2geo(inps.pix_box)
    # Out message
    inps.data_box = (0, 0, inps.width, inps.length)
    vprint('data   coverage in y/x: '+str(inps.data_box))
    vprint('subset coverage in y/x: '+str(inps.pix_box))
    vprint('data   coverage in lat/lon: '+str(coord.box_pixel2geo(inps.data_box)))
    vprint('subset coverage in lat/lon: '+str(inps.geo_box))
    vprint('------------------------------------------------------------------------')

    # DEM contour display
    if max(inps.pix_box[3] - inps.pix_box[1],
           inps.pix_box[2] - inps.pix_box[0]) > 2e3:
        inps.disp_dem_contour = False
        if inps.dem_file:
            vprint('area exceed 2000 pixels, turn off default DEM contour display')

    # Multilook, if too many subplots in one figure for less memory and faster speed
    if inps.multilook_num > 1:
        inps.multilook = True

    # Colormap
    inps.colormap = pp.auto_colormap_name(
        metadata,
        inps.colormap,
        datasetName=inps.dset[0],
        print_msg=inps.print_msg,
    )

    # Reference Point
    # Convert ref_lalo if existed, to ref_yx, and use ref_yx for the following
    # ref_yx is referenced to input data coverage, not subseted area for display
    if inps.ref_lalo:
        vprint(f'input reference point in lat/lon: {inps.ref_lalo}')
        if not inps.geo_box and not coord.lookup_file:
            print('WARNING: --ref-lalo is NOT supported when 1) file is radar-coded AND 2) no lookup table file found')
            print('    --> ignore the --ref-lalo input and continue.')
            inps.ref_lalo = []

        else:
            inps.ref_yx = coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1])[:2]
            vprint(f'input reference point in y  /x  : {inps.ref_yx}')

    # ref_lalo
    if inps.ref_yx and inps.geo_box:
        inps.ref_lalo = coord.radar2geo(inps.ref_yx[0], inps.ref_yx[1])[:2]
    elif 'REF_LAT' in metadata.keys():
        inps.ref_lalo = [float(metadata['REF_LAT']), float(metadata['REF_LON'])]
    else:
        inps.ref_lalo = None

    # Points of interest
    inps = pp.read_pts2inps(inps, coord)

    # Unit and Wrap
    inps.disp_unit, inps.wrap = pp.check_disp_unit_and_wrap(
        metadata,
        disp_unit=inps.disp_unit,
        wrap=inps.wrap,
        wrap_range=inps.wrap_range,
        print_msg=inps.print_msg,
    )

    # Map Projection via cartopy
    inps = check_map_projection(inps, metadata, print_msg=inps.print_msg)

    # Min / Max - Display
    if not inps.vlim:
        if (any(i in inps.key.lower() for i in ['coherence', '.cor'])
                or (inps.key == 'ifgramStack' and inps.dset[0].split('-')[0] in ['coherence'])
                or inps.dset[0] == 'cmask'):
            inps.vlim = [0.0, 1.0]
        elif inps.key in ['.int'] or inps.wrap:
            inps.vlim = inps.wrap_range

    # Transparency - Alpha
    if not inps.transparency:
        # Auto adjust transparency value when showing shaded relief DEM
        if inps.dem_file and inps.disp_dem_shade:
            inps.transparency = 0.8
        else:
            inps.transparency = 1.0

    # Flip Left-Right / Up-Down
    if inps.auto_flip:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(metadata, print_msg=inps.print_msg)

    # Figure Title
    if not inps.fig_title:
        inps.fig_title = pp.auto_figure_title(
            metadata['FILE_PATH'],
            datasetNames=inps.dset,
            inps_dict=vars(inps),
        )
    vprint('figure title: '+inps.fig_title)

    # Figure output file name
    if not inps.outfile:
        # ignore whitespaces in the filename
        fbase = inps.fig_title.replace(' ', '')
        if not inps.disp_whitespace:
            fbase += '_nws'
        inps.outfile = [f'{fbase}{inps.fig_ext}']

    inps = update_figure_setting(inps)
    return inps


def check_map_projection(inps, metadata, print_msg=True):
    """Check/initiate the map projection object via cartopy based on inps and metadata.

    Cartopy requires that:
      1. file is geocoded AND
      2. file coordinates are in the unit of degrees / meters AND
      3. set to display in geo-coordinates

    Use cartopy (by initiating inps.map_proj_obj) ONLY IF:
      1. show fancy lat/lon label via --lalo-label OR
      2. show coastline via --coastline

    This function will update the following variables:
        inps.map_proj_obj  # cartopy.crs.* object or None
        inps.coord_unit    # degree or meter
        inps.fig_coord     # geo or yx
    """

    inps.map_proj_obj = None
    inps.coord_unit = metadata.get('Y_UNIT', 'degrees').lower()
    if (inps.geo_box
            and inps.coord_unit.startswith(('deg', 'meter'))
            and inps.fig_coord == 'geo'
            and (inps.lalo_label or inps.coastline)):

        # get projection name from the data coord unit
        # https://scitools.org.uk/cartopy/docs/latest/crs/projections.html
        msg = 'initiate cartopy map projection: '
        if inps.coord_unit.startswith('deg'):
            inps.map_proj_obj = ccrs.PlateCarree()
            vprint(msg + 'PlateCarree')

        elif inps.coord_unit.startswith('meter'):
            if 'UTM_ZONE' in metadata.keys():
                utm_zone = metadata['UTM_ZONE']
                inps.map_proj_obj = ccrs.UTM(utm_zone)
                vprint(msg + f'UTM zone {utm_zone}')

                # check --lalo-label (works for PlateCarree only)
                if inps.lalo_label:
                    raise ValueError('--lalo-label is NOT supported for projection: UTM')

            else:
                print(f'WARNING: Un-recognized coordinate unit: {inps.coord_unit}')
                print('    Switch to the native Y/X and continue to plot')
                inps.fig_coord = 'yx'

    return inps


##################################################################################################
def update_data_with_plot_inps(data, metadata, inps):
    # 1. spatial referencing with respect to the seed point
    if inps.ref_yx:
        inps.ref_box = [inps.ref_yx[1],     inps.ref_yx[0],
                        inps.ref_yx[1] + 1, inps.ref_yx[0] + 1]

        # update ref_y/x to subset
        ref_y = inps.ref_yx[0] - inps.pix_box[1]
        ref_x = inps.ref_yx[1] - inps.pix_box[0]

        # update ref_y/x for multilooking
        if inps.multilook_num > 1:
            ref_y = int((ref_y - int(inps.multilook_num / 2)) / inps.multilook_num)
            ref_x = int((ref_x - int(inps.multilook_num / 2)) / inps.multilook_num)

        if len(data.shape) == 2:
            # read ref_val
            if 0 <= ref_y < data.shape[-2] and 0 <= ref_x < data.shape[-1]:
                ref_val = data[ref_y, ref_x]
            else:
                ref_val = readfile.read(
                    inps.file,
                    datasetName=inps.dset[0],
                    box=inps.ref_box,
                    print_msg=False,
                )[0]

            # applying spatial referencing
            if not np.ma.is_masked(ref_val) and not np.isnan(ref_val):
                data -= ref_val
                vprint(f'set reference pixel to: {inps.ref_yx}')
            else:
                msg = f'WARNING: input reference pixel ({ref_y}, {ref_x}) has either masked or NaN value!'
                msg += ' -> skip re-referencing.'
                print(msg)
                inps.ref_yx = None

        elif len(data.shape) == 3:
            # read ref_val
            if 0 <= ref_y < data.shape[-2] and 0 <= ref_x < data.shape[-1]:
                ref_val = np.squeeze(data[:, ref_y, ref_x])
            elif inps.key == 'timeseries':
                ref_val = readfile.read(
                    inps.file,
                    datasetName=inps.dset,
                    box=inps.ref_box,
                    print_msg=False,
                )[0]
            else:
                raise ValueError(f'input reference point {inps.ref_yx} is out of data coverage!')

            # apply spatial referencing
            if not np.ma.is_masked(ref_val) and np.all(~np.isnan(ref_val)):
                data -= np.tile(ref_val.reshape(-1, 1, 1), (1, data.shape[1], data.shape[2]))
                vprint(f'set reference pixel to: {inps.ref_yx}')
            else:
                msg = f'WARNING: input reference pixel ({ref_y}, {ref_x}) has either masked or NaN value!'
                msg += ' -> skip re-referencing.'
                print(msg)
                inps.ref_yx = None
    else:
        inps.ref_yx = None

    # 2. scale data based on the display unit and re-wrap
    (data,
     inps.disp_unit,
     inps.disp_scale,
     inps.wrap) = pp.scale_data4disp_unit_and_rewrap(
        data,
        metadata=metadata,
        disp_unit=inps.disp_unit,
        wrap=inps.wrap,
        wrap_range=inps.wrap_range,
        print_msg=inps.print_msg,
    )
    if inps.wrap:
        inps.vlim = inps.wrap_range

    # math operation
    if inps.math_operation:
        vprint(f'Apply math operation: {inps.math_operation}')
        if   inps.math_operation == 'square' :  data = np.square(data)
        elif inps.math_operation == 'sqrt'   :  data = np.sqrt(data)
        elif inps.math_operation == 'reverse':  data *= -1.
        elif inps.math_operation == 'inverse':  data = 1. / data
        elif inps.math_operation == 'rad2deg':  data *= 180. / np.pi
        elif inps.math_operation == 'deg2rad':  data *= np.pi / 180.
        else: raise ValueError(f'un-recognized math operation: {inps.math_operation}')

    # 4. update display min/max
    inps.dlim = [np.nanmin(data), np.nanmax(data)]
    if not inps.vlim: # and data.ndim < 3:
        (inps.cmap_lut,
         inps.vlim,
         inps.unique_values) = pp.auto_adjust_colormap_lut_and_disp_limit(data, print_msg=inps.print_msg)
    vprint(f'data    range: {inps.dlim} {inps.disp_unit}')
    vprint(f'display range: {inps.vlim} {inps.disp_unit}')

    return data, inps


##################################################################################################
def prep_slice(cmd, auto_fig=False):
    """Prepare data from command line as input, for easy call plot_slice() externally.

    Parameters: cmd  - string, command to be run in terminal
    Returns:    data - 2D np.ndarray, data to be plotted
                atr  - dict, metadata
                inps - namespace, input argument for plot setup
    Example:
        from cartopy import crs as ccrs
        from mintpy.view import prep_slice, plot_slice

        # initiate matplotlib figure/axes
        subplot_kw = dict(projection=ccrs.PlateCarree())
        fig, ax = plt.subplots(figsize=[4, 3], subplot_kw=subplot_kw)

        # compose view.py command
        cmd = 'view.py geo_velocity.h5 velocity --mask geo_maskTempCoh.h5 --dem srtm1.dem --dem-nocontour '
        cmd += f'--sub-lon -91.7 -91.4 --sub-lat -0.5 -0.3 -c jet -v -3 10 '
        cmd += '--cbar-loc bottom --cbar-nbins 3 --cbar-ext both --cbar-label "LOS velocity [cm/year]" '
        cmd += '--lalo-step 0.2 --lalo-loc 1 0 1 0 --scalebar 0.3 0.80 0.05 --notitle'

        # call prep/plot_slice()
        data, atr, inps = prep_slice(cmd)
        ax, inps, im, cbar = plot_slice(ax, data, atr, inps)
        plt.show()
    """
    # parse
    from mintpy.cli.view import cmd_line_parse
    iargs = cmd.split()[1:]

    # support option inputs of a list of characters (separated by whitespaces but quoted)
    # e.g.: --cbar-label "LOS velocity [cm/yar]" --title "S1 asc. velocity"
    # to be consistent with the behavior in command line parsing
    if any(x.startswith(('"', '\'')) for x in iargs):
        # backup and reset
        temp_iargs = list(iargs)
        iargs = []

        # get index of quoted list of characters
        ind0s = np.where([x.startswith(('"', '\'')) for x in temp_iargs])[0]
        ind1s = np.where([x.endswith(('"', '\'')) for x in temp_iargs])[0]
        for i, temp_iarg in enumerate(temp_iargs):
            if any(ind0 <= i <= ind1 for ind0, ind1 in zip(ind0s, ind1s)):
                # quoted list of characters
                for ind0, ind1 in zip(ind0s, ind1s):
                    if i == ind0:
                        temp = temp_iarg[1:]
                    elif ind0 < i < ind1:
                        temp += ' ' + temp_iarg
                    elif i == ind1:
                        temp += ' ' + temp_iarg[:-1]
                        iargs.append(temp)
                        break
            else:
                # regular unquoted list of characters
                iargs.append(temp_iarg)

    # run parse
    inps = cmd_line_parse(iargs)

    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    # read input args
    inps, atr = read_input_file_info(inps)
    inps = update_inps_with_file_metadata(inps, atr)

    inps.msk, inps.mask_file = pp.read_mask(
        inps.file,
        mask_file=inps.mask_file,
        datasetName=inps.dset[0],
        box=inps.pix_box,
        vmin=inps.mask_vmin,
        vmax=inps.mask_vmax,
        print_msg=inps.print_msg,
    )

    # read data
    data, atr = readfile.read(
        inps.file,
        datasetName=inps.dset[0],
        box=inps.pix_box,
        print_msg=inps.print_msg,
    )

    # reference in time
    if inps.ref_date:
        data -= readfile.read(
            inps.file,
            datasetName=inps.ref_date,
            box=inps.pix_box,
            print_msg=False,
        )[0]

    # reference in space for unwrapPhase
    if (inps.key in ['ifgramStack']
            and inps.dset[0].split('-')[0].startswith('unwrapPhase')
            and 'REF_Y' in atr.keys()):
        ref_y, ref_x = int(atr['REF_Y']), int(atr['REF_X'])
        ref_data = readfile.read(
            inps.file,
            datasetName=inps.dset[0],
            box=(ref_x, ref_y, ref_x+1, ref_y+1),
            print_msg=False,
        )[0]
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
    if auto_fig:
        figsize = [i/2.0 for i in inps.fig_size]
        subplot_kw = dict(projection=inps.map_proj_obj) if inps.map_proj_obj is not None else {}
        ax = plt.subplots(figsize=figsize, subplot_kw=subplot_kw)[1]
        return data, atr, inps, ax
    else:
        return data, atr, inps


##################################################################################################
def plot_slice(ax, data, metadata, inps):
    """Plot one slice of matrix
    Parameters: ax       : matplot.pyplot axes object
                data     : 2D np.ndarray,
                metadata : dictionary, attributes of data
                inps     : Namespace, optional, input options for display
    Returns:    ax       : matplot.pyplot axes object
                inps     : Namespace for input options
                im       : matplotlib.image.AxesImage object
                cbar     : matplotlib.colorbar.Colorbar object
    Example: See prep_slice() above for example usage.
    """
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    # colormap: str -> object
    if isinstance(inps.colormap, str):
        inps.colormap = pp.ColormapExt(
            inps.colormap,
            cmap_lut=inps.cmap_lut,
            vlist=inps.cmap_vlist,
        ).colormap

    # read DEM
    if inps.dem_file:
        dem, dem_metadata, dem_pix_box = pp.read_dem(
            inps.dem_file,
            pix_box=inps.pix_box,
            geo_box=inps.geo_box,
            print_msg=inps.print_msg,
        )

    vprint('display data in transparency: '+str(inps.transparency))
    num_row, num_col = data.shape
    lalo_digit = ut.get_lalo_digit4display(metadata, coord_unit=inps.coord_unit)

    #----------------------- Plot in Geo-coordinate --------------------------------------------#
    if (inps.geo_box
            and inps.coord_unit.startswith(('deg', 'meter'))
            and inps.fig_coord == 'geo'):
        vprint('plot in geo-coordinate')

        # extent info for matplotlib.imshow and other functions
        inps.extent = (inps.geo_box[0], inps.geo_box[2], inps.geo_box[3], inps.geo_box[1])  # (W, E, S, N)
        SNWE = (inps.geo_box[3], inps.geo_box[1], inps.geo_box[0], inps.geo_box[2])

        # Draw coastline using cartopy resolution parameters
        if inps.coastline:
            vprint(f'draw coast line with resolution: {inps.coastline}')
            ax.coastlines(resolution=inps.coastline, linewidth=inps.coastline_linewidth)

        # Plot DEM
        if inps.dem_file:
            vprint('plotting DEM background ...')
            pp.plot_dem_background(
                ax=ax,
                geo_box=inps.geo_box,
                dem=dem,
                inps=inps,
                print_msg=inps.print_msg,
            )

        # Reference (InSAR) data to a GNSS site
        coord = ut.coordinate(metadata)
        if inps.disp_gps and inps.gps_component and inps.ref_gps_site:
            ref_site_lalo = GPS(site=inps.ref_gps_site).get_stat_lat_lon(print_msg=False)
            y, x = coord.geo2radar(ref_site_lalo[0], ref_site_lalo[1])[0:2]
            ref_data = data[y - inps.pix_box[1], x - inps.pix_box[0]]
            data -= ref_data
            vprint('referencing InSAR data to the pixel nearest to GNSS station: '
                   f'{inps.ref_gps_site} at [{ref_site_lalo[0]:.6f}, {ref_site_lalo[1]:.6f}] '
                   f'by substrating {ref_data:.3f} {inps.disp_unit}')
            # do not show the original InSAR reference point
            inps.disp_ref_pixel = False

        # Plot data
        if inps.disp_dem_blend:
            im = pp.plot_blend_image(ax, data, dem, inps, print_msg=inps.print_msg)
        else:
            vprint('plotting data ...')
            im = ax.imshow(data, cmap=inps.colormap, vmin=inps.vlim[0], vmax=inps.vlim[1],
                           extent=inps.extent, origin='upper', interpolation=inps.interpolation,
                           alpha=inps.transparency, animated=inps.animation, zorder=1)

        # Draw faultline using GMT lonlat file
        if inps.faultline_file:
            pp.plot_faultline(
                ax=ax,
                faultline_file=inps.faultline_file,
                SNWE=SNWE,
                linewidth=inps.faultline_linewidth,
                min_dist=inps.faultline_min_dist,
                print_msg=inps.print_msg,
            )

        # Scale Bar
        if inps.coord_unit.startswith('deg') and (inps.geo_box[2] - inps.geo_box[0]) > 30:
            # do not plot scalebar if the longitude span > 30 deg
            inps.disp_scalebar = False

        if inps.disp_scalebar:
            vprint(f'plot scale bar: {inps.scalebar}')
            pp.draw_scalebar(
                ax=ax,
                geo_box=inps.geo_box,
                unit=inps.coord_unit,
                loc=inps.scalebar,
                labelpad=inps.scalebar_pad,
                font_size=inps.font_size,
                linewidth=inps.scalebar_linewidth,
            )

        # Lat Lon labels
        if inps.lalo_label:
            pp.draw_lalo_label(
                ax=ax,
                geo_box=inps.geo_box,
                lalo_step=inps.lalo_step,
                lalo_loc=inps.lalo_loc,
                lalo_max_num=inps.lalo_max_num,
                lalo_offset=inps.lalo_offset,
                font_size=inps.lalo_font_size if inps.lalo_font_size else inps.font_size,
                projection=inps.map_proj_obj,
                print_msg=inps.print_msg)
        else:
            ax.tick_params(which='both', direction='in', labelsize=inps.font_size,
                           left=True, right=True, top=True, bottom=True)

        # Plot Reference Point
        if inps.disp_ref_pixel and inps.ref_lalo:
            ax.plot(inps.ref_lalo[1], inps.ref_lalo[0],
                    inps.ref_marker, ms=inps.ref_marker_size)
            vprint('plot reference point')

        # Plot points of interest
        if inps.pts_lalo is not None:
            ax.plot(inps.pts_lalo[:, 1], inps.pts_lalo[:, 0],
                    inps.pts_marker, ms=inps.pts_marker_size,
                    mec='k', mew=1.)
            vprint('plot points of interest')

        # Show UNR GPS stations
        if inps.disp_gps:
            ax = pp.plot_gps(ax, SNWE, inps, metadata, print_msg=inps.print_msg)

        # Status bar
        if inps.dem_file:
            coord_dem = ut.coordinate(dem_metadata)
            dem_len, dem_wid = dem.shape

        def format_coord(x, y):
            # lat/lon
            msg = f'E={x:.{lalo_digit}f}, N={y:.{lalo_digit}f}'
            # value
            col = coord.lalo2yx(x, coord_type='lon') - inps.pix_box[0]
            row = coord.lalo2yx(y, coord_type='lat') - inps.pix_box[1]
            if 0 <= col < num_col and 0 <= row < num_row:
                v = data[row, col]
                msg += ', v=[]' if np.isnan(v) or np.ma.is_masked(v) else f', v={v:.3f}'
                # DEM
                if inps.dem_file:
                    dem_col = coord_dem.lalo2yx(x, coord_type='lon') - dem_pix_box[0]
                    dem_row = coord_dem.lalo2yx(y, coord_type='lat') - dem_pix_box[1]
                    if 0 <= dem_col < dem_wid and 0 <= dem_row < dem_len:
                        h = dem[dem_row, dem_col]
                        msg += ', h=[]' if np.isnan(h) else f', h={h:.1f}'
                # x/y
                msg += f', x={col+inps.pix_box[0]:.0f}'
                msg += f', y={row+inps.pix_box[1]:.0f}'
            return msg

        ax.format_coord = format_coord

    #------------------------ Plot in Y/X-coordinate ------------------------------------------------#
    else:
        inps.fig_coord = 'yx'
        vprint('plotting in Y/X coordinate ...')

        # Plot DEM
        if inps.dem_file:
            vprint('plotting DEM background ...')
            pp.plot_dem_background(
                ax=ax,
                geo_box=None,
                dem=dem,
                inps=inps,
                print_msg=inps.print_msg,
            )

        # extent = (left, right, bottom, top) in data coordinates
        inps.extent = (inps.pix_box[0]-0.5, inps.pix_box[2]-0.5,
                       inps.pix_box[3]-0.5, inps.pix_box[1]-0.5)

        # Plot Data
        if inps.disp_dem_blend:
            im = pp.plot_blend_image(ax, data, dem, inps, print_msg=inps.print_msg)
        else:
            vprint('plotting data ...')
            im = ax.imshow(data, cmap=inps.colormap, vmin=inps.vlim[0], vmax=inps.vlim[1],
                           extent=inps.extent, interpolation=inps.interpolation,
                           alpha=inps.transparency, zorder=1)
        ax.tick_params(labelsize=inps.font_size)

        # Plot Seed Point
        if inps.disp_ref_pixel:
            ref_y, ref_x = None, None
            if inps.ref_yx:
                ref_y, ref_x = inps.ref_yx[0], inps.ref_yx[1]
            elif 'REF_Y' in metadata.keys():
                ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])

            if ref_y and ref_x:
                ax.plot(ref_x, ref_y, inps.ref_marker, ms=inps.ref_marker_size)
                vprint('plot reference point')

        # Plot points of interest
        if inps.pts_yx is not None:
            ax.plot(inps.pts_yx[:, 1], inps.pts_yx[:, 0],
                    inps.pts_marker, ms=inps.ref_marker_size,
                    mec='black', mew=1.)
            vprint('plot points of interest')

        # temporary test code
        temp_test = False
        if temp_test:
            # Champlain Towers South AOI
            pts_yx = np.array([
                [929,1456],
                [933,1457],
                [933,1436],
                [930,1431],
                [929,1456],
            ])
            ax.plot(pts_yx[:, 1], pts_yx[:, 0], '-', ms=inps.ref_marker_size, mec='black', mew=1.)

        ax.set_xlim(inps.extent[0:2])
        ax.set_ylim(inps.extent[2:4])

        # Status bar

        # read lats/lons if exist
        geom_file = os.path.join(os.path.dirname(metadata['FILE_PATH']), 'inputs/geometryRadar.h5')
        if os.path.isfile(geom_file):
            geom_ds_list = readfile.get_dataset_list(geom_file)
            if all(x in geom_ds_list for x in ['latitude', 'longitude']):
                lats = readfile.read(geom_file, datasetName='latitude',  box=inps.pix_box, print_msg=False)[0]
                lons = readfile.read(geom_file, datasetName='longitude', box=inps.pix_box, print_msg=False)[0]
            else:
                msg = f'WARNING: no latitude / longitude found in file: {os.path.basename(geom_file)}, '
                msg += 'skip showing lat/lon in the status bar.'
                vprint(msg)
                geom_file = None
        else:
            geom_file = None

        def format_coord(x, y):
            # y/x
            msg = f'x={x:.1f}, y={y:.1f}'
            # value
            col = int(np.rint(x - inps.pix_box[0]))
            row = int(np.rint(y - inps.pix_box[1]))
            if 0 <= col < num_col and 0 <= row < num_row:
                v = data[row, col]
                msg += ', v=[]' if np.isnan(v) or np.ma.is_masked(v) else f', v={v:.3f}'
                # DEM
                if inps.dem_file:
                    h = dem[row, col]
                    msg += ', h=[]' if np.isnan(h) else f', h={h:.1f} m'
                # lat/lon
                if geom_file:
                    msg += f', E={lons[row, col]:.{lalo_digit}f}'
                    msg += f', N={lats[row, col]:.{lalo_digit}f}'
            return msg

        ax.format_coord = format_coord


    #---------------------- Figure Setting ----------------------------------------#

    # 3.1 Colorbar
    cbar = None
    if inps.disp_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(inps.cbar_loc, inps.cbar_size, pad=inps.cbar_size, axes_class=plt.Axes)
        inps, cbar = pp.plot_colorbar(inps, im, cax)

    # 3.2 Title
    if inps.disp_title:
        ax.set_title(inps.fig_title, fontsize=inps.font_size, color=inps.font_color)

    # 3.3 Flip Left-Right / Up-Down
    if inps.flip_lr:
        vprint('flip figure left and right')
        ax.invert_xaxis()

    if inps.flip_ud:
        vprint('flip figure up and down')
        ax.invert_yaxis()

    # 3.4 Turn off axis
    if not inps.disp_axis:
        ax.axis('off')
        vprint('turn off axis display')

    # 3.5 Tick labels
    if inps.disp_tick:
        # move x-axis tick label to the top if colorbar is at the bottom
        if inps.cbar_loc == 'bottom':
            ax.tick_params(bottom=False, top=True, labelbottom=False, labeltop=True)

        # manually turn ON to enable tick labels for UTM with cartopy
        # link: https://github.com/SciTools/cartopy/issues/491
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True)
    else:
        # turn off tick labels
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    # rotate Y-axis tick labels
    # link: https://stackoverflow.com/questions/10998621
    if inps.ylabel_rot:
        kwargs = dict(rotation=inps.ylabel_rot)
        # center the vertical alignment for vertical tick labels
        if inps.ylabel_rot % 90 == 0:
            kwargs['va'] = 'center'
        plt.setp(ax.get_yticklabels(), **kwargs)
        vprint(f'rotate Y-axis tick labels by {inps.ylabel_rot} deg')

    return ax, inps, im, cbar


##################################################################################################
def read_input_file_info(inps):
    # File Basic Info
    atr = readfile.read_attribute(inps.file)
    msg = 'input file is '
    if not inps.file.endswith(('.h5', '.he5')):
        msg += '{} '.format(atr['PROCESSOR'])
    msg += '{} file: {}'.format(atr['FILE_TYPE'], os.path.abspath(inps.file))
    if 'DATA_TYPE' in atr.keys():
        msg += ' in {} format'.format(atr['DATA_TYPE'])

    vprint(f'run {os.path.basename(__file__)} in {version.version_description}')
    vprint(msg)

    ## size and name
    inps.length = int(atr['LENGTH'])
    inps.width = int(atr['WIDTH'])
    inps.key = atr['FILE_TYPE']
    inps.fileBase = os.path.splitext(os.path.basename(inps.file))[0]
    inps.fileExt = os.path.splitext(inps.file)[1]
    vprint(f'file size in y/x: {(inps.length, inps.width)}')

    # File dataset List
    inps.sliceList = readfile.get_slice_list(inps.file, no_complex=True)

    # Read input list of dataset to display
    inps, atr = read_dataset_input(inps)

    return inps, atr


def search_dataset_input(all_list, in_list=[], in_num_list=[], search_dset=True):
    """Get dataset(es) from input dataset / dataset_num"""
    # make a copy to avoid weird variable behavior
    in_num_list = [x for x in in_num_list]

    # in_list --> in_num_list --> outNumList --> outList
    if in_list:
        if isinstance(in_list, str):
            in_list = [in_list]

        tempList = []
        if search_dset:
            for ds in in_list:
                # style of regular expression
                if '*' not in ds:
                    ds = f'*{ds}*'
                ds = ds.replace('*','.*')

                # search
                tempList += [e for e in all_list if re.match(ds, e) is not None]

        else:
            tempList += [i for i in in_list if i in all_list]
        tempList = sorted(list(set(tempList)))
        in_num_list += [all_list.index(e) for e in tempList]

    # in_num_list --> outNumList
    outNumList = sorted(list(set(in_num_list)))

    # outNumList --> outList
    outList = [all_list[i] for i in outNumList]

    return outList, outNumList


def read_dataset_input(inps):
    """Check input / exclude / reference dataset input with file dataset list"""
    # read inps.dset + inps.dsetNumList --> inps.dsetNumList
    if len(inps.dset) > 0 or len(inps.dsetNumList) > 0:
        # message
        if len(inps.dset) > 0:
            vprint(f'input dataset: "{inps.dset}"')

        # special rule for special file types
        if inps.key == 'velocity':
            inps.search_dset = False
            vprint(f'turning glob search OFF for {inps.key} file')

        elif inps.key == 'timeseries' and len(inps.dset) == 1 and '_' in inps.dset[0]:
            date1, date2 = inps.dset[0].split('_')
            inps.dset = [date2]
            inps.ref_date = date1

        # search
        inps.dsetNumList = search_dataset_input(
            all_list=inps.sliceList,
            in_list=inps.dset,
            in_num_list=inps.dsetNumList,
            search_dset=inps.search_dset)[1]

    else:
        # default dataset to display for certain type of files
        if inps.key == 'ifgramStack':
            inps.dset = ['unwrapPhase']
        elif inps.key == 'HDFEOS':
            inps.dset = ['displacement']
        elif inps.key == 'giantTimeseries':
            inps.dset = 'recons'
        elif inps.key == 'giantIfgramStack':
            obj = giantIfgramStack(inps.file)
            obj.open(print_msg=False)
            inps.dset = [obj.sliceList[0].split('-')[0]]
        else:
            inps.dset = inps.sliceList
            # do not plot 3D-bperp by default
            if inps.key == 'geometry':
                inps.dset = [x for x in inps.dset if not x.startswith('bperp')]

        inps.dsetNumList = search_dataset_input(
            all_list=inps.sliceList,
            in_list=inps.dset,
            in_num_list=inps.dsetNumList,
            search_dset=inps.search_dset)[1]

    # read inps.exDsetList
    inps.exDsetList, inps.exDsetNumList = search_dataset_input(
        all_list=inps.sliceList,
        in_list=inps.exDsetList,
        in_num_list=[],
        search_dset=inps.search_dset)

    # read inps.plot_drop_ifgram
    drop_num_list = []
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']

    if not inps.plot_drop_ifgram:
        if ftype == 'ifgramStack':
            vprint('do not show the dropped interferograms')
            date12_drop_list = ifgramStack(inps.file).get_drop_date12_list()
            drop_slice_list = [x for x in inps.sliceList if x.split('-')[1] in date12_drop_list]
            drop_num_list = [inps.sliceList.index(x) for x in drop_slice_list]
        else:
            print(f'WARNING: --show-kept option does not apply to file type: {ftype}, ignore and continue.')
            inps.plot_drop_ifgram = True

    # get inps.dset
    inps.dsetNumList = sorted(list(set(inps.dsetNumList) - set(inps.exDsetNumList) - set(drop_num_list)))
    inps.dset = [inps.sliceList[i] for i in inps.dsetNumList]
    inps.dsetNum = len(inps.dset)

    if inps.ref_date:
        if inps.key not in TIMESERIES_KEY_NAMES:
            inps.ref_date = None

        ref_date = search_dataset_input(
            all_list=inps.sliceList,
            in_list=[inps.ref_date],
            in_num_list=[],
            search_dset=inps.search_dset)[0][0]

        if not ref_date:
            msg = f'WARNING: input reference date {inps.ref_date} is not included in input file!'
            msg += 'Ignore it and continue'
            print(msg)
            inps.ref_date = None
        else:
            inps.ref_date = ref_date

    if inps.key in ['ifgramStack']:
        vprint(f'num of datasets in file {os.path.basename(inps.file)}: {len(inps.sliceList)}')
        vprint(f'num of datasets to exclude: {len(inps.exDsetList)}')
        vprint(f'num of datasets to display: {len(inps.dset)}')
    else:
        vprint(f'num of datasets in file {os.path.basename(inps.file)}: {len(inps.sliceList)}')
        vprint(f'datasets to exclude ({len(inps.exDsetList)}):\n{inps.exDsetList}')
        vprint(f'datasets to display ({len(inps.dset)}):\n{inps.dset}')
    if inps.ref_date and inps.key in TIMESERIES_KEY_NAMES:
        vprint(f'input reference date: {inps.ref_date}')

    if inps.dsetNum == 0:
        msg = 'No input dataset found!'
        msg += f'\navailable datasets:\n{inps.sliceList}'
        raise Exception(msg)

    atr = readfile.read_attribute(inps.file, datasetName=inps.dset[0].split('-')[0])

    return inps, atr


def update_figure_setting(inps):
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
            # update length/width based on lat/lon
            if inps.geo_box and inps.fig_coord == 'geo':
                length = abs(inps.geo_box[3] - inps.geo_box[1])
                width  = abs(inps.geo_box[2] - inps.geo_box[0])
            # auto figure size
            inps.fig_size = pp.auto_figure_size(
                ds_shape=(length, width),
                disp_cbar=inps.disp_cbar,
                print_msg=inps.print_msg)

    # Multiple Plots
    else:
        if not inps.fig_size:
            inps.fig_size = pp.default_figsize_multi
        vprint(f'figure size : [{inps.fig_size[0]:.2f}, {inps.fig_size[1]:.2f}]')

        # Figure number (<= 200 subplots per figure)
        if not inps.fig_num:
            inps.fig_num = 1
            while inps.dsetNum/float(inps.fig_num) > 160.0:
                inps.fig_num += 1

        # Row/Column number
        if (inps.fig_row_num == 1 and inps.fig_col_num == 1
                and all(i not in inps.argv for i in ['--nrows', '--ncols'])):

            # calculate row and col number based on input info
            data_shape = [length*1.1, width]
            fig_size4plot = [inps.fig_size[0]*0.95, inps.fig_size[1]]
            inps.fig_row_num, inps.fig_col_num = pp.auto_row_col_num(
                inps.dsetNum,
                data_shape,
                fig_size4plot,
                inps.fig_num)
        inps.fig_num = float(inps.dsetNum) / float(inps.fig_row_num * inps.fig_col_num)
        inps.fig_num = np.ceil(inps.fig_num).astype(int)
        vprint('dataset number: '+str(inps.dsetNum))
        vprint('row     number: '+str(inps.fig_row_num))
        vprint('column  number: '+str(inps.fig_col_num))
        vprint('figure  number: '+str(inps.fig_num))

        # adjust figure size for single row plot, to avoid extra whitespace
        if inps.fig_row_num == 1 and '--figsize' not in inps.argv:
            inps.fig_size = pp.auto_figure_size(
                ds_shape=(length, width*inps.fig_col_num),
                disp_cbar=inps.disp_cbar,
                print_msg=False)
            vprint(f'row number is 1, adjust figure size to {inps.fig_size}')

        if not inps.font_size:
            inps.font_size = 12
            if inps.fig_row_num * inps.fig_col_num > 50:
                inps.font_size = 8

        # Output File Name
        if inps.outfile:
            inps.outdir = os.path.dirname(inps.outfile[0])
            inps.outfile_base, inps.fig_ext = os.path.splitext(os.path.basename(inps.outfile[0]))
            inps.fig_ext = inps.fig_ext.lower()

        else:
            inps.outdir = os.path.dirname(inps.file)
            inps.outfile_base = os.path.splitext(os.path.basename(inps.file))[0]
            # suffix
            if (inps.pix_box[2]-inps.pix_box[0])*(inps.pix_box[3]-inps.pix_box[1]) < width*length:
                inps.outfile_base += '_sub'
            if inps.wrap:
                inps.outfile_base += '_wrap'
                if (inps.wrap_range[1] - inps.wrap_range[0]) != 2*np.pi:
                    inps.outfile_base += str(inps.wrap_range[1] - inps.wrap_range[0])
            if inps.ref_date:
                inps.outfile_base += '_ref'+inps.ref_date
            if inps.exDsetList:
                inps.outfile_base += '_ex'

        # output file name list
        if inps.fig_num == 1:
            inps.outfile = [f'{inps.outfile_base}{inps.fig_ext}']
        else:
            inps.outfile = [f'{inps.outfile_base}_{str(j)}{inps.fig_ext}'
                            for j in range(1, inps.fig_num+1)]
        inps.outfile = [os.path.join(inps.outdir, outfile) for outfile in inps.outfile]

    return inps


def read_data4figure(i_start, i_end, inps, metadata):
    """Read multiple datasets for one figure into 3D matrix based on i_start/end"""

    # initiate output matrix
    data = np.zeros((
        i_end - i_start,
        int((inps.pix_box[3] - inps.pix_box[1]) / inps.multilook_num),
        int((inps.pix_box[2] - inps.pix_box[0]) / inps.multilook_num),
    ), dtype=np.float32)

    # common args for readfile.read()
    kwargs = dict(
        box=inps.pix_box,
        xstep=inps.multilook_num,
        ystep=inps.multilook_num,
        print_msg=inps.print_msg,
    )

    # reference pixel info for unwrapPhase
    if inps.file_ref_yx:
        ref_y, ref_x = inps.file_ref_yx
        ref_kwargs = dict(box=(ref_x, ref_y, ref_x+1, ref_y+1), print_msg=False)

    # fast reading for single dataset type
    if (len(inps.dsetFamilyList) == 1
            and inps.key in ['timeseries', 'ifgramStack', 'geometry', 'HDFEOS', 'giantTimeseries']):

        vprint('reading data as a 3D matrix ...')
        dset_list = [inps.dset[i] for i in range(i_start, i_end)]
        if not metadata.get('DATA_TYPE', 'float32').startswith('complex'):
            data[:] = readfile.read(inps.file, datasetName=dset_list, **kwargs)[0]
        else:
            # calculate amplitude time series in dB, e.g. slcStack.h5
            vprint('input data is complex, calculate its amplitude and continue')
            data[:] = np.abs(readfile.read(inps.file, datasetName=dset_list, **kwargs)[0])

        # reference pixel info for unwrapPhase
        if inps.dsetFamilyList[0].startswith('unwrapPhase') and inps.file_ref_yx:
            ref_data = readfile.read(inps.file, datasetName=dset_list, **ref_kwargs)[0]
            for i in range(data.shape[0]):
                mask = data[i, :, :] != 0.
                data[i, mask] -= ref_data[i]

    # slow reading with one 2D matrix at a time
    else:
        vprint('reading data as a list of 2D matrices ...')
        kwargs['print_msg'] = False
        prog_bar = ptime.progressBar(maxValue=i_end-i_start, print_msg=inps.print_msg)
        for i in range(i_start, i_end):
            prog_bar.update(i - i_start + 1, suffix=inps.dset[i].split('/')[-1])

            # read 2D matrix
            d = readfile.read(inps.file, datasetName=inps.dset[i], **kwargs)[0]

            # reference pixel info for unwrapPhase
            if inps.dset[i].startswith('unwrapPhase') and inps.file_ref_yx:
                ref_d = readfile.read(inps.file, datasetName=inps.dset[i], **ref_kwargs)[0]
                d[d!=0] -= ref_d

            # save the matrix
            data[i - i_start, :, :] = d
        prog_bar.close()

    # ref_date for timeseries
    if inps.ref_date:
        vprint('consider input reference date: '+inps.ref_date)
        ref_data = readfile.read(inps.file, datasetName=inps.ref_date, **kwargs)[0]
        data -= ref_data

    # check if all subplots share the same data unit, they could have/be:
    # 1) the same type OR
    # 2) velocity or timeseries OR
    # 3) horizontal/vertical output from asc_desc2horz_vert.py
    # 4) data/model output from load_gbis.py OR
    # 5) binary files with multiple undefined datasets, as band1, band2, etc.
    if (len(inps.dsetFamilyList) == 1
            or inps.key in ['timeseries', 'inversion']
            or all(d in inps.dsetFamilyList for d in ['horizontal', 'vertical'])
            or inps.dsetFamilyList == ['data','model','residual']
            or inps.dsetFamilyList == [f'band{i+1}' for i in range(len(inps.dsetFamilyList))]
            or all(d.endswith('Amplitude') for d in inps.dsetFamilyList)
            or all(d.endswith('Phase') for d in inps.dsetFamilyList)):
        same_unit4all_subplots = True
    else:
        same_unit4all_subplots = False

    # adjust data due to spatial referencing and unit related scaling
    if same_unit4all_subplots:
        data, inps = update_data_with_plot_inps(data, metadata, inps)
    else:
        if any(x in inps.argv for x in ['-u', '--unit']):
            msg = 'WARNING: -u/--unit option is disabled for multi-subplots with different units!'
            msg += 'Ignore it and continue.'
            print(msg)
        inps.disp_unit = None

    # mask
    if inps.zero_mask:
        vprint('masking pixels with zero value')
        data = np.ma.masked_where(data == 0., data)
    if inps.msk is not None:
        vprint('masking data')
        msk = np.tile(inps.msk, (data.shape[0], 1, 1))
        data = np.ma.masked_where(msk == 0., data)

    # update display min/max
    inps.dlim = [np.nanmin(data), np.nanmax(data)]
    if (same_unit4all_subplots
            and all(arg not in inps.argv for arg in ['-v', '--vlim', '--wrap'])
            and not (inps.dsetFamilyList[0].startswith('unwrap') and not inps.file_ref_yx)
            and inps.dsetFamilyList[0] not in ['bperp']):

        inps.cmap_lut, inps.vlim, inps.unique_values = pp.auto_adjust_colormap_lut_and_disp_limit(
            data,
            num_multilook=10,
            print_msg=False,
        )

    return data


def plot_subplot4figure(i, inps, ax, data, metadata):
    """Plot one subplot for one 3D array
    1) Plot DEM, data and reference pixel
    2) axes setting: tick, ticklabel, title, axis etc.
    """
    # Plot DEM
    if inps.dem_file:
        pp.plot_dem_background(
            ax=ax,
            geo_box=None,
            dem_shade=inps.dem_shade,
            dem_contour=inps.dem_contour,
            dem_contour_seq=inps.dem_contour_seq,
            inps=inps,
            print_msg=inps.print_msg)

    # Plot Data
    vlim = inps.vlim if inps.vlim is not None else [np.nanmin(data), np.nanmax(data)]
    inps.extent = (inps.pix_box[0]-0.5, inps.pix_box[2]-0.5,
                   inps.pix_box[3]-0.5, inps.pix_box[1]-0.5)
    im = ax.imshow(data, cmap=inps.colormap, vmin=vlim[0], vmax=vlim[1],
                   interpolation=inps.interpolation, alpha=inps.transparency,
                   extent=inps.extent, zorder=1)

    # Plot Seed Point
    if inps.disp_ref_pixel:
        ref_y, ref_x = None, None
        if inps.ref_yx:
            ref_y, ref_x = inps.ref_yx[0], inps.ref_yx[1]
        elif 'REF_Y' in metadata.keys():
            ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])

        if ref_y and ref_x:
            ax.plot(ref_x, ref_y, inps.ref_marker, ms=inps.ref_marker_size)

    ax.set_xlim(inps.extent[0:2])
    ax.set_ylim(inps.extent[2:4])

    ## Subplot Setting
    # Tick and Label
    if not inps.disp_tick or inps.fig_row_num * inps.fig_col_num > 10:
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    # status bar
    def format_coord(x, y):
        return f'x={x:.1f}, y={y:.1f}, v ='
    ax.format_coord = format_coord

    # Title
    if inps.disp_title:
        # get title
        subplot_title = None
        if inps.key in TIMESERIES_KEY_NAMES or inps.dset[0].startswith('bperp'):
            # support "/" in the dataset names, e.g. HDF-EOS5 and py2-mintpy formats
            date_str = inps.dset[i].replace('/','-').split('-')[-1]
            try:
                subplot_title = dt.datetime.strptime(date_str, '%Y%m%d').isoformat()[0:10]
            except:
                subplot_title = date_str

        else:
            # dset info - name & index
            title_str = inps.dset[i]
            title_ind = inps.sliceList.index(title_str)

            # ignore dataset family info if there is only one type
            if len(inps.dsetFamilyList) == 1 and '-' in title_str:
                title_str = title_str.split('-')[1]
                # for ifgramStack, show index in the date12 list to facilitate the network modification
                if inps.atr['FILE_TYPE'] == 'ifgramStack':
                    title_ind = inps.date12List.index(title_str)

            # title style depending on the number of subplots
            num_subplot = inps.fig_row_num * inps.fig_col_num
            if num_subplot <= 6:
                subplot_title = title_str
            elif 6 < num_subplot <= 20:
                subplot_title = f'{title_ind}\n{title_str}'
            elif 20 < num_subplot <= 50:
                subplot_title = title_str.replace('_','\n').replace('-','\n')
            else:
                subplot_title = f'{title_ind}'

        # plot title
        if subplot_title:
            if inps.fig_title_in:
                prop = dict(size=inps.font_size)
                pp.add_inner_title(ax, subplot_title, loc=1, prop=prop)
            else:
                kwarg = dict(fontsize=inps.font_size)
                # mark dropped interferograms in bold crimson
                if inps.dset[i] in inps.dropDatasetList:
                    kwarg['color'] = 'crimson'
                    kwarg['fontweight'] = 'bold'
                else:
                    # special titles for Sentinel-1 data
                    if metadata.get('PLATFORM', None) == 'Sen' and inps.disp_title4sentinel1:
                        # display S1A/B in different colors
                        s1_sensor = metadata['SENTINEL1_SENSOR'].split()[i]
                        kwarg['color'] = 'C0' if s1_sensor == 'A' else 'C1'
                        # display IPF in subplot title
                        s1_IPF = metadata['SENTINEL1_IPF'].split()[i]
                        subplot_title += f' : {s1_IPF}'
                ax.set_title(subplot_title, **kwarg)

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
    """Plot one figure with multiple subplots
    1) create figure
    2) read all data into 3D array
    3) loop to plot each subplot using plot_subplot4figure()
    4) common colorbar and save
    """
    fig_title = f'Figure {str(j)} - {inps.outfile[j-1]}'
    vprint('----------------------------------------')
    vprint(fig_title)

    # Open a new figure object
    fig, axs = plt.subplots(
        num=j, figsize=inps.fig_size,
        nrows=inps.fig_row_num,
        ncols=inps.fig_col_num,
        sharex=True, sharey=True,
    )
    fig.canvas.manager.set_window_title(fig_title)
    axs = axs.flatten()

    # Read all data for the current figure into 3D np.array
    i_start = (j - 1) * inps.fig_row_num * inps.fig_col_num
    i_end = min([inps.dsetNum, i_start + inps.fig_row_num * inps.fig_col_num])
    data = read_data4figure(i_start, i_end, inps, metadata)

    if isinstance(inps.colormap, str):
        inps.colormap = pp.ColormapExt(
            inps.colormap, cmap_lut=inps.cmap_lut, vlist=inps.cmap_vlist,
        ).colormap

    # Loop - Subplots
    vprint('plotting ...')
    prog_bar = ptime.progressBar(maxValue=i_end-i_start, print_msg=inps.print_msg)
    for i in range(i_start, i_end):
        idx = i - i_start
        prog_bar.update(idx+1, suffix=inps.dset[i].split('/')[-1])

        # plot subplot
        im = plot_subplot4figure(
            i, inps,
            ax=axs[idx],
            data=data[idx, :, :],
            metadata=metadata)

        # add colorbar for each subplot
        if inps.disp_cbar and not inps.vlim:
            cbar = fig.colorbar(im, ax=axs[idx], pad=0.03, shrink=0.5, aspect=30, orientation='vertical')

            # display unit as colorbar label
            data_unit = readfile.read_attribute(inps.file, datasetName=inps.dset[i]).get('UNIT', None)
            if data_unit:
                cbar.set_label(data_unit)
    prog_bar.close()
    del data

    # delete empty axes
    for i in range(i_end-i_start, len(axs)):
        fig.delaxes(axs[i])

    # Min and Max for this figure
    inps.dlim_all = [np.nanmin([inps.dlim_all[0], inps.dlim[0]]),
                     np.nanmax([inps.dlim_all[1], inps.dlim[1]])]
    vprint(f'data    range: {inps.dlim} {inps.disp_unit}')
    if inps.vlim:
        vprint(f'display range: {inps.vlim} {inps.disp_unit}')

    # NOTE: For plt.subplots(), fig.tight_layout() should be run
    # before fig.add_axes(), which is the case of common colorbar
    # after fig.colorbar() and fig.set_size_inches(), which is the case of individual/multiple colorbars
    def adjust_subplots_layout(fig, inps):
        fig.subplots_adjust(
            left=0.02, right=0.98,
            bottom=0.02, top=0.98,
            wspace=0.05, hspace=0.05,
        )
        if inps.fig_wid_space or inps.fig_hei_space:
            fig.subplots_adjust(hspace=inps.fig_hei_space,
                                wspace=inps.fig_wid_space)
        elif inps.fig_tight_layout:
            fig.tight_layout()
        return

    # Colorbar
    if inps.disp_cbar:
        if not inps.vlim:
            vprint('Note: different color scale for EACH subplot!')
            vprint('Adjust figsize for the colorbar of each subplot.')
            fig.set_size_inches(inps.fig_size[0] * 1.1, inps.fig_size[1])
            adjust_subplots_layout(fig, inps)

        else:
            adjust_subplots_layout(fig, inps)
            # plot common colorbar
            cbar_length = 0.4
            if inps.fig_size[1] > 8.0:
                cbar_length /= 2
            vprint('show colorbar')
            fig.subplots_adjust(right=0.93)
            cax = fig.add_axes([0.94, (1.0-cbar_length)/2, 0.005, cbar_length])
            inps, cbar = pp.plot_colorbar(inps, im, cax)

    else:
        adjust_subplots_layout(fig, inps)

    # Save Figure
    if inps.save_fig:
        vprint(f'save figure to {os.path.abspath(inps.outfile[j-1])} with dpi={inps.fig_dpi}')
        fig.savefig(inps.outfile[j-1], bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
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
    inps.dsetFamilyList = sorted(list({x.split('-')[0] for x in inps.dset}))
    inps.dsetFamilyList = sorted(list({x.replace('Std','') for x in inps.dsetFamilyList}))
    if len(inps.dsetFamilyList) == 1 and inps.atr['FILE_TYPE'] == 'ifgramStack':
        inps.date12List = sorted(list({x.split('-')[1] for x in inps.sliceList}))

    if inps.multilook_num > 1 and inps.print_msg:
        print('multilook {0} by {0} with nearest interpolation'.format(inps.multilook_num))

    elif inps.multilook and inps.multilook_num == 1:
        ## calculate multilook_num
        # ONLY IF:
        #   inps.multilook is True (no --nomultilook input) AND
        #   inps.multilook_num ==1 (no --multilook-num input)
        # inps.multilook is used for this check ONLY
        inps.multilook_num = pp.auto_multilook_num(
            inps.pix_box, inps.fig_row_num * inps.fig_col_num,
            max_memory=inps.maxMemory,
            print_msg=inps.print_msg,
        )

    # multilook mask
    if inps.msk is not None and inps.multilook_num > 1:
        inps.msk = multilook_data(
            inps.msk,
            inps.multilook_num,
            inps.multilook_num,
            method='nearest',
        )

    # Reference pixel for timeseries and ifgramStack
    inps.file_ref_yx = None
    if inps.key in ['ifgramStack'] and 'REF_Y' in metadata.keys():
        ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])
        length, width = int(metadata['LENGTH']), int(metadata['WIDTH'])
        if 0 <= ref_y < length and 0 <= ref_x < width:
            inps.file_ref_yx = [ref_y, ref_x]
            vprint(f'consider reference pixel in y/x: {inps.file_ref_yx}')

    if inps.dsetNum > 10:
        inps.ref_marker_size /= 10.
    elif inps.dsetNum > 100:
        inps.ref_marker_size /= 20.

    # Check dropped interferograms
    inps.dropDatasetList = []
    if inps.key == 'ifgramStack' and inps.disp_title:
        obj = ifgramStack(inps.file)
        obj.open(print_msg=False)
        dropDate12List = obj.get_drop_date12_list()
        for i in inps.dsetFamilyList:
            inps.dropDatasetList += [f'{i}-{j}' for j in dropDate12List]
        vprint("mark interferograms with 'dropIfgram=False' in red colored title")

    # Read DEM
    if inps.dem_file:
        dem_metadata = readfile.read_attribute(inps.dem_file)
        if all(dem_metadata[i] == metadata[i] for i in ['LENGTH', 'WIDTH']):
            vprint(f'reading DEM: {os.path.basename(inps.dem_file)} ... ')
            dem = readfile.read(
                inps.dem_file,
                datasetName='height',
                box=inps.pix_box,
                xstep=inps.multilook_num,
                ystep=inps.multilook_num,
                print_msg=False,
            )[0]

            inps.dem_shade, inps.dem_contour, inps.dem_contour_seq = pp.prep_dem_background(
                dem=dem,
                inps=inps,
                print_msg=inps.print_msg,
            )

        else:
            inps.dem_file = None
            inps.transparency = 1.0
            msg = 'WARNING: DEM file has a different size from the data file. '
            msg += 'This feature is only supported for single subplot, and not for multi-subplots.'
            msg += '\n    --> Ignore it and continue.'
            print(msg)
    return inps


##################################################################################################
class viewer():
    """viewer class definition.

    Example:
        cmd = 'view.py timeseries.h5'
        inps = cmd_line_parse(cmd.split()[1:])
        from mintpy.view import viewer
        obj = viewer()
        obj.configure(inps)
        obj.plot()
    """

    def __init__(self, cmd=None, iargs=None):
        if cmd:
            iargs = cmd.split()[1:]
        self.cmd = cmd
        self.iargs = iargs

    def configure(self, inps):
        global vprint
        vprint = print if inps.print_msg else lambda *args, **kwargs: None

        # matplotlib backend setting
        if not inps.disp_fig:
            plt.switch_backend('Agg')

        inps, self.atr = read_input_file_info(inps)
        inps = update_inps_with_file_metadata(inps, self.atr)

        # --update option
        self.flag = 'run'
        if inps.update_mode and run_or_skip(inps) == 'skip':
            self.flag = 'skip'

        # copy inps to self object
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

        # read mask
        self.msk, self.mask_file = pp.read_mask(
            self.file,
            mask_file=self.mask_file,
            datasetName=self.dset[0],
            box=self.pix_box,
            vmin=self.mask_vmin,
            vmax=self.mask_vmax,
            print_msg=self.print_msg)

        return self.flag


    def plot(self):
        # One Subplot
        if self.dsetNum == 1:
            vprint('reading data ...')
            # read data
            data, self.atr = readfile.read(
                self.file,
                datasetName=self.dset[0],
                box=self.pix_box,
                print_msg=False)

            # reference in time
            if self.ref_date:
                data -= readfile.read(
                    self.file,
                    datasetName=self.ref_date,
                    box=self.pix_box,
                    print_msg=False)[0]

            # reference in space for unwrapPhase
            if (self.key in ['ifgramStack']
                    and self.dset[0].split('-')[0].startswith('unwrapPhase')
                    and 'REF_Y' in self.atr.keys()):
                ref_y, ref_x = int(self.atr['REF_Y']), int(self.atr['REF_X'])
                ref_data = readfile.read(
                    self.file,
                    datasetName=self.dset[0],
                    box=(ref_x, ref_y, ref_x+1, ref_y+1),
                    print_msg=False)[0]
                data[data != 0.] -= ref_data

            # masking - input options
            if self.zero_mask:
                vprint('masking pixels with zero value')
                data = np.ma.masked_where(data == 0., data)
            if self.msk is not None:
                vprint('masking data')
                data = np.ma.masked_where(self.msk == 0., data)
            else:
                self.msk = np.ones(data.shape, dtype=np.bool_)

            # masking - NO_DATA_VALUE
            no_data_val = readfile.get_no_data_value(self.file)
            if self.no_data_value is not None:
                vprint(f'masking pixels with NO_DATA_VALUE of {self.no_data_value}')
                data = np.ma.masked_where(data == self.no_data_value, data)
            elif no_data_val is not None and not np.isnan(no_data_val):
                vprint(f'masking pixels with NO_DATA_VALUE of {no_data_val}')
                data = np.ma.masked_where(data == no_data_val, data)

            # update/save mask info
            if np.ma.is_masked(data):
                self.msk *= ~data.mask
                self.msk *= ~np.isnan(data.data)
            else:
                self.msk *= ~np.isnan(data)

            # update data
            data, self = update_data_with_plot_inps(data, self.atr, self)

            # prepare figure
            subplot_kw = dict(projection=self.map_proj_obj) if self.map_proj_obj is not None else {}
            fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw=subplot_kw)
            if not self.disp_whitespace:
                fig.subplots_adjust(left=0,right=1,bottom=0,top=1)

            # plot
            self = plot_slice(ax, data, self.atr, self)[1]

            # save figure
            if self.save_fig:
                vprint(f'save figure to {os.path.abspath(self.outfile[0])} with dpi={self.fig_dpi}')
                if not self.disp_whitespace:
                    fig.savefig(self.outfile[0], transparent=True, dpi=self.fig_dpi, pad_inches=0.0)
                else:
                    fig.savefig(self.outfile[0], transparent=True, dpi=self.fig_dpi, bbox_inches='tight')
                if not self.disp_fig:
                    fig.clf()

        # Multiple Subplots
        else:
            # warn single-subplot options
            opt_names = ['--show-gps', '--coastline', '--lalo-label', '--lalo-step', '--scalebar',
                         '--pts-yx', '--pts-lalo', '--pts-file']
            opt_names = list(set(opt_names) & set(self.argv))
            for opt_name in opt_names:
                print(f'WARNING: {opt_name} is NOT supported for multi-subplots, ignore it and continue.')

            # prepare
            self = prepare4multi_subplots(self, metadata=self.atr)

            # plot
            self.dlim_all = [0., 0.]
            for j in range(1, self.fig_num + 1):
                plot_figure(j, self, metadata=self.atr)

            # stat
            if self.fig_num > 1:
                vprint('----------------------------------------')
                vprint(f'all data range: {self.dlim_all} {self.disp_unit}')
                if self.vlim:
                    vprint(f'display  range: {self.vlim} {self.disp_unit}')

        # Display Figure
        if self.disp_fig:
            vprint('showing ...')
            plt.show()
        return
