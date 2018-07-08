#!/usr/bin/env python3
#######################################################################
# Program is part of PySAR                                            #
# Copyright(c) 2013-2018, Zhang Yunjun, Joshua Zahner, Heresh Fattahi #
# Author:  Zhang Yunjun, Joshua Zahner, Heresh Fattahi                #
#######################################################################


import os
import sys
import argparse
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.widgets import Slider
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pysar.objects import timeseries, giantTimeseries, HDFEOS
from pysar.utils import readfile, ptime, plot as pp, utils as ut
from pysar.multilook import multilook_data


###########################################################################################
EXAMPLE = """example:
  tsview.py timeseries.h5
  tsview.py timeseries.h5  --wrap
  tsview.py timeseries.h5  --yx 300 400 --zero-first  --nodisplay
  tsview.py geo_timeseries.h5  --lalo 33.250 131.665  --nodisplay
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Interactive Time-series Viewer',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('timeseries_file', help='time series file to display')
    parser.add_argument('-n', dest='init_idx', metavar='NUM', type=int,
                        help='Epoch/slice number to display.')
    parser.add_argument('-m', '--mask', dest='mask_file',
                        help='mask to use. Default: geo_maskTempCoh.h5 for geocoded file and maskTempCoh.h5 for radar file')
    parser.add_argument('--error', dest='error_file',
                        help='txt file with error for each date.')
    parser.add_argument('--dem', dest='dem_file',
                        help='DEM file for background shaed relief')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                       help='initial pixel to plot in lat/lon coord')

    ref = parser.add_argument_group('Reference Pixel')
    ref.add_argument('--ref-yx', dest='ref_yx', type=int, metavar=('Y', 'X'), nargs=2,
                     help='change reference pixel to input location')
    ref.add_argument('--ref-lalo', dest='ref_lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                     help='change reference pixel to input location')
    ref.add_argument('--ref-color', dest='seed_color', metavar='COLOR', default='k',
                     help='marker color of reference point')
    ref.add_argument('--ref-symbol', dest='seed_symbol', metavar='SYMBOL', default='s',
                     help='marker symbol of reference point')
    ref.add_argument('--ref-size', dest='seed_size', metavar='SIZE_NUM', type=int, default=6,
                     help='marker size of reference point, default: 10')

    output = parser.add_argument_group('Output Setting')
    output.add_argument('-o', '--output', dest='fig_base',
                        help='Figure base name for output files')
    output.add_argument('--save', action='store_true', dest='save_fig',
                        help='save data and plot to files')
    output.add_argument('--nodisplay', action='store_false', dest='disp_fig',
                        help='save data and plot to files and do not display figures\n')
    output.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=150,
                        help='DPI - dot per inch - for display/write')

    disp = parser.add_argument_group('Display Setting')
    disp.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'),
                      type=float, nargs=2, default=[8.0, 4.5],
                      help='Figure size in inches - width and length. Default: 10.0 5.0\n' +
                           'i.e. 3.5 2 for ppt; ')
    disp.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                      help='Y limits for point plotting.')
    disp.add_argument('--vlim', dest='vlim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                      help='Display limits for matrix plotting.')
    disp.add_argument('--ref-date', dest='ref_date',
                      help='Change reference date for display')
    disp.add_argument('--exclude', '--ex', dest='ex_date_list',
                      nargs='*', help='Exclude date shown as gray.')
    disp.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                      help='Set displacement at first acquisition to zero.')
    disp.add_argument('--wrap', action='store_true',
                      help='re-wrap data to display data in fringes, for map display only.')
    disp.add_argument('-u', dest='disp_unit', metavar='UNIT', default='cm',
                      help='unit for display. Default: cm')
    disp.add_argument('-c', '--colormap', dest='colormap', default='jet',
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r etc.\n'
                           'Support colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html')
    disp.add_argument('-s', '--fontsize', dest='font_size',
                      type=int, default=10, help='Font size for display')
    disp.add_argument('--notitle', dest='disp_title',
                      action='store_false', help='Do not display title in TS plot.')
    disp.add_argument('--flip-lr', dest='flip_lr',
                      action='store_true', help='flip left-right')
    disp.add_argument('--flip-ud', dest='flip_ud',
                      action='store_true', help='flip up-down')
    disp.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=8.0,
                      help='Point marker size. Default: 12.0')
    # disp.add_argument('--mc','--markercolor', dest='marker_color', default='crimson',\
    #                  help='Point marker color. Default: crimson')
    disp.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0,
                      help='Edge width. Default: 1.0')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if (not inps.disp_fig or inps.fig_base) and not inps.save_fig:
        inps.save_fig = True
    if inps.ylim:
        inps.ylim = sorted(inps.ylim)
    return inps


###########################################################################################
def read_exclude_date(input_ex_date, dateListAll):
    ex_date_list = []
    ex_dates = []
    ex_flag = np.ones((len(dateListAll)), np.bool_)

    if input_ex_date:
        input_ex_date = list(input_ex_date)
        if input_ex_date:
            for ex_date in input_ex_date:
                if os.path.isfile(ex_date):
                    ex_date = ptime.read_date_list(ex_date)
                else:
                    ex_date = [ptime.yyyymmdd(ex_date)]
                ex_date_list += list(set(ex_date) - set(ex_date_list))

            # delete dates not existed in input file
            ex_date_list = sorted(list(set(ex_date_list).intersection(dateListAll)))
            ex_dates = ptime.date_list2vector(ex_date_list)[0]
            for i in ex_date_list:
                ex_flag[dateListAll.index(i)] = False
            print('exclude date:'+str(ex_date_list))
    return ex_date_list, ex_dates, ex_flag


def read_init_info(inps):
    # Time Series Info
    atr = readfile.read_attribute(inps.timeseries_file)
    inps.key = atr['FILE_TYPE']
    if inps.key == 'timeseries':
        obj = timeseries(inps.timeseries_file)
    elif inps.key == 'giantTimeseries':
        obj = giantTimeseries(inps.timeseries_file)
    elif inps.key == 'HDFEOS':
        obj = HDFEOS(inps.timeseries_file)
    else:
        raise ValueError('input file is {}, not timeseries.'.format(inps.key))
    obj.open()

    # default mask file
    if not inps.mask_file and 'masked' not in inps.timeseries_file:
        dir_name = os.path.dirname(inps.timeseries_file)
        if 'Y_FIRST' in atr.keys():
            inps.mask_file = os.path.join(dir_name, 'geo_maskTempCoh.h5')
        else:
            inps.mask_file = os.path.join(dir_name, 'maskTempCoh.h5')
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None

    # date info
    inps.date_list = obj.dateList
    inps.num_date = len(inps.date_list)
    inps.dates, inps.yearList = ptime.date_list2vector(inps.date_list)
    (inps.ex_date_list,
     inps.ex_dates,
     inps.ex_flag) = read_exclude_date(inps.ex_date_list, inps.date_list)

    # initial display index
    inps.ref_idx = obj.refIndex
    if inps.ref_date:
        inps.ref_idx = inps.date_list.index(inps.ref_date)
    if not inps.init_idx:
        if inps.ref_idx < inps.num_date / 2.:
            inps.init_idx = -2
        else:
            inps.init_idx = 2

    # Read Error List
    inps.error_ts = None
    inps.ex_error_ts = None
    if inps.error_file:
        error_fileContent = np.loadtxt(inps.error_file, dtype=bytes).astype(str)
        inps.error_ts = error_fileContent[:, 1].astype(np.float)*inps.unit_fac
        if inps.ex_date_list:
            e_ts = inps.error_ts[:]
            inps.ex_error_ts = e_ts[inps.ex_flag == 0]
            inps.error_ts = e_ts[inps.ex_flag == 1]

    # Zero displacement for 1st acquisition
    if inps.zero_first:
        inps.zero_idx = min(0, np.min(np.where(inps.ex_flag)[0]))

    inps.length, inps.width = obj.length, obj.width
    print('data size in (y0, y1, x0, x1): {}'.format((0, inps.length, 0, inps.width)))
    try:
        inps.lon0 = float(atr['X_FIRST'])
        inps.lat0 = float(atr['Y_FIRST'])
        inps.lon_step = float(atr['X_STEP'])
        inps.lat_step = float(atr['Y_STEP'])
        inps.lon1 = inps.lon0 + inps.width *inps.lon_step
        inps.lat1 = inps.lat0 + inps.length*inps.lat_step
        inps.geocoded = True
        print(('data size in (N, S, W, E): '
               '({%.4f}, {%.4f}, {%.4f}, {%.4f})'.format(inps.lat1,
                                                         inps.lat0,
                                                         inps.lon0,
                                                         inps.lon1)))
    except:
        inps.geocoded = False
        lookup_file = ut.get_lookup_file('./INPUTS/geometryRadar.h5')
        if lookup_file is not None:
            inps.lats = readfile.read(lookup_file, datasetName='latitude')[0]
            inps.lons = readfile.read(lookup_file, datasetName='longitude')[0]
        else:
            inps.lats = None
            inps.lons = None

    # reference pixel
    if inps.ref_lalo and 'Y_FIRST' in atr.keys():
        y = int((inps.ref_lalo[0] - inps.lat0) / inps.lat_step + 0.5)
        x = int((inps.ref_lalo[1] - inps.lon0) / inps.lon_step + 0.5)
        inps.ref_yx = [y, x]
    if not inps.ref_yx:
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]

    # Initial Pixel Coord
    if inps.lalo:
        if 'Y_FIRST' in atr.keys():
            y = int((inps.lalo[0] - inps.lat0) / inps.lat_step + 0.5)
            x = int((inps.lalo[1] - inps.lon0) / inps.lon_step + 0.5)
            inps.yx = [y, x]
        elif lookup_file is not None:
            y, x = ut.glob2radar(inps.lalo[0], inps.lalo[1],
                                 lookup_file, atr, print_msg=False)[0:2]
            inps.yx = [y, x]
    if not inps.yx:
        inps.yx = inps.ref_yx

    # Flip up-down / left-right
    if not inps.flip_lr and not inps.flip_ud:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(atr)

    inps.disp_unit_v = inps.disp_unit
    if inps.wrap:
        inps.range2phase = -4. * np.pi / float(atr['WAVELENGTH'])
        if   'cm' in inps.disp_unit:   inps.range2phase /= 100.
        elif 'mm' in inps.disp_unit:   inps.range2phase /= 1000.
        inps.disp_unit_v = 'radian'
        inps.vlim = (-np.pi, np.pi)

    return inps, atr


def read_timeseries_data(inps):
    ts_data, atr = readfile.read(inps.timeseries_file)
    ref_phase = ts_data[:, inps.ref_yx[0], inps.ref_yx[1]]
    ts_data -= np.tile(ref_phase.reshape(-1, 1, 1), (1, inps.length, inps.width))
    ts_data -= np.tile(ts_data[inps.ref_idx, :, :], (inps.num_date, 1, 1))

    # Display Unit
    (ts_data,
     inps.disp_unit,
     inps.unit_fac) = pp.scale_data2disp_unit(ts_data,
                                              metadata=atr,
                                              disp_unit=inps.disp_unit)

    # Mask file: input mask file + non-zero ts pixels
    mask = np.ones((inps.length, inps.width), np.bool_)
    msk = pp.read_mask(inps.timeseries_file,
                       mask_file=inps.mask_file,
                       datasetName='displacement')[0]
    mask[msk == 0.] = False
    del msk

    ts_stack = np.sum(ts_data, axis=0)
    mask[ts_stack == 0.] = False
    mask[np.isnan(ts_stack)] = False
    del ts_stack

    print('masking data')
    ts_mask = np.tile(mask, (inps.num_date, 1, 1))
    ts_data[ts_mask == 0] = np.nan
    ts_data[:, inps.ref_yx[0], inps.ref_yx[1]] = 0.   # keep value on reference pixel
    del ts_mask

    # default vlim
    data_lim = (np.nanmin(ts_data), np.nanmax(ts_data))
    ts_data_mli = multilook_data(ts_data, 10, 10)
    if not inps.vlim:
        inps.vlim = (np.nanmin(ts_data_mli), np.nanmax(ts_data_mli))
    print('data    range: {} {}'.format(data_lim, inps.disp_unit))
    print('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    # default ylim
    if not inps.ylim:
        if inps.zero_first:
            ts_data_mli -= np.tile(ts_data_mli[inps.zero_idx, :, :], (inps.num_date, 1, 1))
        ymin, ymax = np.nanmin(ts_data_mli), np.nanmax(ts_data_mli)
        ybuffer = (ymax - ymin) * 0.05
        inps.ylim = (ymin - ybuffer, ymax + ybuffer)
    del ts_data_mli

    return ts_data, mask, inps


def plot_init_map(ax, ts_data, inps):
    if inps.dem_file:
        print('plotting DEM background ...')
        dem = readfile.read(inps.dem_file, datasetName='height')[0]
        ax = pp.plot_dem_background(ax=ax,
                                    geo_box=None,
                                    dem=dem,
                                    inps_dict=vars(inps))
        del dem

    idx = inps.init_idx
    d_v = np.array(ts_data[idx, :, :])
    if inps.wrap:
        d_v *= inps.range2phase
        d_v -= np.round(d_v/(2*np.pi)) * (2*np.pi)
    im = ax.imshow(d_v,
                   cmap=inps.colormap,
                   clim=inps.vlim,
                   interpolation='nearest')

    # Reference Pixel
    ax.plot(inps.ref_yx[1],
            inps.ref_yx[0],
            inps.seed_color+inps.seed_symbol,
            ms=inps.seed_size)

    # Initial Pixel
    if inps.yx != inps.ref_yx:
        ax.plot(inps.yx[1], inps.yx[0], 'ro', markeredgecolor='black')

    ax.set_xlim(-0.5, inps.width-0.5)
    ax.set_ylim(inps.length-0.5, -0.5)

    # Status Bar
    def format_coord(x, y):
        col = int(x+0.5)
        row = int(y+0.5)
        msg = 'x={:.1f}, y={:.1f}'.format(x, y)
        if 0 <= col < inps.width and 0 <= row < inps.length:
            try:
                z = d_v[row, col]
                msg += ', value={:.4f}'.format(z)
            except:
                msg += ', value=[]'
        try:
            lon = inps.lon0 + x*inps.lon_step
            lat = inps.lat0 + y*inps.lat_step
            msg += ', lon={:.4f}, lat={:.4f}'.format(lon, lat)
        except:
            pass
        return msg
    ax.format_coord = format_coord

    # Title and Axis Label
    disp_date = inps.dates[idx].strftime('%Y-%m-%d')
    ax.set_title('N = {}, Time = {}'.format(idx, disp_date))

    # Flip up-down / left-right
    if inps.flip_lr:
        ax.invert_xaxis()
        print('flip map left and right')
    if inps.flip_ud:
        ax.invert_yaxis()
        print('flip map up and down')

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "2%", pad="2%")
    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=inps.font_size)
    if inps.wrap:
        cbar.set_ticks([-np.pi, 0, np.pi])
        cbar.ax.set_yticklabels([r'-$\pi$', '0', r'$\pi$'])
    cbar.set_label('Displacement ({})'.format(inps.disp_unit_v))
    return ax, im


def plot_init_time_slider(ax, year_list, init_idx=-1, ref_idx=0):
    val_step = np.min(np.diff(year_list))
    tslider = Slider(ax, label='Years',
                     valmin=year_list[0]-val_step,
                     valmax=year_list[-1]+val_step,
                     valinit=year_list[init_idx],
                     valstep=val_step)

    tslider.ax.bar(year_list, np.ones(len(year_list)), facecolor='black', width=0.01, ecolor=None)
    tslider.ax.bar(year_list[init_idx], 1., facecolor='black', width=0.01, ecolor=None)
    tslider.ax.bar(year_list[ref_idx], 1., facecolor='crimson', width=0.03, ecolor=None)

    # xaxis tick format
    if np.floor(year_list[-1]) == np.floor(year_list[0]):
        digit = 10.
    else:
        digit = 1.
    tslider.ax.set_xticks(np.round(np.linspace(year_list[0], year_list[-1], num=5) * digit) / digit)
    tslider.ax.xaxis.set_minor_locator(MultipleLocator(1./12.))
    tslider.ax.set_yticks([])
    tslider.ax.set_facecolor('lightgoldenrodyellow')
    return tslider


def plot_timeseries_errorbar(ax, dis_ts, inps):
    dates = np.array(inps.dates)
    d_ts = dis_ts[:]
    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        (_, caps, _) = ax.errorbar(inps.ex_dates, ex_d_ts, yerr=inps.ex_error_ts,
                                   fmt='-o', color='gray', ms=inps.marker_size,
                                   lw=0, alpha=1, mfc='gray', elinewidth=inps.edge_width,
                                   ecolor='black', capsize=inps.marker_size*0.5)
        for cap in caps:
            cap.set_markeredgewidth(inps.edge_width)

    # Plot kept dates
    (_, caps, _) = ax.errorbar(dates, d_ts, yerr=inps.error_ts,
                               fmt='-o', ms=inps.marker_size,
                               lw=0, alpha=1, elinewidth=inps.edge_width,
                               ecolor='black', capsize=inps.marker_size*0.5)
    for cap in caps:
        cap.set_markeredgewidth(inps.edge_width)
    return ax


def plot_timeseries_scatter(ax, dis_ts, inps):
    dates = np.array(inps.dates)
    d_ts = dis_ts[:]
    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        ax.scatter(inps.ex_dates, ex_d_ts, s=inps.marker_size**2, color='gray')

    # Plot kept dates
    ax.scatter(dates, d_ts, s=inps.marker_size**2)
    return ax


def plot_point_timeseries(yx, fig, ax, ts_data, inps):
    """Plot point time series displacement at pixel [y, x]"""
    d_ts = ts_data[:, yx[0], yx[1]]
    if inps.zero_first:
        d_ts -= d_ts[inps.zero_idx]

    # plot
    ax.cla()
    if inps.error_file:
        ax = plot_timeseries_errorbar(ax, d_ts, inps)
    else:
        ax = plot_timeseries_scatter(ax, d_ts, inps)

    # format
    ax = _adjust_ts_axis(ax, inps)
    title_ts = _get_ts_title(yx[0], yx[1], inps)
    if inps.disp_title:
        ax.set_title(title_ts)

    fig.canvas.draw()

    # Print to terminal
    print('\n---------------------------------------')
    print(title_ts)
    print(d_ts)

    # Slope estimation
    estimate_slope(d_ts, inps.yearList,
                   ex_flag=inps.ex_flag,
                   disp_unit=inps.disp_unit)

    return ax, d_ts


def _adjust_ts_axis(ax, inps):
    ax.tick_params(which='both', direction='in', labelsize=inps.font_size,
                   bottom=True, top=True, left=True, right=True)
    ax = pp.auto_adjust_xaxis_date(ax, inps.yearList, fontsize=inps.font_size)[0]
    ax.set_xlabel('Time', fontsize=inps.font_size)
    ax.set_ylabel('Displacement ({})'.format(inps.disp_unit), fontsize=inps.font_size)
    ax.set_ylim(inps.ylim)
    return ax


def _get_ts_title(y, x, inps):
    title = 'Y = {}, X = {}'.format(y, x)
    if inps.geocoded:
        lat = inps.lat0 + y*inps.lat_step
        lon = inps.lon0 + x*inps.lon_step
        title += ', lat = {:.4f}, lon = {:.4f}'.format(lat, lon)
    elif inps.lats is not None:
        lat = inps.lats[y, x]
        lon = inps.lons[y, x]
        title += ', lat = {:.4f}, lon = {:.4f}'.format(lat, lon)
    return title


def estimate_slope(d_ts, year_list, ex_flag=None, disp_unit='cm', print_msg=True):
    """Estimate linear velocity / STD of the crrent displacement time-series"""
    d_ts = np.array(d_ts)
    years = np.array(year_list)
    if ex_flag is not None:
        d_ts = d_ts[ex_flag == 1]
        years = years[ex_flag == 1]

    d_fit = stats.linregress(years, d_ts)
    vel = d_fit[0]
    std = d_fit[4]

    if print_msg:
        print('linear velocity: {v:.2f} +/- {s:.2f} [{u}/yr]'.format(
            v=vel, s=std, u=disp_unit))
    return vel, std


def save_ts_plot(yx, fig_v, fig_ts, ts_data, inps):
    print('save info on pixel ({}, {})'.format(yx[0], yx[1]))
    if not inps.fig_base:
        inps.fig_base = 'y%d_x%d' % (yx[0], yx[1])

    # read data
    d_ts = ts_data[:, yx[0], yx[1]]
    vel, std = estimate_slope(d_ts, inps.yearList,
                              ex_flag=inps.ex_flag,
                              disp_unit=inps.disp_unit,
                              print_msg=False)

    # TXT - point time series
    outName = '{}_ts.txt'.format(inps.fig_base)
    header_info = 'timeseries_file={}\n'.format(os.path.abspath(inps.timeseries_file))
    header_info += '{}\n'.format(_get_ts_title(yx[0], yx[1], inps))
    header_info += 'reference pixel: y={}, x={}\n'.format(inps.ref_yx[0], inps.ref_yx[1])
    header_info += 'reference date: {}\n'.format(inps.date_list[inps.ref_idx])
    header_info += 'unit: {}\n'.format(inps.disp_unit)
    header_info += 'slope: {:.2f} +/- {:.2f} [{}/yr]'.format(vel, std, inps.disp_unit)

    np.savetxt(outName,
               list(zip(np.array(inps.date_list), d_ts)),
               fmt='%s',
               delimiter='    ',
               header=header_info)
    print('save time series displacement in meter to '+outName)

    # Figure - point time series
    outName = '{}_ts.pdf'.format(inps.fig_base)
    fig_ts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    print('save time series plot to '+outName)

    # Figure - map
    outName = '{}_{}.png'.format(inps.fig_base, inps.date_list[inps.init_idx])
    fig_v.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    print('save map plot to '+outName)
    return


###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    inps, atr = read_init_info(inps)

    ts_data, mask, inps = read_timeseries_data(inps)

    # Figure 1 - Cumulative Displacement Map"""
    fig_v = plt.figure('Cumulative Displacement Map')
    # Axes 1
    ax_v = fig_v.add_axes([0.125, 0.25, 0.75, 0.65])
    ax_v, im = plot_init_map(ax_v, ts_data, inps)
    # Axes 2 - Time Slider
    ax_time = fig_v.add_axes([0.2, 0.1, 0.6, 0.07])
    tslider = plot_init_time_slider(ax=ax_time,
                                    year_list=inps.yearList,
                                    init_idx=inps.init_idx,
                                    ref_idx=inps.ref_idx)

    def update_time_slider(val):
        """Update Displacement Map using Slider"""
        idx = np.argmin(np.abs(np.array(inps.yearList) - tslider.val))
        disp_date = inps.dates[idx].strftime('%Y-%m-%d')
        ax_v.set_title('N = {n}, Time = {t}'.format(n=idx, t=disp_date))
        d_v = np.array(ts_data[idx, :, :])
        if inps.wrap:
            d_v *= inps.range2phase
            d_v -= np.round(d_v/(2*np.pi)) * (2*np.pi)
        im.set_data(d_v)
        fig_v.canvas.draw()
    tslider.on_changed(update_time_slider)

    # Figure 2 - Time Series Displacement - Point
    fig_ts, ax_ts = plt.subplots(num='Point Displacement Time-series', figsize=inps.fig_size)
    ax_ts, d_ts = plot_point_timeseries(inps.yx, fig_ts, ax_ts, ts_data, inps)

    def plot_timeseries_event(event):
        """Event function to get y/x from button press"""
        if event.inaxes == ax_v:
            y, x = int(event.ydata+0.5), int(event.xdata+0.5)
            if mask[y, x] != 0:
                d_ts = plot_point_timeseries((y, x), fig_ts, ax_ts, ts_data, inps)
            else:
                print('pixel is in masked out area.')

    # Output
    if inps.save_fig:
        save_ts_plot(inps.yx, fig_v, fig_ts, ts_data, inps)

    # Final linking of the canvas to the plots.
    cid = fig_v.canvas.mpl_connect('button_press_event', plot_timeseries_event)

    if inps.disp_fig:
        plt.show()

    fig_v.canvas.mpl_disconnect(cid)
    return


#########################################################################################
if __name__ == '__main__':
    main()
