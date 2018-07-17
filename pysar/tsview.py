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
from pysar import view


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
    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                        help='Y limits for point plotting.')

    parser.add_argument('-l','--lookup', dest='lookup_file', type=str,
                        help='lookup table file')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                       help='initial pixel to plot in lat/lon coord')

    pixel.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=8.0,
                       help='Point marker size. Default: 12.0')
    pixel.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0,
                       help='Edge width. Default: 1.0')

    parser.add_argument('-n', dest='init_idx', metavar='NUM', type=int,
                        help='Epoch/slice number to display.')
    parser.add_argument('--error', dest='error_file',
                        help='txt file with error for each date.')

    parser.add_argument('--exclude', '--ex', dest='ex_date_list', nargs='*',
                        help='Exclude date shown as gray.')
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                        help='Set displacement at first acquisition to zero.')

    parser = pp.add_data_disp_argument(parser)
    parser = pp.add_dem_argument(parser)
    parser = pp.add_figure_argument(parser)
    parser = pp.add_gps_argument(parser)
    parser = pp.add_mask_argument(parser)
    parser = pp.add_map_argument(parser)
    parser = pp.add_point_argument(parser)
    parser = pp.add_reference_argument(parser)
    parser = pp.add_save_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if (not inps.disp_fig or inps.outfile) and not inps.save_fig:
        inps.save_fig = True
    if inps.ylim:
        inps.ylim = sorted(inps.ylim)

    # default value
    if not inps.disp_unit:
        inps.disp_unit = 'cm'
    if not inps.colormap:
        inps.colormap = 'jet'
    if not inps.fig_size:
        inps.fig_size = [8.0, 4.5]

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

    # default lookup table file
    if not inps.lookup_file:
        inps.lookup_file = ut.get_lookup_file('./INPUTS/geometryRadar.h5')
    inps.coord = ut.coordinate(atr, inps.lookup_file)

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
            inps.init_idx = -3
        else:
            inps.init_idx = 3

    # Display Unit
    (inps.disp_unit,
     inps.unit_fac) = pp.scale_data2disp_unit(metadata=atr,
                                              disp_unit=inps.disp_unit)[1:3]

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

    # size and lalo info
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

    inps.pix_box = (0, 0, inps.width, inps.length)
    inps.geo_box = inps.coord.box_pixel2geo(inps.pix_box)

    # reference pixel
    if not inps.ref_lalo and 'REF_LAT' in atr.keys():
        inps.ref_lalo = (float(atr['REF_LAT']), float(atr['REF_LON']))
    if inps.ref_lalo:
        inps.ref_yx = inps.coord.geo2radar(inps.ref_lalo[0],
                                           inps.ref_lalo[1],
                                           print_msg=False)[0:2]
    if not inps.ref_yx:
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]

    # Initial Pixel Coord
    if inps.lalo:
        inps.yx = inps.coord.geo2radar(inps.lalo[0],
                                       inps.lalo[1],
                                       print_msg=False)[0:2]
    if not inps.yx:
        inps.yx = inps.ref_yx
    inps.lalo = inps.coord.radar2geo(inps.yx[0],
                                     inps.yx[1],
                                     print_msg=False)[0:2]

    # Flip up-down / left-right
    if not inps.flip_lr and not inps.flip_ud:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(atr)

    inps.disp_unit_v = inps.disp_unit
    if inps.wrap:
        inps.range2phase = -4. * np.pi / float(atr['WAVELENGTH'])
        if   'cm' in inps.disp_unit:   inps.range2phase /= 100.
        elif 'mm' in inps.disp_unit:   inps.range2phase /= 1000.
        inps.disp_unit_v = 'radian'
        inps.vlim = [-np.pi, np.pi]
    inps.cbar_label = 'Displacement ({})'.format(inps.disp_unit_v)

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
    inps.dlim = [np.nanmin(ts_data), np.nanmax(ts_data)]
    ts_data_mli = multilook_data(ts_data, 10, 10)
    if not inps.vlim:
        inps.vlim = [np.nanmin(ts_data_mli), np.nanmax(ts_data_mli)]
    print('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
    print('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    # default ylim
    if not inps.ylim:
        if inps.zero_first:
            ts_data_mli -= np.tile(ts_data_mli[inps.zero_idx, :, :], (inps.num_date, 1, 1))
        ymin, ymax = np.nanmin(ts_data_mli), np.nanmax(ts_data_mli)
        ybuffer = (ymax - ymin) * 0.05
        inps.ylim = [ymin - ybuffer, ymax + ybuffer]
    del ts_data_mli

    return ts_data, mask, inps


def plot_init_map(ax, d_v, inps, metadata):

    # prepare data
    if inps.wrap:
        d_v *= inps.range2phase
        d_v -= np.round(d_v/(2*np.pi)) * (2*np.pi)

    # Title and Axis Label
    disp_date = inps.dates[inps.init_idx].strftime('%Y-%m-%d')
    inps.fig_title = 'N = {}, Time = {}'.format(inps.init_idx, disp_date)

    # Initial Pixel
    if inps.yx != inps.ref_yx:
        inps.pts_yx = np.array(inps.yx).reshape(-1, 2)
        inps.pts_lalo = np.array(inps.lalo).reshape(-1, 2)
        inps.pts_marker = 'ro'

    # call view.py to plot
    ax, inps, im, cbar = view.plot_2d_matrix(ax, d_v, metadata, inps)

    return ax, im


def plot_init_time_slider(ax, year_list, init_idx=-1, ref_idx=0):
    val_step = np.min(np.diff(year_list))
    tslider = Slider(ax, label='Years',
                     valmin=year_list[0]-val_step,
                     valmax=year_list[-1]+val_step,
                     valinit=year_list[init_idx],
                     valstep=val_step)

    bar_width = val_step / 4.
    datex = np.array(year_list) - bar_width / 2.
    tslider.ax.bar(datex, np.ones(len(datex)), bar_width, facecolor='black', ecolor=None)
    tslider.ax.bar(datex[ref_idx], 1., bar_width*3, facecolor='crimson', ecolor=None)

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
        ax.errorbar(inps.ex_dates, ex_d_ts, yerr=inps.ex_error_ts,
                    fmt='-o', color='gray', ms=inps.marker_size,
                    lw=0, alpha=1, mfc='gray',
                    elinewidth=inps.edge_width, ecolor='black',
                    capsize=inps.marker_size*0.5, mew=inps.edge_width)

    # Plot kept dates
    ax.errorbar(dates, d_ts, yerr=inps.error_ts,
                fmt='-o', ms=inps.marker_size,
                lw=0, alpha=1,
                elinewidth=inps.edge_width, ecolor='black',
                capsize=inps.marker_size*0.5, mew=inps.edge_width)
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
    title_ts = _get_ts_title(yx[0], yx[1], inps.coord)
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


def _get_ts_title(y, x, coord):
    title = 'Y = {}, X = {}'.format(y, x)
    try:
        lat, lon = coord.radar2geo(y, x, print_msg=False)[0:2]
        title += ', lat = {:.4f}, lon = {:.4f}'.format(lat, lon)
    except:
        pass
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
    if not inps.outfile:
        inps.outfile = 'y%d_x%d' % (yx[0], yx[1])

    # read data
    d_ts = ts_data[:, yx[0], yx[1]]
    vel, std = estimate_slope(d_ts, inps.yearList,
                              ex_flag=inps.ex_flag,
                              disp_unit=inps.disp_unit,
                              print_msg=False)

    # TXT - point time series
    outName = '{}_ts.txt'.format(inps.outfile)
    header_info = 'timeseries_file={}\n'.format(os.path.abspath(inps.timeseries_file))
    header_info += '{}\n'.format(_get_ts_title(yx[0], yx[1], inps.coord))
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
    outName = '{}_ts.pdf'.format(inps.outfile)
    fig_ts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    print('save time series plot to '+outName)

    # Figure - map
    outName = '{}_{}.png'.format(inps.outfile, inps.date_list[inps.init_idx])
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
    d_v = np.array(ts_data[inps.init_idx, :, :])
    ax_v, im = plot_init_map(ax_v, d_v, inps, atr)

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
            # get row/col number
            if inps.fig_coord == 'geo':
                y, x = inps.coord.geo2radar(event.ydata,
                                            event.xdata,
                                            print_msg=False)[0:2]
            else:
                y, x = int(event.ydata+0.5), int(event.xdata+0.5)

            # plot time-series displacement if selected pixel is valid
            if mask[y, x] != 0:
                d_ts = plot_point_timeseries((y, x), fig_ts, ax_ts, ts_data, inps)
            else:
                print('pixel ({}, {}) is in masked out area.'.format(y, x))

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
