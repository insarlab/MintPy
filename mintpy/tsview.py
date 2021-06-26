#!/usr/bin/env python3
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Zhang Yunjun, Joshua Zahner, Heresh Fattahi, 2013 #
#############################################################


import os
import sys
import argparse
import numpy as np
from scipy import linalg
from datetime import datetime
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.widgets import Slider

from mintpy.objects import timeseries, giantTimeseries, HDFEOS
from mintpy.utils import arg_group, readfile, ptime, plot as pp, utils as ut
from mintpy.multilook import multilook_data
from mintpy import subset, view, timeseries2velocity as ts2vel


###########################################################################################
EXAMPLE = """example:
  tsview.py timeseries.h5
  tsview.py timeseries.h5  --wrap
  tsview.py timeseries.h5  --yx 300 400 --zero-first  --nodisplay
  tsview.py geo_timeseries.h5  --lalo 33.250 131.665  --nodisplay
  tsview.py timeseries_ERA5_ramp_demErr.h5  --sub-x 900 1400 --sub-y 0 500

  # press left / right key to slide images

  # multiple time-series files
  tsview.py timeseries_ERA5_ramp_demErr.h5 timeseries_ERA5_ramp.h5 timeseries_ERA5.h5 timeseries.h5 --off 5
  tsview.py timeseries_ERA5_ramp_demErr.h5 ../GIANT/Stack/LS-PARAMS.h5 --off 5 --label mintpy giant
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Interactive time-series viewer',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',
                        help='time-series file to display\n'
                             'i.e.: timeseries_ERA5_ramp_demErr.h5 (MintPy)\n'
                             '      LS-PARAMS.h5 (GIAnT)\n'
                             '      S1_IW12_128_0593_0597_20141213_20180619.he5 (HDF-EOS5)')
    parser.add_argument('--label', dest='file_label', nargs='*', help='labels to display for multiple input files')
    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float, help='Y limits for point plotting.')
    parser.add_argument('--tick-right', dest='tick_right', action='store_true', help='set tick and tick label to the right')
    parser.add_argument('-l','--lookup', dest='lookup_file', type=str, help='lookup table file')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2, help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2, help='initial pixel to plot in lat/lon coord')

    pixel.add_argument('--marker', type=str, default='o', help='marker style (default: %(default)s).')
    pixel.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=6.0, help='marker size (default: %(default)s).')
    pixel.add_argument('--lw', '--linewidth', dest='linewidth', type=float, default=0, help='line width (default: %(default)s).')
    pixel.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0, help='Edge width for the error bar (default: %(default)s)')

    parser.add_argument('-n', dest='idx', metavar='NUM', type=int, help='Epoch/slice number for initial display.')
    parser.add_argument('--error', dest='error_file', help='txt file with error for each date.')

    parser.add_argument('--start-date', dest='start_date', type=str, help='start date of displacement to display')
    parser.add_argument('--end-date', dest='end_date', type=str, help='end date of displacement to display')
    parser.add_argument('--exclude', '--ex', dest='ex_date_list', nargs='*', default=['exclude_date.txt'], help='Exclude date shown as gray.')
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true', help='Set displacement at first acquisition to zero.')
    parser.add_argument('--off','--offset', dest='offset', type=float, help='Offset for each timeseries file.')

    parser.add_argument('--noverbose', dest='print_msg', action='store_false', help='Disable the verbose message printing.')
    parser.add_argument('--savetxt','--txt', dest='savetxt', type=str, default=False, help='Filename of a text file to save timeseries (and model params) (default: %(default)s).')

    parser = arg_group.add_data_disp_argument(parser)
    parser = arg_group.add_dem_argument(parser)
    parser = arg_group.add_figure_argument(parser)
    parser = arg_group.add_gps_argument(parser)
    parser = arg_group.add_mask_argument(parser)
    parser = arg_group.add_map_argument(parser)
    parser = arg_group.add_reference_argument(parser)
    parser = arg_group.add_save_argument(parser)
    parser = arg_group.add_subset_argument(parser)

    # temporal model fitting
    parser.add_argument('--show-model', '--model', dest='show_model', action='store_true', help='Show the ts2vel model fitting.')
    parser = arg_group.add_timefunc_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if '--gps-comp' in iargs:
        msg = '--gps-comp is not supported for {}'.format(os.path.basename(__file__))
        raise NotImplementedError(msg)

    if inps.file_label:
        if len(inps.file_label) != len(inps.file):
            raise Exception('input number of labels != number of files.')

    if (not inps.disp_fig or inps.outfile) and not inps.save_fig:
        inps.save_fig = True
    if inps.ylim:
        inps.ylim = sorted(inps.ylim)
    if inps.zero_mask:
        inps.mask_file = 'no'

    # default value
    if not inps.disp_unit:
        inps.disp_unit = 'cm'
    if not inps.colormap:
        inps.colormap = 'jet'
    if not inps.fig_size:
        inps.fig_size = [8.0, 4.5]

    # temporal model fitting, initialize the dicts of exp and log funcs
    inps = ts2vel.init_explog_dicts(inps)

    # verbose print using --noverbose option
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


###########################################################################################
def read_init_info(inps):
    # Time Series Info
    ts_file0 = inps.file[0]
    atr = readfile.read_attribute(ts_file0)
    inps.key = atr['FILE_TYPE']
    if inps.key == 'timeseries':
        obj = timeseries(ts_file0)
    elif inps.key == 'giantTimeseries':
        obj = giantTimeseries(ts_file0)
    elif inps.key == 'HDFEOS':
        obj = HDFEOS(ts_file0)
    else:
        raise ValueError('input file is {}, not timeseries.'.format(inps.key))
    obj.open(print_msg=inps.print_msg)

    if not inps.file_label:
        #inps.file_label = [str(i) for i in list(range(len(inps.file)))]
        inps.file_label = []
        for fname in inps.file:
            fbase = os.path.splitext(os.path.basename(fname))[0]
            fbase = fbase.replace('timeseries', '')
            inps.file_label.append(fbase)

    # default mask file
    if not inps.mask_file and 'masked' not in ts_file0:
        dir_name = os.path.dirname(ts_file0)
        if 'Y_FIRST' in atr.keys():
            inps.mask_file = os.path.join(dir_name, 'geo_maskTempCoh.h5')
        else:
            inps.mask_file = os.path.join(dir_name, 'maskTempCoh.h5')
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None

    # date info
    inps.date_list = obj.dateList
    inps.num_date = len(inps.date_list)
    if inps.start_date:
        inps.date_list = [i for i in inps.date_list if int(i) >= int(inps.start_date)]
    if inps.end_date:
        inps.date_list = [i for i in inps.date_list if int(i) <= int(inps.end_date)]
    inps.num_date = len(inps.date_list)
    inps.dates, inps.yearList = ptime.date_list2vector(inps.date_list)
    (inps.ex_date_list,
     inps.ex_dates,
     inps.ex_flag) = read_exclude_date(inps.ex_date_list, inps.date_list)

    # reference date/index
    if not inps.ref_date:
        inps.ref_date = atr.get('REF_DATE', None)
    if inps.ref_date:
        inps.ref_idx = inps.date_list.index(inps.ref_date)
    else:
        inps.ref_idx = None

    # date/index of interest for initial display
    if not inps.idx:
        if (not inps.ref_idx) or (inps.ref_idx < inps.num_date / 2.):
            inps.idx = inps.num_date - 2
        else:
            inps.idx = 2

    # Display Unit
    (inps.disp_unit,
     inps.unit_fac) = pp.scale_data2disp_unit(metadata=atr, disp_unit=inps.disp_unit)[1:3]

    # Map info - coordinate unit
    inps.coord_unit = atr.get('Y_UNIT', 'degrees').lower()

    # Read Error List
    inps.ts_plot_func = plot_ts_scatter
    inps.error_ts = None
    inps.ex_error_ts = None
    if inps.error_file:
        # assign plot function
        inps.ts_plot_func = plot_ts_errorbar
        # read error file
        error_fc = np.loadtxt(inps.error_file, dtype=bytes).astype(str)
        inps.error_ts = error_fc[:, 1].astype(np.float)*inps.unit_fac
        # update error file with exlcude date
        if inps.ex_date_list:
            e_ts = inps.error_ts[:]
            inps.ex_error_ts = e_ts[inps.ex_flag == 0]
            inps.error_ts = e_ts[inps.ex_flag == 1]

    # Zero displacement for 1st acquisition
    if inps.zero_first:
        inps.zero_idx = min(0, np.min(np.where(inps.ex_flag)[0]))

    # default lookup table file
    if not inps.lookup_file:
        inps.lookup_file = ut.get_lookup_file('./inputs/geometryRadar.h5')
    inps.coord = ut.coordinate(atr, inps.lookup_file)

    # size and lalo info
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), atr)
    inps.pix_box = inps.coord.check_box_within_data_coverage(inps.pix_box)
    inps.geo_box = inps.coord.box_pixel2geo(inps.pix_box)
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    vprint('data   coverage in y/x: '+str(data_box))
    vprint('subset coverage in y/x: '+str(inps.pix_box))
    vprint('data   coverage in lat/lon: '+str(inps.coord.box_pixel2geo(data_box)))
    vprint('subset coverage in lat/lon: '+str(inps.geo_box))
    vprint('------------------------------------------------------------------------')

    # reference pixel
    if not inps.ref_lalo and 'REF_LAT' in atr.keys():
        inps.ref_lalo = (float(atr['REF_LAT']), float(atr['REF_LON']))
    if inps.ref_lalo:
        # set longitude to [-180, 180)
        if inps.coord_unit.lower().startswith('deg') and inps.ref_lalo[1] >= 180.:
            inps.ref_lalo[1] -= 360.
        # ref_lalo --> ref_yx if not set in cmd
        if not inps.ref_yx:
            inps.ref_yx = inps.coord.geo2radar(inps.ref_lalo[0],
                                               inps.ref_lalo[1],
                                               print_msg=False)[0:2]
    # use REF_Y/X if ref_yx not set in cmd
    if not inps.ref_yx and 'REF_Y' in atr.keys():
        inps.ref_yx = (int(atr['REF_Y']), int(atr['REF_X']))

    # ref_yx --> ref_lalo if in geo-coord
    if inps.ref_yx and 'Y_FIRST' in atr.keys():
        inps.ref_lalo = inps.coord.radar2geo(inps.ref_yx[0],
                                             inps.ref_yx[1],
                                             print_msg=False)[0:2]

    # do not plot native reference point if it's out of the coverage due to subset
    if (inps.ref_yx and 'Y_FIRST' in atr.keys()
        and inps.ref_yx == (int(atr['REF_Y']), int(atr['REF_X']))
        and not (inps.pix_box[0] <= inps.ref_yx[1] < inps.pix_box[2]
                 and inps.pix_box[1] <= inps.ref_yx[0] < inps.pix_box[3])):
        inps.disp_ref_pixel = False
        print('the native REF_Y/X is out of subset box, thus do not display')

    # Initial Pixel Coord
    if inps.lalo:
        inps.yx = inps.coord.geo2radar(inps.lalo[0], inps.lalo[1], print_msg=False)[0:2]
    try:
        inps.lalo = inps.coord.radar2geo(inps.yx[0], inps.yx[1], print_msg=False)[0:2]
    except:
        inps.lalo = None

    # Flip up-down / left-right
    if inps.auto_flip:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(atr, print_msg=inps.print_msg)

    # Transparency - Alpha
    if not inps.transparency:
        # Auto adjust transparency value when showing shaded relief DEM
        if inps.dem_file and inps.disp_dem_shade:
            inps.transparency = 0.7
        else:
            inps.transparency = 1.0

    # display unit ans wrap
    # if wrap_step == 2*np.pi (default value), set disp_unit_img = radian;
    # otherwise set disp_unit_img = disp_unit
    inps.disp_unit_img = inps.disp_unit
    if inps.wrap:
        inps.range2phase = -4. * np.pi / float(atr['WAVELENGTH'])
        if   'cm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 100.
        elif 'mm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1000.
        elif 'm'  == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1.
        else:
            raise ValueError('un-recognized display unit: {}'.format(inps.disp_unit))

        if (inps.wrap_range[1] - inps.wrap_range[0]) == 2*np.pi:
            inps.disp_unit_img = 'radian'
        inps.vlim = inps.wrap_range
    inps.cbar_label = 'Displacement [{}]'.format(inps.disp_unit_img)

    # whether to fit a velocity model to the time series
    if inps.show_model:
        inps.dateList = inps.date_list
        inps.model, inps.num_param = ts2vel.read_inps2model(inps)
    return inps, atr


def read_exclude_date(input_ex_date, dateListAll):
    """
    Parameters: input_ex_date : list of string in YYYYMMDD or filenames for excluded dates
                dateListAll   : list of string in YYYYMMDD for all dates
    Returns:    ex_date_list  : list of string in YYYYMMDD for excluded dates
                ex_dates      : list of datetime.datetime objects for excluded dates
                ex_flag       : 1D np.ndarray in size of (num_date),
                                1/True for kept, 0/False for excluded
    """
    # default value
    ex_date_list = []
    ex_dates = []
    ex_flag = np.ones((len(dateListAll)), np.bool_)

    ex_date_list = ptime.read_date_list(input_ex_date, date_list_all=dateListAll)
    if ex_date_list:
        ex_dates = ptime.date_list2vector(ex_date_list)[0]
        for i in ex_date_list:
            ex_flag[dateListAll.index(i)] = False
        vprint('exclude date:'+str(ex_date_list))
    return ex_date_list, ex_dates, ex_flag


def read_timeseries_data(inps):
    """Read data of time-series files
    Parameters: inps : Namespace of input arguments
    Returns:    ts_data : list of 3D np.array in size of (num_date, length, width)
                mask : 2D np.array in size of (length, width)
                inps : Namespace of input arguments
    """
    # read list of 3D time-series
    ts_data = []
    for fname in inps.file:
        vprint('reading timeseries from file {} ...'.format(fname))
        data, atr = readfile.read(fname, datasetName=inps.date_list, box=inps.pix_box)
        if inps.ref_yx and inps.ref_yx != (int(atr.get('REF_Y', -1)), int(atr.get('REF_X', -1))):
            ref_phase = data[:, inps.ref_yx[0]-inps.pix_box[1], inps.ref_yx[1]-inps.pix_box[0]]
            data -= np.tile(ref_phase.reshape(-1, 1, 1), (1, data.shape[-2], data.shape[-1]))
            vprint('reference to pixel: {}'.format(inps.ref_yx))
        if inps.ref_idx is not None:
            vprint('reference to date: {}'.format(inps.date_list[inps.ref_idx]))
            data -= np.tile(data[inps.ref_idx, :, :], (inps.num_date, 1, 1))

        # Display Unit
        (data,
         inps.disp_unit,
         inps.unit_fac) = pp.scale_data2disp_unit(data,
                                                  metadata=atr,
                                                  disp_unit=inps.disp_unit)
        ts_data.append(data)

    # Mask file: input mask file + non-zero ts pixels - ref_point
    mask = np.ones(ts_data[0].shape[-2:], np.bool_)
    msk = pp.read_mask(inps.file[0],
                       mask_file=inps.mask_file,
                       datasetName='displacement',
                       box=inps.pix_box,
                       print_msg=inps.print_msg)[0]
    mask[msk == 0.] = False
    del msk

    ts_stack = np.nansum(ts_data[0], axis=0)
    mask[np.isnan(ts_stack)] = False
    # keep all-zero value for unwrapError time-series
    if atr['UNIT'] not in ['cycle']:
        mask[ts_stack == 0.] = False
    del ts_stack

    #do not mask the reference point
    try:
        mask[inps.ref_yx[0]-inps.pix_box[1],
             inps.ref_yx[1]-inps.pix_box[0]] = True
    except:
        pass

    # default vlim
    inps.dlim = [np.nanmin(ts_data[0]), np.nanmax(ts_data[0])]
    if not inps.vlim:
        inps.cmap_lut, inps.vlim = pp.auto_adjust_colormap_lut_and_disp_limit(ts_data[0],
                                                                              num_multilook=10,
                                                                              print_msg=inps.print_msg)
    vprint('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
    vprint('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    # default ylim
    num_file = len(inps.file)
    if not inps.ylim:
        ts_data_mli = multilook_data(np.squeeze(ts_data[-1]), 4, 4)
        if inps.zero_first:
            ts_data_mli -= np.tile(ts_data_mli[inps.zero_idx, :, :], (inps.num_date, 1, 1))
        ymin, ymax = (np.nanmin(ts_data_mli[inps.ex_flag != 0]),
                      np.nanmax(ts_data_mli[inps.ex_flag != 0]))
        ybuffer = (ymax - ymin) * 0.05
        inps.ylim = [ymin - ybuffer, ymax + ybuffer]
        if inps.offset:
            inps.ylim[1] += inps.offset * (num_file - 1)
        del ts_data_mli

    return ts_data, mask, inps


def plot_ts_errorbar(ax, dis_ts, inps, ppar):
    # make a local copy
    dates = np.array(inps.dates)
    d_ts = dis_ts[:]

    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        ax.errorbar(inps.ex_dates, ex_d_ts, yerr=inps.ex_error_ts,
                    fmt='-o', ms=ppar.ms, color='gray',
                    lw=0, alpha=1, mfc='gray',
                    elinewidth=inps.edge_width, ecolor='black',
                    capsize=ppar.ms*0.5, mew=inps.edge_width)

    # Plot kept dates
    ax.errorbar(dates, d_ts, yerr=inps.error_ts,
                fmt='-o', ms=ppar.ms, color=ppar.mfc,
                lw=0, alpha=1,
                elinewidth=inps.edge_width, ecolor='black',
                capsize=ppar.ms*0.5, mew=inps.edge_width,
                label=ppar.label)
    return ax


def plot_ts_scatter(ax, dis_ts, inps, ppar):
    # make a local copy
    dates = np.array(inps.dates)
    d_ts = dis_ts[:]

    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        ax.plot(inps.ex_dates, ex_d_ts, color='gray', lw=0, ms=ppar.ms, marker=inps.marker)

    # Plot kept dates
    ax.plot(dates, d_ts, color=ppar.mfc, label=ppar.label, lw=inps.linewidth, ms=ppar.ms, marker=inps.marker)
    return ax

def plot_ts_pred_curve(ax, dis_ts, inps, ppar):
    # make a local copy
    dates = np.array(inps.dates)
    dates_str = []
    for dt in dates:
        dates_str.append(dt.strftime("%Y%m%d"))
    dates_str = ptime.date2arange(dates_str[0], dates_str[-1])
    dates = [datetime.strptime(dd, '%Y%m%d') for dd in dates_str]
    d_ts = dis_ts[:]

    # plot the model prediction curve in a fine date_lists
    ax.plot(dates, d_ts, color=ppar.color, label=ppar.label, lw=ppar.linewidth, alpha=0.6)

    # print model parameters on the plot
    msg = '\n'.join(inps.model_table)
    ax.text(0.04, 0.88, msg, size=10, ha='left', va='center', transform=ax.transAxes)
    return ax

def _adjust_ts_axis(ax, inps):
    ax.tick_params(which='both', direction='in', labelsize=inps.font_size, bottom=True, top=True, left=True, right=True)
    ax = pp.auto_adjust_xaxis_date(ax, inps.yearList, fontsize=inps.font_size)[0]
    #ax.set_xlabel('Time [years]', fontsize=inps.font_size)
    ax.set_ylabel('Displacement [{}]'.format(inps.disp_unit), fontsize=inps.font_size)
    ax.set_ylim(inps.ylim)
    return ax


def _get_ts_title(y, x, coord):
    title = 'Y/X = {}, {}'.format(y, x)
    try:
        lat, lon = coord.radar2geo(y, x, print_msg=False)[0:2]
        title += ', lat/lon = {:.4f}, {:.4f}'.format(lat, lon)
    except:
        pass
    return title


def estimate_slope(d_ts, year_list, ex_flag=None, disp_unit='cm'):
    """Estimate linear velocity / STD of the crrent displacement time-series"""
    d_ts = np.array(d_ts)
    years = np.array(year_list)
    if ex_flag is not None:
        d_ts = d_ts[ex_flag == 1]
        years = years[ex_flag == 1]

    d_fit = stats.linregress(years, d_ts)
    vel = d_fit[0]
    std = d_fit[4]

    vprint('linear velocity: {v:.2f} +/- {s:.2f} [{u}/yr]'.format(v=vel, s=std, u=disp_unit))
    return vel, std


def save_ts_plot(yx, fig_img, fig_pts, d_ts, inps):
    vprint('save info on pixel ({}, {})'.format(yx[0], yx[1]))
    # output file name
    if inps.outfile:
        inps.outfile_base, ext = os.path.splitext(inps.outfile[0])
        if ext != '.pdf':
            vprint(('Output file extension is fixed to .pdf,'
                    ' input extension {} is ignored.').format(ext))
    else:
        inps.outfile_base = 'y{}x{}'.format(yx[0], yx[1])

    # get aux info
    vel, std = estimate_slope(d_ts[0], inps.yearList,
                              ex_flag=inps.ex_flag,
                              disp_unit=inps.disp_unit)

    # TXT - point time-series
    outName = '{}_ts.txt'.format(inps.outfile_base)
    header_info = 'time-series file={}\n'.format(inps.file)
    header_info += '{}\n'.format(_get_ts_title(yx[0], yx[1], inps.coord))
    try:
        header_info += 'reference pixel: y={}, x={}\n'.format(inps.ref_yx[0], inps.ref_yx[1])
        header_info += 'reference date: {}\n'.format(inps.date_list[inps.ref_idx])
    except:
        pass
    header_info += 'unit: {}\n'.format(inps.disp_unit)
    header_info += 'slope: {:.2f} +/- {:.2f} [{}/yr]'.format(vel, std, inps.disp_unit)

    # prepare data
    data = np.array(inps.date_list).reshape(-1, 1)
    for i in range(len(d_ts)):
        data = np.hstack((data, d_ts[i].reshape(-1, 1)))
    # write
    np.savetxt(outName,
               data,
               fmt='%s',
               delimiter='\t',
               header=header_info)
    vprint('save displacement time-series in meter to '+outName)

    # Figure - point time-series
    outName = '{}_ts.pdf'.format(inps.outfile_base)
    fig_pts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    vprint('save time-series plot to '+outName)

    # Figure - map
    outName = '{}_{}.png'.format(inps.outfile_base, inps.date_list[inps.idx])
    fig_img.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    vprint('save map plot to '+outName)
    return


class timeseriesViewer():
    """Class for tsview.py

    Example:
        cmd = 'tsview.py timeseries_ERA5_ramp_demErr.h5'
        obj = timeseriesViewer(cmd)
        obj.configure()
        obj.plot()
    """

    def __init__(self, cmd=None, iargs=None):
        if cmd:
            iargs = cmd.split()[1:]
        self.cmd = cmd
        self.iargs = iargs
        # print command line
        cmd = '{} '.format(os.path.basename(__file__))
        cmd += ' '.join(iargs)
        print(cmd)

        # figure variables
        self.figname_img = 'Cumulative Displacement Map'
        self.figsize_img = None
        self.fig_img = None
        self.ax_img = None
        self.cbar_img = None
        self.img = None

        self.ax_tslider = None
        self.tslider = None

        self.figname_pts = 'Point Displacement Time-series'
        self.figsize_pts = None
        self.fig_pts = None
        self.ax_pts = None
        return


    def configure(self):
        inps = cmd_line_parse(self.iargs)
        inps, self.atr = read_init_info(inps)
        # copy inps to self object
        for key, value in inps.__dict__.items():
            setattr(self, key, value)
        # input figsize for the point time-series plot
        self.figsize_pts = self.fig_size
        self.pts_marker = 'r^'
        self.pts_marker_size = 6.
        return


    def plot(self):
        # read 3D time-series
        self.ts_data, self.mask = read_timeseries_data(self)[0:2]

        # fit a velocity model to time series
        if self.show_model:
            if self.yx:
                y = self.yx[0] - self.pix_box[1]
                x = self.yx[1] - self.pix_box[0]
                self.m, self.m_std, self.ts_pred, self.ts_fit = self.fit2model(self.ts_data[0][:,y,x])
                self.model_table = self.model_param_table()

        # save plotted data to text file
        if self.savetxt:
            f = open(self.savetxt,'w')
            print('open a txt file: {}'.format(self.savetxt))
            f.close()

        # Figure 1 - Cumulative Displacement Map
        self.figsize_img = pp.auto_figure_size(ds_shape=self.ts_data[0].shape[-2:],
                                               disp_cbar=True,
                                               disp_slider=True,
                                               print_msg=self.print_msg)
        self.fig_img = plt.figure(self.figname_img, figsize=self.figsize_img)

        # Figure 1 - Axes 1 - Displacement Map
        self.ax_img = self.fig_img.add_axes([0.125, 0.25, 0.75, 0.65])
        img_data = np.array(self.ts_data[0][self.idx, :, :])  ####################
        img_data[self.mask == 0] = np.nan
        self.plot_init_image(img_data)

        # Figure 1 - Axes 2 - Time Slider
        self.ax_tslider = self.fig_img.add_axes([0.125, 0.1, 0.75, 0.07])
        self.plot_init_time_slider(init_idx=self.idx, ref_idx=self.ref_idx)
        self.tslider.on_changed(self.update_time_slider)

        # Figure 2 - Time Series Displacement - Point
        self.fig_pts, self.ax_pts = plt.subplots(num=self.figname_pts, figsize=self.figsize_pts)
        if self.yx:
            d_ts = self.plot_point_timeseries(self.yx)

        # Output
        if self.save_fig:
            save_ts_plot(self.yx, self.fig_img, self.fig_pts, d_ts, self)

        # Final linking of the canvas to the plots.
        self.fig_img.canvas.mpl_connect('button_press_event', self.update_plot_timeseries)
        self.fig_img.canvas.mpl_connect('key_press_event', self.on_key_event)
        if self.disp_fig:
            vprint('showing ...')
            msg = '\n------------------------------------------------------------------------'
            msg += '\nTo scroll through the image sequence:'
            msg += '\n1) Move the slider, OR'
            msg += '\n2) Press left or right arrow key (if not responding, click the image and try again).'
            msg += '\n------------------------------------------------------------------------'
            vprint(msg)
            plt.show()
        return


    def plot_init_image(self, img_data):
        # prepare data
        if self.wrap:
            if self.disp_unit_img == 'radian':
                img_data *= self.range2phase
            img_data = ut.wrap(img_data, wrap_range=self.wrap_range)

        # Title and Axis Label
        disp_date = self.dates[self.idx].strftime('%Y-%m-%d')
        self.fig_title = 'N = {}, Time = {}'.format(self.idx, disp_date)

        # Initial Pixel of interest
        self.pts_yx = None
        self.pts_lalo = None
        if self.yx and self.yx != self.ref_yx:
            self.pts_yx = np.array(self.yx).reshape(-1, 2)
            if self.lalo:
                self.pts_lalo = np.array(self.lalo).reshape(-1, 2)

        # call view.py to plot
        self.img, self.cbar_img = view.plot_slice(self.ax_img, img_data, self.atr, self)[2:4]
        return self.img, self.cbar_img


    def plot_init_time_slider(self, init_idx=-1, ref_idx=None):
        val_step = np.min(np.diff(self.yearList))
        val_min = self.yearList[0]
        val_max = self.yearList[-1]

        self.tslider = Slider(self.ax_tslider, label='Time',
                              valinit=self.yearList[init_idx],
                              valmin=val_min,
                              valmax=val_max,
                              valstep=val_step)

        bar_width = val_step / 4.
        datex = np.array(self.yearList) - bar_width / 2.
        self.tslider.ax.bar(datex, np.ones(len(datex)), bar_width, facecolor='black', ecolor=None)
        if ref_idx is not None:
            self.tslider.ax.bar(datex[ref_idx], 1., bar_width*3, facecolor='crimson', ecolor=None)

        # xaxis tick format
        if np.floor(val_max) == np.floor(val_min):
            digit = 10.
        else:
            digit = 1.
        self.tslider.ax.set_xticks(np.round(np.linspace(val_min, val_max, num=5) * digit) / digit)
        self.tslider.ax.xaxis.set_minor_locator(MultipleLocator(1./12.))
        self.tslider.ax.set_xlim([val_min, val_max])
        self.tslider.ax.set_yticks([])
        self.tslider.valtext.set_visible(False)   #hide slider values
        return self.tslider


    def update_time_slider(self, val):
        """Update Displacement Map using Slider"""
        idx = np.argmin(np.abs(np.array(self.yearList) - self.tslider.val))
        # update title
        disp_date = self.dates[idx].strftime('%Y-%m-%d')
        sub_title = 'N = {n}, Time = {t}'.format(n=idx, t=disp_date)
        self.ax_img.set_title(sub_title, fontsize=self.font_size)

        # read data
        data_img = np.array(self.ts_data[0][idx, :, :])
        data_img[self.mask == 0] = np.nan
        if self.wrap:
            if self.disp_unit_img == 'radian':
                data_img *= self.range2phase
            data_img = ut.wrap(data_img, wrap_range=self.wrap_range)

        # update data
        self.img.set_data(data_img)
        self.idx = idx
        self.fig_img.canvas.draw()
        return

    def plot_point_timeseries(self, yx):
        """Plot point displacement time-series at pixel [y, x]
        Parameters: yx : list of 2 int
        Returns:    d_ts : 2D np.array in size of (num_date, num_file)
        """
        self.ax_pts.cla()

        # plot scatter in different size for different files
        num_file = len(self.ts_data)
        if   num_file <= 2: ms_step = 4
        elif num_file == 3: ms_step = 3
        elif num_file == 4: ms_step = 2
        elif num_file >= 5: ms_step = 1

        d_ts = []
        y = yx[0] - self.pix_box[1]
        x = yx[1] - self.pix_box[0]
        for i in range(num_file-1, -1, -1):
            # get displacement data
            d_tsi = self.ts_data[i][:, y, x]
            # plot fitting
            if self.show_model:
                d_tsi_p = self.ts_pred
                d_tsi_f = self.ts_fit
                if self.zero_first:
                    d_tsi -= d_tsi[self.zero_idx]
                    d_tsi_p -= d_tsi_p[self.zero_idx]
                    d_tsi_f -= d_tsi_f[self.zero_idx]
            d_ts.append(d_tsi)

            # get plot parameter - namespace ppar
            ppar = argparse.Namespace()
            ppar.label = self.file_label[i]
            ppar.ms = self.marker_size - ms_step * (num_file - 1 - i)
            ppar.mfc = pp.mplColors[num_file - 1 - i]
            if self.mask[y, x] == 0:
                ppar.mfc = 'gray'
            if self.offset:
                d_tsi += self.offset * (num_file - 1 - i)
                # offset the model prediction
                if self.show_model:
                    d_tsi_p += self.offset * (num_file - 1 - i)
                    d_tsi_f += self.offset * (num_file - 1 - i)

            # plot
            if not np.all(np.isnan(d_tsi)):
                self.ax_pts = self.ts_plot_func(self.ax_pts, d_tsi, self, ppar)
                # plot model prediction
                if self.show_model:
                    import copy
                    ppar_p     = copy.deepcopy(ppar)
                    ppar_p.mfc = 'lightgrey'
                    ppar_p.ms  = ppar.ms * 0.7
                    ppar_f     = copy.deepcopy(ppar)
                    ppar_f.linewidth  = 2
                    ppar_f.color      = 'lightgrey'
                    plot_ts_pred_curve(self.ax_pts, d_tsi_f, self, ppar_f)
                    self.ts_plot_func(self.ax_pts, d_tsi_p, self, ppar_p)

        # axis format
        self.ax_pts = _adjust_ts_axis(self.ax_pts, self)
        title_ts = _get_ts_title(yx[0], yx[1], self.coord)
        if self.mask[y, x] == 0:
            title_ts += ' (masked out)'
        if self.disp_title:
            self.ax_pts.set_title(title_ts, fontsize=self.font_size)
        if self.tick_right:
            self.ax_pts.yaxis.tick_right()
            self.ax_pts.yaxis.set_label_position("right")

        # legend
        if len(self.ts_data) > 1:
            self.ax_pts.legend()

        # Print to terminal
        vprint('\n---------------------------------------')
        vprint(title_ts)
        float_formatter = lambda x: [float('{:.2f}'.format(i)) for i in x]
        vprint(float_formatter(d_ts[0]))

        if not np.all(np.isnan(d_ts[0])):
            # stat info
            vprint('time-series range: [{:.2f}, {:.2f}] {}'.format(np.nanmin(d_ts[0]),
                                                                   np.nanmax(d_ts[0]),
                                                                   self.disp_unit))

            # estimate (print) slope
            estimate_slope(d_ts[0], self.yearList, ex_flag=self.ex_flag, disp_unit=self.disp_unit)

            # update figure
            self.fig_pts.canvas.draw()

            # print model parameters
            if self.show_model:
                for table_row in self.model_table: print(table_row)

            # write output to a txt file
            if self.savetxt is not False:
                print('write outputs to file: {}'.format(self.savetxt))
                self.write_file(title_ts, d_ts)
        return d_ts


    def update_plot_timeseries(self, event):
        """Event function to get y/x from button press"""
        if event.inaxes == self.ax_img:
            # get row/col number
            if self.fig_coord == 'geo':
                y, x = self.coord.geo2radar(event.ydata, event.xdata, print_msg=False)[0:2]
            else:
                y, x = int(event.ydata+0.5), int(event.xdata+0.5)

            # model fitting
            if self.show_model:
                self.m, self.m_std, self.ts_pred, self.ts_fit = self.fit2model(self.ts_data[0][:,y,x])
                self.model_table = self.model_param_table()

            # plot time-series displacement
            self.plot_point_timeseries((y, x))
        return


    def on_key_event(self, event):
        """Slide images with left/right key on keyboard"""
        if event.inaxes and event.inaxes.figure == self.fig_img:
            idx = None
            if event.key == 'left':
                idx = max(self.idx - 1, 0)
            elif event.key == 'right':
                idx = min(self.idx + 1, self.num_date - 1)

            if idx is not None and idx != self.idx:
                # update title
                disp_date = self.dates[idx].strftime('%Y-%m-%d')
                sub_title = 'N = {n}, Time = {t}'.format(n=idx, t=disp_date)
                self.ax_img.set_title(sub_title, fontsize=self.font_size)

                # read data
                data_img = np.array(self.ts_data[0][idx, :, :])
                data_img[self.mask == 0] = np.nan
                if self.wrap:
                    if self.disp_unit_img == 'radian':
                        data_img *= self.range2phase
                    data_img = ut.wrap(data_img, wrap_range=self.wrap_range)

                # update
                self.tslider.set_val(self.yearList[idx]) # update slider
                self.img.set_data(data_img)              # update image
                self.idx = idx
                self.fig_img.canvas.draw()
        return


    def fit2model(self, d_ts):
        """Fir a temporal model to the time series"""
        # convert datetime obj to str YYYYMMDD
        disp_dates = []
        for dt in self.dates:
            disp_dates.append(dt.strftime("%Y%m%d"))
        num_date = len(disp_dates)

        # initiate output
        m = np.zeros(self.num_param, dtype=np.float32)
        m_std = np.zeros(self.num_param, dtype=np.float32)

        ## only option 2 - least squares with uncertainty propagation
        print('\nestimate time functions via linalg.lstsq ...')
        G, m, e2 = ts2vel.estimate_time_func(model=self.model, date_list=disp_dates, dis_ts=d_ts)

        # reconstruction of predicted time series
        ts_pred = np.matmul(G,m)

        # reconstruct the fine resolution function
        disp_dates_fine = ptime.date2arange(disp_dates[0], disp_dates[-1])
        G_fine = timeseries.get_design_matrix4time_func(date_list=disp_dates_fine, model=self.model)
        ts_fit = np.matmul(G_fine,m)

        ## Compute the covariance matrix for model parameters: Gm = d
        # C_m_hat = (G.T * C_d^-1, * G)^-1  # linear propagation from the TS covariance matrix. (option 2.1)
        #         = sigma^2 * (G.T * G)^-1  # assuming obs errors are normally dist. in time.   (option 2.2a)
        # Based on the law of integrated expectation, we estimate the obs sigma^2 using
        # the OLS estimation residual e_hat_i = d_i - d_hat_i
        # option 2.2a - assume obs errors following normal dist. in time
        G_inv = linalg.inv(np.dot(G.T, G))
        m_var = e2.flatten() / (num_date - self.num_param)
        m_std = np.sqrt(np.dot(np.diag(G_inv).reshape(-1, 1), m_var))
        return m, m_std, ts_pred, ts_fit


    def model_param_table(self):
        """Make a table of the estimated model parameters"""
        model_table = ['Complex model param'+' '*15+'Value'+' '*8+'Std']

        str_fmt = '{:s} {:10.3f} {:10.3f}'
        str_lj  = 28

        # deformation model info
        poly_deg   = self.model['polynomial']
        num_period = len(self.model['periodic'])
        num_step   = len(self.model['step'])
        num_exp    = sum([len(val) for key, val in self.model['exp'].items()])

        # time func 1 - polynomial
        for i in range(1, poly_deg+1):
            # dataset name
            if i == 1:
                dsName = 'velocity'
            elif i == 2:
                dsName = 'acceleration'
            else:
                dsName = 'poly{}'.format(i)
            # write
            model_table.append(str_fmt.format(dsName.ljust(str_lj), self.m[i], self.m_std[i]))

        # time func 2 - periodic
        p0 = poly_deg + 1
        for i in range(num_period):
            # calculate the amplitude and phase of the periodic signal
            # following equation (9-10) in Minchew et al. (2017, JGR)
            coef_cos = self.m[p0 + 2*i]
            coef_sin = self.m[p0 + 2*i + 1]
            period_amp = np.sqrt(coef_cos**2 + coef_sin**2)

            # dataset name
            period = self.model['periodic'][i]
            dsNameSuffixes = ['Amplitude', 'Phase']
            if period == 1:
                dsNames = [f'annual{x}' for x in dsNameSuffixes]
            elif period == 0.5:
                dsNames = [f'semiAnnual{x}' for x in dsNameSuffixes]
            else:
                dsNames = [f'periodY{period}{x}' for x in dsNameSuffixes]
            # write
            for dsName in dsNames:
                model_table.append(str_fmt.format(dsName.ljust(str_lj), period_amp, np.nan))

        # time func 3 - step
        p0 = (poly_deg + 1) + (2 * num_period)
        for i in range(num_step):
            # dataset name
            dsName = 'step{}'.format(self.model['step'][i])
            # write
            model_table.append(str_fmt.format(dsName.ljust(str_lj), self.m[p0+i], self.m_std[p0+i]))

        # time func 4 - exponential
        p0 = (poly_deg + 1) + (2 * num_period) + (num_step)
        i = 0
        for exp_onset in self.model['exp'].keys():
            for exp_tau in self.model['exp'][exp_onset]:
                # dataset name
                dsName = 'exp{}Tau{}'.format(exp_onset, exp_tau)
                # write
                model_table.append(str_fmt.format(dsName.ljust(str_lj), self.m[p0+i], self.m_std[p0+i]))
                i += 1

        # time func 5 - logarithmic
        p0 = (poly_deg + 1) + (2 * num_period) + (num_step) + (num_exp)
        i = 0
        for log_onset in self.model['log'].keys():
            for log_tau in self.model['log'][log_onset]:
                # dataset name
                dsName = 'log{}Tau{}'.format(log_onset, log_tau)
                # write
                model_table.append(str_fmt.format(dsName.ljust(str_lj), self.m[p0+i], self.m_std[p0+i]))
                i += 1
        return model_table


    def write_file(self, title_ts, d_ts):
        # get dates string
        dates = []
        for dt in self.dates:
            dates.append(dt.strftime("%Y%m%d"))

        # get a table of ts for output
        float_formatter = lambda x: [float('{:.3f}'.format(i)) for i in x]
        obs = list(map(str, float_formatter(d_ts[0])))
        if self.show_model:
            preds    = list(map(str, float_formatter(self.ts_pred)))
            ts_table = list(zip(dates, obs, preds))
        else:
            ts_table = list(zip(dates, obs))

        # write to file
        with open(self.savetxt,'a') as f:
            f.write('# '+title_ts)
            if not self.show_model:
                f.write('\n')
                f.write('# Date\t\tObs\n')
                for line in ts_table:
                    f.write('{}\t{}\n'.format(line[0],line[1]))
            else:
                f.write('\n')
                for table_row in self.model_table:
                    f.write('# '+table_row+'\n')
                f.write('# Date\t\tObs\tPred\n')
                for line in ts_table:
                    f.write('{}\t{}\t{}\n'.format(line[0],line[1],line[2]))
            f.write('\n')
        return

###########################################################################################
def main(iargs=None):
    obj = timeseriesViewer(iargs=iargs)
    obj.configure()
    obj.plot()
    #obj.fig_img.canvas.mpl_disconnect(obj.cid)
    return


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
