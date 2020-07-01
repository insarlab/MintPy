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
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.widgets import Slider
from mintpy.objects import timeseries, giantTimeseries, HDFEOS
from mintpy.utils import readfile, ptime, plot as pp, utils as ut
from mintpy.multilook import multilook_data
from mintpy import subset, view


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
    parser.add_argument('timeseries_file', nargs='+',
                        help='time-series file to display\n'
                             'i.e.: timeseries_ERA5_ramp_demErr.h5 (MintPy)\n'
                             '      LS-PARAMS.h5 (GIAnT)\n'
                             '      S1_IW12_128_0593_0597_20141213_20180619.he5 (HDF-EOS5)')
    parser.add_argument('--label', dest='file_label', nargs='*',
                        help='labels to display for multiple input files')
    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float,
                        help='Y limits for point plotting.')
    parser.add_argument('--tick-right', dest='tick_right', action='store_true',
                        help='set tick and tick label to the right')

    parser.add_argument('-l','--lookup', dest='lookup_file', type=str,
                        help='lookup table file')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2,
                       help='initial pixel to plot in lat/lon coord')

    pixel.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=6.0,
                       help='Point marker size. Default: 6')
    pixel.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0,
                       help='Edge width. Default: 1.0')

    parser.add_argument('-n', dest='idx', metavar='NUM', type=int,
                        help='Epoch/slice number to display.')
    parser.add_argument('--error', dest='error_file',
                        help='txt file with error for each date.')

    parser.add_argument('--start-date', dest='start_date', type=str,
                        help='start date of displacement to display')
    parser.add_argument('--end-date', dest='end_date', type=str,
                        help='end date of displacement to display')
    parser.add_argument('--exclude', '--ex', dest='ex_date_list',
                        nargs='*', default=['exclude_date.txt'],
                        help='Exclude date shown as gray.')
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                        help='Set displacement at first acquisition to zero.')
    parser.add_argument('--off','--offset', dest='offset', type=float,
                        help='Offset for each timeseries file.')

    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')

    parser = pp.add_data_disp_argument(parser)
    parser = pp.add_dem_argument(parser)
    parser = pp.add_figure_argument(parser)
    parser = pp.add_gps_argument(parser)
    parser = pp.add_mask_argument(parser)
    parser = pp.add_map_argument(parser)
    parser = pp.add_reference_argument(parser)
    parser = pp.add_save_argument(parser)
    parser = pp.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if '--gps-comp' in iargs:
        msg = '--gps-comp is not supported for {}'.format(os.path.basename(__file__))
        raise NotImplementedError(msg)

    if inps.file_label:
        if len(inps.file_label) != len(inps.timeseries_file):
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

    # verbose print using --noverbose option
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


###########################################################################################
def read_init_info(inps):
    # Time Series Info
    ts_file0 = inps.timeseries_file[0]
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
        inps.file_label = [str(i) for i in list(range(len(inps.timeseries_file)))]

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

    # initial display index
    #if atr['REF_DATE'] in inps.date_list:
    #    inps.ref_idx = inps.date_list.index(atr['REF_DATE'])
    #else:
    #    inps.ref_idx = 0

    if inps.ref_date:
        inps.ref_idx = inps.date_list.index(inps.ref_date)
    else:
        inps.ref_idx = 3

    if not inps.idx:
        if inps.ref_idx < inps.num_date / 2.:
            inps.idx = inps.num_date - 3
        else:
            inps.idx = 3

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
        if inps.ref_lalo[1] > 180.:
            inps.ref_lalo[1] -= 360.
        inps.ref_yx = inps.coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1], print_msg=False)[0:2]
    if not inps.ref_yx:
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]

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
    for fname in inps.timeseries_file:
        vprint('reading timeseries from file {} ...'.format(fname))
        data, atr = readfile.read(fname, datasetName=inps.date_list, box=inps.pix_box)
        try:
            ref_phase = data[:, inps.ref_yx[0]-inps.pix_box[1], inps.ref_yx[1]-inps.pix_box[0]]
            data -= np.tile(ref_phase.reshape(-1, 1, 1), (1, data.shape[-2], data.shape[-1]))
            vprint('reference to pixel: {}'.format(inps.ref_yx))
        except:
            pass
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
    msk = pp.read_mask(inps.timeseries_file[0],
                       mask_file=inps.mask_file,
                       datasetName='displacement',
                       box=inps.pix_box,
                       print_msg=inps.print_msg)[0]
    mask[msk == 0.] = False
    del msk

    ts_stack = np.sum(ts_data[0], axis=0)
    mask[ts_stack == 0.] = False
    mask[np.isnan(ts_stack)] = False
    del ts_stack

    #do not mask the reference point
    try:
        mask[inps.ref_yx[0]-inps.pix_box[1],
             inps.ref_yx[1]-inps.pix_box[0]] = True
    except:
        pass

    #vprint('masking data')
    #ts_mask = np.tile(mask, (inps.num_date, 1, 1))
    #for i in range(len(ts_data)):
    #    ts_data[i][ts_mask == 0] = np.nan
    #    try:
    #        ts_data[i][:, inps.ref_yx[0], inps.ref_yx[1]] = 0.   # keep value on reference pixel
    #    except:
    #        pass
    #del ts_mask

    # default vlim
    inps.dlim = [np.nanmin(ts_data[0]), np.nanmax(ts_data[0])]
    ts_data_mli = multilook_data(np.squeeze(ts_data[0]), 10, 10)
    if not inps.vlim:
        inps.vlim = [np.nanmin(ts_data_mli[inps.ex_flag != 0]),
                     np.nanmax(ts_data_mli[inps.ex_flag != 0])]
    vprint('data    range: {} {}'.format(inps.dlim, inps.disp_unit))
    vprint('display range: {} {}'.format(inps.vlim, inps.disp_unit))

    # default ylim
    num_file = len(inps.timeseries_file)
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
        ax.scatter(inps.ex_dates, ex_d_ts, s=ppar.ms**2, color='gray')

    # Plot kept dates
    ax.scatter(dates, d_ts, s=ppar.ms**2, label=ppar.label, color=ppar.mfc)
    return ax


def _adjust_ts_axis(ax, inps):
    ax.tick_params(which='both', direction='in', labelsize=inps.font_size, bottom=True, top=True, left=True, right=True)
    ax = pp.auto_adjust_xaxis_date(ax, inps.yearList, fontsize=inps.font_size)[0]
    ax.set_xlabel('Time [years]', fontsize=inps.font_size)
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
        inps.outfile_base = 'y{}_x{}'.format(yx[0], yx[1])

    # get aux info
    vel, std = estimate_slope(d_ts[0], inps.yearList,
                              ex_flag=inps.ex_flag,
                              disp_unit=inps.disp_unit)

    # TXT - point time-series
    outName = '{}_ts.txt'.format(inps.outfile_base)
    header_info = 'timeseries_file={}\n'.format(inps.timeseries_file)
    header_info += '{}\n'.format(_get_ts_title(yx[0], yx[1], inps.coord))
    header_info += 'reference pixel: y={}, x={}\n'.format(inps.ref_yx[0], inps.ref_yx[1])
    header_info += 'reference date: {}\n'.format(inps.date_list[inps.ref_idx])
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

        # Figure 1 - Cumulative Displacement Map
        self.fig_img = plt.figure(self.figname_img, figsize=self.figsize_img)

        # Figure 1 - Axes 1 - Displacement Map
        self.ax_img = self.fig_img.add_axes([0.125, 0.25, 0.75, 0.65])
        img_data = np.array(self.ts_data[0][self.idx, :, :])  ####################
        img_data[self.mask == 0] = np.nan
        self.plot_init_image(img_data)

        # Figure 1 - Axes 2 - Time Slider
        self.ax_tslider = self.fig_img.add_axes([0.2, 0.1, 0.6, 0.07])
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


    def plot_init_time_slider(self, init_idx=-1, ref_idx=0):
        val_step = np.min(np.diff(self.yearList))
        val_min = self.yearList[0]
        val_max = self.yearList[-1]

        self.tslider = Slider(self.ax_tslider, label='Years',
                              valinit=self.yearList[init_idx],
                              valmin=val_min,
                              valmax=val_max,
                              valstep=val_step)

        bar_width = val_step / 4.
        datex = np.array(self.yearList) - bar_width / 2.
        self.tslider.ax.bar(datex, np.ones(len(datex)), bar_width, facecolor='black', ecolor=None)
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
        return self.tslider


    def update_time_slider(self, val):
        """Update Displacement Map using Slider"""
        idx = np.argmin(np.abs(np.array(self.yearList) - self.tslider.val))
        # update title
        disp_date = self.dates[idx].strftime('%Y-%m-%d')
        self.ax_img.set_title('N = {n}, Time = {t}'.format(n=idx, t=disp_date),
                              fontsize=self.font_size)
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
            if self.zero_first:
                d_tsi -= d_tsi[self.zero_idx]
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

            # plot
            if not np.all(np.isnan(d_tsi)):
                self.ax_pts = self.ts_plot_func(self.ax_pts, d_tsi, self, ppar)

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
            vprint('displacement range: [{:.2f}, {:.2f}] {}'.format(np.nanmin(d_ts[0]),
                                                                    np.nanmax(d_ts[0]),
                                                                    self.disp_unit))

            # estimate (print) slope
            estimate_slope(d_ts[0], self.yearList, ex_flag=self.ex_flag, disp_unit=self.disp_unit)

            # update figure
            self.fig_pts.canvas.draw()
        return d_ts


    def update_plot_timeseries(self, event):
        """Event function to get y/x from button press"""
        if event.inaxes == self.ax_img:
            # get row/col number
            if self.fig_coord == 'geo':
                y, x = self.coord.geo2radar(event.ydata, event.xdata, print_msg=False)[0:2]
            else:
                y, x = int(event.ydata+0.5), int(event.xdata+0.5)

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
                self.ax_img.set_title('N = {n}, Time = {t}'.format(n=idx, t=disp_date),
                                      fontsize=self.font_size)

                # read data
                data_img = np.array(self.ts_data[0][idx, :, :])
                data_img[self.mask == 0] = np.nan
                if self.wrap:
                    if self.disp_unit_img == 'radian':
                        data_img *= self.range2phase
                    data_img = ut.wrap(data_img, wrap_range=self.wrap_range)

                # update
                self.img.set_data(data_img)              # update image
                self.tslider.set_val(self.yearList[idx]) # update slider
                self.idx = idx
                self.fig_img.canvas.draw()
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
