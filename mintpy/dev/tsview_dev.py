#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, Zhang Yunjun, 2017                #
############################################################


import os
import sys
import argparse
from datetime import datetime as dt
import h5py
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mintpy.objects import timeseries
from mintpy.utils import readfile, ptime, utils as ut, plot as pp
from mintpy.mask import mask_matrix


############# Global Variables ################
tims, inps, img, mask, d_v, d_ts = None, None, None, None, None, None
ax_v, fig_ts, fig_v, ax_ts, tslider, second_plot_axis, new_axes = None, None, None, None, None, None, None
h5, k, dateList, atr, date_num = None, None, None, None, None
lat, lon, ullat, ullon, lat_step, lon_step = None, None, None, None, None, None
width, length = None, None

plot_figure, p1_scatter, p2_scatter, scatts = None, None, None, None
p1_scatter_point, p2_scatter_point = None, None
p1_x, p1_y, p2_x, p2_y = None, None, None, None
annot = None
second_plot_axis_visible = False


###########################################################################################
EXAMPLE='''example:
  tsview.py timeseries.h5 --ylim -10 10
  tsview.py timeseries_demErr_ramp.h5 -n 5 -m maskTempCoh.h5
  tsview.py timeseries_demErr_ramp.h5 --yx 300 400 --nodisplay --zero-first
  tsview.py geo_timeseries_demErr_ramp.h5 --lalo 33.250 131.665 --nodisplay
'''

'''
    Creates command line argument parser and sets inps default values
'''
def create_parser():
    parser = argparse.ArgumentParser(description='Interactive Time-series Viewer',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('timeseries_file', help='time series file to display')
    parser.add_argument('-n', dest='epoch_num', metavar='NUM', type=int, default='-2',\
                        help='Epoch/slice number to display, default: the 2nd last.')
    parser.add_argument('-m','--mask', dest='mask_file',\
                        help='mask to use. Default: geo_maskTempCoh.h5 for geocoded file and maskTempCoh.h5 for radar file')
    parser.add_argument('--error', dest='error_file', help='txt file with error for each date.')
    parser.add_argument('--dem', dest='dem_file', help='DEM file for background shaed relief')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, metavar=('Y','X'), nargs=2,\
                       help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, metavar=('LAT','LON'), nargs=2,\
                       help='initial pixel to plot in lat/lon coord')
    pixel.add_argument('--ref-yx', dest='ref_yx', type=int, metavar=('Y','X'), nargs=2,\
                       help='change reference pixel to input location')
    pixel.add_argument('--ref-lalo', dest='ref_lalo', type=float, metavar=('LAT','LON'), nargs=2,\
                       help='change reference pixel to input location')

    output = parser.add_argument_group('Output Setting')
    output.add_argument('-o','--output', dest='fig_base', help='Figure base name for output files')
    output.add_argument('--save', action='store_true', dest='save_fig',\
                        help='save data and plot to files')
    output.add_argument('--nodisplay', action='store_false', dest='disp_fig',\
                        help='save data and plot to files and do not display figures\n')
    output.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=150,\
                        help='DPI - dot per inch - for display/write')

    disp = parser.add_argument_group('Display Setting')
    disp.add_argument('--figsize', dest='fig_size', metavar=('WID','LEN'), type=float, nargs=2, default=[12.0,5.0],\
                      help='Figure size in inches - width and length. Default: 12.0 5.0\n'+\
                           'i.e. 3.5 2 for ppt; ')
    disp.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN','YMAX'), type=float,\
                      help='Y limits for point plotting.')
    disp.add_argument('--ylim-mat', dest='ylim_mat', nargs=2, metavar=('YMIN','YMAX'), type=float,\
                      help='Display limits for matrix plotting.')
    disp.add_argument('--ref-date', dest='ref_date', help='Change reference date for display')
    disp.add_argument('--exclude','--ex', dest='ex_date_list', nargs='*', help='Exclude date shown as gray.')
    disp.add_argument('--zf','--zero-first', dest='zero_first', action='store_true',\
                      help='Set displacement at first acquisition to zero.')
    disp.add_argument('-u', dest='disp_unit', metavar='UNIT', default='cm',\
                      help='unit for display. Default: cm')
    disp.add_argument('-c','--colormap', dest='colormap', default='jet',\
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r etc.\n'
                           'Support colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html')
    disp.add_argument('-s','--fontsize', dest='font_size', type=int, default=10, help='Font size for display')
    disp.add_argument('--notitle', dest='disp_title', action='store_false', help='Do not display title in TS plot.')
    disp.add_argument('--no-flip', dest='auto_flip', action='store_false',\
                      help='Turn off auto flip based on orbit direction.\n'+\
                           'Default: flip left-right for descending data in radar coord\n'+\
                           '         flip up-down    for ascending  data in radar coord\n'+\
                           '         no flip for data in geo coord')
    disp.add_argument('--ms','--markersize', dest='marker_size', type=float, default=8.0,\
                      help='Point marker size. Default: 8.0')
    #disp.add_argument('--mc','--markercolor', dest='marker_color', default='crimson',\
    #                  help='Point marker color. Default: crimson')
    disp.add_argument('--ew','--edgewidth', dest='edge_width', type=float, default=1.0,\
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
def read_timeseries_yx(timeseries_file, y, x, ref_yx=None):
    '''Read time-series displacement on point (y,x) from timeseries_file
    Inputs:
        timeseries_file : string, name/path of timeseries hdf5 file
        y/x : int, row/column number of point of interest
    Output:
        dis_ts : list of float, displacement time-series of point of interest
    '''
    atr = readfile.read_attribute(timeseries_file)
    k = atr['FILE_TYPE']
    dis_ts = []

    if k in ['GIANT_TS']:
        h5 = h5py.File(timeseries_file, 'r')
        date_list = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in h5['dates'][:].tolist()]
        dname = [i for i in ['rawts','recons'] if i in list(h5.keys())][0]
        dis_ts = h5[dname][:,y,x]
        if ref_yx is not None:
            dis_ts = h5[dname][:,ref_yx[0],ref_yx[1]]
        h5.close()
    else:
        box = [x, y, x+1, y+1]
        dis_ts = timeseries(timeseries_file).read(box=box, print_msg=False)
        #date_list = list(h5[k].keys())
        #for date in date_list:
        #    dis = h5[k].get(date)[y,x]
        #    if inps.ref_yx:
        #        dis -= h5[k].get(date)[ref_yx[0], ref_yx[1]]
        #    dis_ts.append(dis)
        #dis_ts = np.array(dis_ts)

    return dis_ts


def read_timeseries_lalo(timeseries_file, lat, lon):
    '''Read time-series displacement on point (y,x) from timeseries_file
    Inputs:
        timeseries_file : string, name/path of timeseries hdf5 file
        lat/lon : float, latitude/longitude of point of interest
    Output:
        dis_ts : list of float, displacement time-series of point of interest
    '''

    atr = readfile.read_attribute(timeseries_file)
    if 'X_FIRST' not in list(atr.keys()):
        print('ERROR: input file is not geocoded')
        return None

    lat0 = float(atr['Y_FIRST'])
    lat_step = float(atr['Y_STEP'])
    lon0 = float(atr['X_FIRST'])
    lon_step = float(atr['X_STEP'])
    y = int(np.rint((lat-lat0)/lat_step))
    x = int(np.rint((lon-lon0)/lon_step))
    dis_ts = read_timeseries_yx(timeseries_file, y, x)
    return dis_ts


################### HELPER FUNCTIONS ##########################

###################### EXTRANEOUS HELPER FUNCTIONS ######################

'''
    Handles displaying of the plot and figure
'''
def display_figure():
    global inps

    if inps.disp_fig:
        plt.show()

'''
    Plots initial point on map and sets timeseries data points on scatter plot
    to appropriate values
'''
def plot_data_from_inital_point():
    global ax_ts, inps, tims, d_ts

    if inps.yx:
        d_ts = update_timeseries(inps.yx[0], inps.yx[1], 1)
    else:
        d_ts = np.zeros(len(tims))
        ax_ts, scatter = plot_timeseries_scatter(ax_ts, d_ts, inps)

'''
    Reads list of errors from error file
'''
def read_error_list():
    global inps, date_num

    inps.error_ts = None
    if inps.error_file:
        error_file_content = np.loadtxt(inps.error_file, dtype=str)
        inps.error_ts = error_file_content[:, 1].astype(np.float) * inps.unit_fac
        if inps.ex_date_list:
            e_ts = inps.error_ts[:]
            inps.ex_error_ts = np.array([e_ts[i] for i in inps.ex_idx_list])
            inps.error_ts = np.array([e_ts[i] for i in range(date_num) if i not in inps.ex_idx_list])

'''
    Saves figure and data to output file
'''
def save_output():
    global inps, lat, lon, ullat, lat_step, ullon, lon_step, atr, fig_ts, dateList

    if inps.save_fig and inps.yx:
        print(('save info for pixel ' + str(inps.yx)))
        if not inps.fig_base:
            inps.fig_base = 'y%d_x%d' % (inps.yx[0], inps.yx[1])

        # TXT - point time series
        outName = inps.fig_base + '_ts.txt'
        header_info = 'timeseries_file=' + inps.timeseries_file
        header_info += '\ny=%d, x=%d' % (inps.yx[0], inps.yx[1])

        try:
            lat = ullat + inps.yx[0] * lat_step
            lon = ullon + inps.yx[1] * lon_step
            header_info += '\nlat=%.6f, lon=%.6f' % (lat, lon)
        except:
            pass

        if inps.ref_yx:
            header_info += '\nreference pixel: y=%d, x=%d' % (inps.ref_yx[0], inps.ref_yx[1])
        else:
            header_info += '\nreference pixel: y=%s, x=%s' % (atr['REF_Y'], atr['REF_X'])

        header_info += '\nunit=m/yr'
        np.savetxt(outName, list(zip(np.array(dateList), np.array(d_ts) / inps.unit_fac)), fmt='%s', \
                   delimiter='    ', header=header_info)
        print(('save time series displacement in meter to ' + outName))


        # Figure - point time series
        #outName = inps.fig_base + '_ts.pdf'
        #fig_ts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        #print(('save time series plot to ' + outName))

        # Figure - map
        outName = '{}_ts.png'.format(inps.fig_base)   # , dateList[inps.epoch_num])
        fig_v.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        print(('save plot to ' + outName))



################### PLOT SETUP HELPER FUCNTIONS ###################

'''
    Reads basic information about timeseries file being viewed
'''
def read_timeseries_info():
    global atr, k, h5, dateList, tims, date_num, inps

    atr = readfile.read_attribute(inps.timeseries_file)
    k = atr['FILE_TYPE']
    print(('input file is '+k+': '+inps.timeseries_file))
    if not k in ['timeseries','GIANT_TS']:
        raise ValueError('Only timeseries file is supported!')

    h5 = h5py.File(inps.timeseries_file,'r')
    if k in ['GIANT_TS']:
        dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in h5['dates'][:].tolist()]
    else:
        dateList = timeseries(inps.timeseries_file).get_date_list()
    date_num = len(dateList)
    inps.dates, tims = ptime.date_list2vector(dateList)


'''
    Sets dates to be excluded from timeseries data points
'''
def exclude_dates():
    global inps, dateList

    if inps.ex_date_list:
        input_ex_date = list(inps.ex_date_list)
        inps.ex_date_list = []

        if input_ex_date:
            for ex_date in input_ex_date:

                if os.path.isfile(ex_date):
                    ex_date = ptime.read_date_list(ex_date)
                else:
                    ex_date = [ptime.yyyymmdd(ex_date)]

                inps.ex_date_list += list(set(ex_date) - set(inps.ex_date_list))

            # delete dates not existed in input file
            inps.ex_date_list = sorted(list(set(inps.ex_date_list).intersection(dateList)))
            inps.ex_dates = ptime.date_list2vector(inps.ex_date_list)[0]
            inps.ex_idx_list = sorted([dateList.index(i) for i in inps.ex_date_list])
            print(('exclude date:' + str(inps.ex_date_list)))


'''
    Sets the 'zero' value for the plot
'''
def set_zero_displacement():
    global inps, date_num

    if inps.zero_first:
        if inps.ex_date_list:
            inps.zero_idx = min(list(set(range(date_num)) - set(inps.ex_idx_list)))
        else:
            inps.zero_idx = 0


'''
    Computed dimensions of file (width, length)
'''
def compute_file_size():
    global atr, width, length

    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    print(('data size in [y0,y1,x0,x1]: [%d, %d, %d, %d]' % (0, length, 0, width)))


'''
    Computes parameters needed for setting lat/lon values
'''
def compute_lat_lon_params():
    global ullon, ullat, lon_step, lat_step, atr, width, length

    try:
        ullon = float(atr['X_FIRST'])
        ullat = float(atr['Y_FIRST'])
        lon_step = float(atr['X_STEP'])
        lat_step = float(atr['Y_STEP'])
        lrlon = ullon + width * lon_step
        lrlat = ullat + length * lat_step
        print(('data size in [lat0,lat1,lon0,lon1]: [%.4f, %.4f, %.4f, %.4f]' % (lrlat, ullat, ullon, lrlon)))
        return ullon, ullat, lon_step, lat_step, lrlon, lrlat
    except:
        pass


'''
    Sets coordinates of initial pixel
'''
def set_inital_pixel_coords():
    global inps, atr

    if inps.lalo and 'Y_FIRST' in list(atr.keys()):
        y, x = set_yx_coords(inps.lalo[0], inps.lalo[1])
        inps.yx = [y, x]
    if inps.ref_lalo and 'Y_FIRST' in list(atr.keys()):
        y, x = set_yx_coords(inps.ref_lalo[0], inps.ref_lalo[1])
        inps.ref_yx = [y, x]

    # Display Unit
    if inps.disp_unit == 'cm': inps.unit_fac = 100.0
    elif inps.disp_unit == 'm': inps.unit_fac = 1.0
    elif inps.disp_unit == 'dm': inps.unit_fac = 10.0
    elif inps.disp_unit == 'mm': inps.unit_fac = 1000.0
    elif inps.disp_unit == 'km': inps.unit_fac = 0.001
    else: raise ValueError('Un-recognized unit: '+inps.disp_unit)
    if k in ['GIANT_TS']:
        print('data    unit: mm')
        inps.unit_fac *= 0.001
    else:
        print('data    unit: m')
    print(('display unit: '+inps.disp_unit))


'''
    Sets x/y coordinates from lalo valyes
'''
def set_yx_coords(y_input, x_input):
    global ullat, ullon, lat_step, lon_step

    y = int((y_input - ullat) / lat_step + 0.5)
    x = int((x_input - ullon) / lon_step + 0.5)

    return y, x


'''
    Sets multiplier for data based on display unit
'''
def set_unit_fraction():
    global inps

    if inps.disp_unit == 'cm':
        inps.unit_fac = 100.0
    elif inps.disp_unit == 'm':
        inps.unit_fac = 1.0
    elif inps.disp_unit == 'dm':
        inps.unit_fac = 10.0
    elif inps.disp_unit == 'mm':
        inps.unit_fac = 1000.0
    elif inps.disp_unit == 'km':
        inps.unit_fac = 0.001
    else:
        raise ValueError('Un-recognized unit: ' + inps.disp_unit)
    print('data    unit: m')
    print(('display unit: ' + inps.disp_unit))


'''
    Flips displat map ud/lr
'''
def flip_map():
    global inps, atr

    if inps.auto_flip:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(atr)
    else:
        inps.flip_ud = False
        inps.flip_lr = False


'''
    Sets mask for map
'''
def set_mask():
    global mask, inps, atr

    if not inps.mask_file:
        if os.path.basename(inps.timeseries_file).startswith('geo_'):
            file_list = ['geo_maskTempCoh.h5']
        else:
            file_list = ['maskTempCoh.h5', 'maskConnComp.h5']

        try:
            inps.mask_file = ut.get_file_list(file_list)[0]
        except:
            inps.mask_file = None

    try:
        mask = readfile.read(inps.mask_file, datasetName='mask')[0]
        mask[mask!=0] = 1
        print(('load mask from file: '+inps.mask_file))
    except:
        mask = None
        print('No mask used.')


'''
    Sets the initial map plot
'''
def set_initial_map():
    global d_v, h5, k, dateList, inps, data_lim

    d_v = h5['timeseries'][inps.epoch_num][:] * inps.unit_fac
    # Initial Map
    print(str(dateList))
    d_v = readfile.read(inps.timeseries_file, datasetName=dateList[inps.epoch_num])[0] * inps.unit_fac
    #d_v = h5[k].get(dateList[inps.epoch_num])[:]*inps.unit_fac
    if inps.ref_date:
        inps.ref_d_v = readfile.read(inps.timeseries_file, datasetName=inps.ref_date)[0]*inps.unit_fac
        d_v -= inps.ref_d_v

    if mask is not None:
        d_v = mask_matrix(d_v, mask)

    if inps.ref_yx:
        d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]

    data_lim = [np.nanmin(d_v), np.nanmax(d_v)]

    if not inps.ylim_mat:
        inps.ylim_mat = data_lim

    print(('Initial data range: '+str(data_lim)))
    print(('Display data range: '+str(inps.ylim_mat)))

    print(('Initial data range: ' + str(data_lim)))
    print(('Display data range: ' + str(inps.ylim)))


def setup_plot():
    # Time Series Info
    read_timeseries_info()
    # Read exclude dates
    exclude_dates()
    # Zero displacement for 1st acquisition
    set_zero_displacement()
    # File Size
    compute_file_size()
    # Latitude Longitude Parameters
    compute_lat_lon_params()
    # Initial Pixel Coordinates
    set_inital_pixel_coords()
    # Display Unit
    set_unit_fraction()
    # Flip up-down / left-right
    flip_map()
    # Mask file
    set_mask()
    # Initial Map
    set_initial_map()



################# PLOT CONFIGURATION HELPER METHODS #######################

'''
    Sets DEM topography file for the map
'''
def set_dem_file():
    global ax_v, inps, img

    if inps.dem_file:
        dem = readfile.read(inps.dem_file, datasetName='height')[0]
        ax_v = pp.plot_dem_yx(ax_v, dem)

    img = ax_v.imshow(d_v, cmap=inps.colormap, clim=inps.ylim_mat, interpolation='nearest')


'''
    Sets reference pixel on map
'''
def set_map_reference_pixel():
    global d_v, inps, ax_v, atr

    if inps.ref_yx:
        d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]
        ax_v.plot(inps.ref_yx[1], inps.ref_yx[0], 'ks', ms=6)
    else:
        try:
            ax_v.plot(int(atr['REF_X']), int(atr['REF_Y']), 'ks', ms=6)
        except:
            pass


'''
    Sets axis limits, labels, and titles
'''
def set_plot_axis_params():
    global inps, d_v, ax_v, atr

    if inps.yx:
        ax_v.plot(inps.yx[1], inps.yx[0], 'ro', markeredgecolor='black')

    ax_v.set_xlim(0, np.shape(d_v)[1])
    ax_v.set_ylim(np.shape(d_v)[0], 0)
    ax_v.format_coord = format_coord

    # Title and Axis Label
    ax_v.set_title('N = %d, Time = %s' % (inps.epoch_num, inps.dates[inps.epoch_num].strftime('%Y-%m-%d')))
    ax_v.yaxis.set_label_position("right")
    ax_v.yaxis.tick_right()

    #if not 'Y_FIRST' in list(atr.keys()):
    #    ax_v.set_xlabel('Range')
    #    ax_v.set_ylabel('Azimuth')


'''
    Flips axis lr/ud
'''
def flip_axis():
    global inps, ax_v

    if inps.flip_lr:
        ax_v.invert_xaxis()
        print('flip map left and right')
    if inps.flip_ud:
        ax_v.invert_yaxis()
        print('flip map up and down')


'''
    Creates colorbar for figure
'''
def make_color_bar():
    global fig_v, ax_v, img, inps
    #divider = make_axes_locatable(ax_v)
    #cax = divider.append_axes("left", "3%", pad="3%")
    cax = fig_v.add_axes([0.03, 0.4, 0.01, 0.4])
    cbar = fig_v.colorbar(img, cax=cax, orientation='vertical')
    #cbar = fig_v.colorbar(img, orientation='vertical')
    cbar.set_label('Displacement [%s]' % inps.disp_unit)
    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.yaxis.tick_right()

'''
    Creates timeseries slider for figure
'''
def make_time_slider():
    global tslider, fig_v, tims, inps

    ax_time = fig_v.add_axes([0.03, 0.10, 0.30, 0.07], facecolor='lightgoldenrodyellow', yticks=[])
    tslider = Slider(ax_time, '', tims[0], tims[-1], valinit=tims[inps.epoch_num])
    tslider.ax.bar(tims, np.ones(len(tims)), facecolor='black', width=0.01, ecolor=None)

    # xaxis tick format
    if np.floor(tims[-1]) == np.floor(tims[0]):
        digit = 10.
    else:
        digit = 1.
    tslider.ax.set_xticks(np.round(np.linspace(tims[0], tims[-1], num=5) * digit) / digit)
    tslider.ax.xaxis.set_minor_locator(MultipleLocator(1./12.))

    tslider.on_changed(time_slider_update)


def configure_plot():
    # DEM File
    set_dem_file()
    # Reference Pixel
    set_map_reference_pixel()
    # Initial Pixel
    set_plot_axis_params()
    # Flip Axis
    flip_axis()
    # Construct Color Bar
    make_color_bar()
    # Construct Time Slider
    make_time_slider()


################### PLOTTING HELPER FUNCTIONS #######################
def time_slider_update(val):
    '''Update Displacement Map using Slider'''
    global tims, tslider, ax_v, d_v, inps, img, fig_v, h5, k, dateList
    timein = tslider.val
    idx_nearest = np.argmin(np.abs(np.array(tims) - timein))
    ax_v.set_title('N = %d, Time = %s' % (idx_nearest, inps.dates[idx_nearest].strftime('%Y-%m-%d')))
    d_v = h5[k][idx_nearest][:] * inps.unit_fac
    if inps.ref_date:
        d_v -= inps.ref_d_v
    if mask is not None:
        d_v = mask_matrix(d_v, mask)
    if inps.ref_yx:
        d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]
    img.set_data(d_v)
    fig_v.canvas.draw()


def format_coord(x, y):
    '''Formats x, y coordinates into useful output string (used for creating plot titles)'''
    global width, length, ullat, lat_step, ullon, lon_step, d_v, lat, lon

    col = int(x + 0.5)
    row = int(y + 0.5)
    if 0 <= col < width and 0 <= row < length:
        z = d_v[row, col]
        try:
            lon = ullon + x * lon_step
            lat = ullat + y * lat_step
            return 'x=%.0f, y=%.0f, value=%.4f, lon=%.4f, lat=%.4f' % (x, y, z, lon, lat)
        except:
            return 'x=%.0f, y=%.0f, value=%.4f' % (x, y, z)


def plot_timeseries_errorbar(ax, dis_ts, inps):
    '''Plots errorbars for timeseries data'''
    global date_num
    dates = list(inps.dates)
    d_ts = dis_ts[:]
    if inps.ex_date_list:
        # Update displacement time-series
        dates = sorted(list(set(inps.dates) - set(inps.ex_dates)))
        ex_d_ts = np.array([dis_ts[i] for i in inps.ex_idx_list])
        d_ts = np.array([dis_ts[i] for i in range(date_num) if i not in inps.ex_idx_list])
        # Plot excluded dates
        (_, caps, _) = ax.errorbar(inps.ex_dates, ex_d_ts, yerr=inps.ex_error_ts, fmt='-o', color='gray', \
                                   ms=inps.marker_size, lw=0, alpha=1, mfc='gray', \
                                   elinewidth=inps.edge_width, ecolor='black', capsize=inps.marker_size * 0.5)
        for cap in caps:  cap.set_markeredgewidth(inps.edge_width)
    # Plot kept dates
    (_, caps, _) = ax.errorbar(dates, d_ts, yerr=inps.error_ts, fmt='-o', \
                               ms=inps.marker_size, lw=0, alpha=1, \
                               elinewidth=inps.edge_width, ecolor='black', capsize=inps.marker_size * 0.5)
    for cap in caps:  cap.set_markeredgewidth(inps.edge_width)
    return ax


def plot_timeseries_scatter(ax, dis_ts, inps, plot_num=1):
    '''Plots scatter points on provioded axis
        Inputs:
            ax      : mpl axis, matplotlib axis object on which to plot scattaer data
            dis_ts  : [float], data points for timeseries
            inps    : [Object], plot settings
            plot_num: int, plot delimiter to determine which scatter plot (1/2) is being updated
        Output:
            ax      : mpl axis, matplotlib axis object on which to plot scattaer data
            scatter : mpl scatter, matplotlib scatter object
    '''
    global date_num

    dates = list(inps.dates)
    d_ts = dis_ts[:]
    if inps.ex_date_list:
        # Update displacement time-series
        dates = sorted(list(set(inps.dates) - set(inps.ex_dates)))
        ex_d_ts = np.array([dis_ts[i] for i in inps.ex_idx_list])
        d_ts = np.array([dis_ts[i] for i in range(date_num) if i not in inps.ex_idx_list])
        # Plot excluded dates
        ax.scatter(inps.ex_dates, ex_d_ts, s=inps.marker_size ** 2, color='gray')  # color='crimson'
    # Plot kept dates
    color = pp.mplColors[0]
    if plot_num == 2:
        #color = 'crimson'
        color = pp.mplColors[1]
    #print(('Color is ' + color))
    scatter = ax.scatter(dates, d_ts, s=inps.marker_size ** 2, label='1', color=color)

    return ax, scatter


def update_timeseries(y, x, plot_number, data_only=False):
    '''Plot point time series displacement at pixel [y, x]
        Inputs:
            y           : int, y coordinate to update
            x           : int, x coordinate to update
            plot_number : int, plot number (1/2) to update
            data_only   : bool, compute and return data only, or set remainder of plot variables
        Outputs:
            d_ts        : [float], timeseries data at x, y point
    '''
    global fig_ts, ax_ts, second_plot_axis, inps, dateList, h5, k, inps, tims, fig_v, date_num, d_ts

    set_scatter_coords(plot_number, x, y)

    if plot_number == 1:
        axis = ax_ts
    else:
        axis = second_plot_axis

    d_ts = []
    for i, date in enumerate(dateList):
        d = h5['timeseries'][i][y, x]
        if inps.ref_yx:
            d -= h5['timeseries'][i][inps.ref_yx[0], inps.ref_yx[1]]
        d_ts.append(d * inps.unit_fac)

    if inps.zero_first:
        d_ts -= d_ts[inps.zero_idx]

    # Returns computed data without setting any plot or figure parameters
    if data_only:
        return d_ts

    axis.cla()
    if inps.error_file:
        axis = plot_timeseries_errorbar(ax_ts, d_ts, inps)
    else:
        axis, scatter = plot_timeseries_scatter(axis, d_ts, inps, plot_number)
        scatter.set_label('2')

    if inps.ylim:
        axis.set_ylim(inps.ylim)
    for tick in axis.yaxis.get_major_ticks():
        tick.label.set_fontsize(inps.font_size)

    # Title
    title_ts = set_axis_title(x, y)
    if inps.disp_title:
        axis.set_title(title_ts)

    axis = pp.auto_adjust_xaxis_date(axis, tims, fontsize=inps.font_size)[0]
    axis.set_xlabel('Time', fontsize=inps.font_size)
    axis.set_ylabel('Displacement [%s]' % inps.disp_unit, fontsize=inps.font_size)

    fig_v.canvas.draw()

    # Print to terminal
    print('\n---------------------------------------')
    print(title_ts)
    print(d_ts)

    # Slope estimation
    estimate_slope()

    return d_ts


def set_axis_title(x, y):
    '''Sets title of axis for a given X, Y Point
        Inputs:
            x   : int, x coordinate
            y   : int, y coordinate
        Outputs:
            title_ts    : string, computed axis title
    '''
    global lat, lon, ullon, ullat, lat_step, lon_step

    if x is None:
        title_ts = 'No Point Selected'
    else:

        title_ts = 'Y = %d, X = %d' % (y, x)
        try:
            lat, lon = xy_to_lat_lon(x, y)
            title_ts += ', lat = %.4f, lon = %.4f' % (lat, lon)
        except:
            pass

    return title_ts


def xy_to_lat_lon(x, y):
    '''Converst x,y coordinated to lat/lon coordinates
        Inputs:
            x   : int, x coordinate
            y`  : int, y coordinate
        Outputs:
            latitude    : double, computed latitude coordinate
            longitude   : double, computed longitude coordinate
    '''
    global ullat, ullon, lat_step, lon_step

    latitude = ullat + y * lat_step
    longitude = ullon + x * lon_step

    return latitude, longitude


def estimate_slope():
    '''Estimates slope of timeseries scatter data'''
    global inps, tims, d_ts, date_num

    if inps.ex_date_list:
        tims_kept = [tims[i] for i in range(date_num) if i not in inps.ex_idx_list]
        d_ts_kept = [d_ts[i] for i in range(date_num) if i not in inps.ex_idx_list]
        d_slope = stats.linregress(np.array(tims_kept), np.array(d_ts_kept))
    else:
        d_slope = stats.linregress(np.array(tims), np.array(d_ts))

    print(('linear velocity: %.2f +/- %.2f [%s/yr]' % (d_slope[0], d_slope[4], inps.disp_unit)))


def set_scatter_coords(plot_number, x, y):
    '''Sets the coordinates or the starting scatter point
        Inputs:
            plot_number     : int, the plot number (1 or 2)
            x               : int, x coordinate
            y               : int, y coordinate
    '''
    global p1_x, p1_y, p2_x, p2_y

    if plot_number == 1:        # Set scatter point 1 coordinates
        p1_x, p1_y = x, y
    else:                       # Set scatter point 2 coordinates
        p2_x, p2_y = x, y


def plot_timeseries_event(event):
    '''Event function to get y/x from button press'''
    global ax_v, d_ts, p1_scatter_point, p2_scatter_point, second_plot_axis, p1_x, p1_y, p2_x, p2_y

    if event.inaxes != ax_v:
        return

    ii = int(event.ydata + 0.5)
    jj = int(event.xdata + 0.5)

    if event.button == 1:       # Compute and update plot 1 data on left mouse-click

        if p1_scatter_point is not None:
            p1_scatter_point.remove()   # remove previous scatter point

        p1_scatter_point = ax_v.scatter(event.xdata, event.ydata, s=50, c='red', marker='o') # place new sactter point

        d_ts = update_timeseries(ii, jj, 1)     # update timeseries scatter plot for plot 1

    elif event.button == 3 and second_plot_axis_visible:    # COmpute and update plot 2 on right mouse-click

        if p2_scatter_point is not None:
            p2_scatter_point.remove()   # remove previous scatter point

        p2_scatter_point = ax_v.scatter(event.xdata, event.ydata, s=50, c='blue', marker='o') # place new scatter point

        d_ts = update_timeseries(ii, jj, 2)     # update timeseries scatter plot for plot 2


'''Displays second data plot to screen'''
def show_second_plot(event):

    global fig_v, second_plot_axis, second_plot_axis_visible

    second_plot_axis = fig_v.add_axes([0.45, 0.18, 0.52, 0.3])
    second_plot_axis_visible = True

    fig_v.canvas.draw()


'''Hides second data plot from screen'''
def hide_second_plot(event):
    global second_plot_axis, fig_v, p2_scatter_point, second_plot_axis_visible

    if p2_scatter_point is not None:
        p2_scatter_point.remove()
        p2_scatter_point = None

    second_plot_axis.remove()

    second_plot_axis_visible = False

    fig_v.canvas.draw()


'''Displays Scatter Plot Data from one or both data axes in separate figure for anlaysis'''
def show_data_as_fig(event):
    global second_plot_axis, ax_ts, second_plot_axis_visible

    if ax_ts == event.inaxes or second_plot_axis == event.inaxes:
        show_figure(1)
        if second_plot_axis_visible:
            show_figure(2)


# Configures and Shows Data Plot as Separate Figure Window
def show_figure(plot_number):
    '''Configures and shows timeseries scatter plot as separate figure window
        Inputs:
            plot_number : int, plot number to show data of
    '''
    global p2_x, p2_y, p1_x, p1_y, ax_ts, inps, plot_figure, p1_scatter, p2_scatter, new_axes, annot

    # Set up new plot figure, window, and axes
    plot_figure = plt.figure("PLOT!!", figsize=(10, 5))
    new_axes = plot_figure.add_subplot(111)
    #new_axes.set_ylim(inps.ylim)

    # Set annotations for new acis
    annot = new_axes.annotate("",
                              xy=(0, 0),
                              xytext=(445, 10),
                              textcoords="axes points",
                              bbox=dict(boxstyle="round", fc="w"))
    annot.set_visible(False)

    # Compute timeseries data
    d_ts_n = set_timeseries_data(plot_number)

    # Compute and plot scatter points on acis
    scatter = plot_timeseries_scatter(new_axes, d_ts_n, inps, plot_number)

    # Set appropriate scatter variable
    if plot_number == 1:
        _, p1_scatter = scatter
    elif plot_number == 2:
        _, p2_scatter = scatter

    set_title_and_legend(new_axes)

    # Connect events to canvas
    plot_figure.canvas.mpl_connect('pick_event', hide_scatter)
    plot_figure.canvas.mpl_connect('motion_notify_event', on_hover)

    plot_figure.show()
    plot_figure.canvas.draw()

''' Defines behavior when hovering over a given data point (ie. showing annotation)'''
def on_hover(event):
    global plot_figure, annot, p1_scatter, p2_scatter, new_axes

    vis = annot.get_visible()
    if event.inaxes == new_axes:
        cont, ind = p1_scatter.contains(event)
        if cont:
            update_annot(ind, p1_scatter)
            annot.set_visible(True)
            plot_figure.canvas.draw_idle()
        else:
            cont, ind = p2_scatter.contains(event) if p2_scatter is not None else (False, 0)
            if cont:
                update_annot(ind, p2_scatter)
                annot.set_visible(True)
                plot_figure.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    plot_figure.canvas.draw_idle()

''' Defins annotation styles when hovering over a data point '''
def update_annot(ind, sc):
    global p1_x, p1_y, p2_x, p2_y, annot, p1_scatter, p2_scatter, tims, lat, lon

    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos

    if sc is p1_scatter and p1_x is not None:
        data = update_timeseries(p1_y, p1_x, 1, True)
        latitude, longitude = xy_to_lat_lon(p1_x, p1_y)
    elif sc is p2_scatter and p2_x is not None:
        data = update_timeseries(p2_y, p2_x, 2, True)
        latitude, longitude = xy_to_lat_lon(p2_x, p2_y)
    else:
        data = np.zeros(len(tims))
        latitude, longitude = None, None

    raw_date = str(dateList[ind["ind"][0]])
    date = list(raw_date)
    date.insert(4, '-')
    date.insert(7, '-')
    date = "".join(date)
    datum = str(data[ind["ind"][0]])

    text = "(%.4f , %.4f)" % (latitude, longitude)
    text += "\nDate: "+date+"\n"+datum
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor('b')
    annot.get_bbox_patch().set_alpha(0.4)


'''Hides Scatter Plot Data on Data Point Figure on Legend Item Click'''
def hide_scatter(event):
    global scatts, plot_figure

    legline = event.artist
    origline = scatts[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)

    # Change the alpha on the line in the legend so we can see what lines
    # have been toggled
    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)

    plot_figure.canvas.draw_idle()


'''Sets title and legend information in Data Point Figure'''
def set_title_and_legend(axis):
    global p1_x, p1_y, p2_x, p2_y, inps, p1_scatter, p2_scatter, scatts

    # Compute title based off lat/lon coords
    series_label_1 = set_axis_title(p1_x, p1_y)
    series_label_2 = None

    title = series_label_1

    if p2_x is not None:
        series_label_2 = set_axis_title(p2_x, p2_y)
        title += " vs " + series_label_2

    # Display title
    if inps.disp_title:
        axis.set_title(title)

    # Set Legend
    legend = axis.legend((p1_scatter, p2_scatter), (series_label_1, series_label_2), fancybox=True)
    legend.get_frame().set_alpha(0.4)
    scatters = [p1_scatter, p2_scatter]
    scatts = dict()

    for legline, scatter in zip(legend.legendHandles, scatters):
        if legline is not None:
            legline.set_picker(5)  # 5 pts tolerance
            scatts[legline] = scatter

''' Sets timeseries data (x/y points) prior to computing timeseries data'''
def set_timeseries_data(plot_number):
    global p1_y, p1_x, p2_y, p2_x

    x_point, y_point = p1_x, p1_y

    if plot_number == 2:
        x_point = p2_x
        y_point = p2_y

    return compute_timeseries_data(plot_number, x_point, y_point)

''' Computes timeseries data for a given x, y points '''
def compute_timeseries_data(plot_number, x_point, y_point):
    global tims

    if x_point is not None:
        d_ts_n = update_timeseries(y_point, x_point, plot_number)
    else:
        d_ts_n = np.zeros(len(tims))

    return d_ts_n


######################## MAIN FUNCTION ########################
def main(iargs=None):
    global fig_v, ax_v, inps, ax_ts, fig_ts, second_plot_axis

    inps = cmd_line_parse(iargs)

    setup_plot()

    ########## Main Figure- Cumulative Displacement Map
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    fig_v = plt.figure('Cumulative Time-series Displacement', figsize=inps.fig_size)

    ######### Map Axis - Displacement Map Axis
    ax_v = fig_v.add_axes([0.08, 0.25, 0.25, 0.70])

    configure_plot()

    ########## Plot Axes - Time Series Displacement - Points
    ax_ts = fig_v.add_axes([0.45, 0.62, 0.53, 0.3])
    second_plot_axis = fig_v.add_axes([0.45, 0.18, 0.53, 0.3])
    second_plot_axis.remove()

    # Read Error List
    read_error_list()
    # Plot Data from Initial Point on Map
    plot_data_from_inital_point()

    ######### Second Plot Axis Buttons - Show/Hide Axis and Data
    ax_button_show = fig_v.add_axes([0.8, 0.03, 0.18, 0.045])
    show_button = Button(ax_button_show, "Display Second Plot")
    show_button.on_clicked(show_second_plot)

    ax_button_hide = fig_v.add_axes([0.61, 0.03, 0.18, 0.045])
    hide_button = Button(ax_button_hide, "Hide Second Plot")
    hide_button.on_clicked(hide_second_plot)

    ########## Output
    save_output()

    ########## MPL Connection Actions
    first_data_point = fig_v.canvas.mpl_connect('button_press_event', plot_timeseries_event)
    show_data_figure = fig_v.canvas.mpl_connect('button_press_event', show_data_as_fig)

    display_figure()

    ########## MPL Disconnect Actions
    fig_v.canvas.mpl_disconnect(first_data_point)
    fig_v.canvas.mpl_disconnect(show_data_figure)

###########################################################################################
if __name__ == '__main__':
    main()
