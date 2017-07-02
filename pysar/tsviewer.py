#! /usr/bin/env python
# Adopted from plotts.py from GIAnT v1.0 for PySAR products

import os
import sys
import argparse
from datetime import datetime as dt

import h5py
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar.view as view
from pysar.mask import mask_matrix


###########################################################################################
EXAMPLE='''example:
  tsviewer.py timeseries.h5 --ylim -10 10
  tsviewer.py timeseries_demErr_plane.h5 -n 5 -m maskTempCoh.h5
  tsviewer.py timeseries_demErr_plane.h5 --yx 300 400 --nodisplay --zero-first
  tsviewer.py geo_timeseries_demErr_plane.h5 --lalo 33.250 131.665 --nodisplay
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Interactive time-series viewer',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('timeseries_file', help='time series file to display')
    parser.add_argument('-n', dest='epoch_num', metavar='NUM', type=int, default='-2',\
                        help='Epoch/slice number to display, default: the 2nd last.')
    parser.add_argument('-m','--mask', dest='mask_file', default='mask.h5', help='mask to use. Default: mask.h5')

    pixel = parser.add_argument_group('Pixel Input')
    pixel.add_argument('--yx', type=int, nargs=2, help='initial pixel to plot in Y/X coord')
    pixel.add_argument('--lalo', type=float, nargs=2, help='initial pixel to plot in lat/lon coord')
    pixel.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='change reference pixel to input location')
    pixel.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2, help='change reference pixel to input location')

    output = parser.add_argument_group('Output Setting')
    output.add_argument('--save', action='store_true', dest='save_fig',\
                        help='save data and plot to files')
    output.add_argument('--nodisplay', action='store_false', dest='disp_fig',\
                        help='save data and plot to files and do not display figures\n')
    output.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=150,\
                        help='DPI - dot per inch - for display/write')

    disp = parser.add_argument_group('Display Setting')
    disp.add_argument('--figsize', dest='fig_size', metavar=('WID','LEN'), type=float, nargs=2, default=[10.0,5.0],\
                      help='Figure size in inches - width and length. Default: 10.0 5.0\n'+\
                           'i.e. 3.5 2 for ppt; ')
    disp.add_argument('--ylim', dest='ylim', nargs=2, type=float, help='Y Limits for plotting.')
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


    inps = parser.parse_args()
    if not inps.disp_fig and not inps.save_fig:
        inps.save_fig = True
    if inps.ylim:
        inps.ylim = sorted(inps.ylim)
    return inps


###########################################################################################
if __name__ == '__main__':
    #######Actual code.
    inps = cmdLineParse()

    # Time Series Info
    atr = readfile.read_attribute(inps.timeseries_file)
    k = atr['FILE_TYPE']
    print 'input file is '+k+': '+inps.timeseries_file
    if not k == 'timeseries':
        raise ValueError('Only timeseries file is supported!')

    h5 = h5py.File(inps.timeseries_file,'r')
    dateList = sorted(h5[k].keys())
    dates, tims = ptime.date_list2vector(dateList)

    # Read exclude dates
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
            inps.ex_date_list = list(set(inps.ex_date_list).intersection(dateList))
            inps.ex_dates = ptime.date_list2vector(inps.ex_date_list)[0]
            inps.ex_idx_list = [dateList.index(i) for i in inps.ex_date_list]
            print 'exclude date:'+str(inps.ex_date_list)

    # File Size
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])
    print 'data size in [y0,y1,x0,x1]: [%d, %d, %d, %d]' % (0, length, 0, width)
    try:
        ullon = float(atr['X_FIRST'])
        ullat = float(atr['Y_FIRST'])
        lon_step = float(atr['X_STEP'])
        lat_step = float(atr['Y_STEP'])
        lrlon = ullon + width*lon_step
        lrlat = ullat + length*lat_step
        print 'data size in [lat0,lat1,lon0,lon1]: [%.4f, %.4f, %.4f, %.4f]' % (lrlat, ullat, ullon, lrlon)
    except:
        pass
    
    # Initial Pixel Coord
    if inps.lalo and 'Y_FIRST' in atr.keys():
        y = int((inps.lalo[0] - ullat) / lat_step + 0.5)
        x = int((inps.lalo[1] - ullon) / lon_step + 0.5)
        inps.yx = [y, x]

    if inps.ref_lalo and 'Y_FIRST' in atr.keys():
        y = int((inps.ref_lalo[0] - ullat) / lat_step + 0.5)
        x = int((inps.ref_lalo[1] - ullon) / lon_step + 0.5)
        inps.ref_yx = [y, x]

    # Display Unit
    if inps.disp_unit == 'cm': inps.unit_fac = 100.0
    elif inps.disp_unit == 'm': inps.unit_fac = 1.0
    elif inps.disp_unit == 'dm': inps.unit_fac = 10.0
    elif inps.disp_unit == 'mm': inps.unit_fac = 1000.0
    elif inps.disp_unit == 'km': inps.unit_fac = 0.001
    else: raise ValueError('Un-recognized unit: '+inps.disp_unit)
    print 'data    unit: m'
    print 'display unit: '+inps.disp_unit

    # Flip up-down / left-right
    if inps.auto_flip:
        inps.flip_lr, inps.flip_ud = view.auto_flip_direction(atr)
    else:
        inps.flip_ud = False
        inps.left_lr = False

    # Mask file
    try:
        mask = readfile.read(inps.mask_file)[0]
        print 'load mask from file: '+inps.mask_file
        mask[mask!=0] = 1
    except:
        mask = None
        print 'No mask used.'

    # Initial Map
    d_v = h5[k].get(dateList[inps.epoch_num])[:]*inps.unit_fac
    if inps.ref_date:
        inps.ref_d_v = h5[k].get(inps.ref_date)[:]*inps.unit_fac
        d_v -= inps.ref_d_v
    if mask is not None:
        d_v = mask_matrix(d_v, mask)
    if inps.ref_yx:
        d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]
    data_lim = [np.nanmin(d_v), np.nanmax(d_v)]

    if not inps.ylim:
        inps.ylim = data_lim
    print 'Initial data range: '+str(data_lim)
    print 'Display data range: '+str(inps.ylim)

    ########## Fig 1 - Cumulative Displacement Map
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    fig_v = plt.figure('Cumulative Displacement')

    # Axes 1 
    #ax_v = fig_v.add_subplot(111)
    #ax_v.set_position([0.125,0.25,0.75,0.65])
    #This works on OSX. Original worked on Linux.
    # rect[left, bottom, width, height]
    ax_v = fig_v.add_axes([0.125,0.25,0.75,0.65])
    img = ax_v.imshow(d_v, cmap=inps.colormap, clim=inps.ylim)

    # Reference Pixel
    if inps.ref_yx:
        d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]
        ax_v.plot(inps.ref_yx[1], inps.ref_yx[0], 'ks', ms=6)
    else:
        try:
            ax_v.plot(int(atr['ref_x']), int(atr['ref_y']), 'ks', ms=6)
        except:
            pass

    # Initial Pixel
    if inps.yx:
        ax_v.plot(inps.yx[1], inps.yx[0], 'ro', markeredgecolor='black')

    ax_v.set_xlim(0, np.shape(d_v)[1])
    ax_v.set_ylim(np.shape(d_v)[0], 0)

    # Status Bar
    def format_coord(x,y):
        global d_v
        col = int(x+0.5)
        row = int(y+0.5)
        if 0<=col<width and 0<=row<length:
            z = d_v[row,col]
            try:
                lon = ullon + x*lon_step
                lat = ullat + y*lat_step
                return 'x=%.0f, y=%.0f, value=%.4f, lon=%.4f, lat=%.4f'%(x,y,z,lon,lat)
            except:
                return 'x=%.0f, y=%.0f, value=%.4f'%(x,y,z)
    ax_v.format_coord = format_coord

    # Title and Axis Label
    ax_v.set_title('N = %d, Time = %s' % (inps.epoch_num, dates[inps.epoch_num].strftime('%Y-%m-%d')))
    if not 'Y_FIRST' in atr.keys():
        ax_v.set_xlabel('Range')
        ax_v.set_ylabel('Azimuth')

    # Flip axis
    if inps.flip_lr:
        ax_v.invert_xaxis()
        print 'flip map left and right'
    if inps.flip_ud:
        ax_v.invert_yaxis()
        print 'flip map up and down'

    # Colorbar
    cbar = fig_v.colorbar(img, orientation='vertical')
    cbar.set_label('Displacement [%s]' % inps.disp_unit)

    # Axes 2 - Time Slider
    ax_time = fig_v.add_axes([0.125,0.1,0.6,0.07], axisbg='lightgoldenrodyellow', yticks=[])
    tslider = Slider(ax_time, 'Years', tims[0], tims[-1], valinit=tims[inps.epoch_num])
    tslider.ax.bar(tims, np.ones(len(tims)), facecolor='black', width=0.01, ecolor=None)
    tslider.ax.set_xticks(np.round(np.linspace(tims[0],tims[-1],num=5)*100)/100)

    def time_slider_update(val):
        '''Update Displacement Map using Slider'''
        global fig_v,ax_v,img,mask,inps,dates,tims
        timein = tslider.val
        idx_nearest = np.argmin(np.abs(np.array(tims)-timein))
        ax_v.set_title('N = %d, Time = %s' % (idx_nearest, dates[idx_nearest].strftime('%Y-%m-%d')))
        d_v = h5[k].get(dateList[idx_nearest])[:]*inps.unit_fac
        if inps.ref_date:
            d_v -= inps.ref_d_v
        if mask is not None:
            d_v = mask_matrix(d_v, mask)
        if inps.ref_yx:
            d_v -= d_v[inps.ref_yx[0], inps.ref_yx[1]]
        img.set_data(d_v)
        fig_v.canvas.draw()

    tslider.on_changed(time_slider_update)


    ########## Fig 2 - Time Series Displacement - Point
    fig_ts = plt.figure('Time series - point', figsize=inps.fig_size)
    ax_ts = fig_ts.add_subplot(111)

    def update_timeseries(y, x):
        '''Plot point time series displacement at pixel [y, x]'''
        global fig_ts,ax_ts,inps,dates
        d_ts = []
        for date in dateList:
            d = h5[k].get(date)[y,x]
            if inps.ref_yx:
                d -= h5[k].get(date)[inps.ref_yx[0], inps.ref_yx[1]]
            d_ts.append(d*inps.unit_fac)
        
        if inps.zero_first:
            d_ts -= d_ts[0]

        ax_ts.cla()
        ax_ts.scatter(dates, d_ts)
        if inps.ex_date_list:
            ex_d_ts = [d_ts[i] for i in inps.ex_idx_list]
            ax_ts.scatter(inps.ex_dates, ex_d_ts, color='gray')   # color='crimson'
        ax_ts.set_ylim(inps.ylim)

        # Title
        title_ts = 'Y = %d, X = %d'%(y,x)
        try:
            lat = ullat + y*lat_step
            lon = ullon + x*lon_step
            title_ts += ', lat = %.4f, lon = %.4f' % (lat, lon)
        except:
            pass
        if inps.disp_title:
            ax_ts.set_title(title_ts)

        ax_ts = ptime.auto_adjust_xaxis_date(ax_ts, tims, fontSize=inps.font_size)[0]
        ax_ts.set_xlabel('Time', fontsize=inps.font_size)
        ax_ts.set_ylabel('Displacement [%s]' % inps.disp_unit, fontsize=inps.font_size)

        fig_ts.canvas.draw()

        # Print to terminal
        print '\n---------------------------------------'
        print title_ts
        print d_ts

        return d_ts

    # Initial point time series plot
    if inps.yx:
        d_ts = update_timeseries(inps.yx[0], inps.yx[1])
    else:
        ax_ts.scatter(dates, np.zeros(len(tims)))
        if inps.ex_date_list:
            ax_ts.scatter(inps.ex_dates, np.zeros(len(inps.ex_date_list)), color='gray')  # color='crimson'

    def plot_timeseries_event(event):
        '''Event function to get y/x from button press'''
        global ax_v
        if event.inaxes != ax_v:
            return

        ii = int(event.ydata+0.5)
        jj = int(event.xdata+0.5)
        d_ts = update_timeseries(ii, jj)


    ########## Output
    if inps.save_fig and inps.yx:
        print 'save info for pixel '+str(inps.yx)
        fig_base = 'y%d_x%d' % (inps.yx[0], inps.yx[1])

        # TXT - point time series
        outName = fig_base+'_ts.txt'
        header_info = 'timeseries_file='+inps.timeseries_file
        header_info += '\ny=%d, x=%d' % (inps.yx[0], inps.yx[1])
        try:
            lat = ullat + inps.yx[0]*lat_step
            lon = ullon + inps.yx[1]*lon_step
            header_info += '\nlat=%.6f, lon=%.6f' % (lat, lon)
        except:
            pass
        if inps.ref_yx:
            header_info += '\nreference pixel: y=%d, x=%d' % (inps.ref_yx[0], inps.ref_yx[1])
        else:
            header_info += '\nreference pixel: y=%s, x=%s' % (atr['ref_y'], atr['ref_x'])
        header_info += '\nunit=m/yr'
        np.savetxt(outName, zip(np.array(dateList), np.array(d_ts)/inps.unit_fac), fmt='%s',\
                   delimiter='    ', header=header_info)
        print 'save time series displacement in meter to '+outName

        # Figure - point time series
        outName = fig_base+'_ts.pdf'
        fig_ts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        print 'save time series plot to '+outName

        # Figure - map
        outName = fig_base+'_'+dateList[inps.epoch_num]+'.png'
        fig_v.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        print 'save map plot to '+outName


    ########## Final linking of the canvas to the plots.
    cid = fig_v.canvas.mpl_connect('button_press_event', plot_timeseries_event)
    if inps.disp_fig:
        plt.show()
    fig_v.canvas.mpl_disconnect(cid)

