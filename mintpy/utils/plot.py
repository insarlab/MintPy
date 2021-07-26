############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#     from mintpy.utils import plot as pp


import os
import argparse
import warnings
import datetime as dt
import h5py
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt, ticker, dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyproj
from cartopy import crs as ccrs
from cartopy.mpl import geoaxes, ticker as cticker

from mintpy.objects.coord import coordinate
from mintpy.objects.colors import ColormapExt
from mintpy.objects import timeseriesKeyNames, timeseriesDatasetNames
from mintpy.utils import (
    arg_group,
    ptime,
    readfile,
    network as pnet,
    utils0 as ut0,
)


min_figsize_single = 6.0       # default min size in inch, for single plot
max_figsize_single = 10.0      # default min size in inch, for single plot
# default size in inch, for multiple subplots
default_figsize_multi = [15.0, 8.0]
max_figsize_height = 8.0       # max figure size in vertical direction in inch


# default color names in matplotlib
# ref: https://matplotlib.org/users/dflt_style_changes.html
mplColors = ['#1f77b4',
             '#ff7f0e',
             '#2ca02c',
             '#d62728',
             '#9467bd',
             '#8c564b',
             '#e377c2',
             '#7f7f7f',
             '#bcbd22',
             '#17becf']


########################################### Parser utilities ##############################################
def cmd_line_parse(iargs=''):
    parser = argparse.ArgumentParser(description='Ploting Parser')
    parser = arg_group.add_data_disp_argument(parser)
    parser = arg_group.add_dem_argument(parser)
    parser = arg_group.add_figure_argument(parser)
    parser = arg_group.add_gps_argument(parser)
    parser = arg_group.add_mask_argument(parser)
    parser = arg_group.add_map_argument(parser)
    parser = arg_group.add_point_argument(parser)
    parser = arg_group.add_reference_argument(parser)
    parser = arg_group.add_save_argument(parser)
    parser = arg_group.add_subset_argument(parser)

    inps = parser.parse_args(args=iargs)
    return inps


def read_pts2inps(inps, coord_obj):
    """Read pts_* options"""
    ## 1. merge pts_file/lalo/yx into pts_yx
    # pts_file --> pts_lalo
    if inps.pts_file and os.path.isfile(inps.pts_file):
        print('read points lat/lon from text file: {}'.format(inps.pts_file))
        inps.pts_lalo = np.loadtxt(inps.pts_file, usecols=(0,1), dtype=bytes).astype(float)

    # pts_lalo --> pts_yx
    if inps.pts_lalo is not None:
        # format pts_lalo to 2D array in size of (num_pts, 2)
        inps.pts_lalo = np.array(inps.pts_lalo).reshape(-1, 2)
        # pts_lalo --> pts_yx
        inps.pts_yx = coord_obj.geo2radar(inps.pts_lalo[:, 0],
                                          inps.pts_lalo[:, 1],
                                          print_msg=False)[:2]
        inps.pts_yx = np.array(inps.pts_yx).T.reshape(-1, 2)

    ## 2. pts_yx --> pts_yx/lalo
    if inps.pts_yx is not None:
        # format pts_yx to 2D array in size of (num_pts, 2)
        inps.pts_yx = np.array(inps.pts_yx).reshape(-1, 2)

        # pts_yx --> pts_lalo
        try:
            inps.pts_lalo = coord_obj.radar2geo(inps.pts_yx[:, 0],
                                                inps.pts_yx[:, 1],
                                                print_msg=False)[:2]
            inps.pts_lalo = np.array(inps.pts_lalo).T.reshape(-1, 2)
        except:
            pass

    return inps


############################################ Plot Utilities #############################################
def add_inner_title(ax, title, loc, prop=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if prop is None:
        prop = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=prop,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return at


def auto_figure_size(ds_shape, scale=1.0, disp_cbar=False, disp_slider=False,
                     cbar_ratio=0.25, slider_ratio=0.30, print_msg=True):
    """Get auto figure size based on input data shape
    Adjust if display colobar on the right and/or slider on the bottom

    Parameters: ds_shape          - tuple/list of 2 int for the 2D matrix shape in [length, width]
                scale             - floag, scale the final figure size
                disp_cbar/slider  - bool, plot colorbar on the right / slider on the bottom
                cbar/slider_ratio - float, size ratio of the additional colobar / slider
    Returns:    figsize           - list of 2 float for the figure size in [width, lenght] in inches
    """
    # figure shape
    fig_shape = list(ds_shape)[::-1]
    if disp_cbar:
        fig_shape[0] *= (1 + cbar_ratio)
    if disp_slider:
        fig_shape[1] *= (1 + slider_ratio)

    # get scale to meet the min/max figure size constrain
    fig_scale = min(min_figsize_single / min(fig_shape),
                    max_figsize_single / max(fig_shape),
                    max_figsize_height / fig_shape[1])

    # fig_shape/scale --> fig_size
    fig_size = [i*fig_scale*scale for i in fig_shape]
    if print_msg:
        print('figure size : [{:.2f}, {:.2f}]'.format(fig_size[0], fig_size[1]))

    return fig_size


def auto_figure_title(fname, datasetNames=[], inps_dict=None):
    """Get auto figure title from meta dict and input options
    Parameters: fname : str, input file name
                datasetNames : list of str, optional, dataset to read for multi dataset/group files
                inps_dict : dict, optional, processing attributes, including:
                    ref_date
                    pix_box
                    wrap
    Returns:    fig_title : str, output figure title
    Example:    'geo_velocity.h5' = auto_figure_title('geo_velocity.h5', None, vars(inps))
                '101020-110220_ERA5_ramp_demErr' = auto_figure_title('timeseries_ERA5_ramp_demErr.h5', '110220')
    """
    if not datasetNames:
        datasetNames = []
    if isinstance(datasetNames, str):
        datasetNames = [datasetNames]

    fbase, fext = os.path.splitext(os.path.basename(fname))
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    num_pixel = int(atr['WIDTH']) * int(atr['LENGTH'])

    if k == 'ifgramStack':
        if len(datasetNames) == 1:
            fig_title = datasetNames[0]
            if 'unwCor' in fname:
                fig_title += '_unwCor'
        else:
            fig_title = datasetNames[0].split('-')[0]

    elif k in timeseriesKeyNames and len(datasetNames) == 1:
        if 'ref_date' in inps_dict.keys():
            ref_date = inps_dict['ref_date']
        elif 'REF_DATE' in atr.keys():
            ref_date = atr['REF_DATE']
        else:
            ref_date = None

        if not ref_date:
            fig_title = datasetNames[0]
        else:
            fig_title = '{}_{}'.format(ref_date, datasetNames[0])

        # grab info of applied phase corrections from filename
        try:
            processMark = os.path.basename(fname).split('timeseries')[1].split(fext)[0]
            fig_title += processMark
        except:
            pass

    elif k == 'geometry':
        if len(datasetNames) == 1:
            fig_title = datasetNames[0]
        elif datasetNames[0].startswith('bperp'):
            fig_title = 'bperp'
        else:
            fig_title = fbase

    elif fext in ['.h5','.he5']:
        fig_title = fbase

    else:
        fig_title = os.path.basename(fname)
        # show dataset name for multi-band binry files
        num_band = int(atr.get('number_bands', '1'))
        if num_band > 1 and len(datasetNames) == 1:
            fig_title += ' - {}'.format(datasetNames[0])

    # suffix for subset
    if inps_dict.get('pix_box', None):
        box = inps_dict['pix_box']
        if (box[2] - box[0]) * (box[3] - box[1]) < num_pixel:
            fig_title += '_sub'

    # suffix for re-wrapping
    if inps_dict.get('wrap', False):
        fig_title += '_wrap'
        wrap_range = inps_dict.get('wrap_range', [-1.*np.pi, np.pi])
        wrap_range = wrap_range[1] - wrap_range[0]
        if wrap_range != 2*np.pi:
            if wrap_range == int(wrap_range):
                wrap_range = int(wrap_range)
            fig_title += str(wrap_range)

    return fig_title


def auto_flip_direction(metadata, ax=None, print_msg=True):
    """Check flip left-right and up-down based on attribute dict, for radar-coded file only"""
    # default value
    flip_lr = False
    flip_ud = False

    # auto flip for file in radar coordinates
    if 'Y_FIRST' not in metadata.keys() and 'ORBIT_DIRECTION' in metadata.keys():
        msg = '{} orbit'.format(metadata['ORBIT_DIRECTION'])
        if metadata['ORBIT_DIRECTION'].lower().startswith('a'):
            flip_ud = True
            msg += ' -> flip up-down'
        else:
            flip_lr = True
            msg += ' -> flip left-right'
        if print_msg:
            print(msg)

    if ax is not None:
        if flip_lr:
            ax.invert_xaxis()
        if flip_ud:
            ax.invert_yaxis()
        return ax

    return flip_lr, flip_ud


def auto_multilook_num(box, num_time, print_msg=True):
    """Calcualte the default/auto multilook number based on the input 3D shape.
    Parameters: box           - tuple of 4 int in (x0, y0, x1, y1) for the spatial bounding box
                num_time      - int, the 3rd / time dimension size
    Returns:    multilook_num - int, multilook number
    """
    # calc total number of pixels
    num_pixel = (box[2] - box[0]) * (box[3] - box[1]) * num_time

    # calc auto multilook_num
    if   num_pixel > (64e6*320):  multilook_num = 32;      # 16k * 4k image with 320 subplots
    elif num_pixel > (32e6*160):  multilook_num = 16;      #  8k * 4k image with 160 subplots
    elif num_pixel > ( 8e6*160):  multilook_num = 8;       #  4k * 2k image with 160 subplots
    elif num_pixel > ( 4e6*80) :  multilook_num = 4;       #  2k * 2k image with 80  subplots
    elif num_pixel > ( 4e6*20) :  multilook_num = 2;       #  2k * 2k image with 20  subplots
    else:                         multilook_num = 1;

    # print out msg
    if multilook_num > 1 and print_msg:
        print('total number of pixels: {:.1E}'.format(num_pixel))
        print('* multilook {0} by {0} with nearest interpolation for display to save memory'.format(multilook_num))

    return multilook_num


def auto_row_col_num(subplot_num, data_shape, fig_size, fig_num=1):
    """Get optimal row and column number given figure size number of subplots
    Parameters: subplot_num : int, total number of subplots
                data_shape : list of 2 float, data size in pixel in row and column direction of each plot
                fig_size : list of 2 float, figure window size in inches
                fig_num : int, number of figure windows, optional, default = 1.
    Returns:    row_num : number of subplots in row    direction per figure
                col_num : number of subplots in column direction per figure
    """
    subplot_num_per_fig = int(np.ceil(float(subplot_num) / float(fig_num)))

    data_shape_ratio = float(data_shape[0]) / float(data_shape[1])
    num_ratio = fig_size[1] / fig_size[0] / data_shape_ratio
    row_num = max(np.sqrt(subplot_num_per_fig * num_ratio), 1.)
    col_num = max(np.sqrt(subplot_num_per_fig / num_ratio), 1.)

    if row_num == 1.:
        col_num = subplot_num_per_fig
    elif col_num == 1.:
        row_num = subplot_num_per_fig

    while np.rint(row_num) * np.rint(col_num) < subplot_num_per_fig:
        if row_num % 1 > col_num % 1:
            row_num += 0.5
        else:
            col_num += 0.5

    row_num = int(np.rint(row_num))
    col_num = int(np.rint(col_num))
    return row_num, col_num


def auto_shared_lalo_location(axs, loc=(1,0,0,1), flatten=False):
    """Return the auto lat/lon label location of subplots
    Parameters: axs : 2D np.ndarray of matplotlib.axes._subplots.AxesSubplot object
                loc : tuple of 4 bool, for (left, right, top, bottom)
                flatten : bool, return variable in 2D np.ndarray or in list of flattened array
    Returns:    locs : 2D np.ndarray of tuple of 4 bool.
    """
    nrows, ncols = axs.shape
    locs = np.zeros([nrows, ncols, 4], dtype=int)
    locs[ :, 0,0] = loc[0]
    locs[ :,-1,1] = loc[1]
    locs[ 0, :,2] = loc[2]
    locs[-1, :,3] = loc[3]

    if flatten:
        loc_list = list(locs.tolist())
        locs = []
        for i in range(nrows):
            for j in range(ncols):
                locs.append(loc_list[i][j])
    return locs


def auto_colormap_name(metadata, cmap_name=None, datasetName=None, print_msg=True):
    gray_dataset_key_words = ['coherence', 'temporalCoherence',
                              'waterMask', 'shadowMask',
                              '.cor', '.mli', '.slc', '.amp', '.ramp']
    if not cmap_name:
        if any(i in gray_dataset_key_words for i in [metadata['FILE_TYPE'],
                                                     str(datasetName).split('-')[0]]):
            cmap_name = 'gray'
        else:
            cmap_name = 'jet'
    if print_msg:
        print('colormap:', cmap_name)

    return cmap_name


def auto_adjust_colormap_lut_and_disp_limit(data, num_multilook=1, max_discrete_num_step=20, print_msg=True):

    # max step size / min step number for a uniform colormap
    unique_values = np.unique(data[~np.isnan(data) * np.isfinite(data)])
    min_val = np.min(unique_values).astype(float)
    max_val = np.max(unique_values).astype(float)

    if min_val == max_val:
        cmap_lut = 256
        vlim = [min_val, max_val]

    else:
        vstep = np.min(np.abs(np.diff(unique_values))).astype(float)
        min_num_step = int((max_val - min_val) / vstep + 1)

        # use discrete colromap for data with uniform AND limited unique values
        # e.g. unwrapError time-series
        if min_num_step <= max_discrete_num_step:
            cmap_lut = min_num_step
            vlim = [min_val - vstep/2, max_val + vstep/2]

            if print_msg:
                msg = 'data has uniform and limited number ({} <= {})'.format(unique_values.size, max_discrete_num_step)
                msg += ' of unique values --> discretize colormap'
                print(msg)

        else:
            from mintpy.multilook import multilook_data
            data_mli = multilook_data(data, num_multilook, num_multilook)

            cmap_lut = 256
            vlim = [np.nanmin(data_mli), np.nanmax(data_mli)]

    return cmap_lut, vlim


def auto_adjust_xaxis_date(ax, datevector, fontsize=12, every_year=1, buffer_year=0.2):
    """Adjust X axis
    Input:
        ax          - matplotlib figure axes object
        datevector  - list of float, date in years
                         i.e. [2007.013698630137, 2007.521917808219, 2007.6463470319634]
                      OR list of datetime.datetime objects
        every_year  - int, number of years per major locator
        buffer_year - float in years, None for keep the original xlim range.
    Output:
        ax  - matplotlib figure axes object
        dss - datetime.datetime object, xmin
        dee - datetime.datetime object, xmax
    """
    # convert datetime.datetime format into date in years
    if isinstance(datevector[0], dt.datetime):
        datevector = [i.year + (i.timetuple().tm_yday-1)/365.25 for i in datevector]

    # Min/Max
    if buffer_year is not None:
        ts = datevector[0]  - buffer_year;        ys=int(ts);  ms=int((ts - ys) * 12.0)
        te = datevector[-1] + buffer_year + 0.1;  ye=int(te);  me=int((te - ye) * 12.0)
        if ms > 12:   ys = ys + 1;   ms = 1
        if me > 12:   ye = ye + 1;   me = 1
        if ms < 1:    ys = ys - 1;   ms = 12
        if me < 1:    ye = ye - 1;   me = 12
        dss = dt.datetime(ys, ms, 1, 0, 0)
        dee = dt.datetime(ye, me, 1, 0, 0)
    else:
        (dss, dee) = ax.get_xlim()
    ax.set_xlim(dss, dee)

    # Label/Tick format
    ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
    ax.xaxis.set_major_locator(mdates.YearLocator(every_year))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    # Label font size
    ax.tick_params(labelsize=fontsize)
    # fig2.autofmt_xdate()     #adjust x overlap by rorating, may enble again
    return ax, dss, dee


def auto_adjust_yaxis(ax, dataList, fontsize=12, ymin=None, ymax=None):
    """Adjust Y axis
    Input:
        ax       : matplot figure axes object
        dataList : list of float, value in y axis
        fontsize : float, font size
        ymin     : float, lower y axis limit
        ymax     : float, upper y axis limit
    Output:
        ax
    """
    # Min/Max
    dataRange = max(dataList) - min(dataList)
    if ymin is None:
        ymin = min(dataList) - 0.1*dataRange
    if ymax is None:
        ymax = max(dataList) + 0.1*dataRange
    ax.set_ylim([ymin, ymax])
    # Tick/Label setting
    #xticklabels = plt.getp(ax, 'xticklabels')
    #yticklabels = plt.getp(ax, 'yticklabels')
    #plt.setp(yticklabels, 'color', 'k', fontsize=fontsize)
    #plt.setp(xticklabels, 'color', 'k', fontsize=fontsize)

    return ax


####################################### Plot ################################################
def plot_coherence_history(ax, date12List, cohList, p_dict={}):
    """Plot min/max Coherence of all interferograms for each date"""
    # Figure Setting
    if 'ds_name'     not in p_dict.keys():   p_dict['ds_name']     = 'Coherence'
    if 'fontsize'    not in p_dict.keys():   p_dict['fontsize']    = 12
    if 'linewidth'   not in p_dict.keys():   p_dict['linewidth']   = 2
    if 'markercolor' not in p_dict.keys():   p_dict['markercolor'] = 'orange'
    if 'markersize'  not in p_dict.keys():   p_dict['markersize']  = 12
    if 'disp_title'  not in p_dict.keys():   p_dict['disp_title']  = True
    if 'every_year'  not in p_dict.keys():   p_dict['every_year']  = 1
    if 'vlim'        not in p_dict.keys():   p_dict['vlim']        = [0.2, 1.0]

    # Get date list
    date12List = ptime.yyyymmdd_date12(date12List)
    m_dates = [date12.split('_')[0] for date12 in date12List]
    s_dates = [date12.split('_')[1] for date12 in date12List]
    dateList = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))

    dates, datevector = ptime.date_list2vector(dateList)
    bar_width = ut0.most_common(np.diff(dates).tolist())*3/4
    x_list = [i-bar_width/2 for i in dates]

    coh_mat = pnet.coherence_matrix(date12List, cohList)

    ax.bar(x_list, np.nanmax(coh_mat, axis=0), bar_width.days, label='Max {}'.format(p_dict['ds_name']))
    ax.bar(x_list, np.nanmin(coh_mat, axis=0), bar_width.days, label='Min {}'.format(p_dict['ds_name']))

    if p_dict['disp_title']:
        ax.set_title('{} History of All Related Pairs'.format(p_dict['ds_name']))

    ax = auto_adjust_xaxis_date(ax, datevector, fontsize=p_dict['fontsize'],
                                every_year=p_dict['every_year'])[0]
    ax.set_ylim([p_dict['vlim'][0], p_dict['vlim'][1]])

    ax.set_xlabel('Time [years]', fontsize=p_dict['fontsize'])
    ax.set_ylabel(p_dict['ds_name'], fontsize=p_dict['fontsize'])
    ax.legend(loc='lower right')

    return ax


def plot_network(ax, date12List, dateList, pbaseList, p_dict={}, date12List_drop=[], print_msg=True):
    """Plot Temporal-Perp baseline Network
    Inputs
        ax : matplotlib axes object
        date12List : list of string for date12 in YYYYMMDD_YYYYMMDD format
        dateList   : list of string, for date in YYYYMMDD format
        pbaseList  : list of float, perp baseline, len=number of acquisition
        p_dict   : dictionary with the following items:
                      fontsize
                      linewidth
                      markercolor
                      markersize

                      cohList : list of float, coherence value of each interferogram, len = number of ifgrams
                      colormap : string, colormap name
                      disp_title : bool, show figure title or not, default: True
                      disp_drop: bool, show dropped interferograms or not, default: True
    Output
        ax : matplotlib axes object
    """

    # Figure Setting
    if 'fontsize'    not in p_dict.keys():  p_dict['fontsize']    = 12
    if 'linewidth'   not in p_dict.keys():  p_dict['linewidth']   = 2
    if 'markercolor' not in p_dict.keys():  p_dict['markercolor'] = 'orange'
    if 'markersize'  not in p_dict.keys():  p_dict['markersize']  = 12

    # For colorful display of coherence
    if 'cohList'     not in p_dict.keys():  p_dict['cohList']     = None
    if 'xlabel'      not in p_dict.keys():  p_dict['xlabel']      = 'Time [years]'
    if 'ylabel'      not in p_dict.keys():  p_dict['ylabel']      = 'Perp Baseline [m]'
    if 'cbar_label'  not in p_dict.keys():  p_dict['cbar_label']  = 'Average Spatial Coherence'
    if 'cbar_size'   not in p_dict.keys():  p_dict['cbar_size']   = '3%'
    if 'disp_cbar'   not in p_dict.keys():  p_dict['disp_cbar']   = True
    if 'colormap'    not in p_dict.keys():  p_dict['colormap']    = 'RdBu'
    if 'vlim'        not in p_dict.keys():  p_dict['vlim']        = [0.2, 1.0]
    if 'disp_title'  not in p_dict.keys():  p_dict['disp_title']  = True
    if 'disp_drop'   not in p_dict.keys():  p_dict['disp_drop']   = True
    if 'disp_legend' not in p_dict.keys():  p_dict['disp_legend'] = True
    if 'every_year'  not in p_dict.keys():  p_dict['every_year']  = 1
    if 'number'      not in p_dict.keys():  p_dict['number']      = None

    # support input colormap: string for colormap name, or colormap object directly
    if isinstance(p_dict['colormap'], str):
        cmap = ColormapExt(p_dict['colormap']).colormap
    elif isinstance(p_dict['colormap'], mpl.colors.LinearSegmentedColormap):
        cmap = p_dict['colormap']
    else:
        raise ValueError('unrecognized colormap input: {}'.format(p_dict['colormap']))

    cohList = p_dict['cohList']
    transparency = 0.7

    # Date Convert
    dateList = ptime.yyyymmdd(sorted(dateList))
    dates, datevector = ptime.date_list2vector(dateList)
    tbaseList = ptime.date_list2tbase(dateList)[0]

    ## maxBperp and maxBtemp
    date12List = ptime.yyyymmdd_date12(date12List)
    ifgram_num = len(date12List)
    pbase12 = np.zeros(ifgram_num)
    tbase12 = np.zeros(ifgram_num)
    for i in range(ifgram_num):
        m_date, s_date = date12List[i].split('_')
        m_idx = dateList.index(m_date)
        s_idx = dateList.index(s_date)
        pbase12[i] = pbaseList[s_idx] - pbaseList[m_idx]
        tbase12[i] = tbaseList[s_idx] - tbaseList[m_idx]
    if print_msg:
        print('max perpendicular baseline: {:.2f} m'.format(np.max(np.abs(pbase12))))
        print('max temporal      baseline: {} days'.format(np.max(tbase12)))

    ## Keep/Drop - date12
    date12List_keep = sorted(list(set(date12List) - set(date12List_drop)))
    idx_date12_keep = [date12List.index(i) for i in date12List_keep]
    idx_date12_drop = [date12List.index(i) for i in date12List_drop]
    if not date12List_drop:
        p_dict['disp_drop'] = False

    ## Keep/Drop - date
    m_dates = [i.split('_')[0] for i in date12List_keep]
    s_dates = [i.split('_')[1] for i in date12List_keep]
    dateList_keep = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    dateList_drop = sorted(list(set(dateList) - set(dateList_keep)))
    idx_date_keep = [dateList.index(i) for i in dateList_keep]
    idx_date_drop = [dateList.index(i) for i in dateList_drop]

    # Ploting
    if cohList is not None:
        data_min = min(cohList)
        data_max = max(cohList)
        disp_min = p_dict['vlim'][0]
        disp_max = p_dict['vlim'][1]
        if print_msg:
            print('showing coherence')
            print('data range: {}'.format([data_min, data_max]))
            print('display range: {}'.format(p_dict['vlim']))

        if p_dict['disp_cbar']:
            cax = make_axes_locatable(ax).append_axes("right", p_dict['cbar_size'], pad=p_dict['cbar_size'])
            norm = mpl.colors.Normalize(vmin=disp_min, vmax=disp_max)
            cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
            cbar.ax.tick_params(labelsize=p_dict['fontsize'])
            cbar.set_label(p_dict['cbar_label'], fontsize=p_dict['fontsize'])

        # plot low coherent ifgram first and high coherence ifgram later
        cohList_keep = [cohList[date12List.index(i)] for i in date12List_keep]
        date12List_keep = [x for _, x in sorted(zip(cohList_keep, date12List_keep))]

    # Dot - SAR Acquisition
    if idx_date_keep:
        x_list = [dates[i] for i in idx_date_keep]
        y_list = [pbaseList[i] for i in idx_date_keep]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=p_dict['markersize'], mfc=p_dict['markercolor'])
    if idx_date_drop:
        x_list = [dates[i] for i in idx_date_drop]
        y_list = [pbaseList[i] for i in idx_date_drop]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=p_dict['markersize'], mfc='gray')

    ## Line - Pair/Interferogram
    # interferograms dropped
    if p_dict['disp_drop']:
        for date12 in date12List_drop:
            date1, date2 = date12.split('_')
            idx1 = dateList.index(date1)
            idx2 = dateList.index(date2)
            x = np.array([dates[idx1], dates[idx2]])
            y = np.array([pbaseList[idx1], pbaseList[idx2]])
            if cohList is not None:
                coh = cohList[date12List.index(date12)]
                coh_norm = (coh - disp_min) / (disp_max - disp_min)
                ax.plot(x, y, '--', lw=p_dict['linewidth'], alpha=transparency, c=cmap(coh_norm))
            else:
                ax.plot(x, y, '--', lw=p_dict['linewidth'], alpha=transparency, c='k')

    # interferograms kept
    for date12 in date12List_keep:
        date1, date2 = date12.split('_')
        idx1 = dateList.index(date1)
        idx2 = dateList.index(date2)
        x = np.array([dates[idx1], dates[idx2]])
        y = np.array([pbaseList[idx1], pbaseList[idx2]])
        if cohList is not None:
            coh = cohList[date12List.index(date12)]
            coh_norm = (coh - disp_min) / (disp_max - disp_min)
            ax.plot(x, y, '-', lw=p_dict['linewidth'], alpha=transparency, c=cmap(coh_norm))
        else:
            ax.plot(x, y, '-', lw=p_dict['linewidth'], alpha=transparency, c='k')

    if p_dict['disp_title']:
        ax.set_title('Interferogram Network', fontsize=p_dict['fontsize'])

    # axis format
    ax = auto_adjust_xaxis_date(ax, datevector, fontsize=p_dict['fontsize'],
                                every_year=p_dict['every_year'])[0]
    ax = auto_adjust_yaxis(ax, pbaseList, fontsize=p_dict['fontsize'])
    ax.set_xlabel(p_dict['xlabel'], fontsize=p_dict['fontsize'])
    ax.set_ylabel(p_dict['ylabel'], fontsize=p_dict['fontsize'])
    ax.tick_params(which='both', direction='in', labelsize=p_dict['fontsize'],
                   bottom=True, top=True, left=True, right=True)

    if p_dict['number'] is not None:
        ax.annotate(p_dict['number'], xy=(0.03, 0.92), color='k',
                    xycoords='axes fraction', fontsize=p_dict['fontsize'])

    # Legend
    if p_dict['disp_drop'] and p_dict['disp_legend']:
        solid_line = mpl.lines.Line2D([], [], color='k', ls='solid',  label='Ifgram used')
        dash_line  = mpl.lines.Line2D([], [], color='k', ls='dashed', label='Ifgram dropped')
        ax.legend(handles=[solid_line, dash_line])

    return ax


def plot_perp_baseline_hist(ax, dateList, pbaseList, p_dict={}, dateList_drop=[]):
    """ Plot Perpendicular Spatial Baseline History
    Inputs
        ax : matplotlib axes object
        dateList : list of string, date in YYYYMMDD format
        pbaseList : list of float, perp baseline
        p_dict : dictionary with the following items:
                    fontsize
                    linewidth
                    markercolor
                    markersize
                    disp_title : bool, show figure title or not, default: True
                    every_year : int, number of years for the major tick on xaxis
        dateList_drop : list of string, date dropped in YYYYMMDD format
                          e.g. ['20080711', '20081011']
    Output:
        ax : matplotlib axes object
    """
    # Figure Setting
    if 'fontsize'    not in p_dict.keys():   p_dict['fontsize']    = 12
    if 'linewidth'   not in p_dict.keys():   p_dict['linewidth']   = 2
    if 'markercolor' not in p_dict.keys():   p_dict['markercolor'] = 'orange'
    if 'markersize'  not in p_dict.keys():   p_dict['markersize']  = 12
    if 'disp_title'  not in p_dict.keys():   p_dict['disp_title']  = True
    if 'every_year'  not in p_dict.keys():   p_dict['every_year']  = 1
    transparency = 0.7

    # Date Convert
    dateList = ptime.yyyymmdd(dateList)
    dates, datevector = ptime.date_list2vector(dateList)

    # Get index of date used and dropped
    # dateList_drop = ['20080711', '20081011']  # for debug
    idx_keep = list(range(len(dateList)))
    idx_drop = []
    for i in dateList_drop:
        idx = dateList.index(i)
        idx_keep.remove(idx)
        idx_drop.append(idx)

    # Plot
    # ax=fig.add_subplot(111)

    # Plot date used
    if idx_keep:
        x_list = [dates[i] for i in idx_keep]
        y_list = [pbaseList[i] for i in idx_keep]
        ax.plot(x_list, y_list, '-ko', alpha=transparency, lw=p_dict['linewidth'],
                ms=p_dict['markersize'], mfc=p_dict['markercolor'])

    # Plot date dropped
    if idx_drop:
        x_list = [dates[i] for i in idx_drop]
        y_list = [pbaseList[i] for i in idx_drop]
        ax.plot(x_list, y_list, 'ko', alpha=transparency,
                ms=p_dict['markersize'], mfc='gray')

    if p_dict['disp_title']:
        ax.set_title('Perpendicular Baseline History', fontsize=p_dict['fontsize'])

    # axis format
    ax = auto_adjust_xaxis_date(ax, datevector, fontsize=p_dict['fontsize'],
                                every_year=p_dict['every_year'])[0]
    ax = auto_adjust_yaxis(ax, pbaseList, fontsize=p_dict['fontsize'])
    ax.set_xlabel('Time [years]', fontsize=p_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]', fontsize=p_dict['fontsize'])

    return ax

def plot_rotate_diag_coherence_matrix(ax, coh_list, date12_list, date12_list_drop=[],
                                      rotate_deg=-45., cmap='RdBu', disp_half=False, disp_min=0.2):
    """Plot Rotated Coherence Matrix, suitable for Sentinel-1 data with sequential network"""
    # support input colormap: string for colormap name, or colormap object directly
    if isinstance(cmap, str):
        cmap = ColormapExt(cmap).colormap
    elif isinstance(cmap, mpl.colors.LinearSegmentedColormap):
        pass
    else:
        raise ValueError('unrecognized colormap input: {}'.format(cmap))

    #calculate coherence matrix
    coh_mat = pnet.coherence_matrix(date12_list, coh_list)
    #for date12_list_drop, set value to nan in upper triangle
    if date12_list_drop:
        m_dates = [i.split('_')[0] for i in date12_list]
        s_dates = [i.split('_')[1] for i in date12_list]
        date_list = sorted(list(set(m_dates + s_dates)))
        for date12 in date12_list_drop:
            idx1, idx2 = [date_list.index(i) for i in date12.split('_')]
            coh_mat[idx2, idx1] = np.nan

    #aux info
    num_img = coh_mat.shape[0]
    idx1, idx2 = np.where(~np.isnan(coh_mat))
    num_conn = np.max(np.abs(idx1 - idx2))

    #plot diagonal - black
    diag_mat = np.diag(np.ones(num_img))
    diag_mat[diag_mat == 0.] = np.nan
    im = ax.imshow(diag_mat, cmap='gray_r', vmin=0.0, vmax=1.0)
    im.set_transform(mpl.transforms.Affine2D().rotate_deg(rotate_deg) + ax.transData)

    im = ax.imshow(coh_mat, vmin=disp_min, vmax=1, cmap=cmap)
    im.set_transform(mpl.transforms.Affine2D().rotate_deg(rotate_deg) + ax.transData)

    #axis format
    ymax = np.sqrt(num_conn**2 / 2.) + 0.9
    ax.set_xlim(-1, np.sqrt(num_img**2 * 2)-0.7)
    if disp_half:
        ax.set_ylim(0, ymax)
    else:
        ax.set_ylim(-1.*ymax, ymax)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.axis('off')
    return ax, im


def plot_coherence_matrix(ax, date12List, cohList, date12List_drop=[], p_dict={}):
    """Plot Coherence Matrix of input network
    if date12List_drop is not empty, plot KEPT pairs in the upper triangle and
                                           ALL  pairs in the lower triangle.
    Parameters: ax : matplotlib.pyplot.Axes,
                date12List : list of date12 in YYYYMMDD_YYYYMMDD format
                cohList    : list of float, coherence value
                date12List_drop : list of date12 for date12 marked as dropped
                p_dict  : dict of plot settting
    Returns:    ax : matplotlib.pyplot.Axes
                coh_mat : 2D np.array in size of [num_date, num_date]
                im : mappable
    """
    # Figure Setting
    if 'ds_name'     not in p_dict.keys():   p_dict['ds_name']     = 'Coherence'
    if 'fontsize'    not in p_dict.keys():   p_dict['fontsize']    = 12
    if 'linewidth'   not in p_dict.keys():   p_dict['linewidth']   = 2
    if 'markercolor' not in p_dict.keys():   p_dict['markercolor'] = 'orange'
    if 'markersize'  not in p_dict.keys():   p_dict['markersize']  = 12
    if 'disp_title'  not in p_dict.keys():   p_dict['disp_title']  = True
    if 'fig_title'   not in p_dict.keys():   p_dict['fig_title']   = '{} Matrix'.format(p_dict['ds_name'])
    if 'colormap'    not in p_dict.keys():   p_dict['colormap']    = 'jet'
    if 'cbar_label'  not in p_dict.keys():   p_dict['cbar_label']  = p_dict['ds_name']
    if 'vlim'        not in p_dict.keys():   p_dict['vlim']        = (0.2, 1.0)
    if 'disp_cbar'   not in p_dict.keys():   p_dict['disp_cbar']   = True
    if 'legend_loc'  not in p_dict.keys():   p_dict['legend_loc']  = 'best'
    if 'disp_legend' not in p_dict.keys():   p_dict['disp_legend'] = True

    # support input colormap: string for colormap name, or colormap object directly
    if isinstance(p_dict['colormap'], str):
        cmap = ColormapExt(p_dict['colormap']).colormap
    elif isinstance(p_dict['colormap'], mpl.colors.LinearSegmentedColormap):
        cmap = p_dict['colormap']
    else:
        raise ValueError('unrecognized colormap input: {}'.format(p_dict['colormap']))

    date12List = ptime.yyyymmdd_date12(date12List)
    coh_mat = pnet.coherence_matrix(date12List, cohList)

    if date12List_drop:
        # Date Convert
        m_dates = [i.split('_')[0] for i in date12List]
        s_dates = [i.split('_')[1] for i in date12List]
        dateList = sorted(list(set(m_dates + s_dates)))
        # Set dropped pairs' value to nan, in upper triangle only.
        for date12 in date12List_drop:
            idx1, idx2 = [dateList.index(i) for i in date12.split('_')]
            coh_mat[idx1, idx2] = np.nan

    # Show diagonal value as black, to be distinguished from un-selected interferograms
    diag_mat = np.diag(np.ones(coh_mat.shape[0]))
    diag_mat[diag_mat == 0.] = np.nan
    im = ax.imshow(diag_mat, cmap='gray_r', vmin=0.0, vmax=1.0, interpolation='nearest')
    im = ax.imshow(coh_mat, cmap=cmap,
                   vmin=p_dict['vlim'][0],
                   vmax=p_dict['vlim'][1],
                   interpolation='nearest')

    date_num = coh_mat.shape[0]
    if date_num < 30:    tick_step = 5
    elif date_num < 50:  tick_step = 10
    elif date_num < 100: tick_step = 20
    else:                tick_step = 30
    tick_list = list(range(0, date_num, tick_step))
    ax.get_xaxis().set_ticks(tick_list)
    ax.get_yaxis().set_ticks(tick_list)
    ax.set_xlabel('Image Number', fontsize=p_dict['fontsize'])
    ax.set_ylabel('Image Number', fontsize=p_dict['fontsize'])
    ax.tick_params(which='both', direction='in', labelsize=p_dict['fontsize'],
                   bottom=True, top=True, left=True, right=True)

    if p_dict['disp_title']:
        ax.set_title(p_dict['fig_title'])

    # Colorbar
    if p_dict['disp_cbar']:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(p_dict['cbar_label'], fontsize=p_dict['fontsize'])

    # Legend
    if date12List_drop and p_dict['disp_legend']:
        ax.plot([], [], label='Upper: Ifgrams used')
        ax.plot([], [], label='Lower: Ifgrams all')
        ax.legend(loc=p_dict['legend_loc'], handlelength=0)

    return ax, coh_mat, im


def read_dem(dem_file, pix_box=None, geo_box=None, print_msg=True):
    if print_msg:
        print('reading DEM: {} ...'.format(os.path.basename(dem_file)))

    dem_metadata = readfile.read_attribute(dem_file)
    # read dem data
    if dem_metadata['FILE_TYPE'] == 'geometry':
        dsName = 'height'
    else:
        dsName = None

    # get dem_pix_box
    coord = coordinate(dem_metadata)
    if pix_box is None:
        pix_box = (0, 0, int(dem_metadata['WIDTH']), int(dem_metadata['LENGTH']))

    # Support DEM with different Resolution and Coverage
    if geo_box:
        dem_pix_box = coord.box_geo2pixel(geo_box)
    else:
        dem_pix_box = pix_box
    box2read = coord.check_box_within_data_coverage(dem_pix_box, print_msg=False)

    dem, dem_metadata = readfile.read(dem_file,
                                      datasetName=dsName,
                                      box=box2read,
                                      print_msg=print_msg)

    # if input DEM does not cover the entire AOI, fill with NaN
    if pix_box is not None and box2read != dem_pix_box:
        if print_msg:
            print('align DEM to the input data file')
        dem_tmp = np.zeros((dem_pix_box[3] - dem_pix_box[1],
                            dem_pix_box[2] - dem_pix_box[0]), dtype=dem.dtype) * np.nan
        dem_tmp[box2read[1]-dem_pix_box[1]:box2read[3]-dem_pix_box[1],
                box2read[0]-dem_pix_box[0]:box2read[2]-dem_pix_box[0]] = dem
        dem = np.array(dem_tmp)
    return dem, dem_metadata, dem_pix_box


def prepare_dem_background(dem, inps=None, print_msg=True):
    """Prepare to plot DEM on background
    Parameters: dem : 2D np.int16 matrix, dem data
                inps : Namespace with the following 4 items:
                    'disp_dem_shade'    : bool,  True/False
                    'disp_dem_contour'  : bool,  True/False
                    'dem_contour_step'  : float, 200.0
                    'dem_contour_smooth': float, 3.0
    Returns:    dem_shade : 3D np.array in size of (length, width, 4)
                dem_contour : 2D np.array in size of (length, width)
                dem_contour_sequence : 1D np.array
    Examples:   dem = readfile.read('inputs/geometryRadar.h5')[0]
                dem_shade, dem_contour, dem_contour_seq = pp.prepare_dem_background(dem=dem)
    """
    # default returns
    dem_shade = None
    dem_contour = None
    dem_contour_sequence = None

    # default inputs
    if inps is None:
        inps = cmd_line_parse()
    if inps.shade_max == 999.:
        inps.shade_max = np.nanmax(dem) + 2000

    # prepare shade relief
    if inps.disp_dem_shade:
        from matplotlib.colors import LightSource
        ls = LightSource(azdeg=inps.shade_azdeg, altdeg=inps.shade_altdeg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            dem_shade = ls.shade(dem, vert_exag=inps.shade_exag,
                                 cmap=ColormapExt('gray').colormap,
                                 vmin=inps.shade_min,
                                 vmax=inps.shade_max)
        dem_shade[np.isnan(dem_shade[:, :, 0])] = np.nan
        if print_msg:
            print('show shaded relief DEM')

    # prepare contour
    if inps.disp_dem_contour:
        from scipy import ndimage
        dem_contour = ndimage.gaussian_filter(dem, sigma=inps.dem_contour_smooth, order=0)
        dem_contour_sequence = np.arange(inps.dem_contour_step, 9000, step=inps.dem_contour_step)
        if print_msg:
            print(('show contour in step of {} m '
                   'with smoothing factor of {}').format(inps.dem_contour_step,
                                                         inps.dem_contour_smooth))
    return dem_shade, dem_contour, dem_contour_sequence


def plot_dem_background(ax, geo_box=None, dem_shade=None, dem_contour=None, dem_contour_seq=None,
                        dem=None, inps=None, print_msg=True):
    """Plot DEM as background.
    Parameters: ax : matplotlib.pyplot.Axes or BasemapExt object
                geo_box : tuple of 4 float in order of (E, N, W, S), geo bounding box
                dem_shade : 3D np.array in size of (length, width, 4)
                dem_contour : 2D np.array in size of (length, width)
                dem_contour_sequence : 1D np.array
                dem : 2D np.array of DEM data
                inps : Namespace with the following 4 items:
                    'disp_dem_shade'    : bool,  True/False
                    'disp_dem_contour'  : bool,  True/False
                    'dem_contour_step'  : float, 200.0
                    'dem_contour_smooth': float, 3.0
                    'pix_box'           : 4-tuple of int, (x0, y0, x1, y1)
    Returns:    ax : matplotlib.pyplot.Axes or BasemapExt object
    Examples:   m = pp.plot_dem_background(m, geo_box=inps.geo_box, dem=dem, inps=inps)
                ax = pp.plot_dem_background(ax=ax, geo_box=None, dem_shade=dem_shade,
                                            dem_contour=dem_contour, dem_contour_seq=dem_contour_seq)
    """
    # default inputs
    if inps is None:
        inps = cmd_line_parse()

    if all(i is None for i in [dem_shade, dem_contour, dem_contour_seq]) and dem is not None:
        (dem_shade,
         dem_contour,
         dem_contour_seq) = prepare_dem_background(dem, inps=inps, print_msg=print_msg)

    # get extent - (left, right, bottom, top) in data coordinates
    if geo_box is not None:
        geo_extent = (geo_box[0], geo_box[2],
                      geo_box[3], geo_box[1])
    else:
        if hasattr(inps, 'pix_box'):
            pix_box = tuple(inps.pix_box)
        else:
            data = [i for i in [dem, dem_shade, dem_contour] if i is not None][0]
            pix_box = (0, 0, data.shape[1], data.shape[0])
        rdr_extent = (pix_box[0]-0.5, pix_box[2]-0.5,
                      pix_box[3]-0.5, pix_box[1]-0.5)

    # plot shaded relief
    if dem_shade is not None:
        # config
        kwargs = dict(interpolation='spline16', zorder=0, origin='upper')

        # geo coordinates
        if geo_box is not None:
            ax.imshow(dem_shade, extent=geo_extent, **kwargs)

        # radar coordinates
        elif isinstance(ax, plt.Axes):
            ax.imshow(dem_shade, extent=rdr_extent, **kwargs)

    # plot topo contour
    if dem_contour is not None and dem_contour_seq is not None:
        # config
        kwargs = dict(origin='upper', colors='black',
                      linewidths=inps.dem_contour_linewidth,
                      alpha=0.5, zorder=1)
        # plot contour line above data (zorder=1) if no DEM shade
        if dem_shade is None:
            kwargs['zorder'] = 2

        # geo coordinates
        if geo_box is not None:
            yy, xx = np.mgrid[geo_box[1]:geo_box[3]:dem_contour.shape[0]*1j,
                              geo_box[0]:geo_box[2]:dem_contour.shape[1]*1j]

            ax.contour(xx, yy, dem_contour, dem_contour_seq, extent=geo_extent, **kwargs)

        # radar coordinates
        elif isinstance(ax, plt.Axes):
            ax.contour(dem_contour, dem_contour_seq, extent=rdr_extent, **kwargs)

    return ax


def plot_gps(ax, SNWE, inps, metadata=dict(), print_msg=True):
    from mintpy.objects import gps
    vprint = print if print_msg else lambda *args, **kwargs: None

    marker_size = 7
    vmin, vmax = inps.vlim
    cmap = ColormapExt(inps.colormap).colormap if isinstance(inps.colormap, str) else inps.colormap

    atr = dict(metadata)
    atr['UNIT'] = 'm'
    unit_fac = scale_data2disp_unit(metadata=atr, disp_unit=inps.disp_unit)[2]

    if not inps.gps_start_date:
        inps.gps_start_date = metadata.get('START_DATE', None)

    if not inps.gps_end_date:
        inps.gps_end_date = metadata.get('END_DATE', None)

    site_names, site_lats, site_lons = gps.search_gps(SNWE, inps.gps_start_date, inps.gps_end_date)
    if inps.ref_gps_site and inps.ref_gps_site not in site_names:
        raise ValueError('input reference GPS site "{}" not available!'.format(inps.ref_gps_site))

    k = metadata['FILE_TYPE']
    if inps.gps_component and k not in ['velocity', 'timeseries']:
        inps.gps_component = None
        vprint('WARNING: --gps-comp is not implemented for {} file yet, set --gps-comp = None and continue'.format(k))

    if inps.gps_component:
        # plot GPS velocity/displacement along LOS direction
        vprint('-'*30)
        msg = 'plotting GPS '
        msg += 'velocity' if k == 'velocity' else 'displacement'
        msg += ' in LOS direction'
        vprint(msg)
        vprint('number of available GPS stations: {}'.format(len(site_names)))
        vprint('start date: {}'.format(inps.gps_start_date))
        vprint('end   date: {}'.format(inps.gps_end_date))
        vprint('components projection: {}'.format(inps.gps_component))

        # get GPS LOS observations
        site_obs = gps.get_gps_los_obs(
            insar_file=inps.file,
            site_names=site_names,
            start_date=inps.gps_start_date,
            end_date=inps.gps_end_date,
            gps_comp=inps.gps_component,
            print_msg=print_msg,
            redo=inps.gps_redo,
        )

        # reference GPS
        if inps.ref_gps_site:
            vprint('referencing all GPS LOS observations to site: {}'.format(inps.ref_gps_site))
            ref_ind = site_names.tolist().index(inps.ref_gps_site)
            # plot label of the reference site
            ax.annotate(site_names[ref_ind], xy=(site_lons[ref_ind], site_lats[ref_ind]),
                        fontsize=inps.font_size)
            # update value
            ref_val = site_obs[ref_ind]
            if not np.isnan(ref_val):
                site_obs -= ref_val

        # scale to the same unit as InSAR
        site_obs *= unit_fac

        # plot
        for lat, lon, obs in zip(site_lats, site_lons, site_obs):
            color = cmap( (obs - vmin) / (vmax - vmin) ) if obs else 'none'
            ax.scatter(lon, lat, color=color, s=marker_size**2, edgecolors='k', zorder=10)

    else:
        # plot GPS locations only
        vprint('showing GPS locations')
        ax.scatter(site_lons, site_lats, s=marker_size**2, color='w', edgecolors='k', zorder=10)

    # plot GPS label
    if inps.disp_gps_label:
        for i in range(len(site_names)):
            ax.annotate(site_names[i], xy=(site_lons[i], site_lats[i]),
                        fontsize=inps.font_size)

    return ax


def plot_colorbar(inps, im, cax):
    # extend
    if not inps.cbar_ext:
        if   inps.vlim[0] <= inps.dlim[0] and inps.vlim[1] >= inps.dlim[1]: inps.cbar_ext='neither'
        elif inps.vlim[0] >  inps.dlim[0] and inps.vlim[1] >= inps.dlim[1]: inps.cbar_ext='min'
        elif inps.vlim[0] <= inps.dlim[0] and inps.vlim[1] <  inps.dlim[1]: inps.cbar_ext='max'
        else:  inps.cbar_ext='both'

    # orientation
    if inps.cbar_loc in ['left', 'right']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    # plot colorbar
    if inps.wrap and (inps.wrap_range[1] - inps.wrap_range[0]) == 2.*np.pi:
        cbar = plt.colorbar(im, cax=cax, orientation=orientation, ticks=[-np.pi, 0, np.pi])
        cbar.ax.set_yticklabels([r'-$\pi$', '0', r'$\pi$'])
    else:
        cbar = plt.colorbar(im, cax=cax, orientation=orientation, extend=inps.cbar_ext)

    # update ticks
    if inps.cmap_lut <= 5:
        inps.cbar_nbins = inps.cmap_lut

    if inps.cbar_nbins:
        if inps.cbar_nbins <= 2:
            # manually set tick for better positions when the color step is not a common number
            # e.g. for numInvIfgram.h5
            cbar.set_ticks(inps.dlim)
        else:
            cbar.locator = ticker.MaxNLocator(nbins=inps.cbar_nbins)
            cbar.update_ticks()

    cbar.ax.tick_params(which='both', direction='out', labelsize=inps.font_size, colors=inps.font_color)

    # add label
    if inps.cbar_label:
        cbar.set_label(inps.cbar_label, fontsize=inps.font_size, color=inps.font_color)
    elif inps.disp_unit != '1':
        cbar.set_label(inps.disp_unit,  fontsize=inps.font_size, color=inps.font_color)

    return inps, cbar



##################### Data Scale based on Unit and Wrap Range ##################
def check_disp_unit_and_wrap(metadata, disp_unit=None, wrap=False, wrap_range=[-1.*np.pi, np.pi], print_msg=True):
    """Get auto disp_unit for input dataset
    Example:
        if not inps.disp_unit:
            inps.disp_unit = pp.auto_disp_unit(atr)
    """
    # default display unit if not given
    if not disp_unit:
        k = metadata['FILE_TYPE']
        k = k.replace('.','')
        disp_unit = metadata['UNIT'].lower()
        if (k in ['timeseries', 'giantTimeseries', 'velocity', 'HDFEOS']
                and disp_unit.split('/')[0].endswith('m')):
            disp_unit = 'cm'
        elif k in ['mli', 'slc', 'amp']:
            disp_unit = 'dB'

    if wrap:
        # wrap is supported for displacement file types only
        if disp_unit.split('/')[0] not in ['radian', 'm', 'cm', 'mm', '1']:
            wrap = False
            print('WARNING: re-wrap is disabled for disp_unit = {}'.format(disp_unit))
        elif disp_unit.split('/')[0] != 'radian' and (wrap_range[1] - wrap_range[0]) == 2.*np.pi:
            disp_unit = 'radian'
            if print_msg:
                print('change disp_unit = radian due to rewrapping')

    return disp_unit, wrap


def scale_data2disp_unit(data=None, metadata=dict(), disp_unit=None):
    """Scale data based on data unit and display unit
    Inputs:
        data    : 2D np.array
        metadata  : dictionary, meta data
        disp_unit : str, display unit
    Outputs:
        data    : 2D np.array, data after scaling
        disp_unit : str, display unit
    Default data file units in MintPy are:  m, m/yr, radian, 1
    """
    if not metadata:
        metadata['UNIT'] = 'm'

    # Initial
    scale = 1.0
    data_unit = metadata['UNIT'].lower().split('/')
    disp_unit = disp_unit.lower().split('/')

    # if data and display unit is the same
    if disp_unit == data_unit:
        return data, metadata['UNIT'], scale

    # Calculate scaling factor  - 1
    # phase unit - length / angle
    if data_unit[0].endswith('m'):
        if   disp_unit[0] == 'mm': scale *= 1000.0
        elif disp_unit[0] == 'cm': scale *= 100.0
        elif disp_unit[0] == 'dm': scale *= 10.0
        elif disp_unit[0] == 'm' : scale *= 1.0
        elif disp_unit[0] == 'km': scale *= 1/1000.0
        elif disp_unit[0] in ['radians','radian','rad','r']:
            range2phase = -(4*np.pi) / float(metadata['WAVELENGTH'])
            scale *= range2phase
        else:
            print('Unrecognized display phase/length unit:', disp_unit[0])
            pass

        if   data_unit[0] == 'mm': scale *= 0.001
        elif data_unit[0] == 'cm': scale *= 0.01
        elif data_unit[0] == 'dm': scale *= 0.1
        elif data_unit[0] == 'km': scale *= 1000.

    elif data_unit[0] == 'radian':
        phase2range = -float(metadata['WAVELENGTH']) / (4*np.pi)
        if   disp_unit[0] == 'mm': scale *= phase2range * 1000.0
        elif disp_unit[0] == 'cm': scale *= phase2range * 100.0
        elif disp_unit[0] == 'dm': scale *= phase2range * 10.0
        elif disp_unit[0] == 'm' : scale *= phase2range * 1.0
        elif disp_unit[0] == 'km': scale *= phase2range * 1/1000.0
        elif disp_unit[0] in ['radians','radian','rad','r']:
            pass
        else:
            print('Unrecognized phase/length unit:', disp_unit[0])
            pass

    # amplitude/coherence unit - 1
    elif data_unit[0] == '1':
        if disp_unit[0] == 'db' and data is not None:
            disp_unit[0] = 'dB'

            if metadata['FILE_TYPE'] in ['.cor', '.int', '.unw']:
                # dB for power quantities
                data = 10 * np.log10(np.clip(data, a_min=1e-1, a_max=None))

            else:
                # dB for field quantities, e.g. amp, slc
                data = 20 * np.log10(np.clip(data, a_min=1e-1, a_max=None))

        else:
            try:
                scale /= float(disp_unit[0])
            except:
                print('Un-scalable display unit:', disp_unit[0])
    else:
        print('Un-scalable data unit:', data_unit)

    # Calculate scaling factor  - 2
    if len(data_unit) == 2:
        try:
            disp_unit[1]
            if   disp_unit[1] in ['y','yr','year'  ]: disp_unit[1] = 'year'
            elif disp_unit[1] in ['m','mon','month']: disp_unit[1] = 'mon'; scale *= 12.0
            elif disp_unit[1] in ['d','day'        ]: disp_unit[1] = 'day'; scale *= 365.25
            else: print('Unrecognized time unit for display:', disp_unit[1])
        except:
            disp_unit.append('year')
        disp_unit = disp_unit[0]+'/'+disp_unit[1]
    else:
        disp_unit = disp_unit[0]

    # scale input data
    if data is not None and scale != 1.0:
        data *= np.array(scale, dtype=data.dtype)

    return data, disp_unit, scale


def scale_data4disp_unit_and_rewrap(data, metadata, disp_unit=None, wrap=False, wrap_range=[-1.*np.pi, np.pi],
                                    print_msg=True):
    """Scale 2D matrix value according to display unit and re-wrapping flag
    Inputs:
        data - 2D np.array
        metadata  - dict, including the following attributes:
               UNIT
               FILE_TYPE
               WAVELENGTH
        disp_unit  - string, optional
        wrap - bool, optional
    Outputs:
        data
        disp_unit
        wrap
    """
    if not disp_unit:
        disp_unit, wrap = check_disp_unit_and_wrap(metadata,
                                                   disp_unit=None,
                                                   wrap=wrap,
                                                   wrap_range=wrap_range,
                                                   print_msg=print_msg)

    # Data Operation - Scale to display unit
    disp_scale = 1.0
    if not disp_unit == metadata['UNIT']:
        data, disp_unit, disp_scale = scale_data2disp_unit(data,
                                                           metadata=metadata,
                                                           disp_unit=disp_unit)

    # Data Operation - wrap
    if wrap:
        data = wrap_range[0] + np.mod(data - wrap_range[0], wrap_range[1] - wrap_range[0])
        if print_msg:
            print('re-wrapping data to {}'.format(wrap_range))
    return data, disp_unit, disp_scale, wrap


def read_mask(fname, mask_file=None, datasetName=None, box=None, xstep=1, ystep=1,
              vmin=None, vmax=None, print_msg=True):
    """Find and read mask for input data file fname
    Parameters: fname       : string, data file name/path
                mask_file   : string, optional, mask file name
                datasetName : string, optional, dataset name for HDFEOS file type
                box         : tuple of 4 int, for reading part of data
    Returns:    mask        : 2D np.ndarray in bool, mask data
                mask_file   : string, file name of mask data
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']

    # default mask file
    if (not mask_file
            and k in ['velocity', 'timeseries']
            and 'msk' not in fname):

        for mask_file in ['maskTempCoh.h5', 'maskResInv.h5']:
            if 'PhaseVelocity' in fname:
                mask_file = None #'maskSpatialCoh.h5'

            # check coordinate
            if os.path.basename(fname).startswith('geo_'):
                mask_file = 'geo_{}'.format(mask_file)

            # absolute path and file existence
            mask_file = os.path.join(os.path.dirname(fname), str(mask_file))
            if os.path.isfile(mask_file):
                break
            else:
                mask_file = None

    # Read mask file if inputed
    mask = None
    if os.path.isfile(str(mask_file)):
        try:
            atr_msk = readfile.read_attribute(mask_file)
            if all(int(atr_msk[key]) == int(atr[key]) for key in ['LENGTH','WIDTH']):
                # grab dsname for conn comp mask [None for the others]
                dsName=None
                if all(meta['FILE_TYPE'] == 'ifgramStack' for meta in [atr, atr_msk]):
                    date12 = datasetName.split('-')[1]
                    dsName = 'connectComponent-{}'.format(date12)

                # read mask data
                mask = readfile.read(mask_file,
                                     box=box,
                                     datasetName=dsName,
                                     print_msg=print_msg)[0]
                vprint('read mask from file:', os.path.basename(mask_file))

            else:
                mask_file = None
                msg = 'WARNING: input file has different size from mask file: {}'.format(mask_file)
                msg += '\n    data file {} row/column number: {} / {}'.format(fname, atr['LENGTH'], atr['WIDTH'])
                msg += '\n    mask file {} row/column number: {} / {}'.format(mask_file, atr_msk['LENGTH'], atr_msk['WIDTH'])
                msg += '\n    Continue without mask.'
                vprint(msg)

        except:
            mask_file = None
            vprint('Can not open mask file:', mask_file)

    elif k in ['HDFEOS']:
        if datasetName.split('-')[0] in timeseriesDatasetNames:
            mask_file = fname
            mask = readfile.read(fname, datasetName='mask', print_msg=print_msg)[0]
            vprint('read {} contained mask dataset.'.format(k))

    elif fname.endswith('PARAMS.h5'):
        mask_file = fname
        with h5py.File(fname, 'r') as f:
            mask = f['cmask'][:] == 1.
        vprint('read {} contained cmask dataset'.format(os.path.basename(fname)))

    # multilook
    if mask is not None and xstep * ystep > 1:
        # output size if x/ystep > 1
        xsize = int(mask.shape[1] / xstep)
        ysize = int(mask.shape[0] / ystep)

        # sampling
        mask = mask[int(ystep/2)::ystep,
                    int(xstep/2)::xstep]
        mask = mask[:ysize, :xsize]

    # set to bool type
    if mask is not None:
        mask[np.isnan(mask)] = 0

        # vmin/max
        if vmin is not None:
            mask[mask < vmin] = 0
            vprint(f'hide pixels with mask value < {vmin}')
        if vmax is not None:
            mask[mask > vmax] = 0
            vprint(f'hide pixels with mask value > {vmax}')

        mask = mask != 0

    return mask, mask_file



###############################################  Maps  ###############################################
def auto_lalo_sequence(geo_box, lalo_step=None, lalo_max_num=4, step_candidate=[1, 2, 3, 4, 5]):
    """Auto calculate lat/lon label sequence based on input geo_box
    Parameters: geo_box        : 4-tuple of float, defining UL_lon, UL_lat, LR_lon, LR_lat coordinate
                lalo_step      : float
                lalo_max_num   : int, rough major tick number along the longer axis
                step_candidate : list of int, candidate list for the significant number of step
    Returns:    lats/lons : np.array of float, sequence of lat/lon auto calculated from input geo_box
                lalo_step : float, lat/lon label step
    Example:    geo_box = (128.0, 37.0, 138.0, 30.0)
                lats, lons, step = m.auto_lalo_sequence(geo_box)
    """
    max_lalo_dist = max([geo_box[1]-geo_box[3], geo_box[2]-geo_box[0]])

    if not lalo_step:
        # Initial tick step
        lalo_step = ut0.round_to_1(max_lalo_dist/lalo_max_num)

        # reduce decimal if it ends with 8/9
        digit = np.int(np.floor(np.log10(lalo_step)))
        if str(lalo_step)[-1] in ['8','9']:
            digit += 1
            lalo_step = round(lalo_step, digit)

        # Final tick step - choose from candidate list
        lalo_step_candidate = [i*10**digit for i in step_candidate]
        distance = [(i - max_lalo_dist/lalo_max_num) ** 2 for i in lalo_step_candidate]
        lalo_step = lalo_step_candidate[distance.index(min(distance))]

    digit = np.int(np.floor(np.log10(lalo_step)))

    # Auto tick sequence
    lat_major = np.ceil(geo_box[3]/10**(digit+1))*10**(digit+1)
    lats = np.unique(np.hstack((np.arange(lat_major, lat_major-10.*max_lalo_dist, -lalo_step),
                                np.arange(lat_major, lat_major+10.*max_lalo_dist, lalo_step))))
    lats = np.sort(lats[np.where(np.logical_and(lats >= geo_box[3], lats <= geo_box[1]))])

    lon_major = np.ceil(geo_box[0]/10**(digit+1))*10**(digit+1)
    lons = np.unique(np.hstack((np.arange(lon_major, lon_major-10.*max_lalo_dist, -lalo_step),
                                np.arange(lon_major, lon_major+10.*max_lalo_dist, lalo_step))))
    lons = np.sort(lons[np.where(np.logical_and(lons >= geo_box[0], lons <= geo_box[2]))])
    return lats, lons, lalo_step, digit


def draw_lalo_label(geo_box, ax=None, lalo_step=None, lalo_loc=[1, 0, 0, 1], lalo_max_num=4,
                    font_size=12, xoffset=None, yoffset=None, projection=ccrs.PlateCarree(), print_msg=True):
    """Auto draw lat/lon label/tick based on coverage from geo_box
    Parameters: geo_box   : 4-tuple of float, (W, N, E, S) in degree
                ax        : CartoPy axes.
                lalo_step : float
                lalo_loc  : list of 4 bool, positions where the labels are drawn as in [left, right, top, bottom]
                            default: [1,0,0,1]
                lalo_max_num : int
                ...
    Example:    geo_box = (128.0, 37.0, 138.0, 30.0)
                m.draw_lalo_label(geo_box)
    """
    # default ax
    if not ax:
        ax = plt.gca()

    # default lat/lon sequences
    lats, lons, lalo_step, digit = auto_lalo_sequence(geo_box, lalo_step=lalo_step, lalo_max_num=lalo_max_num)
    if print_msg:
        print('plot lat/lon label in step of {} and location of {}'.format(lalo_step, lalo_loc))

    # ticklabel/tick style
    ax.tick_params(which='both', direction='in', labelsize=font_size,
                   left=True, right=True, top=True, bottom=True,
                   labelleft=lalo_loc[0], labelright=lalo_loc[1],
                   labeltop=lalo_loc[2], labelbottom=lalo_loc[3])
    if xoffset is not None:
        ax.tick_params(axis='x', which='major', pad=xoffset)
    if yoffset is not None:
        ax.tick_params(axis='y', which='major', pad=yoffset)

    # ticklabel symbol style
    decimal_digit = max(0, 0-digit)
    lon_formatter = cticker.LongitudeFormatter(number_format='.{}f'.format(decimal_digit))
    lat_formatter = cticker.LatitudeFormatter(number_format='.{}f'.format(decimal_digit))
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.set_xticks(lons, crs=projection)
    ax.set_yticks(lats, crs=projection)
    return ax


def draw_scalebar(ax, geo_box, unit='degrees', loc=[0.2, 0.2, 0.1], labelpad=0.05, font_size=12, color='k'):
    """draw a simple map scale from x1,y to x2,y in map projection coordinates, label it with actual distance
    ref_link: http://matplotlib.1069221.n5.nabble.com/basemap-scalebar-td14133.html
    Parameters: ax       : matplotlib.pyplot.axes object
                geo_box  : tuple of 4 float in (x0, y0, x1, y1) for (W, N, E, S) in degrees / meters
                unit     : str, coordinate unit - degrees or meters
                loc      : list of 3 float in (length, lat, lon) of scale bar center in ratio of width, relative coord
                labelpad : float
    Returns:    ax
    Example:    from mintpy.utils import plot as pp
                pp.draw_scale_bar(ax, geo_box)
    """
    geod = pyproj.Geod(ellps='WGS84')
    if not ax:
        ax = plt.gca()

    ## location - center
    lon_c = geo_box[0] + loc[1] * (geo_box[2] - geo_box[0])
    lat_c = geo_box[3] + loc[2] * (geo_box[1] - geo_box[3])

    ## length
    # 1. calc scene width in meters
    if unit.startswith('deg'):
        scene_width = geod.inv(geo_box[0], geo_box[3],
                               geo_box[2], geo_box[3])[2]
    elif unit.startswith('meter'):
        scene_width = geo_box[2] - geo_box[0]

    # 2. convert length ratio to length in meters
    length_meter = ut0.round_to_1(scene_width * loc[0])
    if length_meter > 1000.0:
        # round to the nearest km
        length_meter = np.rint(length_meter/1000.0)*1000.0

    # 3. convert length in meters to length in display coord unit
    if unit.startswith('deg'):
        lon_c2 = geod.fwd(lon_c, lat_c, 90, length_meter)[0]
        length_disp = np.abs(lon_c - lon_c2)
    elif unit.startswith('meter'):
        length_disp = length_meter

    ## starting/ending longitude
    lon0 = lon_c - length_disp / 2.0
    lon1 = lon_c + length_disp / 2.0

    ## plot scale bar
    ax.plot([lon0, lon1], [lat_c, lat_c], color=color)
    ax.plot([lon0, lon0], [lat_c, lat_c + 0.1*length_disp], color=color)
    ax.plot([lon1, lon1], [lat_c, lat_c + 0.1*length_disp], color=color)

    ## plot scale bar label
    unit = 'm'
    if length_meter >= 1000.0:
        unit = 'km'
        length_meter *= 0.001
    label = '{:.0f} {}'.format(length_meter, unit)
    txt_offset = (geo_box[1] - geo_box[3]) * labelpad

    ax.text(lon0 + length_disp / 2.0,
            lat_c + txt_offset,
            label,
            verticalalignment='center', horizontalalignment='center',
            fontsize=font_size, color=color)

    return ax
