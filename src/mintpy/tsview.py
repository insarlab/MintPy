#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Zhang Yunjun, Joshua Zahner, Heresh Fattahi, 2013 #
#############################################################


import argparse
import os
import re

import numpy as np
from matplotlib import patches, pyplot as plt, widgets
from scipy import linalg, stats

from mintpy import subset, timeseries2velocity as ts2vel, view
from mintpy.multilook import multilook_data
from mintpy.objects import HDFEOS, giantTimeseries, timeseries
from mintpy.utils import plot as pp, ptime, readfile, time_func, utils as ut


###########################################################################################
def read_init_info(inps):
    # Time Series Info
    atr = readfile.read_attribute(inps.file[0])
    atr['DATA_TYPE'] = atr.get('DATA_TYPE', 'float32')

    inps.key = atr['FILE_TYPE']
    if inps.key == 'timeseries':
        obj = timeseries(inps.file[0])
    elif inps.key == 'giantTimeseries':
        obj = giantTimeseries(inps.file[0])
    elif inps.key == 'HDFEOS':
        obj = HDFEOS(inps.file[0])
    else:
        raise ValueError(f'input file is {inps.key}, not timeseries.')
    obj.open(print_msg=inps.print_msg)
    inps.seconds = atr.get('CENTER_LINE_UTC', 0)

    if not inps.file_label:
        inps.file_label = []
        for fname in inps.file:
            fbase = os.path.splitext(os.path.basename(fname))[0]
            fbase = fbase.replace('timeseries', 'TS')
            inps.file_label.append(fbase)

    # default mask file
    if not inps.mask_file and 'msk' not in inps.file[0]:
        dir_name = os.path.dirname(inps.file[0])
        if 'Y_FIRST' in atr.keys():
            inps.mask_file = os.path.join(dir_name, 'geo_maskTempCoh.h5')
        else:
            inps.mask_file = os.path.join(dir_name, 'maskTempCoh.h5')
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None

    ## date info
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
    inps.disp_unit, inps.unit_fac = pp.scale_data2disp_unit(
        metadata=atr, disp_unit=inps.disp_unit)[1:3]

    # Read Error List
    inps.ts_plot_func = plot_ts_scatter
    inps.error_ts = None
    inps.ex_error_ts = None
    if inps.error_file:
        # assign plot function
        inps.ts_plot_func = plot_ts_errorbar

        # read error file
        error_fc = np.loadtxt(inps.error_file, dtype=bytes).astype(str)
        inps.error_ts = error_fc[:, 1].astype(np.float32)*inps.unit_fac

        # update error file with exclude date
        if inps.ex_date_list:
            e_ts = inps.error_ts[:]
            inps.ex_error_ts = e_ts[inps.ex_flag == 0]
            inps.error_ts = e_ts[inps.ex_flag == 1]

    # Zero displacement for 1st acquisition
    if inps.zero_first:
        inps.zero_idx = min(0, np.min(np.where(inps.ex_flag)[0]))

    # default lookup table file and coordinate object
    if not inps.lookup_file:
        inps.lookup_file = ut.get_lookup_file('./inputs/geometryRadar.h5')
    inps.coord = ut.coordinate(atr, inps.lookup_file)

    ## size and lalo info
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), atr)
    inps.pix_box = inps.coord.check_box_within_data_coverage(inps.pix_box)
    inps.geo_box = inps.coord.box_pixel2geo(inps.pix_box)
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    vprint('data   coverage in y/x: '+str(data_box))
    vprint('subset coverage in y/x: '+str(inps.pix_box))
    vprint('data   coverage in lat/lon: '+str(inps.coord.box_pixel2geo(data_box)))
    vprint('subset coverage in lat/lon: '+str(inps.geo_box))
    vprint('------------------------------------------------------------------------')

    # Map info - coordinate unit
    inps.coord_unit = atr.get('Y_UNIT', 'degrees').lower()
    inps.lalo_digit = ut.get_lalo_digit4display(atr, coord_unit=inps.coord_unit)
    inps = view.check_map_projection(inps, metadata=atr, print_msg=inps.print_msg)

    # calculate multilook_num
    # ONLY IF:
    #   inps.multilook is True (no --nomultilook input) AND
    #   inps.multilook_num ==1 (no --multilook-num input)
    # Note: inps.multilook is used for this check ONLY
    # Note: multilooking is only applied to the 3D data cubes and their related operations:
    # e.g. spatial indexing, referencing, etc. All the other variables are in the original grid
    # so that users get the same result as the non-multilooked version.
    if inps.multilook and inps.multilook_num == 1:
        inps.multilook_num = pp.auto_multilook_num(
            inps.pix_box, inps.num_date,
            max_memory=inps.maxMemory,
            print_msg=inps.print_msg,
        )

    ## reference pixel
    if not inps.ref_lalo and 'REF_LAT' in atr.keys():
        inps.ref_lalo = (float(atr['REF_LAT']), float(atr['REF_LON']))
    if inps.ref_lalo:
        # set longitude to [-180, 180)
        if inps.coord_unit.lower().startswith('deg') and inps.ref_lalo[1] >= 180.:
            inps.ref_lalo[1] -= 360.
        # ref_lalo --> ref_yx if not set in cmd
        if not inps.ref_yx:
            inps.ref_yx = inps.coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1], print_msg=False)[0:2]

    # use REF_Y/X if ref_yx not set in cmd
    if not inps.ref_yx and 'REF_Y' in atr.keys():
        inps.ref_yx = (int(atr['REF_Y']), int(atr['REF_X']))

    # ref_yx --> ref_lalo if in geo-coord
    # for plotting purpose only
    if inps.ref_yx and 'Y_FIRST' in atr.keys():
        inps.ref_lalo = inps.coord.radar2geo(inps.ref_yx[0], inps.ref_yx[1], print_msg=False)[0:2]

    # do not plot native reference point if it's out of the coverage due to subset
    if (inps.ref_yx and 'Y_FIRST' in atr.keys()
        and inps.ref_yx == (int(atr.get('REF_Y',-999)), int(atr.get('REF_X',-999)))
        and not (    inps.pix_box[0] <= inps.ref_yx[1] < inps.pix_box[2]
                 and inps.pix_box[1] <= inps.ref_yx[0] < inps.pix_box[3])):
        inps.disp_ref_pixel = False
        print('the native REF_Y/X is out of subset box, thus do not display')

    ## initial pixel coord
    if inps.lalo:
        inps.yx = inps.coord.geo2radar(inps.lalo[0], inps.lalo[1], print_msg=False)[0:2]
    if inps.yx:
        try:
            inps.lalo = inps.coord.radar2geo(inps.yx[0], inps.yx[1], print_msg=False)[0:2]
        except FileNotFoundError:
            inps.lalo = None

    ## figure settings
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

    ## display unit and wrap
    # if wrap_step == 2*np.pi (default value), set disp_unit_img = radian;
    # otherwise set disp_unit_img = disp_unit
    inps.disp_unit_img = inps.disp_unit
    if inps.wrap:
        inps.vlim = inps.wrap_range

        if (inps.wrap_range[1] - inps.wrap_range[0]) == 2*np.pi:
            inps.disp_unit_img = 'radian'

        if inps.disp_unit_img == 'radian':
            inps.range2phase = -4. * np.pi / float(atr['WAVELENGTH'])
            if   'cm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 100.
            elif 'mm' == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1000.
            elif 'm'  == inps.disp_unit.split('/')[0]:   inps.range2phase /= 1.
            else:
                raise ValueError(f'un-recognized display unit: {inps.disp_unit}')

    inps.cbar_label = 'Amplitude' if atr['DATA_TYPE'].startswith('complex') else 'Displacement'
    inps.cbar_label += f' [{inps.disp_unit_img}]'

    ## fit a suite of time func to the time series
    inps.model = time_func.inps2model(inps, date_list=inps.date_list, print_msg=inps.print_msg)

    # dense TS for plotting
    inps.date_list_fit = ptime.get_date_range(inps.date_list[0], inps.date_list[-1])
    inps.dates_fit = ptime.date_list2vector(inps.date_list_fit)[0]
    inps.G_fit = time_func.get_design_matrix4time_func(
        date_list=inps.date_list_fit,
        model=inps.model,
        seconds=inps.seconds)

    return inps, atr


def subset_and_multilook_yx(yx, pix_box=None, multilook_num=1):
    """Update row/col number of one pixel due to subset and multilooking."""
    y, x = yx
    if pix_box is not None:
        y -= pix_box[1]
        x -= pix_box[0]
    if multilook_num > 1:
        y = int((y - int(multilook_num / 2)) / multilook_num)
        x = int((x - int(multilook_num / 2)) / multilook_num)
    return (y, x)


def read_exclude_date(input_ex_date, dateListAll):
    """Read exclude list of dates
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
    ## read list of 3D time-series
    ts_data = []
    for fname in inps.file:
        msg = f'reading timeseries from file {fname}'
        msg += f' with step of {inps.multilook_num} by {inps.multilook_num}' if inps.multilook_num > 1 else ''
        vprint(msg)

        data, atr = readfile.read(
            fname,
            datasetName=inps.date_list,
            box=inps.pix_box,
            xstep=inps.multilook_num,
            ystep=inps.multilook_num)

        if atr['DATA_TYPE'].startswith('complex'):
            vprint('input data is complex, calculate its amplitude and continue')
            data = np.abs(data)

        if inps.ref_yx and inps.ref_yx != (int(atr.get('REF_Y', -1)), int(atr.get('REF_X', -1))):
            (ry, rx) = subset_and_multilook_yx(inps.ref_yx, inps.pix_box, inps.multilook_num)
            ref_phase = data[:, ry, rx]
            data -= np.tile(ref_phase.reshape(-1, 1, 1), (1, data.shape[-2], data.shape[-1]))
            vprint(f'reference to pixel: {inps.ref_yx}')

        if inps.ref_idx is not None:
            vprint(f'reference to date: {inps.date_list[inps.ref_idx]}')
            data -= np.tile(data[inps.ref_idx, :, :], (inps.num_date, 1, 1))

        # Display Unit
        data, inps.disp_unit, inps.unit_fac = pp.scale_data2disp_unit(
            data,
            metadata=atr,
            disp_unit=inps.disp_unit)
        ts_data.append(data)

    ## mask file: input mask file + non-zero ts pixels - ref_point
    mask = pp.read_mask(
        inps.file[0],
        mask_file=inps.mask_file,
        datasetName='displacement',
        box=inps.pix_box,
        xstep=inps.multilook_num,
        ystep=inps.multilook_num,
        print_msg=inps.print_msg,
    )[0]

    if mask is None:
        mask = np.ones(ts_data[0].shape[-2:], np.bool_)

    ts_stack = np.nansum(ts_data[0], axis=0)
    mask[np.isnan(ts_stack)] = False
    # keep all-zero value for unwrapError time-series
    if atr['UNIT'] not in ['cycle']:
        mask[ts_stack == 0.] = False
    del ts_stack

    # do not mask the reference point
    x0, y0, x1, y1 = inps.pix_box
    if inps.ref_yx and (x0 <= inps.ref_yx[1] < x1) and (y0 <= inps.ref_yx[0] < y1):
        (ry, rx) = subset_and_multilook_yx(inps.ref_yx, inps.pix_box, inps.multilook_num)
        mask[ry, rx] = True

    ## default vlim
    inps.dlim = [np.nanmin(ts_data[0]), np.nanmax(ts_data[0])]
    if not inps.vlim:
        inps.cmap_lut, inps.vlim = pp.auto_adjust_colormap_lut_and_disp_limit(
            ts_data[0], num_multilook=10, print_msg=inps.print_msg,
        )[:2]
    vprint(f'data    range: {inps.dlim} {inps.disp_unit}')
    vprint(f'display range: {inps.vlim} {inps.disp_unit}')

    ## default ylim
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

    kwargs = dict(
        fmt='-o', ms=ppar.ms, lw=0, alpha=1,
        elinewidth=inps.edge_width,
        capsize=ppar.ms*0.5, mew=inps.edge_width,
    )

    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        ax.errorbar(
            inps.ex_dates, ex_d_ts,
            yerr=inps.ex_error_ts,
            color='gray',
            ecolor='gray',
            **kwargs,
        )

    # Plot kept dates
    ax.errorbar(
        dates, d_ts,
        yerr=inps.error_ts,
        label=ppar.label,
        color=ppar.mfc,
        ecolor=ppar.mfc,
        **kwargs
    )

    handles = ax.get_legend_handles_labels()[0]

    return ax, handles[-1]


def plot_ts_scatter(ax, dis_ts, inps, ppar):
    # make a local copy
    dates = np.array(inps.dates)
    d_ts = dis_ts[:]

    kwargs = dict(ms=ppar.ms, marker=inps.marker)

    if inps.ex_date_list:
        # Update displacement time-series
        d_ts = dis_ts[inps.ex_flag == 1]
        dates = dates[inps.ex_flag == 1]

        # Plot excluded dates
        ex_d_ts = dis_ts[inps.ex_flag == 0]
        ax.plot(inps.ex_dates, ex_d_ts, color='gray', lw=0, **kwargs)

    # Plot kept dates
    handle, = ax.plot(dates, d_ts, color=ppar.mfc, label=ppar.label, lw=inps.linewidth, **kwargs)
    return ax, handle


def plot_ts_fit(ax, ts_fit, inps, ppar, m_strs=None, ts_fit_lim=None):
    """Plot time series fitting results."""
    # plot the model prediction uncertainty boundaries
    h0 = None
    if ts_fit_lim is not None:
        h0 = ax.fill_between(
            inps.dates_fit, ts_fit_lim[0], ts_fit_lim[1],
            fc=ppar.color, ec='none', alpha=0.2,
        )

    # plot the model prediction curve in a fine date_lists
    h1, = ax.plot(inps.dates_fit, ts_fit, color=ppar.color, lw=ppar.linewidth, alpha=0.8)

    # legend
    handles = []
    labels = []
    kwargs = dict(loc='best', fontsize=ppar.fontsize)
    if h0 is not None:
        handles.append((h0, h1))
        labels.append('time func. pred. (w. 95% conf. interval.)')

    # print model parameters on the plot
    # link for auto text loc: https://stackoverflow.com/questions/7045729
    if m_strs:
        if len(labels) == 0:
            kwargs['handlelength'] = 0
            kwargs['handletextpad'] = 0
            kwargs['borderaxespad'] = 1.0

        # remove multiple spaces for better display since matplotlib does not give good alignment as in the terminal
        handles += [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * len(m_strs)
        labels += [re.sub(' +', ' ', x) for x in m_strs]

    # do not print model parameters for multi files [temporarily]
    # solutions:
    # 1. plot multiple legends, need to find out how to layout them, otherwise they would overlap on top of each other
    #    link: https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#multiple-legends-on-the-same-axes
    # 2. update the same legend, by saving the handles and labels for all files
    if len(labels) > 0 and len(inps.file) == 1:
        ax.legend(handles, labels, **kwargs)

    return ax


def get_model_param_str(model, ds_dict, disp_unit='cm'):
    """Summary model parameters in a str paragraph.
    Parameters: model     - dict, model dictionary
                ds_dict   - dict, est. time func. param
                disp_unit - float, display unit for length, which can be scaled
    Returns:    ds_strs   - list of strings, each to summary the result of one time func
    """

    # dataset unit dict
    ds_unit_dict = ts2vel.model2hdf5_dataset(model)[2]
    ds_names = list(ds_unit_dict.keys())

    # update ds_unit_dict based on disp_unit
    for ds_name in ds_names:
        units = ds_unit_dict[ds_name].split('/')
        if units[0] == 'm' and disp_unit != 'm':
            units[0] = disp_unit
            ds_unit_dict[ds_name] = '/'.join(units)

    # list of dataset names
    ds_names = [x for x in ds_dict.keys() if not x.endswith('Std') and x not in ['intercept']]
    w_key = max(len(x) for x in ds_names)
    w_val = max(len(f'{x[0]:.2f}') for x in ds_dict.values())

    ds_strs = []
    for ds_name in ds_names:
        # get param info
        ds_value = ds_dict[ds_name]
        ds_unit = ds_unit_dict[ds_name]
        ds_std_value = ds_dict.get(ds_name+'Std', None)

        # compose string
        ds_str = f'{ds_name:<{w_key}}: {ds_value[0]:>{w_val}.2f}'
        ds_str += f' +/- {ds_std_value[0]:>{w_val}.2f}' if ds_std_value is not None else ''
        ds_str += f' {ds_unit}'
        ds_strs.append(ds_str)

    return ds_strs


def fit_time_func(model, date_list, ts_dis, disp_unit='cm', G_fit=None, conf_level=0.95, seconds=0):
    """Fit a suite of time functions to the time series.
    Equations:  Gm = d
    Parameters: model      - dict of time functions, check utils.time_func.estimate_time_func() for details.
                date_list  - list of dates in YYYYMMDD format
                ts_dis     - 1D np.ndarray, displacement time series
                disp_unit  - str, display unit for length, which can be scaled
                G_fit      - 2D np.ndarray, design matrix for the dense time series prediction plot
                conf_level - float in [0,1], confidence level of the plotted confidence intervals
    Returns:    m_strs     - dict, dictionary in {ds_name: ds_value}
                ts_fit     - 1D np.ndarray, dense time series fit for plotting
                ts_fit_lim - list of 1D np.ndarray, the lower and upper
                             boundaries of dense time series fit for plotting
    """
    # init output
    m_strs = []
    ts_fit = None
    ts_fit_lim = None

    if np.all(np.isnan(ts_dis)):
        return m_strs, ts_fit, ts_fit_lim

    # 1.1 estimate time func parameter via least squares (OLS)
    G, m, e2 = time_func.estimate_time_func(
        model=model,
        date_list=date_list,
        dis_ts=ts_dis,
        seconds=seconds)

    # 1.2 calc the precision of time func parameters
    # using the OLS estimation residues e2 = sum((d - Gm) ** 2)
    # assuming obs errors following normal distribution in time
    num_obs = len(date_list)
    num_param = G.shape[1]
    G_inv = linalg.inv(np.dot(G.T, G))
    m_var_sum = e2.flatten() / (num_obs - num_param)
    m_std = np.sqrt(np.dot(np.diag(G_inv).reshape(-1, 1), m_var_sum))

    # 1.3 translate estimation result into HDF5 ready datasets
    # AND compose list of strings for printout
    m_dict = ts2vel.model2hdf5_dataset(model, m, m_std)[0]
    m_strs = get_model_param_str(model, m_dict, disp_unit=disp_unit)

    # 2. reconstruct the fine resolution function
    if G_fit is not None:
        ts_fit = np.matmul(G_fit, m)
        ts_fit_std = np.sqrt(np.diag(G_fit.dot(np.diag(m_std**2)).dot(G_fit.T)))

        # calc confidence interval
        # references:
        # 1. Exercise 6.4 OMT: Interpretation from Hanssen et al. (2017) EdX online course.
        #    Hanssen, R., Verhagen, S. and Samiei-Esfahany, S., (2017) Observation Theory: Estimating the Unknown,
        #    Available at: https://www.edx.org/course/observation-theory-estimating-the-unknown
        # 2. https://stackoverflow.com/questions/20626994
        alpha = 1 - conf_level                                # level of significance
        conf_int_scale = stats.norm.ppf(1 - alpha / 2)        # scaling factor for confidence interval
        ts_fit_lim = [ts_fit - conf_int_scale * ts_fit_std,
                      ts_fit + conf_int_scale * ts_fit_std]

    return m_strs, ts_fit, ts_fit_lim


def get_point_coord_str(y, x, coord_obj, lalo_digit=5):
    """Get the string of the point coordinates.

    Parameters: y / x      - int, row / column number
                coord_obj  - mintpy.objects.coordinate object
                lalo_digit - int, digit of the decimal place for lat/lon
    Returns:    pts_str    - str, point coordinate string
    """
    coord_str = f'Y/X = {y}, {x}'
    try:
        lat, lon = coord_obj.radar2geo(y, x, print_msg=False)[0:2]
        coord_str += f', lat/lon = {lat:.{lalo_digit}f}, {lon:.{lalo_digit}f}'
    except FileNotFoundError:
        pass
    return coord_str


def save_ts_data_and_plot(yx, d_ts, m_strs, inps):
    """Save TS data and plots into files."""
    y, x = yx
    vprint(f'save info on pixel ({y}, {x})')

    # output file name
    if inps.outfile:
        inps.outfile_base, fext = os.path.splitext(inps.outfile[0])
        if fext != '.pdf':
            msg = 'Output file extension is fixed to .pdf,'
            msg += f' input extension {fext} is ignored.'
            vprint(msg)
    else:
        inps.outfile_base = f'y{y}x{x}'

    # TXT - point time-series and time func param
    outName = f'{inps.outfile_base}_ts.txt'
    header = f'time-series file = {inps.file[0]}\n'
    header += f'{get_point_coord_str(y, x, inps.coord, inps.lalo_digit)}\n'
    header += f'reference pixel: y={inps.ref_yx[0]}, x={inps.ref_yx[1]}\n' if inps.ref_yx else ''
    header += f'reference date: {inps.date_list[inps.ref_idx]}\n' if inps.ref_idx else ''
    header += 'estimated time function parameters:\n'
    for m_str in m_strs:
        header += f'    {m_str}\n'
    header += f'unit: {inps.disp_unit}'

    # prepare data
    data = np.hstack((np.array(inps.date_list).reshape(-1, 1), d_ts.reshape(-1, 1)))

    # write
    np.savetxt(outName, data, fmt='%s', delimiter='\t', header=header)
    print('save displacement time-series to file: '+outName)

    # Figure - point time-series
    outName = f'{inps.outfile_base}_ts.pdf'
    inps.fig_pts.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    print('save time-series plot to file: '+outName)

    # Figure - map
    outName = f'{inps.outfile_base}_{inps.date_list[inps.idx]}.png'
    inps.fig_img.savefig(outName, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    print('save map plot to file: '+outName)
    return


class timeseriesViewer():
    """Class for tsview.py

    Example:
        from mintpy.cli.tsview import cmd_line_parse
        from mintpy.tsview import timeseriesViewer

        cmd = 'timeseries.h5 --yx 273 271 --figsize 8 4'
        inps = cmd_line_parse(cmd.split())
        obj = timeseriesViewer(inps)
        obj.open()
        obj.plot()
    """

    def __init__(self, inps):

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

        # copy inps to self object
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

    def open(self):
        global vprint
        vprint = print if self.print_msg else lambda *args, **kwargs: None

        # print command line
        if self.argv is not None:
            print(f'{os.path.basename(__file__)} ' + ' '.join(self.argv))

        # matplotlib backend setting
        if not self.disp_fig:
            plt.switch_backend('Agg')

        self, self.atr = read_init_info(self)

        # input figsize for the image/point time-series plot
        self.figsize_img = self.fig_size_img
        self.figsize_pts = self.fig_size
        self.pts_marker = 'r^'
        self.pts_marker_size = 6.

    def plot(self):
        # read 3D time-series
        self.ts_data, self.mask = read_timeseries_data(self)[0:2]

        # Figure 1 - Cumulative Displacement Map
        if not self.figsize_img:
            if self.geo_box and self.fig_coord == 'geo':
                w, n, e, s = self.geo_box
                ds_shape = (e - w, n - s)
            else:
                ds_shape = self.ts_data[0].shape[-2:]
            self.figsize_img = pp.auto_figure_size(
                ds_shape=ds_shape,
                disp_cbar=True,
                disp_slider=True,
                print_msg=False)
        vprint(f'create figure for map in size of [{self.figsize_img[0]:.1f}, {self.figsize_img[1]:.1f}]')
        subplot_kw = dict(projection=self.map_proj_obj) if self.map_proj_obj is not None else {}
        self.fig_img, self.ax_img = plt.subplots(figsize=self.figsize_img, subplot_kw=subplot_kw)

        # Figure 1 - Axes 1 - Displacement Map
        img_data = np.array(self.ts_data[0][self.idx, :, :])
        img_data[self.mask == 0] = np.nan
        self.plot_init_image(img_data)

        # Figure 1 - Axes 2 - Time Slider
        self.plot_init_time_slider(init_idx=self.idx, ref_idx=self.ref_idx)
        self.tslider.on_changed(self.update_time_slider)

        # Figure 2 - Time Series Displacement - Point
        vprint(f'create figure for point in size of [{self.figsize_pts[0]:.1f}, {self.figsize_pts[1]:.1f}]')
        self.fig_pts, self.ax_pts = plt.subplots(num=self.figname_pts, figsize=self.figsize_pts)
        if self.yx:
            d_ts, m_strs = self.plot_point_timeseries(self.yx)

            # save figures and data to files
            if self.save_fig:
                save_ts_data_and_plot(self.yx, d_ts, m_strs, self)

        # Final linking of the canvas to the plots.
        self.fig_img.canvas.mpl_connect('button_press_event', self.update_point_timeseries)
        self.fig_img.canvas.mpl_connect('key_press_event', self.update_image)
        if self.disp_fig:
            vprint('showing ...')
            msg = '\n------------------------------------------------------------------------'
            msg += '\nTo scroll through the image sequence:'
            msg += '\n1) Move the slider, OR'
            msg += '\n2) Press left or right arrow key (if not responding, click the image and try again).'
            msg += '\n------------------------------------------------------------------------'
            vprint(msg)

            # --no-show-map option
            # requires --yx/lalo input
            if self.yx and not self.disp_fig_img:
                plt.close(self.fig_img)

            plt.show()
        return


    ##---------- event functions
    def update_point_timeseries(self, event):
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


    def update_image(self, event):
        """Slide images with left/right key on keyboard"""
        if event.inaxes and event.inaxes.figure == self.fig_img:
            idx = None
            if event.key == 'left':
                idx = max(self.idx - 1, 0)
            elif event.key == 'right':
                idx = min(self.idx + 1, self.num_date - 1)

            if idx is not None and idx != self.idx:
                # update title
                disp_date = self.dates[idx].strftime(self.disp_date_format)
                sub_title = f'N = {idx}, Time = {disp_date}'
                self.ax_img.set_title(sub_title, fontsize=self.font_size)

                # read data
                data_img = np.array(self.ts_data[0][idx, :, :])
                data_img[self.mask == 0] = np.nan
                if self.wrap:
                    if self.disp_unit_img == 'radian':
                        data_img *= self.range2phase
                    data_img = ut.wrap(data_img, wrap_range=self.wrap_range)

                # update
                self.tslider.set_val(idx)         # update slider
                self.img.set_data(data_img)       # update image
                self.idx = idx
                self.fig_img.canvas.draw_idle()
                self.fig_img.canvas.flush_events()
        return


    def update_time_slider(self, val):
        """Update Displacement Map using Slider"""
        self.idx = self.tslider.val

        # update title
        disp_date = self.dates[self.idx].strftime(self.disp_date_format)
        sub_title = f'N = {self.idx}, Time = {disp_date}'
        self.ax_img.set_title(sub_title, fontsize=self.font_size)

        # read/update 2D image data
        data_img = self.ts_data[0][self.idx, :, :]
        data_img[self.mask == 0] = np.nan
        if self.wrap:
            if self.disp_unit_img == 'radian':
                data_img *= self.range2phase
            data_img = ut.wrap(data_img, wrap_range=self.wrap_range)
        self.img.set_data(data_img)

        # update figure
        self.fig_img.canvas.draw_idle()
        self.fig_img.canvas.flush_events()
        return


    ##---------- plot functions
    def plot_init_image(self, img_data):
        """Plot the initial 2D image."""
        # prepare data
        if self.wrap:
            if self.disp_unit_img == 'radian':
                img_data *= self.range2phase
            img_data = ut.wrap(img_data, wrap_range=self.wrap_range)

        # Title and Axis Label
        self.disp_date_format = ptime.get_compact_isoformat(self.date_list[0])
        disp_date = self.dates[self.idx].strftime(self.disp_date_format)
        self.fig_title = f'N = {self.idx}, Time = {disp_date}'

        # Initial Pixel of interest
        self.pts_yx = None
        self.pts_lalo = None
        if self.yx and self.yx != self.ref_yx:
            self.pts_yx = np.array(self.yx).reshape(-1, 2)
            if self.lalo:
                self.pts_lalo = np.array(self.lalo).reshape(-1, 2)

        # call view.py to plot
        self.img, self.cbar_img = view.plot_slice(self.ax_img, img_data, self.atr, self)[2:4]
        self.fig_img.canvas.manager.set_window_title(self.figname_img)
        self.fig_img.tight_layout(rect=(0, 0.16, 1, 0.97))

        return self.img, self.cbar_img


    def plot_init_time_slider(self, init_idx=-1, ref_idx=None):
        """Plot the initial slider."""
        # initiate axes
        #self.fig_img.subplots_adjust(bottom=0.16)
        self.ax_tslider = self.fig_img.add_axes([0.125, 0.05, 0.75, 0.03])

        # plot slider
        self.tslider = widgets.Slider(
            ax=self.ax_tslider,
            label='Image',
            valinit=init_idx,
            valmin=0,
            valmax=self.num_date-1,
            valstep=1)

        # plot reference date:
        # as a gray dot on the slider AND
        # as x-axis label
        if ref_idx is not None:
            self.tslider.ax.scatter(ref_idx, 0.5, s=8**2, marker='o', color='gray', edgecolors='w')
            disp_date = self.dates[ref_idx].strftime(self.disp_date_format)
            self.ax_tslider.set_title(f'Reference: N = {ref_idx}, Time = {disp_date}', fontsize=self.font_size)

        return self.tslider


    def plot_point_timeseries(self, yx):
        """Plot point displacement time-series at pixel [y, x]
        Parameters: yx     : list of 2 int
        Returns:    ts_dis : 1D np.array in size of (num_date,) for the 1st file
                    m_strs : list of strings for the est. time func. param. for the 1st file
        """
        ax = self.ax_pts
        ax.cla()

        # plot scatter in different size for different files
        num_file = len(self.ts_data)
        if   num_file <= 2: ms_step = 4
        elif num_file == 3: ms_step = 3
        elif num_file == 4: ms_step = 2
        elif num_file >= 5: ms_step = 1

        # get local Y/X coord for the subsetted and multilooked 3D data cube
        (y, x) = subset_and_multilook_yx(yx, self.pix_box, self.multilook_num)

        handles, labels = [], []
        for i in range(num_file-1, -1, -1):
            # get displacement data
            ts_dis = self.ts_data[i][:, y, x]

            # fit time func
            m_strs, ts_fit, ts_fit_lim = fit_time_func(
                model=self.model,
                date_list=np.array(self.date_list)[self.ex_flag].tolist(),
                ts_dis=ts_dis[self.ex_flag],
                disp_unit=self.disp_unit,
                G_fit=self.G_fit,
                seconds=self.seconds)

            if self.zero_first:
                off = ts_dis[self.zero_idx]
                ts_dis -= off
                ts_fit -= off

            if self.offset:
                ts_dis += self.offset * (num_file - 1 - i)
                ts_fit += self.offset * (num_file - 1 - i)

            # plot
            if not np.all(np.isnan(ts_dis)):
                ppar = argparse.Namespace()
                ppar.label = self.file_label[i]
                ppar.mfc = f'C{num_file-1-i}' if self.mask[y, x] != 0 else 'gray'
                ppar.ms = self.marker_size - ms_step * (num_file - 1 - i)
                # use smaller marker size for very long time series
                ppar.ms /= 10 if self.num_date > 1e3 else 1

                handle = self.ts_plot_func(ax, ts_dis, self, ppar)[1]
                handles.append(handle)
                labels.append(ppar.label)

                # plot model prediction
                if self.plot_model:
                    fpar = argparse.Namespace()
                    fpar.linewidth = 3
                    fpar.color = 'C1' if num_file == 1 else ppar.mfc
                    fpar.fontsize = self.font_size

                    if not self.plot_model_conf_int:
                        ts_fit_lim = None

                    plot_ts_fit(
                        ax, ts_fit, self, fpar,
                        m_strs=m_strs,
                        ts_fit_lim=ts_fit_lim,
                    )

        # axis format
        ax.tick_params(which='both', direction='in', labelsize=self.font_size,
                       bottom=True, top=True, left=True, right=True)
        pp.auto_adjust_xaxis_date(ax, self.yearList, fontsize=self.font_size)
        ax.set_ylabel(self.cbar_label, fontsize=self.font_size)
        ax.set_ylim(self.ylim)
        if self.tick_right:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")

        # title
        title = get_point_coord_str(yx[0], yx[1], self.coord, self.lalo_digit)
        title += ' (masked out)' if self.mask[y, x] == 0 else ''
        if self.disp_title:
            ax.set_title(title, fontsize=self.font_size)

        # legend
        if len(self.ts_data) > 1:
            ax.legend(handles, labels)

        # Print to terminal
        vprint('\n---------------------------------------')
        vprint(title)
        float_formatter = lambda x: [float(f'{i:.2f}') for i in x]
        if self.num_date <= 1e3:
            vprint(float_formatter(ts_dis))

        if not np.all(np.isnan(ts_dis)):
            # min/max displacement
            ts_min, ts_max = np.nanmin(ts_dis), np.nanmax(ts_dis)
            vprint(f'time-series range: [{ts_min:.2f}, {ts_max:.2f}] {self.disp_unit}')

            # time func param
            vprint('time function parameters:')
            for m_str in m_strs:
                vprint(f'    {m_str}')

            # update figure
            # use fig.canvas.draw_idel() instead of fig.canvas.draw()
            # reference: https://stackoverflow.com/questions/64789437
            self.fig_pts.canvas.draw_idle()
            self.fig_pts.canvas.flush_events()

        return ts_dis, m_strs
