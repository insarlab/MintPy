"""Utilities wrapped around cartopy/matplotlib for maps."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Joshua Zahner, Jun 2022            #
############################################################
# Recommend import:
#   from mintpy.utils.map import draw_lalo_label, draw_scalebar
#   from mintpy.utils import plot as pp


import numpy as np
import pyproj
from cartopy import crs as ccrs
from cartopy.mpl import ticker as cticker
from matplotlib import pyplot as plt

##########################################  Lat/Lon Labels  ##########################################

def draw_lalo_label(ax, geo_box, lalo_step=None, lalo_loc=[1, 0, 0, 1], lalo_max_num=4, lalo_offset=None,
                    projection=ccrs.PlateCarree(), font_size=None, print_msg=True):
    """Auto draw lat/lon label/tick based on coverage from geo_box
    Parameters: ax           : cartopy axes.
                geo_box      : 4-tuple of float, (W, N, E, S) in degree
                lalo_step    : float, major tick interval for X/Y axes
                lalo_loc     : list of 4 bool, positions where the labels are drawn as in [left, right, top, bottom]
                               default: [1,0,0,1]
                lalo_max_num : int, maximum number of major ticks for X/Y axes
                lalo_offset  : float, distance in points between tick and label.
                               Set to negative value to move the ticklabel inside the plot.
                projection   : cartopy.crs object
                ...
    Example:    geo_box = (128.0, 37.0, 138.0, 30.0)
                m.draw_lalo_label(geo_box)
    """

    # default lat/lon sequences
    lats, lons, lalo_step, digit = auto_lalo_sequence(geo_box,
                                                      lalo_step=lalo_step,
                                                      lalo_max_num=lalo_max_num)
    if print_msg:
        print(f'plot lat/lon label in step of {lalo_step} and location of {lalo_loc}')

    # ticklabel/tick style
    ax.tick_params(which='both', direction='in', labelsize=font_size,
                   left=True, right=True, top=True, bottom=True,
                   labelleft=lalo_loc[0], labelright=lalo_loc[1],
                   labeltop=lalo_loc[2], labelbottom=lalo_loc[3])

    # custom padding to move the tick labels inside axis
    if lalo_offset:
        ax.tick_params(axis='x', which='major', pad=lalo_offset[1])
        ax.tick_params(axis='y', which='major', pad=lalo_offset[0])

    # ticklabel symbol style
    dec_digit = max(0, 0-digit)
    lon_formatter = cticker.LongitudeFormatter(number_format=f'.{dec_digit}f')
    lat_formatter = cticker.LatitudeFormatter(number_format=f'.{dec_digit}f')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # plot ticks & tick labels
    ax.set_xticks(lons, crs=projection)
    ax.set_yticks(lats, crs=projection)

    return ax


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

    max_lalo_dist = max([geo_box[1] - geo_box[3], geo_box[2] - geo_box[0]])

    if not lalo_step:
        # Initial tick step
        lalo_step = round_to_1(max_lalo_dist / lalo_max_num)

        # reduce decimal if it ends with 8/9
        digit = int(np.floor(np.log10(lalo_step)))
        if str(lalo_step)[-1] in ['8','9']:
            digit += 1
            lalo_step = round(lalo_step, digit)

        # Final tick step - choose from candidate list
        lalo_step_candidate = [i*10**digit for i in step_candidate]
        distance = [(i - max_lalo_dist / lalo_max_num) ** 2 for i in lalo_step_candidate]
        lalo_step = lalo_step_candidate[distance.index(min(distance))]

    digit = int(np.floor(np.log10(lalo_step)))

    # Auto tick sequence
    lat_major = np.ceil( geo_box[3] / 10**(digit+1) ) * 10**(digit+1)
    lon_major = np.ceil( geo_box[0] / 10**(digit+1) ) * 10**(digit+1)
    lats = np.hstack((np.arange(lat_major, lat_major - 10.*max_lalo_dist, -lalo_step),
                      np.arange(lat_major, lat_major + 10.*max_lalo_dist,  lalo_step)))
    lons = np.hstack((np.arange(lon_major, lon_major - 10.*max_lalo_dist, -lalo_step),
                      np.arange(lon_major, lon_major + 10.*max_lalo_dist,  lalo_step)))
    lats = np.sort(np.unique(lats))
    lons = np.sort(np.unique(lons))
    lats = lats[np.where(np.logical_and(lats >= geo_box[3], lats <= geo_box[1]))]
    lons = lons[np.where(np.logical_and(lons >= geo_box[0], lons <= geo_box[2]))]

    return lats, lons, lalo_step, digit


############################################  Scale Bar  #############################################

def draw_scalebar(ax, geo_box, unit='degrees', loc=[0.2, 0.2, 0.1], labelpad=0.05, font_size=12,
                  color='k', linewidth=2):
    """draw a simple map scale from x1,y to x2,y in map projection coordinates, label it with actual distance
    ref_link: http://matplotlib.1069221.n5.nabble.com/basemap-scalebar-td14133.html
    Parameters: ax       : matplotlib.pyplot.axes object
                geo_box  : tuple of 4 float in (x0, y0, x1, y1) for (W, N, E, S) in degrees / meters
                unit     : str, coordinate unit - degrees or meters
                loc      : list of 3 float in (length, y, x) of scale bar location:
                           length = ratio of the total width
                           y / x = axes fraction of the scale bar center
                labelpad : float
    Returns:    ax       : matplotlib.pyplot.axes object
    Example:    from mintpy.utils import plot as pp
                pp.draw_scale_bar(ax, geo_box)
    """
    geod = pyproj.Geod(ellps='WGS84')
    ax = ax if ax else plt.gca()

    ## location - center
    lon_c = geo_box[0] + loc[1] * (geo_box[2] - geo_box[0])
    lat_c = geo_box[3] + loc[2] * (geo_box[1] - geo_box[3])

    ## length
    # 1. calc scene width in meters
    if unit.startswith('meter'):
        scene_width = geo_box[2] - geo_box[0]
    elif unit.startswith('deg'):
        scene_width = geod.inv(geo_box[0], geo_box[3],
                               geo_box[2], geo_box[3])[2]

        # do not plot scalebar if the longitude span > 30 deg
        if (geo_box[2] - geo_box[0]) > 30:
            return ax

    # 2. convert length ratio to length in meters
    length_meter = round_to_1(scene_width * loc[0])
    # round to the nearest km
    if length_meter > 1000.0:
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
    kwargs = dict(color=color, linewidth=linewidth)
    ax.plot([lon0, lon1], [lat_c, lat_c], **kwargs)
    ax.plot([lon0, lon0], [lat_c, lat_c + 0.1*length_disp], **kwargs)
    ax.plot([lon1, lon1], [lat_c, lat_c + 0.1*length_disp], **kwargs)

    ## plot scale bar label
    unit = 'm'
    if length_meter >= 1000.0:
        unit = 'km'
        length_meter *= 0.001
    txt_offset = (geo_box[1] - geo_box[3]) * labelpad

    ax.text(x=lon0+length_disp/2.0,
            y=lat_c+txt_offset,
            s=f'{length_meter:.0f} {unit}',
            va='center',
            ha='center',
            fontsize=font_size,
            color=color)

    return ax


##############################################  Utils  ###############################################
# copied from utils.util0.py to simplify the module dependency

def round_to_1(x):
    """Return the most significant digit of input number"""
    digit = int(np.floor(np.log10(abs(x))))
    return round(x, -digit)
