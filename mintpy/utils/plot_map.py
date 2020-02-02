############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, Zhang Yunjun, 2020                #
############################################################
# Recommend import:
#     from mintpy.utils import plot as pp OR
#     from mintpy.utils.plot_map import draw_lalo_label, draw_scalebar


import pyproj
import numpy as np
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl import ticker
from mintpy.utils import utils0 as ut0


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

    # Auto tick sequence
    digit = np.int(np.floor(np.log10(lalo_step)))
    lat_major = np.ceil(geo_box[3]/10**(digit+1))*10**(digit+1)
    lats = np.unique(np.hstack((np.arange(lat_major, lat_major-10.*max_lalo_dist, -lalo_step),
                                np.arange(lat_major, lat_major+10.*max_lalo_dist, lalo_step))))
    lats = np.sort(lats[np.where(np.logical_and(lats >= geo_box[3], lats <= geo_box[1]))])

    lon_major = np.ceil(geo_box[0]/10**(digit+1))*10**(digit+1)
    lons = np.unique(np.hstack((np.arange(lon_major, lon_major-10.*max_lalo_dist, -lalo_step),
                                np.arange(lon_major, lon_major+10.*max_lalo_dist, lalo_step))))
    lons = np.sort(lons[np.where(np.logical_and(lons >= geo_box[0], lons <= geo_box[2]))])

    ## Need to artificially add the geo_box boundaries in order to get cartopy to draw gridlines
    ## all the way to the map edge. Can be removed if drawing of grid lines is not needed
    #lons = np.insert(lons, [0, -1], [geo_box[0], geo_box[2]])
    #lats = np.insert(lats, [0, -1], [geo_box[3], geo_box[1]])
    return lats, lons, lalo_step, digit


def draw_lalo_label(geo_box, ax=None, lalo_step=None, lalo_loc=[1, 0, 0, 1], lalo_max_num=4,
                    font_size=12, xoffset=None, yoffset=None, yrotate='horizontal',
                    projection=ccrs.PlateCarree(), print_msg=True):
    """Auto draw lat/lon label/tick based on coverage from geo_box
    Parameters: geo_box   : 4-tuple of float, defining UL_lon, UL_lat, LR_lon, LR_lat coordinate
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
    lon_formatter = ticker.LongitudeFormatter(number_format='.{}f'.format(decimal_digit))
    lat_formatter = ticker.LatitudeFormatter(number_format='.{}f'.format(decimal_digit))
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.set_xticks(lons, crs=projection)
    ax.set_yticks(lats, crs=projection)
    return ax


def draw_scalebar(ax, geo_box, loc=[0.2, 0.2, 0.1], labelpad=0.05, font_size=12, color='k'):
    """draw a simple map scale from x1,y to x2,y in map projection coordinates, label it with actual distance
    ref_link: http://matplotlib.1069221.n5.nabble.com/basemap-scalebar-td14133.html
    Parameters: ax       : matplotlib.pyplot.axes object
                geo_box  : tuple of 4 float in (x0, y0, x1, y1) for (W, N, E, S) in degrees
                loc      : list of 3 float, distance, lat/lon of scale bar center in ratio of width, relative coord
                labelpad : float
    Returns:    ax
    Example:    from mintpy.utils import plot as pp
                pp.draw_scale_bar(ax, geo_box)
    """
    if not ax:
        ax = plt.gca()

    geod = pyproj.Geod(ellps='WGS84')

    # length in meter
    scene_width = geod.inv(geo_box[0], geo_box[3], geo_box[2], geo_box[3])[2]
    distance = ut0.round_to_1(scene_width * loc[0])
    lon_c = geo_box[0] + loc[1] * (geo_box[2] - geo_box[0])
    lat_c = geo_box[3] + loc[2] * (geo_box[1] - geo_box[3])

    # plot scale bar
    if distance > 1000.0:
        distance = np.rint(distance/1000.0)*1000.0
    lon_c2, lat_c2 = geod.fwd(lon_c, lat_c, 90, distance)[0:2]
    length = np.abs(lon_c - lon_c2)
    lon0 = lon_c - length/2.0
    lon1 = lon_c + length/2.0

    ax.plot([lon0, lon1], [lat_c, lat_c], color=color)
    ax.plot([lon0, lon0], [lat_c, lat_c + 0.1*length], color=color)
    ax.plot([lon1, lon1], [lat_c, lat_c + 0.1*length], color=color)

    # plot scale bar label
    unit = 'm'
    if distance >= 1000.0:
        unit = 'km'
        distance *= 0.001
    label = '{:.0f} {}'.format(distance, unit)
    txt_offset = (geo_box[1] - geo_box[3]) * labelpad

    ax.text(lon0+0.5*length, lat_c+txt_offset, label,
            verticalalignment='center', horizontalalignment='center',
            fontsize=font_size, color=color)
    return ax


def draw_scalebar_cartopy(ax, img_extent, location=[0.1, 0.1], length=0.2):
    """Draws scalebar on plot at location with a input proportional length
    Parameters: ax          : CartoPy axes.
                img_extent  : 4-tuple of floats, img extent in UL lon, UL lat, LR lon, LR lat
                location    : list of float, location on plot to draw lower left edge of scalebar
                length      : int, proportional length to draw scalebar
    Returns:    scalebar    : CartopyScalebar
    Example:    scalebar = pp.draw_scalebar(ax, img_extent, location=[0.2, 0.2], length=0.5)
    """
    scalebar = CartopyScalebar()

    # Compute scalebar length as a proportion of dataset width (km)
    geod = pyproj.Geod(ellps='WGS84')
    dist = geod.inv(img_extent[0], img_extent[2], img_extent[1], img_extent[3])[2] / 1000  # convert to kilometers
    length = length * dist  # scale by length proportion

    scalebar.draw(ax=ax, location=location, length=length)
    return scalebar



######################################### Cartopy Scalebar class begin ####################################
class CartopyScalebar:
    """From: https://stackoverflow.com/a/50674451/2676166"""

    def _axes_to_lonlat(self, ax, coords):
        """(lon, lat) from axes coordinates."""
        display = ax.transAxes.transform(coords)
        data = ax.transData.inverted().transform(display)
        lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)
        return lonlat


    def _upper_bound(self, start, direction, distance, dist_func):
        """A point farther than distance from start, in the given direction.

        It doesn't matter which coordinate system start is given in, as long
        as dist_func takes points in that coordinate system.

        Args:
            start:     Starting point for the line.
            direction  Nonzero (2, 1)-shaped array, a direction vector.
            distance:  Positive distance to go past.
            dist_func: A two-argument function which returns distance.

        Returns:
            Coordinates of a point (a (2, 1)-shaped NumPy array).
        """
        if distance <= 0:
            raise ValueError("Minimum distance is not positive: {}".format(distance))

        if np.linalg.norm(direction) == 0:
            raise ValueError("Direction vector must not be zero.")

        # Exponential search until the distance between start and end is
        # greater than the given limit.
        length = 0.1
        end = start + length * direction

        while dist_func(start, end) < distance:
            length *= 2
            end = start + length * direction

        return end


    def _distance_along_line(self, start, end, distance, dist_func, tol):
        """Point at a distance from start on the segment  from start to end.

        It doesn't matter which coordinate system start is given in, as long
        as dist_func takes points in that coordinate system.

        Args:
            start:     Starting point for the line.
            end:       Outer bound on point's location.
            distance:  Positive distance to travel.
            dist_func: Two-argument function which returns distance.
            tol:       Relative error in distance to allow.

        Returns:
            Coordinates of a point (a (2, 1)-shaped NumPy array).
        """
        initial_distance = dist_func(start, end)
        if initial_distance < distance:
            raise ValueError("End is closer to start ({}) than given distance ({}).".format(
                initial_distance, distance))

        if tol <= 0:
            raise ValueError("Tolerance is not positive: {}".format(tol))

        # Binary search for a point at the given distance.
        left = start
        right = end

        while not np.isclose(dist_func(start, right), distance, rtol=tol):
            midpoint = (left + right) / 2

            # If midpoint is too close, search in second half.
            if dist_func(start, midpoint) < distance:
                left = midpoint
            # Otherwise the midpoint is too far, so search in first half.
            else:
                right = midpoint

        return right


    def _point_along_line(self, ax, start, distance, angle=0, tol=0.01):
        """Point at a given distance from start at a given angle.

        Args:
            ax:       CartoPy axes.
            start:    Starting point for the line in axes coordinates.
            distance: Positive physical distance to travel.
            angle:    Anti-clockwise angle for the bar, in radians. Default: 0
            tol:      Relative error in distance to allow. Default: 0.01

        Returns:
            Coordinates of a point (a (2, 1)-shaped NumPy array).
        """
        from cartopy import geodesic as cgeo
        # Physical distance between points.
        def dist_func(a_axes, b_axes):
            a_phys = self._axes_to_lonlat(ax, a_axes)
            b_phys = self._axes_to_lonlat(ax, b_axes)
            dist = cgeo.Geodesic().inverse(a_phys, b_phys).base[0, 0]
            return dist

        # Direction vector of the line in axes coordinates.
        direction = np.array([np.cos(angle), np.sin(angle)])
        end = self._upper_bound(start, direction, distance, dist_func)
        return self._distance_along_line(start, end, distance, dist_func, tol)


    def draw(self, ax, location, length=0.2, metres_per_unit=1000, unit_name='km',
                  tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
                  ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
                  **kwargs):
        """Add a scale bar to CartoPy axes.

        For angles between 0 and 90 the text and line may be plotted at
        slightly different angles for unknown reasons. To work around this,
        override the 'rotation' keyword argument with text_kwargs.

        Args:
            ax:              CartoPy axes.
            location:        Position of left-side of bar in axes coordinates.
            length:          Geodesic length of the scale bar.
            metres_per_unit: Number of metres in the given unit. Default: 1000
            unit_name:       Name of the given unit. Default: 'km'
            tol:             Allowed relative error in length of bar. Default: 0.01
            angle:           Anti-clockwise rotation of the bar.
            color:           Color of the bar and text. Default: 'black'
            linewidth:       Same argument as for plot.
            text_offset:     Perpendicular offset for text in axes coordinates.
                             Default: 0.005
            ha:              Horizontal alignment. Default: 'center'
            va:              Vertical alignment. Default: 'bottom'
            **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
            **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
            **kwargs:        Keyword arguments for both plot and text.
        """
        # Setup kwargs, update plot_kwargs and text_kwargs.
        if plot_kwargs is None:
            plot_kwargs = {}
        if text_kwargs is None:
            text_kwargs = {}

        plot_kwargs = {'linewidth': linewidth, 'color': color,
                       **plot_kwargs, **kwargs}
        text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                       **text_kwargs, **kwargs}

        # Convert all units and types.
        location = np.asarray(location)  # For vector addition.
        length_metres = round(length) * metres_per_unit
        angle_rad = angle * np.pi / 180

        # End-point of bar.
        end = self._point_along_line(ax, location, length_metres,
                                     angle=angle_rad,
                                     tol=tol)

        # Coordinates are currently in axes coordinates, so use transAxes to
        # put into data coordinates. *zip(a, b) produces a list of x-coords,
        # then a list of y-coords.
        ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)

        # Push text away from bar in the perpendicular direction.
        midpoint = (location + end) / 2
        offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
        text_location = midpoint + offset

        # 'rotation' keyword argument is in text_kwargs.
        ax.text(*text_location, "{} {}".format(round(length), unit_name),
                rotation_mode='anchor', transform=ax.transAxes, **text_kwargs)
        return ax
######################################### Cartopy Scalebar class end ####################################
