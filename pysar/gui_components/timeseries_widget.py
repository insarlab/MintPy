# Major Package Imports
import os
import sys
import argparse
from datetime import datetime as dt

# PyQT Package Imports
from PyQt5 import QtWidgets as qt
from PyQt5 import Qt as Qt

# Matplotlib Package Imports
from matplotlib.widgets import Slider
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# Scientific Computing Package Imports
import h5py
import numpy as np
import scipy.stats as stats

# Pysar Package Imports
from pysar.objects import timeseries
from pysar.utils import readfile, ptime, utils as ut, plot as pp
from pysar.mask import mask_matrix


class TimeSeriesWidget(qt.QWidget):

    def __init__(self, parent=None, iargs=None):

        super(TimeSeriesWidget, self).__init__(parent)

        # Declaration of All Instance Variables
        self.canvas, self.toolbar                               =       None, None
        self.fig_v, self.fig_ts, self.ax_v, self.ax_ts          =       None, None, None, None
        self.second_plot_axis, self.second_plot_axis_visible    =       None, False

        self.atr, self.k, self.h5                               =       None, None, None
        self.dateList, self.date_num, self.tims                 =       None, None, None

        self.width, self.length                                 =       None, None
        self.ullat, self.lat_step, self.lat                     =       None, None, None
        self.ullon, self.lon_step, self.lon                     =       None, None, None

        self.d_v, self.d_ts, self.data_lim                      =       None, None, None

        self.p1_scatter_point, self.p1_x, self.p1_y             =       None, None, None
        self.p2_scatter_point, self.p2_x, self.p2_y             =       None, None, None

        self.inps, self.mask, self.img, self.tslider            =       None, None, None, None
        self.first_data_point = None


        # Setup Matplotlib Widgets
        self.fig_v, self.canvas, self.toolbar = self.main(iargs)

        # Set layout parameters
        layout = qt.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)

        self.resize(1000, 600)

    def main(self, args=None):

        # Parse command line arguments into inps
        self.inps = self.cmd_line_parse(args)

        # Do inital setup of plot axes and elements
        self.setup_plot()

        # Create new figure object
        self.fig_v = Figure(self.inps.fig_size)
        # Create new figure canvas object
        self.canvas = FigureCanvas(self.fig_v)
        # Create new MPL toolbar object
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Set click focus to be on the figure canvas
        self.canvas.setFocusPolicy(Qt.Qt.ClickFocus)
        self.canvas.setFocus()

        # Set image plotting axis
        self.ax_v = self.fig_v.add_axes([0.035, 0.42, 0.5, 0.5])

        # Do initial plot configuration
        self.configure_plot()

        # Set timeseries data axis
        self.ax_ts = self.fig_v.add_axes([0.55, 0.62, 0.42, 0.3])

        # Set second timeseries data axis
        self.second_plot_axis = self.fig_v.add_axes([0.55, 0.18, 0.42, 0.3])
        # Remove second timeseries data axis from screen
        self.second_plot_axis.remove()

        # Read in error file and compute error list
        self.read_error_list()

        # Plot data from initial point on timeseries axis
        self.plot_data_from_inital_point()

        # Save output to file if necesarry
        self.save_output()

        # Setup button press event to a new plot_timeseries_event
        self.first_data_point = self.fig_v.canvas.mpl_connect('button_press_event', self.plot_timeseries_event)

        # Refresh the canvas to display changes
        if self.inps.disp_fig:
            self.canvas.draw()

        # Disconnect the button press listener
        #self.fig_v.canvas.mpl_disconnect(self.first_data_point)

        return self.fig_v, self.canvas, self.toolbar


    ###########################################   Setting Up Plot     ###########################################

    '''
        Sets up plotting interface
    '''
    def setup_plot(self):

        # Time Series Info
        self.read_timeseries_info()

        # Read exclude dates
        self.exclude_dates()

        # Zero displacement for 1st acquisition
        self.set_zero_displacement()

        # File Size
        self.compute_file_size()

        # Latitude Longitude Parameters
        self.compute_lat_lon_params()

        # Initial Pixel Coordinates
        self.set_inital_pixel_coords()

        # Display Unit
        self.set_unit_fraction()

        # Flip up-down / left-right
        self.flip_map()

        # Mask file
        self.set_mask()

        # Initial Map
        self.set_initial_map()

    '''
        Reads basic information about timeseries file being viewed
    '''
    def read_timeseries_info(self):

        # Read file attributes from timeseries file
        self.atr = readfile.read_attribute(self.inps.timeseries_file)

        # Determine file type of input file
        self.k = self.atr['FILE_TYPE']
        print(('input file is ' + self.k + ': ' + self.inps.timeseries_file))

        # Raise error if file is not a timeseries file
        if not self.k in ['timeseries', 'GIANT_TS']:
            raise ValueError('Only timeseries file is supported!')

        # Read data out of h5 input file
        self.h5 = h5py.File(self.inps.timeseries_file, 'r')

        # Sets datelist for each valid file type
        if self.k in ['GIANT_TS']:
            self.dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in self.h5['dates'][:].tolist()]
        else:
            self.dateList = timeseries(self.inps.timeseries_file).get_date_list()

        # Determine number of dates in datelist
        self.date_num = len(self.dateList)

        # Formats dates into useful objects
        self.inps.dates, self.tims = ptime.date_list2vector(self.dateList)

    '''
        Sets dates to be excluded from timeseries data points
    '''
    def exclude_dates(self):

        # Determines if any dates are to be excluded
        if self.inps.ex_date_list:

            # Gets list of dates to be excluded
            input_ex_date = list(self.inps.ex_date_list)
            self.inps.ex_date_list = []

            if input_ex_date:

                # For each date to be excluded, format the date properly and add it
                # to the excluded date list
                for ex_date in input_ex_date:

                    if os.path.isfile(ex_date):
                        ex_date = ptime.read_date_list(ex_date)
                    else:
                        ex_date = [ptime.yyyymmdd(ex_date)]

                    self.inps.ex_date_list += list(set(ex_date) - set(self.inps.ex_date_list))

                # Remove dates that don't exist in the input file
                self.inps.ex_date_list = sorted(list(set(self.inps.ex_date_list).intersection(self.dateList)))
                self.inps.ex_dates = ptime.date_list2vector(self.inps.ex_date_list)[0]
                self.inps.ex_idx_list = sorted([self.dateList.index(i) for i in self.inps.ex_date_list])
                print(('exclude date:' + str(self.inps.ex_date_list)))

    '''
        Sets the 'zero' value for the plot
    '''
    def set_zero_displacement(self):

        if self.inps.zero_first:

            # Set 'zero-value' for the plot to be true '0' or smallest value in excluded date list
            if self.inps.ex_date_list:
                self.inps.zero_idx = min(list(set(range(self.date_num)) - set(self.inps.ex_idx_list)))
            else:
                self.inps.zero_idx = 0

    '''
        Computed dimensions of file (width, length)
    '''
    def compute_file_size(self):

        # Computes and assigns the length and width of the file
        self.length = int(self.atr['LENGTH'])
        self.width = int(self.atr['WIDTH'])
        print(('data size in [y0,y1,x0,x1]: [%d, %d, %d, %d]' % (0, self.length, 0, self.width)))

    '''
        Computes parameters needed for setting lat/lon values
    '''
    def compute_lat_lon_params(self):

        try:

            # Sets series of parameters needed for computing visual length of a
            # latitude/longitude parameter
            self.ullon = float(self.atr['X_FIRST'])
            self.ullat = float(self.atr['Y_FIRST'])
            self.lon_step = float(self.atr['X_STEP'])
            self.lat_step = float(self.atr['Y_STEP'])

            # Computes visual length of latitude/longitude on screen
            lrlon = self.ullon + self.width * self.lon_step
            lrlat = self.ullat + self.length * self.lat_step
            
            print(('data size in [lat0,lat1,lon0,lon1]: [%.4f, %.4f, %.4f, %.4f]' % (self.lrlat, 
                                                                                     self.ullat,
                                                                                     self.ullon,
                                                                                     self.lrlon)
                   ))
            return self.ullon, self.ullat, self.lon_step, self.lat_step, lrlon, lrlat
        except:
            pass

    '''
        Sets coordinates of initial pixel
    '''
    def set_inital_pixel_coords(self):

        # Sets (x, y) coordinates for provided (lat, lon) coordinates of initial pixel
        if self.inps.lalo and 'Y_FIRST' in list(self.atr.keys()):
            y, x = self.set_yx_coords(self.inps.lalo[0], self.inps.lalo[1])
            self.inps.yx = [y, x]

        # Sets (x, y) coordinates for provided (lat, lon) coordinates of reference pixel
        if self.inps.ref_lalo and 'Y_FIRST' in list(self.atr.keys()):
            y, x = self.set_yx_coords(self.inps.ref_lalo[0], self.inps.ref_lalo[1])
            self.inps.ref_yx = [y, x]

        # Set unit scaling fraction based on display unit
        self.set_unit_fraction()

    '''
        Sets x/y coordinates from lalo values
    '''
    def set_yx_coords(self, y_input, x_input):

        # Computes (x, y) coords from provided (lat, lon) values
        y = int((y_input - self.ullat) / self.lat_step + 0.5)
        x = int((x_input - self.ullon) / self.lon_step + 0.5)

        return y, x

    '''
        Sets multiplier for data based on display unit
    '''
    def set_unit_fraction(self):

        # Set unit scaling fraction based on display unit
        if self.inps.disp_unit == 'cm':
            self.inps.unit_fac = 100.0
        elif self.inps.disp_unit == 'm':
            self.inps.unit_fac = 1.0
        elif self.inps.disp_unit == 'dm':
            self.inps.unit_fac = 10.0
        elif self.inps.disp_unit == 'mm':
            self.inps.unit_fac = 1000.0
        elif self.inps.disp_unit == 'km':
            self.inps.unit_fac = 0.001
        else:
            raise ValueError('Un-recognized unit: ' + self.inps.disp_unit)

        # Defines default unit for each file type
        if self.k in ['GIANT_TS']:
            print('data    unit: mm')
            self.inps.unit_fac *= 0.001
        else:
            print('data    unit: m')
        print(('display unit: ' + self.inps.disp_unit))

    '''
        Flips display map ud/lr
    '''
    def flip_map(self):

        # Flips map left/right, or up/down, as specified
        if self.inps.auto_flip:
            self.inps.flip_lr, self.inps.flip_ud = pp.auto_flip_direction(self.atr)
        else:
            self.inps.flip_ud = False
            self.inps.flip_lr = False

    '''
        Sets mask for map
    '''
    def set_mask(self):

        # If mask has not been explcitly set
        if not self.inps.mask_file:

            # Create file list of potential mask files
            if os.path.basename(self.inps.timeseries_file).startswith('geo_'):
                file_list = ['geo_maskTempCoh.h5']
            else:
                file_list = ['maskTempCoh.h5', 'mask.h5']

            # Sets mask file from potential list if exists
            try:
                self.inps.mask_file = ut.get_file_list(file_list)[0]
            except:
                self.inps.mask_file = None

        # Attempt to read in mask file and create masking matrix
        try:
            self.mask = readfile.read(self.inps.mask_file, datasetName='mask')[0]
            self.mask[self.mask != 0] = 1
            print(('load mask from file: ' + self.inps.mask_file))
        except:
            self.mask = None
            print('No mask used.')

    '''
        Sets the initial map plot
    '''
    def set_initial_map(self):

        # Set inital data values for map
        self.d_v = self.h5['timeseries'][self.inps.epoch_num][:] * self.inps.unit_fac

        # Read in initial data values from file, provided the starting epoch
        print(str(self.dateList))
        self.d_v = readfile.read(self.inps.timeseries_file, datasetName=self.dateList[self.inps.epoch_num])[0] * self.inps.unit_fac

        # If reference data is specified, read in data values for that date, and subtract
        # values from inital values to create data value relative to the provided date
        if self.inps.ref_date:
            self.inps.ref_d_v = readfile.read(self.inps.timeseries_file, datasetName=self.inps.ref_date)[0] * self.inps.unit_fac
            self.d_v -= self.inps.ref_d_v

        # If mask if set, mask the data values using the computed masking file and matrix
        if self.mask is not None:
            self.d_v = mask_matrix(self.d_v, self.mask)

        # If a reference pixel was set, subtract data values of that (x, y) point from the
        # initial data values matrix
        if self.inps.ref_yx:
            self.d_v -= self.d_v[self.inps.ref_yx[0], self.inps.ref_yx[1]]

        # Determine minimum and maximum data values in matrix
        self.data_lim = [np.nanmin(self.d_v), np.nanmax(self.d_v)]

        # Assign y-limits to min and max values, if not set explicitly
        if not self.inps.ylim_mat:
            self.inps.ylim_mat = self.data_lim

        print(('Initial data range: ' + str(self.data_lim)))
        print(('Display data range: ' + str(self.inps.ylim_mat)))

        print(('Initial data range: ' + str(self.data_lim)))
        print(('Display data range: ' + str(self.inps.ylim)))

    ###########################################   Configuring Plot  ###########################################

    def configure_plot(self):

        # DEM File
        self.set_dem_file()

        # Reference Pixel
        self.set_map_reference_pixel()

        # Initial Pixel
        self.set_plot_axis_params()

        # Flip Axis
        self.flip_axis()

        # Construct Color Bar
        self.make_color_bar()

        # Construct Time Slider
        self.make_time_slider()

    '''
        Sets DEM topography file for the map
    '''
    def set_dem_file(self):

        if self.inps.dem_file:
            dem = readfile.read(self.inps.dem_file, datasetName='height')[0]
            self.ax_v = pp.plot_dem_yx(self.ax_v, dem)

        self.img = self.ax_v.imshow(self.d_v, cmap=self.inps.colormap, clim=self.inps.ylim_mat, interpolation='nearest')

    '''
        Sets reference pixel on map
    '''
    def set_map_reference_pixel(self):

        if self.inps.ref_yx:
            self.d_v -= self.d_v[self.inps.ref_yx[0], self.inps.ref_yx[1]]
            self.ax_v.plot(self.inps.ref_yx[1], self.inps.ref_yx[0], 'ks', ms=6)
        else:
            try:
                self.ax_v.plot(int(self.atr['ref_x']), int(self.atr['ref_y']), 'ks', ms=6)
            except:
                pass

    '''
        Sets axis limits, labels, and titles
    '''
    def set_plot_axis_params(self):

        if self.inps.yx:
            self.ax_v.plot(self.inps.yx[1], self.inps.yx[0], 'ro', markeredgecolor='black')

        self.ax_v.set_xlim(0, np.shape(self.d_v)[1])
        self.ax_v.set_ylim(np.shape(self.d_v)[0], 0)
        self.ax_v.format_coord = self.format_coord

        # Title and Axis Label
        self.ax_v.set_title('N = %d, Time = %s' % (self.inps.epoch_num, self.inps.dates[self.inps.epoch_num].strftime('%Y-%m-%d')))

        if not 'Y_FIRST' in list(self.atr.keys()):
            self.ax_v.set_xlabel('Range')
            self.ax_v.set_ylabel('Azimuth')

    '''
        Flips axis lr/ud
    '''
    def flip_axis(self):

        if self.inps.flip_lr:
            self.ax_v.invert_xaxis()
            print('flip map left and right')
        if self.inps.flip_ud:
            self.ax_v.invert_yaxis()
            print('flip map up and down')

    '''
        Creates colorbar for figure
    '''
    def make_color_bar(self):

        # Colorbar
        cbar_axes = self.fig_v.add_axes([0.075, 0.28, 0.40, 0.03])
        cbar = self.fig_v.colorbar(self.img, cax=cbar_axes, orientation='horizontal')
        cbar.set_label('Displacement [%s]' % self.inps.disp_unit)

    '''
        Creates timeseries slider for figure
    '''
    def make_time_slider(self):

        ax_time = self.fig_v.add_axes([0.07, 0.10, 0.37, 0.07], facecolor='lightgoldenrodyellow', yticks=[])
        self.tslider = Slider(ax_time, '', self.tims[0], self.tims[-1], valinit=self.tims[self.inps.epoch_num])
        self.tslider.ax.bar(self.tims, np.ones(len(self.tims)), facecolor='black', width=0.01, ecolor=None)
        self.tslider.ax.set_xticks(np.round(np.linspace(self.tims[0], self.tims[-1], num=5) * 100) / 100)
        self.tslider.on_changed(self.time_slider_update)

    def time_slider_update(self, val):
        '''Update Displacement Map using Slider'''

        timein = self.tslider.val
        idx_nearest = np.argmin(np.abs(np.array(self.tims) - timein))

        self.ax_v.set_title('N = %d, Time = %s' % (idx_nearest, self.inps.dates[idx_nearest].strftime('%Y-%m-%d')))
        self.d_v = self.h5[self.k][idx_nearest][:] * self.inps.unit_fac

        if self.inps.ref_date:
            self.d_v -= self.inps.ref_d_v

        if self.mask is not None:
            self.d_v = mask_matrix(self.d_v, self.mask)

        if self.inps.ref_yx:
            self.d_v -= self.d_v[self.inps.ref_yx[0], self.inps.ref_yx[1]]

        self.img.set_data(self.d_v)
        self.fig_v.canvas.draw()

    def format_coord(self, x, y):
        '''Formats x, y coordinates into useful output string (used for creating plot titles)'''

        col = int(x + 0.5)
        row = int(y + 0.5)
        if 0 <= col < self.width and 0 <= row < self.length:
            z = self.d_v[row, col]
            try:
                self.lon = self.ullon + x * self.lon_step
                self.lat = self.ullat + y * self.lat_step
                return 'x=%.0f, y=%.0f, value=%.4f, lon=%.4f, lat=%.4f' % (x, y, z, self.lon, self.lat)
            except:
                return 'x=%.0f, y=%.0f, value=%.4f' % (x, y, z)

    ###########################################   Miscellaneous Functions  ###########################################

    def read_error_list(self):

        self.inps.error_ts = None
        if self.inps.error_file:
            error_file_content = np.loadtxt(self.inps.error_file, dtype=str)
            self.inps.error_ts = error_file_content[:, 1].astype(np.float) * self.inps.unit_fac
            if self.inps.ex_date_list:
                e_ts = self.inps.error_ts[:]
                self.inps.ex_error_ts = np.array([e_ts[i] for i in self.inps.ex_idx_list])
                self.inps.error_ts = np.array([e_ts[i] for i in range(self.date_num) if i not in self.inps.ex_idx_list])

    '''
        Saves figure and data to output file
    '''
    def save_output(self):

        if self.inps.save_fig and self.inps.yx:
            print(('save info for pixel ' + str(self.inps.yx)))
            if not self.inps.fig_base:
                self.inps.fig_base = 'y%d_x%d' % (self.inps.yx[0], self.inps.yx[1])

            # TXT - point time series
            outName = self.inps.fig_base + '_ts.txt'
            header_info = 'timeseries_file=' + self.inps.timeseries_file
            header_info += '\ny=%d, x=%d' % (self.inps.yx[0], self.inps.yx[1])

            try:
                lat = self.ullat + self.inps.yx[0] * self.lat_step
                lon = self.ullon + self.inps.yx[1] * self.lon_step
                header_info += '\nlat=%.6f, lon=%.6f' % (lat, lon)
            except:
                pass

            if self.inps.ref_yx:
                header_info += '\nreference pixel: y=%d, x=%d' % (self.inps.ref_yx[0], self.inps.ref_yx[1])
            else:
                header_info += '\nreference pixel: y=%s, x=%s' % (self.atr['ref_y'], self.atr['ref_x'])

            header_info += '\nunit=m/yr'
            np.savetxt(outName, list(zip(np.array(self.dateList), np.array(self.d_ts) / self.inps.unit_fac)), fmt='%s', \
                       delimiter='    ', header=header_info)
            print(('save time series displacement in meter to ' + outName))

            # Figure - point time series
            outName = self.inps.fig_base + '_ts.pdf'
            self.fig_ts.savefig(outName, bbox_inches='tight', transparent=True, dpi=self.inps.fig_dpi)
            print(('save time series plot to ' + outName))

            # Figure - map
            outName = self.inps.fig_base + '_' + self.dateList[self.inps.epoch_num] + '.png'
            self.fig_v.savefig(outName, bbox_inches='tight', transparent=True, dpi=self.inps.fig_dpi)
            print(('save map plot to ' + outName))

    '''
        Plots initial point on map and sets timeseries data points on scatter plot
        to appropriate values
    '''
    def plot_data_from_inital_point(self):

        if self.inps.yx:
            self.d_ts = self.update_timeseries(self.inps.yx[0], self.inps.yx[1], 1)
        else:
            self.d_ts = np.zeros(len(self.tims))
            self.ax_ts, scatter = self.plot_timeseries_scatter(self.ax_ts, self.d_ts, self.inps)

    ###########################################   Timeseries Data Functions #########################################

    def plot_timeseries_errorbar(self, ax, dis_ts, inps):
        '''Plots errorbars for timeseries data'''

        dates = list(self.inps.dates)
        self.d_ts = dis_ts[:]
        if self.inps.ex_date_list:
            # Update displacement time-series
            dates = sorted(list(set(self.inps.dates) - set(self.inps.ex_dates)))
            ex_d_ts = np.array([dis_ts[i] for i in self.inps.ex_idx_list])
            self.d_ts = np.array([dis_ts[i] for i in range(self.date_num) if i not in inps.ex_idx_list])
            # Plot excluded dates
            (_, caps, _) = ax.errorbar(self.inps.ex_dates, 
                                       ex_d_ts, 
                                       yerr=self.nps.ex_error_ts, 
                                       fmt='-o', color='gray', 
                                       ms=self.inps.marker_size, 
                                       lw=0, alpha=1, mfc='gray',
                                       elinewidth=self.inps.edge_width, 
                                       ecolor='black', 
                                       capsize=self.inps.marker_size * 0.5)
            
            for cap in caps:  cap.set_markeredgewidth(self.inps.edge_width)
        
        # Plot kept dates
        (_, caps, _) = ax.errorbar(dates, 
                                   self.d_ts, 
                                   yerr=self.inps.error_ts, 
                                   fmt='-o',
                                   ms=self.inps.marker_size, 
                                   lw=0, alpha=1,
                                   elinewidth=self.inps.edge_width, 
                                   ecolor='black', 
                                   capsize=self.inps.marker_size * 0.5)
        
        for cap in caps:  cap.set_markeredgewidth(self.inps.edge_width)
        return ax

    def plot_timeseries_scatter(self, ax, dis_ts, inps, plot_num=1):
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

        dates = list(self.inps.dates)

        print(dis_ts)

        self.d_ts = dis_ts[:]
        if self.inps.ex_date_list:
            # Update displacement time-series
            dates = sorted(list(set(self.inps.dates) - set(self.inps.ex_dates)))
            ex_d_ts = np.array([dis_ts[i] for i in self.inps.ex_idx_list])
            self.d_ts = np.array([dis_ts[i] for i in range(self.date_num) if i not in self.inps.ex_idx_list])
            # Plot excluded dates
            ax.scatter(self.inps.ex_dates, ex_d_ts, s=self.inps.marker_size ** 2, color='gray')  # color='crimson'

        # Plot kept dates
        color = 'blue'
        if plot_num == 2:
            color = 'crimson'
        print(('Color is ' + color))

        print(dates)
        print(self.d_ts)

        scatter = ax.scatter(dates, self.d_ts, s=self.inps.marker_size ** 2, label='1', color=color)

        return ax, scatter

    def update_timeseries(self, y, x, plot_number, data_only=False):
        '''Plot point time series displacement at pixel [y, x]
            Inputs:
                y           : int, y coordinate to update
                x           : int, x coordinate to update
                plot_number : int, plot number (1/2) to update
                data_only   : bool, compute and return data only, or set remainder of plot variables
            Outputs:
                d_ts        : [float], timeseries data at x, y point
        '''

        self.set_scatter_coords(plot_number, x, y)

        if plot_number == 1:
            axis = self.ax_ts
        #else:
            #axis = second_plot_axis

            self.d_ts = []
        for i, date in enumerate(self.dateList):
            d = self.h5['timeseries'][i][y, x]
            if self.inps.ref_yx:
                d -= self.h5['timeseries'][i][self.inps.ref_yx[0], self.inps.ref_yx[1]]

            self.d_ts.append(d * self.inps.unit_fac)

        if self.inps.zero_first:
            self.d_ts -= self.d_ts[self.inps.zero_idx]

        # Returns computed data without setting any plot or figure parameters
        if data_only:
            return self.d_ts

        axis.cla()
        if self.inps.error_file:
            axis = self.plot_timeseries_errorbar(self.ax_ts, self.d_ts, self.inps)
        else:
            axis, scatter = self.plot_timeseries_scatter(axis, self.d_ts, self.inps, plot_number)
            scatter.set_label('2')

        axis.set_ylim(self.inps.ylim_mat[0] * 2, self.inps.ylim_mat[1] * 2)
        for tick in axis.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.inps.font_size)

        # Title
        title_ts = self.set_axis_title(x, y)
        if self.inps.disp_title:
            axis.set_title(title_ts)

        axis = pp.auto_adjust_xaxis_date(axis, self.tims, fontSize=self.inps.font_size)[0]
        axis.set_xlabel('Time', fontsize=self.inps.font_size)
        axis.set_ylabel('Displacement [%s]' % self.inps.disp_unit, fontsize=self.inps.font_size)

        self.fig_v.canvas.draw()

        # Print to terminal
        print('\n---------------------------------------')
        print(title_ts)
        print(self.d_ts)

        # Slope estimation
        self.estimate_slope()

        return self.d_ts

    def set_axis_title(self, x, y):
        '''Sets title of axis for a given X, Y Point
            Inputs:
                x   : int, x coordinate
                y   : int, y coordinate
            Outputs:
                title_ts    : string, computed axis title
        '''

        if x is None:
            title_ts = 'No Point Selected'
        else:

            title_ts = 'Y = %d, X = %d' % (y, x)
            try:
                lat, lon = self.xy_to_lat_lon(x, y)
                title_ts += ', lat = %.4f, lon = %.4f' % (lat, lon)
            except:
                pass

        return title_ts

    def xy_to_lat_lon(self, x, y):
        '''Converst x,y coordinated to lat/lon coordinates
            Inputs:
                x   : int, x coordinate
                y`  : int, y coordinate
            Outputs:
                latitude    : double, computed latitude coordinate
                longitude   : double, computed longitude coordinate
        '''

        latitude = self.ullat + y * self.lat_step
        longitude = self.ullon + x * self.lon_step

        return latitude, longitude

    def estimate_slope(self):
        '''Estimates slope of timeseries scatter data'''

        if self.inps.ex_date_list:
            tims_kept = [self.tims[i] for i in range(self.date_num) if i not in self.inps.ex_idx_list]
            d_ts_kept = [self.d_ts[i] for i in range(self.date_num) if i not in self.inps.ex_idx_list]
            d_slope = stats.linregress(np.array(tims_kept), np.array(d_ts_kept))
        else:
            d_slope = stats.linregress(np.array(self.tims), np.array(self.d_ts))

        print(('linear velocity: %.2f +/- %.2f [%s/yr]' % (d_slope[0], d_slope[4], self.inps.disp_unit)))

    def set_scatter_coords(self, plot_number, x, y):
        '''Sets the coordinates or the starting scatter point
            Inputs:
                plot_number     : int, the plot number (1 or 2)
                x               : int, x coordinate
                y               : int, y coordinate
        '''

        if plot_number == 1:  # Set scatter point 1 coordinates
            self.p1_x, self.p1_y = x, y
        else:  # Set scatter point 2 coordinates
            self.p2_x, self.p2_y = x, y

    def plot_timeseries_event(self, event):
        '''Event function to get y/x from button press'''

        if event.inaxes != self.ax_v:
            return

        ii = int(event.ydata + 0.5)
        jj = int(event.xdata + 0.5)

        if event.button == 1:  # Compute and update plot 1 data on left mouse-click

            if self.p1_scatter_point is not None:
                self.p1_scatter_point.remove()  # remove previous scatter point

            self.p1_scatter_point = self.ax_v.scatter(event.xdata, event.ydata, s=50, c='red',
                                            marker='o')  # place new sactter point

            self.d_ts = self.update_timeseries(ii, jj, 1)  # update timeseries scatter plot for plot 1

        elif event.button == 3 and self.second_plot_axis_visible:  # Compute and update plot 2 on right mouse-click

            if self.p2_scatter_point is not None:
                self.p2_scatter_point.remove()  # remove previous scatter point

            self.p2_scatter_point = self.ax_v.scatter(event.xdata, event.ydata, s=50, c='blue',
                                            marker='o')  # place new scatter point

            self.d_ts = self.update_timeseries(ii, jj, 2)  # update timeseries scatter plot for plot 2

    ###########################################   Parser Information  ###########################################

    def create_parser(self):

        EXAMPLE = '''example:
                      tsview.py timeseries.h5 --ylim -10 10
                      tsview.py timeseries_demErr_plane.h5 -n 5 -m maskTempCoh.h5
                      tsview.py timeseries_demErr_plane.h5 --yx 300 400 --nodisplay --zero-first
                      tsview.py geo_timeseries_demErr_plane.h5 --lalo 33.250 131.665 --nodisplay
                   '''

        parser = argparse.ArgumentParser(description='Interactive Time-series Viewer', \
                                         formatter_class=argparse.RawTextHelpFormatter, \
                                         epilog=EXAMPLE)
        parser.add_argument('timeseries_file', help='time series file to display')
        parser.add_argument('-n', dest='epoch_num', metavar='NUM', type=int, default='-2', \
                            help='Epoch/slice number to display, default: the 2nd last.')
        parser.add_argument('-m', '--mask', dest='mask_file', \
                            help='mask to use. Default: geo_maskTempCoh.h5 for geocoded file and maskTempCoh.h5 for radar file')
        parser.add_argument('--error', dest='error_file', help='txt file with error for each date.')
        parser.add_argument('--dem', dest='dem_file', help='DEM file for background shaed relief')

        pixel = parser.add_argument_group('Pixel Input')
        pixel.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2, \
                           help='initial pixel to plot in Y/X coord')
        pixel.add_argument('--lalo', type=float, metavar=('LAT', 'LON'), nargs=2, \
                           help='initial pixel to plot in lat/lon coord')
        pixel.add_argument('--ref-yx', dest='ref_yx', type=int, metavar=('Y', 'X'), nargs=2, \
                           help='change reference pixel to input location')
        pixel.add_argument('--ref-lalo', dest='ref_lalo', type=float, metavar=('LAT', 'LON'), nargs=2, \
                           help='change reference pixel to input location')

        output = parser.add_argument_group('Output Setting')
        output.add_argument('-o', '--output', dest='fig_base', help='Figure base name for output files')
        output.add_argument('--save', action='store_true', dest='save_fig', \
                            help='save data and plot to files')
        output.add_argument('--nodisplay', action='store_false', dest='disp_fig', \
                            help='save data and plot to files and do not display figures\n')
        output.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=150, \
                            help='DPI - dot per inch - for display/write')

        disp = parser.add_argument_group('Display Setting')
        disp.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                          default=[10.0, 5.0], \
                          help='Figure size in inches - width and length. Default: 10.0 5.0\n' + \
                               'i.e. 3.5 2 for ppt; ')
        disp.add_argument('--ylim', dest='ylim', nargs=2, metavar=('YMIN', 'YMAX'), type=float, \
                          help='Y limits for point plotting.')
        disp.add_argument('--ylim-mat', dest='ylim_mat', nargs=2, metavar=('YMIN', 'YMAX'), type=float, \
                          help='Display limits for matrix plotting.')
        disp.add_argument('--ref-date', dest='ref_date', help='Change reference date for display')
        disp.add_argument('--exclude', '--ex', dest='ex_date_list', nargs='*', help='Exclude date shown as gray.')
        disp.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true', \
                          help='Set displacement at first acquisition to zero.')
        disp.add_argument('-u', dest='disp_unit', metavar='UNIT', default='cm', \
                          help='unit for display. Default: cm')
        disp.add_argument('-c', '--colormap', dest='colormap', default='jet', \
                          help='colormap used for display, i.e. jet, RdBu, hsv, jet_r etc.\n'
                               'Support colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html')
        disp.add_argument('-s', '--fontsize', dest='font_size', type=int, default=10, help='Font size for display')
        disp.add_argument('--notitle', dest='disp_title', action='store_false', help='Do not display title in TS plot.')
        disp.add_argument('--no-flip', dest='auto_flip', action='store_false', \
                          help='Turn off auto flip based on orbit direction.\n' + \
                               'Default: flip left-right for descending data in radar coord\n' + \
                               '         flip up-down    for ascending  data in radar coord\n' + \
                               '         no flip for data in geo coord')
        disp.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=12.0, \
                          help='Point marker size. Default: 12.0')
        # disp.add_argument('--mc','--markercolor', dest='marker_color', default='crimson',\
        #                  help='Point marker color. Default: crimson')
        disp.add_argument('--ew', '--edgewidth', dest='edge_width', type=float, default=1.0, \
                          help='Edge width. Default: 1.0')
        return parser

    def cmd_line_parse(self, iargs=None):
        parser = self.create_parser()
        inps = parser.parse_args(args=iargs)

        if (not inps.disp_fig or inps.fig_base) and not inps.save_fig:
            inps.save_fig = True
        if inps.ylim:
            inps.ylim = sorted(inps.ylim)
        return inps


if __name__ == '__main__':
    app = qt.QApplication(sys.argv)

    main = TimeSeriesWidget(iargs=sys.argv[1:])
    main.show()

    sys.exit(app.exec_())
