# Major Package Imports
import os
import sys
import argparse
from datetime import datetime as dt

# Matplotlib Package Imports
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# PyQT Package Imports
from PyQt5 import QtWidgets as qt
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


class TimeSerieWidget(qt.QWidget):

    def __init__(self, parent=None, iargs=None):
        super(TimeSerieWidget, self).__init__(parent)

        # fig_v, ax_v, inps, ax_ts, fig_ts, second_plot_axis
        self.canvas = None
        self.toolbar = None
        self.fig_v = None
        self.ax_v = None
        self.inps = None
        self.ax_ts = None
        self.fig_ts = None
        self.second_plot_axis = None


        self.EXAMPLE = '''example:
                          tsview.py timeseries.h5 --ylim -10 10
                          tsview.py timeseries_demErr_plane.h5 -n 5 -m maskTempCoh.h5
                          tsview.py timeseries_demErr_plane.h5 --yx 300 400 --nodisplay --zero-first
                          tsview.py geo_timeseries_demErr_plane.h5 --lalo 33.250 131.665 --nodisplay
                       '''

        self.fig_v, self.canvas, self.toolbar = self.main(['/Users/joshua/Desktop/pysar/test_data/new_data/timeseries_ECMWF_demErr_plane.h5', '--figsize', '5', '3'])

        layout = qt.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)

        self.resize(1000, 600)

    def main(self, args=None):

        self.inps = self.cmd_line_parse(args)

        self.fig_v = Figure(self.inps.fig_size)
        self.canvas = FigureCanvas(self.fig_v)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.ax_v = self.fig_v.add_axes([0.035, 0.42, 0.5, 0.5])
        self.ax_ts = self.fig_v.add_axes([0.55, 0.62, 0.42, 0.3])
        self.second_plot_axis = self.fig_v.add_axes([0.55, 0.18, 0.42, 0.3])
        self.second_plot_axis.remove()

        if self.inps.disp_fig:
            self.canvas.draw()

        return self.fig_v, self.canvas, self.toolbar


    def create_parser(self):
        parser = argparse.ArgumentParser(description='Interactive Time-series Viewer', \
                                         formatter_class=argparse.RawTextHelpFormatter, \
                                         epilog=self.EXAMPLE)
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

    main = TimeSerieWidget()
    main.show()

    sys.exit(app.exec_())
