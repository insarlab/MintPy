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

    def __init__(self, parent=None):
        super(TimeSerieWidget, self).__init__(parent)

        # a figure instance to plot on
        self.figure = Figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = qt.QPushButton('Plot')
        #self.button.clicked.connect(self.plot)

        # set the layout
        layout = qt.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)



if __name__ == '__main__':
    app = qt.QApplication(sys.argv)

    main = TimeSerieWidget()
    main.show()

    sys.exit(app.exec_())
