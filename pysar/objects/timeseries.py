# Author: Zhang Yunjun, 2018 Mar.

import os, sys
import h5py
import numpy as np


chunk_shape = (128,128)
dataType = np.float32


class timeseries:
    '''
    Time-series object for displacement of a set of SAR images from the same platform and track.
    Attributes are saved in the root level.
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.
    
    /                    Root level
    Attributes           Dictionary for metadata
    timeseries           group name
        /date            1D array of string in size of (n,) in YYYYMMDD format
        /bperp           1D array of float32 in size of (n,) in meter.
        /timeseries      3D array of float32 in size of (n, l, w) in meter.
    '''
    def __init__(self, file=None):
        self.file = file
    


    def close(self):
        print('close timeseries: {}'.format(os.path.basename(self.file)))
        self.h5.close()