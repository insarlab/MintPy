############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2018, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi, 2018              #
############################################################


import os, sys, glob
import time
from datetime import datetime as dt
import h5py
import numpy as np
from pysar.utils import readfile, datetime as ptime, utils as ut
from pysar.objects import ifgramDatasetNames, geometryDatasetNames, timeseriesDatasetNames

dataType = np.float32


########################################################################################
class timeseries:
    '''
    Time-series object for displacement of a set of SAR images from the same platform and track.
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.

    /                    Root level
    Attributes           Dictionary for metadata
    timeseries           group name
        /date            1D array of string  in size of (n,) in YYYYMMDD format
        /bperp           1D array of float32 in size of (n,) in meter.
        /timeseries      3D array of float32 in size of (n, l, w) in meter.
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'timeseries'

    def close(self):
        print('close timeseries: {}'.format(os.path.basename(self.file)))
        self.f.close()

    def open(self):
        print('openning timeseries file: {} ...'.format(self.file))
        self.f = h5py.File(self.file,'r')

        self.dateList = list(self.f['timeseries'].keys())
        dset = self.f['timeseries'].get(self.dateList[0])
        self.rows,self.cols = dset.shape
        self.numPixels = self.rows*self.cols
        self.numDates = len(self.dateList)



class ifgramStack:
    ''' Interferograms Stack object.
    /ifgramStack           Root level group name
        Attributes         Dictionary for metadata
        /date              2D array of string  in size of (m, 2   ) in YYYYMMDD format for master and slave date
        /bperp             1D array of float32 in size of (m,     ) in meter.
        /unwrapPhase       3D array of float32 in size of (m, l, w) in radian.
        /coherence         3D array of float32 in size of (m, l, w).
        /connectComponent  3D array of int16   in size of (m, l, w).           (optional)
        /wrapPhase         3D array of float32 in size of (m, l, w) in radian. (optional)
        /rangeOffset       3D array of float32 in size of (m, l, w).           (optional)
        /azimuthOffset     3D array of float32 in size of (m, l, w).           (optional)
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'ifgramStack'

    def open(self):
        print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f = h5py.File(self.file, 'r')
        dset = self.f[self.key].get(ifgramDatasetNames[0])
        self.numIfgram, self.length, self.width = dset.shape  #Size
        self.read_datetime()                                  #Time
        self.metadata = readfile.read_attribute(self.file)    #Metadata

    def read_datetime(self):
        self.f = h5py.File(self.file, 'r')
        d = self.f[self.key].get('date')
        self.masterDates = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in list(d[:,0])])
        self.slaveDates = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in list(d[:,1])])
        return self.masterDates, self.slaveDates

    def close(self):
        print('close {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f.close()

    def get_temporal_baseline(self):
        '''Return array of temporal baseline in days'''
        self.read_datetime()
        self.temporalBaseline = np.array([i.days for i in self.slaveDates - self.masterDates], dtype=np.int)
        return self.temporalBaseline

    def get_perp_baseline(self):
        '''Return array of spatial perpendicular baseline in meters'''
        self.perpBaseline = self.f[self.key].get('bperp')[:]
        return self.perpBaseline

    def read(self, datasetName=ifgramDatasetNames[0], box=None):
        dset = self.f[self.key].get(datasetName)
        if box:
            data = dset[:, box[1]:box[3], box[0]:box[2]]
        else:
            data = dset[:]
        return data















