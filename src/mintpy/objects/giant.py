"""Classes for HDF5/GIAnT file operations"""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#   from mintpy.objects import giantTimeseries, giantIfgramStack

import os
from datetime import datetime as dt

import h5py
import numpy as np

GIANT_DSET_NAMES = [
    'recons',        #Reconstructed filtered time-series in mm
    'rawts',         #Raw time-series in mm
    'ifgcnt',        #Number of interferograms used for every pixel.
    'figram',        #Deramped + atmosphere corrected interferograms. in mm
    'igram',         #Unwrapped IFGs read straight from files in mm
    'cmask',         #Common mask for pixels
    'igram_aps',     #Atmosphere corrected interferogram stack in mm
    'sar_aps',       #Atmospheric phase screen for each of the SAR scenes in mm
]


########################################################################################
FILE_STRUCTURE_GIANT_TIMESERIES = """
/                Root level
/cmask           2D array of float32 in size of (   l, w).
/dates           1D array of float32 in size of (n,     ) in ordinal date format
/recons          3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - filtered
/rawts           3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - un-filtered
/ifgcnt          2D array of int32   in size of (   l, w), number of interferograms used for every pixel
"""

class giantTimeseries:
    """
    Time-series object for displacement of a set of SAR images from the same platform and track.
    """

    def __init__(self, file=None):
        self.file = file
        self.name = 'giantTimeseries'

    def open(self, print_msg=True):
        if print_msg:
            print(f'open {self.name} file: {os.path.basename(self.file)}')
        self.get_size()
        self.get_metadata()
        self.numPixel = self.length * self.width

        # Time Info
        self.times = np.array([dt.strptime(i, "%Y%m%d") for i in self.dateList])
        self.tbase = np.array([i.days for i in self.times - self.times[self.refIndex]], dtype=np.float32)
        self.yearList = [i.year + (i.timetuple().tm_yday-1)/365.25 for i in self.times]  #e.g. 2014.95

        # Dataset Info
        with h5py.File(self.file, 'r') as f:
            # get existed datasetNames in the order of GIANT_DSET_NAMES
            dsNames = [i for i in f.keys()
                       if (isinstance(f[i], h5py.Dataset)
                           and f[i].shape[-2:] == (self.length, self.width))]
            self.datasetNames = [i for i in GIANT_DSET_NAMES if i in dsNames]
            self.datasetNames += [i for i in dsNames if i not in GIANT_DSET_NAMES]

            self.sliceList = []
            for dsName in self.datasetNames:
                ds = f[dsName]
                if len(ds.shape) == 3:
                    self.sliceList += [f'{dsName}-{i}' for i in self.dateList]
                elif len(ds.shape) == 2:
                    self.sliceList.append(dsName)
                else:
                    raise ValueError(('un-recognized dataset dimension for {}:'
                                      ' {}').format(dsName, ds.shape))

    def get_size(self, dsName='recons'):
        with h5py.File(self.file, 'r') as f:
            self.numDate, self.length, self.width = f[dsName].shape
        return self.numDate, self.length, self.width

    def get_date_list(self):
        with h5py.File(self.file, 'r') as f:
            self.dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d')
                             for i in f['dates'][:].tolist()]
        return self.dateList

    def get_metadata(self):
        # read existing metadata
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value

        # size
        self.get_size()
        self.metadata['LENGTH'] = str(self.length)
        self.metadata['WIDTH'] = str(self.width)

        # ref_date/index
        dateList = self.get_date_list()
        if 'REF_DATE' not in self.metadata.keys():
            self.metadata['REF_DATE'] = dateList[0]
        self.refIndex = dateList.index(self.metadata['REF_DATE'])
        return self.metadata




########################################################################################
FILE_STRUCTURE_GIANT_IFGRAMSTACK = """
/                Root level
/cmask           2D array of float32 in size of (   l, w).
/dates           1D array of float32 in size of (n,     ) in ordinal date format
/recons          3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - filtered
/rawts           3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - un-filtered
/ifgcnt          2D array of int32   in size of (   l, w), number of interferograms used for every pixel
"""

class giantIfgramStack:
    """
    Time-series object for displacement of a set of SAR images from the same platform and track.
    """

    def __init__(self, file=None):
        self.file = file
        self.name = 'giantIfgramStack'

    def open(self, print_msg=True):
        if print_msg:
            print(f'open {self.name} file: {os.path.basename(self.file)}')
        self.get_size()
        self.get_date12_list()
        self.get_metadata()
        self.numPixel = self.length * self.width

        # Dataset Info
        with h5py.File(self.file, 'r') as f:
            self.pbaseIfgram = f['bperp'][:]
            # get existed datasetNames in the order of GIANT_DSET_NAMES
            dsNames = [i for i in f.keys()
                       if (isinstance(f[i], h5py.Dataset)
                           and f[i].shape[-2:] == (self.length, self.width))]
            self.datasetNames = [i for i in GIANT_DSET_NAMES if i in dsNames]
            self.datasetNames += [i for i in dsNames if i not in GIANT_DSET_NAMES]

            self.sliceList = []
            for dsName in self.datasetNames:
                ds = f[dsName]
                if len(ds.shape) == 3:
                    if ds.shape[0] == self.numIfgram:
                        self.sliceList += [f'{dsName}-{i}' for i in self.date12List]
                    elif ds.shape[0] == self.numDate:
                        self.sliceList += [f'{dsName}-{i}' for i in self.dateList]
                elif len(ds.shape) == 2:
                    self.sliceList.append(dsName)
                else:
                    raise ValueError(('un-recognized dataset dimension for {}:'
                                      ' {}').format(dsName, ds.shape))

    def get_size(self, dsName='igram'):
        with h5py.File(self.file, 'r') as f:
            dsName = [i for i in f.keys() if dsName in i][0]
            self.numIfgram, self.length, self.width = f[dsName].shape
        return self.numIfgram, self.length, self.width

    def get_date12_list(self):
        with h5py.File(self.file, 'r') as f:
            self.dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d')
                             for i in f['dates'][:].tolist()]
            self.numDate = len(self.dateList)
            # grab date12 from Jmat
            dates = np.array(self.dateList)
            Jmat = f['Jmat'][:]
            mDates = []
            sDates = []
            for i in range(Jmat.shape[0]):
                mDates.append(dates[Jmat[i, :] ==  1][0])
                sDates.append(dates[Jmat[i, :] == -1][0])
            self.date12List = [f'{m}_{s}' for m, s in zip(mDates, sDates)]
            self.mDates = mDates
            self.sDates = sDates
        return self.date12List

    def get_metadata(self):
        # metadata
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
            dateList = [dt.fromordinal(int(i)).strftime('%Y%m%d')
                        for i in f['dates'][:].tolist()]
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value
        self.metadata['START_DATE'] = dateList[0]
        self.metadata['END_DATE'] = dateList[-1]
        # size
        self.get_size()
        self.metadata['LENGTH'] = str(self.length)
        self.metadata['WIDTH'] = str(self.width)
        return self.metadata
