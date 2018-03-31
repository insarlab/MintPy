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
from pysar.utils import readfile, datetime as ptime
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
        /timeseries      3D array of float32 in size of (n, l, w) in meter.
        /date            1D array of string  in size of (n,) in YYYYMMDD format
        /bperp           1D array of float32 in size of (n,) in meter.
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'timeseries'

    def close(self):
        print('close timeseries: {}'.format(os.path.basename(self.file)))
        self.f.close()

    def open(self):
        print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.metadata = readfile.read_attribute(self.file)
        self.f = h5py.File(self.file,'r')
        self.numDate, self.length, self.width = self.f[self.key].get(self.key).shape
        dates = self.f[self.key].get('date')[:]
        self.times = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in dates])
        self.dateList = [i.decode('utf8') for i in dates]

        self.refIndex = self.dateList.index(self.metadata['REF_DATE'])
        self.btemp = np.array([i.days for i in self.times - self.times[self.refIndex]], dtype=np.int16)
        self.bperp = self.f[self.key].get('bperp')[:]
        self.bperp -= self.bperp[self.refIndex]

    def write2hdf5(self, data, outFile=None, dates=None, bperp=None, metadata=None, refFile=None):
        '''
        Parameters: data  : 3D array of float32
                    dates : 1D array/list of string in YYYYMMDD format
                    bperp : 1D array/list of float32
                    metadata : dict
                    outFile : string
                    refFile : string
        Returns: outFile : string
        Examples:
            from pysar.objects import timeseries
            tsobj = timeseries('timeseries.h5')
            timeseries.write(data, outFile='timeseries.h5', dates=dateList, bperp=bperp, metadata=atr)
            tsobj = timeseries('timeseries_demErr.h5')
            timeseries.write(data, outFile='timeseries_demErr.h5', refFile='timeseries.h5')
        '''
        if not outFile:
            outFile = self.file
        if refFile:
            refobj = timeseries(refFile)
            refobj.open()
            metadata = refobj.metadata
            dates = refobj.dateList
            bperp = refobj.bperp
            refobj.close()
        data = np.array(data, dtype=np.float32)
        dates = np.array(dates, dtype=np.string_)
        bperp = np.array(bperp, dtype=np.float32)

        ##### 3D and 1D dataset
        gName = 'timeseries'
        print('create timeseries HDF5 file: {} with w mode'.format(outFile))
        f = h5py.File(outFile,'w')
        group = f.create_group(gName)
        print('create dataset /{}/timeseries of {:<10} in size of {}'.format(gName, str(data.dtype), data.shape))
        dset = group.create_dataset('timeseries', data=data, chunks=True)

        print('create dataset /{}/dates      of {:<10} in size of {}'.format(gName, str(dates.dtype), dates.shape))
        dset = group.create_dataset('date', data=dates, chunks=True)

        print('create dataset /{}/bperp      of {:<10} in size of {}'.format(gName, str(bperp.dtype), bperp.shape))
        dset = group.create_dataset('bperp', data=bperp, chunks=True)

        #### Attributes
        for key, value in metadata.items():
            group.attrs[key] = str(value)

        f.close()
        print('finished writing to {}'.format(outFile))
        return outFile






########################################################################################
class ifgramStack:
    ''' Interferograms Stack object.
    /ifgramStack           Root level group name
        Attributes         Dictionary for metadata
        /date              2D array of string  in size of (m, 2   ) in YYYYMMDD format for master and slave date
        /bperp             1D array of float32 in size of (m,     ) in meter.
        /dropIfgram        1D array of bool    in size of (m,     ).
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
        '''
        Time format/rules:
            All datetime.datetime objects named with time
            All string in YYYYMMDD        named with date (following roipac)
        '''
        print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.metadata = readfile.read_attribute(self.file)
        self.f = h5py.File(self.file, 'r')
        self.get_size()
        self.read_datetimes()

        self.mDateList = [i.strftime("%Y%m%d") for i in self.mTimes]
        self.sDateList = [i.strftime("%Y%m%d") for i in self.sTimes]
        self.date12List = ['{}_{}'.format(i,j) for i,j in zip(self.mDateList,self.sDateList)]
        self.btemp = np.array([i.days for i in self.sTimes - self.mTimes], dtype=np.int16)
        self.bperp = self.f[self.key].get('bperp')[:][self.dropIfgram]

        #Time in timeseries domain
        self.dateList = sorted(list(set(np.hstack((self.mDateList, self.sDateList)))))
        self.numDate = len(self.dateList)
        self.btempHist = []
        d1 = dt(*time.strptime(self.dateList[0],"%Y%m%d")[0:5])
        for i in range(self.numDate):
            di = dt(*time.strptime(self.dateList[i],"%Y%m%d")[0:5])
            self.btempHist.append((di-d1).days)
        self.btempHistDiff = np.diff(self.btempHist)

    def close(self):
        print('close {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f.close()

    def get_size(self):
        self.length, self.width = self.f[self.key].get(ifgramDatasetNames[0]).shape[1:3]
        try:    self.dropIfgram = self.f[self.key].get('dropIfgram')[:]
        except: self.dropIfgram = np.ones(numIfgram, dtype=np.bool_)
        self.numIfgram = np.sum(self.dropIfgram)
        return self.numIfgram, self.length, self.width

    def read_datetimes(self):
        '''Read master/slave dates into array of datetime.datetime objects'''
        dates = self.f[self.key].get('date')
        self.mTimes = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in list(dates[:,0])])[self.dropIfgram]
        self.sTimes = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in list(dates[:,1])])[self.dropIfgram]
        return self.mTimes, self.sTimes     

    def read(self, datasetName=ifgramDatasetNames[0], box=None):
        '''Read 3D dataset with bounding box in space'''
        dset = self.f[self.key].get(datasetName)
        if box is None:
            box = (0,0,self.width,self.length)
        data = dset[self.dropIfgram, box[1]:box[3], box[0]:box[2]]
        return data


    ##### Functions for Network Inversion
    def get_design_matrix(self, refDate=None):
        '''Return design matrix of the input ifgramStack'''
        if not refDate:
            refDate = self.dateList[0]
        refIndex = self.dateList.index(refDate)

        ## calculate design matrix
        A = np.zeros((self.numIfgram, self.numDate))
        B = np.zeros(A.shape)
        for i in range(self.numIfgram):
            m_idx, s_idx = [self.dateList.index(j) for j in self.date12List[i].split('_')]
            A[i, m_idx] = -1
            A[i, s_idx] = 1
            B[i, m_idx:s_idx] = self.btemp[m_idx+1:s_idx+1] - self.btemp[m_idx:s_idx]
        #Remove reference date as it can not be resolved
        A = np.hstack((A[:,0:refIndex], A[:,(refIndex+1):]))
        B = B[:,:-1]
        return A, B

    def get_perp_baseline_timeseries(self):
        '''Get spatial perpendicular baseline in timeseries from input '''
        B = self.get_design_matrix()[1]
        B_inv = np.linalg.pinv(B)
        bperpRate = np.dot(B_inv, self.bperp)
        zero = np.array([0.],np.float32)
        pbaseTimeseries = np.concatenate((zero, np.cumsum([bperpRate*self.btempHistDiff])))
        return pbaseTimeseries












