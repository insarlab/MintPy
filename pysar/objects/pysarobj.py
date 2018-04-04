############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2018, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi, 2018              #
############################################################
# class used for file operation within PySAR
# Recommended usage:
#     from pysar.objects import timeseries, ifgramStack, geometry

import os, sys, glob
import time
from datetime import datetime as dt
import h5py
import numpy as np

dataType = np.float32

##------------------ Variables ---------------------##
timeseriesKeyNames = ['timeseries','HDFEOS','GIANT_TS']

timeseriesDatasetNames = ['raw',
                          'troposphericDelay',
                          'topographicResidual',
                          'ramp',
                          'displacement']

geometryDatasetNames = ['height',
                        'latitude',
                        'longitude',
                        'rangeCoord',
                        'azimuthCoord',
                        'incidenceAngle',
                        'headingAngle',
                        'slantRangeDistance',
                        'shadowMask',
                        'waterMask',
                        'commonMask',
                        'bperp']

ifgramDatasetNames = ['unwrapPhase',
                      'coherence',
                      'connectComponent',
                      'wrapPhase',
                      'iono',
                      'rangeOffset',
                      'azimuthOffset']



########################################################################################
class timeseries:
    '''
    Time-series object for displacement of a set of SAR images from the same platform and track.
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.

    /                    Root level
    timeseries           group name
        Attributes       Dictionary for metadata
        /timeseries      3D array of float32 in size of (n, l, w) in meter.
        /date            1D array of string  in size of (n,     ) in YYYYMMDD format
        /bperp           1D array of float32 in size of (n,     ) in meter. (optional)
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'timeseries'

    def close(self, printMsg=True):
        try:
            self.f.close()
            if printMsg:
                print('close timeseries file: {}'.format(os.path.basename(self.file)))
        except:
            pass

    def open(self, printMsg=True):
        if printMsg:
            print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f = h5py.File(self.file,'r')
        self.get_metadata()
        self.numDate, self.length, self.width = self.f[self.key].get(self.key).shape
        self.numPixel = self.length * self.width

        dates = self.f[self.key].get('date')[:]
        self.times = np.array([dt(*time.strptime(i.decode('utf8'),"%Y%m%d")[0:5]) for i in dates])
        self.dateList = [i.decode('utf8') for i in dates]
        self.datasetList = list(self.dateList)

        #Temporal baseline in days
        if 'REF_DATE' in self.metadata.keys():
            self.refIndex = self.dateList.index(self.metadata['REF_DATE'])
            self.btemp = np.array([i.days for i in self.times - self.times[self.refIndex]], dtype=np.int16)
        else:
            self.refIndex = None

        #Perpendicular baseline in meters
        if 'bperp' in self.f[self.key].keys():
            self.bperp = self.f[self.key].get('bperp')[:]
            if self.refIndex:
                self.bperp -= self.bperp[self.refIndex]
        else:
            self.bperp = None

    def get_metadata(self):
        self.f = h5py.File(self.file, 'r')
        self.metadata = dict(self.f[self.key].attrs)
        for key, value in self.metadata.items():
            try:     self.metadata[key] = value.decode('utf8')
            except:  self.metadata[key] = value
        return self.metadata

    def read(self, datasetName=None, box=None, printMsg=True):
        '''Read dataset from timeseries file
        Parameters: self : timeseries object
                    datasetName : (list of) string in YYYYMMDD format
                    box : tuple of 4 int, indicating x0,y0,x1,y1 of range
        Returns:    data : 2D or 3D dataset
        Examples:   from pysar.objects import timeseries
                    tsobj = timeseries('timeseries_ECMWF_demErr.h5')
                    data = tsobj.read(datasetName='20161020')
                    data = tsobj.read(datasetName='20161020', box=(100,300,500,800))
                    data = tsobj.read(datasetName=['20161020','20161026','20161101'])
                    data = tsobj.read(box=(100,300,500,800))
        '''
        try: self.f
        except: self.open(printMsg=printMsg)
        ds = self.f[self.key].get(self.key)
        ##Index in time/1st dimension
        if not datasetName:
            dsIndex = range(self.numDate)
        elif isinstance(datasetName, str):
            dsIndex = self.dateList.index(datasetName)
        elif isinstance(datasetName, list):
            dsIndex = []
            for e in datasetName:
                dsIndex.append(self.dateList.index(e))
        ##Index in space/2_3 dimension
        if box is None:
            box = [0,0,self.width,self.length]
        data = ds[dsIndex, box[1]:box[3], box[0]:box[2]]
        data = np.squeeze(data)
        return data


    def write2hdf5(self, data, outFile=None, dates=None, bperp=None, metadata=None, refFile=None):
        '''
        Parameters: data  : 3D array of float32
                    dates : 1D array/list of string in YYYYMMDD format
                    bperp : 1D array/list of float32 (optional)
                    metadata : dict
                    outFile : string
                    refFile : string
        Returns: outFile : string
        Examples:
            from pysar.objects import timeseries

            ##Generate a new timeseries file
            tsobj = timeseries('timeseries.h5')
            timeseries.write(data, dates=dateList, bperp=bperp, metadata=atr)

            ##Generate a timeseries with same attributes and same date/bperp info
            tsobj = timeseries('timeseries_demErr.h5')
            timeseries.write(data, refFile='timeseries.h5')
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

        if bperp.shape != ():
            print('create dataset /{}/bperp      of {:<10} in size of {}'.format(gName, str(bperp.dtype), bperp.shape))
            dset = group.create_dataset('bperp', data=bperp, chunks=True)

        #### Attributes
        for key, value in metadata.items():
            group.attrs[key] = str(value)

        f.close()
        print('finished writing to {}'.format(outFile))
        return outFile




########################################################################################
class geometry:
    ''' Geometry object.
    /geometry                    Root level group name
        Attributes               Dictionary for metadata. 'X/Y_FIRST/STEP' attribute for geocoded.
        /height                  2D array of float32 in size of (l, w   ) in meter.
        /latitude (azimuthCoord) 2D array of float32 in size of (l, w   ) in degree.
        /longitude (rangeCoord)  2D array of float32 in size of (l, w   ) in degree.
        /incidenceAngle          2D array of float32 in size of (l, w   ) in degree.
        /slantRangeDistance      2D array of float32 in size of (l, w   ) in meter.
        /headingAngle            2D array of float32 in size of (l, w   ) in degree. (optional)
        /shadowMask              2D array of bool    in size of (l, w   ).           (optional)
        /waterMask               2D array of bool    in size of (l, w   ).           (optional)
        /bperp                   3D array of float32 in size of (n, l, w) in meter   (optional)
        /date                    1D array of string  in size of (n,     ) in YYYYMMDD(optional)
        ...
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'geometry'

    def close(self, printMsg=True):
        try:
            self.f.close()
            if printMsg:
                print('close geometry file: {}'.format(os.path.basename(self.file)))
        except: 
            pass

    def open(self, printMsg=True):
        if printMsg:
            print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f = h5py.File(self.file,'r')
        self.get_metadata()
        self.length, self.width = self.f[self.key].get(geometryDatasetNames[0]).shape
        self.numPixel = self.length * self.width

        self.datasetList = list(set(self.f[self.key].keys()) & set(geometryDatasetNames))
        if 'bperp' in self.f[self.key].keys():
            self.dateList = [i.decode('utf8') for i in self.f[self.key].get('date')[:]]
            self.numDate = len(self.dateList)
            ##Update bperp datasetNames
            try: self.datasetList.remove('bperp')
            except: pass
            self.datasetList += ['bperp-'+d for d in self.dateList]
        else:
            self.dateList = None

        self.geocoded = False
        if 'Y_FIRST' in self.metadata.keys():
            self.geocoded = True

    def get_metadata(self):
        self.f = h5py.File(self.file, 'r')
        self.metadata = dict(self.f[self.key].attrs)
        for key, value in self.metadata.items():
            try:     self.metadata[key] = value.decode('utf8')
            except:  self.metadata[key] = value
        return self.metadata

    def read(self, datasetName=geometryDatasetNames[0], box=None, printMsg=True):
        '''Read 2D / 3D dataset with bounding box in space
        Parameters: datasetName : string, to point to specific 2D dataset, e.g.:
                        height
                        incidenceAngle
                        bperp
                        ...
                        bperp-20161020
                        bperp-20161026
                        bperp-...
                    box : tuple of 4 int, for (x0,y0,x1,y1)
                    printMsg : bool
        Returns: data : 2D or 3D array
        Example:
            obj = geometry('./INPUTS/geometryRadar.h5')
            obj.read(datasetName='height')
            obj.read(datasetName='incidenceAngle')
            obj.read(datasetName='bperp')
            obj.read(datasetName='bperp-20161020')
        '''
        try: self.f
        except: self.open(printMsg=printMsg)
        if box is None:
            box = (0,0,self.width,self.length)
        if datasetName is None:
            datasetName = geometryDatasetNames[0]

        datasetName = datasetName.split('-')
        dset = self.f[self.key].get(datasetName[0])
        if len(dset.shape) == 2:
            data = dset[box[1]:box[3], box[0]:box[2]]

        ## bperp
        elif len(dset.shape) == 3:
            if len(datasetName) == 1:
                data = dset[:, box[1]:box[3], box[0]:box[2]]
            else:
                data = dset[self.dateList.index(datasetName[1]), box[1]:box[3], box[0]:box[2]]
                data = np.squeeze(data)
        return data







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

    def close(self, printMsg=True):
        try:
            self.f.close()
            if printMsg:
                print('close {} file: {}'.format(self.key, os.path.basename(self.file)))
        except:
            pass

    def open(self, printMsg=True):
        '''
        Time format/rules:
            All datetime.datetime objects named with time
            All string in YYYYMMDD        named with date (following roipac)
        '''
        if printMsg:
            print('open {} file: {}'.format(self.key, os.path.basename(self.file)))
        self.f = h5py.File(self.file, 'r')
        self.get_metadata()
        self.get_size()
        self.read_datetimes()
        self.numPixel = self.length * self.width

        #Get datasetList for self.read()
        self.datasetList = []
        mDateList = [i.decode('utf8') for i in self.f[self.key].get('date')[:,0]]
        sDateList = [i.decode('utf8') for i in self.f[self.key].get('date')[:,1]]
        for dsName in list(set(self.f[self.key].keys()) & set(ifgramDatasetNames)):
            self.datasetList += ['{}-{}_{}'.format(dsName,m,s) for m,s in zip(mDateList, sDateList)]

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

    def get_metadata(self):
        self.f = h5py.File(self.file, 'r')
        self.metadata = dict(self.f[self.key].attrs)
        for key, value in self.metadata.items():
            try:     self.metadata[key] = value.decode('utf8')
            except:  self.metadata[key] = value
        return self.metadata

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

    def read(self, datasetName=ifgramDatasetNames[0], box=None, printMsg=True):
        '''Read 3D dataset with bounding box in space
        Parameters: datasetName : string, to point to specific 2D dataset, e.g.:
                        unwrapPhase
                        coherence
                        connectComponent
                        ...
                        unwrapPhase-20161020_20161026
                        unwrapPhase-...
                        coherence-20161020_20161026
                        ...
                    box : tuple of 4 int, for (x0,y0,x1,y1)
                    printMsg : bool
        Returns: data : 2D or 3D array
        Example:
            obj = ifgramStack('./INPUTS/ifgramStack.h5')
            obj.read(datasetName='unwrapPhase')
            obj.read(datasetName='coherence')
            obj.read(datasetName='unwrapPhase-20161020_20161026')
        '''
        try: self.f
        except: self.open(printMsg=printMsg)
        if box is None:
            box = (0,0,self.width,self.length)
        if datasetName is None:
            datasetName = ifgramDatasetNames[0]

        datasetName = datasetName.split('-')
        dset = self.f[self.key].get(datasetName[0])
        if len(datasetName) == 1:
            data = dset[self.dropIfgram, box[1]:box[3], box[0]:box[2]]
        else:
            data = dset[self.date12List.index(datasetName[1]), box[1]:box[3], box[0]:box[2]]
            data = np.squeeze(data)
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






########################################################################################
class hdfEos5:
    '''
    Time-series object in HDF-EOS5 format for Univ of Miami's InSAR Time-series Web Viewer (http://insarmaps.miami.edu)
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.

    /HDFEOS                           Root level group name
        /GRIDS                        2nd level group name for products in geo coordinates
            Attributes                metadata in dict.
            /timeseries               3rd level group name for time-series InSAR product
                /date                 1D array of string  in size of (n,     ) in YYYYMMDD format.
                /bperp                1D array of float32 in size of (n,     ) in meter
                /temporalCoherence    2D array of float32 in size of (   l, w).
                /mask                 2D array of bool_   in size of (   l, w).
                /raw                  3D array of float32 in size of (n, l, w) in meter
                /troposphericDelay    3D array of float32 in size of (n, l, w) in meter
                /topographicResidual  3D array of float32 in size of (n, l, w) in meter
                /ramp                 3D array of float32 in size of (n, l, w) in meter
                /displacement         3D array of float32 in size of (n, l, w) in meter
            /geometry                 3rd level group name for geometry data
                /height               2D array of float32 in size of (   l, w) in meter.
                /incidenceAngle       2D array of float32 in size of (   l, w) in degree.
                /slantRangeDistance   2D array of float32 in size of (   l, w) in meter.
                /azimuthCoord         2D array of float32 in size of (   l, w) in degree.
                /rangeCoord           2D array of float32 in size of (   l, w) in degree.
                /headingAngle         2D array of float32 in size of (   l, w) in degree. (optional)
                /shadowMask           2D array of bool    in size of (   l, w).           (optional)
                /waterMask            2D array of bool    in size of (   l, w).           (optional)
                /bperp                3D array of float32 in size of (n, l, w) in meter.  (optional)
        /SWATHS
            Attributes
            /ifgramStack
                ...
            /geometry
                ...
    '''
    def __init__(self, file=None):
        self.file = file
        self.key = 'HDFEOS'

    def close(self, printMsg=True):
        try:
            self.f.close()
            if printMsg:
                print('close timeseries file: {}'.format(os.path.basename(self.file)))
        except:
            pass









