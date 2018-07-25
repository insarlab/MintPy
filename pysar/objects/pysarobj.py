############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi, 2018              #
############################################################
# class used for file operation within PySAR
# Recommend import:
#     from pysar.objects import timeseries, ifgramStack, geometry


import os
import sys
import glob
import time
from datetime import datetime as dt
import h5py
import numpy as np

BOOL_ZERO = np.bool_(0)
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)

dataType = np.float32

dataTypeDict = {'bool': np.bool_, 'byte': np.bool_, 'flag': np.bool_,
                'int': np.int16, 'int16': np.int16, 'short': np.int16, 'int32': np.int32,
                'int64': np.int64, 'long': np.int64,
                'float': np.float32, 'float32': np.float32,
                'float_': np.float64, 'float64': np.float64,
                'complex': np.complex64, 'complex64': np.complex64, 'cpx_float32': np.complex64,
                'cfloat': np.complex64, 'cfloat32': np.complex64,
                'complex128': np.complex128, 'complex_': np.complex128, 'cpx_float64': np.complex128
                }

##------------------ Variables ---------------------##
timeseriesKeyNames = ['timeseries', 'HDFEOS', 'giantTimeseries']

timeseriesDatasetNames = ['timeseries',
                          'raw',
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
                        'azimuthAngle',
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
                      'azimuthOffset',
                      'refPhase']

datasetUnitDict = {'unwrapPhase'        :'radian',
                   'coherence'          :'1',
                   'connectComponent'   :'1',
                   'wrapPhase'          :'radian',
                   'iono'               :'radian',  #not sure
                   'rangeOffset'        :'1',
                   'azimuthOffset'      :'1',

                   'height'             :'m',
                   'latitude'           :'degree',
                   'longitude'          :'degree',
                   'rangeCoord'         :'1',
                   'azimuthCoord'       :'1',
                   'incidenceAngle'     :'degree',
                   'azimuthAngle'       :'degree',
                   'slantRangeDistance' :'m',
                   'shadowMask'         :'1',
                   'waterMask'          :'1',
                   'commonMask'         :'1',
                   'bperp'              :'m',

                   'timeseries'         :'m',
                   'giantTimeseries'    :'mm',
                   'raw'                :'m',
                   'troposphericDelay'  :'m',
                   'topographicResidual':'m',
                   'ramp'               :'m',
                   'displacement'       :'m',

                   'velocity'           :'m/year',
                   'acceleration'       :'m/year^2',
                   'mask'               :'1',

                   '.unw'               :'radian',
                   '.int'               :'radian',
                   '.flat'              :'radian',
                   '.cor'               :'1',
                   '.dem'               :'m',
                   '.hgt'               :'m',
                   '.hgt_sim'           :'m',
                   }


########################################################################################
FILE_STRUCTURE_TIMESERIES = """
/                Root level
Attributes       Dictionary for metadata
/timeseries      3D array of float32 in size of (n, l, w) in meter.
/date            1D array of string  in size of (n,     ) in YYYYMMDD format
/bperp           1D array of float32 in size of (n,     ) in meter. (optional)
"""

class timeseries:
    """
    Time-series object for displacement of a set of SAR images from the same platform and track.
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.
    """

    def __init__(self, file=None):
        self.file = file
        self.name = 'timeseries'
        self.file_structure = FILE_STRUCTURE_TIMESERIES

    def close(self, print_msg=True):
        try:
            self.f.close()
            if print_msg:
                print('close timeseries file: {}'.format(os.path.basename(self.file)))
        except:
            pass

    def open(self, print_msg=True):
        if print_msg:
            print('open {} file: {}'.format(self.name, os.path.basename(self.file)))
        self.get_metadata()
        self.get_size()
        self.get_date_list()
        self.numPixel = self.length * self.width

        with h5py.File(self.file, 'r') as f:
            try:
                self.pbase = f['bperp'][:]
                self.pbase -= self.pbase[self.refIndex]
            except:
                self.pbase = None
        self.times = np.array([dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in self.dateList])
        self.tbase = np.array([i.days for i in self.times - self.times[self.refIndex]],
                              dtype=np.float32)
        # list of float for year, 2014.95
        self.yearList = [i.year + (i.timetuple().tm_yday-1)/365.25 for i in self.times]
        self.datasetList = ['{}-{}'.format(self.name, i) for i in self.dateList]

    def get_metadata(self):
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
            dates = f['date'][:]
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value
        # ref_date/index
        dateList = [i.decode('utf8') for i in dates]
        if 'REF_DATE' not in self.metadata.keys():
            self.metadata['REF_DATE'] = dateList[0]
        self.refIndex = dateList.index(self.metadata['REF_DATE'])
        self.metadata['START_DATE'] = dateList[0]
        self.metadata['END_DATE'] = dateList[-1]
        return self.metadata

    def get_size(self):
        with h5py.File(self.file, 'r') as f:
            self.numDate, self.length, self.width = f[self.name].shape
        return self.numDate, self.length, self.width

    def get_date_list(self):
        with h5py.File(self.file, 'r') as f:
            self.dateList = [i.decode('utf8') for i in f['date'][:]]
        return self.dateList

    def read(self, datasetName=None, box=None, print_msg=True):
        """Read dataset from timeseries file
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
        """
        if print_msg:
            print('reading {} data from file: {} ...'.format(self.name, self.file))
        self.open(print_msg=False)

        # convert input datasetName into list of dates
        if not datasetName or datasetName == 'timeseries':
            datasetName = []
        elif isinstance(datasetName, str):
            datasetName = [datasetName]
        datasetName = [i.replace('timeseries', '').replace('-', '') for i in datasetName]

        with h5py.File(self.file, 'r') as f:
            ds = f[self.name]
            if isinstance(ds, h5py.Group):  # support for old pysar files
                ds = ds[self.name]

            # Get dateFlag - mark in time/1st dimension
            dateFlag = np.zeros((self.numDate), dtype=np.bool_)
            if not datasetName:
                dateFlag[:] = True
            else:
                for e in datasetName:
                    dateFlag[self.dateList.index(e)] = True

            # Get Index in space/2_3 dimension
            if box is None:
                box = [0, 0, self.width, self.length]

            data = ds[dateFlag, box[1]:box[3], box[0]:box[2]]
            data = np.squeeze(data)
        return data

    def write2hdf5(self, data, outFile=None, dates=None, bperp=None, metadata=None, refFile=None, compression=None):
        """
        Parameters: data  : 3D array of float32
                    dates : 1D array/list of string in YYYYMMDD format
                    bperp : 1D array/list of float32 (optional)
                    metadata : dict
                    outFile : string
                    refFile : string
                    compression : string or None
        Returns: outFile : string
        Examples:
            from pysar.objects import timeseries

            ##Generate a new timeseries file
            tsobj = timeseries('timeseries.h5')
            timeseries.write(data, dates=dateList, bperp=bperp, metadata=atr)

            ##Generate a timeseries with same attributes and same date/bperp info
            tsobj = timeseries('timeseries_demErr.h5')
            timeseries.write(data, refFile='timeseries.h5')
        """
        if not outFile:
            outFile = self.file
        if refFile:
            refobj = timeseries(refFile)
            refobj.open(print_msg=False)
            if metadata is None:
                metadata = refobj.metadata
            if dates is None:
                dates = refobj.dateList
            if bperp is None:
                bperp = refobj.pbase
            refobj.close(print_msg=False)
        data = np.array(data, dtype=np.float32)
        dates = np.array(dates, dtype=np.string_)
        bperp = np.array(bperp, dtype=np.float32)
        metadata['FILE_TYPE'] = self.name

        # 3D dataset - timeseries
        print('create timeseries HDF5 file: {} with w mode'.format(outFile))
        f = h5py.File(outFile, 'w')
        print(('create dataset /timeseries of {t:<10} in size of {s} '
               'with compression={c}').format(t=str(data.dtype),
                                              s=data.shape,
                                              c=compression))
        dset = f.create_dataset('timeseries', data=data, chunks=True, compression=compression)

        # 1D dataset - date / bperp
        print('create dataset /dates      of {:<10} in size of {}'.format(str(dates.dtype), dates.shape))
        dset = f.create_dataset('date', data=dates)

        if bperp.shape != ():
            print('create dataset /bperp      of {:<10} in size of {}'.format(str(bperp.dtype), bperp.shape))
            dset = f.create_dataset('bperp', data=bperp)

        # Attributes
        for key, value in metadata.items():
            f.attrs[key] = str(value)

        f.close()
        print('finished writing to {}'.format(outFile))
        return outFile

    def timeseries_std(self, maskFile=None, outFile=None):
        """Calculate the standard deviation (STD) for each epoch,
           output result to a text file.
        """
        # Calculate STD
        data = self.read()
        if maskFile:
            mask = singleDataset(maskFile).read()
            print('read mask from file: '+maskFile)
            data[:, mask == 0] = np.nan
        self.std = np.nanstd(data, axis=(1, 2))

        # Write text file
        header = 'Standard Deviation in space for each epoch of timeseries\n'
        header += 'Timeseries file: {}\n'.format(self.file)
        header += 'Mask file: {}\n'.format(maskFile)
        header += 'Date\t\tSTD (m)'
        if not outFile:
            outFile = os.path.join(os.path.dirname(os.path.abspath(self.file)),
                                   'std_{}.txt'.format(os.path.splitext(os.path.basename(self.file))[0]))
        np.savetxt(outFile, np.hstack((np.array(self.dateList).reshape(-1, 1), self.std.reshape(-1, 1))),
                   fmt='%s', delimiter='\t', header=header)
        print('save timeseries STD to text file: {}'.format(outFile))
        return outFile

    def timeseries_rms(self, maskFile=None, outFile=None):
        """Calculate the Root Mean Square for each epoch of input timeseries file
            and output result to a text file.
        """
        # Calculate RMS
        data = self.read()
        if maskFile and os.path.isfile(maskFile):
            print('read mask from file: '+maskFile)
            mask = singleDataset(maskFile).read()
            data[:, mask == 0] = np.nan
        self.rms = np.sqrt(np.nanmean(np.square(data), axis=(1, 2)))

        # Write text file
        header = 'Root Mean Square in space for each epoch of timeseries\n'
        header += 'Timeseries file: {}\n'.format(self.file)
        header += 'Mask file: {}\n'.format(maskFile)
        header += 'Date\t\tRMS (m)'
        if not outFile:
            outFile = os.path.join(os.path.dirname(os.path.abspath(self.file)),
                                   'rms_{}.txt'.format(os.path.splitext(os.path.basename(self.file))[0]))
        np.savetxt(outFile, np.hstack((np.array(self.dateList).reshape(-1, 1), self.rms.reshape(-1, 1))),
                   fmt='%s', delimiter='\t', header=header)
        print('save timeseries RMS to text file: {}'.format(outFile))
        return outFile

    def spatial_average(self, maskFile=None, box=None):
        self.open(print_msg=False)
        data = self.read(box=box)
        if maskFile and os.path.isfile(maskFile):
            print('read mask from file: '+maskFile)
            mask = singleDataset(maskFile).read(box=box)
            data[:, mask == 0] = np.nan
        dmean = np.nanmean(data, axis=(1, 2))
        return dmean, self.dateList

    def temporal_average(self):
        print('calculating the temporal average of timeseries file: {}'.format(self.file))
        self.open(print_msg=False)
        data = self.read()
        dmean = np.nanmean(data, axis=0)
        return dmean

    def save2bl_list_file(self, out_file='bl_list.txt'):
        """Generate bl_list.txt file from timeseries h5 file."""
        self.open(print_msg=False)
        date6_list = [i[2:8] for i in self.dateList]
        pbase_list = self.pbase.tolist()
        print('write baseline list info to file: {}'.format(out_file))
        with open(out_file, 'w') as f:
            for d, pbase in zip(date6_list, pbase_list):
                f.write('{}\t{}\n'.format(d, pbase))
        return out_file


########################################################################################
FILE_STRUCTURE_GEOMETRY = """
/                        Root level
Attributes               Dictionary for metadata. 'X/Y_FIRST/STEP' attribute for geocoded.
/height                  2D array of float32 in size of (l, w   ) in meter.
/latitude (azimuthCoord) 2D array of float32 in size of (l, w   ) in degree.
/longitude (rangeCoord)  2D array of float32 in size of (l, w   ) in degree.
/incidenceAngle          2D array of float32 in size of (l, w   ) in degree.
/slantRangeDistance      2D array of float32 in size of (l, w   ) in meter.
/azimuthAngle            2D array of float32 in size of (l, w   ) in degree. (optional)
/shadowMask              2D array of bool    in size of (l, w   ).           (optional)
/waterMask               2D array of bool    in size of (l, w   ).           (optional)
/bperp                   3D array of float32 in size of (n, l, w) in meter   (optional)
/date                    1D array of string  in size of (n,     ) in YYYYMMDD(optional)
...
"""

class geometry:
    """ Geometry object."""

    def __init__(self, file=None):
        self.file = file
        self.name = 'geometry'
        self.file_structure = FILE_STRUCTURE_GEOMETRY

    def close(self, print_msg=True):
        try:
            self.f.close()
            if print_msg:
                print('close geometry file: {}'.format(os.path.basename(self.file)))
        except:
            pass

    def open(self, print_msg=True):
        if print_msg:
            print('open {} file: {}'.format(self.name, os.path.basename(self.file)))
        self.get_metadata()
        self.get_size()
        self.numPixel = self.length * self.width
        self.geocoded = False
        if 'Y_FIRST' in self.metadata.keys():
            self.geocoded = True

        with h5py.File(self.file, 'r') as f:
            self.datasetNames = [i for i in geometryDatasetNames if i in f.keys()]
            self.datasetList = list(self.datasetNames)
            if 'bperp' in f.keys():
                self.dateList = [i.decode('utf8') for i in f['date'][:]]
                self.numDate = len(self.dateList)
                # Update bperp datasetNames
                try:
                    self.datasetList.remove('bperp')
                except:
                    pass
                self.datasetList += ['bperp-'+d for d in self.dateList]
            else:
                self.dateList = None

    def get_size(self):
        with h5py.File(self.file, 'r') as f:
            dsName = [i for i in f.keys() if i in geometryDatasetNames][0]
            dsShape = f[dsName].shape
            if len(dsShape) == 3:
                self.length, self.width = dsShape[1:3]
            else:
                self.length, self.width = dsShape
        return self.length, self.width

    def get_metadata(self):
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value
        return self.metadata

    def read(self, datasetName=geometryDatasetNames[0], box=None, print_msg=True):
        """Read 2D / 3D dataset with bounding box in space
        Parameters: datasetName : (list of) string, to point to specific 2D dataset, e.g.:
                        height
                        incidenceAngle
                        bperp
                        ...
                        bperp-20161020
                        bperp-20161026
                        bperp-...
                    box : tuple of 4 int, for (x0,y0,x1,y1)
                    print_msg : bool
        Returns: data : 2D or 3D array
        Example:
            obj = geometry('./INPUTS/geometryRadar.h5')
            obj.read(datasetName='height')
            obj.read(datasetName='incidenceAngle')
            obj.read(datasetName='bperp')
            obj.read(datasetName='bperp-20161020')
            obj.read(datasetName=['bperp-20161020',
                                  'bperp-20161026'])
        """
        self.open(print_msg=False)
        if box is None:
            box = (0, 0, self.width, self.length)

        if datasetName is None:
            datasetName = geometryDatasetNames[0]
        elif isinstance(datasetName, str):
            datasetName = [datasetName]

        with h5py.File(self.file, 'r') as f:
            familyName = datasetName[0].split('-')[0]
            ds = f[familyName]
            if print_msg:
                print('reading {} data from file: {} ...'.format(familyName, self.file))

            if len(ds.shape) == 1:
                data = ds[:]
            elif len(ds.shape) == 2:
                data = ds[box[1]:box[3], box[0]:box[2]]
            else:
                # get dateFlag - mark in time/1st dimension
                dateFlag = np.zeros((ds.shape[0]), dtype=np.bool_)
                datasetName = [i.replace(familyName, '').replace('-', '') for i in datasetName]
                if any(not i for i in datasetName):
                    dateFlag[:] = True
                else:
                    for e in datasetName:
                        dateFlag[self.dateList.index(e)] = True
                data = ds[dateFlag, box[1]:box[3], box[0]:box[2]]
                data = np.squeeze(data)
        return data


########################################################################################
FILE_STRUCTURE_IFGRMA_STACK = """
/                  Root level group name
Attributes         Dictionary for metadata
/date              2D array of string  in size of (m, 2   ) in YYYYMMDD format for master and slave date
/bperp             1D array of float32 in size of (m,     ) in meter.
/dropIfgram        1D array of bool    in size of (m,     ) with 0 for drop and 1 for keep
/unwrapPhase       3D array of float32 in size of (m, l, w) in radian.
/coherence         3D array of float32 in size of (m, l, w).
/connectComponent  3D array of int16   in size of (m, l, w).           (optional)
/wrapPhase         3D array of float32 in size of (m, l, w) in radian. (optional)
/rangeOffset       3D array of float32 in size of (m, l, w).           (optional)
/azimuthOffset     3D array of float32 in size of (m, l, w).           (optional)
"""

class ifgramStack:
    """ Interferograms Stack object."""

    def __init__(self, file=None):
        self.file = file
        self.name = 'ifgramStack'
        self.file_structure = FILE_STRUCTURE_IFGRMA_STACK

    def close(self, print_msg=True):
        try:
            self.f.close()
            if print_msg:
                print('close {} file: {}'.format(self.name, os.path.basename(self.file)))
        except:
            pass

    def open(self, print_msg=True):
        """
        Time format/rules:
            All datetime.datetime objects named with time
            All string in YYYYMMDD        named with date (following roipac)
        """
        if print_msg:
            print('open {} file: {}'.format(self.name, os.path.basename(self.file)))
        self.get_metadata()
        self.get_size()
        self.read_datetimes()
        self.numPixel = self.length * self.width

        with h5py.File(self.file, 'r') as f:
            self.dropIfgram = f['dropIfgram'][:]
            self.pbaseIfgram = f['bperp'][:]
            self.datasetNames = [key for key in list(f.keys())
                                 if (isinstance(f[key], h5py.Dataset)
                                     and f[key].shape[-2:] == (self.length, self.width))]
        self.date12List = ['{}_{}'.format(i, j) for i, j in zip(self.mDates, self.sDates)]
        self.tbaseIfgram = np.array([i.days for i in self.sTimes - self.mTimes], dtype=np.float32)

        # Get datasetList for self.read()
        self.datasetList = []
        for dsName in self.datasetNames:
            self.datasetList += ['{}-{}'.format(dsName, i) for i in self.date12List]

        # Time in timeseries domain
        self.dateList = self.get_date_list(dropIfgram=False)
        self.numDate = len(self.dateList)

        # Reference pixel
        try:
            self.refY = int(self.metadata['REF_Y'])
            self.refX = int(self.metadata['REF_X'])
        except:
            self.refY = None
            self.refX = None
        try:
            self.refLat = float(self.metadata['REF_LAT'])
            self.refLon = float(self.metadata['REF_LON'])
        except:
            self.refLat = None
            self.refLon = None

    def get_metadata(self):
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
            dates = f['date'][:].flatten()
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value
        dateList = sorted([i.decode('utf8') for i in dates])
        self.metadata['START_DATE'] = dateList[0]
        self.metadata['END_DATE'] = dateList[-1]
        return self.metadata

    def get_size(self, dropIfgram=False):
        with h5py.File(self.file, 'r') as f:
            self.numIfgram, self.length, self.width = f[ifgramDatasetNames[0]].shape
        return self.numIfgram, self.length, self.width

    def read_datetimes(self):
        """Read master/slave dates into array of datetime.datetime objects"""
        with h5py.File(self.file, 'r') as f:
            dates = f['date'][:]
        self.mDates = np.array([i.decode('utf8') for i in dates[:, 0]])
        self.sDates = np.array([i.decode('utf8') for i in dates[:, 1]])
        self.mTimes = np.array([dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in self.mDates])
        self.sTimes = np.array([dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in self.sDates])

    def read(self, datasetName=ifgramDatasetNames[0], box=None, print_msg=True, dropIfgram=False):
        """Read 3D dataset with bounding box in space
        Parameters: datasetName : string, to point to specific 2D dataset, e.g.:
                        unwrapPhase
                        coherence
                        connectComponent
                        ...
                        unwrapPhase-20161020_20161026
                        unwrapPhase-...
                        coherence-20161020_20161026
                        ...
                        ['unwrapPhase-20161020_20161026',
                         'unwrapPhase-20161020_20161101',
                         ...]
                    box : tuple of 4 int, for (x0,y0,x1,y1)
                    print_msg : bool
        Returns: data : 2D or 3D array
        Example:
            obj = ifgramStack('./INPUTS/ifgramStack.h5')
            obj.read(datasetName='unwrapPhase')
            obj.read(datasetName='coherence')
            obj.read(datasetName='unwrapPhase-20161020_20161026')
            obj.read(datasetName=['unwrapPhase-20161020_20161026',
                                  'unwrapPhase-20161020_20161101'])
        """
        self.get_size()
        date12List = self.get_date12_list(dropIfgram=False)

        # convert input datasetName into list
        if datasetName is None:
            datasetName = [ifgramDatasetNames[0]]
        elif isinstance(datasetName, str):
            datasetName = [datasetName]

        with h5py.File(self.file, 'r') as f:
            familyName = datasetName[0].split('-')[0]
            ds = f[familyName]
            if print_msg:
                print('reading {} data from file: {} ...'.format(familyName, self.file))

            # get dateFlag - mark in time/1st dimension
            dateFlag = np.zeros((self.numIfgram), dtype=np.bool_)
            datasetName = [i.replace(familyName, '').replace('-', '') for i in datasetName]
            if any(not i for i in datasetName):
                if dropIfgram:
                    dateFlag = f['dropIfgram'][:]
                else:
                    dateFlag[:] = True
            else:
                for e in datasetName:
                    dateFlag[date12List.index(e)] = True

            # get index in space/2-3 dimension
            if box is None:
                box = (0, 0, self.width, self.length)

            data = ds[dateFlag, box[1]:box[3], box[0]:box[2]]
            data = np.squeeze(data)
        return data

    def spatial_average(self, datasetName=ifgramDatasetNames[1], maskFile=None, box=None):
        if datasetName is None:
            datasetName = ifgramDatasetNames[1]
        print('calculating spatial average of {} in file {} ...'.format(datasetName, self.file))
        if maskFile and os.path.isfile(maskFile):
            print('read mask from file: '+maskFile)
            mask = singleDataset(maskFile).read(box=box)
        else:
            maskFile = None

        with h5py.File(self.file, 'r') as f:
            dset = f[datasetName]
            numIfgram = dset.shape[0]
            dmean = np.zeros((numIfgram), dtype=np.float32)
            for i in range(numIfgram):
                data = dset[i, box[1]:box[3], box[0]:box[2]]
                if maskFile:
                    data[mask == 0] = np.nan
                dmean[i] = np.nanmean(data)
                sys.stdout.write('\rreading interferogram {}/{} ...'.format(i+1, numIfgram))
                sys.stdout.flush()
            print('')
        return dmean, self.date12List

    # Functions considering dropIfgram value

    def get_date12_list(self, dropIfgram=True):
        with h5py.File(self.file, 'r') as f:
            dates = f['date'][:]
            if dropIfgram:
                dates = dates[f['dropIfgram'][:], :]
        mDates = np.array([i.decode('utf8') for i in dates[:, 0]])
        sDates = np.array([i.decode('utf8') for i in dates[:, 1]])
        date12List = ['{}_{}'.format(i, j) for i, j in zip(mDates, sDates)]
        return date12List

    def get_drop_date12_list(self):
        with h5py.File(self.file, 'r') as f:
            dates = f['date'][:]
            dates = dates[~f['dropIfgram'][:], :]
        mDates = np.array([i.decode('utf8') for i in dates[:, 0]])
        sDates = np.array([i.decode('utf8') for i in dates[:, 1]])
        date12List = ['{}_{}'.format(i, j) for i, j in zip(mDates, sDates)]
        return date12List

    def get_date_list(self, dropIfgram=False):
        with h5py.File(self.file, 'r') as f:
            dates = f['date'][:]
            if dropIfgram:
                dates = dates[f['dropIfgram'][:], :]
        mDates = [i.decode('utf8') for i in dates[:, 0]]
        sDates = [i.decode('utf8') for i in dates[:, 1]]
        dateList = sorted(list(set(mDates + sDates)))
        return dateList

    def nonzero_mask(self, datasetName=None, print_msg=True, dropIfgram=True):
        """Return the common mask of pixels with non-zero value in dataset of all ifgrams.
           Ignoring dropped ifgrams
        """
        self.open(print_msg=False)
        with h5py.File(self.file, 'r') as f:
            if datasetName is None:
                datasetName = [i for i in ['connectComponent', 'unwrapPhase']
                               if i in f.keys()][0]
            print('calculate the common mask of pixels with non-zero {} value'.format(datasetName))

            dset = f[datasetName]
            mask = np.ones(dset.shape[1:3], dtype=np.bool_)
            dropIfgramFlag = np.ones(dset.shape[0], dtype=np.bool_)
            if dropIfgram:
                dropIfgramFlag = self.dropIfgram
            num2read = np.sum(dropIfgramFlag)
            idx2read = np.where(dropIfgramFlag)[0]
            for i in range(num2read):  # Loop to save memory usage
                data = dset[idx2read[i], :, :]
                mask[data == 0.] = 0
                mask[np.isnan(data)] = 0
                sys.stdout.write('\rreading interferogram {}/{} ...'.format(i+1, num2read))
                sys.stdout.flush()
            print('')
        return mask

    def temporal_average(self, datasetName=ifgramDatasetNames[1], dropIfgram=True):
        self.open(print_msg=False)
        if datasetName is None:
            datasetName = ifgramDatasetNames[0]
        print('calculate the temporal average of {} in file {} ...'.format(datasetName, self.file))
        if datasetName == 'unwrapPhase':
            phase2range = -1 * float(self.metadata['WAVELENGTH']) / (4.0 * np.pi)
            tbaseIfgram = self.tbaseIfgram / 365.25

        with h5py.File(self.file, 'r') as f:
            dset = f[datasetName]
            num_ifgram, length, width = dset.shape
            dmean = np.zeros((length, width), dtype=np.float32)
            drop_ifgram_flag = np.ones(num_ifgram, dtype=np.bool_)
            if dropIfgram:
                drop_ifgram_flag = self.dropIfgram
                if np.all(drop_ifgram_flag == 0.):
                    raise Exception(('ALL interferograms are marked as dropped, '
                                     'can not calculate temporal average.'))

            #r_step = 100
            #for i in range(int(np.ceil(length / r_step))):
            #    r0 = i * r_step
            #    r1 = min(r0 + r_step, length)
            #    data = dset[drop_ifgram_flag, r0:r1, :]
            #    dmean[r0:r1, :] = np.nanmean(data, axis=0)
            #    sys.stdout.write('\rreading lines {}/{} ...'.format(r0, length))
            #    sys.stdout.flush()

            num2read = np.sum(drop_ifgram_flag)
            idx2read = np.where(drop_ifgram_flag)[0]
            for i in range(num2read):
                idx = idx2read[i]
                data = dset[idx, :, :]
                if datasetName == 'unwrapPhase':
                    data *= (phase2range * (1./tbaseIfgram[idx]))
                dmean += data
                sys.stdout.write('\rreading interferogram {}/{} ...'.format(i+1, num2read))
                sys.stdout.flush()
            dmean *= 1./np.sum(self.dropIfgram)
            print('')
        return dmean

    # Functions for Unwrap error correction
    def get_design_matrix4ifgram_triangle(self, date12_list=None):
        """Generate the design matrix of ifgram triangle for unwrap error correction using phase closure
        Parameters: date12_list : list of string in YYYYMMDD_YYYYMMDD format
        Returns:    C : 2D np.array in size of (num_tri, num_ifgram) consisting 0, 1, -1
                        for 3 SAR acquisition in t1, t2 and t3 in time order,
                        ifg1 for (t1, t2) with 1
                        ifg2 for (t1, t3) with -1
                        ifg3 for (t2, t3) with 1
        Examples:   stack_obj = ifgramStack('./INPUTS/ifgramStack.h5')
                    C = stack_obj.get_design_matrix4ifgram_triangle
        """
        if not date12_list:
            date12_list = self.get_date12_list(dropIfgram=None)

        # calculate triangle_idx
        triangle_idx = []
        for ifgram1 in date12_list:
            # ifgram1 (date1, date2)
            date1, date2 = ifgram1.split('_')

            # ifgram2 candidates (date1, date3)
            date3_list = []
            for ifgram2 in date12_list:
                if date1 == ifgram2.split('_')[0] and ifgram2 != ifgram1:
                    date3_list.append(ifgram2.split('_')[1])

            # ifgram2/3
            if len(date3_list) > 0:
                for date3 in date3_list:
                    ifgram3 = '{}_{}'.format(date2, date3)
                    if ifgram3 in date12_list:
                        ifgram1 = '{}_{}'.format(date1, date2)
                        ifgram2 = '{}_{}'.format(date1, date3)
                        ifgram3 = '{}_{}'.format(date2, date3)
                        triangle_idx.append([date12_list.index(ifgram1),
                                             date12_list.index(ifgram2),
                                             date12_list.index(ifgram3)])

        triangle_idx = np.array(triangle_idx, np.int16)
        triangle_idx = np.unique(triangle_idx, axis=0)

        # triangle_idx to C
        num_triangle = triangle_idx.shape[0]
        C = np.zeros((num_triangle, len(date12_list)), np.float32)
        for i in range(num_triangle):
            C[i, triangle_idx[i, 0]] = 1
            C[i, triangle_idx[i, 1]] = -1
            C[i, triangle_idx[i, 2]] = 1
        return C


    # Functions for Network Inversion
    def get_design_matrix4timeseries_estimation(self, refDate=None, dropIfgram=True, date12_list=None):
        """Return design matrix of the input ifgramStack for timeseries estimation
        Ignoring dropped ifgrams by default
        Parameters: refDate : str, date in YYYYMMDD format
                    dropIfgram : bool, use dropped ifgram info or not
        Returns:    A : 2D array of float32 in size of (num_ifgram, num_date-1)
                    B : 2D array of float32 in size of (num_ifgram, num_date-1)
        Examples:   stack_obj = ifgramStack('./INPUTS/ifgramStack.h5')
                    A, B = stack_obj.get_design_matrix4timeseries_estimation()
                    A, B = stack_obj.get_design_matrix4timeseries_estimation(date12_list=date12_list)
        """
        # Date info
        if date12_list:
            date12List = list(date12_list)
        else:
            date12List = self.get_date12_list(dropIfgram=dropIfgram)
        mDates = [i.split('_')[0] for i in date12List]
        sDates = [i.split('_')[1] for i in date12List]
        dateList = sorted(list(set(mDates + sDates)))
        dates = [dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in dateList]
        tbase = np.array([(i - dates[0]).days for i in dates], np.float32) / 365.25
        numIfgram = len(date12List)
        numDate = len(dateList)

        # calculate design matrix
        A = np.zeros((numIfgram, numDate), np.float32)
        B = np.zeros(A.shape, np.float32)
        for i in range(numIfgram):
            m_idx, s_idx = [dateList.index(j) for j in date12List[i].split('_')]
            A[i, m_idx] = -1
            A[i, s_idx] = 1
            B[i, m_idx:s_idx] = tbase[m_idx+1:s_idx+1] - tbase[m_idx:s_idx]

        # Remove reference date as it can not be resolved
        if not refDate:
            refDate = dateList[0]
        refIndex = dateList.index(refDate)
        A = np.hstack((A[:, 0:refIndex], A[:, (refIndex+1):]))
        B = B[:, :-1]
        return A, B

    def get_perp_baseline_timeseries(self, dropIfgram=True):
        """Get spatial perpendicular baseline in timeseries from ifgramStack, ignoring dropped ifgrams"""
        # Get tbase_diff
        date12List = self.get_date12_list(dropIfgram=dropIfgram)
        mDates = [i.split('_')[0] for i in date12List]
        sDates = [i.split('_')[1] for i in date12List]
        dateList = sorted(list(set(mDates + sDates)))
        dates = [dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in dateList]
        tbase = np.array([(i - dates[0]).days for i in dates], np.float32) / 365.25
        tbase_diff = np.diff(tbase).flatten()

        B = self.get_design_matrix4timeseries_estimation(dropIfgram=dropIfgram)[1]
        B_inv = np.linalg.pinv(B)
        with h5py.File(self.file, 'r') as f:
            pbaseIfgram = f['bperp'][:]
            if dropIfgram:
                pbaseIfgram = pbaseIfgram[f['dropIfgram'][:]]
        pbaseRate = np.dot(B_inv, pbaseIfgram)
        pbaseTimeseries = np.concatenate((np.array([0.], dtype=np.float32),
                                          np.cumsum([pbaseRate * tbase_diff])))
        return pbaseTimeseries

    def update_drop_ifgram(self, date12List_to_drop):
        """Update dropIfgram dataset based on input date12List_to_drop"""
        if date12List_to_drop is None:
            return
        date12ListAll = self.get_date12_list(dropIfgram=False)
        with h5py.File(self.file, 'r+') as f:
            print('open file {} with r+ mode'.format(self.file))
            print('update HDF5 dataset "/dropIfgram".')
            f['dropIfgram'][:] = np.array([i not in date12List_to_drop for i in date12ListAll],
                                          dtype=np.bool_)


########################################################################################
class singleDataset:
    def __init__(self, file=None):
        self.file = file

    def close(self, print_msg=True):
        try:
            self.f.close()
            if print_msg:
                print('close file: {}'.format(os.path.basename(self.file)))
        except:
            pass

    def read(self, box=None):
        self.f = h5py.File(self.file, 'r')
        k = list(self.f.keys())[0]
        data = self.f[k][:]
        if box is not None:
            data = data[box[1]:box[3], box[0]:box[2]]
        return data


########################################################################################
FILE_STRUCTURE_HDFEOS = """
/                             Root level group
Attributes                    metadata in dict.
/HDFEOS/GRIDS/timeseries      timeseries group
    /observation
        /displacement         3D array of float32 in size of (n, l, w) in meter
        /date                 1D array of string  in size of (n,     ) in YYYYMMDD format.
        /bperp                1D array of float32 in size of (n,     ) in meter
    /quality
        /temporalCoherence    2D array of float32 in size of (   l, w).
        /mask                 2D array of bool_   in size of (   l, w).
    /geometry
        /height               2D array of float32 in size of (   l, w) in meter.
        /incidenceAngle       2D array of float32 in size of (   l, w) in degree.
        /slantRangeDistance   2D array of float32 in size of (   l, w) in meter.
        /azimuthAngle         2D array of float32 in size of (   l, w) in degree. (optional)
        /shadowMask           2D array of bool    in size of (   l, w).           (optional)
        /waterMask            2D array of bool    in size of (   l, w).           (optional)
        /bperp                3D array of float32 in size of (n, l, w) in meter.  (optional)
"""

class HDFEOS:
    """
    Time-series object in HDF-EOS5 format for Univ of Miami's InSAR Time-series Web Viewer
        Link: http://insarmaps.miami.edu
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.
    """

    def __init__(self, file=None):
        self.file = file
        self.name = 'HDFEOS'
        self.file_structure = FILE_STRUCTURE_HDFEOS
        self.datasetGroupNameDict = {'displacement'       : 'observation',
                                     'raw'                : 'observation',
                                     'troposphericDelay'  : 'observation',
                                     'topographicResidual': 'observation',
                                     'ramp'               : 'observation',
                                     'date'               : 'observation', 
                                     'temporalCoherence'  : 'quality', 
                                     'mask'               : 'quality', 
                                     'coherence'          : 'quality', 
                                     'variance'           : 'quality', 
                                     'uncertainty'        : 'quality', 
                                     'height'             : 'geometry',
                                     'incidenceAngle'     : 'geometry',
                                     'slantRangeDistance' : 'geometry',
                                     'azimuthAngle'       : 'geometry',
                                     'shadowMask'         : 'geometry',
                                     'waterMask'          : 'geometry',
                                     'bperp'              : 'geometry'
                                    }

    def close(self, print_msg=True):
        try:
            self.f.close()
            if print_msg:
                print('close timeseries file: {}'.format(os.path.basename(self.file)))
        except:
            pass

    def open(self, print_msg=True):
        if print_msg:
            print('open {} file: {}'.format(self.name, os.path.basename(self.file)))
        self.get_metadata()
        self.length = int(self.metadata['LENGTH'])
        self.width = int(self.metadata['WIDTH'])

        self.datasetList = []
        with h5py.File(self.file, 'r') as f:
            g = f['HDFEOS/GRIDS/timeseries/observation']
            self.dateList = [i.decode('utf8') for i in g['date'][:]]
            self.pbase = g['bperp'][:]
            self.datasetList += ['displacement-{}'.format(i) for i in self.dateList]
            self.numDate = len(self.dateList)

            g = f['HDFEOS/GRIDS/timeseries/quality']
            for key in g.keys():
                self.datasetList.append(key)

            g = f['HDFEOS/GRIDS/timeseries/geometry']
            for key in g.keys():
                obj = g[key]
                if isinstance(obj, h5py.Dataset):
                    if len(obj.shape) == 2:
                        self.datasetList.append(key)
                    # elif len(obj.shape) == 3:
                    #    self.datesetList += ['{}-{}'.format(key, i) for i in self.dateList]

    def get_metadata(self):
        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)
            dates = f['HDFEOS/GRIDS/timeseries/observation/date'][:]
        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except:
                self.metadata[key] = value
        self.metadata['FILE_TYPE'] = self.name

        # ref_date/index
        dateList = [i.decode('utf8') for i in dates]
        if 'REF_DATE' not in self.metadata.keys():
            self.metadata['REF_DATE'] = dateList[0]
        self.refIndex = dateList.index(self.metadata['REF_DATE'])
        return self.metadata

    def read(self, datasetName=None, box=None, print_msg=True):
        """Read dataset from HDF-EOS5 file
        Parameters: self : HDFEOS object
                    datasetName : (list of) str
                    box : tuple of 4 int, for (x0, y0, x1, y1)
                    print_msg : bool
        Returns:    data: 2D or 3D array
        Example:    obj = HDFEOS('S1_IW1_128_0593_0597_20141213_20171221.he5')
                    obj.read(datasetName='displacement')
                    obj.read(datasetName='displacement-20150915')
                    obj.read(datasetName=['displacement-20150915',
                                          'displacement-20150921'])
                    obj.read(datasetName='incidenceAngle')
        """
        self.open(print_msg=False)
        if box is None:
            box = [0, 0, self.width, self.length]

        if datasetName is None:
            datasetName = [timeseriesDatasetNames[-1]]
        elif isinstance(datasetName, str):
            datasetName = [datasetName]

        with h5py.File(self.file, 'r') as f:
            familyName = datasetName[0].split('-')[0]
            groupName = self.datasetGroupNameDict[familyName]
            ds = f['HDFEOS/GRIDS/timeseries/{}/{}'.format(groupName, familyName)]
            if print_msg:
                print('reading {} data from file: {} ...'.format(familyName, self.file))

            if len(ds.shape) == 1:
                data = ds[:]
            elif len(ds.shape) == 2:
                data = ds[box[1]:box[3], box[0]:box[2]]
            else:
                # get dateFlag - mark in time/1st dimension
                dateFlag = np.zeros((ds.shape[0]), dtype=np.bool_)
                datasetName = [i.replace(familyName, '').replace('-', '') for i in datasetName]
                if any(not i for i in datasetName):
                    dateFlag[:] = True
                else:
                    for e in datasetName:
                        dateFlag[self.dateList.index(e)] = True
                data = ds[dateFlag, box[1]:box[3], box[0]:box[2]]
                data = np.squeeze(data)
        return data


########################################################################################
FILE_STRUCTURE_GIANT_TIMESERIES = """
/                Root level
/cmask           2D array of float32 in size of (   l, w).
/dates           1D array of float32 in size of (n,     ) in ordinal date format
/recons          3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - filtered
/rawts           3D array of float32 in size of (n, l, w) in mm, reconstructed timeseries - un-filtered
"""

class giantTimeseries:
    """
    Time-series object for displacement of a set of SAR images from the same platform and track.
    It contains a "timeseries" group and three datasets: date, bperp and timeseries.
    """

    def __init__(self, file=None):
        self.file = file
        self.name = 'giantTimeseries'
        self.file_structure = FILE_STRUCTURE_GIANT_TIMESERIES

    def open(self, print_msg=True):
        if print_msg:
            print('open {} file: {}'.format(self.name, os.path.basename(self.file)))
        self.get_size()
        self.get_metadata()
        self.numPixel = self.length * self.width

        self.times = np.array([dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in self.dateList])
        self.tbase = np.array([i.days for i in self.times - self.times[self.refIndex]],
                              dtype=np.float32)
        # list of float for year, 2014.95
        self.yearList = [i.year + (i.timetuple().tm_yday-1)/365.25 for i in self.times]
        self.datasetList = ['{}-{}'.format(self.name, i) for i in self.dateList]

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

        # ref_date/index
        dateList = self.get_date_list()
        if 'REF_DATE' not in self.metadata.keys():
            self.metadata['REF_DATE'] = dateList[0]
        self.refIndex = dateList.index(self.metadata['REF_DATE'])
        return self.metadata

    def read(self, datasetName=None, box=None, print_msg=True):
        """Read dataset from timeseries file
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
        """
        if print_msg:
            print('reading {} data from file: {} ...'.format(self.name, self.file))
        self.open(print_msg=False)

        # convert input datasetName into list of dates
        if not datasetName or datasetName == 'timeseries':
            datasetName = []
        elif isinstance(datasetName, str):
            datasetName = [datasetName]
        datasetName = [i.replace('timeseries', '').replace('-', '') for i in datasetName]

        with h5py.File(self.file, 'r') as f:
            ds = f[self.name]
            if isinstance(ds, h5py.Group):  # support for old pysar files
                ds = ds[self.name]

            # Get dateFlag - mark in time/1st dimension
            dateFlag = np.zeros((self.numDate), dtype=np.bool_)
            if not datasetName:
                dateFlag[:] = True
            else:
                for e in datasetName:
                    dateFlag[self.dateList.index(e)] = True

            # Get Index in space/2_3 dimension
            if box is None:
                box = [0, 0, self.width, self.length]

            data = ds[dateFlag, box[1]:box[3], box[0]:box[2]]
            data = np.squeeze(data)
        return data

