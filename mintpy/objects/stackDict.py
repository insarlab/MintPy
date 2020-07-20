############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2017               #
############################################################
# class used for data loading from InSAR stack to MintPy timeseries
# Recommend import:
#     from mintpy.objects.stackDict import (geometryDict,
#                                           ifgramStackDict,
#                                           ifgramDict)


import os
import time
import warnings
import h5py
import numpy as np
try:
    from skimage.transform import resize
except ImportError:
    raise ImportError('Could not import skimage!')

from mintpy.objects import (
    dataTypeDict,
    geometryDatasetNames,
    ifgramDatasetNames,
)
from mintpy.utils import readfile, ptime, utils0 as ut


########################################################################################
class ifgramStackDict:
    '''
    IfgramStack object for a set of InSAR pairs from the same platform and track.

    Example:
        from mintpy.objects.insarobj import ifgramStackDict
        pairsDict = {('20160524','20160530'):ifgramObj1,
                     ('20160524','20160605'):ifgramObj2,
                     ('20160524','20160611'):ifgramObj3,
                     ('20160530','20160605'):ifgramObj4,
                     ...
                     }
        stackObj = ifgramStackDict(pairsDict=pairsDict)
        stackObj.write2hdf5(outputFile='ifgramStack.h5', box=(200,500,300,600))
    '''

    def __init__(self, name='ifgramStack', pairsDict=None, dsName0=ifgramDatasetNames[0]):
        self.name = name
        self.pairsDict = pairsDict
        self.dsName0 = dsName0        #reference dataset name, unwrapPhase OR azimuthOffset

    def get_size(self, box=None):
        self.numIfgram = len(self.pairsDict)
        ifgramObj = [v for v in self.pairsDict.values()][0]
        self.length, ifgramObj.width = ifgramObj.get_size(family=self.dsName0)
        if box:
            self.length = box[3] - box[1]
            self.width = box[2] - box[0]
        else:
            self.length = ifgramObj.length
            self.width = ifgramObj.width
        return self.numIfgram, self.length, self.width

    def get_date12_list(self):
        pairs = [pair for pair in self.pairsDict.keys()]
        self.date12List = ['{}_{}'.format(i[0], i[1]) for i in pairs]
        return self.date12List

    def get_metadata(self):
        ifgramObj = [v for v in self.pairsDict.values()][0]
        self.metadata = ifgramObj.get_metadata(family=self.dsName0)
        if 'UNIT' in self.metadata.keys():
            self.metadata.pop('UNIT')
        return self.metadata

    def get_dataset_data_type(self, dsName):
        ifgramObj = [v for v in self.pairsDict.values()][0]
        dsFile = ifgramObj.datasetDict[dsName]
        metadata = readfile.read_attribute(dsFile)
        dataType = dataTypeDict[metadata.get('DATA_TYPE', 'float32').lower()]
        return dataType

    def write2hdf5(self, outputFile='ifgramStack.h5', access_mode='w', box=None, compression=None, extra_metadata=None):
        '''Save/write an ifgramStackDict object into an HDF5 file with the structure defined in:

        https://mintpy.readthedocs.io/en/latest/api/data_structure/#ifgramstack

        Parameters: outputFile : str, Name of the HDF5 file for the InSAR stack
                    access_mode : str, access mode of output File, e.g. w, r+
                    box : tuple, subset range in (x0, y0, x1, y1)
                    extra_metadata : dict, extra metadata to be added into output file
        Returns:    outputFile
        '''

        self.outputFile = outputFile
        f = h5py.File(self.outputFile, access_mode)
        print('create HDF5 file {} with {} mode'.format(self.outputFile, access_mode))

        self.pairs = sorted([pair for pair in self.pairsDict.keys()])
        self.dsNames = list(self.pairsDict[self.pairs[0]].datasetDict.keys())
        self.dsNames = [i for i in ifgramDatasetNames if i in self.dsNames]
        maxDigit = max([len(i) for i in self.dsNames])
        self.get_size(box)

        self.bperp = np.zeros(self.numIfgram)
        ###############################
        # 3D datasets containing unwrapPhase, coherence, connectComponent, wrapPhase, etc.
        for dsName in self.dsNames:
            dsShape = (self.numIfgram, self.length, self.width)
            dsDataType = np.float32
            dsCompression = compression
            if dsName in ['connectComponent']:
                dsDataType = np.int16
                dsCompression = 'lzf'

            print(('create dataset /{d:<{w}} of {t:<25} in size of {s}'
                   ' with compression = {c}').format(d=dsName,
                                                     w=maxDigit,
                                                     t=str(dsDataType),
                                                     s=dsShape,
                                                     c=dsCompression))
            ds = f.create_dataset(dsName,
                                  shape=dsShape,
                                  maxshape=(None, dsShape[1], dsShape[2]),
                                  dtype=dsDataType,
                                  chunks=True,
                                  compression=dsCompression)

            prog_bar = ptime.progressBar(maxValue=self.numIfgram)
            for i in range(self.numIfgram):
                ifgramObj = self.pairsDict[self.pairs[i]]
                data = ifgramObj.read(dsName, box=box)[0]
                ds[i, :, :] = data
                self.bperp[i] = ifgramObj.get_perp_baseline(family=self.dsName0)
                prog_bar.update(i+1, suffix='{}_{}'.format(self.pairs[i][0],
                                                           self.pairs[i][1]))
            prog_bar.close()
            ds.attrs['MODIFICATION_TIME'] = str(time.time())

        ###############################
        # 2D dataset containing reference and secondary dates of all pairs
        dsName = 'date'
        dsDataType = np.string_
        dsShape = (self.numIfgram, 2)
        print('create dataset /{d:<{w}} of {t:<25} in size of {s}'.format(d=dsName,
                                                                          w=maxDigit,
                                                                          t=str(dsDataType),
                                                                          s=dsShape))
        data = np.array(self.pairs, dtype=dsDataType)
        f.create_dataset(dsName, data=data)

        ###############################
        # 1D dataset containing perpendicular baseline of all pairs
        dsName = 'bperp'
        dsDataType = np.float32
        dsShape = (self.numIfgram,)
        print('create dataset /{d:<{w}} of {t:<25} in size of {s}'.format(d=dsName,
                                                                          w=maxDigit,
                                                                          t=str(dsDataType),
                                                                          s=dsShape))
        data = np.array(self.bperp, dtype=dsDataType)
        f.create_dataset(dsName, data=data)

        ###############################
        # 1D dataset containing bool value of dropping the interferograms or not
        dsName = 'dropIfgram'
        dsDataType = np.bool_
        dsShape = (self.numIfgram,)
        print('create dataset /{d:<{w}} of {t:<25} in size of {s}'.format(d=dsName,
                                                                          w=maxDigit,
                                                                          t=str(dsDataType),
                                                                          s=dsShape))
        data = np.ones(dsShape, dtype=dsDataType)
        f.create_dataset(dsName, data=data)

        ###############################
        # Attributes
        self.get_metadata()
        if extra_metadata:
            self.metadata.update(extra_metadata)
            print('add extra metadata: {}'.format(extra_metadata))
        self.metadata = ut.subset_attribute(self.metadata, box)
        self.metadata['FILE_TYPE'] = self.name
        for key, value in self.metadata.items():
            f.attrs[key] = value

        f.close()
        print('Finished writing to {}'.format(self.outputFile))
        return self.outputFile


########################################################################################
class ifgramDict:
    """
    Ifgram object for a single InSAR pair of interferogram. It includes dataset name (family) of:
        'unwrapPhase','coherence','connectComponent','wrapPhase','ionoPhase','rangeOffset','azimuthOffset', etc.

    Example:
        from mintpy.objects.insarobj import ifgramDict
        datasetDict = {'unwrapPhase'     :'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.unw',
                       'coherence'       :'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.cor',
                       'connectComponent':'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.unw.conncomp',
                       'wrapPhase'       :'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.int',
                       'ionoPhase'       :'$PROJECT_DIR/merged/ionosphere/20151220_20160206/iono.bil.unwCor.filt',
                       ...
                      }
        ifgramObj = ifgramDict(datasetDict=datasetDict)
        data, atr = ifgramObj.read('unwrapPhase')
    """

    def __init__(self, name='ifgram', datasetDict={}, metadata=None):
        self.name = name
        self.datasetDict = datasetDict

        self.platform = None
        self.track = None
        self.processor = None
        # platform, track and processor can get values from metadat if they exist
        if metadata is not None:
            for key, value in metadata.items():
                setattr(self, key, value)

    def read(self, family, box=None):
        self.file = self.datasetDict[family]
        data, metadata = readfile.read(self.file, box=box)
        return data, metadata

    def get_size(self, family=ifgramDatasetNames[0]):
        self.file = self.datasetDict[family]
        metadata = readfile.read_attribute(self.file)
        self.length = int(metadata['LENGTH'])
        self.width = int(metadata['WIDTH'])
        return self.length, self.width

    def get_perp_baseline(self, family=ifgramDatasetNames[0]):
        self.file = self.datasetDict[family]
        metadata = readfile.read_attribute(self.file)
        self.bperp_top = float(metadata['P_BASELINE_TOP_HDR'])
        self.bperp_bottom = float(metadata['P_BASELINE_BOTTOM_HDR'])
        self.bperp = (self.bperp_top + self.bperp_bottom) / 2.0
        return self.bperp

    def get_metadata(self, family=ifgramDatasetNames[0]):
        self.file = self.datasetDict[family]
        self.metadata = readfile.read_attribute(self.file)
        self.length = int(self.metadata['LENGTH'])
        self.width = int(self.metadata['WIDTH'])

        # if self.processor is None:
        #    ext = self.file.split('.')[-1]
        #    if 'PROCESSOR' in self.metadata.keys():
        #        self.processor = self.metadata['PROCESSOR']
        #    elif os.path.exists(self.file+'.xml'):
        #        self.processor = 'isce'
        #    elif os.path.exists(self.file+'.rsc'):
        #        self.processor = 'roipac'
        #    elif os.path.exists(self.file+'.par'):
        #        self.processor = 'gamma'
        #    elif ext == 'grd':
        #        self.processor = 'gmtsar'
        #    #what for DORIS/SNAP
        #    else:
        #        self.processor = 'isce'
        #self.metadata['PROCESSOR'] = self.processor

        if self.track:
            self.metadata['TRACK'] = self.track

        if self.platform:
            self.metadata['PLATFORM'] = self.platform

        return self.metadata


########################################################################################
class geometryDict:
    '''
    Geometry object for Lat, Lon, Heigt, Incidence, Heading, Bperp, ... from the same platform and track.

    Example:
        from mintpy.utils import readfile
        from mintpy.utils.insarobj import geometryDict
        datasetDict = {'height'        :'$PROJECT_DIR/merged/geom_reference/hgt.rdr',
                       'latitude'      :'$PROJECT_DIR/merged/geom_reference/lat.rdr',
                       'longitude'     :'$PROJECT_DIR/merged/geom_reference/lon.rdr',
                       'incidenceAngle':'$PROJECT_DIR/merged/geom_reference/los.rdr',
                       'heandingAngle' :'$PROJECT_DIR/merged/geom_reference/los.rdr',
                       'shadowMask'    :'$PROJECT_DIR/merged/geom_reference/shadowMask.rdr',
                       'bperp'         :bperpDict
                       ...
                      }
        bperpDict = {'20160406':'$PROJECT_DIR/merged/baselines/20160406/bperp',
                     '20160418':'$PROJECT_DIR/merged/baselines/20160418/bperp',
                     ...
                    }
        metadata = readfile.read_attribute('$PROJECT_DIR/merged/interferograms/20160629_20160723/filt_fine.unw')
        geomObj = geometryDict(processor='isce', datasetDict=datasetDict, extraMetadata=metadata)
        geomObj.write2hdf5(outputFile='geometryRadar.h5', access_mode='w', box=(200,500,300,600))
    '''

    def __init__(self, name='geometry', processor=None, datasetDict={}, extraMetadata=None):
        self.name = name
        self.processor = processor
        self.datasetDict = datasetDict
        self.extraMetadata = extraMetadata

        # get extra metadata from geometry file if possible
        self.dsNames = list(self.datasetDict.keys())
        if not self.extraMetadata:
            dsFile = self.datasetDict[self.dsNames[0]]
            metadata = readfile.read_attribute(dsFile)
            if all(i in metadata.keys() for i in ['STARTING_RANGE', 'RANGE_PIXEL_SIZE']):
                self.extraMetadata = metadata

    def read(self, family, box=None):
        self.file = self.datasetDict[family]
        data, metadata = readfile.read(self.file,
                                       datasetName=family,
                                       box=box)
        return data, metadata

    def get_slant_range_distance(self, box=None):
        """Generate 2D slant range distance if missing from input template file"""
        if not self.extraMetadata:
            return None

        if 'Y_FIRST' in self.extraMetadata.keys():
            # for dataset in geo-coordinates, use contant value from SLANT_RANGE_DISTANCE.
            key = 'SLANT_RANGE_DISTANCE'
            print('geocoded input, use contant value from metadata {}'.format(key))
            if key in self.extraMetadata.keys():
                length = int(self.extraMetadata['LENGTH'])
                width = int(self.extraMetadata['WIDTH'])
                range_dist = float(self.extraMetadata[key])
                data = np.ones((length, width), dtype=np.float32) * range_dist

            else:
                return None

        else:
            # for dataset in radar-coordinates, calculate 2D pixel-wise value from geometry
            data = ut.range_distance(self.extraMetadata,
                                     dimension=2,
                                     print_msg=False)

        # subset
        if box is not None:
            data = data[box[1]:box[3], box[0]:box[2]]
        return data

    def get_incidence_angle(self, box=None):
        """Generate 2D slant range distance if missing from input template file"""
        if not self.extraMetadata:
            return None

        if 'Y_FIRST' in self.extraMetadata.keys():
            # for dataset in geo-coordinates, use contant value from INCIDENCE_ANGLE.
            key = 'INCIDENCE_ANGLE'
            print('geocoded input, use contant value from metadata {}'.format(key))
            if key in self.extraMetadata.keys():
                length = int(self.extraMetadata['LENGTH'])
                width = int(self.extraMetadata['WIDTH'])
                inc_angle = float(self.extraMetadata[key])
                data = np.ones((length, width), dtype=np.float32) * inc_angle

            else:
                return None

        else:
            # read DEM if available for more previse calculation
            if 'height' in self.dsNames:
                dem = readfile.read(self.datasetDict['height'], datasetName='height')[0]
            else:
                dem = None

            # for dataset in radar-coordinates, calculate 2D pixel-wise value from geometry
            data = ut.incidence_angle(self.extraMetadata,
                                      dem=dem,
                                      dimension=2,
                                      print_msg=False)

        # subset
        if box is not None:
            data = data[box[1]:box[3], box[0]:box[2]]
        return data

    def get_size(self, family=None, box=None):
        if not family:
            family = [i for i in self.datasetDict.keys() if i != 'bperp'][0]
        self.file = self.datasetDict[family]
        metadata = readfile.read_attribute(self.file)
        if box:
            length = box[3] - box[1]
            width = box[2] - box[0]
        else:
            length = int(metadata['LENGTH'])
            width = int(metadata['WIDTH'])
        return length, width

    def get_dataset_list(self):
        self.datasetList = list(self.datasetDict.keys())
        return self.datasetList

    def get_metadata(self, family=None):
        if not family:
            family = [i for i in self.datasetDict.keys() if i != 'bperp'][0]
        self.file = self.datasetDict[family]
        self.metadata = readfile.read_attribute(self.file)
        self.length = int(self.metadata['LENGTH'])
        self.width = int(self.metadata['WIDTH'])
        if 'UNIT' in self.metadata.keys():
            self.metadata.pop('UNIT')

        # if self.processor is None:
        #    ext = self.file.split('.')[-1]
        #    if 'PROCESSOR' in self.metadata.keys():
        #        self.processor = self.metadata['PROCESSOR']
        #    elif os.path.exists(self.file+'.xml'):
        #        self.processor = 'isce'
        #    elif os.path.exists(self.file+'.rsc'):
        #        self.processor = 'roipac'
        #    elif os.path.exists(self.file+'.par'):
        #        self.processor = 'gamma'
        #    elif ext == 'grd':
        #        self.processor = 'gmtsar'
        #    #what for DORIS/SNAP
        #    else:
        #        self.processor = 'isce'
        #self.metadata['PROCESSOR'] = self.processor
        return self.metadata

    def write2hdf5(self, outputFile='geometryRadar.h5', access_mode='w', box=None, compression='lzf', extra_metadata=None):
        ''' Save/write to HDF5 file with structure defined in:

        https://mintpy.readthedocs.io/en/latest/api/data_structure/#geometry
        '''
        if len(self.datasetDict) == 0:
            print('No dataset file path in the object, skip HDF5 file writing.')
            return None

        self.outputFile = outputFile
        f = h5py.File(self.outputFile, access_mode)
        print('create HDF5 file {} with {} mode'.format(self.outputFile, access_mode))

        #groupName = self.name
        #group = f.create_group(groupName)
        #print('create group   /{}'.format(groupName))

        maxDigit = max([len(i) for i in geometryDatasetNames])
        length, width = self.get_size(box=box)
        self.length, self.width = self.get_size()

        ###############################
        for dsName in self.dsNames:
            # 3D datasets containing bperp
            if dsName == 'bperp':
                self.dateList = list(self.datasetDict[dsName].keys())
                dsDataType = np.float32
                self.numDate = len(self.dateList)
                dsShape = (self.numDate, length, width)
                ds = f.create_dataset(dsName,
                                      shape=dsShape,
                                      maxshape=(None, dsShape[1], dsShape[2]),
                                      dtype=dsDataType,
                                      chunks=True,
                                      compression=compression)
                print(('create dataset /{d:<{w}} of {t:<25} in size of {s}'
                       ' with compression = {c}').format(d=dsName,
                                                         w=maxDigit,
                                                         t=str(dsDataType),
                                                         s=dsShape,
                                                         c=str(compression)))

                print('read coarse grid baseline files and linear interpolate into full resolution ...')
                prog_bar = ptime.progressBar(maxValue=self.numDate)
                for i in range(self.numDate):
                    fname = self.datasetDict[dsName][self.dateList[i]]
                    data = read_isce_bperp_file(fname=fname,
                                                out_shape=(self.length, self.width),
                                                box=box)
                    ds[i, :, :] = data
                    prog_bar.update(i+1, suffix=self.dateList[i])
                prog_bar.close()

                # Write 1D dataset date
                dsName = 'date'
                dsShape = (self.numDate,)
                dsDataType = np.string_
                print(('create dataset /{d:<{w}} of {t:<25}'
                       ' in size of {s}').format(d=dsName,
                                                 w=maxDigit,
                                                 t=str(dsDataType),
                                                 s=dsShape))
                data = np.array(self.dateList, dtype=dsDataType)
                ds = f.create_dataset(dsName, data=data)

            # 2D datasets containing height, latitude, incidenceAngle, shadowMask, etc.
            else:
                dsDataType = np.float32
                if dsName.lower().endswith('mask'):
                    dsDataType = np.bool_
                dsShape = (length, width)
                print(('create dataset /{d:<{w}} of {t:<25} in size of {s}'
                       ' with compression = {c}').format(d=dsName,
                                                         w=maxDigit,
                                                         t=str(dsDataType),
                                                         s=dsShape,
                                                         c=str(compression)))
                data = np.array(self.read(family=dsName, box=box)[0], dtype=dsDataType)
                ds = f.create_dataset(dsName,
                                      data=data,
                                      chunks=True,
                                      compression=compression)

        ###############################
        # Generate Dataset if not existed in binary file: incidenceAngle, slantRangeDistance
        for dsName in [i for i in ['incidenceAngle', 'slantRangeDistance']
                       if i not in self.dsNames]:
            # Calculate data
            data = None
            if dsName == 'incidenceAngle':
                data = self.get_incidence_angle(box=box)
            elif dsName == 'slantRangeDistance':
                data = self.get_slant_range_distance(box=box)

            # Write dataset
            if data is not None:
                dsShape = data.shape
                dsDataType = np.float32
                print(('create dataset /{d:<{w}} of {t:<25} in size of {s}'
                       ' with compression = {c}').format(d=dsName,
                                                         w=maxDigit,
                                                         t=str(dsDataType),
                                                         s=dsShape,
                                                         c=str(compression)))
                ds = f.create_dataset(dsName,
                                      data=data,
                                      dtype=dsDataType,
                                      chunks=True,
                                      compression=compression)

        ###############################
        # Attributes
        self.get_metadata()
        if extra_metadata:
            self.metadata.update(extra_metadata)
            print('add extra metadata: {}'.format(extra_metadata))
        self.metadata = ut.subset_attribute(self.metadata, box)
        self.metadata['FILE_TYPE'] = self.name
        for key, value in self.metadata.items():
            f.attrs[key] = value

        f.close()
        print('Finished writing to {}'.format(self.outputFile))
        return self.outputFile


########################################################################################
def read_isce_bperp_file(fname, out_shape, box=None):
    '''Read ISCE coarse grid perpendicular baseline file, and project it to full size
    Parameters: self : geometry object,
                fname : str, bperp file name
                outShape : tuple of 2int, shape of file in full resolution
                box : tuple of 4 int, subset range in (x0, y0, x1, y1) with respect to full resolution
    Returns:    data : 2D array of float32
    Example:    fname = '$PROJECT_DIR/merged/baselines/20160418/bperp'
                data = self.read_sice_bperp_file(fname, (3600,2200), box=(200,400,1000,1000))
    '''
    # read original data
    data_c = readfile.read(fname)[0]

    # resize to full resolution
    data_min, data_max = np.nanmin(data_c), np.nanmax(data_c)
    if data_max != data_min:
        data_c = (data_c - data_min) / (data_max - data_min)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        data = resize(data_c, out_shape, order=1, mode='edge', preserve_range=True)
    if data_max != data_min:
        data = data * (data_max - data_min) + data_min

    # for debug
    debug_mode=False
    if debug_mode:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8,6));
        im = ax1.imshow(readfile.read(fname)[0]); fig.colorbar(im, ax=ax1)
        im = ax2.imshow(data);  fig.colorbar(im, ax=ax2)
        plt.show()

    if box is not None:
        data = data[box[1]:box[3], box[0]:box[2]]
    return data


########################################################################################
class platformTrack:

    def __init__(self, name='platformTrack'):  # , pairDict = None):
        self.pairs = None

    def getPairs(self, pairDict, platTrack):
        pairs = pairDict.keys()
        self.pairs = {}
        for pair in pairs:
            if pairDict[pair].platform_track == platTrack:
                self.pairs[pair] = pairDict[pair]

    def getSize_geometry(self, dsName):
        pairs = self.pairs.keys()
        pairs2 = []
        width = []
        length = []
        files = []
        for pair in pairs:
            self.pairs[pair].get_metadata(dsName)
            if self.pairs[pair].length != 0 and self.pairs[pair].file not in files:
                files.append(self.pairs[pair].file)
                pairs2.append(pair)
                width.append(self.pairs[pair].width)
                length.append(self.pairs[pair].length)

        length = median(length)
        width = median(width)
        return pairs2, length, width

    def getSize(self):
        pairs = self.pairs.keys()
        self.numPairs = len(pairs)
        width = []
        length = []
        for pair in pairs:
            length.append(self.pairs[pair].length)
            width.append(self.pairs[pair].width)
        self.length = median(length)
        self.width = median(width)

    def getDatasetNames(self):
        # extract the name of the datasets which are actually the keys of
        # observations, quality and geometry dictionaries.

        pairs = [pair for pair in self.pairs.keys()]
        # Assuming all pairs of a given platform-track have the same observations
        # let's extract the keys of the observations of the first pair.

        if self.pairs[pairs[0]].observationsDict is not None:
            self.dsetObservationNames = [k for k in self.pairs[pairs[0]].observationsDict.keys()]
        else:
            self.dsetObservationNames = []

        # Assuming all pairs of a given platform-track have the same quality files
        # let's extract the keys of the quality dictionary of the first pair.
        if self.pairs[pairs[0]].qualityDict is not None:
            self.dsetQualityNames = [k for k in self.pairs[pairs[0]].qualityDict.keys()]
        else:
            self.dsetQualityNames = []

        ##################
        # Despite the observation and quality files, the geometry may not exist
        # for all pairs. Therfore we need to look at all pairs and get possible
        # dataset names.
        self.dsetGeometryNames = []
        for pair in pairs:
            if self.pairs[pair].geometryDict is not None:
                keys = [k for k in self.pairs[pair].geometryDict.keys()]
                self.dsetGeometryNames = list(set(self.dsetGeometryNames) | set(keys))
