############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun, 2017              #
############################################################


import os, sys, glob
import h5py
import numpy as np
from pysar.utils import readfile, datetime as ptime

dataType = np.float32

timeseriesDatasetNames = ['raw',\
                          'troposphericDelay',\
                          'topographicResidual',\
                          'ramp',\
                          'displacement']

geometryDatasetNames = ['height',\
                        'latitude',\
                        'longitude',\
                        'rangeCoord',\
                        'azimuthCoord',\
                        'incidenceAngle',\
                        'headingAngle',\
                        'slantRangeDistance',\
                        'shadowMask',\
                        'waterMask',\
                        'commonMask',\
                        'bperp']

ifgramDatasetNames = ['unwrapPhase',
                      'coherence',\
                      'connectComponent',\
                      'wrapPhase',\
                      'iono',\
                      'rangeOffset',\
                      'azimuthOffset']

########################################################################################
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


########################################################################################
class geometry:
    '''
    Geometry object for Lat, Lon, Heigt, Incidence, Heading, Bperp, ... from the same platform and track.
    Attributes are saved in the root level.
    It contains a "geometry" group and some datasets: date12, bperp, unwrapIfgram, coherence, ...

    /                    Root level
    Attributes           Dictionary for metadata
    geometry             group name
        /height          2D array of float32 in size of (l, w) in meter.
        /incidenceAngle  2D array of float32 in size of (l, w) in degree.
        /latitude        2D array of float32 in size of (l, w) in degree.
        /longitude       2D array of float32 in size of (l, w) in degree.
        /rangeCoord      2D array of float32 in size of (l, w).
        /azimuthCoord    2D array of float32 in size of (l, w).
        /headingAngle    2D array of float32 in size of (l, w) in degree. (optional) 
        /shadowMask      2D array of bool in size of (l, w). (optional)
        /waterMask       2D array of bool in size of (l, w). (optional)
        /bperp           3D array of float32 in size of (n, lc, wc) in meter (optional)
        ...
    '''
    def __init__(self, name='geometry', datasetDict={}, metadata=None):
        self.datasetDict = datasetDict

        if metadata is not None:
            for key , value in metadata.items():
                setattr(self, key, value)

    def get_size(self):
        return None
    def get_metadata(self):
        return None
    def read(self):
        return None
    def save2h5(self):
        return None

########################################################################################
class ifgramStack:
    '''
    IfgramStack object for a set of InSAR pairs from the same platform and track.

    Example:
        from pysar.objects.ifgramStack import ifgramStack
        pairsDict = {('20160524','20160530'):ifgramObj1,
                     ('20160524','20160605'):ifgramObj2,
                     ('20160524','20160611'):ifgramObj3,
                     ('20160530','20160605'):ifgramObj4,
                     ...
                     }
        stackObj = ifgramStack(pairsDict=pairsDict)
        stackObj.save2h5(outputFile='ifgramStack.h5', box=(200,500,300,600))
    '''

    def __init__(self, name='ifgramStack', pairsDict=None):
        self.pairsDict = pairsDict

    def get_size(self, box=None):
        self.numIfgram = len(self.pairsDict)
        ifgramObj = [v for v in self.pairsDict.values()][0]
        self.length, ifgramObj.width = ifgramObj.get_size()
        if box:
            self.length = box[3] - box[1]
            self.width = box[2] - box[0]
        else:
            self.length = ifgramObj.length
            self.width = ifgramObj.width
        return self.numIfgram, self.length, self.width

    def get_metadata(self):
        ifgramObj = [v for v in self.pairsDict.values()][0]
        self.metadata = ifgramObj.get_metadata()
        return self.metadata

    def save2h5(self, outputFile='ifgramStack.h5', access_mode='w', box=None):
        '''Save/write an ifgramStack object into an HDF5 file with the structure below:

        /interferograms        Root level group name
            Attributes         Dictionary for metadata
            /date              2D array of string  in size of (m, 2   ) in YYYYMMDD format for master and slave date
            /bperp             1D array of float32 in size of (m,     ) in meter.
            /unwrapPhase       3D array of float32 in size of (m, l, w) in radian.
            /coherence         3D array of float32 in size of (m, l, w).
            /connectComponent  3D array of int16   in size of (m, l, w). (optional)
            /wrapPhase         3D array of float32 in size of (m, l, w) in radian. (optional)
            /rangeOffset       3D array of float32 in size of (m, l, w). (optional)
            /azimuthOffset     3D array of float32 in size of (m, l, w). (optional)

        Parameters: outputFile : string
                        Name of the HDF5 file for the InSAR stack
                    access_mode : string
                        Access mode of output File, e.g. w, r+
                    box : tuple
                        Subset range in (x0, y0, x1, y1)
        Returns:    outputFile
        '''

        self.outputFile = outputFile
        f = h5py.File(self.outputFile, access_mode)
        print('create HDF5 file {} with {} mode'.format(self.outputFile, access_mode))

        groupName = 'ifgramStack'
        group = f.create_group(groupName)
        print('create group "/{}"'.format(groupName))

        self.pairs = [pair for pair in self.pairsDict.keys()]
        self.dsNames = [dsName for dsName in self.pairsDict[self.pairs[0]].datasetDict.keys()]
        self.get_size(box)

        self.bperp = np.zeros(self.numIfgram)
        ###############################
        # 3D datasets containing unwrapPhase, coherence, connectComponent, wrapPhase, etc.
        for dsName in self.dsNames:
            print('create dataset "/{}/{}"'.format(groupName, dsName))
            ds = group.create_dataset(dsName, shape=(self.numIfgram, self.length, self.width),\
                                      maxshape=(None, self.length, self.width), dtype=dataType, chunks=True)

            progBar = ptime.progress_bar(maxValue=self.numIfgram)
            for i in range(self.numIfgram):
                ifgramObj = self.pairsDict[self.pairs[i]]
                data, metadata = ifgramObj.read(dsName, box=box)
                ds[i,:,:] = data
                self.bperp[i] = ifgramObj.get_perp_baseline()
                progBar.update(i+1, suffix='{}-{}'.format(self.pairs[i][0],self.pairs[i][1]))
            progBar.close()

        ###############################
        # 2D dataset containing master and slave dates of all pairs
        dsDateName = 'date'
        print('create dataset "/{}/{}"'.format(groupName, dsDateName))
        dsDate = group.create_dataset(dsDateName, data=np.array(self.pairs, dtype=np.string_))

        ###############################
        # 1D dataset containing perpendicular baseline of all pairs
        # practice resizable matrix here for update mode
        dsBperpName = 'bperp'
        print('create dataset "/{}/{}"'.format(groupName, dsBperpName))
        if dsBperpName not in group.keys():
            dsBperp = group.create_dataset(dsBperpName, shape=(self.numIfgram,), maxshape=(None,), dtype=dataType)
        else:
            dsBperp = group.get(dsBperpName)
            dsBperp.resize(self.numIfgram, axis=0)
        dsBperp[:] = self.bperp

        ###############################
        # Attributes
        self.get_metadata()
        for key,value in self.metadata.items():
            group.attrs[key] = value

        ###################################
        # 3D datasets for geometry (Lat, Lon, Heigt, Incidence, 
        # Heading, Bperp, ...). For a given platform from a specific 
        # track, a common viewing geometry is assumed. Therfore each 
        # of Lat, Lon, Height, Incidence and Heading can be stored as
        # 2D dataset. Baselines if provided should be 3D. 

        #for dsName in platTrackObj.dsetGeometryNames:
        #    print ('Create dataset for ', dsName)
        #    pairs, length, width = platTrackObj.getSize_geometry(dsName)
        #    numPairs = len(pairs)
        #    dsg = geometryGroup.create_dataset(dsName, shape=(numPairs, length, width),
        #            dtype=dataType) #, chunks=chunk_shape)
        #    for i in range(numPairs):
        #        data, metadata = platTrackObj.pairs[pairs[i]].read(dsName)
        #        dsg[i,:,:] = data

        #for key,value in self.metadata.items():
        #    group.attrs[key] = value

        f.close()
        print('Finished writing to {}'.format(self.outputFile))
        return self.outputFile

    def addDatasets(self, platTrack, fileList, nameList, bands):
        # appends a list of 2D or 3D datsets to the geometry group.
        # Can be lat.rdr, lon.rdr, z.rdr, los.rdr, a mask file, etc 
        import reader
        if fileList is not None:
            f = h5py.File(self.outputFile, 'a')
         
            numDataSets = len(fileList)
            for i in range(numDataSets):
                print('adding ',fileList[i])
                if bands is None:
                    data = reader.read(fileList[i])
                else:
                    data = reader.read(fileList[i] , bands=[bands[i]])
                dsg = f['/geometry'].create_dataset(nameList[i], data=data, shape=data.shape, dtype=data.dtype)
            f.close()


########################################################################################
class ifgram:
    """
    Ifgram object for a single InSAR pair of interferogram. It includes dataset name (family) of:
        'unwrapPhase','coherence','connectComponent','wrapPhase','iono','rangeOffset','azimuthOffset', etc.

    Example:
        from pysar.objects.ifgramStack import ifgram
        datasetDict = {'unwrapPhase':'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.unw',
                       'coherence':'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.cor',
                       'connectComponent':'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.unw.conncomp',
                       'wrapPhase':'$PROJECT_DIR/merged/interferograms/20151220_20160206/filt_fine.int',
                       ...
                      }
        ifgramObj = ifgram(dates=('20160524','20160530'), datasetDict=datasetDict)
        data, atr = ifgramObj.read('unwrapPhase')
    """
    def __init__(self, name='ifgram', dates=None, datasetDict={}, metadata=None):


        self.masterDate, self.slaveDate = dates
        self.datasetDict = datasetDict

        self.platform = None
        self.track = None
        self.processor = None
        # platform, track and processor can get values from metadat if they exist   
        if metadata is not None:
            for key , value in metadata.items():
                setattr(self, key, value)

    def read(self, family=ifgramDatasetNames[0], box=None):
        self.get_metadata()
        data = readfile.read(self.file, box=box)[0]
        return data, self.metadata

    def get_size(self):
        self.file = self.datasetDict[ifgramDatasetNames[0]]
        self.metadata = readfile.read_attribute(self.file)
        self.length = int(self.metadata['LENGTH'])
        self.width = int(self.metadata['WIDTH'])
        return self.length, self.width

    def get_perp_baseline(self):
        self.file = self.datasetDict[ifgramDatasetNames[0]]
        self.metadata = readfile.read_attribute(self.file)
        self.bperp_top = float(self.metadata['P_BASELINE_TOP_HDR'])
        self.bperp_bottom = float(self.metadata['P_BASELINE_BOTTOM_HDR'])
        self.bperp = (self.bperp_top + self.bperp_bottom) / 2.0
        return self.bperp

    def get_metadata(self, family=ifgramDatasetNames[0]):
        self.file = self.datasetDict[family]
        self.metadata = readfile.read_attribute(self.file)
        self.length = int(self.metadata['LENGTH'])
        self.width = int(self.metadata['WIDTH'])

        if self.processor is None:
            ext = self.file.split('.')[-1]
            if os.path.exists(self.file+'.xml'):
                self.processor = 'isce' 
            elif os.path.exists(self.file+'.rsc'):
                self.processor = 'roipac'
            elif os.path.exists(self.file+'.par'):
                self.processor = 'gamma'
            elif ext == 'grd':
                self.processor = 'gmtsar'
            #what for DORIS/SNAP
            elif 'PROCESSOR' in self.metadata.keys():
                self.processor = self.metadata['PROCESSOR']               
            else:
                self.processor = 'isce'
        return self.metadata



########################################################################################
class platformTrack:

    def __init__(self, name='platformTrack'): #, pairDict = None):
        self.pairs = None
         
    def getPairs(self, pairDict, platTrack):
        pairs = pairDict.keys()
        self.pairs = {}
        for pair in pairs:
            if pairDict[pair].platform_track == platTrack:
                self.pairs[pair]=pairDict[pair]

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
        width  = median(width)
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
        self.width  = median(width)

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
            if self.pairs[pair].geometryDict  is not None:
                keys = [k for k in self.pairs[pair].geometryDict.keys()]       
                self.dsetGeometryNames = list(set(self.dsetGeometryNames) | set(keys))




