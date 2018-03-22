# Author: Heresh Fattahi, Zhang Yunjun

import os, sys, glob
import h5py
import numpy as np
from pysar.utils import readfile, writefile, datetime as ptime

dataType = np.float32

class ifgramStack:
    '''
    Stack of Interferograms object for a set of interferograms and coherence from the same platform and track.
    Attributes are saved in the root level.
    It contains a "interferograms" group and some datasets: date12, bperp, unwrapIfgram, coherence, ...
    
    /                    Root level
    Attributes           Dictionary for metadata
    interferograms       group name
        /date            2D array of string in size of (m, 2) in YYYYMMDD format, 1st column for master, 2nd for slave.
        /bperp           1D array of float32 in size of (m,) in meter.
        /unwrapIfgram    3D array of float32 in size of (m, l, w) in radian.
        /coherence       3D array of float32 in size of (m, l, w).
        /connComp        3D array of int16 in size of (m, l, w). (optional)
        /wrapIfgram      3D array of float32 in size of (m, l, w) in radian. (optional)
        /rangeOffset     3D array of float32 in size of (m, l, w). (optional)
        /azimuthOffset   3D array of float32 in size of (m, l, w). (optional)
        ...
    All pairs for a given platform-track are required to have the same size.
    Pairs of different platforms-tracks may have different size.
    '''

    def __init__(self, name='ifgramStack', pairsDict=None):
        self.pairsDict = pairsDict
        self.datasetNameDict = {}

    def get_size(self, box=None):
        self.numIfgram = len(self.pairsDict)
        ifgramObj = [v for v in self.pairsDict.values()][0]
        self.metadata = ifgramObj.get_metadata()
        if box:
            self.length = box[3] - box[1]
            self.width = box[2] - box[0]
        else:
            self.length = ifgramObj.length
            self.width = ifgramObj.width
        return self.numIfgram, self.length, self.width

    def save2h5(self, outputFile='ifgramStack.h5', access_mode='w', box=None):
        '''
        h5OutName : Name of the HDF5 file for the InSAR stack

        platform_tracks : A list containing the platform_tracks to be stored 
                          in the HDF5 file. If None all platform_tracks are
                          exctracted from pairs. If pairs does not contain information
                          about the platform_tracks, then all pairs in the pairsDict are 
                          considered from a single platform single track.

        ref_pixel       : A dictionary containing refernce pixels for each platform_track.
                          eg: ref_pixDict = {'Sentinel1A/Track144':(500,300) , 'Sentinel1A/Track100':(440,700)}     
                              first pixel is in y direction (lines) and second pixel in x direction (columns)
        '''
        self.get_size(box)
        groupName = 'ifgramStack'

        pairsTime = [pair for pair in self.pairsDict.keys()]
        self.outputFile = outputFile
        self.h5file = h5py.File(outputFile, access_mode)
        group = self.h5file.create_group(groupName)

        self.dsNames = [dsName for dsName in self.pairsDict[pairsTime[0]].datasetDict.keys()]
        #import pdb; pdb.set_trace()
        for dsName in self.dsNames:
            print ('Create dataset "/{}/{}"'.format(groupName, dsName))
            ds = group.create_dataset(dsName, shape=(self.numIfgram, self.length, self.width),\
                                      maxshape=(None, self.length, self.width), dtype=dataType, chunks=True)

            progBar = ptime.progress_bar(maxValue=self.numIfgram)
            for i in range(self.numIfgram):
                ifgramObj = self.pairsDict[pairsTime[i]]
                data, metadata = readfile.read(ifgramObj.datasetDict[dsName], box=box)
                ds[i,:,:] = data
                #import pdb; pdb.set_trace()
                mDate, sDate = [i.strftime('%Y%m%d') for i in pairsTime[i]]
                progBar.update(i+1, suffix='{}-{}'.format(mDate, sDate))
            progBar.close()
        
        ###############################
        # A 2D dataset containing a 2D array of strings. First column 
        # is the master time and second column the slave time of pairs.
        #if len(platTrackObj.dsetObservationNames) > 0:
        #    piars_idx = vstack((masterTimes,slaveTimes)).T
        #    dspairs = group.create_dataset('pairs_idx', data=piars_idx, dtype=piars_idx.dtype)

        ###################################
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
#
        #    for i in range(numPairs):
        #        data, metadata = platTrackObj.pairs[pairs[i]].read(dsName)
        #        dsg[i,:,:] = data

        #for key,value in self.metadata.items():
        #    group.attrs[key] = value

        self.h5file.close()

    def addDatasets(self, platTrack, fileList, nameList, bands):
        # appends a list of 2D or 3D datsets to the geometry group. Can be lat.rdr, lon.rdr, z.rdr, los.rdr, a mask file, etc 
        import reader
        if fileList is not None:
            self.h5file = h5py.File(self.outputFile, 'a')
         
            numDataSets = len(fileList)
            for i in range(numDataSets):
                print('adding ',fileList[i])
                if bands is None:
                    data = reader.read(fileList[i])
                else:
                    data = reader.read(fileList[i] , bands=[bands[i]])
                dsg = self.h5file['/geometry'].create_dataset(nameList[i], data=data, shape=data.shape, dtype=data.dtype)
            self.h5file.close()


    


# A dictionary to help with reading the data when more than one band exists
bandsDict = {

             'isce-unwrapPhase' : [2] , 
             'isce-unwrapAmplitude' : [1], 
             'isce-wrapPhase' : [1], 
             'isce-azimuthOffset':[1],
             'isce-rangeOffset':[2],             
             'isce-offsetSnr':[1],
             'isce-iono':[1]
}

class ifgram:
    """
     ifgram object for a single InSAR pair of interferogram, including unwrapped/wrapped phase, amplitude,
     coherence, offset in az/rg direction, iono, etc.
    """
    def __init__(self, name='ifgram', dates=None, datasetDict={}, metadata=None):

        self.masterDate, self.slaveDate = dates
        #######################################
        self.datasetDict = datasetDict
        self.platform = None
        self.track = None
        self.processor = None

        # platform, track and processor can get values from metadat if they exist   
        if metadata is not None:
            for key , value in metadata.items():
                setattr(self, key, value)

    def read(self, family):
        self.get_metadata(family)
        data = readfile.read(self.file)[0]
        return data, self.metadata

    def get_metadata(self, family=None):
        if not family:
            family = [k for k in self.datasetDict.keys()][0]
        self.file = self.datasetDict[family]

        # if the processor is not known, find the processor based on the file extension
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

        if os.path.exists(self.file):
            self.metadata = readfile.read_attribute(self.file)
            self.length = int(self.metadata['LENGTH'])
            self.width = int(self.metadata['WIDTH'])
        else:
            self.metadata = {}
            self.length = 0
            self.width = 0

        if self.processor is None:
            if 'PROCESSOR' in self.metadata.keys():
                self.processor = self.metadata['PROCESSOR']               
            else:
                self.processor = 'isce'
        return self.metadata
        
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




