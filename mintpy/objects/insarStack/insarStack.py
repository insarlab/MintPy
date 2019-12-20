# Author: Heresh Fattahi

import os
import sys
import h5py
import insarPair as insarPair
from numpy import median, float32, complex64, vstack

chunk_shape = (128, 128)
dataType = float32
#dataType = complex64


class insarStack:
    """    
    InsarStack object for a stack of InSAR data (multi-platform multi-track data).
    Method save2h5 creates a HDF5 file to store the stack data.
    Each platform-track is a group with three sub groups : observation, quality, geometry.
    Each sub-group may have different datasets. Some common datasets are:

    observations (3D) :  unwrapped interferograms, RangeOffset, AzimuthOffset, ...
    quality  (3D)    :  coherence, uncertainty, ...
    geometry (2D or 3D)    :  incidence, heading angle, latitute, longitude, ...              

    All pairs for a given platform-track are required to have the same size.
    Pairs of different platforms-tracks may have different size.
    """

    def __init__(self, name='insarStack', pairsDict=None):
        self.pairs = pairsDict

    def save2h5(self, output='data.h5', access_mode='w', platform_tracks=None,
                ref_pixels=None, ref_pixel_method='average_coherence'):
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
        self.output = output
        self.h5file = h5py.File(output, access_mode)
        self.platform_tracks = platform_tracks
        self.ref_pixels = ref_pixels
        if self.platform_tracks is None:
            self.get_platform_tracks()

        for platTrack in self.platform_tracks:
            print('platform-track : ', platTrack)
            group = self.h5file.create_group(platTrack)
            obsGroup = group.create_group('observations')
            qualityGroup = group.create_group('quality')
            geometryGroup = group.create_group('geometry')
            ################################
            # A class object for the platformTrack
            platTrackObj = platformTrack()
            platTrackObj.getPairs(self.pairs, platTrack)
            platTrackObj.getSize()
            platTrackObj.getDatasetNames()
            ###############################
            # Get the reference pixel for a given platform/track
            if self.ref_pixels is not None:
                platTrackObj.ref_pixel = self.ref_pixels[platTrack]
            else:
                platTrackObj.ref_pixel = None
            print('refernce pixel: ', platTrackObj.ref_pixel)
            ###############################
            # 3D datasets for observation quality (Coherence, uncertainty, ...)
            pairs = [pair for pair in platTrackObj.pairs.keys()]
            for dsName in platTrackObj.dsetQualityNames:
                print('Create dataset for ', dsName)
                dsq = qualityGroup.create_dataset(dsName, shape=(
                    platTrackObj.numPairs, platTrackObj.length, platTrackObj.width), dtype=dataType)

                masterTimes = [None]*platTrackObj.numPairs
                slaveTimes = [None]*platTrackObj.numPairs
                for i in range(platTrackObj.numPairs):
                    data, metadata = platTrackObj.pairs[pairs[i]].read(dsName)
                    dsq[i, :, :] = data
                    master, slave = pairs[i]
                    masterTimes[i] = master.strftime(
                        '%Y-%m-%d %H:%M:%S').encode('utf8')
                    slaveTimes[i] = slave.strftime(
                        '%Y-%m-%d %H:%M:%S').encode('utf8')

            ###############################
            # store the pair times as a 2D dataset
            if len(platTrackObj.dsetQualityNames) > 0:
                piars_idx = vstack((masterTimes, slaveTimes)).T
                dsq = qualityGroup.create_dataset(
                    'pairs_idx', data=piars_idx, dtype=piars_idx.dtype)
            ###############################
            # if the reference pixel is not given let's choose a pixel with maximum average coherence
            # if platTrackObj.ref_pixel is None:
            #    platTrackObj.ref_pixel = self.choose_ref_pixel(platTrack , method == 'average_coherence')

            ###############################
            # 3D datasets for observations (possible datasets: unwrapped-phase, RangeOffset, AzimuthOffset, unwrapped-amplitude, etc)
            # There should be no limitation for storing any other possible observations.

            pairs = [pair for pair in platTrackObj.pairs.keys()]

            for dsName in platTrackObj.dsetObservationNames:
                print('Create dataset for ', dsName)
                dso = obsGroup.create_dataset(dsName, shape=(platTrackObj.numPairs, platTrackObj.length, platTrackObj.width),
                                              dtype=dataType)  # , chunks=chunk_shape)

                masterTimes = [None]*platTrackObj.numPairs
                slaveTimes = [None]*platTrackObj.numPairs
                for i in range(platTrackObj.numPairs):
                    print(pairs[i])
                    data, metadata = platTrackObj.pairs[pairs[i]].read(dsName)
                    if platTrackObj.ref_pixel:
                        dso[i, :, :] = data - data[0, platTrackObj.ref_pixel[0], platTrackObj.ref_pixel[1]]
                    else:
                        dso[i, :, :] = data
                    master, slave = pairs[i]
                    masterTimes[i] = master.strftime('%Y-%m-%d %H:%M:%S').encode("ascii", "ignore")
                    slaveTimes[i] = slave.strftime('%Y-%m-%d %H:%M:%S').encode("ascii", "ignore")

            ###############################
            # A 2D dataset containing a 2D array of strings. First column
            # is the master time and second column the slave time of pairs.
            if len(platTrackObj.dsetObservationNames) > 0:
                piars_idx = vstack((masterTimes, slaveTimes)).T
                dspairs = group.create_dataset('pairs_idx', data=piars_idx, dtype=piars_idx.dtype)
            ###################################
            for key, value in metadata.items():
                obsGroup.attrs[key] = value

            ###################################
            # 3D datasets for geometry (Lat, Lon, Heigt, Incidence,
            # Heading, Bperp, ...). For a given platform from a specific
            # track, a common viewing geometry is assumed. Therfore each
            # of Lat, Lon, Height, Incidence and Heading can be stored as
            # 2D dataset. Baselines if provided should be 3D.

            for dsName in platTrackObj.dsetGeometryNames:
                print('Create dataset for ', dsName)
                pairs, length, width = platTrackObj.getSize_geometry(dsName)
                numPairs = len(pairs)
                dsg = geometryGroup.create_dataset(dsName, shape=(numPairs, length, width),
                                                   dtype=dataType)  # , chunks=chunk_shape)

                for i in range(numPairs):
                    data, metadata = platTrackObj.pairs[pairs[i]].read(dsName)
                    dsg[i, :, :] = data

            for key, value in metadata.items():
                geometryGroup.attrs[key] = value

        self.h5file.close()

    def addDatasets(self, platTrack, fileList, nameList, bands):
        # appends a list of 2D or 3D datsets to the geometry group. Can be lat.rdr, lon.rdr, z.rdr, los.rdr, a mask file, etc
        import reader
        if fileList is not None:
            self.h5file = h5py.File(self.output, 'a')

            numDataSets = len(fileList)
            for i in range(numDataSets):
                print('adding ', fileList[i])
                if bands is None:
                    data = reader.read(fileList[i])
                else:
                    data = reader.read(fileList[i], bands=[bands[i]])

                dsg = self.h5file['/'+platTrack+'/geometry'].create_dataset(
                    nameList[i], data=data, shape=data.shape, dtype=data.dtype)
            self.h5file.close()

    def get_platform_tracks(self):

        self.platform_tracks = []
        for pair in self.pairs.keys():
            if self.pairs[pair].platform_track not in self.platform_tracks:
                self.platform_tracks.append(self.pairs[pair].platform_track)

    # def loadh5(self, platform_track , groupName='observation', datasetName='unwrapped', method = , method_par, )

    #     method     :  chunck       , block       , all
    #     method_par :  Chunck_size  , block_size  ,

#    def choose_reference_pixel(self, platTrack , method):

        # compute average coherence of the 3D dataset
        # find the pixel with maximum value


#    def time_baseline_timeseries():


##################################

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
