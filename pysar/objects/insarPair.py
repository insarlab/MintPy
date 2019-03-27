# Author: Heresh Fattahi, Zhang Yunjun

import os
from pysar.utils import readfile


# A dictionary to help with reading the data when more than one band exists
bandsDict = {

    'isce-unwrapPhase': [2],
    'isce-unwrapAmplitude': [1],
    'isce-wrapPhase': [1],
    'isce-azimuthOffset': [1],
    'isce-rangeOffset': [2],
    'isce-offsetSnr': [1],
    'isce-iono': [1]
}


class ifgram:
    """
     ifgram object for a single InSAR pair of interferogram, including unwrapped/wrapped phase, amplitude,
     coherence, offset in az/rg direction, iono, etc.
    """

    def __init__(self, name='insarPair', dates=None, observation={}, quality={}, geometry={}, metadata=None):

        self.masterDate, self.slaveDate = dates
        #######################################

        self.observationsDict = observation
        self.qualityDict = quality
        self.geometryDict = geometry

        self.platform = 'platform'
        self.track = 'track'
        self.processor = None

        # platform, track and processor can get values from metadat if they exist
        if metadata is not None:
            for key, value in metadata.items():
                setattr(self, key, value)

    def read(self, family):

        self.get_metadata(family)
        bands_key = self.processor + '-' + family
        if bands_key in bandsDict.keys():
            bands = bandsDict[bands_key]
        else:
            bands = None
        print(self.file)
        data = reader.read(self.file, self.processor, bands=bands)
        return data, self.metadata

    def get_metadata(self, family):

        if family in self.observationsDict.keys():
            self.file = self.observationsDict[family]
        elif family in self.qualityDict.keys():
            self.file = self.qualityDict[family]
        elif family in self.geometryDict.keys():
            self.file = self.geometryDict[family]
        else:
            self.file = ''
        # else:
        #   '''error message'''
        ############################
        # if the processor is not known, find the processor based on the file extension
        if self.processor is None:

            ext = self.file.split('.')[-1]
            if os.path.exists(self.file+'.xml'):
                self.processor = 'isce'
            elif os.path.exists(self.file+'.rsc'):
                self.processor = 'roipac'
            elif os.path.exists(self.file+'.par'):
                self.processor = 'gamma'
            elif os.path.exists(self.file+'.par'):
                self.processor = 'GAMMA'
            elif ext == 'grd':
                self.processor = 'GMT'
            # what for DORIS

        if os.path.exists(self.file):
            self.metadata = reader.read_metadata(self.file, self.processor)
            self.length = int(self.metadata['LENGTH'])
            self.width = int(self.metadata['WIDTH'])
        else:
            self.metadata = {}
            self.length = 0
            self.width = 0

        if self.platform is None and 'platform' in self.metadata.keys():
            self.platform = self.metadata['platform']

        if self.track is None and 'track' in self.metadata.keys():
            self.track = self.metadata['track']

        self.platform_track = self.platform + '-' + self.track

        if self.processor is None:
            if 'processor' in self.metadata.keys():
                self.processor = self.metadata['processor']
            else:
                self.processor = 'ISCE'
