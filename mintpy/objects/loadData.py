#!/usr/bin/env python3
# Author: Heresh Fattahi

import os
import imp
import sys
import glob
import argparse
import configparser
import datetime
import time
from insarPair import insarPair
from insarStack import insarStack


#################################################################
def create_parser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser(
        description='saving a stack of InSAR pairs to an HDF5 file')
    parser.add_argument('-i', '--input', type=str, nargs='+', dest='input', required=True,
                        help='Directories with the pair directories that includes data for each pair')
    parser.add_argument('-I', '--input_quality', type=str, nargs='+', dest='inputQuality', default=None,
                        help='Directories with the pair directories that includes quality for each pair')
    parser.add_argument('--input_geometry', type=str, nargs='+', dest='inputGeometry', default=None,
                        help='Directories with the pair directories that includes geometry for each pair')

    parser.add_argument('-p', '--name_pattern', type=str, nargs='+', dest='namePattern', required=True,
                        help='name pattern of observation files to be included in the HDF5 file')
    parser.add_argument('-P', '--name_pattern_quality', type=str, nargs='+', dest='namePatternQuality',
                        help='name pattern of quality files to be included in the HDF5 file')
    parser.add_argument('--name_pattern_geometry', type=str, nargs='+', dest='namePatternGeometry',
                        help='name pattern of geometry files to be included in the HDF5 file')

    parser.add_argument('-d', '--observation', type=str, nargs='+', dest='observation', default='unwrapped-phase',
                        help='name of the observation 3D dataset to be stored in the HDF5 file')
    parser.add_argument('-q', '--quality', type=str, nargs='+', dest='quality', default='coherence',
                        help='name of the quality 3D dataset to be stored in the HDF5 file')
    parser.add_argument('-g', '--geometry', type=str, nargs='+', dest='geometry',
                        help='name of the geometry 3D dataset to be stored in the HDF5 file')

    parser.add_argument('-f', '--file_list', type=str, nargs='+', dest='fileList', default=None,
                        help='A list of files to be added to the HDF5 file')
    parser.add_argument('-n', '--name_list', type=str, nargs='+', dest='nameList', default=None,
                        help='A list of datset names corresponding to the fileList to be added to the HDF5 file')
    parser.add_argument('-b', '--band_list', type=int, nargs='+', dest='bandList', default=None,
                        help='A list of band numbers corresponding to the fileList to be added to the HDF5 file. If not specified all bands are added')

    parser.add_argument('-o', '--output', type=str, dest='output', required=True,
                        help='output H5 file that inculdes all datasets')
    parser.add_argument('-r', '--reference_pixel', type=str, dest='refPixel', default=None,
                        help='reference pixel')

    return parser


def cmd_line_parse(iargs=None):
    '''
    Command line parser.
    '''

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps


def write2h5(inps):
    # dumping all pairs to an h5 file
    numObservations = len(inps.observation)
    if inps.inputQuality is not None:
        numQuality = len(inps.inputQuality)
    if inps.inputGeometry is not None:
        numGeometry = len(inps.inputGeometry)

    dirs = glob.glob(os.path.join(inps.input[0], '*'))
    pairsDict = {}
    for dir in dirs:
        dates = os.path.basename(dir).split('_')
        try:
            t1 = time.strptime(dates[0], '%Y%m%d')
        except:
            continue
        Time1 = datetime.datetime(t1.tm_year, t1.tm_mon, t1.tm_mday)

        try:
            t2 = time.strptime(dates[1], '%Y%m%d')
        except:
            continue
        Time2 = datetime.datetime(t2.tm_year, t2.tm_mon, t2.tm_mday)
        metadataDict = {'platform': 'platform', 'processor': 'ISCE'}

        #####################################
        # a dictionary of observartions for a given pair. One pair may
        # have several types of observations.
        # example obsDict = {'unwrapped phase': /pathToFile/filt.unw, 'iono':/PathToFile/iono.bil}

        obsDict = {}
        for i in range(numObservations):
            diri = os.path.join(inps.input[i], dates[0] + '_' + dates[1])
            file = glob.glob(os.path.join(diri, inps.namePattern[i]))

            if len(file) > 0:
                if os.path.exists(file[0]):
                    file = file[0]
                    obsDict[inps.observation[i]] = file
        #####################################
        # a dictionary of quality measures for a given pair.
        qualityDict = {}
        if inps.inputQuality is not None:
            for i in range(numQuality):
                diri = os.path.join(inps.inputQuality[i], dates[0] + '_' + dates[1])
                file = glob.glob(os.path.join(diri, inps.namePatternQuality[i]))
                if len(file) > 0:
                    file = file[0]
                    qualityDict[inps.quality[i]] = file

        #####################################
        # a dictionary of geometry file for a given pair.
        geometryDict = {}
        if inps.inputGeometry is not None:

            for i in range(numGeometry):
                diri = os.path.join(inps.inputGeometry[i], dates[0] + '_' + dates[1])
                file = glob.glob(os.path.join(diri, inps.namePatternGeometry[i]))

                if len(file) > 0:
                    file = file[0]
                    geometryDict[inps.geometry[i]] = file

        #pairObj = insarPair(dates=(Time1 , Time2) ,observation = obsDict, metadata=metadataDict)
        pairObj = insarPair(dates=(Time1, Time2), observation=obsDict,
                            quality=qualityDict, geometry=geometryDict, metadata=metadataDict)
        pairObj.get_metadata(inps.observation[0])
        pairsDict[(Time1, Time2)] = pairObj

    ############################################
    #stackObj = insarStack(pairsDict = pairsDict, nonPairsDict = nonPairsDict)
    stackObj = insarStack(pairsDict=pairsDict)
    stackObj.get_platform_tracks()
    outFile = inps.output
    stackObj.save2h5(output=outFile, access_mode='w', ref_pixels={'platform-track': inps.refPixel})

    stackObj.addDatasets('platform-track', inps.fileList, inps.nameList, inps.bandList)
    return outFile


def main(iargs=None):

    inps = cmd_line_parse(iargs)
    if inps.refPixel:
        inps.refPixel = [int(i) for i in inps.refPixel.split(',')]
        inps.refPixel = (inps.refPixel[0], inps.refPixel[1])

    h5File = write2h5(inps)


if __name__ == '__main__':
    '''
    loading a stack of InSAR pairs to and HDF5 file
    '''

    main()
