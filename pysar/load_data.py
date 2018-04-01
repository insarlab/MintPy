#!/usr/bin/env python3
# Author: Heresh Fattahi, Zhang Yunjun

import os, sys, glob
import argparse
import datetime, time

import pysar
from pysar.utils import readfile, datetime as ptime
from pysar.utils.sensor import sensors
from pysar.objects import ifgramDatasetNames, geometryDatasetNames
from pysar.objects.insarobj import ifgram, ifgramStack, geometry
from pysar.defaults import default_path


#################################################################
datasetName2templateKey={'unwrapPhase'     :'pysar.load.unwFile',
                         'coherence'       :'pysar.load.corFile',
                         'connectComponent':'pysar.load.connCompFile',
                         'wrapPhase'       :'pysar.load.intFile',
                         'height'          :'pysar.load.demFile',
                         'latitude'        :'pysar.load.lookupYFile',
                         'longitude'       :'pysar.load.lookupXFile',
                         'azimuthCoord'    :'pysar.load.lookupYFile',
                         'rangeCoord'      :'pysar.load.lookupXFile',
                         'incidenceAngle'  :'pysar.load.incAngleFile',
                         'headingAngle'    :'pysar.load.headAngleFile',
                         'shadowMask'      :'pysar.load.shadowMaskFile'
                         }

DEFAULT_TEMPLATE='''template:
## 1. Load Data (--load to exit after this step)
{}\n
{}\n
{}\n
##----------Subset Range Setting
pysar.subset.yx   = auto
pysar.subset.lalo = auto
'''.format(default_path.isceAutoPath,\
           default_path.roipacAutoPath,\
           default_path.gammaAutoPath)

TEMPLATE='''template:
pysar.load.processor      = auto  #[isce,            roipac,              gamma,             ], auto for isce
pysar.load.unwFile        = auto  #[filt*.unw,       filt*rlks_c*.unw,    diff*rlks.unw      ]
pysar.load.corFile        = auto  #[filt*.cor,       filt*rlks.cor,       *filt*rlks.cor     ]
pysar.load.connCompFile   = auto  #[filt*.unw.conn*, filt*rlks_snap*.byt, None               ]
pysar.load.intFile        = auto  #[filt*.int,       filt*rlks.int,       diff*rlks.int      ]

pysar.load.demFile        = auto  #[hgt.rdr,         radar*.hgt,          sim_*rlks.rdc.hgt  ]
pysar.load.lookupYFile    = auto  #[lat.rdr,         geomap*.trans,       sim*rlks.UTM_TO_RDC]
pysar.load.lookupXFile    = auto  #[lon.rdr,         geomap*.trans,       sim*rlks.UTM_TO_RDC]
pysar.load.incAngleFile   = auto  #[los.rdr,         None,                None               ]
pysar.load.headAngleFile  = auto  #[los.rdr,         None,                None               ]
pysar.load.shadowMaskFile = auto  #[shadowMask.rdr,  None,                None               ]

pysar.subset.yx   = auto
'''

NOTE='''NOTE:
  unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the master and slave date, either in file name or folder name.
'''

EXAMPLE='''example:
  load_data.py -t GalapagosSenDT128.tempalte
  load_data.py -H #Show example input template for ISCE/ROI_PAC/GAMMA products
'''

def createParser():
    '''Create command line parser.'''
    parser = argparse.ArgumentParser(description='Saving a stack of Interferograms to an HDF5 file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+NOTE+'\n'+EXAMPLE)

    parser.add_argument('-H', dest='print_example_template', action='store_true',\
                        help='Print/Show the example template file for loading.')

    parser.add_argument('-t','--template', type=str, nargs='+', dest='template_file', help='template file with path info.')

    parser.add_argument('--project', type=str, dest='project_name', help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor', type=str, dest='processor', choices={'isce','roipac','gamma','doris','gmtsar'},\
                        help='InSAR processor/software of the file')
    #parser.add_argument('-x', type=int, nargs=2, dest='subset_x', metavar=('X_MIN','X_MAX'),\
    #                    help='Subset range in x/range direction')
    #parser.add_argument('-y', type=int, nargs=2, dest='subset_y', metavar=('Y_MIN','Y_MAX'),\
    #                    help='Subset range in y/azimuth direction')

    parser.add_argument('--enforce', dest='update_mode', action='store_false',\
                        help='Disable the update mode, or skip checking dataset already loaded. [not implemented yet]')
    parser.add_argument('-o','--output', type=str, nargs=3, dest='outfile',\
                        default=['./INPUTS/ifgramStack.h5','./INPUTS/geometryRadar.h5','./INPUTS/geometryGeo.h5'],\
                        help='output HDF5 file')
    return parser


def cmdLineParse(iargs = None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)

    if inps.print_example_template:
        print(DEFAULT_TEMPLATE)
        sys.exit(1)

    if not inps.template_file:
        parser.print_usage()
        print('{}: error: the following arguments are required: -t/--template'.format(os.path.basename(__file__)))
        print('{} -H to show the example template file'.format(os.path.basename(__file__)))
        sys.exit(1)

    return inps


#################################################################
def read_inps2dict(inps):
    '''Read input Namespace object info into inpsDict'''
    ## Read input info into inpsDict
    inpsDict = dict()
    inpsDict['PROJECT_NAME'] = None
    inpsDict['processor'] = None
    inpsDict['platform'] = None

    for file in inps.template_file:
        inpsDict.update(readfile.read_template(file))
        if any(i in file for i in sensors):
            inpsDict['PROJECT_NAME'] = os.path.splitext(os.path.basename(file))[1]

    ##Here to insert code to check default file path for miami user
    # Check 1) SCRATCHDIR exists, 2) pysar.auto_path is True and 3) template['PROJECT_NAME'] is not None

    if len(glob.glob(inpsDict['pysar.load.unwFile'])) == 0:
        print('ERROR: No required .unw file found with pattern: {}'.format(inpsDict['pysar.load.unwFile']))

    if not inpsDict['processor']:
        inpsDict['processor'] = inpsDict['pysar.load.processor']
    if inps.processor:
        inpsDict['processor'] = inps.processor.lower()
    if not inpsDict['processor']:
        inpsDict['processor'] = 'isce'

    if inps.project_name:
        inpsDict['PROJECT_NAME'] = inps.project_name
    if inpsDict['PROJECT_NAME']:
        try:
            inpsDict['platform'] = [i for i in sensors if i in inpsDict['PROJECT_NAME']][0].upper()
        except:
            pass

    return inpsDict


def read_subset_box(inpsDict):
    file0 = glob.glob(inpsDict['pysar.load.unwFile'])[0]
    atr = readfile.read_attribute(file0)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    box = (0,0,width,length)

    key = 'pysar.subset.yx'
    if inpsDict[key].lower() not in ['auto','no']:
        try:
            sub = [i.strip() for i in inpsDict[key].split(',')]
            y0, y1 = sorted([int(i.strip()) for i in sub[0].split(':')])
            x0, x1 = sorted([int(i.strip()) for i in sub[1].split(':')])
            box = (x0, y0, x1, y1)
        except:
            print('ERROR: input subset info {} is not in required format x0:x1,y0:y1'.format(inpsDict[key]))
            sys.exit(1)        

    print('box of input  files: {}'.format((0,0,width,length)))
    print('box of data to read: {}'.format(box))
    return box


def read_inps_dict2ifgram_stack_object(inpsDict):
    '''Read input arguments into dict of ifgramStack object'''
    ########## inpsDict --> dsPathDict
    print('-'*50)
    print('searching interferometric pairs info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(datasetName2templateKey.keys())])
    dsPathDict = {}
    dsNumDict = {}
    for dsName in ifgramDatasetNames:
        if dsName in datasetName2templateKey.keys():
            key = datasetName2templateKey[dsName]
            files = glob.glob(inpsDict[key])
            if len(files) > 0:
                dsPathDict[dsName] = files
                dsNumDict[dsName] = len(files)
                print('{:<{width}}: {path}'.format(dsName, width=maxDigit, path=inpsDict[key]))

    ##Check required dataset
    dsName0 = ifgramDatasetNames[0]
    if dsName0 not in dsPathDict.keys():
        print('ERROR: No reqired {} data files found!'.format(dsName0))
        sys.exit(1)

    ##Check number of files for all dataset types
    dsNumList = list(dsNumDict.values())
    if any(i != dsNumList[0] for i in dsNumList):
        print('ERROR: Not all types of dataset have the same number of files:')
        for key, value in dsNumDict.items():
            print('number of {:<{width}}: {num}'.format(key, width=maxDigit, num=value))
        sys.exit(1)
    print('number of files per type: {}'.format(dsNumList[0]))

    ##Check data dimension for all files
    

    ########## dsPathDict --> pairsDict --> stackObj
    dsNameList = list(dsPathDict.keys())
    pairsDict = {}
    for dsPath in dsPathDict[dsName0]:
        dates = ptime.yyyymmdd(readfile.read_attribute(dsPath)['DATE12'].split('-'))

        #####################################
        # A dictionary of data files for a given pair.
        # One pair may have several types of dataset.
        # example ifgramPathDict = {'unwrapPhase': /pathToFile/filt.unw, 'iono':/PathToFile/iono.bil}
        # All path of data file must contain the master and slave date, either in file name or folder name.

        ifgramPathDict = {}
        for i in range(len(dsNameList)):
            dsName = dsNameList[i]
            dsPath1 = dsPathDict[dsName][i]
            if all(d[2:8] in dsPath1 for d in dates):
                ifgramPathDict[dsName] = dsPath1
            else:
                dsPath2 = [i for i in dsPathDict[dsName] if all(d[2:8] in i for d in dates)]
                if len(dsPath2)>0:
                    ifgramPathDict[dsName] = dsPath2[0]
                else:
                    print('WARNING: {} file missing for pair {}'.format(dsName, dates))
        ifgramObj = ifgram(dates=tuple(dates), datasetDict=ifgramPathDict)
        pairsDict[tuple(dates)] = ifgramObj

    if len(pairsDict)>0:
        stackObj = ifgramStack(pairsDict=pairsDict)
    else:
        stackObj = None
    return stackObj


def read_inps_dict2geometry_object(inpsDict):

    ########## eliminate dsName by processor
    if inpsDict['processor'] in ['isce','doris']:
        datasetName2templateKey.pop('azimuthCoord')
        datasetName2templateKey.pop('rangeCoord')
    elif inpsDict['processor'] in ['roipac','gamma']:
        datasetName2templateKey.pop('latitude')
        datasetName2templateKey.pop('longitude')
    else:
        print('Un-recognized InSAR processor: {}'.format(inpsDict['processor']))

    ########## inpsDict --> dsPathDict
    print('-'*50)
    print('searching geometry files info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(datasetName2templateKey.keys())])
    dsPathDict = {}
    for dsName in geometryDatasetNames:
        if dsName in datasetName2templateKey.keys():
            key = datasetName2templateKey[dsName]
            files = glob.glob(inpsDict[key])
            if len(files) > 0:
                dsPathDict[dsName] = files[0]
                print('{:<{width}}: {path}'.format(dsName, width=maxDigit, path=files[0]))

    ##Check required dataset
    dsName0 = geometryDatasetNames[0]
    if dsName0 not in dsPathDict.keys():
        print('ERROR: No reqired {} data files found!'.format(dsName0))
        sys.exit(1)

    ########## metadata
    ifgramRadarMetadata = None
    key = datasetName2templateKey[ifgramDatasetNames[0]]
    files = glob.glob(inpsDict[key])
    if len(files)>0:
        atr = readfile.read_attribute(files[0])
        if 'Y_FIRST' not in atr.keys():
            ifgramRadarMetadata = atr.copy()

    ########## dsPathDict --> dsGeoPathDict + dsRadarPathDict
    dsNameList = list(dsPathDict.keys())
    dsGeoPathDict = {}
    dsRadarPathDict = {}
    for dsName in dsNameList:
        atr = readfile.read_attribute(dsPathDict[dsName])
        if 'Y_FIRST' in atr.keys():
            dsGeoPathDict[dsName] = dsPathDict[dsName]
        else:
            dsRadarPathDict[dsName] = dsPathDict[dsName]

    geomRadarObj = None
    geomGeoObj = None
    if len(dsRadarPathDict)>0:
        geomRadarObj = geometry(processor=inpsDict['processor'], datasetDict=dsRadarPathDict, ifgramMetadata=ifgramRadarMetadata)
    if len(dsGeoPathDict)>0:
        geomGeoObj = geometry(processor=inpsDict['processor'], datasetDict=dsGeoPathDict, ifgramMetadata=None)
    return geomRadarObj, geomGeoObj


#################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    inps.outfile = [os.path.abspath(i) for i in inps.outfile]
    inps.outdir = os.path.dirname(inps.outfile[0])
    if not os.path.isdir(inps.outdir):
        os.makedirs(inps.outdir)
        print('create directory: {}'.format(inps.outdir))

    inpsDict = read_inps2dict(inps)
    box = read_subset_box(inpsDict)

    stackObj = read_inps_dict2ifgram_stack_object(inpsDict)
    geomRadarObj, geomGeoObj = read_inps_dict2geometry_object(inpsDict)

    if stackObj:
        print('-'*50)
        stackObj.save2h5(outputFile=inps.outfile[0], access_mode='w', box=box)
    if geomRadarObj:
        print('-'*50)
        geomRadarObj.save2h5(outputFile=inps.outfile[1], access_mode='w', box=box)
    if geomGeoObj:
        print('-'*50)
        geomGeoObj.save2h5(outputFile=inps.outfile[2], access_mode='w')

    return inps.outfile


#################################################################
if __name__ == '__main__' :
    '''
    loading a stack of InSAR pairs to and HDF5 file
    '''
    main()
