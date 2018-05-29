#!/usr/bin/env python3
# Author: Heresh Fattahi, Zhang Yunjun

import os
import sys
import glob
import argparse
import datetime
import time

from pysar.defaults import auto_path
from pysar.objects import ifgramDatasetNames, geometryDatasetNames, ifgramStack, geometry, sensor
from pysar.objects.insarobj import ifgramDict, ifgramStackDict, geometryDict
from pysar.utils import readfile, ptime, utils as ut
from pysar import subset


#################################################################
datasetName2templateKey = {'unwrapPhase'     : 'pysar.load.unwFile',
                           'coherence'       : 'pysar.load.corFile',
                           'connectComponent': 'pysar.load.connCompFile',
                           'wrapPhase'       : 'pysar.load.intFile',
                           'height'          : 'pysar.load.demFile',
                           'latitude'        : 'pysar.load.lookupYFile',
                           'longitude'       : 'pysar.load.lookupXFile',
                           'azimuthCoord'    : 'pysar.load.lookupYFile',
                           'rangeCoord'      : 'pysar.load.lookupXFile',
                           'incidenceAngle'  : 'pysar.load.incAngleFile',
                           'headingAngle'    : 'pysar.load.headAngleFile',
                           'shadowMask'      : 'pysar.load.shadowMaskFile',
                           'bperp'           : 'pysar.load.bperpFile'
                           }

DEFAULT_TEMPLATE = """template:
########## 1. Load Data (--load to exit after this step)
{}\n
{}\n
{}\n
""".format(auto_path.isceAutoPath,
           auto_path.roipacAutoPath,
           auto_path.gammaAutoPath)

TEMPLATE = """template:
########## 1. Load Data (--load to exit after this step)
## auto - automatic path pattern for Univ of Miami file structure
## load_data.py -H to check more details and example inputs.
pysar.load.processor      = auto  #[isce,roipac,gamma,], auto for isce
pysar.load.updateMode     = auto  #[yes / no], auto for yes, skip re-loading if HDF5 files are complete
##---------interferogram datasets:
pysar.load.unwFile        = auto  #[path2unw_file]
pysar.load.corFile        = auto  #[path2cor_file]
pysar.load.connCompFile   = auto  #[path2conn_file]
pysar.load.intFile        = auto  #[path2int_file]
##---------geometry datasets:
pysar.load.demFile        = auto  #[path2hgt_file]
pysar.load.lookupYFile    = auto  #[path2lat_file]]
pysar.load.lookupXFile    = auto  #[path2lon_file]
pysar.load.incAngleFile   = auto  #[path2los_file]
pysar.load.headAngleFile  = auto  #[path2los_file]
pysar.load.shadowMaskFile = auto  #[path2shadow_file]
pysar.load.bperpFile      = auto  #[path2bperp_file]

## 1.1 Subset (optional)
## if both yx and lalo are specified, use lalo option unless a) no lookup file AND b) dataset is in radar coord
pysar.subset.lalo     = auto    #[31.5:32.5,130.5:131.0 / no], auto for no
pysar.subset.yx       = auto    #[1800:2000,700:800 / no], auto for no
pysar.subset.tightBox = auto    #[yes / no], auto for yes, tight bounding box for files in geo coord
"""

NOTE = """NOTE:
  unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the master and slave date, either in file name or folder name.
"""

EXAMPLE = """example:
  load_data.py -t GalapagosSenDT128.tempalte
  load_data.py -t pysarApp_template.txt
  load_data.py -H #Show example input template for ISCE/ROI_PAC/GAMMA products
"""


def create_parser():
    """Create command line parser."""
    parser = argparse.ArgumentParser(description='Saving a stack of Interferograms to an HDF5 file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+NOTE+'\n'+EXAMPLE)
    parser.add_argument('-H', dest='print_example_template', action='store_true',
                        help='Print/Show the example template file for loading.')
    parser.add_argument('-t', '--template', type=str,
                        dest='template_file', help='template file with path info.')

    parser.add_argument('--project', type=str, dest='PROJECT_NAME',
                        help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor', type=str, dest='processor',
                        choices={'isce', 'roipac', 'gamma', 'doris', 'gmtsar'},
                        help='InSAR processor/software of the file', default='isce')
    parser.add_argument('--enforce', '-f', dest='updateMode', action='store_false',
                        help='Disable the update mode, or skip checking dataset already loaded.')
    parser.add_argument('--compression', choices={'gzip', 'lzf', None}, default=None,
                        help='compress loaded geometry while writing HDF5 file, default: None.')

    parser.add_argument('-o', '--output', type=str, nargs=3, dest='outfile',
                        default=['./INPUTS/ifgramStack.h5',
                                 './INPUTS/geometryRadar.h5', './INPUTS/geometryGeo.h5'],
                        help='output HDF5 file')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.print_example_template:
        sys.exit(DEFAULT_TEMPLATE)

    if not inps.template_file:
        parser.print_usage()
        print(('{}: error: the following arguments are required:'
               ' -t/--template'.format(os.path.basename(__file__))))
        print('{} -H to show the example template file'.format(os.path.basename(__file__)))
        sys.exit(1)

    inps.outfile = [os.path.abspath(i) for i in inps.outfile]
    inps.outdir = os.path.dirname(inps.outfile[0])

    return inps


#################################################################
def read_inps2dict(inps):
    """Read input Namespace object info into inpsDict"""
    # Read input info into inpsDict
    inpsDict = vars(inps)
    inpsDict['PLATFORM'] = None

    # Read template file
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)
    # FA 1/18: insert general options from template file as pysar.* options
    if 'processor' in template.keys():
        template['pysar.load.processor'] = template['processor']

    prefix = 'pysar.load.'
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        value = template[prefix+key]
        if key in ['processor', 'updateMode', 'compression']:
            inpsDict[key] = template[prefix+key]
        elif value:
            inpsDict[prefix+key] = template[prefix+key]

    if inpsDict['compression'] == False:
        inpsDict['compression'] = None

    # PROJECT_NAME --> PLATFORM
    (inpsDict['PLATFORM'],
     inpsDict['PROJECT_NAME']) = sensor.project_name2sensor([inpsDict['PROJECT_NAME'],
                                                             inps.template_file])
    if inpsDict['PLATFORM']:
        print('Find PLATFORM from PROJECT_NAME as: {}'.format(inpsDict['PLATFORM']))

    # Here to insert code to check default file path for miami user
    if (auto_path.autoPath
        and 'SCRATCHDIR' in os.environ
        and inpsDict['PROJECT_NAME'] is not None):
        print(('check auto path setting for Univ of Miami users'
               ' for processor: {}'.format(inpsDict['processor'])))
        if inpsDict['processor'] == 'isce':
            inpsDict = auto_path.get_auto_path4isce(inpsDict['PROJECT_NAME'], inpsDict)
        elif inpsDict['processor'] == 'roipac':
            inpsDict = auto_path.get_auto_path4roipac(inpsDict['PROJECT_NAME'], inpsDict)
        elif inpsDict['processor'] == 'gamma':
            inpsDict = auto_path.get_auto_path4gamma(inpsDict['PROJECT_NAME'], inpsDict)
    return inpsDict


def read_subset_box(inpsDict):
    # Read subset info from template
    inpsDict['box'] = None
    inpsDict['box4geo_lut'] = None
    pix_box, geo_box = subset.read_subset_template2box(inpsDict['template_file'])

    # Grab required info to read input geo_box into pix_box
    try:
        lookupFile = [glob.glob(str(inpsDict['pysar.load.lookupYFile']))[0],
                      glob.glob(str(inpsDict['pysar.load.lookupXFile']))[0]]
    except:
        lookupFile = None

    try:
        pathKey = [i for i in datasetName2templateKey.values()
                   if i in inpsDict.keys()][0]
        file = glob.glob(str(inpsDict[pathKey]))[0]
        atr = readfile.read_attribute(file)
    except:
        atr = dict()

    geocoded = None
    if 'Y_FIRST' in atr.keys():
        geocoded = True
    else:
        geocoded = False

    # Check conflict
    if geo_box and not geocoded and lookupFile is None:
        geo_box = None
        print(('WARNING: pysar.subset.lalo is not supported'
               ' if 1) no lookup file AND'
               '    2) radar/unkonwn coded dataset'))
        print('\tignore it and continue.')
    if not geo_box and not pix_box:
        return inpsDict

    # geo_box --> pix_box
    if geo_box is not None:
        if geocoded:
            pix_box = (0, 0, 0, 0)
            pix_box[0:3:2] = ut.coord_geo2radar(geo_box[0:3:2], atr, 'lon')
            pix_box[1:4:2] = ut.coord_geo2radar(geo_box[1:4:2], atr, 'lat')
        else:
            pix_box = subset.bbox_geo2radar(geo_box, atr, lookupFile)
        print('input bounding box of interest in lalo: {}'.format(geo_box))
    print('box to read for datasets in y/x: {}'.format(pix_box))

    # Get box for geocoded lookup table (for gamma/roipac)
    box4geo_lut = None
    if lookupFile is not None:
        atrLut = readfile.read_attribute(lookupFile[0])
        if not geocoded and 'Y_FIRST' in atrLut.keys():
            geo_box = subset.bbox_radar2geo(pix_box, atr, lookupFile)
            box4geo_lut = subset.bbox_geo2radar(geo_box, atrLut)
            print('box to read for geocoded lookup file in y/x: {}'.format(box4geo_lut))

    inpsDict['box'] = pix_box
    inpsDict['box4geo_lut'] = box4geo_lut
    return inpsDict


def read_inps_dict2ifgram_stack_dict_object(inpsDict):
    """Read input arguments into dict of ifgramStackDict object"""
    # inpsDict --> dsPathDict
    print('-'*50)
    print('searching interferometric pairs info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(datasetName2templateKey.keys())])
    dsPathDict = {}
    dsNumDict = {}
    for dsName in [i for i in ifgramDatasetNames
                   if i in datasetName2templateKey.keys()]:
        key = datasetName2templateKey[dsName]
        if key in inpsDict.keys():
            files = sorted(glob.glob(str(inpsDict[key])))
            if len(files) > 0:
                dsPathDict[dsName] = files
                dsNumDict[dsName] = len(files)
                print('{:<{width}}: {path}'.format(dsName,
                                                   width=maxDigit,
                                                   path=inpsDict[key]))

    # Check required dataset
    dsName0 = ifgramDatasetNames[0]
    if dsName0 not in dsPathDict.keys():
        print('WARNING: No reqired {} data files found!'.format(dsName0))
        return None

    # Check number of files for all dataset types
    dsNumList = list(dsNumDict.values())
    if any(i != dsNumList[0] for i in dsNumList):
        print('WARNING: Not all types of dataset have the same number of files:')
        for key, value in dsNumDict.items():
            print('number of {:<{width}}: {num}'.format(key,
                                                        width=maxDigit,
                                                        num=value))
    print('number of files per type: {}'.format(dsNumList[0]))

    # Check data dimension for all files

    # dsPathDict --> pairsDict --> stackObj
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
                dsPath2 = [i for i in dsPathDict[dsName]
                           if all(d[2:8] in i for d in dates)]
                if len(dsPath2) > 0:
                    ifgramPathDict[dsName] = dsPath2[0]
                else:
                    print('WARNING: {} file missing for pair {}'.format(dsName, dates))
        ifgramObj = ifgramDict(dates=tuple(dates),
                               datasetDict=ifgramPathDict)
        pairsDict[tuple(dates)] = ifgramObj

    if len(pairsDict) > 0:
        stackObj = ifgramStackDict(pairsDict=pairsDict)
    else:
        stackObj = None
    return stackObj


def read_inps_dict2geometry_dict_object(inpsDict):

    # eliminate dsName by processor
    if inpsDict['processor'] in ['isce', 'doris']:
        datasetName2templateKey.pop('azimuthCoord')
        datasetName2templateKey.pop('rangeCoord')
    elif inpsDict['processor'] in ['roipac', 'gamma']:
        datasetName2templateKey.pop('latitude')
        datasetName2templateKey.pop('longitude')
    else:
        print('Un-recognized InSAR processor: {}'.format(inpsDict['processor']))

    # inpsDict --> dsPathDict
    print('-'*50)
    print('searching geometry files info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(datasetName2templateKey.keys())])
    dsPathDict = {}
    for dsName in [i for i in geometryDatasetNames
                   if i in datasetName2templateKey.keys()]:
        key = datasetName2templateKey[dsName]
        if key in inpsDict.keys():
            files = sorted(glob.glob(str(inpsDict[key])))
            if len(files) > 0:
                if dsName == 'bperp':
                    bperpDict = {}
                    for file in files:
                        date = ptime.yyyymmdd(os.path.basename(os.path.dirname(file)))
                        bperpDict[date] = file
                    dsPathDict[dsName] = bperpDict
                    print('{:<{width}}: {path}'.format(dsName,
                                                       width=maxDigit,
                                                       path=inpsDict[key]))
                    print('number of bperp files: {}'.format(len(list(bperpDict.keys()))))
                else:
                    dsPathDict[dsName] = files[0]
                    print('{:<{width}}: {path}'.format(dsName,
                                                       width=maxDigit,
                                                       path=files[0]))

    # Check required dataset
    dsName0 = geometryDatasetNames[0]
    if dsName0 not in dsPathDict.keys():
        print('WARNING: No reqired {} data files found!'.format(dsName0))
        # sys.exit(1)

    # metadata
    ifgramRadarMetadata = None
    ifgramKey = datasetName2templateKey[ifgramDatasetNames[0]]
    if ifgramKey in inpsDict.keys():
        ifgramFiles = glob.glob(str(inpsDict[ifgramKey]))
        if len(ifgramFiles) > 0:
            atr = readfile.read_attribute(ifgramFiles[0])
            if 'Y_FIRST' not in atr.keys():
                ifgramRadarMetadata = atr.copy()

    # dsPathDict --> dsGeoPathDict + dsRadarPathDict
    dsNameList = list(dsPathDict.keys())
    dsGeoPathDict = {}
    dsRadarPathDict = {}
    for dsName in dsNameList:
        if dsName == 'bperp':
            atr = readfile.read_attribute(next(iter(dsPathDict[dsName].values())))
        else:
            atr = readfile.read_attribute(dsPathDict[dsName])
        if 'Y_FIRST' in atr.keys():
            dsGeoPathDict[dsName] = dsPathDict[dsName]
        else:
            dsRadarPathDict[dsName] = dsPathDict[dsName]

    geomRadarObj = None
    geomGeoObj = None
    if len(dsRadarPathDict) > 0:
        geomRadarObj = geometryDict(processor=inpsDict['processor'],
                                    datasetDict=dsRadarPathDict,
                                    extraMetadata=ifgramRadarMetadata)
    if len(dsGeoPathDict) > 0:
        geomGeoObj = geometryDict(processor=inpsDict['processor'],
                                  datasetDict=dsGeoPathDict,
                                  extraMetadata=None)
    return geomRadarObj, geomGeoObj


def update_object(outFile, inObj, box, updateMode=True):
    """Do not write h5 file if: 1) h5 exists and readable,
                                2) it contains all date12 from ifgramStackDict,
                                            or all datasets from geometryDict"""
    updateFile = True
    if updateMode and not ut.update_file(outFile, check_readable=True):
        if inObj.name == 'ifgramStack':
            outObj = ifgramStack(outFile)
            if (outObj.get_size() == inObj.get_size(box=box) 
                    and sorted(outObj.get_date12_list(dropIfgram=False)) == sorted(inObj.get_date12_list())):
                print(('All date12   exists in file {} with same size as required,'
                       ' no need to re-load.'.format(outFile)))
                updateFile = False

        elif inObj.name == 'geometry':
            outObj = geometry(outFile)
            outObj.open(print_msg=False)
            if (outObj.get_size() == inObj.get_size(box=box)
                    and all(i in outObj.datasetNames for i in inObj.get_dataset_list())):
                print(('All datasets exists in file{} with same size as required,'
                       ' no need to re-load.'.format(outFile)))
                updateFile = False
    return updateFile


def prepare_metadata(inpsDict):
    prepCmd0 = None
    if inpsDict['processor'] == 'gamma':
        prepCmd0 = 'prep_gamma.py'
    elif inpsDict['processor'] == 'roipac':
        prepCmd0 = 'prep_roipac.py'

    if prepCmd0:
        print('-'*50)
        print(('prepare metadata files for {} products before loading'
               ' into PySAR'.format(inpsDict['processor'])))
        for key in inpsDict.keys():
            if (key.startswith('pysar.load.') and key.endswith('File')
                    and len(glob.glob(str(inpsDict[key]))) > 0):
                prepCmd = '{} {}'.format(prepCmd0, inpsDict[key])
                print(prepCmd)
                os.system(prepCmd)
    return


def print_write_setting(inpsDict):
    updateMode = inpsDict['updateMode']
    comp = inpsDict['compression']
    print('-'*50)
    print('updateMode : {}'.format(updateMode))
    print('compression: {}'.format(comp))
    box = inpsDict['box']
    boxGeo = inpsDict['box4geo_lut']

    return updateMode, comp, box, boxGeo


def get_extra_metadata(inpsDict):
    """Extra metadata to be written into stack file"""
    extraDict = {}
    for key in ['PROJECT_NAME', 'PLATFORM']:
        if inpsDict[key]:
            extraDict[key] = inpsDict[key]
    return extraDict


#################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not os.path.isdir(inps.outdir):
        os.makedirs(inps.outdir)
        print('create directory: {}'.format(inps.outdir))

    inpsDict = read_inps2dict(inps)
    inpsDict = read_subset_box(inpsDict)
    prepare_metadata(inpsDict)

    stackObj = read_inps_dict2ifgram_stack_dict_object(inpsDict)
    geomRadarObj, geomGeoObj = read_inps_dict2geometry_dict_object(inpsDict)

    updateMode, comp, box, boxGeo = print_write_setting(inpsDict)
    if stackObj and update_object(inps.outfile[0], stackObj, box, updateMode=updateMode):
        print('-'*50)
        extraDict = get_extra_metadata(inpsDict)
        stackObj.write2hdf5(outputFile=inps.outfile[0],
                            access_mode='w',
                            box=box,
                            compression=comp,
                            extra_metadata=extraDict)

    if geomRadarObj and update_object(inps.outfile[1], geomRadarObj, box, updateMode=updateMode):
        print('-'*50)
        geomRadarObj.write2hdf5(outputFile=inps.outfile[1],
                                access_mode='w',
                                box=box,
                                compression=comp)

    if geomGeoObj and update_object(inps.outfile[2], geomGeoObj, boxGeo, updateMode=updateMode):
        print('-'*50)
        geomGeoObj.write2hdf5(outputFile=inps.outfile[2],
                              access_mode='w',
                              box=boxGeo,
                              compression=comp)

    return inps.outfile


#################################################################
if __name__ == '__main__':
    """
    loading a stack of InSAR pairs to and HDF5 file
    """
    main()
