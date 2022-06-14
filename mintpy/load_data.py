#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import glob
import time
import argparse
import warnings

from mintpy.defaults import auto_path
from mintpy.defaults.template import get_template_content
from mintpy.objects import (
    geometryDatasetNames,
    geometry,
    ifgramDatasetNames,
    ifgramStack,
    sensor,
)
from mintpy.objects.stackDict import (
    geometryDict,
    ifgramStackDict,
    ifgramDict,
)
from mintpy.utils import readfile, ptime, utils as ut
from mintpy import subset


#################################################################
PROCESSOR_LIST = ['isce', 'aria', 'hyp3', 'gmtsar', 'snap', 'gamma', 'roipac', 'cosicorr']

IFGRAM_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.unwFile',
    'coherence'       : 'mintpy.load.corFile',
    'connectComponent': 'mintpy.load.connCompFile',
    'wrapPhase'       : 'mintpy.load.intFile',
    'magnitude'       : 'mintpy.load.magFile',
}

ION_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.ionUnwFile',
    'coherence'       : 'mintpy.load.ionCorFile',
    'connectComponent': 'mintpy.load.ionConnCompFile',
}

OFFSET_DSET_NAME2TEMPLATE_KEY = {
    'azimuthOffset'   : 'mintpy.load.azOffFile',
    'azimuthOffsetStd': 'mintpy.load.azOffStdFile',
    'rangeOffset'     : 'mintpy.load.rgOffFile',
    'rangeOffsetStd'  : 'mintpy.load.rgOffStdFile',
    'offsetSNR'       : 'mintpy.load.offSnrFile',
}

GEOMETRY_DSET_NAME2TEMPLATE_KEY = {
    'height'          : 'mintpy.load.demFile',
    'latitude'        : 'mintpy.load.lookupYFile',
    'longitude'       : 'mintpy.load.lookupXFile',
    'azimuthCoord'    : 'mintpy.load.lookupYFile',
    'rangeCoord'      : 'mintpy.load.lookupXFile',
    'incidenceAngle'  : 'mintpy.load.incAngleFile',
    'azimuthAngle'    : 'mintpy.load.azAngleFile',
    'shadowMask'      : 'mintpy.load.shadowMaskFile',
    'waterMask'       : 'mintpy.load.waterMaskFile',
    'bperp'           : 'mintpy.load.bperpFile',
}


DEFAULT_TEMPLATE = """template:
########## 1. Load Data (--load to exit after this step)
{}\n
{}\n
{}""".format(
    auto_path.AUTO_PATH_GAMMA,
    auto_path.AUTO_PATH_ISCE_STRIPMAP,
    auto_path.AUTO_PATH_ISCE_TOPS,
)

TEMPLATE = get_template_content('load_data')

NOTE = """NOTE:
  For interferogram, unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the reference and secondary date, either in file name or folder name.
"""

EXAMPLE = """example:
  # show example template file for ISCE/ROI_PAC/GAMMA products
  load_data.py -H

  load_data.py -t smallbaselineApp.cfg
  load_data.py -t smallbaselineApp.cfg GalapagosSenDT128.txt --project GalapagosSenDT128

  # load ionosphere stack
  smallbaselineApp.py SaltonSeaSenDT173.txt -g
  load_data.py -t smallbaselineApp.cfg --ion

  # load geometry ONLY
  smallbaselineApp.py SaltonSeaSenDT173.txt -g
  load_data.py -t smallbaselineApp.cfg --geom
"""


def create_parser():
    """Create command line parser."""
    parser = argparse.ArgumentParser(description='Saving a stack of Interferograms to an HDF5 file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+NOTE+'\n'+EXAMPLE)
    parser.add_argument('-H', dest='print_example_template', action='store_true',
                        help='Print/Show the example template file for loading.')
    parser.add_argument('-t', '--template', dest='template_file', type=str, nargs='+',
                        help='template file(s) with path info.')

    # write single file
    single = parser.add_mutually_exclusive_group(required=False)
    single.add_argument('--ion','--ionosphere', dest='only_load_ionosphere', action='store_true',
                        help='Switch to load the ionospheric pairs into ionStack.h5 file.')
    single.add_argument('--geom','--geometry', dest='only_load_geometry', action='store_true',
                        help='Load the geometry file(s) ONLY.')

    # options from template file name & content
    parser.add_argument('--project', type=str, dest='PROJECT_NAME',
                        help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor', type=str, dest='processor', choices=PROCESSOR_LIST,
                        help='InSAR processor/software of the file', default='isce')
    parser.add_argument('--enforce', '-f', dest='updateMode', action='store_false',
                        help='Disable the update mode, or skip checking dataset already loaded.')
    parser.add_argument('--compression', choices={'gzip', 'lzf', None}, default=None,
                        help='compress loaded geometry while writing HDF5 file, default: None.')

    # output
    parser.add_argument('-o', '--output', type=str, nargs=3, dest='outfile',
                        default=['./inputs/ifgramStack.h5',
                                 './inputs/geometryRadar.h5',
                                 './inputs/geometryGeo.h5'],
                        help='output HDF5 file')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check --template option
    if inps.template_file:
        pass

    elif inps.print_example_template:
        print(DEFAULT_TEMPLATE)
        sys.exit(0)

    else:
        parser.print_usage()
        print(('{}: error: one of the following arguments are required:'
               ' -t/--template, -H'.format(os.path.basename(__file__))))
        print('{} -H to show the example template file'.format(os.path.basename(__file__)))
        sys.exit(1)

    inps.outfile = [os.path.abspath(i) for i in inps.outfile]
    inps.outdir = os.path.dirname(inps.outfile[0])

    if inps.only_load_ionosphere and os.path.basename(inps.outfile[0]) == 'ifgramStack.h5':
        inps.outfile[0] = os.path.join(inps.outdir, 'ionStack.h5')
        print(f'load ionosphere only --> set the default output file name to: {inps.outfile[0]}.')
    return inps


#################################################################
def read_inps2dict(inps):
    """Read input Namespace object info into iDict.

    It grab the following contents into iDict
    1. inps & all template files
    2. configurations: processor, autoPath, updateMode, compression, x/ystep
    3. extra metadata: PLATFORM, PROJECT_NAME,
    4. translate autoPath

    Parameters: inps  - namespace, input arguments from command line & template file
    Returns:    iDict - dict,      input arguments from command line & template file
    """
    # Read input info into iDict
    iDict = vars(inps)
    iDict['PLATFORM'] = None

    # Read template file
    template = {}
    for fname in inps.template_file:
        temp = readfile.read_template(fname)
        temp = ut.check_template_auto_value(temp)
        template.update(temp)
    for key, value in template.items():
        iDict[key] = value
    if 'processor' in template.keys():
        template['mintpy.load.processor'] = template['processor']

    prefix = 'mintpy.load.'
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        value = template[prefix+key]
        if key in ['processor', 'autoPath', 'updateMode', 'compression']:
            iDict[key] = template[prefix+key]
        elif key in ['xstep', 'ystep']:
            iDict[key] = int(template[prefix+key])
        elif value:
            iDict[prefix+key] = template[prefix+key]
    print('processor : {}'.format(iDict['processor']))

    if iDict['compression'] == False:
        iDict['compression'] = None

    iDict['xstep'] = iDict.get('xstep', 1)
    iDict['ystep'] = iDict.get('ystep', 1)

    # PROJECT_NAME --> PLATFORM
    if not iDict['PROJECT_NAME']:
        cfile = [i for i in list(inps.template_file) if os.path.basename(i) != 'smallbaselineApp.cfg']
        iDict['PROJECT_NAME'] = sensor.project_name2sensor_name(cfile)[1]

    msg = 'SAR platform/sensor : '
    sensor_name = sensor.project_name2sensor_name(str(iDict['PROJECT_NAME']))[0]
    if sensor_name:
        msg += str(sensor_name)
        iDict['PLATFORM'] = str(sensor_name)
    else:
        msg += 'unknown from project name "{}"'.format(iDict['PROJECT_NAME'])
    print(msg)

    # update file path with auto
    if iDict.get('autoPath', False):
        print('use auto path defined in mintpy.defaults.auto_path for options in auto')
        iDict = auto_path.get_auto_path(processor=iDict['processor'],
                                        work_dir=os.path.dirname(iDict['outdir']),
                                        template=iDict)

    # copy dset_name2template_key info into iDict
    if inps.only_load_ionosphere:
        iDict['dset_name2template_key'] = {
            **ION_DSET_NAME2TEMPLATE_KEY,
            **GEOMETRY_DSET_NAME2TEMPLATE_KEY}
    else:
        iDict['dset_name2template_key'] = {
            **IFGRAM_DSET_NAME2TEMPLATE_KEY,
            **OFFSET_DSET_NAME2TEMPLATE_KEY,
            **GEOMETRY_DSET_NAME2TEMPLATE_KEY}

    return iDict


def read_subset_box(iDict):
    """read the following items:
    geocoded
    box
    box4geo_lut
    """
    # Read subset info from template
    iDict['box'] = None
    iDict['box4geo_lut'] = None
    pix_box, geo_box = subset.read_subset_template2box(iDict['template_file'][0])

    # Grab required info to read input geo_box into pix_box
    try:
        lookupFile = [glob.glob(str(iDict['mintpy.load.lookupYFile']))[0],
                      glob.glob(str(iDict['mintpy.load.lookupXFile']))[0]]
    except:
        lookupFile = None

    try:
        pathKey = [i for i in iDict['dset_name2template_key'].values()
                   if i in iDict.keys()][0]
        file = glob.glob(str(iDict[pathKey]))[0]
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
        print(('WARNING: mintpy.subset.lalo is not supported'
               ' if 1) no lookup file AND'
               '    2) radar/unkonwn coded dataset'))
        print('\tignore it and continue.')

    if not geo_box and not pix_box:
        # adjust for the size inconsistency problem in SNAP geocoded products
        # ONLY IF there is no input subset
        # Use the min bbox if files size are different
        if iDict['processor'] == 'snap':
            fnames = ut.get_file_list(iDict['mintpy.load.unwFile'])
            pix_box = update_box4files_with_inconsistent_size(fnames)

        if not pix_box:
            return iDict

    # geo_box --> pix_box
    coord = ut.coordinate(atr, lookup_file=lookupFile)
    if geo_box is not None:
        pix_box = coord.bbox_geo2radar(geo_box)
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print('input bounding box of interest in lalo: {}'.format(geo_box))
    print('box to read for datasets in y/x: {}'.format(pix_box))

    # Get box for geocoded lookup table (for gamma/roipac)
    box4geo_lut = None
    if lookupFile is not None:
        atrLut = readfile.read_attribute(lookupFile[0])
        if not geocoded and 'Y_FIRST' in atrLut.keys():
            geo_box = coord.bbox_radar2geo(pix_box)
            box4geo_lut = ut.coordinate(atrLut).bbox_geo2radar(geo_box)
            print('box to read for geocoded lookup file in y/x: {}'.format(box4geo_lut))

    iDict['geocoded'] = geocoded
    iDict['box'] = pix_box
    iDict['box4geo_lut'] = box4geo_lut
    return iDict


#################################################################
def update_box4files_with_inconsistent_size(fnames):
    """Check the size (row / column number) of a list of files
    For SNAP geocoded products has one line missing in some interferograms, Andre, 2019-07-16
    Parameters: fnames  - list of path for interferogram files
    Returns:    pix_box - None if all files are in same size
                          (0, 0, min_width, min_length) if not.
    """
    atr_list = [readfile.read_attribute(fname) for fname in fnames]
    length_list = [int(atr['LENGTH']) for atr in atr_list]
    width_list = [int(atr['WIDTH']) for atr in atr_list]
    if any(len(set(i)) for i in [length_list, width_list]):
        min_length = min(length_list)
        min_width = min(width_list)
        pix_box = (0, 0, min_width, min_length)

        # print out warning message
        msg = '\n'+'*'*80
        msg += '\nWARNING: NOT all input unwrapped interferograms have the same row/column number!'
        msg += '\nMinimum size is: ({}, {})'.format(min_length, min_width)
        msg += '\n'+'-'*30
        msg += '\nThe following dates have different size:'

        for i in range(len(fnames)):
            if length_list[i] != min_length or width_list[i] != min_width:
                msg += '\n\t{}\t({}, {})'.format(atr_list[i]['DATE12'], length_list[i], width_list[i])

        msg += '\n'+'-'*30
        msg += '\nAssuming the interferograms above have:'
        msg += '\n\textra line(s) at the bottom OR'
        msg += '\n\textra column(s) at the right'
        msg += '\nContinue to load data using subset of the minimum size.'
        msg += '\n'+'*'*80+'\n'
        print(msg)
    else:
        pix_box = None
    return pix_box


def skip_files_with_inconsistent_size(dsPathDict, pix_box=None, dsName='unwrapPhase'):
    """Skip files by removing the file path from the input dsPathDict."""
    atr_list = [readfile.read_attribute(fname) for fname in dsPathDict[dsName]]
    length_list = [int(atr['LENGTH']) for atr in atr_list]
    width_list = [int(atr['WIDTH']) for atr in atr_list]

    # Check size requirements
    drop_inconsistent_files = False
    if any(len(set(size_list)) > 1 for size_list in [length_list, width_list]):
        if pix_box is None:
            drop_inconsistent_files = True
        else:
            # if input subset is within the min file sizes: do NOT drop
            max_box_width, max_box_length = pix_box[2:4]
            if max_box_length > min(length_list) or max_box_width > min(width_list):
                drop_inconsistent_files = True

    # update dsPathDict
    if drop_inconsistent_files:
        common_length = ut.most_common(length_list)
        common_width = ut.most_common(width_list)

        # print out warning message
        msg = '\n'+'*'*80
        msg += '\nWARNING: NOT all input unwrapped interferograms have the same row/column number!'
        msg += '\nThe most common size is: ({}, {})'.format(common_length, common_width)
        msg += '\n'+'-'*30
        msg += '\nThe following dates have different size:'

        dsNames = list(dsPathDict.keys())
        date12_list = [atr['DATE12'] for atr in atr_list]
        num_drop = 0
        for i in range(len(date12_list)):
            if length_list[i] != common_length or width_list[i] != common_width:
                date12 = date12_list[i]
                dates = ptime.yyyymmdd(date12.split('-'))
                # update file list for all datasets
                for dsName in dsNames:
                    fnames = [i for i in dsPathDict[dsName]
                              if all(d[2:8] in i for d in dates)]
                    if len(fnames) > 0:
                        dsPathDict[dsName].remove(fnames[0])
                msg += '\n\t{}\t({}, {})'.format(date12, length_list[i], width_list[i])
                num_drop += 1

        msg += '\n'+'-'*30
        msg += '\nSkip loading the above interferograms ({}).'.format(num_drop)
        msg += '\nContinue to load the rest interferograms ({}).'.format(len(date12_list) - num_drop)
        msg += '\n'+'*'*80+'\n'
        print(msg)
    return dsPathDict


def read_inps_dict2ifgram_stack_dict_object(iDict):
    """Read input arguments into ifgramStackDict object.

    Parameters: iDict    - dict, input arguments from command line & template file
    Returns:    stackObj - ifgramStackDict object or None
    """
    if iDict['only_load_geometry']:
        return None
    elif iDict['only_load_ionosphere']:
        ds_type = 'ionospheric'
    else:
        ds_type = 'interferometric'

    # iDict --> dsPathDict
    print('-'*50)
    print(f'searching {ds_type} pairs info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(iDict['dset_name2template_key'].keys())])
    dsPathDict = {}
    for dsName in [i for i in ifgramDatasetNames
                   if i in iDict['dset_name2template_key'].keys()]:
        key = iDict['dset_name2template_key'][dsName]
        if key in iDict.keys():
            files = sorted(glob.glob(str(iDict[key])))
            if len(files) > 0:
                dsPathDict[dsName] = files
                print('{:<{width}}: {path}'.format(dsName,
                                                   width=maxDigit,
                                                   path=iDict[key]))

    # Check 1: required dataset
    dsName0s = ['unwrapPhase', 'rangeOffset', 'azimuthOffset']
    dsName0 = [i for i in dsName0s if i in dsPathDict.keys()]
    if len(dsName0) == 0:
        print('WARNING: No reqired {} data files found!'.format(dsName0s))
        return None
    else:
        dsName0 = dsName0[0]

    # Check 2: data dimension for unwrapPhase files
    dsPathDict = skip_files_with_inconsistent_size(dsPathDict,
                                                   pix_box=iDict['box'],
                                                   dsName=dsName0)

    # Check 3: number of files for all dataset types
    # dsPathDict --> dsNumDict
    dsNumDict = {}
    for key in dsPathDict.keys():
        num_file = len(dsPathDict[key])
        dsNumDict[key] = num_file
        print('number of {:<{width}}: {num}'.format(key, width=maxDigit, num=num_file))

    dsNumList = list(dsNumDict.values())
    if any(i != dsNumList[0] for i in dsNumList):
        msg = 'WARNING: NOT all types of dataset have the same number of files.'
        msg += ' -> skip interferograms with missing files and continue.'
        print(msg)
        #raise Exception(msg)

    # dsPathDict --> pairsDict --> stackObj
    dsNameList = list(dsPathDict.keys())

    #####################################
    # A dictionary of data file paths for a list of pairs, e.g.:
    # pairsDict = {
    #     ('date1', 'date2') : ifgramPathDict1,
    #     ('date1', 'date3') : ifgramPathDict2,
    #     ...,
    # }

    pairsDict = {}
    for i, dsPath0 in enumerate(dsPathDict[dsName0]):
        # date string used in the file/dir path
        # YYYYDDD       for gmtsar [modern Julian date]
        # YYYYMMDDTHHMM for uavsar
        # YYYYMMDD      for all the others
        date6s = readfile.read_attribute(dsPath0)['DATE12'].replace('_','-').split('-')
        if iDict['processor'] == 'gmtsar':
            date12MJD = os.path.basename(os.path.dirname(dsPath0))
        else:
            date12MJD = None

        #####################################
        # A dictionary of data file paths for a given pair.
        # One pair may have several types of dataset, e.g.:
        # ifgramPathDict1 = {
        #     'unwrapPhase': /dirPathToFile/filt_fine.unw,
        #     'coherence'  : /dirPathToFile/filt_fine.cor,
        #     ...
        # }
        # All path of data file must contain the reference and secondary date, in file/dir name.

        ifgramPathDict = {}
        for dsName in dsNameList:
            # search the matching data file for the given date12
            # 1st guess: file in the same order as the one for dsName0
            dsPath1 = dsPathDict[dsName][i]
            if (all(d6 in dsPath1 for d6 in date6s)
                    or (date12MJD and date12MJD in dsPath1)):
                ifgramPathDict[dsName] = dsPath1

            else:
                # 2nd guess: any file in the list
                dsPath2 = [p for p in dsPathDict[dsName]
                           if (all(d6 in p for d6 in date6s)
                                   or (date12MJD and date12MJD in dsPath1))]

                if len(dsPath2) > 0:
                    ifgramPathDict[dsName] = dsPath2[0]
                else:
                    print('WARNING: {:>18} file missing for pair {}'.format(dsName, date6s))

        # initiate ifgramDict object
        ifgramObj = ifgramDict(datasetDict=ifgramPathDict)

        # update pairsDict object
        date8s = ptime.yyyymmdd(date6s)
        pairsDict[tuple(date8s)] = ifgramObj

    if len(pairsDict) > 0:
        stackObj = ifgramStackDict(pairsDict=pairsDict, dsName0=dsName0)
    else:
        stackObj = None
    return stackObj


def read_inps_dict2geometry_dict_object(iDict):
    """Read input arguments into geometryDict object(s).

    Parameters: iDict        - dict, input arguments from command line & template file
    Returns:    geomRadarObj - geometryDict object in radar coordinates or None
                geomGeoObj   - geometryDict object in geo   coordinates or None
    """
    if iDict['only_load_ionosphere']:
        return None, None

    # eliminate lookup table dsName for input files in radar-coordinates
    if iDict['processor'] in ['isce', 'doris']:
        # for processors with lookup table in radar-coordinates, remove azimuth/rangeCoord
        iDict['dset_name2template_key'].pop('azimuthCoord')
        iDict['dset_name2template_key'].pop('rangeCoord')
    elif iDict['processor'] in ['roipac', 'gamma']:
        # for processors with lookup table in geo-coordinates, remove latitude/longitude
        iDict['dset_name2template_key'].pop('latitude')
        iDict['dset_name2template_key'].pop('longitude')
    elif iDict['processor'] in ['aria', 'gmtsar', 'hyp3', 'snap', 'cosicorr']:
        # for processors with geocoded products support only, do nothing for now.
        # check again when adding products support in radar-coordiantes
        pass
    else:
        print('Un-recognized InSAR processor: {}'.format(iDict['processor']))

    # iDict --> dsPathDict
    print('-'*50)
    print('searching geometry files info')
    print('input data files:')
    maxDigit = max([len(i) for i in list(iDict['dset_name2template_key'].keys())])
    dsPathDict = {}
    for dsName in [i for i in geometryDatasetNames
                   if i in iDict['dset_name2template_key'].keys()]:
        key = iDict['dset_name2template_key'][dsName]
        if key in iDict.keys():
            files = sorted(glob.glob(str(iDict[key])))
            if len(files) > 0:
                if dsName == 'bperp':
                    bperpDict = {}
                    for file in files:
                        date = ptime.yyyymmdd(os.path.basename(os.path.dirname(file)))
                        bperpDict[date] = file
                    dsPathDict[dsName] = bperpDict
                    print('{:<{width}}: {path}'.format(dsName,
                                                       width=maxDigit,
                                                       path=iDict[key]))
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

    # metadata
    ifgramMetaGeo = None
    ifgramMetaRadar = None
    ifgramKey = iDict['dset_name2template_key']['unwrapPhase']
    if ifgramKey in iDict.keys():
        ifgramFiles = glob.glob(str(iDict[ifgramKey]))
        if len(ifgramFiles) > 0:
            atr = readfile.read_attribute(ifgramFiles[0])
            if 'Y_FIRST' in atr.keys():
                ifgramMetaGeo = atr.copy()
            else:
                ifgramMetaRadar = atr.copy()

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
        geomRadarObj = geometryDict(processor=iDict['processor'],
                                    datasetDict=dsRadarPathDict,
                                    extraMetadata=ifgramMetaRadar)
    if len(dsGeoPathDict) > 0:
        geomGeoObj = geometryDict(processor=iDict['processor'],
                                  datasetDict=dsGeoPathDict,
                                  extraMetadata=ifgramMetaGeo)
    return geomRadarObj, geomGeoObj


#################################################################
def run_or_skip(outFile, inObj, box, updateMode=True, xstep=1, ystep=1):
    """Check if re-writing is necessary.

    Do not write HDF5 file if ALL the following meet:
        1. HDF5 file exists and is readable,
        2. HDF5 file constains all the datasets and in the same size
        3. For ifgramStackDict, HDF5 file contains all date12.

    Parameters: outFile    - str, path to the output HDF5 file
                inObj      - ifgramStackDict or geometryDict, object to write
                box        - tuple of int, bounding box in (x0, y0, x1, y1)
                updateMode - bool
                x/ystep    - int
    Returns:    flag       - str, run or skip
    """

    flag = 'run'
    if updateMode and ut.run_or_skip(outFile, check_readable=True) == 'skip':
        if inObj.name == 'ifgramStack':
            in_size = inObj.get_size(box=box, xstep=xstep, ystep=ystep)[1:]
            in_dset_list = inObj.get_dataset_list()
            in_date12_list = inObj.get_date12_list()

            outObj = ifgramStack(outFile)
            outObj.open(print_msg=False)
            out_size = (outObj.length, outObj.width)
            out_dset_list = outObj.datasetNames
            out_date12_list = outObj.date12List

            if (out_size == in_size
                    and set(in_dset_list).issubset(set(out_dset_list))
                    and set(in_date12_list).issubset(set(out_date12_list))):
                print(('All date12   exists in file {} with same size as required,'
                       ' no need to re-load.'.format(os.path.basename(outFile))))
                flag = 'skip'

        elif inObj.name == 'geometry':
            in_size = inObj.get_size(box=box, xstep=xstep, ystep=ystep)
            in_dset_list = inObj.get_dataset_list()

            outObj = geometry(outFile)
            outObj.open(print_msg=False)
            out_size = (outObj.length, outObj.width)
            out_dset_list = outObj.datasetNames

            if (out_size == in_size
                    and set(in_dset_list).issubset(set(out_dset_list))):
                print(('All datasets exists in file {} with same size as required,'
                       ' no need to re-load.'.format(os.path.basename(outFile))))
                flag = 'skip'

    return flag


def prepare_metadata(iDict):
    """Prepare metadata via prep_{insar_processor}.py scripts."""

    processor = iDict['processor']
    script_name = 'prep_{}.py'.format(processor)
    print('-'*50)
    print('prepare metadata files for {} products'.format(processor))

    if processor in ['gamma', 'hyp3', 'roipac', 'snap', 'cosicorr']:
        # import prep_module
        if processor == 'gamma':
            from mintpy import prep_gamma as prep_module
        elif processor == 'hyp3':
            from mintpy import prep_hyp3 as prep_module
        elif processor == 'roipac':
            from mintpy import prep_roipac as prep_module
        elif processor == 'snap':
            from mintpy import prep_snap as prep_module
        elif processor == 'cosicorr':
            from mintpy import prep_cosicorr as prep_module

        # run prep_{processor} module
        for key in [i for i in iDict.keys()
                    if (i.startswith('mintpy.load.')
                        and i.endswith('File')
                        and i != 'mintpy.load.metaFile')]:
            if len(glob.glob(str(iDict[key]))) > 0:
                # print command line
                script_name = '{}.py'.format(os.path.basename(prep_module.__name__).split('.')[-1])
                iargs = [iDict[key]]
                if processor == 'gamma' and iDict['PLATFORM']:
                    iargs += ['--sensor', iDict['PLATFORM'].lower()]
                elif processor == 'cosicorr':
                    iargs += ['--metadata', iDict['mintpy.load.metaFile']]
                print(script_name, ' '.join(iargs))
                # run
                prep_module.main(iargs)

    elif processor == 'isce':
        from mintpy import prep_isce
        from mintpy.utils.isce_utils import get_processor

        # --meta-file
        meta_files = sorted(glob.glob(iDict['mintpy.load.metaFile']))
        if len(meta_files) > 0:
            meta_file = meta_files[0]
        else:
            warnings.warn('No input metadata file found: {}'.format(iDict['mintpy.load.metaFile']))
            meta_file = 'auto'

        # --baseline-dir / --geometry-dir
        baseline_dir = iDict['mintpy.load.baselineDir']
        geom_dir = os.path.dirname(iDict['mintpy.load.demFile'])

        # --dset-dir / --file-pattern
        obs_keys = ['mintpy.load.unwFile', 'mintpy.load.ionUnwFile',
                    'mintpy.load.rgOffFile', 'mintpy.load.azOffFile']
        obs_keys = [i for i in obs_keys if i in iDict['dset_name2template_key'].values()]
        obs_paths = [iDict[key] for key in obs_keys if iDict[key].lower() != 'auto']
        if len(obs_paths) > 0:

            # ifgramStack
            processor = get_processor(meta_file) if os.path.isfile(meta_file) else 'topsStack'
            if processor == 'alosStack':
                obs_dir = os.path.dirname(obs_paths[0])
            else:
                obs_dir = os.path.dirname(os.path.dirname(obs_paths[0]))
            obs_file = os.path.basename(obs_paths[0])

            # ionStack
            if 'mintpy.ionUnwFile' in obs_keys:
                parts = obs_paths[0].rsplit(os.sep, 3)
                obs_dir = parts[0]
                obs_file = os.path.join(parts[-2], parts[-1])

        else:
            obs_dir = None
            obs_file = None

        # --geom-files
        geom_names = ['dem', 'lookupY', 'lookupX', 'incAngle', 'azAngle', 'shadowMask', 'waterMask']
        geom_keys = ['mintpy.load.{}File'.format(i) for i in geom_names]
        geom_files = [os.path.basename(iDict[key]) for key in geom_keys
                      if (iDict[key] and iDict[key].lower() != 'auto')]

        # compose list of input arguments
        iargs = ['-m', meta_file, '-g', geom_dir]
        if baseline_dir:
            iargs += ['-b', baseline_dir]
        if obs_dir is not None:
            iargs += ['-d', obs_dir, '-f', obs_file]
        if geom_files:
            iargs += ['--geom-files'] + geom_files

        # run module
        print('prep_isce.py', ' '.join(iargs))
        try:
            prep_isce.main(iargs)
        except:
            warnings.warn('prep_isce.py failed. Assuming its result exists and continue...')

    elif processor == 'aria':
        from mintpy import prep_aria

        ## compose input arguments
        # use the default template file if exists & input
        default_temp_files = [fname for fname in iDict['template_file']
                              if fname.endswith('smallbaselineApp.cfg')]
        if len(default_temp_files) > 0:
            temp_file = default_temp_files[0]
        else:
            temp_file = iDict['template_file'][0]
        iargs = ['--template', temp_file]

        # file name/dir/path
        ARG2OPT_DICT = {
            '--stack-dir'           : 'mintpy.load.unwFile',
            '--unwrap-stack-name'   : 'mintpy.load.unwFile',
            '--coherence-stack-name': 'mintpy.load.corFile',
            '--conn-comp-stack-name': 'mintpy.load.connCompFile',
            '--dem'                 : 'mintpy.load.demFile',
            '--incidence-angle'     : 'mintpy.load.incAngleFile',
            '--azimuth-angle'       : 'mintpy.load.azAngleFile',
            '--water-mask'          : 'mintpy.load.waterMaskFile',
        }

        for arg_name, opt_name in ARG2OPT_DICT.items():
            arg_value = iDict.get(opt_name, 'auto')
            if arg_value.lower() not in ['auto', 'no', 'none']:
                if arg_name.endswith('dir'):
                    iargs += [arg_name, os.path.dirname(arg_value)]
                elif arg_name.endswith('name'):
                    iargs += [arg_name, os.path.basename(arg_value)]
                else:
                    iargs += [arg_name, arg_value]

        # configurations
        if iDict['compression']:
            iargs += ['--compression', iDict['compression']]
        if iDict['updateMode']:
            iargs += ['--update']

        ## run
        print('prep_aria.py', ' '.join(iargs))
        try:
            prep_aria.main(iargs)
        except:
            warnings.warn('prep_aria.py failed. Assuming its result exists and continue...')

    elif processor == 'gmtsar':
        from mintpy import prep_gmtsar

        # use the custom template file if exists & input
        custom_temp_files = [fname for fname in iDict['template_file']
                             if not fname.endswith('smallbaselineApp.cfg')]
        if len(custom_temp_files) == 0:
            raise FileExistsError('Custom template file NOT found and is required for GMTSAR!')

        # run prep_*.py
        iargs = [custom_temp_files[0], '--mintpy-dir', os.path.dirname(iDict['outdir'])]
        print('prep_gmtsar.py', ' '.join(iargs))
        try:
            prep_gmtsar.main(iargs)
        except:
            warnings.warn('prep_gmtsar.py failed. Assuming its result exists and continue...')

    else:
        msg = 'un-recognized InSAR processor: {}'.format(processor)
        msg += '\nsupported processors: {}'.format(PROCESSOR_LIST)
        raise ValueError(msg)

    return


def print_write_setting(iDict):
    updateMode = iDict['updateMode']
    comp = iDict['compression']
    print('-'*50)
    print('updateMode : {}'.format(updateMode))
    print('compression: {}'.format(comp))
    print('x/ystep: {}/{}'.format(iDict['xstep'], iDict['ystep']))

    # box
    box = iDict['box']
    # box for geometry file in geo-coordinates
    if not iDict.get('geocoded', False):
        boxGeo = iDict['box4geo_lut']
    else:
        boxGeo = box

    return updateMode, comp, box, boxGeo


def get_extra_metadata(iDict):
    """Extra metadata with key names in MACRO_CASE to be written into stack file.

    Parameters: iDict     - dict, input arguments from command lines & template file
                extraDict - dict, extra metadata from template file:
                            E.g. PROJECT_NAME, PLATFORM, ORBIT_DIRECTION, SUBSET_X/YMIN, etc.
    """
    extraDict = {}
    # all keys in MACRO_CASE
    upper_keys = [i for i in iDict.keys() if i.isupper()]
    for key in upper_keys:
        value  = iDict[key]
        if key in ['PROJECT_NAME', 'PLATFORM']:
            if value:
                extraDict[key] = value
        else:
            extraDict[key] = value
    return extraDict


#################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # read input options
    iDict = read_inps2dict(inps)

    # prepare metadata
    prepare_metadata(iDict)
    extraDict = get_extra_metadata(iDict)

    # prepare for subset [need the metadata from above]
    iDict = read_subset_box(iDict)

    # skip data writing for aria as it is included in prep_aria
    if iDict['processor'] == 'aria':
        return

    # initiate objects
    stackObj = read_inps_dict2ifgram_stack_dict_object(iDict)
    geomRadarObj, geomGeoObj = read_inps_dict2geometry_dict_object(iDict)

    # prepare write
    updateMode, comp, box, boxGeo = print_write_setting(iDict)
    if any([stackObj, geomRadarObj, geomGeoObj]) and not os.path.isdir(inps.outdir):
        os.makedirs(inps.outdir)
        print('create directory: {}'.format(inps.outdir))

    # write
    if stackObj and run_or_skip(inps.outfile[0], stackObj, box,
                                updateMode=updateMode,
                                xstep=iDict['xstep'],
                                ystep=iDict['ystep']) == 'run':
        print('-'*50)
        stackObj.write2hdf5(outputFile=inps.outfile[0],
                            access_mode='w',
                            box=box,
                            xstep=iDict['xstep'],
                            ystep=iDict['ystep'],
                            compression=comp,
                            extra_metadata=extraDict)

    if geomRadarObj and run_or_skip(inps.outfile[1], geomRadarObj, box,
                                    updateMode=updateMode,
                                    xstep=iDict['xstep'],
                                    ystep=iDict['ystep']) == 'run':
        print('-'*50)
        geomRadarObj.write2hdf5(outputFile=inps.outfile[1],
                                access_mode='w',
                                box=box,
                                xstep=iDict['xstep'],
                                ystep=iDict['ystep'],
                                compression='lzf',
                                extra_metadata=extraDict)

    if geomGeoObj and run_or_skip(inps.outfile[2], geomGeoObj, boxGeo,
                                  updateMode=updateMode,
                                  xstep=iDict['xstep'],
                                  ystep=iDict['ystep']) == 'run':
        print('-'*50)
        geomGeoObj.write2hdf5(outputFile=inps.outfile[2],
                              access_mode='w',
                              box=boxGeo,
                              xstep=iDict['xstep'],
                              ystep=iDict['ystep'],
                              compression='lzf')

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return


#################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
