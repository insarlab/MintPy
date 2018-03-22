#!/usr/bin/env python3
# Author: Heresh Fattahi, Zhang Yunjun

import os, sys, glob
import argparse
import datetime, time

from pysar.utils import readfile, datetime as ptime
from pysar.utils.insarObj import ifgram, ifgramStack, ifgramDatasetNames


#################################################################
TEMPLATE='''
## 1. Load Data (--load to exit after this step)
#------------------------ ISCE/sentinelStack ---------------------:
pysar.load.processor     = isce
pysar.load.unwFile       = $PROJECT_DIR/merged/interferograms/*/filt_*.unw
pysar.load.corFile       = $PROJECT_DIR/merged/interferograms/*/filt_*.cor
pysar.load.connCompFile  = $PROJECT_DIR/merged/interferograms/*/filt_*.unw.conncomp
pysar.load.intFile       = $PROJECT_DIR/merged/interferograms/*/filt_*.int
#------------------------ ROI_PAC --------------------------------:
pysar.load.processor     = roipac
pysar.load.unwFile       = $PROJECT_DIR/PROCESS/DONE/filt_*rlks_c10.unw
pysar.load.corFile       = $PROJECT_DIR/PROCESS/DONE/filt_*rlks.cor
pysar.load.connCompFile  = $PROJECT_DIR/PROCESS/DONE/filt_*_snap_connect.byt
pysar.load.intFile       = $PROJECT_DIR/PROCESS/DONE/filt_*.int
#------------------------ GAMMA ----------------------------------:
pysar.load.processor     = roipac
pysar.load.unwFile       = $PROJECT_DIR/PROCESS/DONE/diff_*rlks.unw
pysar.load.corFile       = $PROJECT_DIR/PROCESS/DONE/filt_*rlks.cor
pysar.load.connCompFile  = $PROJECT_DIR/PROCESS/DONE/filt_*_snap_connect.byt
pysar.load.intFile       = $PROJECT_DIR/PROCESS/DONE/diff_*.int

pysar.subset.yx   = auto
#pysar.subset.lalo = auto
'''

NOTE='''NOTE:
  unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the master and slave date, either in file name or folder name.
'''

EXAMPLE='''example:
  load_ifgram.py -t GalapagosSenDT128.tempalte
  load_ifgram.py -p "$PROJECT_DIR/merged/interferograms/*/filt_*.unw" "$PROJECT_DIR/merged/interferograms/*/filt_*.cor" --dset-name unwrapPhase coherence
'''

def createParser():
    '''Create command line parser.'''
    parser = argparse.ArgumentParser(description='Saving a stack of Interferograms to an HDF5 file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=NOTE+'\n'+EXAMPLE)

    parser.add_argument('-t','--template', type=str, nargs='*', dest='template_file', help='template file with path info.')
    parser.add_argument('-p','--path-pattern', type=str, nargs='*', dest='path_pattern',\
                        help='path patterns of data files to be loaded in the HDF5 file.\n'+\
                             'one pattern for one type of files, e.g.:\n'+\
                             '"$PROJECT_DIR/merged/interferograms/*/filt_*.unw"\n'+\
                             '"$PROJECT_DIR/merged/interferograms/*/filt_*.cor"\n'+\
                             '"..."\n'+\
                             'Note to quote "" the input path/file pattern!')
    parser.add_argument('--dset-name','--file-type', type=str, nargs='+', dest='dset_name', metavar='NAME',\
                        choices=set(ifgramDatasetNames), help='name of the 3D dataset to be stored in HDF5 file')

    parser.add_argument('--project', type=str, dest='project_name', help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor', type=str, dest='processor', choices={'isce','roipac','gamma','doris','gmtsar'},\
                        help='InSAR processor/software of the file')
    parser.add_argument('-x', type=int, nargs=2, dest='subset_x', metavar=('X_MIN','X_MAX'),\
                        help='Subset range in x/range direction')
    parser.add_argument('-y', type=int, nargs=2, dest='subset_y', metavar=('Y_MIN','Y_MAX'),\
                        help='Subset range in y/azimuth direction')

    parser.add_argument('--enforce', dest='update_mode', action='store_false',\
                        help='Disable the update mode, or skip checking dataset already loaded. [not implemented yet]')
    parser.add_argument('-o','--output', type=str, dest='outfile', default='ifgramStack.h5', help='output HDF5 file')

    return parser


def cmdLineParse(iargs = None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    if not inps.path_pattern and not inps.template_file:
        parser.print_usage()
        print('ERROR: at least one input required: -p/--path-pattern or -t/--template')
        sys.exit(1)

    return inps


#################################################################
def read_template2inps(template_files, inps=None):
    if not inps:
        inps = cmdLineParse()

    template = dict()
    for template_file in template_files:
        template.update(readfile.read_template(template_file))

    # Read template option
    prefix = 'pysar.load.'
    key = prefix+'processor'
    if key in template.keys():
        inps.processor = template[key]

    key = prefix+'ifgramDir'
    if key in template.keys():
        inps.input_dir = template[key]

    key = prefix+'ifgramDir'
    if key in template.keys():
        inps.input_dir = template[key]

    prefix = 'pysar.subset.'
    key = prefix+'yx'
    if key in template.keys():
        sub = [i.strip() for i in template[key].split(',')]
        inps.subset_y = [int(i.strip()) for i in sub[0].split(':')]
        inps.subset_x = [int(i.strip()) for i in sub[1].split(':')]

    return inps


def get_dataset_name(fname):
    '''unwrapPhase, coherence, connectComponent, wrapPhase, rangeOffset, azimuthOffset, ...
    '''
    ext = os.path.splitext(fname)[1]
    if ext in ['.unw']:
        dsetName = 'unwrapPhase'
    elif ext in ['.cor','.cc']:
        dsetName = 'coherence'
    elif ext in ['.conncomp']:
        dsetName = 'connectComponent'
    elif ext in ['.int']:
        dsetName = 'wrapPhase'
    else:
        dsetName = None
    return dsetName


def read_subset_box(inps):
    file0 = glob.glob(inps.path_pattern[0])[0]
    atr = readfile.read_attribute(file0)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    if inps.subset_x:
        x0, x1 = sorted(inps.subset_x)
    else:
        x0, x1 = [0, width]
    if inps.subset_y:
        y0, y1 = sorted(inps.subset_y)
    else:
        y0, y1 = [0, length]
    box = (x0, y0, x1, y1)
    print('box of input  files: {}'.format((0,0,width,length)))
    print('box of data to read: {}'.format(box))
    return box


def read_inps2ifgram_stack_obj(inps):
    '''Read input arguments into dict of ifgramStack object'''
    dsNameNum = len(inps.path_pattern)
    print('searching interferometric pairs info')

    ##Get files number and path for each data type
    print('input data files:')
    dsNameList = []
    dsNumDict = {}
    dsPathDict = {}
    for pathPattern in inps.path_pattern:
        files = glob.glob(pathPattern)
        if len(files)>0:
            dsName = get_dataset_name(files[0])
            dsNameList.append(dsName)
            dsNumDict[dsName] = len(files)
            dsPathDict[dsName] = files
            print('{}: {}'.format(dsName, pathPattern))

    key = 'unwrapPhase'
    if key not in dsNumDict.keys():
        print('ERROR: No reqired {} data files found!'.format(key))
        sys.exit(1)

    dsNumList = [i for i in dsNumDict.values()]
    if any(i != dsNumList[0] for i in dsNumList):
        print('ERROR: Not all types of dataset have the same number of files:')
        for key, value in dsNumDict.items():
            print('number of {}: {}'.format(key,value))
        sys.exit(1)
    print('number of files per type: {}'.format(dsNumList[0]))


    pairsDict = {}
    for dsPath in dsPathDict[key]:
        atr = readfile.read_attribute(dsPath)
        dates = ptime.yyyymmdd(atr['DATE12'].split('-'))

        #####################################
        # A dictionary of data files for a given pair.
        # One pair may have several types of dataset.
        # example ifgramPathDict = {'unwrapPhase': /pathToFile/filt.unw, 'iono':/PathToFile/iono.bil}

        ifgramPathDict = {}
        for i in range(dsNameNum):
            dsName = dsNameList[i]
            dsPath1 = dsPathDict[dsName][i]
            if all(d in dsPath1 for d in dates):
                ifgramPathDict[dsName] = dsPath1
            else:
                dsPath2 = [i for i in dsPathDict[dsName] if all(d in i for d in dates)]
                if len(dsPath2)>0:
                    ifgramPathDict[dsName] = dsPath2[0]
                else:
                    print('WARNING: {} file missing for pair {}'.format(dsName, dates))
        ifgramObj = ifgram(dates=tuple(dates), datasetDict=ifgramPathDict)
        pairsDict[tuple(dates)] = ifgramObj

    stackObj = ifgramStack(pairsDict=pairsDict)
    return stackObj


def write2h5(inps):
    '''Write HDF5 file based on input namespace'''
    stackObj = read_inps2ifgram_stack_obj(inps)
    box = read_subset_box(inps)
    inps.outfile = stackObj.save2h5(outputFile=inps.outfile, access_mode='w', box=box)
    return inps.outfile


#################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    inps.path_pattern = [readfile.check_variable_name(i) for i in inps.path_pattern]
    if inps.subset_x:  inps.subset_x = sorted(inps.subset_x)
    if inps.subset_y:  inps.subset_y = sorted(inps.subset_y)

    inps.outfile = write2h5(inps)
    return inps.outfile


#################################################################
if __name__ == '__main__' :
    '''
    loading a stack of InSAR pairs to and HDF5 file
    '''
    main()
