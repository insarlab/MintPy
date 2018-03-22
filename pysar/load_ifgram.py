#!/usr/bin/env python3
# Author: Heresh Fattahi, Zhang Yunjun

import os, sys, glob
import argparse
import datetime, time

from pysar.utils import readfile
from pysar.objects.ifgramStack import ifgram, ifgramStack

datasetNames = ['unwrapPhase','coherence','connectComponent','wrapPhase','rangeOffset','azimuthOffset']


#################################################################
EXAMPLE='''example:
  load_ifgram.py -t GalapagosSenDT128.tempalte
  load_ifgram.py -p isce -i $SC/GalapagosSenDT128/merged/interferograms/ -f filt_*.unw filt_*.cor filt_*.unw.conncomp
  load_ifgram.py -p roipac -i $SC/GalapagosSenDT128/PROCESS/DONE/ -f filt_*c10.unw filt_*.cor filt_*.int
'''

def createParser():
    '''Create command line parser.'''
    parser = argparse.ArgumentParser(description='Saving a stack of Interferograms to an HDF5 file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('-t','--template', type=str, dest='template_file', help='template file with path info.')
    parser.add_argument('--project', type=str, dest='project_name', help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor','-p', type=str, dest='processor', choices={'isce','roipac','gamma','doris','gmtsar'},\
                        help='InSAR processor/software of the file')

    parser.add_argument('-i','--input-dir', type=str, dest='input_dir', help='directory that include data files')
    parser.add_argument('-f','--file-pattern', type=str, nargs='+', dest='file_pattern',\
                        default=['filt_*.unw','filt_*.cor','filt_*.unw.conncomp','filt_*.int'],\
                        help='name patterns of data files to be loaded in the HDF5 file.')
    parser.add_argument('--dset-name','--file-type', type=str, nargs='+', dest='dset_name', metavar='NAME',\
                        choices=set(datasetNames), help='name of the 3D dataset to be stored in HDF5 file')
    parser.add_argument('-x', type=int, nargs=2, dest='subset_x', metavar=('X_MIN','X_MAX'),\
                        help='Subset range in x/range direction')
    parser.add_argument('-y', type=int, nargs=2, dest='subset_y', metavar=('Y_MIN','Y_MAX'),\
                        help='Subset range in y/azimuth direction')

    parser.add_argument('--enforce', dest='update_mode', action='store_false',\
                        help='Disable the update mode, or skip checking dataset already loaded.')
    parser.add_argument('-o','--output', type=str, dest='outfile', default='ifgramStack.h5', help='output HDF5 file')

    return parser


def cmdLineParse(iargs = None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    if not inps.input_dir and not inps.template_file:
        parser.print_usage()
        print('ERROR: at least one input required: -i/--input-dir or -t/--template')
        sys.exit(1)

    return inps


def read_template2inps(template_file, inps=None):
    if not inps:
        inps = cmdLineParse()
    template = readfile.read_template(template_file)

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


def read_subset_box(inps, length, width):
    if inps.subset_x:
        x0, x1 = sorted(inps.subset_x)
    else:
        x0, x1 = [0, width]
    if inps.subset_y:
        y0, y1 = sorted(inps.subset_y)
    else:
        y0, y1 = [0, length]
    box = (x0, y0, x1, y1)
    length = y1 - y0
    width = x1 - x0
    return box


def write2h5(inps):
    numDataset = len(inps.file_pattern)
    pairsDict = {}

    ifgramDirs = os.path.join(inps.input_dir,'*')
    print('searching interferometric pairs info in directory: {}'.format(ifgramDirs))
    ifgramDirs = glob.glob(ifgramDirs)
    for ifgramDir in ifgramDirs:
        ###### Get master/slave datetime
        file0 = glob.glob(os.path.join(ifgramDir, inps.file_pattern[0]))[0]
        atr = readfile.read_attribute(file0)
        dates = atr['DATE12'].split('-')
        t1 = time.strptime(dates[0],'%y%m%d')
        t2 = time.strptime(dates[1],'%y%m%d')
        Time1 = datetime.datetime(t1.tm_year,t1.tm_mon,t1.tm_mday)
        Time2 = datetime.datetime(t2.tm_year,t2.tm_mon,t2.tm_mday)

        #####################################
        # a dictionary of observartions for a given pair. One pair may 
        # have several types of observations.
        # example obsDict = {'unwrapped phase': /pathToFile/filt.unw, 'iono':/PathToFile/iono.bil}
        datasetDict = {}
        for i in range(numDataset):
            file = glob.glob(os.path.join(ifgramDir, inps.file_pattern[i]))
            if len(file) > 0 and os.path.exists(file[0]):
                datasetName = get_dataset_name(file[0])
                datasetDict[datasetName] = file[0]

        ifgramObj = ifgram(dates=(Time1, Time2), datasetDict=datasetDict, metadata=vars(inps))
        ifgramObj.get_metadata()

        pairsDict[(Time1,Time2)] = ifgramObj

    ############################################
    stackObj = ifgramStack(pairsDict=pairsDict)
    #import pdb; pdb.set_trace()

    outFile = inps.outfile
    box = read_subset_box(inps, ifgramObj.length, ifgramObj.width)
    stackObj.save2h5(outputFile=inps.outfile, access_mode='w', box=box)

    #stackObj.addDatasets('platform-track',inps.fileList, inps.nameList, inps.bandList)
    return inps.outfile


def main(iargs=None):
    inps = cmdLineParse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    if inps.subset_x:  inps.subset_x = sorted(inps.subset_x)
    if inps.subset_y:  inps.subset_y = sorted(inps.subset_y)

    inps.outfile = write2h5(inps)

    return inps.outfile


if __name__ == '__main__' :
    '''
    loading a stack of InSAR pairs to and HDF5 file
    '''
    main()
