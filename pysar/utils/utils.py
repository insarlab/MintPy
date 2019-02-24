############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Recommend import:
#   from pysar.utils import utils as ut


import os
import sys
import glob
import time
from datetime import datetime as dt
import shutil
import errno
import subprocess
from argparse import Namespace
import h5py
import numpy as np
from scipy import ndimage, linalg
import matplotlib.pyplot as plt
import multiprocessing

from pysar.objects import (
    deramp,
    geometryDatasetNames,
    geometry,
    ifgramDatasetNames,
    ifgramStack,
    sensor,
    timeseries,
)

from pysar.utils import (
    ptime,
    readfile,
    writefile,
    network as pnet,
)

from pysar.utils.utils0 import *
from pysar.utils.utils1 import *
from pysar.objects.coord import coordinate



#####################################  pysarApp utilities begin ############################################    
def check_loaded_dataset(work_dir='./', print_msg=True):
    """Check the result of loading data for the following two rules:
        1. file existance
        2. file attribute readability

    Parameters: work_dir  : string, PySAR working directory
                print_msg : bool, print out message
    Returns:    True, if all required files and dataset exist; otherwise, ERROR
                    If True, PROCESS, SLC folder could be removed.
                stack_file  : 
                geom_file   :
                lookup_file :
    Example:    work_dir = os.path.expandvars('$SCRATCHDIR/SinabungT495F50AlosA/PYSAR')
                ut.check_loaded_dataset(work_dir)
    """
    load_complete = True

    if not work_dir:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    # 1. interferograms stack file: unwrapPhase, coherence
    flist = [os.path.join(work_dir, 'INPUTS/ifgramStack.h5')]
    stack_file = is_file_exist(flist, abspath=True)
    if stack_file is not None:
        obj = ifgramStack(stack_file)
        obj.open(print_msg=False)
        for dname in ['unwrapPhase', 'coherence']:
            if dname not in obj.datasetNames:
                raise ValueError('required dataset "{}" is missing in file {}'.format(dname, stack_file))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './INPUTS/ifgramStack.h5')

    atr = readfile.read_attribute(stack_file)

    # 2. geom_file: height
    if 'X_FIRST' in atr.keys():
        flist = [os.path.join(work_dir, 'INPUTS/geometryGeo.h5')]
    else:
        flist = [os.path.join(work_dir, 'INPUTS/geometryRadar.h5')]
    geom_file = is_file_exist(flist, abspath=True)
    if geom_file is not None:
        obj = geometry(geom_file)
        obj.open(print_msg=False)
        dname = geometryDatasetNames[0]
        if dname not in obj.datasetNames:
            raise ValueError('required dataset "{}" is missing in file {}'.format(dname, geom_file))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './INPUTS/geometry*.h5')

    # 3. lookup_file: latitude,longitude or rangeCoord,azimuthCoord
    # could be different than geometry file in case of roipac and gamma
    flist = [os.path.join(work_dir, 'INPUTS/geometry*.h5')]
    lookup_file = get_lookup_file(flist, abspath=True, print_msg=print_msg)
    if lookup_file is not None:
        obj = geometry(lookup_file)
        obj.open(print_msg=False)

        if atr['PROCESSOR'] in ['isce', 'doris']:
            dnames = [geometryDatasetNames[1],
                      geometryDatasetNames[2]]
        elif atr['PROCESSOR'] in ['gamma', 'roipac']:
            dnames = [geometryDatasetNames[3],
                      geometryDatasetNames[4]]
        else:
            raise AttributeError('InSAR processor: {}'.format(atr['PROCESSOR']))

        for dname in dnames:
            if dname not in obj.datasetNames:
                load_complete = False
                raise Exception('required dataset "{}" is missing in file {}'.format(dname, lookup_file))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), './INPUTS/geometry*.h5')

    # print message
    if print_msg:
        print(('Loaded dataset are processed by '
               'InSAR software: {}'.format(atr['PROCESSOR'])))
        if 'X_FIRST' in atr.keys():
            print('Loaded dataset is in GEO coordinates')
        else:
            print('Loaded dataset is in RADAR coordinates')
        print('Interferograms Stack: {}'.format(stack_file))
        print('Geometry File       : {}'.format(geom_file))
        print('Lookup Table File   : {}'.format(lookup_file))
        if load_complete:
            print('-'*50)
            print('All data needed found/loaded/copied. Processed 2-pass InSAR data can be removed.')
        print('-'*50)
    return load_complete, stack_file, geom_file, lookup_file


############################################################
def get_temporal_coherence_mask(inps, template):
    """Generate mask from temporal coherence"""
    configKeys = ['pysar.networkInversion.minTempCoh']
    inps.maskFile = 'maskTempCoh.h5'
    inps.minTempCoh = template['pysar.networkInversion.minTempCoh']
    maskCmd = 'generate_mask.py {} -m {} -o {} --shadow {}'.format(inps.tempCohFile,
                                                                   inps.minTempCoh,
                                                                   inps.maskFile,
                                                                   inps.geomFile)
    print(maskCmd)

    # update mode checking
    # run if 1) output file exists; 2) newer than input file and 3) all config keys are the same
    run = False
    if run_or_skip(out_file=inps.maskFile,
                   in_file=inps.tempCohFile,
                   print_msg=False) == 'run':
        run = True
    else:
        print('  1) output file: {} already exists and newer than input file: {}'.format(inps.maskFile,
                                                                                         inps.tempCohFile))
        meta_dict = readfile.read_attribute(inps.maskFile)
        if any(str(template[i]) != meta_dict.get(i, 'False') for i in configKeys):
            run = True
            print('  2) NOT all key configration parameters are the same --> run.\n\t{}'.format(configKeys))
        else:
            print('  2) all key configuration parameters are the same:\n\t{}'.format(configKeys))
    # result
    print('run this step:', run)
    if run:
        status = subprocess.Popen(maskCmd, shell=True).wait()
        if status is not 0:
            raise Exception('Error while generating mask file from temporal coherence.')

        # update configKeys
        meta_dict = {}
        for key in configKeys:
            meta_dict[key] = template[key]
        add_attribute(inps.maskFile, meta_dict)

    # check number of pixels selected in mask file for following analysis
    min_num_pixel = float(template['pysar.networkInversion.minNumPixel'])
    msk = readfile.read(inps.maskFile)[0]
    num_pixel = np.sum(msk != 0.)
    print('number of reliable pixels: {}'.format(num_pixel))
    if num_pixel < min_num_pixel:
        msg = "Not enough reliable pixels (minimum of {}). ".format(int(min_num_pixel))
        msg += "Try the following:\n"
        msg += "1) Check the reference pixel and make sure it's not in areas with unwrapping errors\n"
        msg += "2) Check the network and make sure it's fully connected without subsets"
        raise RuntimeError(msg)
    del msk
    return


def correct_tropospheric_delay(inps, template):
    """Correct tropospheric delay with options from template"""
    inps.tropPolyOrder = template['pysar.troposphericDelay.polyOrder']
    inps.tropModel     = template['pysar.troposphericDelay.weatherModel']
    inps.tropMethod    = template['pysar.troposphericDelay.method']

    # check existing tropospheric delay file
    try:
        fileList = [os.path.join(inps.workDir, 'INPUTS/{}.h5'.format(inps.tropModel))]
        inps.tropFile = get_file_list(fileList)[0]
    except:
        inps.tropFile = None

    # run
    if inps.tropMethod:
        fbase = os.path.splitext(inps.timeseriesFile)[0]

        # Phase/Elevation Ratio (Doin et al., 2009)
        if inps.tropMethod == 'height_correlation':
            outName = '{}_tropHgt.h5'.format(fbase)
            print('tropospheric delay correction with height-correlation approach')
            tropCmd = ('tropcor_phase_elevation.py {t} -g {d} -p {p}'
                       ' -m {m} -o {o}').format(t=inps.timeseriesFile,
                                                d=inps.geomFile,
                                                p=inps.tropPolyOrder,
                                                m=inps.maskFile,
                                                o=outName)
            print(tropCmd)
            if run_or_skip(out_file=outName, in_file=inps.timeseriesFile) == 'run':
                status = subprocess.Popen(tropCmd, shell=True).wait()
                if status is not 0:
                    raise Exception('Error while correcting tropospheric delay.\n')
            inps.timeseriesFile = outName
            inps.timeseriesFiles.append(outName)

        # Weather Re-analysis Data (Jolivet et al., 2011;2014)
        elif inps.tropMethod == 'pyaps':
            inps.weatherDir = template['pysar.troposphericDelay.weatherDir']
            outName = '{}_{}.h5'.format(fbase, inps.tropModel)
            print(('Atmospheric correction using Weather Re-analysis dataset'
                   ' (PyAPS, Jolivet et al., 2011)'))
            print('Weather Re-analysis dataset:', inps.tropModel)
            tropCmd = ('tropcor_pyaps.py -f {t} --model {m} -g {g}'
                       ' -w {w}').format(t=inps.timeseriesFile,
                                         m=inps.tropModel,
                                         g=inps.geomFile,
                                         w=inps.weatherDir)
            print(tropCmd)
            if run_or_skip(out_file=outName, in_file=inps.timeseriesFile) == 'run':
                if inps.tropFile:
                    tropCmd = 'diff.py {} {} -o {} --force'.format(inps.timeseriesFile,
                                                                   inps.tropFile,
                                                                   outName)
                    print('--------------------------------------------')
                    print('Use existed tropospheric delay file: {}'.format(inps.tropFile))
                    print(tropCmd)
                status = subprocess.Popen(tropCmd, shell=True).wait()
                if status is not 0:
                    print('\nError while correcting tropospheric delay, try the following:')
                    print('1) Check the installation of PyAPS')
                    print('   http://earthdef.caltech.edu/projects/pyaps/wiki/Main')
                    print('   Try in command line: python -c "import pyaps"')
                    print('2) Use other tropospheric correction method, height-correlation, for example')
                    print('3) or turn off the option by setting pysar.troposphericDelay.method = no.\n')
                    raise RuntimeError()
            inps.timeseriesFile = outName
            inps.timeseriesFiles.append(outName)
        else:
            print('Un-recognized atmospheric delay correction method: {}'.format(inps.tropMethod))

    # Grab tropospheric delay file
    try:
        fileList = [os.path.join(inps.workDir, 'INPUTS/{}.h5'.format(inps.tropModel))]
        inps.tropFile = get_file_list(fileList)[0]
    except:
        inps.tropFile = None
    return




def save_hdfeos5(inps, customTemplate=None):
    if not inps.geocoded:
        warnings.warn('Dataset is in radar coordinates, skip writting to HDF-EOS5 format.')
    else:
        # Add attributes from custom template to timeseries file
        if customTemplate is not None:
            add_attribute(inps.timeseriesFile, customTemplate)

        # Save to HDF-EOS5 format
        print('--------------------------------------------')
        hdfeos5Cmd = ('save_hdfeos5.py {t} -c {c} -m {m} -g {g}'
                      ' -t {e}').format(t=inps.timeseriesFile,
                                        c=inps.tempCohFile,
                                        m=inps.maskFile,
                                        g=inps.geomFile,
                                        e=inps.templateFile)
        print(hdfeos5Cmd)
        atr = readfile.read_attribute(inps.timeseriesFile)
        SAT = sensor.get_unavco_mission_name(atr)
        try:
            inps.hdfeos5File = get_file_list('{}_*.he5'.format(SAT))[0]
        except:
            inps.hdfeos5File = None
        if run_or_skip(out_file=inps.hdfeos5File, in_file=[inps.timeseriesFile,
                                                           inps.tempCohFile,
                                                           inps.maskFile,
                                                           inps.geomFile]) == 'run':
            status = subprocess.Popen(hdfeos5Cmd, shell=True).wait()
            if status is not 0:
                raise Exception('Error while generating HDF-EOS5 time-series file.\n')
    return



def plot_pysarApp(inps):
    def grab_latest_update_date(fname, prefix='# Latest update:'):
        with open(fname, 'r') as f:
            lines = f.readlines()
        try:
            line = [i for i in lines if prefix in i][0]
            t_update = re.findall('\d{4}-\d{2}-\d{2}', line)[0]
            t_update = dt.strptime(t_update, '%Y-%m-%d')
        except:
            t_update = None
        return t_update
        
    inps.plotShellFile = os.path.join(os.path.dirname(__file__), '../../sh/plot_pysarApp.sh')
    plotCmd = './'+os.path.basename(inps.plotShellFile)
    print('\n**********  Plot Results / Save to PIC  **********')
    # copy to workding directory if not existed yet.
    if not os.path.isfile(plotCmd):
        print('copy {} to work directory: {}'.format(inps.plotShellFile, inps.workDir))
        shutil.copy2(inps.plotShellFile, inps.workDir)
    # rename and copy if obsolete file detected
    else:
        t_exist = grab_latest_update_date(plotCmd)
        t_src = grab_latest_update_date(inps.plotShellFile)
        if not t_exist or t_exist < t_src:
            print('obsolete shell file detected.')
            cmd = 'mv {f} {f}_obsolete'.format(f=os.path.basename(plotCmd))
            print('rename existing file: {}'.format(cmd))
            os.system(cmd)
            print('copy {} to work directory: {}'.format(inps.plotShellFile, inps.workDir))
            shutil.copy2(inps.plotShellFile, inps.workDir)

    if os.path.isfile(plotCmd):
        print(plotCmd)
        status = subprocess.Popen(plotCmd, shell=True).wait()
        msg = '\n'+'-'*50
        msg += '\nUse info.py to check the HDF5 file structure and metadata.'
        msg += '\nUse the following scripts for more visualization options:'
        msg += '\n    view.py                  - 2D map(s) view'
        msg += '\n    tsview.py                - 1D point time-series (interactive)'
        msg += '\n    transect.py              - 1D profile/transection (interactive)'
        msg += '\n    plot_coherence_matrix.py - plot coherence matrix of one pixel (interactive)'
        msg += '\n    plot_network.py          - plot network configuration of the whole dataset'
        print(msg)
        if status is not 0:
            raise Exception('Error while plotting data files using {}'.format(plotCmd))
    return inps
#####################################  pysarApp utilities end ######################################



##################################### Utilities Functions ##########################################
def read_timeseries_lalo(lat, lon, ts_file, lookup_file=None, ref_lat=None, ref_lon=None):
    """ Read time-series of one pixel with input lat/lon
    Parameters: lat/lon     : float, latitude/longitude
                ts_file     : string, filename of time-series HDF5 file
                lookup_file : string, filename of lookup table file
                ref_lat/lon : float, latitude/longitude of reference pixel
    Returns:    dates : 1D np.array of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis   : 1D np.array of float in meter
    """
    # read date
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    dates = ptime.date_list2vector(obj.dateList)[0]
    dates = np.array(dates)

    # read displacement
    coord = coordinate(obj.metadata, lookup_file=lookup_file)
    y, x = coord.geo2radar(lat, lon)[0:2]
    box = (x, y, x+1, y+1)
    dis = readfile.read(ts_file, box=box)[0]
    # reference pixel
    if ref_lat is not None:
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]
    #start at zero
    dis -= dis[0]
    return dates, dis


def read_timeseries_yx(y, x, ts_file, lookup_file=None, ref_y=None, ref_x=None):
    """ Read time-series of one pixel with input y/x
    Parameters: y/x         : int, row/column number of interest
                ts_file     : string, filename of time-series HDF5 file
                lookup_file : string, filename of lookup table file
                ref_y/x     : int, row/column number of reference pixel
    Returns:    dates : 1D np.array of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                dis   : 1D np.array of float in meter
    """
    # read date
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    dates = ptime.date_list2vector(obj.dateList)[0]
    dates = np.array(dates)

    # read displacement
    box = (x, y, x+1, y+1)
    dis = readfile.read(ts_file, box=box)[0]
    # reference pixel
    if ref_y is not None:
        ref_box = (ref_x, ref_y, ref_x+1, ref_y+1)
        dis -= readfile.read(ts_file, box=ref_box)[0]
    #start at zero
    dis -= dis[0]
    return dates, dis
