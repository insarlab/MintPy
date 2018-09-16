#!/usr/bin/env python3
############################################################
# Project: PySAR                                           #
# Purpose: InSAR Time Series Analysis in Python            #
# Author: Zhang Yunjun, Heresh Fattahi                     #
# Created: July 2013                                       #
# Copyright (c) 2013-2018, Zhang Yunjun, Heresh Fattahi    #
############################################################


import os
import re
import glob
import time
from datetime import datetime as dt
import argparse
import warnings
import shutil
import subprocess
import numpy as np
from pysar.objects import sensor
from pysar.defaults.auto_path import autoPath
from pysar.utils import readfile, utils as ut
from pysar import version


##########################################################################
EXAMPLE = """example:
  pysarApp.py                                             #Run / Rerun
  pysarApp.py  SanAndreasT356EnvD.template  --fast        #Fast processing
  pysarApp.py  SanAndreasT356EnvD.template  --load-data   #Exit after loading data into HDF5 files

  # Template options
  pysarApp.py -H                               #Print    default template
  pysarApp.py -g                               #Generate default template
  pysarApp.py -g SanAndreasT356EnvD.template   #Generate default template considering input custom template
"""


def create_parser():
    parser = argparse.ArgumentParser(description='PySAR Routine Time Series Analysis',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "It's equivalent to None if default pysarApp_template.txt is input.")
    parser.add_argument('--dir', dest='workDir',
                        help='PySAR working directory, default is:\n' +
                             'a) current directory, or\n' +
                             'b) $SCRATCHDIR/projectName/PYSAR, if meets the following 3 requirements:\n' +
                             '    1) autoPath = True in pysar/defaults/auto_path.py\n' +
                             '    2) environmental variable $SCRATCHDIR exists\n' +
                             '    3) input custom template with basename same as projectName\n')
    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='Generate default template (and merge with custom template), then exit.')
    parser.add_argument('-H', dest='print_auto_template', action='store_true',
                        help='Print/Show the example template file for routine processing.')
    parser.add_argument('--version', action='store_true', help='print version number')
    parser.add_argument('--fast', action='store_true',
                        help='Fast processing without using pixel-wised network inversion and DEM error correction')

    parser.add_argument('--reset', action='store_true',
                        help='Reset files attributes to re-run pysarApp.py after loading data by:\n' +
                             '    1) removing ref_y/x/lat/lon for unwrapIfgram.h5 and coherence.h5\n' +
                             '    2) set DROP_IFGRAM=no for unwrapIfgram.h5 and coherence.h5')
    parser.add_argument('--load-data', dest='load_dataset', action='store_true',
                        help='Step 1. Load/check dataset, then exit')
    parser.add_argument('--modify-network', dest='modify_network', action='store_true',
                        help='Step 4. Modify the network, then exit')
    parser.add_argument('--invert-network', dest='invert_network', action='store_true',
                        help='Step 5. Inverse network of interferograms into time-series, then exit')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.autoTemplateFile = os.path.join(os.path.dirname(__file__), '../docs/pysarApp_template.txt')

    if inps.print_auto_template:
        with open(inps.autoTemplateFile, 'r') as f:
            print(f.read())
        raise SystemExit()

    if (inps.customTemplateFile 
            and os.path.basename(inps.customTemplateFile) == 'pysarApp_template.txt'):
        inps.customTemplateFile = None
    return inps


def copy_aux_file(inps):
    # for Univ of Miami
    fileList = ['PROCESS/unavco_attributes.txt',
                'PROCESS/bl_list.txt',
                'SLC/summary*slc.jpg']
    try:
        projectDir = os.path.join(os.getenv('SCRATCHDIR'), inps.projectName)
        fileList = ut.get_file_list([os.path.join(projectDir, i) for i in fileList],
                                    abspath=True)
        for file in fileList:
            if ut.run_or_skip(out_file=os.path.basename(file),
                              in_file=file,
                              check_readable=False) == 'run':
                shutil.copy2(file, inps.workDir)
                print('copy {} to work directory'.format(os.path.basename(file)))
    except:
        pass
    return inps


def check_obsolete_default_template(inps):
    """Update pysarApp_template.txt file if it's obsolete, a.k.a. lack new option names"""
    template_file = os.path.join(inps.workDir, 'pysarApp_template.txt')
    obsolete_template = False
    current_dict = readfile.read_template(template_file)
    latest_dict = readfile.read_template(inps.autoTemplateFile)
    for key in latest_dict.keys():
        if key not in current_dict.keys():
            obsolete_template = True

    if obsolete_template:
        print('obsolete default template detected, update to the latest template options.')
        shutil.copy2(inps.autoTemplateFile, inps.workDir)
        template_file = ut.update_template_file(template_file, current_dict)
    else:
        print('latest template file detected:', template_file)
    return template_file


def read_template(inps):
    print('\n**********  Read Template File  **********')
    # default template
    inps.templateFile = os.path.join(inps.workDir, 'pysarApp_template.txt')
    if not os.path.isfile(inps.templateFile):
        print('generate default template file:', inps.templateFile)
        shutil.copy2(inps.autoTemplateFile, inps.workDir)
    else:
        check_obsolete_default_template(inps)

    # custom template
    customTemplate = None
    if inps.customTemplateFile:
        # Copy custom template file to work directory
        if ut.run_or_skip(out_file=os.path.basename(inps.customTemplateFile),
                          in_file=inps.customTemplateFile,
                          check_readable=False) == 'run':
            shutil.copy2(inps.customTemplateFile, inps.workDir)
            print('copy {} to work directory'.format(os.path.basename(inps.customTemplateFile)))

        # Read custom template
        print('read custom template file:', inps.customTemplateFile)
        customTemplate = readfile.read_template(inps.customTemplateFile)
        # correct some loose type errors
        standardValues = {'def':'auto', 'default':'auto',
                          'y':'yes', 'on':'yes', 'true':'yes',
                          'n':'no', 'off':'no', 'false':'no'
                         }
        for key, value in customTemplate.items():
            if value in standardValues.keys():
                customTemplate[key] = standardValues[value]
        for key in ['pysar.deramp', 'pysar.troposphericDelay.method']:
            if key in customTemplate.keys():
                customTemplate[key] = customTemplate[key].lower().replace('-', '_')
        if 'processor' in customTemplate.keys():
            customTemplate['pysar.load.processor'] = customTemplate['processor']

        # Update default template with custom input template
        print('update default template based on input custom template')
        inps.templateFile = ut.update_template_file(inps.templateFile, customTemplate)

    if inps.generate_template:
        raise SystemExit('Exit as planned after template file generation.')

    print('read default template file:', inps.templateFile)
    template = readfile.read_template(inps.templateFile)
    template = ut.check_template_auto_value(template)

    # Get existing files name: unavco_attributes.txt
    try:
        inps.unavcoMetadataFile = ut.get_file_list('unavco_attribute*txt', abspath=True)[0]
    except:
        inps.unavcoMetadataFile = None
        print('No UNAVCO attributes file found.')

    return inps, template, customTemplate


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
    if ut.run_or_skip(out_file=inps.maskFile, in_file=inps.tempCohFile, print_msg=False) == 'run':
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
        ut.add_attribute(inps.maskFile, meta_dict)

    # check number of pixels selected in mask file for following analysis
    min_num_pixel = float(template['pysar.networkInversion.minNumPixel'])
    msk = readfile.read(inps.maskFile)[0]
    num_pixel = np.sum(msk != 0.)
    print('number of pixels selected: {}'.format(num_pixel))
    if num_pixel < min_num_pixel:
        msg = "Not enought coherent pixels selected (minimum of {}). ".format(int(min_num_pixel))
        msg += "Try the following:\n"
        msg += "1) Check the reference pixel and make sure it's not in areas with unwrapping errors\n"
        msg += "2) Check the network and make sure it's fully connected without subsets"
        raise RuntimeError(msg)
    del msk
    return


def correct_unwrap_error(inps, template):
    unw_cor_method = template['pysar.unwrapError.method']
    if unw_cor_method:
        print('\n**********  Unwrapping Error Correction **********')
        if unw_cor_method == 'phase_closure':
            unwCmd = 'unwrap_error_phase_closure.py {} -t {} --update'.format(inps.stackFile,
                                                                              inps.templateFile)
            if inps.fast:
                unwCmd += ' --fast'

        elif unw_cor_method == 'bridging':
            unwCmd = 'unwrap_error_bridging.py {} -t {} --update'.format(inps.stackFile,
                                                                         inps.templateFile)
        elif unw_cor_method == 'bridging+phase_closure':
            unwCmd = ('unwrap_error_bridging.py {} -t {} --update'
                      ' -i unwrapPhase -o unwrapPhase_bridge').format(inps.stackFile,
                                                                      inps.templateFile)
            unwCmd += ('\nunwrap_error_phase_closure.py {} --update --fast'
                      ' -i unwrapPhase_bridge -o unwrapPhase_bridge_closure').format(inps.stackFile)
        else:
            raise ValueError('un-recognized method: {}'.format(unw_cor_method))

        print(unwCmd)
        status = subprocess.Popen(unwCmd, shell=True).wait()
        if status is not 0:
            raise Exception('Error while correcting phase unwrapping errors.\n')
    return inps


def correct_tropospheric_delay(inps, template):
    """Correct tropospheric delay with options from template"""
    inps.tropPolyOrder = template['pysar.troposphericDelay.polyOrder']
    inps.tropModel     = template['pysar.troposphericDelay.weatherModel']
    inps.tropMethod    = template['pysar.troposphericDelay.method']

    # check existing tropospheric delay file
    try:
        fileList = [os.path.join(inps.workDir, 'INPUTS/{}.h5'.format(inps.tropModel))]
        inps.tropFile = ut.get_file_list(fileList)[0]
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
            if ut.run_or_skip(out_file=outName, in_file=inps.timeseriesFile) == 'run':
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
            if ut.run_or_skip(out_file=outName, in_file=inps.timeseriesFile) == 'run':
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
        inps.tropFile = ut.get_file_list(fileList)[0]
    except:
        inps.tropFile = None
    return


def save_hdfeos5(inps, customTemplate=None):
    if not inps.geocoded:
        warnings.warn('Dataset is in radar coordinates, skip writting to HDF-EOS5 format.')
    else:
        # Add attributes from custom template to timeseries file
        if customTemplate is not None:
            ut.add_attribute(inps.timeseriesFile, customTemplate)

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
            inps.hdfeos5File = ut.get_file_list('{}_*.he5'.format(SAT))[0]
        except:
            inps.hdfeos5File = None
        if ut.run_or_skip(out_file=inps.hdfeos5File, in_file=[inps.timeseriesFile,
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
        
    inps.plotShellFile = os.path.join(os.path.dirname(__file__), '../sh/plot_pysarApp.sh')
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
        msg += '\nPlay with the following scripts for more plotting options:'
        msg += '\nview.py, tsview.py, transect.py, plot_network.py'
        print(msg)
        if status is not 0:
            raise Exception('Error while plotting data files using {}'.format(plotCmd))
    return inps


##########################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)
    if inps.version:
        raise SystemExit(version.version_description)

    #########################################
    # Initiation
    #########################################
    print(version.logo)

    # Project Name
    inps.projectName = None
    if inps.customTemplateFile:
        inps.customTemplateFile = os.path.abspath(inps.customTemplateFile)
        inps.projectName = os.path.splitext(os.path.basename(inps.customTemplateFile))[0]
        print('Project name:', inps.projectName)

    # Work directory
    if not inps.workDir:
        if autoPath and 'SCRATCHDIR' in os.environ and inps.projectName:
            inps.workDir = os.path.join(os.getenv('SCRATCHDIR'), inps.projectName, 'PYSAR')
        else:
            inps.workDir = os.getcwd()
    inps.workDir = os.path.abspath(inps.workDir)
    if not os.path.isdir(inps.workDir):
        os.makedirs(inps.workDir)
    os.chdir(inps.workDir)
    print("Go to work directory:", inps.workDir)

    copy_aux_file(inps)

    inps, template, customTemplate = read_template(inps)

    #########################################
    # Loading Data
    #########################################
    print('\n**********  Load Data  **********')
    loadCmd = 'load_data.py --template {}'.format(inps.templateFile)
    if inps.customTemplateFile:
        loadCmd += ' {}'.format(inps.customTemplateFile)
    if inps.projectName:
        loadCmd += ' --project {}'.format(inps.projectName)
    print(loadCmd)
    status = subprocess.Popen(loadCmd, shell=True).wait()
    os.chdir(inps.workDir)

    print('-'*50)
    inps, atr = ut.check_loaded_dataset(inps.workDir, inps)

    # Add template options into HDF5 file metadata
    if inps.customTemplateFile:
        metaCmd = 'add_attribute.py {} {}'.format(inps.stackFile, inps.customTemplateFile)
        print(metaCmd)
        status = subprocess.Popen(metaCmd, shell=True).wait()
    #ut.add_attribute(inps.stackFile, template)

    if inps.load_dataset:
        raise SystemExit('Exit as planned after loading/checking the dataset.')

    if inps.reset:
        print('Reset dataset attributtes for a fresh re-run.\n'+'-'*50)
        # Reset reference pixel
        refPointCmd = 'reference_point.py {} --reset'.format(inps.stackFile)
        print(refPointCmd)
        status = subprocess.Popen(refPointCmd, shell=True).wait()
        # Reset network modification
        networkCmd = 'modify_network.py {} --reset'.format(inps.stackFile)
        print(networkCmd)
        status = subprocess.Popen(networkCmd, shell=True).wait()

    #########################################
    # Generating Aux files
    #########################################
    print('\n**********  Generate Auxiliary Files  **********')
    inps.waterMaskFile = 'waterMask.h5'
    if not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    # Initial mask (pixels with valid unwrapPhase or connectComponent in ALL interferograms)
    inps.maskFile = 'mask.h5'
    maskCmd = 'generate_mask.py {} --nonzero -o {} --update'.format(inps.stackFile, inps.maskFile)
    print(maskCmd)
    status = subprocess.Popen(maskCmd, shell=True).wait()

    # Average phase velocity - Stacking
    inps.avgPhaseVelFile = 'avgPhaseVelocity.h5'
    avgCmd = 'temporal_average.py {i} --dataset unwrapPhase -o {o} --update'.format(i=inps.stackFile,
                                                                                  o=inps.avgPhaseVelFile)
    print(avgCmd)
    status = subprocess.Popen(avgCmd, shell=True).wait()

    # Average spatial coherence
    inps.avgSpatialCohFile = 'avgSpatialCoherence.h5'
    avgCmd = 'temporal_average.py {i} --dataset coherence -o {o} --update'.format(i=inps.stackFile,
                                                                                  o=inps.avgSpatialCohFile)
    print(avgCmd)
    status = subprocess.Popen(avgCmd, shell=True).wait()

    # mask based on average spatial coherence
    inps.maskSpatialCohFile = 'maskSpatialCoh.h5'
    if ut.run_or_skip(out_file=inps.maskSpatialCohFile, in_file=inps.avgSpatialCohFile) == 'run':
        maskCmd = 'generate_mask.py {i} -m 0.7 -o {o}'.format(i=inps.avgSpatialCohFile,
                                                              o=inps.maskSpatialCohFile)
        if inps.waterMaskFile:
            maskCmd += ' --base {}'.format(inps.waterMaskFile)
        print(maskCmd)
        status = subprocess.Popen(maskCmd, shell=True).wait()


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print('\n**********  Select Reference Point  **********')
    refPointCmd = 'reference_point.py {} -t {} -c {}'.format(inps.stackFile,
                                                             inps.templateFile,
                                                             inps.avgSpatialCohFile)
    print(refPointCmd)
    status = subprocess.Popen(refPointCmd, shell=True).wait()
    if status is not 0:
        raise Exception('Error while finding reference pixel in space.\n')

    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    correct_unwrap_error(inps, template)

    #########################################
    # Network Modification (Optional)
    #########################################
    print('\n**********  Modify Network  **********')
    networkCmd = 'modify_network.py {} -t {}'.format(inps.stackFile,
                                                     inps.templateFile)
    print(networkCmd)
    status = subprocess.Popen(networkCmd, shell=True).wait()
    if status is not 0:
        raise Exception('Error while modifying the network of interferograms.\n')

    # Plot network colored in spatial coherence
    print('--------------------------------------------------')
    plotCmd = 'plot_network.py {} --template {} --nodisplay'.format(inps.stackFile,
                                                                    inps.templateFile)
    print(plotCmd)
    inps.cohSpatialAvgFile = '{}_coherence_spatialAverage.txt'.format(
        os.path.splitext(os.path.basename(inps.stackFile))[0])
    try:
        outFile = [i for i in ['Network.pdf', 'PIC/Network.pdf'] if os.path.isfile(i)][0]
    except:
        outFile = None
    if ut.run_or_skip(out_file=outFile,
                      in_file=[inps.stackFile,
                               inps.cohSpatialAvgFile,
                               inps.templateFile],
                      check_readable=False) == 'run':
        status = subprocess.Popen(plotCmd, shell=True).wait()

    if inps.modify_network:
        raise SystemExit('Exit as planned after network modification.')

    #########################################
    # Inversion of Interferograms
    ########################################
    print('\n**********  Invert Network of Interferograms into Time-series  **********')
    invCmd = 'ifgram_inversion.py {} --template {} --update '.format(inps.stackFile,
                                                                     inps.templateFile)
    if inps.fast:
        invCmd += ' --fast'
    if inps.waterMaskFile:
        invCmd += ' -m {}'.format(inps.waterMaskFile)
    print(invCmd)
    inps.timeseriesFile = 'timeseries.h5'
    inps.tempCohFile = 'temporalCoherence.h5'
    inps.timeseriesFiles = ['timeseries.h5']       #all ts files
    status = subprocess.Popen(invCmd, shell=True).wait()
    if status is not 0:
        raise Exception('Error while inverting network interferograms into timeseries')

    print('\n--------------------------------------------')
    print('Update Mask based on Temporal Coherence ...')
    get_temporal_coherence_mask(inps, template)

    if inps.invert_network:
        raise SystemExit('Exit as planned after network inversion.')

    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   for Envisat data in radar coord only
    ##############################################
    if atr['PLATFORM'].lower().startswith('env'):
        print('\n**********  Local Oscillator Drift Correction for Envisat  **********')
        outName = os.path.splitext(inps.timeseriesFile)[0]+'_LODcor.h5'
        lodCmd = 'local_oscilator_drift.py {} {} -o {}'.format(inps.timeseriesFile,
                                                               inps.geomFile,
                                                               outName)
        print(lodCmd)
        if ut.run_or_skip(out_file=outName, in_file=[inps.timeseriesFile, inps.geomFile]) == 'run':
            status = subprocess.Popen(lodCmd, shell=True).wait()
            if status is not 0:
                raise Exception('Error while correcting Local Oscillator Drift.\n')
        inps.timeseriesFile = outName
        inps.timeseriesFiles.append(outName)

    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print('\n**********  Tropospheric Delay Correction  **********')
    correct_tropospheric_delay(inps, template)

    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print('\n**********  Topographic Residual (DEM error) Correction  **********')
    outName = os.path.splitext(inps.timeseriesFile)[0]+'_demErr.h5'
    topoCmd = 'dem_error.py {} -t {} -o {} --update '.format(inps.timeseriesFile,
                                                             inps.templateFile,
                                                             outName)
    if not inps.fast:
        topoCmd += ' -g {}'.format(inps.geomFile)
    print(topoCmd)
    inps.timeseriesResFile = None
    if template['pysar.topographicResidual']:
        status = subprocess.Popen(topoCmd, shell=True).wait()
        if status is not 0:
            raise Exception('Error while correcting topographic phase residual.\n')
        inps.timeseriesFile = outName
        inps.timeseriesResFile = 'timeseriesResidual.h5'
        inps.timeseriesFiles.append(outName)
    else:
        print('No correction for topographic residuals.')

    # Timeseries Residual Standard Deviation
    print('\n**********  Timeseries Residual Root Mean Square  **********')
    if inps.timeseriesResFile:
        rmsCmd = 'timeseries_rms.py {} -t {}'.format(inps.timeseriesResFile,
                                                     inps.templateFile)
        print(rmsCmd)
        status = subprocess.Popen(rmsCmd, shell=True).wait()
        if status is not 0:
            raise Exception('Error while calculating RMS of time series phase residual.\n')
    else:
        print('No timeseries residual file found! Skip residual RMS analysis.')

    # Reference in Time
    print('\n**********  Select Reference Date  **********')
    if template['pysar.reference.date']:
        refCmd = 'reference_date.py -t {} '.format(inps.templateFile)
        for fname in inps.timeseriesFiles:
            refCmd += ' {}'.format(fname)
        print(refCmd)
        status = subprocess.Popen(refCmd, shell=True).wait()
        if status is not 0:
            raise Exception('Error while changing reference date.\n')
    else:
        print('No reference change in time.')

    ##############################################
    # Phase Ramp Correction (Optional)
    ##############################################
    print('\n**********  Remove Phase Ramp  **********')
    inps.derampMaskFile = template['pysar.deramp.maskFile']
    inps.derampMethod = template['pysar.deramp']
    if inps.derampMethod:
        print('Phase Ramp Removal method: {}'.format(inps.derampMethod))
        ramp_list = ['linear', 'quadratic',
                     'linear_range', 'quadratic_range',
                     'linear_azimuth', 'quadratic_azimuth']
        if inps.derampMethod in ramp_list:
            outName = '{}_ramp.h5'.format(os.path.splitext(inps.timeseriesFile)[0])
            derampCmd = 'remove_ramp.py {} -s {} -m {} -o {}'.format(inps.timeseriesFile,
                                                                     inps.derampMethod,
                                                                     inps.derampMaskFile,
                                                                     outName)
            print(derampCmd)
            if ut.run_or_skip(out_file=outName, in_file=inps.timeseriesFile) == 'run':
                status = subprocess.Popen(derampCmd, shell=True).wait()
                if status is not 0:
                    raise Exception('Error while removing phase ramp for time-series.\n')
            inps.timeseriesFile = outName
            inps.timeseriesFiles.append(outName)
        else:
            msg = 'un-recognized phase ramp method: {}'.format(inps.derampMethod)
            msg += '\navailable ramp types:\n{}'.format(ramp_list)
            raise ValueError(msg)
    else:
        print('No phase ramp removal.')

    #############################################
    # Velocity and rmse maps
    #############################################
    print('\n**********  Estimate Velocity  **********')
    inps.velFile = 'velocity.h5'
    velCmd = 'timeseries2velocity.py {} -t {} -o {} --update'.format(inps.timeseriesFile,
                                                                     inps.templateFile,
                                                                     inps.velFile)
    print(velCmd)
    status = subprocess.Popen(velCmd, shell=True).wait()
    if status is not 0:
        raise Exception('Error while estimating linear velocity from time-series.\n')

    # Velocity from Tropospheric delay
    if inps.tropFile:
        suffix = os.path.splitext(os.path.basename(inps.tropFile))[0].title()
        inps.tropVelFile = '{}{}.h5'.format(os.path.splitext(inps.velFile)[0], suffix)
        velCmd = 'timeseries2velocity.py {} -t {} -o {} --update'.format(inps.tropFile,
                                                                         inps.templateFile,
                                                                         inps.tropVelFile)
        print(velCmd)
        status = subprocess.Popen(velCmd, shell=True).wait()

    ############################################
    # Post-processing
    # Geocodeing --> Masking --> KMZ & HDF-EOS5
    ############################################
    print('\n**********  Post-processing  **********')
    if template['pysar.save.hdfEos5'] is True and template['pysar.geocode'] is False:
        print('Turn ON pysar.geocode to be able to save to HDF-EOS5 format.')
        template['pysar.geocode'] = True

    # Geocoding
    if not inps.geocoded:
        if template['pysar.geocode'] is True:
            print('\n--------------------------------------------')
            geo_dir = os.path.abspath('./GEOCODE')
            if not os.path.isdir(geo_dir):
                os.makedirs(geo_dir)
                print('create directory: {}'.format(geo_dir))
            geoCmd = ('geocode.py {v} {c} {t} {g} -l {l} -t {e}'
                      ' --outdir {d} --update').format(v=inps.velFile,
                                                       c=inps.tempCohFile,
                                                       t=inps.timeseriesFile,
                                                       g=inps.geomFile,
                                                       l=inps.lookupFile,
                                                       e=inps.templateFile,
                                                       d=geo_dir)
            print(geoCmd)
            status = subprocess.Popen(geoCmd, shell=True).wait()
            if status is not 0:
                raise Exception('Error while geocoding.\n')
            else:
                inps.velFile        = os.path.join(geo_dir, 'geo_'+os.path.basename(inps.velFile))
                inps.tempCohFile    = os.path.join(geo_dir, 'geo_'+os.path.basename(inps.tempCohFile))
                inps.timeseriesFile = os.path.join(geo_dir, 'geo_'+os.path.basename(inps.timeseriesFile))
                inps.geomFile       = os.path.join(geo_dir, 'geo_'+os.path.basename(inps.geomFile))
                inps.geocoded = True

            # generate mask based on geocoded temporal coherence
            print('\n--------------------------------------------')
            outName = os.path.join(geo_dir, 'geo_maskTempCoh.h5')
            genCmd = 'generate_mask.py {} -m {} -o {}'.format(inps.tempCohFile,
                                                              inps.minTempCoh,
                                                              outName)
            print(genCmd)
            if ut.run_or_skip(out_file=outName, in_file=inps.tempCohFile) == 'run':
                status = subprocess.Popen(genCmd, shell=True).wait()
            inps.maskFile = outName

    # mask velocity file
    if inps.velFile and inps.maskFile:
        outName = '{}_masked.h5'.format(os.path.splitext(inps.velFile)[0])
        maskCmd = 'mask.py {} -m {} -o {}'.format(inps.velFile,
                                                  inps.maskFile,
                                                  outName)
        print(maskCmd)
        if ut.run_or_skip(out_file=outName, in_file=[inps.velFile, inps.maskFile]) == 'run':
            status = subprocess.Popen(maskCmd, shell=True).wait()
        try:
            inps.velFile = glob.glob(outName)[0]
        except:
            inps.velFile = None

    # Save to Google Earth KML file
    if inps.geocoded and inps.velFile and template['pysar.save.kml'] is True:
        print('\n--------------------------------------------')
        print('creating Google Earth KMZ file for geocoded velocity file: ...')
        outName = '{}.kmz'.format(os.path.splitext(os.path.basename(inps.velFile))[0])
        kmlCmd = 'save_kml.py {} -o {}'.format(inps.velFile, outName)
        print(kmlCmd)
        try:
            outFile = [i for i in [outName, 'PIC/{}'.format(outName)] if os.path.isfile(i)][0]
        except:
            outFile = None
        if ut.run_or_skip(out_file=outFile, in_file=inps.velFile, check_readable=False) == 'run':
            status = subprocess.Popen(kmlCmd, shell=True).wait()
            if status is not 0:
                raise Exception('Error while generating Google Earth KMZ file.')

    #############################################
    # Save Timeseries to HDF-EOS5 format
    #############################################
    if template['pysar.save.hdfEos5'] is True:
        print('\n**********  Save Time-series in HDF-EOS5 Format  **********')
        save_hdfeos5(inps, customTemplate)

    #############################################
    # Plot Figures
    #############################################
    if template['pysar.plot']:
        plot_pysarApp(inps)

    #############################################
    # Timing                                    #
    #############################################
    m, s = divmod(time.time()-start_time, 60)
    print('\n###############################################')
    print('End of PySAR processing!')
    print('################################################\n')
    print('time used: {:02.0f} mins {:02.1f} secs'.format(m, s))


###########################################################################################
if __name__ == '__main__':
    main()
