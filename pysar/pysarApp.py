#!/usr/bin/env python3
############################################################
# Project: PySAR                                           #
# Purpose: InSAR Time Series Analysis in Python            #
# Author: Zhang Yunjun, Heresh Fattahi                     #
# Created: July 2013                                       #
# Copyright (c) 2013-2019, Zhang Yunjun, Heresh Fattahi    #
############################################################


import os
import re
import glob
import time
import shutil
import argparse
import warnings
import subprocess

import numpy as np

from pysar.objects import sensor
from pysar.utils import readfile, utils as ut
from pysar.defaults.auto_path import autoPath
from pysar import version


##########################################################################
STEP_LIST = [
    'loadData',
    'refPoint',
    'stacking',
    'unwCor',
    'netModify',
    'netInversion',
    'tropo',
    'deramp',
    'topo',
    'residRms',
    'refDate',
    'ts2vel',
    'geocode',
    'googleEarth',
    'hdfEos5',
]

EXAMPLE = """example:
  pysarApp.py                       #Run / Rerun
  pysarApp.py <template_file>       #Run / Rerun
  pysarApp.py -h / --help           #Help
  pysarApp.py -H                    #Print all template options

  # Run with --start/stop/dostep
  pysarApp.py GalapagosSenDT128.template --dostep startup   #Do generate default_template from custom_template
  pysarApp.py GalapagosSenDT128.template --stop load_data   #End processing after loading data
"""

def create_parser():
    parser = argparse.ArgumentParser(description='PySAR Routine Time Series Analysis',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "Ignored if the default pysarApp_template.txt is input.")
    parser.add_argument('--dir', dest='workDir',
                        help='Working directory, default:\n' +
                             'a) current directory, OR\n' +
                             'b) $SCRATCHDIR/projectName/PYSAR, if:\n' +
                             '    1) autoPath == True in $PYSAR_HOME/pysar/defaults/auto_path.py AND\n' +
                             '    2) environment variable $SCRATCHDIR exists AND\n' +
                             '    3) customTemplateFile is specified (projectName.*)\n')

    parser.add_argument('-H', dest='print_auto_template', action='store_true',
                        help='Print/Show the example template file for routine processing.')
    parser.add_argument('-v','--version', action='store_true', help='print software version')

    step = parser.add_argument_group('Steps', 'Options for steps processing with start/end/dostep')
    step.add_argument('--start','-s', dest='startStep', default=STEP_LIST[0],
                      help='Start processing at the named step, default: {}'.format(STEP_LIST[0]))
    step.add_argument('--end','-e', dest='endStep',  default=STEP_LIST[-1],
                      help='End processing at the named step, default: {}'.format(STEP_LIST[-1]))
    step.add_argument('--dostep', dest='doStep',
                      help='Run processing only at the named step')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # print full template
    if inps.print_auto_template:
        autoTemplateFile = os.path.join(os.path.dirname(__file__), 'defaults/pysarApp_template.txt')
        print(open(autoTemplateFile, 'r').read())
        raise SystemExit()

    # print software version
    print(version.description)
    if inps.version:
        raise SystemExit()

    # ignore if pysarApp_template.txt is input as custom template
    if (inps.customTemplateFile
            and os.path.basename(inps.customTemplateFile) == 'pysarApp_template.txt'):
        inps.customTemplateFile = None

    print('InSAR time series analysis with {}'.format(os.path.basename(__file__)))
    # check input --start/end/dostep
    for key in ['startStep', 'endStep', 'doStep']:
        value = vars(inps)[key]
        if value and value not in STEP_LIST:
            msg = 'Input step not found: {}'.format(value)
            msg += '\nAvailable steps: {}'.format(STEP_LIST)
            raise ValueError(msg)

    # ignore --start/end input if --dostep is specified
    if inps.doStep:
        inps.startStep = inps.doStep
        inps.endStep = inps.doStep

    # get list of steps to run
    idx0 = STEP_LIST.index(inps.startStep)
    idx1 = STEP_LIST.index(inps.endStep)
    inps.runSteps = STEP_LIST[idx0:idx1]
    print('Run processing on steps: {}'.format(inps.runSteps))
    print('-'*50)
    return inps


##########################################################################
class TimeSeriesAnalysis:
    """ Routine workflow object for InSAR time series analysis
    The routine workflow consists a series of hardwired steps, each step has:
        do${stepName} : bool, to mark whether to run this 
    
    """
    def __init__(self, customTemplateFile=None, runSteps=STEP_LIST, workDir=None):
        self.customTemplateFile = customTemplateFile
        self.runSteps = runSteps
        self.workDir = workDir
        self.cwd = os.path.abspath(os.getcwd())
        return

    def startup(self):
        """The starting point of the workflow. It runs everytime. 
        It 1) grab project name if given
           2) grab and go to work directory
           3) get and read template(s) options
        """
        #1. Get projectName
        self.projectName = None
        if self.customTemplateFile:
            self.projectName = os.path.splitext(os.path.basename(self.customTemplateFile))[0]
            print('Project name:', self.projectName)

        #2. Go to the work directory
        #2.1 Get workDir
        if not self.workDir:
            if autoPath and 'SCRATCHDIR' in os.environ and self.projectName:
                self.workDir = os.path.join(os.getenv('SCRATCHDIR'), self.projectName, 'PYSAR')
            else:
                self.workDir = os.getcwd()
        self.workDir = os.path.abspath(inps.workDir)

        #2.2 Go to workDir
        if not os.path.isdir(self.workDir):
            os.makedirs(self.workDir)
            print('create directory:', self.workDir)
        os.chdir(self.workDir)
        print("Go to work directory:", self.workDir)

        #2.3 Create sub-folders
        for sub_folder in ['GEOCODE','INPUTS','PIC']:
            sub_dir = os.path.join(self.workDir, sub_folder)
            if not os.path.isdir(sub_dir):
                os.makedirs(sub_dir)
                print('create sub-directory:', sub_dir)

        #3. Read templates
        #3.1 Get default template file
        lfile = os.path.join(os.path.dirname(__file__), 'defaults/pysarApp_template.txt')  #latest version
        cfile = os.path.join(self.workDir, 'pysarApp_template.txt')                        #current version
        if not os.path.isfile(cfile):
            print('copy default template file {} to work directory'.format(lfile))
            shutil.copy2(lfile, self.workDir)
        else:
            #cfile is obsolete if any key is missing
            ldict = readfile.read_template(lfile)
            cdict = readfile.read_template(cfile)
            if any([key not in cdict.keys() for key in ldict.keys()]):
                print('obsolete default template detected, update to the latest version.')
                shutil.copy2(lfile, self.workDir)
                #keep the existing option value from obsolete template file
                template_file = ut.update_template_file(cfile, cdict)
            else:
                print('latest template file detected:', cfile)
        self.templateFile = cfile

        # 3.2 read (custom) template files into dicts
        self._read_template()
        return


    def _read_template(self):
        # read custom template, to:
        # 1) update default template
        # 2) add metadata to ifgramStack file and HDF-EOS5 file
        self.customTemplate = None
        if self.customTemplateFile:
            cfile = self.customTemplateFile
            # Copy custom template file to INPUTS directory for backup
            inputs_dir = os.path.join(self.workDir, 'INPUTS')
            if ut.run_or_skip(out_file=os.path.join(inputs_dir, os.path.basename(cfile)),
                              in_file=cfile,
                              check_readable=False) == 'run':
                shutil.copy2(cfile, inputs_dir)
                print('copy {} to INPUTS directory'.format(os.path.basename(cfile)))

            # Read custom template
            print('read custom template file:', cfile)
            cdict = readfile.read_template(cfile)

            # correct some loose type errors
            standardValues = {'def':'auto', 'default':'auto',
                              'y':'yes', 'on':'yes', 'true':'yes',
                              'n':'no', 'off':'no', 'false':'no'
                             }
            for key, value in cdict.items():
                if value in standardValues.keys():
                    cdict[key] = standardValues[value]

            for key in ['pysar.deramp', 'pysar.troposphericDelay.method']:
                if key in cdict.keys():
                    cdict[key] = cdict[key].lower().replace('-', '_')

            if 'processor' in cdict.keys():
                cdict['pysar.load.processor'] = cdict['processor']

            # these metadata are used in load_data.py only, not needed afterwards
            # (in order to manually add extra offset when the lookup table is shifted)
            # (seen in ROI_PAC product sometimes)
            for key in ['SUBSET_XMIN', 'SUBSET_YMIN']:
                if key in cdict.keys():
                    cdict.pop(key)

            self.customTemplate = dict(cdict)

            # Update default template file based on custom template
            print('update default template based on input custom template')
            self.templateFile = ut.update_template_file(self.templateFile, self.customTemplate)

        print('read default template file:', self.templateFile)
        self.template = readfile.read_template(self.templateFile)
        self.template = ut.check_template_auto_value(self.template)    
        return


    def configure(self):
        """
        """
        self.config = argparse.Namespace
        for sname in STEP_list:
            step = argparse.Namespace
            if sname == 'loadData':
                step.cmd = ['load_data.py']
                step.input = []
                step.output = ['INPUTS/ifgramStack.h5']

            elif sname == 'refPoint':
                step.cmd = ['reference_point.py']
                step.input = ['INPUTS/ifgramStack.h5']
                step.output = ['INPUTS/ifgramStack.h5', ]

            # flag - run the step or not
            step.run = True
            if sname not in self.runSteps:
                step.run = False
            vars(self.config)[sname] = step
        return self.config

        
    def run_load_data(self):
        """step - loadData
        It 1) copy auxiliary files into PYSAR/PIC directory (for Unvi of Miami only)
           2) load all interferograms stack files into PYSAR/INPUTS directory.
           3) check loading result
           4) add custom metadata (optional, for HDF-EOS5 format only)
        """
        # 1) copy aux files (optional)
        self._copy_aux_file()

        # 2) loading data
        cmd = 'load_data.py --template {}'.format(self.templateFile)
        if self.customTemplateFile:
            cmd += ' {}'.format(self.customTemplateFile)
        if self.projectName:
            cmd += ' --project {}'.format(self.projectName)
        # run
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        os.chdir(self.workDir)

        # 3) check loading result
        status, stack_file = ut.check_loaded_dataset(self.workDir, print_msg=True)[0:2]

        # 4) add custom metadata (optional)
        if self.customTemplateFile:
            print('updating {} metadata based on custom template file: {}'.format(
                os.path.basename(stack_file), inps.customTemplateFile))
            # use ut.add_attribute() instead of add_attribute.py because of
            # better control of special metadata, such as SUBSET_X/YMIN
            ut.add_attribute(stack_file, customTemplate)
        return status


    def _copy_aux_file(self):
        if not self.projectName:
            return

        # for Univ of Miami
        flist = ['PROCESS/unavco_attributes.txt',
                 'PROCESS/bl_list.txt',
                 'SLC/summary*slc.jpg']
        try:
            proj_dir = os.path.join(os.getenv('SCRATCHDIR'), self.projectName)
            flist = get_file_list([os.path.join(proj_dir, i) for i in flist], abspath=True)
            for fname in flist:
                if run_or_skip(out_file=os.path.basename(fname),
                               in_file=fname,
                               check_readable=False) == 'run':
                    shutil.copy2(fname, self.workDir)
                    print('copy {} to work directory'.format(os.path.basename(fname)))
        except:
            pass
        return


    def run_reference_point(self):
        """step - refPoint
        It 1) generate mask file from common conn comp
           2) generate average spatial coherence
           3) add REF_X/Y and/or REF_LAT/LON attribute to stack file
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        mask_file = 'maskConnComp.h5'
        coh_file = 'avgSpatialCoh.h5'

        # 1) generate mask file from the common connected components
        cmd = 'generate_mask.py {} --nonzero -o {} --update'.format(stack_file, mask_file)
        print(cmd)
        subprocess.Popen(cmd, shell=True).wait()

        # 2) generate average spatial coherence
        cmd = 'temporal_average.py {} --dataset coherence -o {} --update'.format(stack_file, coh_file)
        print(cmd)
        subprocess.Popen(cmd, shell=True).wait()

        # 3) select reference point
        cmd = 'reference_point.py {} -t {} -c {}'.format(stack_file, self.templateFile, coh_file)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_ifgram_stacking(self):
        """step - stacking
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        water_mask_file = 'waterMask.h5'
        pha_vel_file = 'avgPhaseVelocity.h5'
        coh_mask_file = 'maskSpatialCoh.h5'

        # 1) stacking - average phase velocity
        cmd = 'temporal_average.py {} --dataset unwrapPhase -o {} --update'.format(stack_file, pha_vel_file)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()

        # 2) generate mask based on average spatial coherence
        if ut.run_or_skip(out_file=coh_mask_file, in_file=pha_vel_file) == 'run':
            cmd = 'generate_mask.py {} -m 0.7 -o {}'.format(pha_vel_file, coh_mask_file)
            if os.path.isfile(water_mask_file):
                cmd += ' --base {}'.format(water_mask_file)
            print(cmd)
            subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_unwrap_error_correction(self):
        """step - unwCor
        """
        method = self.template['pysar.unwrapError.method']
        if not method:
            return True

        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        mask_file = 'maskConnComp.h5'

        cmd1 = 'unwrap_error_bridging.py {} -t {} --update'.format(stack_file, self.templateFile)
        cmd2 = 'unwrap_error_phase_closure.py {} {} -t {} --update'.format(stack_file, mask_file, self.templateFile)

        if method == 'bridging':
            cmd = cmd1
        elif method == 'phase_closure':
            cmd = cmd2
        elif method == 'bridging+phase_closure':
            cmd = cmd1 + ' -i unwrapPhase -o unwrapPhase_bridging\n'
            cmd += cmd2 + ' -i unwrapPhase_bridging -o unwrapPhase_bridging_phaseClosure'
        else:
            raise ValueError('un-recognized method: {}'.format(method))

        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_network_modification(self):
        """step - netModify
        """

        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        coh_txt = '{}_coherence_spatialAvg.txt'.format(os.path.splitext(os.path.basename(stack_file))[0])
        try:
            net_fig = [i for i in ['Network.pdf', 'PIC/Network.pdf'] if os.path.isfile(i)][0]
        except:
            net_fig = None

        # 1) modify network
        cmd = 'modify_network.py {} -t {}'.format(stack_file, self.templateFile)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()

        # 2) plot network
        cmd = 'plot_network.py {} -t {} --nodisplay'.format(stack_file, self.templateFile)
        print(cmd)
        if ut.run_or_skip(out_file=net_fig,
                          in_file=[stack_file, coh_txt, self.templateFile],
                          check_readable=False) == 'run':
            subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_network_inversion(self):
        """step - netInversion
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]

        # 1) invert ifgramStack for time-series
        cmd = 'ifgram_inversion.py {} -t {} --update '.format(stack_file, self.templateFile)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()

        # 2) 
        return status


    def run():

        self.run_load_data()

        self.run_reference_point()

        self.run_ifgram_stacking()

        self.run_unwrap_error_correction()

        self.run_network_modification()

        self.run_network_inversion()

        self.run_tropospheric_delay_correction()

        self.run_phase_deramping()

        self.run_topographic_residual_correction()

        self.run_residual_phase_rms()

        self.run_reference_date()

        self.run_timeseries2velocity()

        self.run_geocode()

        self.run_save2google_earth()

        self.run_save2hdfeos5()

        # plot results before exit
        self.plot()

        # Go back to original directory
        print('Go to directory:', self.cwd)
        os.chdir(self.cwd)
        return


##########################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    app = TimeSeriesAnalysis(inps.customTemplateFile, inps.runSteps, inps.workDir)
    app.startup()
    app.configure()
    app.run()

    # Timing
    m, s = divmod(time.time()-start_time, 60)
    print('\nTotal time: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return

###########################################################################################
if __name__ == '__main__':
    main()
