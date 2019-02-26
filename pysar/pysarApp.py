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

from pysar.objects import sensor, RAMP_LIST
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
        self.workDir = os.path.abspath(self.workDir)

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

        # correct some loose setup conflicts
        if self.template['pysar.geocode'] is False:
            for key in ['pysar.save.hdfEos5', 'pysar.save.kml']:
                if self.template[key] is True:
                    self.template['pysar.geocode'] = True
                    print('Turn ON pysar.geocode in order to run {}.'.format(key))
                    break
        return


    def configure(self):
        """
        """
        self.config = argparse.Namespace
        for sname in STEP_LIST:
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
        1) network inversion --> timeseries.h5, temporalCoherence.h5, numInvIfgram.h5
        2) temporalCoherence.h5 --> maskTempCoh.h5
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]

        # 1) invert ifgramStack for time-series
        cmd = 'ifgram_inversion.py {} -t {} --update '.format(stack_file, self.templateFile)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()

        # 2) get reliable pixel mask: maskTempCoh.h5
        self.generate_temporal_coherence_mask()
        return status


    def generate_temporal_coherence_mask(self):
        """"""
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        tcoh_file = 'temporalCoherence.h5'
        mask_file = 'maskTempCoh.h5'
        tcoh_min = self.template['pysar.networkInversion.minTempCoh']

        cmd = 'generate_mask.py {} -m {} -o {} --shadow {}'.format(tcoh_file, tcoh_min, mask_file, geom_file)
        print(cmd)

        # update mode checking, run if:
        # 1) output file exists and newer than input file, AND
        # 2) all config keys are the same
        config_keys = ['pysar.networkInversion.minTempCoh']
        run = False
        if run_or_skip(out_file=mask_file, in_file=tcoh_file, print_msg=False) == 'run':
            run = True
        else:
            print('  1) output file: {} already exists and newer than input file: {}'.format(mask_file, tcoh_file))
            atr = readfile.read_attribute(mask_file)
            if any(str(self.template[i]) != atr.get(i, 'False') for i in config_keys):
                run = True
                print('  2) NOT all key configration parameters are the same --> run.\n\t{}'.format(config_keys))
            else:
                print('  2) all key configuration parameters are the same:\n\t{}'.format(config_keys))
        print('run this step:', run)

        if run:
            subprocess.Popen(cmd, shell=True).wait()
            # update configKeys
            atr = {}
            for key in config_keys:
                atr[key] = self.template[key]
            add_attribute(mask_file, atr)

        # check number of pixels selected in mask file for following analysis
        num_pixel = np.sum(readfile.read(mask_file)[0] != 0.)
        print('number of reliable pixels: {}'.format(num_pixel))

        min_num_pixel = float(self.template['pysar.networkInversion.minNumPixel'])
        if num_pixel < min_num_pixel:
            msg = "Not enough reliable pixels (minimum of {}). ".format(int(min_num_pixel))
            msg += "Try the following:\n"
            msg += "1) Check the reference pixel and make sure it's not in areas with unwrapping errors\n"
            msg += "2) Check the network and make sure it's fully connected without subsets"
            raise RuntimeError(msg)
        return


    @staticmethod
    def get_timeseries_filename(template):
        """Get input/output time-series filename for each step
        Parameters: template : dict, content of pysarApp_template.txt
        Returns:    steps    : dict of dicts, input/output filenames for each step
        """
        steps = dict()
        fname0 = 'timeseries.h5'
        fname1 = 'timeseries.h5'

        # loop for all steps
        phase_correction_steps = ['LOD', 'tropo', 'deramp', 'topo']
        for sname in phase_correction_steps:
            # fname0 == fname1 if no valid correction method is set.
            fname0 = fname1

            if sname == 'LOD':
                atr = readfile.read_attribute(fname0)
                if atr['PLATFORM'].lower().startswith('env'):
                    fname1 = '{}_LODcor.h5'.format(os.path.splitext(fname0)[0])

            elif sname == 'tropo':
                method = template['pysar.troposphericDelay.method']
                model  = template['pysar.troposphericDelay.weatherModel']
                if method:
                    if method == 'height_correlation':
                        fname1 = '{}_tropHgt.h5'.format(os.path.splitext(fname0)[0])

                    elif method == 'pyaps':
                        fname1 = '{}_{}.h5'.format(os.path.splitext(fname0)[0], model)

                    else:
                        msg = 'Un-recognized tropospheric  crrection method: {}'.format(method)
                        raise ValueError(msg)

            elif sname == 'deramp':
                method = template['pysar.deramp']
                if method:
                    if method in RAMP_LIST:
                        fname1 = '{}_ramp.h5'.format(os.path.splitext(fname0)[0])
                    else:
                        msg = 'un-recognized phase ramp type: {}'.format(method)
                        msg += '\navailable ramp types:\n{}'.format(RAMP_LIST)
                        raise ValueError(msg)

            elif sname == 'topo':
                method = template['pysar.topographicResidual']
                if method:
                    fname1 = '{}_demErr.h5'.format(os.path.splitext(fname0)[0])

            step = dict()
            step['input'] = fname0
            step['output'] = fname1
            steps[sname] = step

        # step - refDate
        fnames = [steps[sname]['output'] for sname in phase_correction_steps]
        fnames += [steps[sname]['input'] for sname in phase_correction_steps]
        fnames = sorted(list(set(fnames)))
        step = dict()
        step['input'] = fnames
        steps['refDate'] = step

        # step - ts2vel / geocode
        step = dict()
        step['input'] = steps['refDate']['input'][-1]
        steps['ts2vel'] = step
        steps['geocode'] = step
        return steps


    def run_local_oscillator_drift_correction(self):
        """Correct local oscillator drift (LOD)
        Automatically applied for Envisat data.
        """
        step = 'LOD'
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        fnames = get_timeseries_filename(self.template)[step]
        in_file = fnames['input']
        out_file = fnames['output']
        status = True
        if in_file != out_file:
            print('\n******************** step - {} (auto for Envisat) ********************'.format(step))
            cmd = 'local_oscilator_drift.py {} {} -o {}'.format(in_file, geom_file, out_file)
            print(cmd)
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                status = subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_tropospheric_delay_correction(self):
        """step - tropo
        """
        step = 'tropo'
        print('\n******************** step - {} ********************'.format(step))
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        mask_file = 'maskTempCoh.h5'

        fnames = get_timeseries_filename(self.template)[step]
        in_file = fnames['input']
        out_file = fnames['output']
        status = True
        if in_file != out_file:
            poly_order  = self.template['pysar.troposphericDelay.polyOrder']
            tropo_model = self.template['pysar.troposphericDelay.weatherModel']
            weather_dir = self.template['pysar.troposphericDelay.weatherDir']
            method      = self.template['pysar.troposphericDelay.method']

            # Phase/Elevation Ratio (Doin et al., 2009)
            if method == 'height_correlation':
                cmd = 'tropcor_phase_elevation.py {f} -g {g} -p {p} -m {m} -o {o}'.format(f=in_file,
                                                                                          g=geom_file,
                                                                                          p=poly_order,
                                                                                          m=mask_file,
                                                                                          o=out_file)
                print('tropospheric delay correction with height-correlation approach')
                print(cmd)
                if run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    status = subprocess.Popen(cmd, shell=True).wait()

            # Weather Re-analysis Data (Jolivet et al., 2011;2014)
            elif method == 'pyaps':
                cmd = 'tropcor_pyaps.py -f {f} --model {m} -g {g} -w {w}'.format(t=in_file,
                                                                                 m=tropo_model,
                                                                                 g=geom_file,
                                                                                 w=weather_dir)
                print('Atmospheric correction using Weather Re-analysis dataset (PyAPS, Jolivet et al., 2011)')
                print('Weather Re-analysis dataset:', tropo_model)
                print(cmd)
                if run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    tropo_file = './INPUTS/{}.h5'.format(tropo_model)
                    if os.path.isfile(tropo_file):
                        cmd = 'diff.py {f} {t} -o {o} --force'.format(f=in_file, t=tropo_file, o=out_file)
                        print('--------------------------------------------')
                        print('Use existed tropospheric delay file: {}'.format(tropo_file))
                        print(cmd)
                    status = subprocess.Popen(cmd, shell=True).wait()
                    if status is not 0:
                        msg = '\nError in step: tropo with {} method'.format(method)
                        msg += '\nTry the following:'
                        msg += '\n1) Check the installation of PyAPS'
                        msg += '\n   http://earthdef.caltech.edu/projects/pyaps/wiki/Main'
                        msg += '\n   Try in command line: python -c "import pyaps"'
                        msg += '\n2) Use other tropospheric correction method, height-correlation, for example'
                        msg += '\n3) or turn off the option by setting pysar.troposphericDelay.method = no.\n'
                        print(msg)
        else:
            print('No tropospheric delay correction.')
        return status


    def run_phase_deramping(self):
        """step - deramp
        """
        step = 'deramp'
        print('\n******************** step - {} ********************'.format(step))
        mask_file = self.template['pysar.deramp.maskFile']
        method    = self.template['pysar.deramp']

        fnames = get_timeseries_filename(self.template)[step]
        in_file = fnames['input']
        out_file = fnames['output']
        status = True
        if in_file != out_file:
            print('Remove for each acquisition a phase ramp: {}'.format(method))
            cmd = 'remove_ramp.py {f} -s {s} -m {m} -o {o}'.format(f=in_file, s=method, m=mask_file, o=out_file)
            print(cmd)
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                status = subprocess.Popen(cmd, shell=True).wait()
        else:
            print('No phase ramp removal.')
        return status


    def run_topographic_residual_correction(self):
        """step - topo
        """
        step = 'topo'
        print('\n******************** step - {} ********************'.format(step))
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        fnames = get_timeseries_filename(self.template)[step]
        in_file = fnames['input']
        out_file = fnames['output']
        status = True
        if in_file != out_file:
            cmd = 'dem_error.py {f} -g {g} -t {t} -o {o} --update '.format(f=in_file,
                                                                           g=geom_file,
                                                                           t=self.templateFile,
                                                                           o=out_file)
            print(cmd)
            status = subprocess.Popen(cmd, shell=True).wait()
        else:
            print('No topographic residual correction.')
        return status


    def run_residual_phase_rms(self):
        """step - residRms"""
        step = 'residRms'
        print('\n******************** step - {} ********************'.format(step))
        res_file = 'timeseriesResidual.h5'
        status = True
        if os.path.isfile(res_file):
            cmd = 'timeseries_rms.py {} -t {}'.format(res_file, self.templateFile)
            print(cmd)
            status = subprocess.Popen(cmd, shell=True).wait()
        else:
            print('No residual phase file found! Skip residual RMS analysis.')
        return status


    def run_reference_date(self):
        """step - refDate"""
        status = True
        step = 'refDate'
        print('\n******************** step - {} ********************'.format(step))
        if self.template['pysar.reference.date']:
            in_files = get_timeseries_filename(self.template)[step]['input']
            cmd = 'reference_date.py -t {} '.format(self.templateFile)
            for in_file in in_files:
                cmd += ' {}'.format(in_file)
            print(cmd)
            status = subprocess.Popen(cmd, shell=True).wait()
        else:
            print('No reference date change.')
        return status


    def run_timeseries2velocity(self):
        """step - ts2vel"""
        status = True
        step = 'ts2vel'
        print('\n******************** step - {} ********************'.format(step))
        ts_file = get_timeseries_filename(self.template)[step]['input']
        vel_file = 'velocity.h5'
        cmd = 'timeseries2velocity.py {f} -t {t} -o {o} --update'.format(f=ts_file,
                                                                         t=self.templateFile,
                                                                         o=vel_file)
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()

        # Velocity from estimated tropospheric delays
        tropo_file = './INPUTS/{}.h5'.format(tropo_model)
        if os.path.isfile(tropo_file):
            suffix = os.path.splitext(os.path.basename(tropo_file))[0].title()
            tropo_vel_file = '{}{}.h5'.format(os.path.splitext(vel_file)[0], suffix)
            cmd = 'timeseries2velocity.py {f} -t {t} -o {o} --update'.format(f=tropo_file,
                                                                             t=self.templateFile,
                                                                             o=tropo_vel_file)
            print(cmd)
            subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_geocode(self):
        """step - geocode"""
        status = True
        step = 'geocode'
        print('\n******************** step - {} ********************'.format(step))
        if self.template['pysar.geocode']:
            ts_file = get_timeseries_filename(self.template)[step]['input']
            atr = readfile.read_attribute(ts_file)
            if 'Y_FIRST' not in atr.keys():
                # 1. geocode
                geom_file, lookup_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2:4]
                in_files = [geom_file, 'temporalCoherence.h5', ts_file, 'velocity.h5']
                out_dir = os.path.join(self.workDir, 'GEOCODE')
                geocode_script = os.path.join(os.path.dirname(__file__), 'geocode.py')

                cmd = '{scp} -l {l} -t {t} --outdir {o} --update '.format(l=lookup_file, t=self.templateFile, o=out_dir)
                for in_file in in_files:
                    cmd += ' {}'.format(in_file)
                print(cmd)
                status = subprocess.Popen(cmd, shell=True).wait()

                # 2. generate reliable pixel mask in geo coordinate
                geom_file = os.path.join(out_dir, 'geo_{}'.format(os.path.basename(geom_file)))
                tcoh_file = os.path.join(out_dir, 'geo_temporalCoherence.h5')
                mask_file = os.path.join(out_dir, 'geo_maskTempCoh.h5')
                tcoh_min = self.template['pysar.networkInversion.minTempCoh']
                cmd = 'generate_mask.py {} -m {} -o {} --shadow {}'.format(tcoh_file, tcoh_min, mask_file, geom_file)
                print(cmd)
                if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file) == 'run':
                    subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_save2google_earth(self):
        """step - googleEarth
        """
        status = True
        step = 'googleEarth'
        if self.template['pysar.save.kml'] is True:
            print('\n******************** step - {} ********************'.format(step))
            print('creating Google Earth KMZ file for geocoded velocity file: ...')
            # input
            vel_file = 'velocity.h5'
            atr = readfile.read_attribute(vel_file)
            if 'Y_FIRST' not in atr.keys():
                vel_file = os.path.join(self.workDir, 'GEOCODE/geo_velocity.h5')

            # output
            kmz_file = '{}.kmz'.format(os.path.splitext(os.path.basename(vel_file))[0])
            try:
                kmz_file = [i for i in [kmz_file, 'PIC/{}'.format(kmz_file)] if os.path.isfile(i)][0]
            except:
                kmz_file = None

            cmd = 'save_kmz.py {} -o {}'.format(vel_file, kmz_file)
            print(cmd)
            if ut.run_or_skip(out_file=kmz_file, in_file=vel_file, check_readable=False) == 'run':
                status = subprocess.Popen(cmd, shell=True).wait()
        return status


    def run_save2hdfeos5(self):
        """step - hdfeos5"""
        status = True
        return status


    def run():

        self.run_load_data()

        self.run_reference_point()

        self.run_ifgram_stacking()

        self.run_unwrap_error_correction()

        self.run_network_modification()

        self.run_network_inversion()

        self.run_local_oscillator_drift_correction()

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
