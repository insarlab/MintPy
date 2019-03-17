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
from datetime import datetime as dt
import shutil
import argparse
import subprocess
import numpy as np
import pysar
from pysar.objects import sensor, RAMP_LIST
from pysar.utils import readfile, utils as ut
from pysar.defaults.auto_path import autoPath


##########################################################################
STEP_LIST = [
    'load_data',
    'reference_point',
    'stack_interferograms',
    'correct_unwrap_error',
    'modify_network',
    'invert_network',
    'correct_LOD',
    'correct_troposphere',
    'deramp',
    'correct_topography',
    'residual_RMS',
    'reference_date',
    'velocity',
    'geocode',
    'google_earth',
    'hdfeos5',
]

STEP_HELP = """Command line options for steps processing with names are chosen from the following list:

{}
{}
{}

In order to use either --start or --dostep, it is necessary that a
previous run was done using one of the steps options to process at least
through the step immediately preceding the starting step of the current run.
""".format(STEP_LIST[0:5], STEP_LIST[5:10], STEP_LIST[10:])

EXAMPLE = """example:
  pysarApp.py                         #run with default template 'pysarApp_template.txt'
  pysarApp.py <custom_template>       #run with default and custom templates
  pysarApp.py -h / --help             #help
  pysarApp.py -H                      #print    default template options
  pysarApp.py -g                      #generate default template if it does not exist
  pysarApp.py -g <custom_template>    #generate/update default template based on custom template

  # Run with --start/stop/dostep options
  pysarApp.py GalapagosSenDT128.template --dostep velocity  #run at step 'velocity' only
  pysarApp.py GalapagosSenDT128.template --end load_data    #end after step 'load_data'
"""

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, F. Amelung (2019), InSAR time series analysis: error
  correction and noise reduction (under review).
"""

def create_parser():
    parser = argparse.ArgumentParser(description='PySAR Routine Time Series Analysis',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "ignored if the default pysarApp_template.txt is input.")
    parser.add_argument('--dir', dest='workDir',
                        help='specify custom working directory. The default is:\n' +
                             'a) current directory, OR\n' +
                             'b) $SCRATCHDIR/projectName/PYSAR, if:\n' +
                             '    1) autoPath == True in $PYSAR_HOME/pysar/defaults/auto_path.py AND\n' +
                             '    2) environment variable $SCRATCHDIR exists AND\n' +
                             '    3) customTemplateFile is specified (projectName.*)\n')

    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='generate default template (if it does not exist) and exit.')
    parser.add_argument('-H', dest='print_template', action='store_true',
                        help='print the default template file and exit.')
    parser.add_argument('-v','--version', action='store_true', help='print software version and exit')

    parser.add_argument('--noplot', dest='plot', action='store_false',
                        help='do not plot results at the end of the processing.')

    step = parser.add_argument_group('steps processing (start/end/dostep)', STEP_HELP)
    step.add_argument('--start', dest='startStep', metavar='STEP', default=STEP_LIST[0],
                      help='start processing at the named step, default: {}'.format(STEP_LIST[0]))
    step.add_argument('--end', dest='endStep', metavar='STEP',  default=STEP_LIST[-1],
                      help='end processing at the named step, default: {}'.format(STEP_LIST[-1]))
    step.add_argument('--dostep', dest='doStep', metavar='STEP',
                      help='run processing at the named step only')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    template_file = os.path.join(os.path.dirname(__file__), 'defaults/pysarApp_template.txt')

    # print default template
    if inps.print_template:
        raise SystemExit(open(template_file, 'r').read())

    # print software version
    if inps.version:
        raise SystemExit(pysar.version.description)

    if not inps.customTemplateFile and not os.path.isfile(os.path.basename(template_file)):
        parser.print_usage()
        print(EXAMPLE)
        msg = "ERROR: no template file found! It requires:"
        msg += "\n  1) input a custom template file, OR"
        msg += "\n  2) there is a default template 'pysarApp_template.txt' in current directory." 
        print(msg)
        raise SystemExit()

    # invalid input of custom template
    if inps.customTemplateFile:
        inps.customTemplateFile = os.path.abspath(inps.customTemplateFile)
        if not os.path.isfile(inps.customTemplateFile):
            raise FileNotFoundError(inps.customTemplateFile)
        elif os.path.basename(inps.customTemplateFile) == os.path.basename(template_file):
            # ignore if pysarApp_template.txt is input as custom template
            inps.customTemplateFile = None

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
    if idx0 > idx1:
        msg = 'input start step "{}" is AFTER input end step "{}"'.format(inps.startStep, inps.endStep)
        raise ValueError(msg)
    inps.runSteps = STEP_LIST[idx0:idx1+1]

    # empty the step list for -g option
    if inps.generate_template:
        inps.runSteps = []

    # message - software version
    if len(inps.runSteps) <= 1:
        print(pysar.version.description)
    else:
        print(pysar.version.logo)

    # mssage - processing steps
    if len(inps.runSteps) > 0:
        print('Run routine processing with {} on steps: {}'.format(os.path.basename(__file__), inps.runSteps))
        print('Time: {}'.format(dt.now()))

    if inps.doStep:
        print('Remaining steps: {}'.format(STEP_LIST[idx0+1:]))
        print('--dostep option enabled, disable the plotting at the end of the processing.')
        inps.plot = False

    print('-'*50)
    return inps


##########################################################################
class TimeSeriesAnalysis:
    """ Routine processing workflow for time series analysis of small baseline InSAR stacks
    """

    def __init__(self, customTemplateFile=None, workDir=None):
        self.customTemplateFile = customTemplateFile
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
            if not os.path.isdir(inputs_dir):
                os.makedirs(inputs_dir)
                print('create directory:', inputs_dir)
            if ut.run_or_skip(out_file=os.path.join(inputs_dir, os.path.basename(cfile)),
                              in_file=cfile,
                              check_readable=False) == 'run':
                shutil.copy2(cfile, inputs_dir)
                print('copy {} to INPUTS directory for backup.'.format(os.path.basename(cfile)))

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
            for key in ['pysar.save.hdfEos5', 'pysar.save.kmz']:
                if self.template[key] is True:
                    self.template['pysar.geocode'] = True
                    print('Turn ON pysar.geocode in order to run {}.'.format(key))
                    break
        return


    def run_load_data(self, step_name):
        """Load InSAR stacks into HDF5 files in ./INPUTS folder.
        It 1) copy auxiliary files into work directory (for Unvi of Miami only)
           2) load all interferograms stack files into PYSAR/INPUTS directory.
           3) check loading result
           4) add custom metadata (optional, for HDF-EOS5 format only)
        """
        # 1) copy aux files (optional)
        self._copy_aux_file()

        # 2) loading data
        scp_args = '--template {}'.format(self.templateFile)
        if self.customTemplateFile:
            scp_args += ' {}'.format(self.customTemplateFile)
        if self.projectName:
            scp_args += ' --project {}'.format(self.projectName)
        # run
        print("load_data.py", scp_args)
        pysar.load_data.main(scp_args.split())
        os.chdir(self.workDir)

        # 3) check loading result
        load_complete, stack_file, geom_file = ut.check_loaded_dataset(self.workDir, print_msg=True)[0:3]

        # 4) add custom metadata (optional)
        if self.customTemplateFile:
            print('updating {}, {} metadata based on custom template file: {}'.format(
                os.path.basename(stack_file),
                os.path.basename(geom_file),
                os.path.basename(self.customTemplateFile)))
            # use ut.add_attribute() instead of add_attribute.py because of
            # better control of special metadata, such as SUBSET_X/YMIN
            ut.add_attribute(stack_file, self.customTemplate)
            ut.add_attribute(geom_file, self.customTemplate)

        # 5) if not load_complete, plot and raise exception
        if not load_complete:
            # plot result if error occured
            self.plot_result(print_aux=False, plot=plot)

            # go back to original directory
            print('Go back to directory:', self.cwd)
            os.chdir(self.cwd)

            # raise error
            msg = 'step {}: NOT all required dataset found, exit.'.format(step_name)
            raise RuntimeError(msg)
        return


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
                if ut.run_or_skip(out_file=os.path.basename(fname), in_file=fname, check_readable=False) == 'run':
                    shutil.copy2(fname, self.workDir)
                    print('copy {} to work directory'.format(os.path.basename(fname)))
        except:
            pass
        return


    def run_reference_point(self, step_name):
        """Select reference point.
        It 1) generate mask file from common conn comp
           2) generate average spatial coherence and its mask
           3) add REF_X/Y and/or REF_LAT/LON attribute to stack file
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        water_mask_file = 'waterMask.h5'
        cc_mask_file = 'maskConnComp.h5'
        coh_mask_file = 'maskSpatialCoh.h5'
        coh_file = 'avgSpatialCoh.h5'

        # 1) generate mask file from the common connected components
        scp_args = '{} --nonzero -o {} --update'.format(stack_file, cc_mask_file)
        print('generate_mask.py', scp_args)
        pysar.generate_mask.main(scp_args.split())

        # 2.1) generate average spatial coherence
        scp_args = '{} --dataset coherence -o {} --update'.format(stack_file, coh_file)
        print('temporal_average.py', scp_args)
        pysar.temporal_average.main(scp_args.split())

        # 2.2) generate mask based on average spatial coherence
        if ut.run_or_skip(out_file=coh_mask_file, in_file=coh_file) == 'run':
            scp_args = '{} -m 0.7 -o {}'.format(coh_file, coh_mask_file)
            if os.path.isfile(water_mask_file):
                scp_args += ' --base {}'.format(water_mask_file)
            print('generate_mask.py', scp_args)
            pysar.generate_mask.main(scp_args.split())

        # 3) select reference point
        scp_args = '{} -t {} -c {}'.format(stack_file, self.templateFile, coh_file)
        print('reference_point.py', scp_args)
        pysar.reference_point.main(scp_args.split())


    def run_ifgram_stacking(self, step_name):
        """Traditional interferograms stacking."""
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        pha_vel_file = 'avgPhaseVelocity.h5'
        scp_args = '{} --dataset unwrapPhase -o {} --update'.format(stack_file, pha_vel_file)
        print('temporal_average.py', scp_args)
        pysar.temporal_average.main(scp_args.split())


    def run_unwrap_error_correction(self, step_name):
        """Correct phase-unwrapping errors"""
        method = self.template['pysar.unwrapError.method']
        if not method:
            print('phase-unwrapping error correction is OFF.')
            return

        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        mask_file = 'maskConnComp.h5'

        scp_args_bridge = '{} -t {} --update'.format(stack_file, self.templateFile)
        scp_args_closure = '{} {} -t {} --update'.format(stack_file, mask_file, self.templateFile)

        if method == 'bridging':
            pysar.unwrap_error_bridging.main(scp_args_bridge.split())
        elif method == 'phase_closure':
            pysar.unwrap_error_phase_closure.main(scp_args_closure.split())
        elif method == 'bridging+phase_closure':
            scp_args_bridge += ' -i unwrapPhase -o unwrapPhase_bridging'
            pysar.unwrap_error_bridging.main(scp_args_bridge.split())
            scp_args_closure += ' -i unwrapPhase_bridging -o unwrapPhase_bridging_phaseClosure'
            pysar.unwrap_error_phase_closure.main(scp_args_closure.split())
        else:
            raise ValueError('un-recognized method: {}'.format(method))


    def run_network_modification(self, step_name):
        """Modify network of interferograms before the network inversion."""
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        coh_txt = '{}_coherence_spatialAvg.txt'.format(os.path.splitext(os.path.basename(stack_file))[0])
        try:
            net_fig = [i for i in ['Network.pdf', 'PIC/Network.pdf'] if os.path.isfile(i)][0]
        except:
            net_fig = None

        # 1) modify network
        scp_args = '{} -t {}'.format(stack_file, self.templateFile)
        print('modify_network.py', scp_args)
        pysar.modify_network.main(scp_args.split())

        # 2) plot network
        scp_args = '{} -t {} --nodisplay'.format(stack_file, self.templateFile)
        print('\nplot_network.py', scp_args)
        if ut.run_or_skip(out_file=net_fig,
                          in_file=[stack_file, coh_txt, self.templateFile],
                          check_readable=False) == 'run':
            pysar.plot_network.main(scp_args.split())


    def run_network_inversion(self, step_name):
        """Invert network of interferograms for raw phase time-series.
        1) network inversion --> timeseries.h5, temporalCoherence.h5, numInvIfgram.h5
        2) temporalCoherence.h5 --> maskTempCoh.h5
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]

        # 1) invert ifgramStack for time-series
        scp_args = '{} -t {} --update '.format(stack_file, self.templateFile)
        print('ifgram_inversion.py', scp_args)
        pysar.ifgram_inversion.main(scp_args.split())

        # 2) get reliable pixel mask: maskTempCoh.h5
        self.generate_temporal_coherence_mask()


    def generate_temporal_coherence_mask(self):
        """Generate reliable pixel mask from temporal coherence"""
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        tcoh_file = 'temporalCoherence.h5'
        mask_file = 'maskTempCoh.h5'
        tcoh_min = self.template['pysar.networkInversion.minTempCoh']

        scp_args = '{} -m {} -o {} --shadow {}'.format(tcoh_file, tcoh_min, mask_file, geom_file)
        print('generate_mask.py', scp_args)

        # update mode: run only if:
        # 1) output file exists and newer than input file, AND
        # 2) all config keys are the same
        config_keys = ['pysar.networkInversion.minTempCoh']
        print('update mode: ON')
        flag = 'skip'
        if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file, print_msg=False) == 'run':
            flag = 'run'
        else:
            print('1) output file: {} already exists and newer than input file: {}'.format(mask_file, tcoh_file))
            atr = readfile.read_attribute(mask_file)
            if any(str(self.template[i]) != atr.get(i, 'False') for i in config_keys):
                flag = 'run'
                print('2) NOT all key configration parameters are the same: {}'.format(config_keys))
            else:
                print('2) all key configuration parameters are the same: {}'.format(config_keys))
        print('run or skip: {}'.format(flag))

        if flag == 'run':
            pysar.generate_mask.main(scp_args.split())
            # update configKeys
            atr = {}
            for key in config_keys:
                atr[key] = self.template[key]
            ut.add_attribute(mask_file, atr)

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
        atr = readfile.read_attribute(fname0)

        # loop for all steps
        phase_correction_steps = ['correct_LOD', 'correct_troposphere', 'deramp', 'correct_topography']
        for sname in phase_correction_steps:
            # fname0 == fname1 if no valid correction method is set.
            fname0 = fname1
            if sname == 'correct_LOD':
                if atr['PLATFORM'].lower().startswith('env'):
                    fname1 = '{}_LODcor.h5'.format(os.path.splitext(fname0)[0])

            elif sname == 'correct_troposphere':
                method = template['pysar.troposphericDelay.method']
                model  = template['pysar.troposphericDelay.weatherModel']
                if method:
                    if method == 'height_correlation':
                        fname1 = '{}_tropHgt.h5'.format(os.path.splitext(fname0)[0])

                    elif method == 'pyaps':
                        fname1 = '{}_{}.h5'.format(os.path.splitext(fname0)[0], model)

                    else:
                        msg = 'un-recognized tropospheric correction method: {}'.format(method)
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

            elif sname == 'correct_topography':
                method = template['pysar.topographicResidual']
                if method:
                    fname1 = '{}_demErr.h5'.format(os.path.splitext(fname0)[0])

            step = dict()
            step['input'] = fname0
            step['output'] = fname1
            steps[sname] = step

        # step - reference_date
        fnames = [steps[sname]['output'] for sname in phase_correction_steps]
        fnames += [steps[sname]['input'] for sname in phase_correction_steps]
        fnames = sorted(list(set(fnames)))
        step = dict()
        step['input'] = fnames
        steps['reference_date'] = step

        # step - velocity / geocode
        step = dict()
        step['input'] = steps['reference_date']['input'][-1]
        steps['velocity'] = step
        steps['geocode'] = step

        # step - hdfeos5
        if 'Y_FIRST' not in atr.keys():
            step = dict()
            step['input'] = './GEOCODE/geo_{}'.format(steps['reference_date']['input'][-1])
        steps['hdfeos5'] = step
        return steps


    def run_local_oscillator_drift_correction(self, step_name):
        """Correct local oscillator drift (LOD).
        Automatically applied for Envisat data.
        Automatically skipped for all the other data.
        """
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            scp_args = '{} {} -o {}'.format(in_file, geom_file, out_file)
            print('local_oscilator_drift.py', scp_args)
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                pysar.local_oscilator_drift.main(scp_args.split())
        else:
            atr = readfile.read_attribute(in_file)
            sat = atr.get('PLATFORM', None)
            print('No local oscillator drift correction is needed for {}.'.format(sat))


    def run_tropospheric_delay_correction(self, step_name):
        """Correct tropospheric delays."""
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        mask_file = 'maskTempCoh.h5'

        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            poly_order  = self.template['pysar.troposphericDelay.polyOrder']
            tropo_model = self.template['pysar.troposphericDelay.weatherModel']
            weather_dir = self.template['pysar.troposphericDelay.weatherDir']
            method      = self.template['pysar.troposphericDelay.method']

            # Phase/Elevation Ratio (Doin et al., 2009)
            if method == 'height_correlation':
                scp_args = '{f} -g {g} -p {p} -m {m} -o {o}'.format(f=in_file,
                                                                    g=geom_file,
                                                                    p=poly_order,
                                                                    m=mask_file,
                                                                    o=out_file)
                print('tropospheric delay correction with height-correlation approach')
                print('tropo_phase_elevation.py', scp_args)
                if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    pysar.tropo_phase_elevation.main(scp_args.split())

            # Weather Re-analysis Data (Jolivet et al., 2011;2014)
            elif method == 'pyaps':
                scp_args = '-f {f} --model {m} -g {g} -w {w}'.format(f=in_file,
                                                                     m=tropo_model,
                                                                     g=geom_file,
                                                                     w=weather_dir)
                print('Atmospheric correction using Weather Re-analysis dataset (PyAPS, Jolivet et al., 2011)')
                print('Weather Re-analysis dataset:', tropo_model)
                print('tropo_pyaps.py ', scp_args)
                if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    tropo_file = './INPUTS/{}.h5'.format(tropo_model)
                    if os.path.isfile(tropo_file):
                        scp_args = '{f} {t} -o {o} --force'.format(f=in_file, t=tropo_file, o=out_file)
                        print('--------------------------------------------')
                        print('Use existed tropospheric delay file: {}'.format(tropo_file))
                        print('diff.py', scp_args)
                        pysar.diff.main(scp_args.split())
                    else:
                        from pysar import tropo_pyaps
                        tropo_pyaps.main(scp_args.split())
        else:
            print('No tropospheric delay correction.')


    def run_phase_deramping(self, step_name):
        """Estimate and remove phase ramp from each acquisition."""
        mask_file = self.template['pysar.deramp.maskFile']
        method    = self.template['pysar.deramp']

        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            print('Remove for each acquisition a phase ramp: {}'.format(method))
            scp_args = '{f} -s {s} -m {m} -o {o}'.format(f=in_file, s=method, m=mask_file, o=out_file)
            print('remove_ramp.py', scp_args)
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                pysar.remove_ramp.main(scp_args.split())
        else:
            print('No phase ramp removal.')


    def run_topographic_residual_correction(self, step_name):
        """step - correct_topography
        Topographic residual (DEM error) correction (optional).
        """
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            scp_args = '{f} -g {g} -t {t} -o {o} --update '.format(f=in_file,
                                                                   g=geom_file,
                                                                   t=self.templateFile,
                                                                   o=out_file)
            print('dem_error.py', scp_args)
            pysar.dem_error.main(scp_args.split())
        else:
            print('No topographic residual correction.')


    def run_residual_phase_rms(self, step_name):
        """Noise evaluation based on the phase residual."""
        res_file = 'timeseriesResidual.h5'
        if os.path.isfile(res_file):
            scp_args = '{} -t {}'.format(res_file, self.templateFile)
            print('timeseries_rms.py', scp_args)
            pysar.timeseries_rms.main(scp_args.split())
        else:
            print('No residual phase file found! Skip residual RMS analysis.')


    def run_reference_date(self, step_name):
        """Change reference date for all time-series files (optional)."""
        if self.template['pysar.reference.date']:
            in_files = self.get_timeseries_filename(self.template)[step_name]['input']
            scp_args = '-t {} '.format(self.templateFile)
            for in_file in in_files:
                scp_args += ' {}'.format(in_file)
            print('reference_date.py', scp_args)
            pysar.reference_date.main(scp_args.split())
        else:
            print('No reference date change.')


    def run_timeseries2velocity(self, step_name):
        """Estimate average velocity from displacement time-series"""
        ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
        vel_file = 'velocity.h5'
        scp_args = '{f} -t {t} -o {o} --update'.format(f=ts_file,
                                                       t=self.templateFile,
                                                       o=vel_file)
        print('timeseries2velocity.py', scp_args)
        pysar.timeseries2velocity.main(scp_args.split())

        # Velocity from estimated tropospheric delays
        tropo_model = self.template['pysar.troposphericDelay.weatherModel']
        tropo_file = './INPUTS/{}.h5'.format(tropo_model)
        if os.path.isfile(tropo_file):
            suffix = os.path.splitext(os.path.basename(tropo_file))[0].title()
            tropo_vel_file = '{}{}.h5'.format(os.path.splitext(vel_file)[0], suffix)
            scp_args= '{f} -t {t} -o {o} --update'.format(f=tropo_file,
                                                          t=self.templateFile,
                                                          o=tropo_vel_file)
            print('timeseries2velocity.py', scp_args)
            pysar.timeseries2velocity.main(scp_args.split())


    def run_geocode(self, step_name):
        """geocode data files in radar coordinates into ./GEOCODE folder."""
        if self.template['pysar.geocode']:
            ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
            atr = readfile.read_attribute(ts_file)
            if 'Y_FIRST' not in atr.keys():
                # 1. geocode
                out_dir = os.path.join(self.workDir, 'GEOCODE')
                if not os.path.isdir(out_dir):
                    os.makedirs(out_dir)
                    print('create directory:', out_dir)

                geom_file, lookup_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2:4]
                in_files = [geom_file, 'temporalCoherence.h5', ts_file, 'velocity.h5']
                scp_args = '-l {l} -t {t} --outdir {o} --update '.format(l=lookup_file,
                                                                         t=self.templateFile,
                                                                         o=out_dir)
                for in_file in in_files:
                    scp_args += ' {}'.format(in_file)
                print('geocode.py', scp_args)
                pysar.geocode.main(scp_args.split())

                # 2. generate reliable pixel mask in geo coordinate
                geom_file = os.path.join(out_dir, 'geo_{}'.format(os.path.basename(geom_file)))
                tcoh_file = os.path.join(out_dir, 'geo_temporalCoherence.h5')
                mask_file = os.path.join(out_dir, 'geo_maskTempCoh.h5')
                tcoh_min = self.template['pysar.networkInversion.minTempCoh']
                scp_args = '{} -m {} -o {} --shadow {}'.format(tcoh_file, tcoh_min, mask_file, geom_file)
                print('generate_mask.py', scp_args)
                if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file) == 'run':
                    pysar.generate_mask.main(scp_args.split())
        else:
            print('geocoding is OFF')
        return


    def run_save2google_earth(self, step_name):
        """Save velocity file in geo coordinates into Google Earth raster image."""
        if self.template['pysar.save.kmz'] is True:
            print('creating Google Earth KMZ file for geocoded velocity file: ...')
            # input
            vel_file = 'velocity.h5'
            atr = readfile.read_attribute(vel_file)
            if 'Y_FIRST' not in atr.keys():
                vel_file = os.path.join(self.workDir, 'GEOCODE/geo_velocity.h5')

            # output
            kmz_file = '{}.kmz'.format(os.path.splitext(vel_file)[0])
            scp_args = '{} -o {}'.format(vel_file, kmz_file)
            print('save_kmz.py', scp_args)

            # update mode
            try:
                kmz_file = [i for i in [kmz_file, 'PIC/{}'.format(kmz_file)] if os.path.isfile(i)][0]
            except:
                kmz_file = None
            if ut.run_or_skip(out_file=kmz_file, in_file=vel_file, check_readable=False) == 'run':
                pysar.save_kmz.main(scp_args.split())
        else:
            print('save velocity to Google Earth format is OFF.')


    def run_save2hdfeos5(self, step_name):
        """Save displacement time-series and its aux data in geo coordinate into HDF-EOS5 format"""
        if self.template['pysar.save.hdfEos5'] is True:
            # input
            ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
            # Add attributes from custom template to timeseries file
            if self.customTemplate is not None:
                ut.add_attribute(ts_file, self.customTemplate)

            tcoh_file = 'temporalCoherence.h5'
            mask_file = 'geo_maskTempCoh.h5'
            geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
            if 'GEOCODE' in ts_file:
                tcoh_file = './GEOCODE/geo_temporalCoherence.h5'
                mask_file = './GEOCODE/geo_maskTempCoh.h5'
                geom_file = './GEOCODE/geo_{}'.format(os.path.basename(geom_file))

            # cmd
            print('--------------------------------------------')
            scp_args = '{f} -c {c} -m {m} -g {g} -t {t}'.format(f=ts_file,
                                                                c=tcoh_file,
                                                                m=mask_file,
                                                                g=geom_file,
                                                                t=self.templateFile)
            print('save_hdfeos5.py', scp_args)

            # output (check existing file)
            atr = readfile.read_attribute(ts_file)
            SAT = sensor.get_unavco_mission_name(atr)
            try:
                hdfeos5_file = get_file_list('{}_*.he5'.format(SAT))[0]
            except:
                hdfeos5_file = None
            if ut.run_or_skip(out_file=hdfeos5_file, in_file=[ts_file, tcoh_file, mask_file, geom_file]) == 'run':
                pysar.save_hdfeos5.main(scp_args.split())
        else:
            print('save time-series to HDF-EOS5 format is OFF.')


    def plot_result(self, print_aux=True, plot=True):
        """Plot data files and save to figures in PIC folder"""
        if not self.template['pysar.plot'] or not plot:
            return

        def grab_latest_update_date(fname, prefix='# Latest update:'):
            try:
                lines = open(fname, 'r').readlines()
                line = [i for i in lines if prefix in i][0]
                t = re.findall('\d{4}-\d{2}-\d{2}', line)[0]
                t = dt.strptime(t, '%Y-%m-%d')
            except:
                t = dt.strptime('2010-01-01', '%Y-%m-%d') #a arbitrary old date
            return t

        print('\n******************** plot & save to PIC ********************')
        sh_file = os.path.join(os.path.dirname(__file__), '../sh/plot_pysarApp.sh')

        # 1) copy to work directory (if not existed yet).
        if not os.path.isfile(os.path.basename(sh_file)):
            print('copy {} to work directory: {}'.format(sh_file, self.workDir))
            shutil.copy2(sh_file, self.workDir)

        # 2) copy to work directory (if obsolete file detected) and rename the existing one
        elif grab_latest_update_date(os.path.basename(sh_file)) < grab_latest_update_date(sh_file):
            os.system('mv {f} {f}_obsolete'.format(f=os.path.basename(sh_file)))
            print('obsolete shell file detected, renamed it to: {}_obsolete'.format(os.path.basename(sh_file)))
            print('copy {} to work directory: {}'.format(sh_file, self.workDir))
            shutil.copy2(sh_file, self.workDir)

        # cmd
        cmd = './'+os.path.basename(sh_file)
        print(cmd)
        subprocess.Popen(cmd, shell=True).wait()

        # message for more visualization scripts
        msg = """Explore more info & visualization options with the following scripts:
        info.py                    #check HDF5 file structure and metadata
        view.py                    #2D map view
        tsview.py                  #1D point time-series (interactive)   
        transect.py                #1D profile (interactive)
        plot_coherence_matrix.py   #plot coherence matrix for one pixel (interactive)
        plot_network.py            #plot network configuration of the dataset    
        save_kmz.py                #generate Google Earth KMZ file in raster image
        save_kmz_timeseries.py     #generate Goodle Earth KMZ file in points for time-series (interactive)
        """
        if print_aux:
            print(msg)
        return


    def run(self, steps=STEP_LIST, plot=True):
        # run the chosen steps
        for sname in steps:
            print('\n\n******************** step - {} ********************'.format(sname))

            if sname == 'load_data':
                self.run_load_data(sname)

            elif sname == 'reference_point':
                self.run_reference_point(sname)

            elif sname == 'stack_interferograms':
                self.run_ifgram_stacking(sname)

            elif sname == 'correct_unwrap_error':
                self.run_unwrap_error_correction(sname)

            elif sname == 'modify_network':
                self.run_network_modification(sname)

            elif sname == 'invert_network':
                self.run_network_inversion(sname)

            elif sname == 'correct_LOD':
                self.run_local_oscillator_drift_correction(sname)

            elif sname == 'correct_troposphere':
                self.run_tropospheric_delay_correction(sname)

            elif sname == 'deramp':
                self.run_phase_deramping(sname)

            elif sname == 'correct_topography':
                self.run_topographic_residual_correction(sname)

            elif sname == 'residual_RMS':
                self.run_residual_phase_rms(sname)

            elif sname == 'reference_date':
                self.run_reference_date(sname)

            elif sname == 'velocity':
                self.run_timeseries2velocity(sname)

            elif sname == 'geocode':
                self.run_geocode(sname)

            elif sname == 'google_earth':
                self.run_save2google_earth(sname)

            elif sname == 'hdfeos5':
                self.run_save2hdfeos5(sname)

        # plot result (show aux visualization message more multiple steps processing)
        print_aux = len(steps) > 1
        self.plot_result(print_aux=print_aux, plot=plot)

        # go back to original directory
        print('Go back to directory:', self.cwd)
        os.chdir(self.cwd)

        # message
        msg = '\n################################################'
        msg += '\nNormal end of PySAR routine processing workflow!'
        msg += '\n################################################'
        print(msg)
        return


##########################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    app = TimeSeriesAnalysis(inps.customTemplateFile, inps.workDir)
    app.startup()
    if len(inps.runSteps) > 0:
        app.run(steps=inps.runSteps, plot=inps.plot)

    # Timing
    m, s = divmod(time.time()-start_time, 60)
    print('Time used: {:02.0f} mins {:02.1f} secs\n'.format(m, s))
    return


###########################################################################################
if __name__ == '__main__':
    main()
