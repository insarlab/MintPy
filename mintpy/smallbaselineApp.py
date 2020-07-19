#!/usr/bin/env python3
############################################################
# Project: MintPy                                          #
# Purpose: Miami InSAR Time-series software in Python      #
# Author: Zhang Yunjun, Heresh Fattahi                     #
# Created: July 2013                                       #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
############################################################


import os
import re
import sys
import time
import datetime
import shutil
import argparse
import subprocess
import numpy as np

import mintpy
from mintpy.objects import sensor, RAMP_LIST
from mintpy.utils import readfile, writefile, utils as ut
from mintpy.defaults.template import STEP_LIST
# dynamic import for modules used by smallbaselineApp workflow
import mintpy.workflow


##########################################################################
STEP_HELP = """Command line options for steps processing with names are chosen from the following list:

{}
{}
{}

In order to use either --start or --dostep, it is necessary that a
previous run was done using one of the steps options to process at least
through the step immediately preceding the starting step of the current run.
""".format(STEP_LIST[0:5], STEP_LIST[5:10], STEP_LIST[10:])

EXAMPLE = """example:
  smallbaselineApp.py                         #run with default template 'smallbaselineApp.cfg'
  smallbaselineApp.py <custom_template>       #run with default and custom templates
  smallbaselineApp.py -h / --help             #help
  smallbaselineApp.py -H                      #print    default template options
  smallbaselineApp.py -g                      #generate default template if it does not exist
  smallbaselineApp.py -g <custom_template>    #generate/update default template based on custom template

  # Run with --start/stop/dostep options
  smallbaselineApp.py GalapagosSenDT128.template --dostep velocity  #run at step 'velocity' only
  smallbaselineApp.py GalapagosSenDT128.template --end load_data    #end after step 'load_data'
"""

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis: 
  Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
  doi:10.1016/j.cageo.2019.104331.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Routine Time Series Analysis for Small Baseline InSAR Stack',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "ignored if the default smallbaselineApp.cfg is input.")
    parser.add_argument('--dir', '--work-dir', dest='workDir', default='./',
                        help='work directory, (default: %(default)s).')

    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='generate default template (if it does not exist) and exit.')
    parser.add_argument('-H', dest='print_template', action='store_true',
                        help='print the default template file and exit.')
    parser.add_argument('-v','--version', action='store_true', help='print software version and exit')

    parser.add_argument('--noplot', dest='plot', action='store_false',
                        help='do not plot results at the end of the processing.')

    step = parser.add_argument_group('steps processing (start/end/dostep)', STEP_HELP)
    step.add_argument('--start', dest='startStep', metavar='STEP', default=STEP_LIST[0],
                      help='start processing at the named step (default: %(default)s).')
    step.add_argument('--end','--stop', dest='endStep', metavar='STEP',  default=STEP_LIST[-1],
                      help='end processing at the named step (default: %(default)s)')
    step.add_argument('--dostep', dest='doStep', metavar='STEP',
                      help='run processing at the named step only')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    template_file = os.path.join(os.path.dirname(__file__), 'defaults/smallbaselineApp.cfg')

    # -H (print default template)
    if inps.print_template:
        print(open(template_file, 'r').read())
        sys.exit(0)

    # -v (print software version)
    if inps.version:
        print(mintpy.version.description)
        sys.exit(0)

    # check all input template files
    if (not inps.customTemplateFile
            and not os.path.isfile(os.path.basename(template_file))
            and not inps.generate_template):
        parser.print_usage()
        print(EXAMPLE)
        msg = "no template file found! It requires:"
        msg += "\n  1) input a custom template file, OR"
        msg += "\n  2) there is a default template 'smallbaselineApp.cfg' in current directory." 
        raise SystemExit('ERROR: '+msg)

    # check custom input template file
    if inps.customTemplateFile:
        if not os.path.isfile(inps.customTemplateFile):
            raise FileNotFoundError(inps.customTemplateFile)

        inps.customTemplateFile = os.path.abspath(inps.customTemplateFile)
        # ignore if smallbaselineApp.cfg is input as custom template
        if os.path.basename(inps.customTemplateFile) == os.path.basename(template_file):
            inps.customTemplateFile = None

    # check input --start/end/dostep
    inps.runSteps = read_inps2run_steps(inps)
    # skip plotting if --dostep is specified
    if inps.doStep:
        inps.plot = False

    return inps


def read_inps2run_steps(inps, step_list=STEP_LIST):
    """read/get run_steps from """
    # check inputs
    for key in ['startStep', 'endStep', 'doStep']:
        value = vars(inps)[key]
        if value and value not in step_list:
            msg = 'Input step not found: {}'.format(value)
            msg += '\nAvailable steps: {}'.format(step_list)
            raise ValueError(msg)

    # ignore --start/end input if --dostep is specified
    if inps.doStep:
        inps.startStep = inps.doStep
        inps.endStep = inps.doStep

    # get list of steps to run
    idx0 = step_list.index(inps.startStep)
    idx1 = step_list.index(inps.endStep)
    if idx0 > idx1:
        msg = 'input start step "{}" is AFTER input end step "{}"'.format(inps.startStep, inps.endStep)
        raise ValueError(msg)
    run_steps = step_list[idx0:idx1+1]

    # empty the step list for -g option
    if inps.generate_template:
        run_steps = []

    # mssage - processing steps
    if len(run_steps) > 0:
        # for single step - compact version info
        if len(run_steps) == 1:
            print(mintpy.version.description)
        else:
            print(mintpy.version.logo)

        print('--RUN-at-{}--'.format(datetime.datetime.now()))
        print('Current directory: {}'.format(os.getcwd()))
        print('Run routine processing with {} on steps: {}'.format(os.path.basename(__file__), run_steps))
        print('Remaining steps: {}'.format(step_list[idx0+1:]))
    print('-'*50)
    return run_steps


def get_the_latest_default_template_file(work_dir):
    """Get the latest version of default template file.
    If an obsolete file exists in the working directory, the existing option values are kept.
    """
    lfile = os.path.join(os.path.dirname(__file__), 'defaults/smallbaselineApp.cfg')  #latest version
    cfile = os.path.join(work_dir, 'smallbaselineApp.cfg')                            #current version
    if not os.path.isfile(cfile):
        print('copy default template file {} to work directory'.format(lfile))
        shutil.copy2(lfile, work_dir)
    else:
        #cfile is obsolete if any key is missing
        ldict = readfile.read_template(lfile)
        cdict = readfile.read_template(cfile)
        if any([key not in cdict.keys() for key in ldict.keys()]):
            print('obsolete default template detected, update to the latest version.')
            shutil.copy2(lfile, work_dir)

            #keep the existing option value from obsolete template file
            ut.update_template_file(cfile, cdict)
    return cfile


##########################################################################
class TimeSeriesAnalysis:
    """ Routine processing workflow for time series analysis of small baseline InSAR stacks
    """

    def __init__(self, customTemplateFile=None, workDir=None):
        self.customTemplateFile = customTemplateFile
        self.workDir = os.path.abspath(workDir)
        self.cwd = os.path.abspath(os.getcwd())
        return

    def startup(self):
        """The starting point of the workflow. It runs everytime. 
        It 1) grab project name if given
           2) go to work directory
           3) get and read template(s) options
           4) get plot shell script to work directory
        """

        #1. Get projectName
        self.projectName = None
        if self.customTemplateFile:
            self.projectName = os.path.splitext(os.path.basename(self.customTemplateFile))[0]
            print('Project name:', self.projectName)

        #2. Go to work directory
        os.makedirs(self.workDir, exist_ok=True)
        os.chdir(self.workDir)
        print('Go to work directory:', self.workDir)

        #3. Read templates
        self.templateFile = get_the_latest_default_template_file(self.workDir)
        self._read_template()

        # 4. Copy the plot shell file
        sh_file = os.path.join(os.path.dirname(__file__), '../sh/plot_smallbaselineApp.sh')

        def grab_latest_update_date(fname, prefix='# Latest update:'):
            try:
                lines = open(fname, 'r').readlines()
                line = [i for i in lines if prefix in i][0]
                t = re.findall('\d{4}-\d{2}-\d{2}', line)[0]
                t = datetime.datetime.strptime(t, '%Y-%m-%d')
            except:
                t = datetime.datetime.strptime('2010-01-01', '%Y-%m-%d') #a arbitrary old date
            return t

        # 1) copy to work directory (if not existed yet).
        if not os.path.isfile(os.path.basename(sh_file)):
            print('copy {} to work directory: {}'.format(sh_file, self.workDir))
            shutil.copy2(sh_file, self.workDir)

        # 2) copy to work directory (if obsolete file detected) and rename the existing one
        elif grab_latest_update_date(os.path.basename(sh_file)) < grab_latest_update_date(sh_file):
            shutil.move(os.path.basename(sh_file), os.path.basename(sh_file)+'_obsolete')
            print('obsolete shell file detected, renamed it to: {}_obsolete'.format(os.path.basename(sh_file)))
            print('copy {} to work directory: {}'.format(sh_file, self.workDir))
            shutil.copy2(sh_file, self.workDir)

        self.plot_sh_cmd = './'+os.path.basename(sh_file)
        return


    def _read_template(self):
        # 1) update default template
        self.customTemplate = None
        if self.customTemplateFile:
            # customTemplateFile --> customTemplate
            print('read custom template file:', self.customTemplateFile)
            cdict = readfile.read_template(self.customTemplateFile)

            # correct some loose type errors
            standardValues = {'def':'auto', 'default':'auto',
                              'y':'yes', 'on':'yes', 'true':'yes',
                              'n':'no', 'off':'no', 'false':'no'
                             }
            for key, value in cdict.items():
                if value in standardValues.keys():
                    cdict[key] = standardValues[value]

            for key in ['mintpy.deramp', 'mintpy.troposphericDelay.method']:
                if key in cdict.keys():
                    cdict[key] = cdict[key].lower().replace('-', '_')

            if 'processor' in cdict.keys():
                cdict['mintpy.load.processor'] = cdict['processor']

            # these metadata are used in load_data.py only, not needed afterwards
            # (in order to manually add extra offset when the lookup table is shifted)
            # (seen in ROI_PAC product sometimes)
            for key in ['SUBSET_XMIN', 'SUBSET_YMIN']:
                if key in cdict.keys():
                    cdict.pop(key)

            self.customTemplate = dict(cdict)

            # customTemplate --> templateFile
            print('update default template based on input custom template')
            self.templateFile = ut.update_template_file(self.templateFile, self.customTemplate)

        # 2) backup custome/default template file in inputs/pic folder
        for backup_dirname in ['inputs', 'pic']:
            backup_dir = os.path.join(self.workDir, backup_dirname)
            # create directory
            os.makedirs(backup_dir, exist_ok=True)

            # back up to the directory
            for tfile in [self.customTemplateFile, self.templateFile]:
                if tfile and  ut.run_or_skip(out_file=os.path.join(backup_dir, os.path.basename(tfile)),
                                             in_file=tfile,
                                             check_readable=False,
                                             print_msg=False) == 'run':
                    shutil.copy2(tfile, backup_dir)
                    print('copy {} to {} directory for backup.'.format(os.path.basename(tfile),
                                                                       os.path.basename(backup_dir)))

        # 3) read default template file
        print('read default template file:', self.templateFile)
        self.template = readfile.read_template(self.templateFile)
        self.template = ut.check_template_auto_value(self.template)

        # correct some loose setup conflicts
        if self.template['mintpy.geocode'] is False:
            for key in ['mintpy.save.hdfEos5', 'mintpy.save.kmz']:
                if self.template[key] is True:
                    self.template['mintpy.geocode'] = True
                    print('Turn ON mintpy.geocode in order to run {}.'.format(key))
                    break
        return


    def run_load_data(self, step_name):
        """Load InSAR stacks into HDF5 files in ./inputs folder.
        It 1) copy auxiliary files into work directory (for Unvi of Miami only)
           2) load all interferograms stack files into mintpy/inputs directory.
           3) check loading result
           4) add custom metadata (optional, for HDF-EOS5 format only)
        """
        # 1) copy aux files (optional)
        self._copy_aux_file()

        # 2) loading data
        stack_processor = self.template['mintpy.load.processor'].lower()
        if stack_processor == 'aria':
            from mintpy import prep_aria
            iargs = ['--template', self.templateFile, '--update']
            prep_aria.main(iargs)

        else:
            # compose list of input arguments
            # instead of using command line then split
            # to support path with whitespace
            iargs = ['--template', self.templateFile]
            if self.customTemplateFile:
                iargs += [self.customTemplateFile]
            if self.projectName:
                iargs += ['--project', self.projectName]

            # run command line
            print('load_data.py', ' '.join(iargs))
            mintpy.load_data.main(iargs)

        # come back to working directory
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
            flist = ut.get_file_list([os.path.join(proj_dir, i) for i in flist], abspath=True)
            for fname in flist:
                if ut.run_or_skip(out_file=os.path.basename(fname), in_file=fname, check_readable=False) == 'run':
                    shutil.copy2(fname, self.workDir)
                    print('copy {} to work directory'.format(os.path.basename(fname)))
        except:
            pass
        return


    def run_network_modification(self, step_name, plot=True):
        """Modify network of interferograms before the network inversion."""
        # check the existence of ifgramStack.h5
        stack_file, geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1:3]
        coh_txt = '{}_coherence_spatialAvg.txt'.format(os.path.splitext(os.path.basename(stack_file))[0])
        try:
            net_fig = [i for i in ['Network.pdf', 'pic/Network.pdf'] if os.path.isfile(i)][0]
        except:
            net_fig = None

        # 1) output waterMask.h5 to simplify the detection/use of waterMask
        water_mask_file = 'waterMask.h5'
        if 'waterMask' in readfile.get_dataset_list(geom_file):
            print('generate {} from {} for conveniency'.format(water_mask_file, geom_file))
            if ut.run_or_skip(out_file=water_mask_file, in_file=geom_file) == 'run':
                water_mask, atr = readfile.read(geom_file, datasetName='waterMask')
                # ignore no-data pixels in geometry files
                ds_name_list = readfile.get_dataset_list(geom_file)
                for ds_name in ['latitude','longitude']:
                    if ds_name in ds_name_list:
                        print('set pixels with 0 in {} to 0 in waterMask'.format(ds_name))
                        ds = readfile.read(geom_file, datasetName=ds_name)[0]
                        water_mask[ds == 0] = 0
                atr['FILE_TYPE'] = 'waterMask'
                writefile.write(water_mask, out_file=water_mask_file, metadata=atr)

        # 2) modify network
        iargs = [stack_file, '-t', self.templateFile]
        print('modify_network.py', ' '.join(iargs))
        mintpy.modify_network.main(iargs)

        # 3) plot network
        if self.template['mintpy.plot'] and plot:
            iargs = [stack_file, '-t', self.templateFile, '--nodisplay']

            dsNames = readfile.get_dataset_list(stack_file)
            if any('phase' in i.lower() for i in dsNames):
                iargs += ['-d', 'coherence', '-v', '0.2', '1.0']
            elif any('offset' in i.lower() for i in dsNames):
                iargs += ['-d', 'offsetSNR', '-v', '0', '20']

            print('\nplot_network.py', ' '.join(iargs))
            if ut.run_or_skip(out_file=net_fig,
                              in_file=[stack_file, coh_txt, self.templateFile],
                              check_readable=False) == 'run':
                mintpy.plot_network.main(iargs)
        return


    def generate_ifgram_aux_file(self):
        """Generate auxiliary files from ifgramStack file"""
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        dsNames = readfile.get_dataset_list(stack_file)
        mask_file = 'maskConnComp.h5'
        coh_file = 'avgSpatialCoh.h5'
        snr_file = 'avgSpatialSnr.h5'

        # 1) generate mask file from the common connected components
        if any('phase' in i.lower() for i in dsNames):
            iargs = [stack_file, '--nonzero', '-o', mask_file, '--update']
            print('\ngenerate_mask.py', ' '.join(iargs))
            mintpy.generate_mask.main(iargs)

        # 2) generate average spatial coherence
        if any('phase' in i.lower() for i in dsNames):
            iargs = [stack_file, '--dataset', 'coherence', '-o', coh_file, '--update']
        elif any('offset' in i.lower() for i in dsNames):
            iargs = [stack_file, '--dataset', 'offsetSNR', '-o', snr_file, '--update']
        print('\ntemporal_average.py', ' '.join(iargs))
        mintpy.temporal_average.main(iargs)
        return


    def run_reference_point(self, step_name):
        """Select reference point.
        It 1) generate mask file from common conn comp
           2) generate average spatial coherence and its mask
           3) add REF_X/Y and/or REF_LAT/LON attribute to stack file
        """
        # 1-2) aux files: maskConnComp and avgSpatialCoh
        self.generate_ifgram_aux_file()

        # 3) add REF_X/Y(/LAT/LON) of the reference point
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        coh_file = 'avgSpatialCoh.h5'

        iargs = [stack_file, '-t', self.templateFile, '-c', coh_file]
        print('reference_point.py', ' '.join(iargs))
        mintpy.reference_point.main(iargs)
        return


    def run_quick_overview(self, step_name):
        """A quick overview on the interferogram stack for:
            1) avgPhaseVelocity.h5: possible ground deformation through interferogram stacking
            2) numNonzeroIntClosure.h5: phase unwrapping errors through the integer ambiguity of phase closure
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]

        # 1) stack interferograms
        pha_vel_file = 'avgPhaseVelocity.h5'
        iargs = [stack_file, '--dataset', 'unwrapPhase', '-o', pha_vel_file, '--update']
        print('temporal_average.py', ' '.join(iargs))
        mintpy.temporal_average.main(iargs)

        # 2) calculate the integer ambiguity of closure phase
        water_mask_file = 'waterMask.h5'
        iargs = [stack_file, '--water-mask', water_mask_file, '--action', 'calculate', '--update']
        print('unwrap_error_phase_closure.py', ' '.join(iargs))
        mintpy.unwrap_error_phase_closure.main(iargs)
        return


    def run_unwrap_error_correction(self, step_name):
        """Correct phase-unwrapping errors"""
        method = self.template['mintpy.unwrapError.method']
        if not method:
            print('phase-unwrapping error correction is OFF.')
            return

        # check required input files
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        mask_file = 'maskConnComp.h5'

        iargs_bridge = [stack_file, '--template', self.templateFile, '--update']
        iargs_closure = iargs_bridge + ['--cc-mask', mask_file]

        if method == 'bridging':
            print('unwrap_error_bridging.py', ' '.join(iargs_bridge))
            mintpy.unwrap_error_bridging.main(iargs_bridge)

        elif method == 'phase_closure':
            print('unwrap_error_phase_closure.py', ' '.join(iargs_closure))
            mintpy.unwrap_error_phase_closure.main(iargs_closure)

        elif method == 'bridging+phase_closure':
            iargs_bridge += ['-i', 'unwrapPhase',
                             '-o', 'unwrapPhase_bridging']
            print('unwrap_error_bridging.py', ' '.join(iargs_bridge))
            mintpy.unwrap_error_bridging.main(iargs_bridge)

            iargs_closure += ['-i', 'unwrapPhase_bridging',
                              '-o', 'unwrapPhase_bridging_phaseClosure']
            print('unwrap_error_phase_closure.py', ' '.join(iargs_closure))
            mintpy.unwrap_error_phase_closure.main(iargs_closure)

        else:
            raise ValueError('un-recognized method: {}'.format(method))
        return


    def run_network_inversion(self, step_name):
        """Invert network of interferograms for raw phase time-series.
        1) network inversion --> timeseries.h5, temporalCoherence.h5, numInvIfgram.h5
        2) temporalCoherence.h5 --> maskTempCoh.h5
        """
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]

        # 1) invert ifgramStack for time-series
        iargs = [stack_file, '-t', self.templateFile, '--update']
        print('ifgram_inversion.py', ' '.join(iargs))
        mintpy.ifgram_inversion.main(iargs)

        # 2) get reliable pixel mask: maskTempCoh.h5
        self.generate_temporal_coherence_mask()
        return


    def generate_temporal_coherence_mask(self):
        """Generate reliable pixel mask from temporal coherence"""
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        tcoh_file = 'temporalCoherence.h5'
        mask_file = 'maskTempCoh.h5'
        tcoh_min = self.template['mintpy.networkInversion.minTempCoh']

        # compose list of arguments
        iargs = [tcoh_file, '-m', tcoh_min, '-o', mask_file]
        # exclude pixels in shadow if shadowMask dataset is available
        if (self.template['mintpy.networkInversion.shadowMask'] is True 
                and 'shadowMask' in readfile.get_dataset_list(geom_file)):
            iargs += ['--base', geom_file, '--base-dataset', 'shadowMask', '--base-value', '1']
        print('generate_mask.py', ' '.join(iargs))

        # update mode: run only if:
        # 1) output file exists and newer than input file, AND
        # 2) all config keys are the same
        config_keys = ['mintpy.networkInversion.{}'.format(i) for i in ['minTempCoh','shadowMask']]
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
            mintpy.generate_mask.main(iargs)
            # update configKeys
            atr = {}
            for key in config_keys:
                atr[key] = self.template[key]
            ut.add_attribute(mask_file, atr)

        # check number of pixels selected in mask file for following analysis
        num_pixel = np.sum(readfile.read(mask_file)[0] != 0.)
        print('number of reliable pixels: {}'.format(num_pixel))

        min_num_pixel = float(self.template['mintpy.networkInversion.minNumPixel'])
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
        Parameters: template : dict, content of smallbaselineApp.cfg
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
                method = template['mintpy.troposphericDelay.method']
                model  = template['mintpy.troposphericDelay.weatherModel']
                if method:
                    if method == 'height_correlation':
                        fname1 = '{}_tropHgt.h5'.format(os.path.splitext(fname0)[0])

                    elif method == 'pyaps':
                        fname1 = '{}_{}.h5'.format(os.path.splitext(fname0)[0], model)

                    else:
                        msg = 'un-recognized tropospheric correction method: {}'.format(method)
                        raise ValueError(msg)

            elif sname == 'deramp':
                method = template['mintpy.deramp']
                if method:
                    if method in RAMP_LIST:
                        fname1 = '{}_ramp.h5'.format(os.path.splitext(fname0)[0])
                    else:
                        msg = 'un-recognized phase ramp type: {}'.format(method)
                        msg += '\navailable ramp types:\n{}'.format(RAMP_LIST)
                        raise ValueError(msg)

            elif sname == 'correct_topography':
                method = template['mintpy.topographicResidual']
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
            step['input'] = './geo/geo_{}'.format(steps['reference_date']['input'][-1])
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
            iargs = [in_file, geom_file, '-o', out_file]
            print('local_oscilator_drift.py', ' '.join(iargs))
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                mintpy.local_oscilator_drift.main(iargs)
        else:
            atr = readfile.read_attribute(in_file)
            sat = atr.get('PLATFORM', None)
            print('No local oscillator drift correction is needed for {}.'.format(sat))
        return



    def run_tropospheric_delay_correction(self, step_name):
        """Correct tropospheric delays."""
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        mask_file = 'maskTempCoh.h5'

        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            poly_order  = self.template['mintpy.troposphericDelay.polyOrder']
            tropo_model = self.template['mintpy.troposphericDelay.weatherModel']
            weather_dir = self.template['mintpy.troposphericDelay.weatherDir']
            method      = self.template['mintpy.troposphericDelay.method']

            def get_dataset_size(fname):
                atr = readfile.read_attribute(fname)
                return (atr['LENGTH'], atr['WIDTH'])

            # Phase/Elevation Ratio (Doin et al., 2009)
            if method == 'height_correlation':
                tropo_look = self.template['mintpy.troposphericDelay.looks']
                tropo_min_cor = self.template['mintpy.troposphericDelay.minCorrelation']
                iargs = [in_file, 
                         '-g', geom_file,
                         '-p', poly_order,
                         '-m', mask_file,
                         '-o', out_file,
                         '-l', tropo_look, 
                         '-t', tropo_min_cor]
                print('tropospheric delay correction with height-correlation approach')
                print('tropo_phase_elevation.py', ' '.join(iargs))
                if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    mintpy.tropo_phase_elevation.main(iargs)

            # Weather Re-analysis Data (Jolivet et al., 2011;2014)
            elif method == 'pyaps':
                iargs = ['-f', in_file, '--model', tropo_model, '-g', geom_file, '-w', weather_dir]
                print('Atmospheric correction using Weather Re-analysis dataset (PyAPS, Jolivet et al., 2011)')
                print('Weather Re-analysis dataset:', tropo_model)
                tropo_file = './inputs/{}.h5'.format(tropo_model)
                if ut.run_or_skip(out_file=out_file, in_file=[in_file, tropo_file]) == 'run':
                    if os.path.isfile(tropo_file) and get_dataset_size(tropo_file) == get_dataset_size(in_file):
                        iargs = [in_file, tropo_file, '-o', out_file, '--force']
                        print('--------------------------------------------')
                        print('Use existed tropospheric delay file: {}'.format(tropo_file))
                        print('diff.py', ' '.join(iargs))
                        mintpy.diff.main(iargs)
                    else:
                        if tropo_model in ['ERA5']:
                            from mintpy import tropo_pyaps3
                            print('tropo_pyaps3.py', ' '.join(iargs))
                            tropo_pyaps3.main(iargs)
                        else:
                            # opt 1 - using tropo_pyaps as python module and call its main function
                            # prefered, disabled for now to make it compatible with python2-pyaps
                            #print('tropo_pyaps.py', ' '.join(iargs))
                            #from mintpy import tropo_pyaps
                            #tropo_pyaps.main(iargs)
                            # opt 2 - using tropo_pyaps as executable script
                            # will be deprecated after python3-pyaps is fully funcational
                            cmd = 'tropo_pyaps.py '+' '.join(iargs)
                            print(cmd)
                            status = subprocess.Popen(cmd, shell=True).wait()

        else:
            print('No tropospheric delay correction.')
        return


    def run_phase_deramping(self, step_name):
        """Estimate and remove phase ramp from each acquisition."""
        mask_file = self.template['mintpy.deramp.maskFile']
        method    = self.template['mintpy.deramp']

        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            print('Remove for each acquisition a phase ramp: {}'.format(method))
            iargs = [in_file, '-s', method, '-m', mask_file, '-o', out_file, '--update']
            print('remove_ramp.py', ' '.join(iargs))
            mintpy.remove_ramp.main(iargs)
        else:
            print('No phase ramp removal.')
        return


    def run_topographic_residual_correction(self, step_name):
        """step - correct_topography
        Topographic residual (DEM error) correction (optional).
        """
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
        fnames = self.get_timeseries_filename(self.template)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']

        if in_file != out_file:
            iargs = [in_file, '-t', self.templateFile, '-o', out_file, '--update']
            if self.template['mintpy.topographicResidual.pixelwiseGeometry']:
                iargs += ['-g', geom_file]
            print('dem_error.py', ' '.join(iargs))
            mintpy.dem_error.main(iargs)

        else:
            print('No topographic residual correction.')
        return


    def run_residual_phase_rms(self, step_name):
        """Noise evaluation based on the phase residual."""
        res_file = 'timeseriesResidual.h5'
        if os.path.isfile(res_file):
            iargs = [res_file, '-t', self.templateFile]
            print('timeseries_rms.py', ' '.join(iargs))
            mintpy.timeseries_rms.main(iargs)
        else:
            print('No residual phase file found! Skip residual RMS analysis.')
        return


    def run_reference_date(self, step_name):
        """Change reference date for all time-series files (optional)."""
        if self.template['mintpy.reference.date']:
            iargs = ['-t', self.templateFile]
            in_files = self.get_timeseries_filename(self.template)[step_name]['input']
            for in_file in in_files:
                iargs += [in_file]
            print('reference_date.py', ' '.join(iargs))
            mintpy.reference_date.main(iargs)
        else:
            print('No reference date change.')
        return


    def run_timeseries2velocity(self, step_name):
        """Estimate average velocity from displacement time-series"""
        ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
        vel_file = 'velocity.h5'

        iargs = [ts_file, '-t', self.templateFile, '-o', vel_file, '--update']
        print('timeseries2velocity.py', ' '.join(iargs))
        mintpy.timeseries2velocity.main(iargs)

        # Velocity from estimated tropospheric delays
        tropo_model = self.template['mintpy.troposphericDelay.weatherModel']
        tropo_file = './inputs/{}.h5'.format(tropo_model)
        if os.path.isfile(tropo_file):
            suffix = os.path.splitext(os.path.basename(tropo_file))[0]
            tropo_vel_file = '{}{}.h5'.format(os.path.splitext(vel_file)[0], suffix)

            iargs = [tropo_file, '-t', self.templateFile, '-o', tropo_vel_file, '--update']
            print('timeseries2velocity.py', ' '.join(iargs))
            mintpy.timeseries2velocity.main(iargs)
        return


    def run_geocode(self, step_name):
        """geocode data files in radar coordinates into ./geo folder."""
        if self.template['mintpy.geocode']:
            ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
            atr = readfile.read_attribute(ts_file)
            if 'Y_FIRST' not in atr.keys():
                # 1. geocode
                out_dir = os.path.join(self.workDir, 'geo')
                os.makedirs(out_dir, exist_ok=True)

                geom_file, lookup_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2:4]
                in_files = [geom_file, 'temporalCoherence.h5', ts_file, 'velocity.h5']
                iargs = ['-l', lookup_file, '-t', self.templateFile, '--outdir', out_dir, '--update']
                for in_file in in_files:
                    iargs += [in_file]
                print('geocode.py', ' '.join(iargs))
                mintpy.geocode.main(iargs)

                # 2. generate reliable pixel mask in geo coordinate
                geom_file = os.path.join(out_dir, 'geo_{}'.format(os.path.basename(geom_file)))
                tcoh_file = os.path.join(out_dir, 'geo_temporalCoherence.h5')
                mask_file = os.path.join(out_dir, 'geo_maskTempCoh.h5')
                tcoh_min = self.template['mintpy.networkInversion.minTempCoh']

                iargs = [tcoh_file, '-m', tcoh_min, '-o', mask_file]
                # exclude pixels in shadow if shadowMask dataset is available
                if (self.template['mintpy.networkInversion.shadowMask'] is True 
                        and 'shadowMask' in readfile.get_dataset_list(geom_file)):
                    iargs += ['--base', geom_file, '--base-dataset', 'shadowMask', '--base-value', '1']
                print('generate_mask.py', ' '.join(iargs))

                if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file) == 'run':
                    mintpy.generate_mask.main(iargs)
            else:
                print('dataset is geocoded, skip geocoding and continue.')
        else:
            print('geocoding is OFF')
        return


    def run_save2google_earth(self, step_name):
        """Save velocity file in geo coordinates into Google Earth raster image."""
        if self.template['mintpy.save.kmz'] is True:
            print('creating Google Earth KMZ file for geocoded velocity file: ...')
            # input
            vel_file = 'velocity.h5'
            atr = readfile.read_attribute(vel_file)
            if 'Y_FIRST' not in atr.keys():
                vel_file = os.path.join(self.workDir, 'geo/geo_velocity.h5')

            # output
            kmz_file = '{}.kmz'.format(os.path.splitext(vel_file)[0])
            iargs = [vel_file, '-o', kmz_file]
            print('save_kmz.py', ' '.join(iargs))

            # update mode
            try:
                fbase = os.path.basename(kmz_file)
                kmz_file = [i for i in [fbase, './geo/{}'.format(fbase), './pic/{}'.format(fbase)] 
                            if os.path.isfile(i)][0]
            except:
                kmz_file = None
            if ut.run_or_skip(out_file=kmz_file, in_file=vel_file, check_readable=False) == 'run':
                mintpy.save_kmz.main(iargs)
        else:
            print('save velocity to Google Earth format is OFF.')
        return


    def run_save2hdfeos5(self, step_name):
        """Save displacement time-series and its aux data in geo coordinate into HDF-EOS5 format"""
        if self.template['mintpy.save.hdfEos5'] is True:
            # input
            ts_file = self.get_timeseries_filename(self.template)[step_name]['input']
            # Add attributes from custom template to timeseries file
            if self.customTemplate is not None:
                ut.add_attribute(ts_file, self.customTemplate)

            tcoh_file = 'temporalCoherence.h5'
            mask_file = 'geo_maskTempCoh.h5'
            geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[2]
            if 'geo' in ts_file:
                tcoh_file = './geo/geo_temporalCoherence.h5'
                mask_file = './geo/geo_maskTempCoh.h5'
                geom_file = './geo/geo_{}'.format(os.path.basename(geom_file))

            # cmd
            print('--------------------------------------------')
            iargs = [ts_file,
                     '-c', tcoh_file, 
                     '-m', mask_file,
                     '-g', geom_file, 
                     '-t', self.templateFile]
            print('save_hdfeos5.py', ' '.join(iargs))

            # output (check existing file)
            atr = readfile.read_attribute(ts_file)
            SAT = sensor.get_unavco_mission_name(atr)
            try:
                hdfeos5_file = ut.get_file_list('{}_*.he5'.format(SAT))[0]
            except:
                hdfeos5_file = None
            if ut.run_or_skip(out_file=hdfeos5_file, in_file=[ts_file, tcoh_file, mask_file, geom_file]) == 'run':
                mintpy.save_hdfeos5.main(iargs)
        else:
            print('save time-series to HDF-EOS5 format is OFF.')
        return


    def plot_result(self, print_aux=True, plot=True):
        """Plot data files and save to figures in pic folder"""
        if self.template['mintpy.plot'] and plot:
            print('\n******************** plot & save to pic ********************')
            print(self.plot_sh_cmd)
            subprocess.Popen(self.plot_sh_cmd, shell=True).wait()

        # message for more visualization scripts
        msg = """Explore more info & visualization options with the following scripts:
        info.py                    #check HDF5 file structure and metadata
        view.py                    #2D map view
        tsview.py                  #1D point time-series (interactive)   
        transect.py                #1D profile (interactive)
        plot_coherence_matrix.py   #plot coherence matrix for one pixel (interactive)
        plot_network.py            #plot network configuration of the dataset    
        plot_transection.py        #plot 1D profile along a line of a 2D matrix (interactive)
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

            elif sname == 'modify_network':
                self.run_network_modification(sname, plot=plot)

            elif sname == 'reference_point':
                self.run_reference_point(sname)

            elif sname == 'quick_overview':
                self.run_quick_overview(sname)

            elif sname == 'correct_unwrap_error':
                self.run_unwrap_error_correction(sname)

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

        # plot result (show aux visualization message for multiple steps processing)
        print_aux = len(steps) > 1
        self.plot_result(print_aux=print_aux, plot=plot)

        # go back to original directory
        print('Go back to directory:', self.cwd)
        os.chdir(self.cwd)

        # message
        msg  = '\n################################################'
        msg += '\n   Normal end of smallbaselineApp processing!'
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
