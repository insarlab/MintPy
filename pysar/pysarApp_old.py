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
import argparse
import warnings
import subprocess
import numpy as np
from pysar.objects import sensor
from pysar.utils import readfile, utils as ut
from pysar import version


##########################################################################
EXAMPLE = """example:
  pysarApp.py                       #Run / Rerun
  pysarApp.py <template_file>       #Run / Rerun
  pysarApp.py -h / --help           #Help on usage
  pysarApp.py -h / --help --steps   #Help on 
  pysarApp.py -H                    #Print all template options

  # Run with --start/stop/dostep
  pysarApp.py GalapagosSenDT128.template --dostep startup   #Do generate default_template from custom_template
  pysarApp.py GalapagosSenDT128.template --stop load_data   #End processing after loading data
"""


STEP_HELP="""
Command line options for steps processing are formed by
combining the following three options as required:

'--start=<step>', '--end=<step>', '--dostep=<step>'

The step names are chosen from the following list:

['startup', 'load_data', 'ref_point', 'stacking']
['unw_cor', 'net_modify', 'net_inversion', 'temp_coh']
['tropo', 'deramp', 'topo', 'resid_rms', 'ref_date', 'ts2vel']
['geocode', 'google_earth', 'hdfeos5', 'plot']

If --start is missing, then processing starts at the first step.
If --end is missing, then processing ends at the final step.
If --dostep is used, then only the named step is processed.

In order to use either --start or --dostep, it is necessary that a
previous run was done using one of the steps options to process at least
through the step immediately preceding the starting step of the current run.
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
                             'a) current directory, OR\n' +
                             'b) $SCRATCHDIR/projectName/PYSAR, if:\n' +
                             '    1) autoPath = True in pysar/defaults/auto_path.py AND\n' +
                             '    2) environment variable $SCRATCHDIR exists AND\n' +
                             '    3) input custom template named as: projectName.*\n')
    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='Generate default template (and merge with custom template), then exit.')
    parser.add_argument('-H', dest='print_auto_template', action='store_true',
                        help='Print/Show the example template file for routine processing.')
    parser.add_argument('-v','--version', action='store_true', help='print software version')

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

    # software version
    if inps.version:
        raise SystemExit(version.description)

    # full template
    if inps.print_auto_template:
        with open(inps.autoTemplateFile, 'r') as f:
            print(f.read())
        raise SystemExit()

    # ignore if pysarApp_template.txt is input as custom template
    if (inps.customTemplateFile 
            and os.path.basename(inps.customTemplateFile) == 'pysarApp_template.txt'):
        inps.customTemplateFile = None
    return inps
##########################################################################


##########################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    #########################################
    # Initiation
    #########################################
    print(version.logo)

    

    #########################################
    # Loading Data
    #########################################
    print('\n**********  Load Data  **********')
    ut.copy_aux_file(inps)

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
        #metaCmd = 'add_attribute.py {} {}'.format(inps.stackFile, inps.customTemplateFile)
        #print(metaCmd)
        #status = subprocess.Popen(metaCmd, shell=True).wait()
        # better control of special metadata, such as SUBSET_X/YMIN
        print('updating {} metadata based on custom template file: {}'.format(
            os.path.basename(inps.stackFile), inps.customTemplateFile))
        ut.add_attribute(inps.stackFile, customTemplate)

    if inps.load_dataset:
        if inps.plot:
            ut.plot_pysarApp(inps)
        raise SystemExit('Exit as planned after loading/checking the dataset.')


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print('\n**********  Select Reference Point  **********')
    # Initial mask (pixels with valid unwrapPhase or connectComponent in ALL interferograms)
    inps.maskFile = 'maskConnComp.h5'
    maskCmd = 'generate_mask.py {} --nonzero -o {} --update'.format(inps.stackFile, inps.maskFile)
    print(maskCmd)
    status = subprocess.Popen(maskCmd, shell=True).wait()

    # Average spatial coherence
    inps.avgSpatialCohFile = 'avgSpatialCoh.h5'
    avgCmd = 'temporal_average.py {i} --dataset coherence -o {o} --update'.format(i=inps.stackFile,
                                                                                  o=inps.avgSpatialCohFile)
    print(avgCmd)
    status = subprocess.Popen(avgCmd, shell=True).wait()

    # Select reference point
    refPointCmd = 'reference_point.py {} -t {} -c {}'.format(inps.stackFile,
                                                             inps.templateFile,
                                                             inps.avgSpatialCohFile)
    print(refPointCmd)
    status = subprocess.Popen(refPointCmd, shell=True).wait()
    if status is not 0:
        if inps.plot:
            plot_pysarApp(inps)
        raise Exception('Error while finding reference pixel in space.\n')

    ###############################################
    # Average velocity from interferogram stacking
    ###############################################
    print('\n**********  Quick assessment with interferogram stacking  **********')
    inps.waterMaskFile = 'waterMask.h5'
    if not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    # Average phase velocity - Stacking
    inps.avgPhaseVelFile = 'avgPhaseVelocity.h5'
    avgCmd = 'temporal_average.py {i} --dataset unwrapPhase -o {o} --update'.format(i=inps.stackFile,
                                                                                    o=inps.avgPhaseVelFile)
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


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    ut.correct_unwrap_error(inps, template)

    #########################################
    # Network Modification (Optional)
    #########################################
    print('\n**********  Modify Network  **********')
    networkCmd = 'modify_network.py {} -t {}'.format(inps.stackFile,
                                                     inps.templateFile)
    print(networkCmd)
    status = subprocess.Popen(networkCmd, shell=True).wait()
    if status is not 0:
        if inps.plot:
            plot_pysarApp(inps)
        raise Exception('Error while modifying the network of interferograms.\n')

    # Plot network colored in spatial coherence
    print('--------------------------------------------------')
    plotCmd = 'plot_network.py {} --template {} --nodisplay'.format(inps.stackFile,
                                                                    inps.templateFile)
    print(plotCmd)
    inps.cohSpatialAvgFile = '{}_coherence_spatialAvg.txt'.format(
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
        if inps.plot:
            plot_pysarApp(inps)
        raise SystemExit('Exit as planned after network modification.')

    #########################################
    # Inversion of Interferograms
    ########################################
    print('\n**********  Invert Network of Interferograms into Time-series  **********')
    invCmd = 'ifgram_inversion.py {} --template {} --update '.format(inps.stackFile,
                                                                     inps.templateFile)
    print(invCmd)
    inps.timeseriesFile = 'timeseries.h5'
    inps.tempCohFile = 'temporalCoherence.h5'
    inps.timeseriesFiles = ['timeseries.h5']       #all ts files
    status = subprocess.Popen(invCmd, shell=True).wait()
    if status is not 0:
        if inps.plot:
            plot_pysarApp(inps)
        raise Exception('Error while inverting network interferograms into timeseries')

    print('\n--------------------------------------------')
    print('Update Mask based on Temporal Coherence ...')
    ut.get_temporal_coherence_mask(inps, template)

    if inps.invert_network:
        if inps.plot:
            plot_pysarApp(inps)
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
                if inps.plot:
                    plot_pysarApp(inps)
                raise Exception('Error while correcting Local Oscillator Drift.\n')
        inps.timeseriesFile = outName
        inps.timeseriesFiles.append(outName)

    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print('\n**********  Tropospheric Delay Correction  **********')
    ut.correct_tropospheric_delay(inps, template)

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
                    if inps.plot:
                        plot_pysarApp(inps)
                    raise Exception('Error while removing phase ramp for time-series.\n')
            inps.timeseriesFile = outName
            inps.timeseriesFiles.append(outName)
        else:
            msg = 'un-recognized phase ramp method: {}'.format(inps.derampMethod)
            msg += '\navailable ramp types:\n{}'.format(ramp_list)
            raise ValueError(msg)
    else:
        print('No phase ramp removal.')

    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print('\n**********  Topographic Residual (DEM error) Correction  **********')
    outName = os.path.splitext(inps.timeseriesFile)[0]+'_demErr.h5'
    topoCmd = 'dem_error.py {i} -t {t} -o {o} --update '.format(i=inps.timeseriesFile,
                                                                t=inps.templateFile,
                                                                o=outName)
    print(topoCmd)
    inps.timeseriesResFile = None
    if template['pysar.topographicResidual']:
        status = subprocess.Popen(topoCmd, shell=True).wait()
        if status is not 0:
            if inps.plot:
                plot_pysarApp(inps)
            raise Exception('Error while correcting topographic phase residual.\n')
        inps.timeseriesFile = outName
        inps.timeseriesResFile = 'timeseriesResidual.h5'
        inps.timeseriesFiles.append(outName)
    else:
        print('No correction for topographic residuals.')

    ##############################################
    # Phase Residual for Noise Evaluation
    ##############################################
    # Timeseries Residual RMS
    print('\n**********  Timeseries Residual Root Mean Square  **********')
    if inps.timeseriesResFile:
        rmsCmd = 'timeseries_rms.py {} -t {}'.format(inps.timeseriesResFile,
                                                     inps.templateFile)
        print(rmsCmd)
        status = subprocess.Popen(rmsCmd, shell=True).wait()
        if status is not 0:
            if inps.plot:
                plot_pysarApp(inps)
            raise Exception('Error while calculating RMS of residual phase time-series.\n')
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
            if inps.plot:
                plot_pysarApp(inps)
            raise Exception('Error while changing reference date.\n')
    else:
        print('No reference change in time.')

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
        if inps.plot:
            plot_pysarApp(inps)
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
            geo_dir = os.path.join(inps.workDir, 'GEOCODE')
            geocode_script = os.path.join(os.path.dirname(__file__), 'geocode.py')
            geoCmd = ('{scp} {v} {c} {t} {g} -l {l} -t {e}'
                      ' --outdir {d} --update').format(scp=geocode_script,
                                                       v=inps.velFile,
                                                       c=inps.tempCohFile,
                                                       t=inps.timeseriesFile,
                                                       g=inps.geomFile,
                                                       l=inps.lookupFile,
                                                       e=inps.templateFile,
                                                       d=geo_dir)
            print(geoCmd)
            status = subprocess.Popen(geoCmd, shell=True).wait()
            if status is not 0:
                if inps.plot:
                    plot_pysarApp(inps)
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
        kmlCmd = 'save_kmz.py {} -o {}'.format(inps.velFile, outName)
        print(kmlCmd)
        try:
            outFile = [i for i in [outName, 'PIC/{}'.format(outName)] if os.path.isfile(i)][0]
        except:
            outFile = None
        if ut.run_or_skip(out_file=outFile, in_file=inps.velFile, check_readable=False) == 'run':
            status = subprocess.Popen(kmlCmd, shell=True).wait()
            if status is not 0:
                if inps.plot:
                    plot_pysarApp(inps)
                raise Exception('Error while generating Google Earth KMZ file.')

    #############################################
    # Save Timeseries to HDF-EOS5 format
    #############################################
    if template['pysar.save.hdfEos5'] is True:
        print('\n**********  Save Time-series in HDF-EOS5 Format  **********')
        ut.save_hdfeos5(inps, customTemplate)

    #############################################
    # Plot Figures
    #############################################
    if inps.plot:
        ut.plot_pysarApp(inps)

    #############################################
    # Timing                                    #
    #############################################
    m, s = divmod(time.time()-start_time, 60)
    print('\n###############################################')
    print('End of PySAR Routine Processing Workflow!')
    print('###############################################\n')
    print('time used: {:02.0f} mins {:02.1f} secs'.format(m, s))


###########################################################################################
if __name__ == '__main__':
    main()
