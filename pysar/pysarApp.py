#!/usr/bin/env python3
###############################################################################
# Project: PySAR 
# Purpose: Python Module for InSAR Time Series Analysis
# Author: Heresh Fattahi, Zhang Yunjun
# Created: July 2013
# Copyright (c) 2013, Heresh Fattahi, Zhang Yunjun
###############################################################################
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
###############################################################################


import os, sys, glob
import time
import argparse
import warnings
import shutil
import subprocess

import h5py
import numpy as np

import pysar
from pysar.utils import readfile, writefile, utils as ut
from pysar.objects import ifgramStack
from pysar.defaults import autoPath
from pysar import geocode, subset, save_hdfeos5 as hdfeos5


##########################################################################
#generate_from: http://patorjk.com/software/taag/
LOGO='''
_________________________________________________
       ____             __     __     ____  
       /    )         /    )   / |    /    )
------/____/----------\-------/__|---/___ /------
     /        /   /    \     /   |  /    |  
____/________(___/_(____/___/____|_/_____|_______
                /                           
            (_ /                            

 A Python package for InSAR time series analysis.
               PySAR v{}, {}
 Geodesy Lab, University of Miami, Maimi FL, USA
_________________________________________________
'''.format(pysar.release_version, pysar.release_date)

TEMPLATE='''# vim: set filetype=cfg:
##------------------------ pysarApp_template.txt ------------------------##
########## 1. Load Data (--load to exit after this step)
## auto - automatic path pattern for Univ of Miami file structure
## load_data.py -H to check more details and example inputs.
pysar.load.processor      = auto  #[isce,roipac,gamma,], auto for isce
##---------interferogram datasets:
pysar.load.unwFile        = auto  #[path2unw_file]
pysar.load.corFile        = auto  #[path2cor_file]
pysar.load.connCompFile   = auto  #[path2conn_file]
pysar.load.intFile        = auto  #[path2int_file]
##---------geometry datasets:
pysar.load.demFile        = auto  #[path2hgt_file]
pysar.load.lookupYFile    = auto  #[path2lat_file]]
pysar.load.lookupXFile    = auto  #[path2lon_file]
pysar.load.incAngleFile   = auto  #[path2los_file]
pysar.load.headAngleFile  = auto  #[path2los_file]
pysar.load.shadowMaskFile = auto  #[path2shadow_file]
pysar.load.bperpFile      = auto  #[path2bperp_file]

## 1.1 Subset (optional)
## if both yx and lalo are specified, use lalo option unless a) no lookup file AND b) dataset is in radar coord
pysar.subset.yx       = auto    #[1800:2000,700:800 / no], auto for no
pysar.subset.lalo     = auto    #[31.5:32.5,130.5:131.0 / no], auto for no
pysar.subset.tightBox = auto    #[yes / no], auto for yes, tight bounding box for files in geo coord

## 1.3 Reference in Space
## reference all interferograms to one common point in space
## auto - randomly select a pixel with coherence > minCoherence
pysar.reference.yx            = auto   #[257,151 / auto]
pysar.reference.lalo          = auto   #[31.8,130.8 / auto]
pysar.reference.coherenceFile = auto   #[filename], auto for avgSpatialCoherence.h5
pysar.reference.minCoherence  = auto   #[0.0-1.0], auto for 0.85, minimum coherence for auto method
pysar.reference.maskFile      = auto   #[filename / no], auto for mask.h5


## 1.4 Unwrapping Error Correction (optional and not recommended)
## unwrapping error correction based on the following two methods:
## a. phase closure (Fattahi, 2015, PhD Thesis)
## b. connecting bridge
pysar.unwrapError.method   = auto   #[bridging / phase_closure / no], auto for no
pysar.unwrapError.maskFile = auto   #[filename / no], auto for no
pysar.unwrapError.ramp     = auto   #[plane / quadratic], auto for plane
pysar.unwrapError.yx       = auto   #[y1_start,x1_start,y1_end,x1_end;y2_start,...], auto for none


########## 2. Network Inversion
## 2.1 Modify Network (optional)
## Coherence-based network modification = MST + Threshold, by default
## 1) calculate a average coherence value for each interferogram using spatial coherence and input mask (with AOI)
## 2) find a minimum spanning tree (MST) network with inverse of average coherence as weight (keepMinSpanTree)
## 3) for all interferograms except for MST's, exclude those with average coherence < minCoherence.
pysar.network.coherenceBased  = auto  #[yes / no], auto for yes, exclude interferograms with coherence < minCoherence
pysar.network.keepMinSpanTree = auto  #[yes / no], auto for yes, keep interferograms in Min Span Tree network
pysar.network.minCoherence    = auto  #[0.0-1.0], auto for 0.5
pysar.network.maskFile        = auto  #[file name, no], auto for mask.h5, no for all pixels
pysar.network.maskAoi.yx      = auto  #[y0:y1,x0:x1 / no], auto for no, area of interest for coherence calculation
pysar.network.maskAoi.lalo    = auto  #[lat0:lat1,lon0:lon1 / no], auto for no - use the whole area

## Network modification based on temporal/perpendicular baselines, date etc.
pysar.network.tempBaseMax     = auto  #[1-inf, no], auto for no, maximum temporal baseline in days
pysar.network.perpBaseMax     = auto  #[1-inf, no], auto for no, maximum perpendicular spatial baseline in meter
pysar.network.referenceFile   = auto  #[date12_list.txt / Modified_unwrapIfgram.h5 / no], auto for no
pysar.network.excludeDate     = auto  #[20080520,20090817 / no], auto for no
pysar.network.excludeIfgIndex = auto  #[1:5,25 / no], auto for no, list of ifg index (start from 0)
pysar.network.startDate       = auto  #[20090101 / no], auto for no
pysar.network.endDate         = auto  #[20110101 / no], auto for no


## 2.2 Invert network of interferograms into time series using weighted least sqaure (WLS) estimator.
## Temporal coherence is calculated using Tazzani et al. (Tizzani et al., 2007, IEEE-TGRS)
## Singular-Value Decomposition (SVD) is applied if network are not fully connected for no weight scenario.
## There are 4 weighting options:
## a. fim       - WLS, use Fisher Information Matrix as weight (Seymour & Cumming, 1994, IGARSS). [Recommended]
## b. variance  - WLS, use inverse of covariance as weight (Guarnieri & Tebaldini, 2008, TGRS)
## c. coherence - WLS, use coherence as weight (Perissin & Wang, 2012, IEEE-TGRS)
## d. no        - LS, no/uniform weight (Berardino et al., 2002, TGRS)
pysar.networkInversion.weightFunc    = auto #[fim / variance / coherence / no], auto for no
pysar.networkInversion.waterMaskFile = auto #[filename / no], auto for no
pysar.networkInversion.residualNorm  = auto #[L2 ], auto for L2, norm minimization solution
pysar.networkInversion.minTempCoh    = auto #[0.0-1.0], auto for 0.7, min temporal coherence for mask


########## Local Oscillator Drift (LOD) Correction (for Envisat only)
## reference: Marinkovic and Larsen, 2013, Proc. LPS
## correct LOD if input dataset comes from Envisat
## skip this step for all the other satellites.


########## 3. Tropospheric Delay Correction (optional and recommended)
## correct tropospheric delay using the following methods:
## a. pyaps - use weather re-analysis data (Jolivet et al., 2011, GRL, need to install PyAPS)
## b. height_correlation - correct stratified tropospheric delay (Doin et al., 2009, J Applied Geop)
## c. base_trop_cor - (not recommend) baseline error and stratified tropo simultaneously (Jo et al., 2010, Geo J)
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.weatherModel = auto  #[ERA / MERRA / NARR], auto for ECMWF, for pyaps method
pysar.troposphericDelay.polyOrder    = auto  #[1 / 2 / 3], auto for 1, for height_correlation method
pysar.troposphericDelay.looks        = auto  #[1-inf], auto for 8, for height_correlation, number of looks applied to
                                             #interferogram for empirical estimation of topography correlated atmosphere.


########## 4. Topographic Residual (DEM Error) Correction (optional and recommended)
## reference: Fattahi and Amelung, 2013, IEEE-TGRS
## Specify stepFuncDate option if you know there are sudden displacement jump in your area,
## i.e. volcanic eruption, or earthquake, and check timeseriesStepModel.h5 afterward for their estimation.
pysar.topographicResidual               = auto  #[yes / no], auto for yes
pysar.topographicResidual.polyOrder     = auto  #[1-inf], auto for 2, poly order of temporal deformation model
pysar.topographicResidual.stepFuncDate  = auto  #[20080529,20100611 / no], auto for no, date of step jump
pysar.topographicResidual.excludeDate   = auto  #[20070321 / txtFile / no], auto for no, date exlcuded for error estimation
pysar.topographicResidual.phaseVelocity = auto  #[yes / no], auto for no - phase, use phase velocity for error estimation

## 4.1 Phase Residual Root Mean Square
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
## Recommendation: quadratic for whole image, plane for local/small area
pysar.residualRms.maskFile        = auto  #[filename / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)

## 4.2 Select Reference Date
## reference all timeseries to one date in time
## minRMS - choose date with minimum residual RMS using value from step 8.1
## no     - do not change the default reference date (1st date)
pysar.reference.date = auto   #[reference_date.txt / 20090214 / minRMS / no], auto for reference_date.txt


########## 5. Phase Ramp Removal (optional)
## remove phase ramp for each epoch, useful to check localized deformation, i.e. volcanic, land subsidence, etc.
## [plane, quadratic, plane_range, plane_azimuth, quadratic_range, quadratic_azimuth, baseline_cor, base_trop_cor]
pysar.deramp          = auto  #[no / plane / quadratic], auto for no - no ramp will be removed
pysar.deramp.maskFile = auto  #[filename / no], auto for maskTempCoh.h5, mask file for ramp estimation


########## 6. Velocity Inversion
## estimate linear velocity from timeseries, and from tropospheric delay file if exists.
pysar.velocity.excludeDate = auto   #[exclude_date.txt / 20080520,20090817 / no], auto for exclude_date.txt
pysar.velocity.startDate   = auto   #[20070101 / no], auto for no
pysar.velocity.endDate     = auto   #[20101230 / no], auto for no


########## 7. Post-processing (geocode, output to Google Earth, HDF-EOS5, etc.)
## 7.1 Geocode
## For data processed by ROI_PAC/Gamma, output resolution for geocoded file is the same as their lookup table file.
## For data processed by ISCE/Doris, output resolution is assign by user with resolution option:
## 1) float number - resolution in degree, 0.001 by default, around 100 m on equator
## 2) file name    - use the resolution from a file in geo coordinates, e.g. demGeo.h5
pysar.geocode            = auto  #[yes / no], auto for yes
pysar.geocode.resolution = auto  #[0.0-inf / filename], auto for 0.001 (~100 m), output resolution for ISCE processor

## 7.2 Export to other formats
pysar.save.hdfEos5         = auto   #[yes / no], auto for no, save timeseries to HDF-EOS5 format
pysar.save.hdfEos5.update  = auto   #[yes / no], auto for no, put XXXXXXXX as endDate in output filename
pysar.save.hdfEos5.subset  = auto   #[yes / no], auto for no, put subset range info   in output filename
pysar.save.kml     = auto   #[yes / no], auto for yes, save geocoded velocity to Google Earth KMZ file
pysar.save.geotiff = auto   #[yes / no], auto for no, save geocoded velocity to Geotiff format [not implemented yet]

## 7.3 Plot
pysar.plot = auto   #[yes / no], auto for yes, plot files generated by pysarApp default processing to PIC folder
'''

EXAMPLE='''example:
  pysarApp.py
  pysarApp.py  SanAndreasT356EnvD.template
  pysarApp.py  SanAndreasT356EnvD.template  --load-data
  pysarApp.py  SanAndreasT356EnvD.template  --dir ~/insarlab/SanAndreasT356EnvD/PYSAR

  # Generate template file:
  pysarApp.py -g
  pysarApp.py SanAndreasT356EnvD.template -g

  -----------------------------------------------------
  Read pysar_template.txt file for more option details.
  -----------------------------------------------------
'''

UM_FILE_STRUCT='''
    scratch/                 # $SCRATCHDIR defined in environmental variable
        SanAndreasT356EnvD/  # my_projectName, same as the basename of template file
            DEM/             # DEM file(s) (for topographic phase and geocode)
            DOWNLOAD/        # (optional) Data downloaded from agencies
            PROCESS/         # Interferograms processed by ROI_PAC, Gamma, ISCE, ... 
            PYSAR/           # PySAR work directory for time series analysis
                subset/      # PySAR subset
            RAW/             # (optional) Raw SAR data untared from DOWNLOAD directory
            SLC/             # (optional) SLC SAR data after focusing from RAW directory
            WEATHER/         # Weather data (e.g. PyAPS products)
                ECMWF/
                MERRA/
'''


def createParser():
    parser = argparse.ArgumentParser(description='Time Series Analysis Routine',
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('templateFileCustom', nargs='?',\
                        help='custom template with option settings.\n'+\
                             "It's equivalent to None, if pysarApp_template.txt is input, as it will be read always.")
    parser.add_argument('--dir', dest='workDir',\
                        help='PySAR working directory, default is:\n'+\
                             'a) current directory, or\n'+\
                             'b) $SCRATCHDIR/projectName/PYSAR, if meets the following 3 requirements:\n'+\
                             '    1) autoPath = True in pysar/defaults/auto_path.py\n'+\
                             '    2) environmental variable $SCRATCHDIR exists\n'+\
                             '    3) input custom template with basename same as projectName\n')
    parser.add_argument('-g', dest='generate_template', action='store_true',\
                        help='Generate default template (and merge with custom template), then exit.')
    parser.add_argument('--reset', action='store_true',\
                        help='Reset files attributes to re-run pysarApp.py after loading data by:\n'+\
                             '    1) removing ref_y/x/lat/lon for unwrapIfgram.h5 and coherence.h5\n'+\
                             '    2) set DROP_IFGRAM=no for unwrapIfgram.h5 and coherence.h5')
    parser.add_argument('--load-data', dest='load_dataset', action='store_true',\
                        help='Step 1. Load/check dataset, then exit')
    parser.add_argument('--modify-network', dest='modify_network', action='store_true',\
                        help='Step 4. Modify the network, then exit')
    parser.add_argument('--invert-network', dest='invert_network', action='store_true',\
                        help='Step 5. Inverse network of interferograms into time-series, then exit')
    return parser


def cmdLineParse(iargs=None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    if inps.templateFileCustom and os.path.basename(inps.templateFileCustom) == 'pysarApp_template.txt':
        inps.templateFileCustom = None
    return inps


###############################################################################
def check_geocode_file(lookupFile, File, templateFile=None, outFile=None):
    '''Geocode input file or use existed geocoded file.'''
    if not File:
        return None
    else:
        atr = readfile.read_attribute(File)
        if 'Y_FIRST' in atr.keys():
            return File
    
    if not lookupFile:
        warnings.warn('No lookup file found! Skip geocoding.')
        return None

    if not outFile:
        outFile = geocode.geocode_output_filename(File)

    if ut.update_file(outFile, File):
        geocodeCmd = 'geocode.py %s -l %s -o %s' % (File, os.path.basename(lookupFile), outFile)
        if templateFile:
            geocodeCmd += ' -t '+templateFile
        print(geocodeCmd)
        status = subprocess.Popen(geocodeCmd, shell=True).wait()

    try:    outFile = glob.glob(outFile)[0]
    except: outFile = None
    return outFile


def copy_aux_file(inps):
    #####for Univ of Miami
    fileList = ['PROCESS/unavco_attributes.txt', 'PROCESS/bl_list.txt','SLC/summary*slc.jpg']
    try:
        projectDir = os.path.join(os.getenv('SCRATCHDIR'), inps.projectName)
        fileList = ut.get_file_list([os.path.join(projectDir, i) for i in fileList], abspath=True)
        for file in fileList:
            if ut.update_file(os.path.basename(file), file, check_readable=False):
                shutil.copy2(file, inps.workDir)
                print('copy {} to work directory'.format(os.path.basename(file)))
    except:
        pass
    return inps


def read_template(inps):
    print('\n**********  Read Template File  **********')
    ##### default template
    inps.templateFile = os.path.join(inps.workDir, 'pysarApp_template.txt')
    if not os.path.isfile(inps.templateFile):
        print('generate default template file: '+inps.templateFile)
        f = open(inps.templateFile, 'w')
        f.write(TEMPLATE)
        f.close()
    else:
        print('default template file exists: '+inps.templateFile)

    ##### custom template
    if inps.templateFileCustom:
        # Copy custom template file to work directory
        if ut.update_file(os.path.basename(inps.templateFileCustom), inps.templateFileCustom, check_readable=False):
            shutil.copy2(inps.templateFileCustom, inps.workDir)

        # Read custom template
        print('read custom template file: '+inps.templateFileCustom)
        templateCustom = readfile.read_template(inps.templateFileCustom)
        # correct some loose type errors
        for key in templateCustom.keys():
            if   templateCustom[key].lower() in ['default']:          templateCustom[key] = 'auto'
            elif templateCustom[key].lower() in ['n','off','false']:  templateCustom[key] = 'no'
            elif templateCustom[key].lower() in ['y','on','true']:    templateCustom[key] = 'yes'
        for key in ['pysar.deramp', 'pysar.troposphericDelay.method']:
            if key in templateCustom.keys():
                templateCustom[key] = templateCustom[key].lower().replace('-','_')

        # FA 1/18: insert general options from template file as pysar.* options
        if 'processor' in templateCustom.keys():
            templateCustom['pysar.load.processor'] = templateCustom['processor']

        # Update default template with custom input template
        print('update default template based on input custom template')
        inps.templateFile = ut.update_template_file(inps.templateFile, templateCustom)

    if inps.generate_template:
        print('Exit as planned after template file generation.')
        sys.exit(0)

    print('read default template file: '+inps.templateFile)
    template = readfile.read_template(inps.templateFile)
    template = ut.check_template_auto_value(template)

    # Get existing files name: unavco_attributes.txt
    try:    inps.unavcoMetadataFile = ut.get_file_list('unavco_attribute*txt', abspath=True)[0]
    except: inps.unavcoMetadataFile = None;   print('No UNAVCO attributes file found.')

    return inps, template


##########################################################################
def main(iargs=None):
    start = time.time()
    inps = cmdLineParse(iargs)
    #########################################
    # Initiation
    #########################################
    print(LOGO)

    # Project Name
    inps.projectName = None
    if inps.templateFileCustom:
        inps.templateFileCustom = os.path.abspath(inps.templateFileCustom)
        inps.projectName = os.path.splitext(os.path.basename(inps.templateFileCustom))[0]
        print('Project name: '+inps.projectName)

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
    print("Go to work directory: "+inps.workDir)

    copy_aux_file(inps)

    inps, template = read_template(inps)


    #########################################
    # Loading Data
    #########################################
    print('\n**********  Load Data  **********')
    loadCmd = 'load_data.py --template {}'.format(inps.templateFile)
    if inps.templateFileCustom:
        loadCmd += ' {} --project {}'.format(inps.templateFileCustom, inps.projectName)
    print(loadCmd)
    status = subprocess.Popen(loadCmd, shell=True).wait()
    os.chdir(inps.workDir)

    print('-'*50)
    inps, atr = ut.check_loaded_dataset(inps.workDir, inps)
    if not inps.stackFile:
        sys.exit('ERROR: No interferograms stack file found!')

    ##Add template options into HDF5 file metadata
    if inps.templateFileCustom:
        atrCmd = 'add_attribute.py {} {}'.format(inps.stackFile, inps.templateFileCustom)
        print(atrCmd)
        status = subprocess.Popen(atrCmd, shell=True).wait()
    #ut.add_attribute(inps.stackFile, template)

    if inps.load_dataset:
        print('Exit as planned after loading/checking the dataset.')
        sys.exit(0)

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


    stackobj = ifgramStack(inps.stackFile)
    stackobj.open()
    #########################################
    # Generating Aux files
    #########################################
    print('\n**********  Generate Auxiliary Files  **********')
    ##### Initial mask (pixels with valid unwrapPhase or connectComponent in ALL interferograms)
    inps.maskFile = 'mask.h5'
    if ut.update_file(inps.maskFile, inps.stackFile):
        maskCmd = 'generate_mask.py {} --nonzero -o {}'.format(inps.stackFile, inps.maskFile)
        print(maskCmd)
        status = subprocess.Popen(maskCmd, shell=True).wait()

    ##### Average spatial coherence
    inps.avgSpatialCohFile = 'avgSpatialCoherence.h5'
    if ut.update_file(inps.avgSpatialCohFile, inps.stackFile):
        avgCmd = 'temporal_average.py {} --dataset coherence -o {}'.format(inps.stackFile, inps.avgSpatialCohFile)
        print(avgCmd)
        status = subprocess.Popen(avgCmd, shell=True).wait()


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print('\n**********  Select Reference Point  **********')
    refPointCmd = 'reference_point.py {} -t {} -c {}'.format(inps.stackFile, inps.templateFile,inps.avgSpatialCohFile)
    print(refPointCmd)
    status = subprocess.Popen(refPointCmd, shell=True).wait()
    if status is not 0:
        print('\nError while finding reference pixel in space.\n')
        sys.exit(-1)


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    if template['pysar.unwrapError.method']:
        print('\n**********  Unwrapping Error Correction  **********')
        outName = os.path.splitext(inps.stackFile)[0]+'_unwCor.h5'
        unwCmd='unwrap_error.py '+inps.stackFile+' --mask '+inps.maskFile+' --template '+inps.templateFile
        print(unwCmd)
        if ut.update_file(outName, inps.stackFile):
            print('This might take a while depending on the size of your data set!')
            status = subprocess.Popen(unwCmd, shell=True).wait()
            if status is not 0:
                print('\nError while correcting phase unwrapping errors.\n')
                sys.exit(-1)
        inps.stackFile = outName


    #########################################
    # Network Modification (Optional)
    #########################################
    print('\n**********  Modify Network  **********')
    networkCmd = 'modify_network.py {} -t {}'.format(inps.stackFile, inps.templateFile)
    print(networkCmd)
    print('-'*50+'\nTo use/restore all interferograms in the file, run modify_network.py with --reset option.\n'+'-'*50)
    status = subprocess.Popen(networkCmd, shell=True).wait()
    if status is not 0:
        print('\nError while modifying the network of interferograms.\n')
        sys.exit(-1)

    # Plot network colored in spatial coherence
    print('--------------------------------------------------')
    plotCmd = 'plot_network.py {} --template {} --nodisplay'.format(inps.stackFile, inps.templateFile)
    print(plotCmd)
    inps.cohSpatialAvgFile = '{}_coherence_spatialAverage.txt'.format(os.path.splitext(os.path.basename(inps.stackFile))[0])
    if ut.update_file('Network.pdf',check_readable=False,inFile=[inps.stackFile,inps.cohSpatialAvgFile,inps.templateFile]):
        status = subprocess.Popen(plotCmd, shell=True).wait()

    if inps.modify_network:
        sys.exit('Exit as planned after network modification.')


    #########################################
    # Inversion of Interferograms
    ########################################
    print('\n**********  Invert Network of Interferograms into Time-series  **********')
    invCmd = 'ifgram_inversion.py {} --template {}'.format(inps.stackFile, inps.templateFile)
    print(invCmd)
    inps.timeseriesFile = 'timeseries.h5'
    inps.tempCohFile = 'temporalCoherence.h5'
    if ut.update_file(inps.timeseriesFile, inps.stackFile):
        status = subprocess.Popen(invCmd, shell=True).wait()
        if status is not 0:
            print('\nError while inverting network of interferograms to time-series.\n')
            sys.exit(-1)

    print('\n--------------------------------------------')
    print('Update Mask based on Temporal Coherence ...')
    inps.maskFile = 'maskTempCoh.h5'
    inps.minTempCoh = template['pysar.networkInversion.minTempCoh']
    maskCmd = 'generate_mask.py {} -m {} -o {}'.format(inps.tempCohFile, inps.minTempCoh, inps.maskFile)
    print(maskCmd)
    if ut.update_file(inps.maskFile, inps.tempCohFile):
        status = subprocess.Popen(maskCmd, shell=True).wait()
        if status is not 0:
            print('\nError while generating mask file from temporal coherence.\n')
            sys.exit(-1)

    if inps.invert_network:
        sys.exit('Exit as planned after network inversion.')


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   for Envisat data in radar coord only
    ############################################## 
    if atr['PLATFORM'].lower().startswith('env'):
        print('\n**********  Local Oscillator Drift Correction for Envisat  **********')
        outName = os.path.splitext(inps.timeseriesFile)[0]+'_LODcor.h5'
        lodCmd = 'local_oscilator_drift.py {} -r {} -o {}'.format(inps.timeseriesFile, inps.geomFile, outName)
        print(lodCmd)
        if ut.update_file(outName, [inps.timeseriesFile, inps.geomFile]):
            status = subprocess.Popen(lodCmd, shell=True).wait()
            if status is not 0:
                print('\nError while correcting Local Oscillator Drift.\n')
                sys.exit(-1)
        inps.timeseriesFile = outName


    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print('\n**********  Tropospheric Delay Correction  **********')
    inps.tropPolyOrder = template['pysar.troposphericDelay.polyOrder']
    inps.tropModel     = template['pysar.troposphericDelay.weatherModel']
    inps.tropMethod    = template['pysar.troposphericDelay.method']
    if inps.tropMethod:
        ##Check Conflict with base_trop_cor
        if template['pysar.deramp'] == 'base_trop_cor':
            msg='''
            Method Conflict: base_trop_cor is in conflict with {} option!
            base_trop_cor applies simultaneous ramp removal AND tropospheric correction.
            IGNORE base_trop_cor input and continue pysarApp.py.
            '''
            warnings.warn(msg)
            template['pysar.deramp'] = False

        fbase = os.path.splitext(inps.timeseriesFile)[0]
        # Call scripts
        if inps.tropMethod == 'height_correlation':
            outName = '{}_tropHgt.h5'.format(fbase)
            print('tropospheric delay correction with height-correlation approach')
            tropCmd = 'tropcor_phase_elevation.py {} -d {} -p {} -m {} -o {}'.format(inps.timeseriesFile,\
                                                                                     inps.geomFile,\
                                                                                     inps.tropPolyOrder,\
                                                                                     inps.maskFile,\
                                                                                     outName)
            print(tropCmd)
            if ut.update_file(outName, inps.timeseriesFile):
                status = subprocess.Popen(tropCmd, shell=True).wait()
                if status is not 0:
                    print('\nError while correcting tropospheric delay.\n')
                    sys.exit(-1)
            inps.timeseriesFile = outName

        elif inps.tropMethod == 'pyaps':
            outName = '{}_{}.h5'.format(fbase, inps.tropModel)
            print('Atmospheric correction using Weather Re-analysis dataset (using PyAPS software)')
            print('Weather Re-analysis dataset: '+inps.tropModel)
            tropCmd = 'tropcor_pyaps.py -f {} --model {} --dem {} -i {} -w {}'.format(inps.timeseriesFile,\
                                                                                      inps.tropModel,\
                                                                                      inps.geomFile,\
                                                                                      inps.geomFile,
                                                                                      os.path.join(inps.workDir,'../WEATHER'))
            print(tropCmd)
            if ut.update_file(outName, inps.timeseriesFile):
                try:
                    fileList = [os.path.join(inps.workDir, 'INPUTS/{}.h5'.format(inps.tropModel))]
                    inps.tropFile = ut.get_file_list(fileList)[0]
                    tropCmd = 'diff.py {} {} -o {}'.format(inps.timeseriesFile, inps.tropFile, outName)
                    print('Use existed tropospheric delay file: {}'.format(inps.tropFile))
                    print(tropCmd)
                except:
                    pass
                status = subprocess.Popen(tropCmd, shell=True).wait()
                if status is not 0:
                    print('\nError while correcting tropospheric delay, try the following:')
                    print('1) Check the installation of PyAPS (http://earthdef.caltech.edu/projects/pyaps/wiki/Main)')
                    print('   Try in command line: python -c "import pyaps"')
                    print('2) Use other tropospheric correction method, height-correlation, for example')
                    print('3) or turn off the option by setting pysar.troposphericDelay.method = no in template file.\n')
                    sys.exit(-1)
            inps.timeseriesFile = outName
        else:
            print('No atmospheric delay correction.')

    ## Grab tropospheric delay file
    try:    inps.tropFile = ut.get_file_list(inps.tropModel+'.h5')[0]
    except: inps.tropFile = None


    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print('\n**********  Topographic Residual (DEM error) Correction  **********')
    outName = os.path.splitext(inps.timeseriesFile)[0]+'_demErr.h5'
    topoCmd = 'dem_error.py {} -g {} -t {} -o {}'.format(inps.timeseriesFile,inps.geomFile,inps.templateFile,outName)
    print(topoCmd)
    inps.timeseriesResFile = None
    if template['pysar.topographicResidual']:
        if ut.update_file(outName, inps.timeseriesFile):
            status = subprocess.Popen(topoCmd, shell=True).wait()
            if status is not 0:
                print('\nError while correcting topographic phase residual.\n')
                sys.exit(-1)
        inps.timeseriesFile = outName
        inps.timeseriesResFile = 'timeseriesResidual.h5'
    else:
        print('No correction for topographic residuals.')


    ##############################################
    # Timeseries Residual Standard Deviation
    ##############################################
    print('\n**********  Timeseries Residual Root Mean Square  **********')
    if inps.timeseriesResFile:
        rmsCmd = 'timeseries_rms.py {} -t {}'.format(inps.timeseriesResFile, inps.templateFile)
        print(rmsCmd)
        status = subprocess.Popen(rmsCmd, shell=True).wait()
        if status is not 0:
            print('\nError while calculating RMS of time series phase residual.\n')
            sys.exit(-1)
    else:
        print('No timeseries residual file found! Skip residual RMS analysis.')


    ##############################################
    # Reference in Time
    ##############################################
    print('\n**********  Select Reference Date  **********')
    if template['pysar.reference.date']:
        outName = '{}_refDate.h5'.format(os.path.splitext(inps.timeseriesFile)[0])
        refCmd = 'reference_date.py {} -t {} -o {}'.format(inps.timeseriesFile, inps.templateFile, outName)
        print(refCmd)
        if ut.update_file(outName, inps.timeseriesFile):
            status = subprocess.Popen(refCmd, shell=True).wait()
            if status is not 0:
                print('\nError while changing reference date.\n')
                sys.exit(-1)
        inps.timeseriesFile = outName
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
        if inps.geocoded and inps.derampMethod in ['baseline_cor','base_trop_cor']:
            warnings.warn('dataset is in geo coordinates, can not apply {} method'.format(inps.derampMethod))
            print('skip deramping and continue.')

        ##Get executable command and output name
        derampCmd = None
        fbase = os.path.splitext(inps.timeseriesFile)[0]
        if inps.derampMethod in ['plane', 'quadratic', 'plane_range', 'quadratic_range', 'plane_azimuth', 'quadratic_azimuth']:
            outName = '{}_{}.h5'.format(fbase, inps.derampMethod)
            derampCmd = 'remove_ramp.py {} -s {} -m {} -o {}'.format(inps.timeseriesFile, inps.derampMethod,\
                                                                     inps.derampMaskFile, outName)

        elif inps.derampMethod == 'baseline_cor':
            outName = '{}_baselineCor.h5'.format(fbase)
            derampCmd = 'baseline_error.py {} {}'.format(inps.timeseriesFile, inps.maskFile)

        elif inps.derampMethod in ['base_trop_cor','basetropcor','baselinetropcor']:
            print('Joint estimation of Baseline error and tropospheric delay [height-correlation approach]')
            outName = '{}_baseTropCor.h5'.format(fbase)
            derampCmd = 'baseline_trop.py {} {} {} range_and_azimuth {}'.format(inps.timeseriesFile,\
                                                                                inps.geomFile,\
                                                                                inps.tropPolyOrder,\
                                                                                inps.maskFile)
        else:
            warnings.warn('Unrecognized phase ramp method: '+template['pysar.deramp'])

        ##Execute command
        if derampCmd:
            print(derampCmd)
            if ut.update_file(outName, inps.timeseriesFile):
                status = subprocess.Popen(derampCmd, shell=True).wait()
                if status is not 0:
                    print('\nError while removing phase ramp for each acquisition of time-series.\n')
                    sys.exit(-1)
            inps.timeseriesFile = outName
    else:
        print('No phase ramp removal.')


    #############################################
    # Velocity and rmse maps
    #############################################
    print('\n**********  Estimate Velocity  **********')
    inps.velFile = 'velocity.h5'
    velCmd = 'timeseries2velocity.py {} -t {} -o {}'.format(inps.timeseriesFile, inps.templateFile, inps.velFile)
    print(velCmd)
    if ut.update_file(inps.velFile, [inps.timeseriesFile, inps.templateFile]):
        status = subprocess.Popen(velCmd, shell=True).wait()
        if status is not 0:
            print('\nError while estimating linear velocity from time-series.\n')
            sys.exit(-1)

    # Velocity from Tropospheric delay
    if inps.tropFile:
        suffix = os.path.splitext(os.path.basename(inps.tropFile))[0].title()
        inps.tropVelFile = '{}{}.h5'.format(os.path.splitext(inps.velFile)[0], suffix)
        velCmd = 'timeseries2velocity.py {} -t {} -o {}'.format(inps.tropFile, inps.templateFile, inps.tropVelFile)
        print(velCmd)
        if ut.update_file(inps.tropVelFile, [inps.tropFile, inps.templateFile]):
            status = subprocess.Popen(velCmd, shell=True).wait()


    ############################################
    # Post-processing
    # Geocodeing, masking and save to KML 
    ############################################
    print('\n**********  Post-processing  **********')
    ###### Geocoding
    if inps.geocoded:
        inps.geoVelFile = inps.velFile
        inps.geoTempCohFile = inps.tempCohFile
        inps.geoTimeseriesFile = inps.timeseriesFile
    else:
        inps.geoVelFile = None
        inps.geoTempCohFile = None
        inps.geoTimeseriesFile = None
        if template['pysar.geocode'] is True:
            print('\n--------------------------------------------')
            inps.geoVelFile        = check_geocode_file(inps.lookupFile, inps.velFile,        inps.templateFile)
            inps.geoTempCohFile    = check_geocode_file(inps.lookupFile, inps.tempCohFile,    inps.templateFile)
            inps.geoTimeseriesFile = check_geocode_file(inps.lookupFile, inps.timeseriesFile, inps.templateFile)


    ##### Mask in geo coord
    inps.geoMaskFile = None
    if inps.geoTempCohFile:
        # Generate mask in geo coord
        print('\n--------------------------------------------')
        outName = 'maskTempCoh.h5'
        if os.path.basename(inps.geoTempCohFile).startswith('geo_'):
            outName = 'geo_'+outName
        maskCmd = 'generate_mask.py {} -m {} -o {}'.format(inps.geoTempCohFile, inps.minTempCoh, outName)
        print(maskCmd)
        if ut.update_file(outName, inps.geoTempCohFile):
            status = subprocess.Popen(maskCmd, shell=True).wait()
        inps.geoMaskFile = outName

        # Mask geo_velocity file
        if inps.geoVelFile and inps.geoMaskFile:
            outName = os.path.splitext(inps.geoVelFile)[0]+'_masked.h5'
            maskCmd = 'mask.py {} -m {} -o {}'.format(inps.geoVelFile, inps.geoMaskFile, outName)
            print(maskCmd)
            if ut.update_file(outName, [inps.geoVelFile, inps.geoMaskFile]):
                status = subprocess.Popen(maskCmd, shell=True).wait()
            try: inps.geoVelFile = glob.glob(outName)[0]
            except: pass

    ##### Save to Google Earth KML file
    if inps.geoVelFile and template['pysar.save.kml'] is True:
        print('\n--------------------------------------------')
        print('creating Google Earth KMZ file for geocoded velocity file: {} ...'.format(inps.geoVelFile))
        outName = os.path.splitext(inps.geoVelFile)[0]+'.kmz'
        kmlCmd = 'save_kml.py {}'.format(inps.geoVelFile)
        print(kmlCmd)
        if ut.update_file(outName, inps.geoVelFile, check_readable=False):
            status = subprocess.Popen(kmlCmd, shell=True).wait()


    #############################################
    # Save Timeseries to HDF-EOS5 format
    #############################################
    if template['pysar.save.hdfEos5'] is True:
        print('\n**********  Save Time-series in HDF-EOS5 Format  **********')
        if 'Y_FIRST' not in atr.keys() and not inps.lookupFile:
            warnings.warn('Dataset is in radar coordinates without lookup table file.'+\
                          'Can not geocode.'+\
                          'Skip saving.')
        else:
            # 1. Time series file
            inps.geoTimeseriesFile = check_geocode_file(inps.lookupFile, inps.timeseriesFile, inps.templateFile)
            # Add HDF-EOS5 attributes
            if inps.unavcoMetadataFile:
                atrCmd = 'add_attribute.py '+inps.geoTimeseriesFile+' '+inps.unavcoMetadataFile
                print(atrCmd)
                status = subprocess.Popen(atrCmd, shell=True).wait()
                if status is not 0:
                    print('\nError while adding HDF-EOS5 attributes to time series file.\n')
                    sys.exit(-1)

            # 2. Temporal Coherence
            inps.geoTempCohFile = check_geocode_file(inps.lookupFile, inps.tempCohFile, inps.templateFile)

            # 3. Mask file
            if not inps.geoMaskFile:
                outName = 'maskTempCoh.h5'
                if os.path.basename(inps.geoTempCohFile).startswith('geo_'):
                    outName = 'geo_'+outName
                maskCmd = 'generate_mask.py '+inps.geoTempCohFile+' -m '+str(inps.minTempCoh)+' -o '+outName
                print(maskCmd)
                if ut.update_file(outName, inps.geoTempCohFile):
                    status = subprocess.Popen(maskCmd, shell=True).wait()
                    if status is not 0:
                        sys.exit('\nError while generating mask file.\n')
                inps.geoMaskFile = outName

            # Save to HDF-EOS5 format
            print('--------------------------------------------')
            SAT = hdfeos5.get_mission_name(atr)
            try:
                inps.hdfeos5_file = ut.get_file_list(SAT+'_*.he5')[0]
                print('Find existed HDF-EOS5 time-series file: '+inps.hdfeos5_file)
            except:
                inps.hdfeos5_file = None
                print('No HDF-EOS5 time-series file exists yet.')
            hdfeos5Cmd = 'save_hdfeos5.py '+inps.geoTimeseriesFile+' -t '+inps.templateFile+\
                         ' -c '+inps.geoTempCohFile+' -m '+inps.geoMaskFile
            print(hdfeos5Cmd)
            if ut.update_file(inps.hdfeos5_file, [inps.geoTimeseriesFile, inps.geoTempCohFile, inps.geoMaskFile,\
                                              inps.inc_angle_geo_file, inps.dem_geo_file], check_readable=False):
                status = subprocess.Popen(hdfeos5Cmd, shell=True).wait()
                #if status is not 0:
                #    sys.exit('\nError while generating HDF-EOS5 time-series file.\n')


    #############################################
    # Plot Figures
    #############################################
    inps.plot = template['pysar.plot']
    if inps.plot is True:
        print('\n**********  Plot Results / Save to PIC  **********')
        inps.plotShellFile = 'plot_pysarApp.sh'

        # Copy to workding directory if not existed yet.
        if not os.path.isfile('./'+inps.plotShellFile):
            print('copy $PYSAR_HOME/bin/{} to work directory: {}'.format(inps.plotShellFile, inps.workDir))
            try:
                shutil.copy2(ut.which(inps.plotShellFile), './')
            except:
                print('WARNING: no {} found in the environment variable path.'.format(inps.plotShellFile))
                print('Check if $PYSAR_HOME/bin is in $PATH, and re-run')
                inps.plot = False

    if inps.plot and ut.update_file('./PIC', [inps.plotShellFile, inps.templateFile], check_readable=False):
        plotCmd = './'+inps.plotShellFile
        print(plotCmd)
        status = subprocess.Popen(plotCmd, shell=True).wait()
        print('\n'+'-'*50)
        msg='''Here are some suggestions for better figures:
        1) Edit parameters in plot_pysarApp.sh and re-run this script.
        2) Play with view.py, tsview.py and save_kml.py for more advanced/customized figures.
        '''
        print(msg)


    #############################################
    # Time                                      #
    #############################################
    s = time.time()-start;  m, s = divmod(s, 60);  h, m = divmod(m, 60)
    print('\nTime used: %02d hours %02d mins %02d secs' % (h, m, s))
    print('\n###############################################')
    print('End of PySAR processing!')
    print('################################################\n')



###########################################################################################
if __name__ == '__main__':
    main()

