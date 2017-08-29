#! /usr/bin/env python2
###############################################################################
# 
# Project: PySAR 
# Purpose: Python Module for InSAR Time-series Analysis
# Author: Heresh Fattahi, Zhang Yunjun
# Created: July 2013
#
###############################################################################
# Copyright (c) 2013, Heresh Fattahi, Zhang Yunjun
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


import os
import sys
import glob
import time
import argparse
import warnings
import shutil

import h5py
import numpy as np

import pysar
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar.subset as subset
import pysar.save_unavco as unavco


def check_subset_file(File, inps_dict, outFile=None, overwrite=False):
    '''Subset input file or use existed subseted file.'''
    if not File:
        return None

    if not outFile:
        if os.getcwd() == inps_dict['work_dir']:
            outFile = 'subset_'+os.path.basename(File)
        else:
            # if current dir is not original timeseries directory (e.g. TIMESERIES/subset)
            # use the same filename (but in different directories)
            outFile = os.path.basename(File)

    if ut.update_file(outFile, File, overwrite):
        outFile = subset.subset_file(File, inps_dict, outFile)

    return outFile


def check_geocode_file(geomapFile, File, outFile=None):
    '''Geocode input file or use existed geocoded file.'''
    if not geomapFile:
        warnings.warn('No geomap*.trans file found! Skip geocoding.')
        return None
    if not File:  return None

    if not outFile:  outFile = 'geo_'+os.path.basename(File)

    if ut.update_file(outFile, File):
        geocodeCmd = 'geocode.py '+File+' -l '+os.path.basename(geomapFile)
        print geocodeCmd
        try: os.system(geocodeCmd)
        except: pass

    try:    outFile = glob.glob(outFile)[0]
    except: outFile = None
    return outFile


def subset_dataset(inps, geo_box4geo, pix_box4rdr):
    '''Subset all file within dataset
    with geo_box4geo for all geocoded file and pix_box4rdr for all files in radar coord.
    Inputs:
        inps - Namespace with all files that needs to be subseted and work_dir
        geo_box4geo - tuple of 4 float, subset range in lat/lon for files in geo coord
        pix_box4rdr - tuple of 4 int, subset range in y/x for files in radar coord
    Output:
        inps - Namespace, update file name/path info
    '''
    # Subset directory
    subset_dir = inps.work_dir+'/SUBSET_Y%d_%dX%d_%d' % (pix_box4rdr[1], pix_box4rdr[3], pix_box4rdr[0], pix_box4rdr[2])
    if not os.path.isdir(subset_dir):
        os.mkdir(subset_dir)
    os.chdir(subset_dir)
    print '\n--------------------------------------------'
    print 'Creating subset datasets ...'
    print "Go to subset directory: "+subset_dir

    # Subset files
    print '--------------------------------------------'
    print 'subseting dataset in radar coord pix_box4rdr: '+str(pix_box4rdr)
    inps = subset.subset_box2inps(inps, pix_box4rdr, None)
    inps.ifgram_file    = check_subset_file(inps.ifgram_file, vars(inps))
    inps.dem_radar_file = check_subset_file(inps.dem_radar_file, vars(inps))
    inps.coherence_file = check_subset_file(inps.coherence_file, vars(inps))
    inps.trop_file      = check_subset_file(inps.trop_file, vars(inps))

    print '--------------------------------------------'
    print 'subseting dataset in geo coord geo_box4geo: '+str(geo_box4geo)
    inps = subset.subset_box2inps(inps, None, geo_box4geo)
    inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps))
    inps.trans_file  = check_subset_file(inps.trans_file, vars(inps))

    # adjust trans file value due to subsetted radar file
    if inps.trans_file:
        [rg_offset, az_offset] = pix_box4rdr[0:2]
        print 'shiftting the mapping transformation file due to subsetted dataset in'+\
              ' radar coord: rg/az - %d/%d' % (rg_offset, az_offset)
        rg, az, trans_atr = readfile.read_float32(inps.trans_file)
        rg -= rg_offset
        az -= az_offset
        print 'writing >>> '+inps.trans_file
        writefile.write(rg, az, trans_atr, inps.trans_file)
        #writefile.write_roipac_rsc(trans_atr, inps.trans_file+'.rsc')

    # Remove subset_x0/y0/x1/y1
    for File in [inps.ifgram_file, inps.coherence_file, inps.trop_file,\
                 inps.dem_radar_file, inps.dem_geo_file, inps.trans_file]:
        if File:
            atrCmd = 'add_attribute.py '+File+' subset_x0=None subset_x1=None subset_y0=None subset_y1=None'
            print atrCmd
            os.system(atrCmd)

    return inps


def create_subset_dataset(inps, pix_box=None, geo_box=None):
    '''Create/prepare subset of datasets in different folder for time series analysis.
    For dataset (unwrapped interferograms) in radar coord, only support subset in row/col or y/x
    For dataset (unwrapped interferograms) in geo coord, lalo has higher priority than yx, if both are specified.
    '''
    atr = readfile.read_attribute(inps.ifgram_file)
    if not 'X_FIRST' in atr.keys():
        print 'Loaded dataset is in radar coordinate.'
        trans_file_orig = inps.trans_file

        # subset.lalo has higher priority than subset.yx, except when no geomap*.trans exists.
        # don't subset if only subset.lalo without geomap*.trans exsits, because glob2radar won't be accurate.
        if geo_box:
            if inps.trans_file:
                pix_box = None
            else:
                geo_box = None
                print 'turn off subset.lalo because no geomap*.trans file exsits.'
        if not geo_box and not pix_box:
            sys.exit('ERROR: no valid subset input!')

        # Calculate subset range in lat/lon for geo files and y/x for radar files
        if geo_box:
            print 'use subset input in lat/lon'
            print 'calculate corresponding bounding box in radar coordinate.'
            pix_box = subset.bbox_geo2radar(geo_box, atr, trans_file_orig)
        else:
            print 'use subset input in y/x'
            print 'calculate corresponding bounding box in geo coordinate.'
            pix_box = subset.check_box_within_data_coverage(pix_box, atr)
            geo_box = subset.bbox_radar2geo(pix_box, atr, trans_file_orig)

        # subset
        inps = subset_dataset(inps, geo_box, pix_box)

    else:
        print 'Loaded dataset is in geo coordinate.'
        
    return inps


def multilook_dataset(inps, lks_y=None, lks_x=None):
    '''Create a multilooked dataset'''
    return inps


##########################################################################
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
               PySAR v1.2, Aug 2017
 Geodesy Lab, University of Miami, Maimi FL, USA
_________________________________________________
'''
#generate_from: http://patorjk.com/software/taag/

TEMPLATE='''##------------------------ pysarApp_template.txt ------------------------##
## 1. Load Data (--load to exit after this step)
## recommended input files for data in radar coordinates:
##     pysar.unwrapFiles         = 'path of all unwrapped interferograms'
##     pysar.corFiles            = 'path of all coherence files'
##     pysar.demFile.geoCoord    = 'path of DEM in geo   coordinates'
##     pysar.demFile.radarCoord  = 'path of DEM in radar coordinates'
##     pysar.transFile           = 'path of mapping transformation file'
## recommended input files for data in geo coordinates:
##     pysar.unwrapFiles 
##     pysar.corFiles    
##     pysar.dem.geoCoord
## auto - automatic path pattern for Univ of Miami file structure, which are:
##     pysar.unwrapFiles        = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*.unw
##     pysar.corFiles           = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*rlks.cor
##     pysar.demFile.geoCoord   = $SCRATCHDIR/$PROJECT_NAME/DEM/*.dem
##     pysar.demFile.radarCoord = $SCRATCHDIR/$PROJECT_NAME/DONE/*master_date12*/radar*.hgt
##     pysar.transFile          = $SCRATCHDIR/$PROJECT_NAME/GEO/*master_date12*/geomap*.trans
pysar.insarProcessor     = auto  #[roipac, gamma, isce, doris], auto for roipac, InSAR processor
pysar.unwrapFiles        = auto  #[filt*.unw, diff_*.unw], path of all unwrapped interferograms
pysar.corFiles           = auto  #[filt*.cor, filt_*.cor], path of all coherence files
pysar.transFile          = auto  #[geomap*.trans, sim*.UTM_TO_RDC], path of mapping transformation file
pysar.demFile.radarCoord = auto  #[radar*.hgt, sim*.hgt_sim], path of DEM in radar coordinate
pysar.demFile.geoCoord   = auto  #[*.dem, sim*.utm.dem],      path of DEM in geo   coordinate


## 1.1 Subset (optional, --subset to exit after this step)
## if both yx and lalo are specified, use lalo option
pysar.subset.tightBox = auto    #[yes / no], auto for yes
pysar.subset.yx       = auto    #[1800:2000,700:800 / no], auto for no - use the whole area
pysar.subset.lalo     = auto    #[31.5:32.5,130.5:131.0 / no], auto for no - use the whole area
pysar.multilook.yx    = auto    #[4,4 / no], auto for no [not implemented yet]


## 2. Modify Network (optional)
## Coherence-based network modification = MST + Threshold
## 1) calculate a average coherence value for each interferogram using spatial coherence and input mask (with AOI)
## 2) find a minimum spanning tree network with inverse of average coherence as weight
## 3) for all interferograms except for MST's, exclude those with average coherence < minCoherence.
pysar.network.coherenceBased  = auto  #[yes / no], auto for yes
pysar.network.coherenceFile   = auto  #[filename], auto for coherence.h5
pysar.network.minCoherence    = auto  #[0.0-1.0], auto for 0.7
pysar.network.maskFile        = auto  #[file name, no], auto for mask.h5, no for all pixels
pysar.network.maskAoi.yx      = auto  #[y0:y1,x0:x1 / no], auto for no, area of interest for coherence calculation
pysar.network.maskAoi.lalo    = auto  #[lat0:lat1,lon0:lon1 / no], auto for no - use the whole area

pysar.network.tempBaseMax     = auto  #[1-inf, no], auto for no, maximum temporal baseline in days
pysar.network.perpBaseMax     = auto  #[1-inf, no], auto for no, maximum perpendicular spatial baseline in meter
pysar.network.referenceFile   = auto  #[date12_list.txt / Modified_unwrapIfgram.h5 / no], auto for no
pysar.network.excludeDate     = auto  #[20080520,20090817 / no], auto for no
pysar.network.excludeIfgIndex = auto  #[1:5,25 / no], auto for no, list of interferogram number starting from 1


## 3. Reference in Space
## reference all interferograms to one common point in space
## auto - randomly select a pixel with coherence > minCoherence
pysar.reference.yx            = auto   #[257,151 / auto]
pysar.reference.lalo          = auto   #[31.8,130.8 / auto]

pysar.reference.coherenceFile = auto   #[file name], auto for averageSpatialCoherence.h5
pysar.reference.minCoherence  = auto   #[0.0-1.0], auto for 0.85, minimum coherence for auto method
pysar.reference.maskFile      = auto   #[file name / no], auto for mask.h5


## 4. Unwrapping Error Correction
## unwrapping error correction based on the following two methods:
## a. phase closure (Fattahi, 2015, PhD Thesis)
## b. connecting bridge
pysar.unwrapError  = auto   #[yes / no], auto for no


## 5. Network Inversion
## invert network of interferograms into time series
## if network are not fully connected (multiple subsets), Singular-Value Decomposition (SVD) is applied.


## 5.1 Temporal Coherence
## calculate temporal coherence based on Tizzani et al., 2007 (IEEE-TGRS)
## and generate a mask file with temporal coherence > minCoherence
pysar.temporalCoherence.threshold  = auto    #[0.0-1.0], auto for 0.7


## 6. Local Oscillator Drift (LOD) Correction (for Envisat only, no need to setup, it runs automatically)
## correct LOD if input dataset comes from Envisat and in radar coordinate
## skip this step for all the other satellites.


## 7. Tropospheric Delay Correction
## correct tropospheric delay using the following methods:
## a. pyaps - use weather re-analysis data (Jolivet et al., 2011, GRL, need to install PyAPS)
## b. height_correlation - correct stratified tropospheric delay (Doin et al., 2009, J Applied Geop)
## c. base_trop_cor - (not recommend) baseline error and stratified tropo simultaneously (Jo et al., 2010, Geo J)
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.polyOrder    = auto  #[1 / 2 / 3], auto for 1, for height_correlation method
pysar.troposphericDelay.weatherModel = auto  #[ERA / MERRA / NARR], auto for ECMWF, for pyaps method


## 8. Topographic (DEM) Residual Correction (Fattahi and Amelung, 2013, IEEE-TGRS)
pysar.topoError              = auto    #[yes / no], auto for yes
pysar.topoError.polyOrder    = auto    #[1 / 2 / 3], auto for 2, polynomial order of temporal deformation model
pysar.topoError.excludeDate  = auto    #[20101120 / txtFile / no], auto for no, date not used for error estimation
pysar.topoError.stepFuncDate = auto    #[20080529 / no], auto for no, date of step jump, i.e. eruption/earthquake date


## 8.1 Phase Residual Root Mean Square
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
## Recommendation: quadratic for whole image, plane for local/small area
pysar.residualRms.maskFile        = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)
pysar.residualRms.saveRefDate     = auto  #[yes / no], auto for yes, save date with min RMS to txt/pdf file.
pysar.residualRms.saveExcludeDate = auto  #[yes / no], auto for yes, save date(s) with RMS > minStd to txt/pdf file.


## 9. Reference in Time
## reference all timeseries to one date in time
## auto - choose date with minimum residual RMS using value from step 8.1
## no   - do not change reference date, keep the defaut one (1st date usually) and skip this step
pysar.reference.date = auto   #[auto / reference_date.txt / 20090214 / no]


## 10. Phase Ramp Removal (optional)
## remove phase ramp for each epoch, useful to check localized deformation, i.e. volcanic, land subsidence, etc.
## [plane, quadratic, plane_range, plane_azimuth, quadratic_range, quadratic_azimuth, baseline_cor, base_trop_cor]
pysar.deramp          = auto  #[no / plane / quadratic], auto for no - no ramp will be removed
pysar.deramp.maskFile = auto  #[file name / no], auto for maskTempCoh.h5, mask file for ramp estimation


## 11. Velocity Inversion
## estimate linear velocity from timeseries, and from tropospheric delay file if exists.
pysar.velocity.excludeDate = auto   #[exclude_date.txt / 20080520,20090817 / no], auto for exclude_date.txt
pysar.velocity.startDate   = auto   #[20070101 / no], auto for no
pysar.velocity.endDate     = auto   #[20101230 / no], auto for no


## 12. Post-processing (geocode, output to Google Earth, UNAVCO, etc.)
pysar.geocode      = auto   #[yes / no], auto for yes
pysar.save.kml     = auto   #[yes / no], auto for yes, save geocoded velocity to Google Earth KMZ file
pysar.save.geotiff = auto   #[yes / no], auto for no, save geocoded velocity to Geotiff format [not implemented yet]
pysar.save.unavco  = auto   #[yes / no], auto for no, save timeseries to UNAVCO InSAR Archive format
pysar.plot         = auto   #[yes / no], auto for yes, plot files generated by pysarApp default processing to PIC folder
'''

EXAMPLE='''example:
  pysarApp.py
  pysarApp.py  SanAndreasT356EnvD.template
  pysarApp.py  SanAndreasT356EnvD.template  --load-data
  pysarApp.py  SanAndreasT356EnvD.template  --subset-data
  pysarApp.py  SanAndreasT356EnvD.template  --modify-network
  pysarApp.py  SanAndreasT356EnvD.template  --dir ~/insarlab/SanAndreasT356EnvD/PYSAR

  # Generate template file:
  pysarApp.py -g
  pysarApp.py SanAndreasT356EnvD.template -g

  # Reset and re-run right after loading data
  pysarApp.py --reset
  pysarApp.py

  --------------------------------------------
  Open pysar_template.txt file for details.
  --------------------------------------------
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


def cmdLineParse():
    parser = argparse.ArgumentParser(description=LOGO,
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
                                     #epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('-v','--version', action='version', version='%(prog)s 1.2')
    parser.add_argument('custom_template_file', nargs='?',\
                        help='custom template with option settings.\n'+\
                             "It's equivalent to None, if pysarApp_template.txt is input, as it will be read always.")
    parser.add_argument('--dir', dest='work_dir',\
                        help='PySAR working directory, default is:\n'+\
                             'a) current directory, or\n'+\
                             'b) $SCRATCHDIR/projectName/PYSAR, if meets the following 3 requirements:\n'+\
                             '     1) miami_path = True in pysar/__init__.py\n'+\
                             '     2) environmental variable $SCRATCHDIR exists\n'+\
                             '     3) input custom template with basename same as projectName\n')
    parser.add_argument('-g', dest='generate_template', action='store_true',\
                        help='Generate default template (and merge with custom template), then exit.')
    parser.add_argument('--reset', action='store_true',\
                        help='Reset files attributes to re-run pysarApp.py after loading data\n'+\
                             '    by removing ref_y/x/lat/lon and set drop_ifgram=no for unwrapIfgram.h5 and coherence.h5')
    parser.add_argument('--load-data', dest='load_dataset', action='store_true',\
                        help='Step 1. Load/check dataset, then exit')
    parser.add_argument('--subset-data', dest='subset_dataset', action='store_true',\
                        help='Step 1.1 Subset the whole dataset with setting in template, then exit')
    parser.add_argument('--modify-network', dest='modify_network', action='store_true',\
                        help='Step 2. Modify the network, then exit')

    inps = parser.parse_args()
    if inps.custom_template_file and os.path.basename(inps.custom_template_file) == 'pysarApp_template.txt':
        inps.custom_template_file = None
    return inps


##########################################################################
def main(argv):
    start = time.time()
    inps = cmdLineParse()


    #########################################
    # Initiation
    #########################################
    print LOGO

    # Project Name
    inps.project_name = None
    if inps.custom_template_file:
        inps.custom_template_file = os.path.abspath(inps.custom_template_file)
        inps.project_name = os.path.splitext(os.path.basename(inps.custom_template_file))[0]
        print 'Project name: '+inps.project_name

    # Work directory
    if not inps.work_dir:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ and inps.project_name:
            inps.work_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PYSAR'
        else:
            inps.work_dir = os.getcwd()
    inps.work_dir = os.path.abspath(inps.work_dir)

    if not os.path.isdir(inps.work_dir):
        os.makedirs(inps.work_dir)
    os.chdir(inps.work_dir)
    print "Go to work directory: "+inps.work_dir


    #####for Univ of Miami
    # Copy bl_list.txt file from PROCESS directory
    try:
        process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
        if ut.update_file('bl_list.txt', process_dir+'/bl_list.txt', check_readable=False):
            shutil.copy2(process_dir+'/bl_list.txt', inps.work_dir)
    except: pass

    # Copy UNAVCO attribute txt file
    file_list = ['unavco_attributes.txt']
    try:
        process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
        if ut.update_file(file_list[0], process_dir+'/'+file_list[0], check_readable=False):
            shutil.copy2(process_dir+'/'+file_list[0], inps.work_dir)
    except: pass
    try:
        inps.unavco_atr_file = ut.get_file_list(file_list, abspath=True)[0]
    except:
        inps.unavco_atr_file = None
        if pysar.miami_path and 'SCRATCHDIR' in os.environ and inps.project_name:
            print 'No UNAVCO attributes file found in PROCESS directory, skip copy'


    #########################################
    # Read Template Options
    #########################################
    print '\n*************** Template Options ****************'
    # default template
    inps.template_file = 'pysarApp_template.txt'
    if not os.path.isfile(inps.template_file):
        print 'generate default template file: '+inps.template_file
        f = open(inps.template_file, 'w')
        f.write(TEMPLATE)
        f.close()
    else:
        print 'default template file exists: '+inps.template_file

    # custom template
    if inps.custom_template_file:
        # Copy custom template file to work directory
        if ut.update_file(os.path.basename(inps.custom_template_file), inps.custom_template_file, check_readable=False):
            shutil.copy2(inps.custom_template_file, inps.work_dir)

        # Read custom template
        print 'read custom template file: '+inps.custom_template_file
        custom_template = readfile.read_template(inps.custom_template_file)
        # correct some loose type errors
        for key in custom_template.keys():
            if   custom_template[key].lower() in ['default']:          custom_template[key] = 'auto'
            elif custom_template[key].lower() in ['n','off','false']:  custom_template[key] = 'no'
            elif custom_template[key].lower() in ['y','on','true']:    custom_template[key] = 'yes'
        for key in ['pysar.deramp', 'pysar.troposphericDelay.method']:
            if key in custom_template.keys():
                custom_template[key] = custom_template[key].lower().replace('-','_')

        # Update default template with custom input template
        print 'update default template based on input custom template'
        inps.template_file = ut.update_template_file(inps.template_file, custom_template)

    if inps.generate_template:
        sys.exit('Exit as planed after template file generation.')

    print 'read default template file: '+inps.template_file
    inps.template_file = os.path.abspath(inps.template_file)
    template = readfile.read_template(inps.template_file)

    # Get existing tropo delay file
    inps.trop_model = 'ECMWF'
    key = 'pysar.troposphericDelay.weatherModel'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.trop_model = 'ECMWF'
        else:
            inps.trop_model = value
    # Grab tropospheric delay file
    try:    inps.trop_file = ut.get_file_list(inps.trop_model+'.h5', abspath=True)[0]
    except: inps.trop_file = None


    #########################################
    # Loading Data
    #########################################
    print '\n*************** Load Data ****************'
    loadCmd = 'load_data.py --dir '+inps.work_dir+' --template '+inps.template_file
    if inps.custom_template_file:
        loadCmd += ' '+inps.custom_template_file+' --project '+inps.project_name
    print loadCmd
    os.system(loadCmd)
    os.chdir(inps.work_dir)

    print '--------------------------------------------'
    inps = ut.check_loaded_dataset(inps.work_dir, inps)
    if not inps.ifgram_file:
        print '\nERROR: No interferograms file found!\n'
        sys.exit('Exit.')

    if inps.load_dataset:
        sys.exit('Exit as planed after loading/checking the dataset')

    if inps.reset:
        print 'Reset dataset attributtes for a fresh re-run with options from %s' % os.path.basename(inps.template_file)
        print '-----------------------------------------------------------------------------------'
        # Reset network
        networkCmd = 'modify_network.py '+inps.ifgram_file
        if inps.coherence_file:
            networkCmd +=  ' '+inps.coherence_file
        networkCmd += ' --reset'
        print networkCmd
        os.system(networkCmd)

        # Reset reference pixel
        seedCmd = 'seed_data.py '+inps.ifgram_file+' --reset'
        print seedCmd
        os.system(seedCmd)

        #print 'Exit as planed after reset.'
        #print 'Now run pysarApp.py to re-processing this dataset.'
        #sys.exit()


    #########################################
    # Check the subset (Optional)
    #########################################
    if inps.trans_file and template['pysar.subset.tightBox'] in ['yes','auto']:
        outName = os.path.splitext(inps.trans_file)[0]+'_tight'+os.path.splitext(inps.trans_file)[1]
        # Get bounding box of valid area in geomap*.trans file
        trans_rg, trans_atr = readfile.read(inps.trans_file, (), 'range')
        rg_unique, rg_pos = np.unique(trans_rg, return_inverse=True)
        idx_row, idx_col = np.where(trans_rg != rg_unique[np.bincount(rg_pos).argmax()])
        pix_box = (np.min(idx_col)-10, np.min(idx_row)-10, np.max(idx_col)+10, np.max(idx_row)+10)
        # Subset geomap_file only if it could save > 20% percent of area
        if abs((pix_box[2]-pix_box[0])*(pix_box[3]-pix_box[1])) < 0.8*(trans_rg.shape[0]*trans_rg.shape[1]):
            print "Get tight subset of geomap*.trans file and/or DEM file in geo coord"
            print '--------------------------------------------'
            inps = subset.subset_box2inps(inps, pix_box, None)
            inps.fill_value = 0.0
            inps.trans_file = check_subset_file(inps.trans_file, vars(inps), outName)

            # Subset DEM in geo coord
            if inps.dem_geo_file:
                outName = os.path.splitext(inps.dem_geo_file)[0]+'_tight'+os.path.splitext(inps.dem_geo_file)[1]
                geomap_atr = readfile.read_attribute(inps.trans_file)
                pix_box, geo_box = subset.get_coverage_box(geomap_atr)
                inps = subset.subset_box2inps(inps, pix_box, geo_box)
                inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps), outName, overwrite=True)


    # Subset based on input template
    if not all(template[key] in ['auto', 'no'] for key in ['pysar.subset.yx','pysar.subset.lalo']):
        print '\n*************** Subset ****************'
        # Read subset option from template file, and return None if lalo/yx is not specified.
        pix_box, geo_box = subset.read_subset_template2box(inps.template_file)
        inps = create_subset_dataset(inps, pix_box, geo_box)

    if inps.subset_dataset:
        sys.exit('Exit as planed after subsetting the dataset')


    #########################################
    # Generating Aux files
    #########################################
    print '\n*************** Generate Initial Mask ****************'
    # Initial mask
    inps.mask_file = 'mask.h5'
    if ut.update_file(inps.mask_file, inps.ifgram_file):
        print 'creating mask file using non-zero pixels from file: '+inps.ifgram_file
        inps.mask_file = ut.nonzero_mask(inps.ifgram_file, inps.mask_file)

    # Average spatial coherence
    if inps.coherence_file:
        inps.spatial_coh_file = 'averageSpatialCoherence.h5'
        print 'creating average spatial coherence from file: '+inps.coherence_file
        if ut.update_file(inps.spatial_coh_file, inps.coherence_file):
            inps.spatial_coh_file = ut.temporal_average(inps.coherence_file, inps.spatial_coh_file)
    else:
        inps.spatial_coh_file = None


    #########################################
    # Network Modification (Optional)
    #########################################
    atr = readfile.read_attribute(inps.ifgram_file)
    h5 = h5py.File(inps.ifgram_file, 'r')
    ifgram_list_all = h5[atr['FILE_TYPE']].keys()
    ifgram_list_keep = ut.check_drop_ifgram(h5, atr, ifgram_list_all)
    h5.close()
    ifgram_num_drop = len(ifgram_list_all) - len(ifgram_list_keep)
    if ifgram_num_drop > 0:
        print '\n*************** Modify Network ****************'
        print 'find number of dropped interferograms: %d, skip updating.' % (ifgram_num_drop)
        print 'To re-modify the network, use modify_network.py --reset option to restore all pairs info.'
        msg_str = 'modify_network.py '+inps.ifgram_file
        if inps.coherence_file:
            msg_str +=  ' '+inps.coherence_file
        msg_str += ' --reset'
        print msg_str
    else:
        print '\n*************** Modify Network ****************'
        networkCmd = 'modify_network.py --template '+inps.template_file+' '+inps.ifgram_file
        if inps.coherence_file:
            networkCmd += ' '+inps.coherence_file
        if inps.trans_file:
            networkCmd += ' --trans '+inps.trans_file
        print networkCmd
        os.system(networkCmd)

    # Plot network colored in spatial coherence
    print '--------------------------------------------'
    plotCmd = 'plot_network.py '+inps.ifgram_file+' --coherence '+inps.coherence_file+' --mask '+inps.mask_file+' --nodisplay'
    print plotCmd
    if ut.update_file('Network.pdf', [inps.ifgram_file, inps.coherence_file, inps.mask_file], check_readable=False):
        os.system(plotCmd)

    if inps.modify_network:
        sys.exit('Exit as planed after network modification.')


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print '\n**********  Reference in space  ***************'
    seedCmd = 'seed_data.py '+inps.ifgram_file+' --template '+inps.template_file+' --mark-attribute'
    if inps.trans_file:
        seedCmd += ' --trans '+inps.trans_file
    run_seedCmd = True

    ## Skip calling seed command only if 1) ref_y/x exists AND 2) pysar.reference.yx/lalo == auto
    atr = readfile.read_attribute(inps.ifgram_file)
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print 'Find reference pixel info from %s in y/x: [%d, %d]' % (os.path.basename(inps.ifgram_file), ref_y, ref_x)
        print '    To remove reference pixel info, use seed_data.py --reset option:'
        print '    seed_data.py '+inps.ifgram_file+' --reset'
        prefix = 'pysar.reference.'
        if all(template[prefix+i] in ['auto','no'] for i in ['yx','lalo']):
            run_seedCmd = False
            print 'No specific coordinates input found, no need to re-select reference pixel'
        else:
            print 'Specific coordinates input found, re-select reference pixel with options from template file'
    except: pass

    if run_seedCmd:
        print 'Call seed_data.py to find reference pixel in space'
        print seedCmd
        os.system(seedCmd)


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    if template['pysar.unwrapError'] not in ['auto','no']:
        print '\n**********  Unwrapping Error Correction  **************'
        outName = os.path.splitext(inps.ifgram_file)[0]+'_unwCor.h5'
        unwCmd='unwrap_error.py '+inps.ifgram_file+' '+inps.mask_file
        print unwCmd
        if ut.update_file(outName, inps.ifgram_file):
            print 'This might take a while depending on the size of your data set!'
            os.system(unwCmd)
        inps.ifgram_file = outName


    #########################################
    # Inversion of Interferograms
    ########################################
    print '\n**********  Network Inversion to Time Series  ********************'
    inps.timeseries_file = 'timeseries.h5'
    invertCmd = 'ifgram_inversion.py '+inps.ifgram_file
    print invertCmd
    if ut.update_file(inps.timeseries_file, inps.ifgram_file):
        os.system(invertCmd)

    ## Check DEM file for tropospheric delay setting
    ## DEM is needed with same coord (radar/geo) as timeseries file
    atr = readfile.read_attribute(inps.timeseries_file)
    if 'X_FIRST' in atr.keys():
        demFile = inps.dem_geo_file
    else:
        demFile = inps.dem_radar_file

    if ut.update_file(demFile):
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print 'ERROR:'
        print '    DEM file was not found!'
        if 'pysar.troposphericDelay.method' in template.keys():
            print '    Continue without tropospheric correction ...'
            template['pysar.troposphericDelay.method'] = 'no'
        if template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']:
            template['pysar.deramp'] = 'no'
        print '++++++++++++++++++++++++++++++++++++++++++++++'


    ##############################################
    # Temporal Coherence: 
    #   A parameter to evaluate the consistency 
    #   of timeseries with the interferograms
    ##############################################
    print '\n********** Temporal Coherence file  *********'
    inps.temp_coh_file = 'temporalCoherence.h5'
    tempCohCmd = 'temporal_coherence.py '+inps.ifgram_file+' '+inps.timeseries_file+' '+inps.temp_coh_file
    print tempCohCmd
    if ut.update_file(inps.temp_coh_file, inps.timeseries_file):
        os.system(tempCohCmd)

    print '\n--------------------------------------------'
    print 'Update Mask based on Temporal Coherence ...'
    # Read template option
    inps.min_temp_coh = 0.7
    key = 'pysar.temporalCoherence.threshold'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.min_temp_coh = 0.7
        else:
            inps.min_temp_coh = float(value)
    outName = 'maskTempCoh.h5'
    maskCmd = 'generate_mask.py '+inps.temp_coh_file+' -m '+str(inps.min_temp_coh)+' -o '+outName
    print maskCmd
    if ut.update_file(outName, inps.temp_coh_file):
        os.system(maskCmd)
    inps.mask_file = outName


    ###############################################
    ## Incident Angle
    ###############################################
    #print '\n********** Incident Angle file  *************'
    #inps.inc_angle_file = 'incidenceAngle.h5'
    #incAngleCmd = 'incidence_angle.py '+inps.timeseries_file+' '+inps.inc_angle_file
    #print incAngleCmd
    #if ut.update_file(inps.inc_angle_file, inps.timeseries_file):
    #    os.system(incAngleCmd)


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   for Envisat data in radar coord only
    ############################################## 
    sar_mission = atr['PLATFORM'].lower()
    if sar_mission.startswith('env'):
        print '\n**********  Local Oscillator Drift correction for Envisat  ********'
        if 'Y_FIRST' not in atr.keys():
            outName = os.path.splitext(inps.timeseries_file)[0]+'_LODcor.h5'
            lodCmd = 'lod.py '+inps.timeseries_file
            print lodCmd
            if ut.update_file(outName, inps.timeseries_file):
                os.system(lodCmd)
            inps.timeseries_file = outName
        else:
            warnings.warn('Can not apply LOD correction for file in radar coord. Skip it for now.')


    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print '\n**********  Tropospheric Delay Correction  ******************'
    key = 'pysar.troposphericDelay.method'
    if (key in template.keys() and template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']):
        message='''
        Orbital error correction was BaseTropCor.
        Tropospheric correction was already applied simultaneous with baseline error correction.
        Tropospheric correction can not be applied again.
        To apply the tropospheric correction separated from baseline error correction, \
           choose other existing options for orbital error correction.
        '''
        warnings.warn(message)
        template[key] = 'no'

    # read template option
    inps.trop_method = 'pyaps'
    key = 'pysar.troposphericDelay.method'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.trop_method = 'pyaps'
        else:
            inps.trop_method = value

    inps.trop_poly_order = '1'
    key = 'pysar.troposphericDelay.polyOrder'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.trop_poly_order = '1'
        else:
            inps.trop_poly_order = value

    # Call scripts
    if inps.trop_method == 'height_correlation':
        print 'tropospheric delay correction with height-correlation approach'
        tropCmd = 'tropcor_phase_elevation.py '+inps.timeseries_file+' -d '+demFile+\
                  ' -p '+inps.trop_poly_order+' -m '+inps.mask_file
        print tropCmd
        outName = os.path.splitext(inps.timeseries_file)[0]+'_tropHgt.h5'
        if ut.update_file(outName, inps.timeseries_file):
            os.system(tropCmd)
        inps.timeseries_file = outName

    elif inps.trop_method == 'pyaps':
        print 'Atmospheric correction using Weather Re-analysis dataset (using PyAPS software)'
        print 'Weather Re-analysis dataset: '+inps.trop_model
        tropCmd = 'tropcor_pyaps.py '+inps.timeseries_file+' -d '+demFile+' -s '+inps.trop_model+\
                  ' --weather-dir '+inps.work_dir+'/../WEATHER'
        print tropCmd
        outName = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.trop_model+'.h5'
        if ut.update_file(outName, inps.timeseries_file):
            try:
                inps.trop_file = ut.get_file_list(inps.trop_model+'.h5')[0]
                diffCmd = 'diff.py '+inps.timeseries_file+' '+inps.trop_file+' -o '+outName
                print 'Use existed tropospheric delay file: '+inps.trop_file
                print diffCmd
                os.system(diffCmd)
            except:
                os.system(tropCmd)
        inps.timeseries_file = outName

    else:
        print 'No atmospheric delay correction.'

    # Grab tropospheric delay file
    try:    inps.trop_file = ut.get_file_list(inps.trop_model+'.h5')[0]
    except: inps.trop_file = None


    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print '\n**********  Topographic Residual (DEM error) correction  *******'
    outName = os.path.splitext(inps.timeseries_file)[0]+'_demErr.h5'
    topoCmd = 'dem_error.py '+inps.timeseries_file+' -o '+outName+' --template '+inps.template_file
    print topoCmd
    inps.timeseries_resid_file = None
    if template['pysar.topoError'] in ['yes','auto']:
        print 'Correcting topographic residuals using method from Fattahi and Amelung, 2013, TGRS ...'
        if ut.update_file(outName, inps.timeseries_file):
            os.system(topoCmd)
        inps.timeseries_file = outName
        inps.timeseries_resid_file = os.path.splitext(outName)[0]+'InvResid.h5'
    else:
        print 'No correction for topographic residuals.'


    ##############################################
    # Timeseries Residual Standard Deviation
    ##############################################
    print '\n**********  Timeseries Residual Root Mean Square  *******'
    if inps.timeseries_resid_file:
        rmsCmd = 'timeseries_rms.py '+inps.timeseries_resid_file+' --template '+inps.template_file
        print rmsCmd
        os.system(rmsCmd)
    else:
        print 'No timeseries residual file found! Skip residual RMS analysis.'


    ##############################################
    # Reference in Time
    ##############################################
    print '\n**********  Reference in Time  *******'
    if template['pysar.reference.date'] != 'no':
        outName = os.path.splitext(inps.timeseries_file)[0]+'_refDate.h5'
        refCmd = 'reference_epoch.py '+inps.timeseries_file+' --template '+inps.template_file
        print refCmd

        if ut.update_file(outName, inps.timeseries_file):
            os.system(refCmd)

        if not ut.update_file(outName):
            inps.timeseries_file = outName
    else:
        print 'No reference change in time.'


    ##############################################
    # Phase Ramp Correction (Optional)
    ##############################################
    # Read template option
    inps.deramp_mask_file = 'maskTempCoh.h5'
    key = 'pysar.deramp.maskFile'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.deramp_mask_file = 'maskTempCoh.h5'
        else:
            inps.deramp_mask_file = value

    print '\n**********  Ramp Removal  ***********************'
    if template['pysar.deramp'] not in ['no','auto']:
        inps.deramp_method = template['pysar.deramp']
        print 'Phase Ramp Removal method : '+inps.deramp_method

        if inps.deramp_method in ['plane', 'quadratic', 'plane_range', 'quadratic_range',\
                                  'plane_azimuth', 'quadratic_azimuth']:
            derampCmd = 'remove_plane.py '+inps.timeseries_file+' -s '+inps.deramp_method+' -m '+inps.mask_file
            print derampCmd

            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.deramp_method+'.h5'
            if ut.update_file(outName, inps.timeseries_file):
                os.system(derampCmd)
            inps.timeseries_file = outName

        elif inps.deramp_method in ['baseline_cor','baselinecor']:
            if not 'X_FIRST' in atr.keys():
                derampCmd = 'baseline_error.py '+inps.timeseries_file+' '+inps.mask_file
                print derampCmd

                outName = os.path.splitext(inps.timeseries_file)[0]+'_baselineCor.h5'
                if ut.update_file(outName, inps.timeseries_file):
                    os.system(derampCmd)
                inps.timeseries_file = outName
            else:
                warnings.warn('BaselineCor method can only be applied in radar coordinate, skipping correction')

        elif inps.deramp_method in ['base_trop_cor','basetropcor','baselinetropcor']:
            if not 'X_FIRST' in atr.keys():
                print 'Joint estimation of Baseline error and tropospheric delay [height-correlation approach]'
                try:    poly_order = template['pysar.troposphericDelay.polyOrder']
                except: poly_order = '1'
                derampCmd = 'baseline_trop.py '+inps.timeseries_file+' '+inps.dem_radar_file+' '+\
                            poly_order+' range_and_azimuth'
                print derampCmd

                outName = os.path.splitext(inps.timeseries_file)[0]+'_baseTropCor.h5'
                if ut.update_file(outName, inps.timeseries_file):
                    os.system(derampCmd)
                inps.timeseries_file = outName
            else:
                warnings.warn('BaselineCor method can only be applied in radar coordinate, skipping correction')
        else:
            warnings.warn('Unrecognized phase ramp method: '+template['pysar.deramp'])
    else:
        print 'No phaes ramp removal.'


    #############################################
    # Velocity and rmse maps
    #############################################
    print '\n**********  Velocity estimation  **********************'
    inps.vel_file = 'velocity.h5'
    velCmd = 'timeseries2velocity.py '+inps.timeseries_file+' --template '+inps.template_file+' -o '+inps.vel_file
    print velCmd
    if ut.update_file(inps.vel_file, [inps.timeseries_file, inps.template_file]):
        os.system(velCmd)

    # Velocity from Tropospheric delay
    if inps.trop_file:
        suffix = os.path.splitext(os.path.basename(inps.trop_file))[0]
        suffix = suffix[0].upper()+suffix[1:].lower()
        inps.trop_vel_file = 'velocity'+suffix+'.h5'
        velCmd = 'timeseries2velocity.py '+inps.trop_file+' --template '+inps.template_file+' -o '+inps.trop_vel_file
        print velCmd
        if ut.update_file(inps.trop_vel_file, [inps.trop_file, inps.template_file]):
            os.system(velCmd)


    ############################################
    # Post-processing
    # Geocodeing, masking and save to KML 
    ############################################
    print '\n**********  Post-processing  ********************************'
    inps.geo_vel_file = None
    inps.geo_temp_coh_file = None
    inps.geo_timeseries_file = None

    # Check geocoding requirement
    key = 'pysar.geocode'
    if template[key] in ['auto','yes']:
        if 'Y_FIRST' in atr.keys():
            print 'dataset is in geo coordinate, no need to geocode.'
            template[key] = 'no'
            inps.geo_vel_file        = inps.vel_file
            inps.geo_temp_coh_file   = inps.temp_coh_file
            inps.geo_timeseries_file = inps.timeseries_file
        #elif not ut.which('geocode.pl'):
        #    print 'Can not find executable geocode.pl from ROI_PAC, skip geocoding.'
        #    template[key] = 'no'

    # Geocoding
    if template[key] in ['yes','auto']: 
        print '\n--------------------------------------------'
        inps.geo_vel_file        = check_geocode_file(inps.trans_file, inps.vel_file)
        inps.geo_temp_coh_file   = check_geocode_file(inps.trans_file, inps.temp_coh_file)
        inps.goe_timeseries_file = check_geocode_file(inps.trans_file, inps.timeseries_file)

    # Mask in geo coord
    if inps.geo_temp_coh_file:
        # Generate mask in geo coord
        print '\n--------------------------------------------'
        outName = 'geo_maskTempCoh.h5'
        maskCmd = 'generate_mask.py '+inps.geo_temp_coh_file+' -m '+str(inps.min_temp_coh)+' -o '+outName
        print maskCmd
        if ut.update_file(outName, inps.geo_temp_coh_file):
            os.system(maskCmd)
        inps.geo_mask_file = outName

        # Mask geo_velocity file
        if inps.geo_vel_file and inps.geo_mask_file:
            outName = os.path.splitext(inps.geo_vel_file)[0]+'_masked.h5'
            maskCmd = 'mask.py '+inps.geo_vel_file+' -m '+inps.geo_mask_file+' -o '+outName
            print maskCmd
            if ut.update_file(outName, [inps.geo_vel_file, inps.geo_mask_file]):
                os.system(maskCmd)
            try:  inps.geo_vel_file = glob.glob(outName)[0]
            except:  pass

    # Save to Google Earth KML file
    if inps.geo_vel_file and template['pysar.save.kml'] in ['auto','yes']:
        print '\n--------------------------------------------'
        print 'creating Google Earth KMZ file for geocoded velocity file: '+inps.geo_vel_file+' ...'
        outName = os.path.splitext(inps.geo_vel_file)[0]+'.kmz'
        kmlCmd = 'save_kml.py '+inps.geo_vel_file
        print kmlCmd
        if ut.update_file(outName, inps.geo_vel_file, check_readable=False):
            os.system(kmlCmd)


    #############################################
    # Save to UNAVCO InSAR Archive format
    #############################################
    if template['pysar.save.unavco'] in ['yes']:
        print '\n*********  Output to UNAVCO InSAR Archive Format  ***********'
        if not inps.trans_file or not inps.dem_geo_file:
            warnings.warn('No geomap*.tran file or DEM in geo coord found! Skip saving.')
        else:
            # 1. Time series file
            inps.geo_timeseries_file = check_geocode_file(inps.trans_file, inps.timeseries_file)
            # Add UNAVCO attributes
            if inps.unavco_atr_file:
                atrCmd = 'add_attribute.py '+inps.geo_timeseries_file+' '+inps.unavco_atr_file
                print atrCmd
                os.system(atrCmd)

            # 2. Temporal Coherence
            inps.geo_temp_coh_file = check_geocode_file(inps.trans_file, inps.temp_coh_file)

            # 3. Mask file
            inps.geo_mask_file = 'geo_maskTempCoh.h5'
            if inps.geo_temp_coh_file and ut.update_file(inps.geo_mask_file, inps.geo_temp_coh_file):
                maskCmd = 'generate_mask.py '+inps.geo_temp_coh_file+' -m 0.7 -o '+inps.geo_mask_file
                print maskCmd
                os.system(maskCmd)

            # 4. Incidence Angle
            inps.inc_angle_file = 'incidenceAngle.h5'
            if ut.update_file(inps.inc_angle_file, inps.timeseries_file):
                incAngleCmd = 'incidence_angle.py '+inps.timeseries_file+' '+inps.inc_angle_file
                print incAngleCmd
                os.system(incAngleCmd)
            inps.geo_inc_angle_file = check_geocode_file(inps.trans_file, inps.inc_angle_file)

            # Save to UNAVCO format
            print '--------------------------------------------'
            SAT = unavco.get_mission_name(atr)
            try:
                inps.unavco_file = ut.get_file_list(SAT+'_*.he5')[0]
                print 'Find existed UNAVCO time-series file: '+inps.unavco_file
            except:
                inps.unavco_file = None
                print 'No UNAVCO time-series file exists yet.'
            #inps.unavco_file = unavco.get_unavco_filename(inps.geo_timeseries_file)
            unavcoCmd = 'save_unavco.py '+inps.geo_timeseries_file+' -d '+inps.dem_geo_file+\
                        ' -i '+inps.geo_inc_angle_file+' -c '+inps.geo_temp_coh_file+' -m '+inps.geo_mask_file
            print unavcoCmd
            if ut.update_file(inps.unavco_file, [inps.geo_timeseries_file, inps.geo_temp_coh_file, inps.geo_mask_file,\
                                                 inps.geo_inc_angle_file, inps.dem_geo_file], check_readable=False):
                os.system(unavcoCmd)


    #############################################
    # Plot Figures
    #############################################
    if template['pysar.plot'] in ['yes','auto']:
        print '\n*********  Plot and Save pysarApp runing results to PIC  ***********'
        inps.plot_sh_file = 'plot_pysarApp.sh'

        # Copy to workding directory if not existed yet.
        if not os.path.isfile(inps.work_dir+'/'+inps.plot_sh_file):
            print 'copy $PYSAR_HOME/shellscripts/'+inps.plot_sh_file+' to working directory'
            try:
                shutil.copy2(ut.which(inps.plot_sh_file), inps.work_dir)
            except:
                print 'WARNING: no '+inps.plot_sh_file+' found in the environment variable path, skip plotting.'

        print 'for better performance, edit the input parameters in '+inps.plot_sh_file+' and re-run this script.'
        plotCmd = 'chmod +x '+inps.plot_sh_file;    print plotCmd;     os.system(plotCmd)
        plotCmd = './'+inps.plot_sh_file;           print plotCmd;     os.system(plotCmd)


    #############################################
    #                PySAR v1.2                 #
    #############################################
    s = time.time()-start;  m, s = divmod(s, 60);  h, m = divmod(m, 60)
    print '\nTime used: %02d hours %02d mins %02d secs' % (h, m, s)
    print '\n###############################################'
    print 'End of PySAR processing!'
    print '################################################\n'



###########################################################################################
if __name__ == '__main__':
    main(sys.argv[:])

