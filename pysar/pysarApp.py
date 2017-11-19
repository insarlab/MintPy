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
import subprocess

import h5py
import numpy as np

import pysar
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar.subset as subset
import pysar.save_unavco as unavco


###############################################################################
def check_geocode_file(geomapFile, File, templateFile=None, outFile=None):
    '''Geocode input file or use existed geocoded file.'''
    if not File:
        return None
    else:
        atr = readfile.read_attribute(File)
        if 'Y_FIRST' in atr.keys():
            return File
    
    if not geomapFile:
        warnings.warn('No lookup file found! Skip geocoding.')
        return None

    if not outFile:  outFile = 'geo_'+os.path.basename(File)

    if ut.update_file(outFile, File):
        geocodeCmd = 'geocode.py '+File+' -l '+os.path.basename(geomapFile)
        if templateFile:
            geocodeCmd += ' -t '+templateFile
        print geocodeCmd
        status = subprocess.Popen(geocodeCmd, shell=True).wait()

    try:    outFile = glob.glob(outFile)[0]
    except: outFile = None
    return outFile


def check_subset_file(File, inps_dict, outFile=None, overwrite=False):
    '''Subset input file or use existed subseted file.'''
    if not File:
        return None

    if not outFile:
        if os.getcwd() == os.path.dirname(os.path.abspath(File)):
            outFile = 'subset_'+os.path.basename(File)
        else:
            outFile = os.path.basename(File)

    if ut.update_file(outFile, File, overwrite):
        outFile = subset.subset_file(File, inps_dict, outFile)
    return outFile


def subset_dataset(inps, template_file):
    '''Create/prepare subset of datasets in different folder for time series analysis.
    1) Read subset info from lat/lon or y/x, and convert into y/x
       where lat/lon > y/x in priority unless a) no lookup file AND b) dataset is in radar coord
       While converting lalo to yx, yx should be the bounding box of lalo.
    2) for geo-coord dataset, use y/x from 1) to subset all the files
       for radar-coord dataset, use y/x from 1) to subset all radar-coord files; then get y/x bounding box
          in lat/lon and use it to subset all geo-coord files.
    '''
    pix_box, geo_box = subset.read_subset_template2box(template_file)
    atr = readfile.read_attribute(inps.ifgram_file)

    #####Step 1: Read subset info into pix_box
    #Check conflict
    if geo_box:
        if not inps.lookup_file and 'Y_FIRST' not in atr.keys():
            print 'WARNING: turn off pysar.subset.lalo because:'
            print '    1) no lookup file AND '
            print '    2) dataset is in radar coord'
            geo_box = None
    if not geo_box and not pix_box:
        sys.exit('ERROR: no valid subset input!\nTry to set pysar.subset.yx/lalo all to auto/no')

    #geo_box --> pix_box
    if geo_box:
        if 'Y_FIRST' not in atr.keys():
            print 'calculate bounding box in y/x from lat/lon input'
            pix_box = subset.bbox_geo2radar(geo_box, atr, inps.lookup_file)
        else:
            pix_box = (0,0,0,0)
            pix_box[0:3:2] = subset.coord_geo2radar(geo_box[0:3:2], atr, 'lon')
            pix_box[1:4:2] = subset.coord_geo2radar(geo_box[1:4:2], atr, 'lat')
        print 'Input subset in (lon0,lat0,lon1,lat1): '+str(geo_box)
    print 'Input subset in (x0,  y0,  x1,  y1  ): '+str(pix_box)


    #####Step 2 - subset
    subset_dir = os.path.join(inps.work_dir, 'SUBSET_Y%d_%d_X%d_%d'%(pix_box[1], pix_box[3], pix_box[0], pix_box[2]))
    if not os.path.isdir(subset_dir):
        os.mkdir(subset_dir)
    os.chdir(subset_dir)
    print '\n--------------------------------------------'
    print 'Creating subset datasets ...'
    print "Go to subset directory: "+subset_dir
    lookupFileOld = inps.lookup_file

    if 'Y_FIRST' in atr.keys():
        print '\n--------------------------------------------'
        print 'dataset is in geo coordinates'
        print '    subseting all files in (x0,y0,x1,y1): '+str(pix_box)
        inps = subset.subset_box2inps(inps, pix_box, None)
        inps.ifgram_file    = check_subset_file(inps.ifgram_file,    vars(inps))
        inps.coherence_file = check_subset_file(inps.coherence_file, vars(inps))
        inps.dem_radar_file = check_subset_file(inps.dem_radar_file, vars(inps))
        inps.dem_geo_file   = check_subset_file(inps.dem_geo_file,   vars(inps))
        inps.lookup_file    = check_subset_file(inps.lookup_file,    vars(inps))
        inps.trop_file      = check_subset_file(inps.trop_file,      vars(inps))
    else:
        atr_lut = readfile.read_attribute(inps.lookup_file)
        print '\n--------------------------------------------'
        print 'dataset is in radar coordinates'
        print '    subseting all radar-coord files in (x0,y0,x1,y1): '+str(pix_box)
        inps = subset.subset_box2inps(inps, pix_box, None)
        inps.ifgram_file    = check_subset_file(inps.ifgram_file,    vars(inps))
        inps.coherence_file = check_subset_file(inps.coherence_file, vars(inps))
        inps.dem_radar_file = check_subset_file(inps.dem_radar_file, vars(inps))
        inps.trop_file      = check_subset_file(inps.trop_file,      vars(inps))
        if 'Y_FIRST' not in atr_lut.keys():
            inps.lookup_file  = check_subset_file(inps.lookup_file,  vars(inps))

        print ''
        geo_box = subset.bbox_radar2geo(pix_box, atr, lookupFile=lookupFileOld)
        print '--------------------------------------------'
        print 'dataset is in radar coordinates'
        print '    subseting all geo-coord files in (lon0,lat0,lon1,lat1): '+str(geo_box)
        inps = subset.subset_box2inps(inps, None, geo_box)
        inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps))
        if 'Y_FIRST' in atr_lut.keys():
            inps.lookup_file = check_subset_file(inps.lookup_file, vars(inps))
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

TEMPLATE='''# vim: set filetype=cfg:
##------------------------ pysarApp_template.txt ------------------------##
## 1. Load Data (--load to exit after this step)
## recommended input files for data in radar coordinates:
##     pysar.insarProcessor     = InSAR processor
##     pysar.unwrapFiles        = 'path of all unwrapped interferograms'
##     pysar.corFiles           = 'path of all coherence files'
##     pysar.lookupFile         = 'path of lookup table / mapping transformation file'
##     pysar.demFile.geoCoord   = 'path of DEM in geo   coordinates'
##     pysar.demFile.radarCoord = 'path of DEM in radar coordinates'
## recommended input files for data in geo coordinates:
##     pysar.insarProcessor
##     pysar.unwrapFiles 
##     pysar.corFiles    
##     pysar.dem.geoCoord
## auto - automatic path pattern for Univ of Miami file structure, which are:
##     pysar.insarProcessor     = roipac
##     pysar.unwrapFiles        = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*.unw
##     pysar.corFiles           = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*rlks.cor
##     pysar.lookupFile         = $SCRATCHDIR/$PROJECT_NAME/GEO/*master_date12*/geomap*.trans
##     pysar.demFile.geoCoord   = $SCRATCHDIR/$PROJECT_NAME/DEM/*.dem
##     pysar.demFile.radarCoord = $SCRATCHDIR/$PROJECT_NAME/DONE/*master_date12*/radar*.hgt
pysar.insarProcessor     = auto  #[roipac,        gamma,           isce], auto for roipac
pysar.unwrapFiles        = auto  #[filt*rlks.unw, diff_*rlks.unw,  filt*.unw]
pysar.corFiles           = auto  #[filt*rlks.cor, filt_*rlks.cor,  filt*.cor]
pysar.lookupFile         = auto  #[geomap*.trans, sim*.UTM_TO_RDC, l*.rdr]
pysar.demFile.radarCoord = auto  #[radar*.hgt,    sim*.hgt_sim,    hgt.rdr]
pysar.demFile.geoCoord   = auto  #[*.dem,         sim*.utm.dem,    demLat*.dem.wgs84]


## 1.1 Subset (optional, --subset to exit after this step)
## if both yx and lalo are specified, use lalo option unless a) no lookup file AND b) dataset is in radar coord
pysar.subset.yx       = auto    #[1800:2000,700:800 / no], auto for no
pysar.subset.lalo     = auto    #[31.5:32.5,130.5:131.0 / no], auto for no
pysar.subset.tightBox = auto    #[yes / no], auto for yes, tight bounding box for files in geo coord
pysar.multilook.yx    = auto    #[4,4 / no], auto for no [not implemented yet]

## 1.2 Prepare geometry files
## Prepare incidenceAngle.h5, rangeDistance.h5 files

## 2. Reference in Space
## reference all interferograms to one common point in space
## auto - randomly select a pixel with coherence > minCoherence
pysar.reference.yx            = auto   #[257,151 / auto]
pysar.reference.lalo          = auto   #[31.8,130.8 / auto]

pysar.reference.coherenceFile = auto   #[filename], auto for averageSpatialCoherence.h5
pysar.reference.minCoherence  = auto   #[0.0-1.0], auto for 0.85, minimum coherence for auto method
pysar.reference.maskFile      = auto   #[filename / no], auto for mask.h5


## 3. Unwrapping Error Correction (optional and not recommended)
## unwrapping error correction based on the following two methods:
## a. phase closure (Fattahi, 2015, PhD Thesis)
## b. connecting bridge
pysar.unwrapError.method   = auto   #[bridging / phase_closure / no], auto for no
pysar.unwrapError.maskFile = auto   #[filename / no], auto for no
pysar.unwrapError.ramp     = auto   #[plane / quadratic], auto for plane
pysar.unwrapError.yx       = auto   #[y1_start,x1_start,y1_end,x1_end;y2_start,...], auto for none


## 4. Modify Network (optional)
## Coherence-based network modification = MST + Threshold, by default
## 1) calculate a average coherence value for each interferogram using spatial coherence and input mask (with AOI)
## 2) find a minimum spanning tree (MST) network with inverse of average coherence as weight (keepMinSpanTree)
## 3) for all interferograms except for MST's, exclude those with average coherence < minCoherence.
pysar.network.coherenceBased  = auto  #[yes / no], auto for yes, exclude interferograms with coherence < minCoherence
pysar.network.keepMinSpanTree = auto  #[yes / no], auto for yes, keep interferograms in Min Span Tree network
pysar.network.coherenceFile   = auto  #[filename], auto for coherence.h5
pysar.network.minCoherence    = auto  #[0.0-1.0], auto for 0.7
pysar.network.maskFile        = auto  #[filename, no], auto for [maskLand.h5, mask.h5][0], no for all pixels
pysar.network.maskAoi.yx      = auto  #[y0:y1,x0:x1 / no], auto for no, area of interest for coherence calculation
pysar.network.maskAoi.lalo    = auto  #[lat0:lat1,lon0:lon1 / no], auto for no - use the whole area

pysar.network.tempBaseMax     = auto  #[1-inf, no], auto for no, maximum temporal baseline in days
pysar.network.perpBaseMax     = auto  #[1-inf, no], auto for no, maximum perpendicular spatial baseline in meter
pysar.network.referenceFile   = auto  #[date12_list.txt / Modified_unwrapIfgram.h5 / no], auto for no
pysar.network.excludeDate     = auto  #[20080520,20090817 / no], auto for no
pysar.network.excludeIfgIndex = auto  #[1:5,25 / no], auto for no, list of interferogram number starting from 1
pysar.network.startDate       = auto  #[20090101 / no], auto for no
pysar.network.endDate         = auto  #[20110101 / no], auto for no


## 5. Network Inversion
## Invert network of interferograms into time series
## Temporal coherence (weighted) is calculated using Tazzani et al. (2007, IEEE-TGRS)
## For no/uniform weight approach, use Singular-Value Decomposition (SVD) if network are not fully connected
## For weighted approach, use weighted least square (WLS) solution with the following weighting functions:
##     variance - phase variance due to temporal decorrelation (Yunjun et al., 2017, in prep)
##     no       - no weight, or ordinal inversion with uniform weight (Berardino et al., 2002, IEEE-TGRS)
##     linear   - uniform distribution CDF function (Tong et al., 2016, RSE)
pysar.timeseriesInv.weightFunc    = auto #[variance / no / linear / normal], auto for no, coherence to weight
pysar.timeseriesInv.coherenceFile = auto #[filename / no], auto for coherence.h5, file to read weight data
pysar.timeseriesInv.residualNorm  = auto #[L2 ], auto for L2, norm minimization solution
pysar.timeseriesInv.minTempCoh    = auto #[0.0-1.0], auto for 0.7, min temporal coherence for mask


## 6. Local Oscillator Drift (LOD) Correction (for Envisat only)
## correct LOD if input dataset comes from Envisat
## skip this step for all the other satellites.
## Reference paper: Marinkovic and Larsen, 2013, Proc. LPS


## 7. Tropospheric Delay Correction (optional and recommended)
## correct tropospheric delay using the following methods:
## a. pyaps - use weather re-analysis data (Jolivet et al., 2011, GRL, need to install PyAPS; Dee et al., 2011)
## b. height_correlation - correct stratified tropospheric delay (Doin et al., 2009, J Applied Geop)
## c. base_trop_cor - (not recommend) baseline error and stratified tropo simultaneously (Jo et al., 2010, Geo J)
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.polyOrder    = auto  #[1 / 2 / 3], auto for 1, for height_correlation method
pysar.troposphericDelay.weatherModel = auto  #[ERA / MERRA / NARR], auto for ECMWF, for pyaps method


## 8. Topographic (DEM) Residual Correction (Fattahi and Amelung, 2013, IEEE-TGRS)
## Specify stepFuncDate option if you know there are sudden displacement jump in your area,
## i.e. volcanic eruption, or earthquake, and check timeseriesStepModel.h5 afterward for their estimation.
pysar.topoError              = auto  #[yes / no], auto for yes
pysar.topoError.polyOrder    = auto  #[1-inf], auto for 2, polynomial order of temporal deformation model
pysar.topoError.stepFuncDate = auto  #[20080529,20100611 / no], auto for no, date of step jump
pysar.topoError.excludeDate  = auto  #[20070321,20101120 / txtFile / no], auto for no, date exlcuded for error estimation


## 8.1 Phase Residual Root Mean Square
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
## Recommendation: quadratic for whole image, plane for local/small area
pysar.residualRms.maskFile        = auto  #[filename / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)
pysar.residualRms.saveRefDate     = auto  #[yes / no], auto for yes, save date with min RMS to txt/pdf file.
pysar.residualRms.saveExcludeDate = auto  #[yes / no], auto for yes, save date(s) with RMS > threshold to txt/pdf file.


## 9. Reference in Time
## reference all timeseries to one date in time
## auto - choose date with minimum residual RMS using value from step 8.1
## no   - do not change reference date, keep the defaut one (1st date usually) and skip this step
pysar.reference.date = auto   #[auto / reference_date.txt / 20090214 / no]


## 10. Phase Ramp Removal (optional)
## remove phase ramp for each epoch, useful to check localized deformation, i.e. volcanic, land subsidence, etc.
## [plane, quadratic, plane_range, plane_azimuth, quadratic_range, quadratic_azimuth, baseline_cor, base_trop_cor]
pysar.deramp          = auto  #[no / plane / quadratic], auto for no - no ramp will be removed
pysar.deramp.maskFile = auto  #[filename / no], auto for maskTempCoh.h5, mask file for ramp estimation


## 11. Velocity Inversion
## estimate linear velocity from timeseries, and from tropospheric delay file if exists.
pysar.velocity.excludeDate = auto   #[exclude_date.txt / 20080520,20090817 / no], auto for exclude_date.txt
pysar.velocity.startDate   = auto   #[20070101 / no], auto for no
pysar.velocity.endDate     = auto   #[20101230 / no], auto for no


## 12. Post-processing (geocode, output to Google Earth, UNAVCO, etc.)
## 12.1 Geocode
## For data processed by ROI_PAC/Gamma, output resolution for geocoded file is the same as their lookup table file.
## For data processed by ISCE/Doris, output resolution is assign by user with resolution option:
## 1) float number - resolution in degree, 0.001 by default, around 100 m on equator
## 2) file name    - use the resolution from a file in geo coordinates, e.g. demGeo.h5
pysar.geocode            = auto  #[yes / no], auto for yes
pysar.geocode.resolution = auto  #[0.0-inf / filename], auto for 0.001 (~100 m), output resolution for ISCE processor

## 12.2 Export to other formats
## To update UNAVCO file with new acquisitions, enabling update mode, a.k.a. put XXXXXXXX as endDate 
## in filename if endDate < 1 year
pysar.save.unavco  = auto   #[yes / update / no], auto for no, save timeseries to UNAVCO InSAR Archive format
pysar.save.kml     = auto   #[yes / no], auto for yes, save geocoded velocity to Google Earth KMZ file
pysar.save.geotiff = auto   #[yes / no], auto for no, save geocoded velocity to Geotiff format [not implemented yet]

## 12.3 Plot
pysar.plot = auto   #[yes / no], auto for yes, plot files generated by pysarApp default processing to PIC folder
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
                             '    1) miami_path = True in pysar/__init__.py\n'+\
                             '    2) environmental variable $SCRATCHDIR exists\n'+\
                             '    3) input custom template with basename same as projectName\n')
    parser.add_argument('-g', dest='generate_template', action='store_true',\
                        help='Generate default template (and merge with custom template), then exit.')
    parser.add_argument('--reset', action='store_true',\
                        help='Reset files attributes to re-run pysarApp.py after loading data by:\n'+\
                             '    1) removing ref_y/x/lat/lon for unwrapIfgram.h5 and coherence.h5\n'+\
                             '    2) set drop_ifgram=no for unwrapIfgram.h5 and coherence.h5')
    parser.add_argument('--load-data', dest='load_dataset', action='store_true',\
                        help='Step 1. Load/check dataset, then exit')
    parser.add_argument('--subset-data', dest='subset_dataset', action='store_true',\
                        help='Step 1.1 Subset the whole dataset with setting in template, then exit')
    parser.add_argument('--modify-network', dest='modify_network', action='store_true',\
                        help='Step 4. Modify the network, then exit')
    parser.add_argument('--inverse-network', dest='inverse_network', action='store_true',\
                        help='Step 5. Inverse network of interferograms into time-series, then exit')

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
    # Copy from PROCESS directory
    file_list = ['unavco_attributes.txt', 'bl_list.txt']
    try:
        process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
        file_list = [process_dir+'/'+i for i in file_list]
        file_list = ut.get_file_list(file_list, abspath=True)
        for file in file_list:
            if ut.update_file(os.path.basename(file), file, check_readable=False):
                shutil.copy2(file, inps.work_dir)
                print 'copy '+os.path.basename(file)+' to work directory'
    except: pass

    # Copy from SLC directory
    file_list = ['summary*slc.jpg']
    try:
        slc_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/SLC'
        file_list = [slc_dir+'/'+i for i in file_list]
        for file in file_list:
            if ut.update_file(os.path.basename(file), file, check_readable=False):
                shutil.copy2(file, inps.work_dir)
                print 'copy '+os.path.basename(file)+' to work directory'
    except: pass


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
        sys.exit('Exit as planned after template file generation.')

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

    # Get unavco_attributes.txt file name
    try:    inps.unavco_atr_file = ut.get_file_list('unavco_attributes.txt', abspath=True)[0]
    except: inps.unavco_atr_file = None;   print 'No UNAVCO attributes file found.'


    #########################################
    # Loading Data
    #########################################
    print '\n*************** Load Data ****************'
    loadCmd = 'load_data.py --dir '+inps.work_dir+' --template '+inps.template_file
    if inps.custom_template_file:
        loadCmd += ' '+inps.custom_template_file+' --project '+inps.project_name
    print loadCmd
    status = subprocess.Popen(loadCmd, shell=True).wait()
    os.chdir(inps.work_dir)

    print '--------------------------------------------'
    inps = ut.check_loaded_dataset(inps.work_dir, inps)
    if not inps.ifgram_file:
        print '\nERROR: No interferograms file found!\n'
        sys.exit('Exit.')

    atr = readfile.read_attribute(inps.ifgram_file)
    inps.coord_type = 'radar'
    if 'Y_FIRST' in atr.keys():
        inps.coord_type = 'geo'

    #if inps.unavco_atr_file:
    #    atrCmd = 'add_attribute.py '+inps.ifgram_file+' '+inps.unavco_atr_file       # FA 11/17: add attributes from template file until unavco_attributes is generated automatically
    #import pdb; pdb.set_trace()
    if inps.custom_template_file:
        atrCmd = 'add_attribute.py '+inps.ifgram_file+' '+inps.custom_template_file

        print atrCmd
        status = subprocess.Popen(atrCmd, shell=True).wait()
        if status is not 0:
            print '\nError while adding UNAVCO attributes to unwrapped interferograms file.\n'
            sys.exit(-1)
    if inps.load_dataset:
        sys.exit('Exit as planned after loading/checking the dataset')

    if inps.reset:
        print 'Reset dataset attributtes for a fresh re-run with options from %s' % os.path.basename(inps.template_file)
        print '-----------------------------------------------------------------------------------'
        # Reset network
        networkCmd = 'modify_network.py '+inps.ifgram_file
        if inps.coherence_file:
            networkCmd +=  ' '+inps.coherence_file
        networkCmd += ' --reset'
        print networkCmd
        status = subprocess.Popen(networkCmd, shell=True).wait()
        if status is not 0:
            print '\nError while resetting the network of interferograms.\n'
            sys.exit(-1)

        # Reset reference pixel
        seedCmd = 'seed_data.py '+inps.ifgram_file+' --reset'
        print seedCmd
        status = subprocess.Popen(seedCmd, shell=True).wait()
        if status is not 0:
            print '\nError while resetting the reference pixel in space.\n'
            sys.exit(-1)


    #########################################
    # Check the subset (Optional)
    #########################################
    if inps.lookup_file and template['pysar.subset.tightBox'] in ['yes','auto']:
        ###Tight subset DEM in geo coord
        #subCmd = 'subset.py '+inps.dem_geo_file+' --tight'
        #print subCmd
        #outName = os.path.splitext(inps.dem_geo_file)[0]+'_tight'+os.path.splitext(inps.dem_geo_file)[1]
        #if ut.update_file(outName, inps.dem_geo_file):
        #    status = subprocess.Popen(subCmd, shell=True).wait()        
        #if status is 0 and os.path.isfile(outName):
        #    inps.dem_geo_file = outName

        ##Tight subset lookup table in geo coord (roipac/gamma)
        atr_lut = readfile.read_attribute(inps.lookup_file)
        if 'Y_FIRST' in atr_lut.keys():
            subCmd = 'subset.py '+inps.lookup_file+' --tight'
            print subCmd
            outName = os.path.splitext(inps.lookup_file)[0]+'_tight'+os.path.splitext(inps.lookup_file)[1]
            if ut.update_file(outName, inps.lookup_file):
                status = subprocess.Popen(subCmd, shell=True).wait()
            if status is 0 and os.path.isfile(outName):
                inps.lookup_file = outName

    # Subset based on input template
    if not all(template[key] in ['auto', 'no'] for key in ['pysar.subset.yx','pysar.subset.lalo']):
        print '\n*************** Subset ****************'
        inps = subset_dataset(inps, inps.template_file)

    if inps.subset_dataset:
        sys.exit('Exit as planned after subsetting the dataset')


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

    # Incidence Angle
    inps.inc_angle_file = ut.get_inc_angle_file(abspath=True, print_msg=False)
    if not inps.inc_angle_file:
        print 'creating incidence angle from file: '+inps.ifgram_file
        inps.inc_angle_file = 'incidenceAngle.h5'
        incAngleCmd = 'incidence_angle.py %s %s' % (inps.ifgram_file, inps.inc_angle_file)
        print incAngleCmd
        status = subprocess.Popen(incAngleCmd, shell=True).wait()
    print 'incidence angle file: %s' % (inps.inc_angle_file)

    # Slant range distance
    inps.range_dist_file = ut.get_range_distance_file(abspath=True, print_msg=False)
    if not inps.range_dist_file:
        print 'creating slant range distance from file: '+inps.ifgram_file
        inps.range_dist_file = 'rangeDistance.h5'
        rangeDistCmd = 'range_distance.py %s %s' % (inps.ifgram_file, inps.range_dist_file)
        print rangeDistCmd
        status = subprocess.Popen(rangeDistCmd, shell=True).wait()
    print 'slant range distance file: %s' % (inps.range_dist_file)


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print '\n**********  Reference in space  ***************'
    seedCmd = 'seed_data.py '+inps.ifgram_file+' --template '+inps.template_file+' --mark-attribute'
    resetCmd = 'seed_data.py '+inps.ifgram_file+' --reset'

    ## Skip calling seed command only if 1) ref_y/x exists AND 2) pysar.reference.yx/lalo == auto
    run_seedCmd = True
    atr = readfile.read_attribute(inps.ifgram_file)
    if 'ref_x' in atr.keys():
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        length = int(atr['FILE_LENGTH'])
        width = int(atr['WIDTH'])
        print 'Find reference pixel info from %s in y/x: [%d, %d]' % (os.path.basename(inps.ifgram_file), ref_y, ref_x)
        if not (0 <= ref_y <= length and 0 <= ref_x <= width):
            print 'existed reference pixel is out of data coverage.'
            print '1) reset reference pixel info'
            print resetCmd
            status = subprocess.Popen(resetCmd, shell=True).wait()
            print '2) re-run seed_data.py to select new refernce pixel.'
        else:
            print '----------------------------------------------------------------------------------------'
            print 'To remove reference pixel info, use seed_data.py --reset option:'
            print resetCmd
            print '----------------------------------------------------------------------------------------'
            prefix = 'pysar.reference.'
            if all(template[prefix+i] in ['auto','no'] for i in ['yx','lalo']):
                run_seedCmd = False
                print 'No specific coordinates input found, no need to re-select reference pixel'
            else:
                print 'Specific coordinates input found, re-select reference pixel with options from template file'

    if run_seedCmd:
        print seedCmd
        status = subprocess.Popen(seedCmd, shell=True).wait()
        if status is not 0:
            print '\nError while finding reference pixel in space.\n'
            sys.exit(-1)


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    if template['pysar.unwrapError.method'] not in ['auto','no']:
        print '\n**********  Unwrapping Error Correction  **************'
        outName = os.path.splitext(inps.ifgram_file)[0]+'_unwCor.h5'
        unwCmd='unwrap_error.py '+inps.ifgram_file+' --mask '+inps.mask_file+' --template '+inps.template_file
        print unwCmd
        #if ut.update_file(outName):
        if ut.update_file(outName, inps.ifgram_file):
            print 'This might take a while depending on the size of your data set!'
            status = subprocess.Popen(unwCmd, shell=True).wait()
            if status is not 0:
                print '\nError while correcting phase unwrapping errors.\n'
                sys.exit(-1)
        inps.ifgram_file = outName


    #########################################
    # Network Modification (Optional)
    #########################################
    print '\n*************** Modify Network ****************'
    networkCmd = 'modify_network.py --template '+inps.template_file+' '+inps.ifgram_file
    if inps.coherence_file:
        networkCmd += ' '+inps.coherence_file
    print networkCmd
    print '----------------------------------------------------------------------------------------'
    print 'To use all interferograms in the file, run modify_network.py with --reset option to restore all pairs info.'
    msg_str = 'modify_network.py '+inps.ifgram_file
    if inps.coherence_file:
        msg_str +=  ' '+inps.coherence_file
    msg_str += ' --reset'
    print msg_str
    print '----------------------------------------------------------------------------------------'
    status = subprocess.Popen(networkCmd, shell=True).wait()
    if status is not 0:
        print '\nError while modifying the network of interferograms.\n'
        sys.exit(-1)

    # Plot network colored in spatial coherence
    print '--------------------------------------------'
    plotCmd = 'plot_network.py '+inps.ifgram_file+' --coherence '+inps.coherence_file+\
              ' --template '+inps.template_file+' --nodisplay'
    print plotCmd
    inps.coh_spatialAverage_file = os.path.splitext(inps.coherence_file)[0]+'_spatialAverage.txt'
    if ut.update_file('Network.pdf', check_readable=False, \
                      inFile=[inps.ifgram_file, inps.coh_spatialAverage_file, inps.template_file]):
        status = subprocess.Popen(plotCmd, shell=True).wait()

    if inps.modify_network:
        sys.exit('Exit as planned after network modification.')


    #########################################
    # Inversion of Interferograms
    ########################################
    print '\n**********  Network Inversion to Time Series  ********************'
    inps.timeseries_file = 'timeseries.h5'
    inps.temp_coh_file = 'temporalCoherence.h5'
    invertCmd = 'ifgram_inversion.py '+inps.ifgram_file+' --template '+inps.template_file
    print invertCmd
    if ut.update_file(inps.timeseries_file, inps.ifgram_file):
        status = subprocess.Popen(invertCmd, shell=True).wait()
        if status is not 0:
            print '\nError while inverting network of interferograms to time-series.\n'
            sys.exit(-1)

    print '\n--------------------------------------------'
    print 'Update Mask based on Temporal Coherence ...'
    # Read template option
    inps.min_temp_coh = 0.7
    key = 'pysar.timeseriesInv.minTempCoh'
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
        status = subprocess.Popen(maskCmd, shell=True).wait()
        if status is not 0:
            print '\nError while generating mask file from temporal coherence.\n'
            sys.exit(-1)
    inps.mask_file = outName

    if inps.inverse_network:
        sys.exit('Exit as planned after network inversion.')


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
        #if 'pysar.troposphericDelay.method' in template.keys():
        #    print '    Continue without tropospheric correction ...'
        #    template['pysar.troposphericDelay.method'] = 'no'
        #if template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']:
        #    template['pysar.deramp'] = 'no'
        print '++++++++++++++++++++++++++++++++++++++++++++++'


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   for Envisat data in radar coord only
    ############################################## 
    sar_mission = atr['PLATFORM'].lower()
    if sar_mission.startswith('env'):
        print '\n**********  Local Oscillator Drift correction for Envisat  ********'
        outName = os.path.splitext(inps.timeseries_file)[0]+'_LODcor.h5'
        lodCmd = 'lod.py %s -r %s' % (inps.timeseries_file, inps.range_dist_file)
        print lodCmd
        if ut.update_file(outName, [inps.timeseries_file, inps.range_dist_file]):
            status = subprocess.Popen(lodCmd, shell=True).wait()
            if status is not 0:
                print '\nError while correcting Local Oscillator Drift.\n'
                sys.exit(-1)
        inps.timeseries_file = outName


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
            status = subprocess.Popen(tropCmd, shell=True).wait()
            if status is not 0:
                print '\nError while correcting tropospheric delay.\n'
                sys.exit(-1)
        inps.timeseries_file = outName

    elif inps.trop_method == 'pyaps':
        print 'Atmospheric correction using Weather Re-analysis dataset (using PyAPS software)'
        print 'Weather Re-analysis dataset: '+inps.trop_model
        tropCmd = 'tropcor_pyaps.py '+inps.timeseries_file+' -d '+demFile+' -s '+inps.trop_model+\
                  ' --weather-dir '+inps.work_dir+'/../WEATHER'+' -i '+inps.inc_angle_file
        print tropCmd
        outName = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.trop_model+'.h5'
        if ut.update_file(outName, inps.timeseries_file):
            try:
                inps.trop_file = ut.get_file_list(inps.trop_model+'.h5')[0]
                tropCmd = 'diff.py '+inps.timeseries_file+' '+inps.trop_file+' -o '+outName
                print 'Use existed tropospheric delay file: '+inps.trop_file
                print tropCmd
            except: pass
            status = subprocess.Popen(tropCmd, shell=True).wait()
            if status is not 0:
                print '\nError while correcting tropospheric delay, try the following:'
                print '1) Check the installation of PyAPS (http://earthdef.caltech.edu/projects/pyaps/wiki/Main)'
                print '2) Use other tropospheric correction method, height-correlation, for example'
                print '3) or turn off the option by setting pysar.troposphericDelay.method = no in template file.\n'
                sys.exit(-1)
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
    topoCmd = 'dem_error.py %s --template %s -i %s -r %s -o %s' %\
              (inps.timeseries_file, inps.template_file, inps.inc_angle_file, inps.range_dist_file, outName)
    print topoCmd
    inps.timeseries_resid_file = None
    if template['pysar.topoError'] in ['yes','auto']:
        print 'Correcting topographic residuals using method from Fattahi and Amelung, 2013, TGRS ...'
        if ut.update_file(outName, inps.timeseries_file):
            status = subprocess.Popen(topoCmd, shell=True).wait()
            if status is not 0:
                print '\nError while correcting topographic phase residual.\n'
                sys.exit(-1)
        inps.timeseries_file = outName
        inps.timeseries_resid_file = 'timeseriesResidual.h5'
    else:
        print 'No correction for topographic residuals.'


    ##############################################
    # Timeseries Residual Standard Deviation
    ##############################################
    print '\n**********  Timeseries Residual Root Mean Square  *******'
    if inps.timeseries_resid_file:
        rmsCmd = 'timeseries_rms.py '+inps.timeseries_resid_file+' --template '+inps.template_file
        print rmsCmd
        status = subprocess.Popen(rmsCmd, shell=True).wait()
        if status is not 0:
            print '\nError while calculating RMS of time series phase residual.\n'
            sys.exit(-1)
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
        status = subprocess.Popen(refCmd, shell=True).wait()
        if status is not 0:
            print '\nError while changing reference date.\n'
            sys.exit(-1)

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
        derampCmd = None

        # Get executable command and output name
        if inps.deramp_method in ['plane', 'quadratic', 'plane_range', 'quadratic_range',\
                                  'plane_azimuth', 'quadratic_azimuth']:
            derampCmd = 'remove_plane.py '+inps.timeseries_file+' -s '+inps.deramp_method+' -m '+inps.mask_file
            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.deramp_method+'.h5'

        elif inps.deramp_method in ['baseline_cor','baselinecor']:
            if not 'X_FIRST' in atr.keys():
                derampCmd = 'baseline_error.py '+inps.timeseries_file+' '+inps.mask_file
                outName = os.path.splitext(inps.timeseries_file)[0]+'_baselineCor.h5'
            else:
                warnings.warn('BaselineCor method can only be applied in radar coordinate, skipping correction')

        elif inps.deramp_method in ['base_trop_cor','basetropcor','baselinetropcor']:
            if not 'X_FIRST' in atr.keys():
                print 'Joint estimation of Baseline error and tropospheric delay [height-correlation approach]'
                try:    poly_order = template['pysar.troposphericDelay.polyOrder']
                except: poly_order = '1'
                derampCmd = 'baseline_trop.py '+inps.timeseries_file+' '+inps.dem_radar_file+' '+\
                            poly_order+' range_and_azimuth'
                outName = os.path.splitext(inps.timeseries_file)[0]+'_baseTropCor.h5'
            else:
                warnings.warn('BaselineTropCor method can only be applied in radar coordinate, skipping correction')
        else:
            warnings.warn('Unrecognized phase ramp method: '+template['pysar.deramp'])

        # Execute command
        if derampCmd:
            print derampCmd
            if ut.update_file(outName, inps.timeseries_file):
                status = subprocess.Popen(derampCmd, shell=True).wait()
                if status is not 0:
                    print '\nError while removing phase ramp for each acquisition of time-series.\n'
                    sys.exit(-1)
            inps.timeseries_file = outName
    else:
        print 'No phase ramp removal.'


    #############################################
    # Velocity and rmse maps
    #############################################
    print '\n**********  Velocity estimation  **********************'
    inps.vel_file = 'velocity.h5'
    velCmd = 'timeseries2velocity.py '+inps.timeseries_file+' --template '+inps.template_file+' -o '+inps.vel_file
    print velCmd
    if ut.update_file(inps.vel_file, [inps.timeseries_file, inps.template_file]):
        status = subprocess.Popen(velCmd, shell=True).wait()
        if status is not 0:
            print '\nError while estimating linear velocity from time-series.\n'
            sys.exit(-1)

    # Velocity from Tropospheric delay
    if inps.trop_file:
        suffix = os.path.splitext(os.path.basename(inps.trop_file))[0]
        suffix = suffix[0].upper()+suffix[1:].lower()
        inps.trop_vel_file = 'velocity'+suffix+'.h5'
        velCmd = 'timeseries2velocity.py '+inps.trop_file+' --template '+inps.template_file+' -o '+inps.trop_vel_file
        print velCmd
        if ut.update_file(inps.trop_vel_file, [inps.trop_file, inps.template_file]):
            status = subprocess.Popen(velCmd, shell=True).wait()


    ############################################
    # Post-processing
    # Geocodeing, masking and save to KML 
    ############################################
    print '\n**********  Post-processing  ********************************'

    # Geocoding
    if 'Y_FIRST' in atr.keys():
        inps.geo_vel_file = inps.vel_file
        inps.geo_temp_coh_file = inps.temp_coh_file
        inps.geo_timeseries_file = inps.timeseries_file
    else:
        inps.geo_vel_file = None
        inps.geo_temp_coh_file = None
        inps.geo_timeseries_file = None

    key = 'pysar.geocode'
    if template[key] in ['auto','yes'] and 'Y_FIRST' not in atr.keys():
        print '\n--------------------------------------------'
        inps.geo_vel_file        = check_geocode_file(inps.lookup_file, inps.vel_file,        inps.template_file)
        inps.geo_temp_coh_file   = check_geocode_file(inps.lookup_file, inps.temp_coh_file,   inps.template_file)
        inps.goe_timeseries_file = check_geocode_file(inps.lookup_file, inps.timeseries_file, inps.template_file)


    # Mask in geo coord
    inps.geo_mask_file = None
    if inps.geo_temp_coh_file:
        # Generate mask in geo coord
        print '\n--------------------------------------------'
        outName = 'maskTempCoh.h5'
        if os.path.basename(inps.geo_temp_coh_file).startswith('geo_'):
            outName = 'geo_'+outName
        maskCmd = 'generate_mask.py '+inps.geo_temp_coh_file+' -m '+str(inps.min_temp_coh)+' -o '+outName
        print maskCmd
        if ut.update_file(outName, inps.geo_temp_coh_file):
            status = subprocess.Popen(maskCmd, shell=True).wait()
        inps.geo_mask_file = outName

        # Mask geo_velocity file
        if inps.geo_vel_file and inps.geo_mask_file:
            outName = os.path.splitext(inps.geo_vel_file)[0]+'_masked.h5'
            maskCmd = 'mask.py '+inps.geo_vel_file+' -m '+inps.geo_mask_file+' -o '+outName
            print maskCmd
            if ut.update_file(outName, [inps.geo_vel_file, inps.geo_mask_file]):
                status = subprocess.Popen(maskCmd, shell=True).wait()
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
            status = subprocess.Popen(kmlCmd, shell=True).wait()


    #############################################
    # Save to UNAVCO InSAR Archive format
    #############################################
    if template['pysar.save.unavco'] in ['yes','update']:
        print '\n*********  Output to UNAVCO InSAR Archive Format  ***********'
        if 'Y_FIRST' not in atr.keys() and not inps.lookup_file:
            warnings.warn('Dataset is in radar coordinates without lookup table file.'+\
                          'Can not geocode.'+\
                          'Skip saving.')
        else:
            # 1. Time series file
            inps.geo_timeseries_file = check_geocode_file(inps.lookup_file, inps.timeseries_file, inps.template_file)
            # Add UNAVCO attributes
            if inps.unavco_atr_file:
                atrCmd = 'add_attribute.py '+inps.geo_timeseries_file+' '+inps.unavco_atr_file
                print atrCmd
                status = subprocess.Popen(atrCmd, shell=True).wait()
                if status is not 0:
                    print '\nError while adding UNAVCO attributes to time series file.\n'
                    sys.exit(-1)

            # 2. Temporal Coherence
            inps.geo_temp_coh_file = check_geocode_file(inps.lookup_file, inps.temp_coh_file, inps.template_file)

            # 3. Mask file
            if not inps.geo_mask_file:
                outName = 'maskTempCoh.h5'
                if os.path.basename(inps.geo_temp_coh_file).startswith('geo_'):
                    outName = 'geo_'+outName
                maskCmd = 'generate_mask.py '+inps.geo_temp_coh_file+' -m '+str(inps.min_temp_coh)+' -o '+outName
                print maskCmd
                if ut.update_file(outName, inps.geo_temp_coh_file):
                    status = subprocess.Popen(maskCmd, shell=True).wait()
                    if status is not 0:
                        print '\nError while generating mask file.\n'
                        sys.exit(-1)
                inps.geo_mask_file = outName

            # 4. Incidence Angle
            inps.inc_angle_file = 'incidenceAngle.h5'
            if ut.update_file(inps.inc_angle_file, inps.timeseries_file):
                incAngleCmd = 'incidence_angle.py '+inps.timeseries_file+' '+inps.inc_angle_file
                print incAngleCmd
                status = subprocess.Popen(incAngleCmd, shell=True).wait()
                if status is not 0:
                    print '\nError while generating incidence angle file.\n'
                    sys.exit(-1)
            inps.geo_inc_angle_file = check_geocode_file(inps.lookup_file, inps.inc_angle_file, inps.template_file)

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
            if template['pysar.save.unavco'] == 'update':
                unavcoCmd += ' --update '
            print unavcoCmd
            if ut.update_file(inps.unavco_file, [inps.geo_timeseries_file, inps.geo_temp_coh_file, inps.geo_mask_file,\
                                                 inps.geo_inc_angle_file, inps.dem_geo_file], check_readable=False):
                status = subprocess.Popen(unavcoCmd, shell=True).wait()
                if status is not 0:
                    print '\nError while generating UNAVCO InSAR arhive time series file.\n'
                    sys.exit(-1)


    #############################################
    # Plot Figures
    #############################################
    if template['pysar.plot'] in ['yes','auto']:
        print '\n*********  Plot and Save pysarApp runing results to PIC  ***********'
        inps.plot_sh_file = 'plot_pysarApp.sh'

        # Copy to workding directory if not existed yet.
        if not os.path.isfile('./'+inps.plot_sh_file):
            print 'copy $PYSAR_HOME/shellscripts/'+inps.plot_sh_file+' to working directory'
            try:
                shutil.copy2(ut.which(inps.plot_sh_file), './')
            except:
                print 'WARNING: no '+inps.plot_sh_file+' found in the environment variable path, skip plotting.'
        print 'for better performance, edit the input parameters in '+inps.plot_sh_file+' and re-run this script.'

        #if ut.update_file('PIC', [inps.plot_sh_file, inps.template_file], check_readable=False):
        plotCmd = './'+inps.plot_sh_file
        print plotCmd
        status = subprocess.Popen(plotCmd, shell=True).wait()


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

