#! /usr/bin/env python
###############################################################################
# 
# Project: PySAR 
# Purpose: Python Module for InSAR Time-series Analysis
# Author: Heresh Fattahi
# Created: July 2013
#
###############################################################################
# Copyright (c) 2013, Heresh Fattahi
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
#
# Yunjun, Feb 2015: Update mask, generate incident angle file
# Yunjun, Oct 2015: Add find_filename(), check_subset()
#                   Move radar_or_geo() to _pysar_utilities.py
#                   Update name for pysar.dem.* and pysar.trop*
#                   Finished pysar.subset.yx option.
#                   Add check_mask(), check_geocode()
# Yunjun, Nov 2015: Add pysar.kml option, workDir input option
# Yunjun, Dec 2016: Add command line parser


import os
import sys
import glob
import time
import argparse

import h5py
import numpy as np

import pysar
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar.subset as subset


def check_isfile(File):
    '''Check if input file exists and readable.'''
    if not File or not os.path.isfile(File):
        return False
     
    try:
        atr = readfile.read_attribute(File)
        return True
    except:
        print File+' exists, but can not read, remove it.'
        rmCmd = 'rm '+File
        print rmCmd
        os.system(rmCmd)
        return False


def check_subset_file(File, inps_dict, outFile=None, overwrite=False):
    '''Subset input file or use existed subseted file.'''
    if not File:     return None
    if not outFile:
        if os.getcwd() == inps_dict['work_dir']:
            outFile = 'subset_'+os.path.basename(File)
        else:
            # if current dir is not original timeseries directory (e.g. TIMESERIES/subset)
            # use the same filename (but in different directories)
            outFile = os.path.basename(File)
    
    if check_isfile(outFile) and not overwrite:
        print outFile+' already exists, no need to re-subset.'
    else:
        outFile = subset.subset_file(File, inps_dict, outFile)
    return outFile


def check_geocode_file(geomapFile, File, outFile=None):
    '''Geocode input file or use existed geocoded file.'''
    if not geomapFile:
        print 'WARNING: No geomap*.trans file found! Skip geocoding.'
        return None
    if not File:  return None
    
    if not outFile:  outFile = 'geo_'+os.path.basename(File)
    if check_isfile(outFile):
        print outFile+' already exists, no need to re-geocode.'
    else:
        geocodeCmd = 'geocode.py '+os.path.basename(geomapFile)+' '+File
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
    
    print '--------------------------------------------'
    print 'subseting dataset in geo coord geo_box4geo: '+str(geo_box4geo)
    inps = subset.subset_box2inps(inps, None, geo_box4geo)
    inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps))
    inps.geomap_file  = check_subset_file(inps.geomap_file, vars(inps))
 
    print '--------------------------------------------'
    print 'subseting dataset in radar coord pix_box4rdr: '+str(pix_box4rdr)
    inps = subset.subset_box2inps(inps, pix_box4rdr, None)
    inps.ifgram_file    = check_subset_file(inps.ifgram_file, vars(inps))
    inps.mask_file      = check_subset_file(inps.mask_file, vars(inps))
    inps.dem_radar_file = check_subset_file(inps.dem_radar_file, vars(inps))
    inps.coherence_file = check_subset_file(inps.coherence_file, vars(inps))
    inps.spatial_coherence_file = check_subset_file(inps.spatial_coherence_file, vars(inps))

    return inps


def create_subset_dataset(inps, pix_box=None, geo_box=None):
    '''Create/prepare subset of datasets in different folder for time series analysis.
    For dataset (unwrapped interferograms) in radar coord, only support subset in row/col or y/x
    For dataset (unwrapped interferograms) in geo coord, lalo has higher priority than yx, if both are specified.
    '''
    # Subset directory
    subset_dir = inps.work_dir+'/subset'
    if not os.path.isdir(subset_dir):  os.mkdir(subset_dir)
    os.chdir(subset_dir)
    print '\n--------------------------------------------'
    print 'Creating subset datasets ...'
    print "Go to directory: "+subset_dir
    
    atr = readfile.read_attribute(inps.ifgram_file)
    if not 'X_FIRST' in atr.keys():
        print 'Loaded dataset is in radar coordinate.'
        geomap_file_orig = inps.geomap_file
        
        # subset.lalo has higher priority than subset.yx, except when no geomap*.trans exists.
        # don't subset if only subset.lalo without geomap*.trans exsits, because glob2radar won't be accurate.
        if geo_box:
            if inps.geomap_file:
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
            pix_box = subset.bbox_geo2radar(geo_box, inps.ifgram_file, geomap_file_orig)
        else:
            print 'use subset input in y/x'
            print 'calculate corresponding bounding box in geo coordinate.'
            geo_box = subset.bbox_radar2geo(pix_box, inps.ifgram_file, geomap_file_orig)
        
        # subset
        inps = subset_dataset(inps, geo_box, pix_box)

    else:
        print 'Loaded dataset is in geo coordinate.'
        
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
               PySAR v1.2, Jan 2017
 Geodesy Lab, University of Miami, Maimi FL, USA
_________________________________________________
'''
#generate_from: http://patorjk.com/software/taag/

TEMPLATE='''template:
# Input Data (not needed for Miami user)
pysar.unwrapFiles    = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*.unw
pysar.corFiles       = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.cor
pysar.wrapFiles      = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.int     #optional
pysar.geomap         = /SanAndreasT356EnvD/PROCESS/GEO/*050102-070809*/geomap*.trans
pysar.dem.radarCoord = /SanAndreasT356EnvD/PROCESS/DONE/*050102-070809*/radar*.hgt
pysar.dem.geoCoord   = /SanAndreasT356EnvD/DEM/srtm1_30m.dem                   #optional

pysar.network.reference       = date12.list         #optional
pysar.network.coherenceBase   = yes                 #optional, auto for yes

pysar.subset.yx          = 1800:2000,700:800        #optional, auto/no/off for whole area
pysar.subset.lalo        = 31.5:32.5,130.5:131.0    #optional, auto/no/off for whole area

pysar.reference.yx       = 257 , 151                #optional, auto for max coherence selection
pysar.reference.lalo     = 31.8, 130.8              #optional, auto for max coherence selection
pysar.reference.date     = 20090120                 #optional, auto for the first date
 
pysar.troposphericDelay.method        = pyaps   #[height_correlation], auto for no tropospheric correction
pysar.troposphericDelay.polyOrder     = 1       #for height_correlation method
pysar.troposphericDelay.weatherModel  = ECMWF   #[ERA, MERRA, NARR], for pyaps method

pysar.topoError = yes               #[no], auto for yes
pysar.deramp    = plane             #[plane, quadratic, baseline_cor, base_trop_cor], auto for no
pysar.geocode   = yes               #[no], auto for yes
'''

EXAMPLE='''example:
  pysarApp.py  SanAndreasT356EnvD.template
  pysarApp.py  SanAndreasT356EnvD.template  --dir ~/insarlab/SanAndreasT356EnvD/TIMESERIES
'''

UM_FILE_STRUCT='''
    scratch/                 # $SCRATCHDIR defined in environmental variable
        SanAndreasT356EnvD/  # my_projectName, same as the basename of template file
            DEM/             # DEM file(s) (for topographic phase and geocode)
            DOWNLOAD/        # (optional) Data downloaded from agencies
            PROCESS/         # Interferograms processed by ROI_PAC, Gamma, ISCE, ... 
            RAW/             # (optional) Raw SAR data untared from DOWNLOAD directory
            SLC/             # (optional) SLC SAR data after focusing from RAW directory
            TIMESERIES/           # PySAR work directory for time series analysis
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description=LOGO,
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)
    
    parser.add_argument('template_file', help='template for PySAR setting\n'+\
                        'If item commented out or left empty, use auto/default value.')
    parser.add_argument('--dir', dest='work_dir',\
                        help='specify working directory; default is:\n'+\
                             '$SCRATCHDIR/my_projectName/TIMESERIES\n'+\
                             '    If using University of Miami file/dir structure\n'+\
                             '    To turn it off, change miami_path value to False in pysar/__init__.py\n'
                             './ - current directory\n'+\
                             '    To set this as your default permanetly,'+\
                             ' change miami_path value to False in pysar/__init__.py')
    parser.add_argument('-v','--version', action='version', version='%(prog)s 2.0')

    inps = parser.parse_args()
    return inps


##########################################################################
def main(argv):
    start = time.time()
    inps = cmdLineParse()


    #########################################
    # Initiation
    #########################################
    print LOGO
    # Read template
    inps.project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print 'Project name: '+inps.project_name
    inps.template_file = os.path.abspath(inps.template_file)
    template = readfile.read_template(inps.template_file)
    for key in template.keys():
        if template[key].lower() == 'default':        template[key] = 'auto'
        if template[key].lower() in ['off','false']:  template[key] = 'no'
        if template[key].lower() in ['on','true']:    template[key] = 'yes'
    if 'pysar.deramp' in template.keys():
        template['pysar.deramp'] = template['pysar.deramp'].lower().replace('-','_')
    if 'pysar.troposphericDelay.method' in template.keys():
        template['pysar.troposphericDelay.method'] = template['pysar.troposphericDelay.method'].lower().replace('-','_')

    # work directory
    if not inps.work_dir:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ:
            inps.work_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+"/PYSAR"
            print 'Use file/dir structure in University of Miami.'+\
                  '(To turn it off, change miami_path value to False in pysar/__init__.py)'
        else:
            inps.work_dir = os.getcwd()
    else:
        inps.work_dir = os.path.abspath(inps.work_dir)
    
    if not os.path.isdir(inps.work_dir):   os.mkdir(inps.work_dir)
    os.chdir(inps.work_dir)
    print "Go to work directory: "+inps.work_dir


    #########################################
    # Loading Data
    #########################################
    print '\n*************** Load Data ****************'
    loadCmd = 'load_data.py '+inps.template_file
    print loadCmd
    os.system(loadCmd)

    print '--------------------------------------------'
    ## Find initial files name/path - required files
    # 1. Unwrapped interferograms
    inps.ifgram_file = 'unwrapIfgram.h5'
    try:    inps.ifgram_file = glob.glob(inps.work_dir+'/'+inps.ifgram_file)[0]
    except: inps.ifgram_file = None
    if inps.ifgram_file:  print 'Unwrapped interferograms: '+inps.ifgram_file
    else:  sys.exit('\nERROR: No interferograms file found!\n')
    
    # 2. Mask
    inps.mask_file = 'Mask.h5'
    try:     inps.mask_file = glob.glob(inps.work_dir+'/'+inps.mask_file)[0]
    except:  inps.mask_file = None
    if not inps.mask_file:
        print 'No mask file found. Creating one using non-zero pixels in file: '+inps.ifgram_file
        inps.mask_file = ut.nonzero_mask(inps.ifgram_file, inps.mask_file)
    print 'Mask: '+inps.mask_file

    ## Find initial files name/path - recommended files (None if not found)
    # 3. Spatial coherence for each interferograms
    inps.coherence_file = 'coherence.h5'
    try:    inps.coherence_file = glob.glob(inps.work_dir+'/'+inps.coherence_file)[0]
    except: inps.coherence_file = None
    if inps.coherence_file:  print 'Coherences: '+inps.coherence_file
    else:  print '\nWARNING: No coherences file found. Cannot use coherence-based network modification without it.\n'
    
    # 4. Average spatial coherence
    inps.spatial_coherence_file = 'average_spatial_coherence.h5'
    try:    inps.spatial_coherence_file = glob.glob(inps.work_dir+'/'+inps.spatial_coherence_file)[0]
    except: inps.spatial_coherence_file = None
    if not inps.coherence_file:
        inps.spatial_coherence_file = ut.temporal_average(inps.coherence_file, inps.spatial_coherence_file)

    # 5. DEM in geo coord
    try:    inps.dem_geo_file = os.path.basename(template['pysar.dem.geoCoord'])
    except: inps.dem_geo_file = '*.dem'
    try:    inps.dem_geo_file = glob.glob(inps.work_dir+'/'+inps.dem_geo_file)[0]
    except: inps.dem_geo_file = None
    if inps.dem_geo_file:  print 'DEM in geo   coord: '+str(inps.dem_geo_file)
    else:  print '\nWARNING: No geo coord DEM found.\n'
    
    # 6. DEM in radar coord
    try:    inps.dem_radar_file = os.path.basename(template['pysar.dem.radarCoord'])
    except: inps.dem_radar_file = 'radar*.hgt'
    try:    inps.dem_radar_file = glob.glob(inps.work_dir+'/'+inps.dem_radar_file)[-1]
    except: inps.dem_radar_file = None
    if inps.dem_radar_file:  print 'DEM in radar coord: '+str(inps.dem_radar_file)
    else:  print '\nWARNING: No radar coord DEM found! Cannot use tropospheric delay correction without it.\n'
    
    # 7. Transform file for geocoding
    try:    inps.geomap_file = os.path.basename(template['pysar.geomap'])
    except: inps.geomap_file = 'geomap*.trans'
    try:    inps.geomap_file = glob.glob(inps.work_dir+'/'+inps.geomap_file)[0]
    except: inps.geomap_file = None
    if inps.geomap_file:     print 'Transform     file: '+str(inps.geomap_file)
    else:  print '\nWARNING: No transform file found! Cannot geocoding without it.\n'


    #########################################
    # Check the subset (Optional)
    #########################################
    print '\n*************** Subset ****************'
    print "Get tight subset of geomap*.trans file and/or DEM file in geo coord"
    print '--------------------------------------------'
    if inps.geomap_file:
        outName = os.path.splitext(inps.geomap_file)[0]+'_tight'+os.path.splitext(inps.geomap_file)[1]
        if check_isfile(outName):
            print '\n'+outName+' already existed.\n'
        else:
            subsetCmd = 'subset.py '+inps.geomap_file+' --footprint '+' -o '+outName
            print subsetCmd
            os.system(subsetCmd)
        inps.geomap_file = outName
        
        # Subset DEM in geo coord
        outName = os.path.splitext(inps.dem_geo_file)[0]+'_tight'+os.path.splitext(inps.dem_geo_file)[1]
        geomap_atr = readfile.read_attribute(inps.geomap_file)
        pix_box, geo_box = subset.get_coverage_box(geomap_atr)
        inps = subset.subset_box2inps(inps, pix_box, geo_box)
        inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps), outName, overwrite=True)

    # Subset based on input template
    if [key for key in template.keys()\
        if ('pysar.subset' in key and template[key].lower() not in ['auto','no'])]:
        # Read subset option from template file, and return None if lalo/yx is not specified.
        pix_box, geo_box = subset.read_subset_template2box(inps.template_file)
        inps = create_subset_dataset(inps, pix_box, geo_box)
    else:
        print '\nNo Subset selected. Processing the whole area.'


    #########################################
    # Network Modification (Optional)
    #########################################
    if [key for key in template.keys() if 'pysar.network' in key]:
        print '\n*************** Modify Network ****************'
        outName = 'Modified_'+os.path.basename(inps.ifgram_file)
        if check_isfile(outName):
            print '\n'+outName+' already existed, no need to re-modify network.\n'
        else:
            networkCmd = 'modify_network.py '+inps.ifgram_file+' '+inps.coherence_file+\
                         ' --template '+inps.template_file+' --mask '+inps.mask_file+' --plot'
            print networkCmd
            os.system(networkCmd)
        
        if check_isfile(outName):  inps.ifgram_file = outName
        outName = 'Modified_'+os.path.basename(inps.mask_file)
        if check_isfile(outName):  inps.mask_file = outName
        outName = 'Modified_'+os.path.basename(inps.coherence_file)
        if check_isfile(outName):  inps.coherence_file = outName
        outName = 'Modified_'+os.path.basename(inps.spatial_coherence_file)
        if check_isfile(outName):  inps.spatial_coherence_file = outName


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print '\n**********  Reference in space  ***************'
    outName = 'Seeded_'+os.path.basename(inps.ifgram_file)
    if check_isfile(outName):
        inps.ifgram_file = outName
        print inps.ifgram_file + ' already exists, no need to re-seed.'
    else:
        print 'referncing all interferograms to the same pixel.'
        seedCmd = 'seed_data.py '+inps.ifgram_file+' -t '+inps.template_file+' -m '+inps.mask_file+\
                  ' -c '+inps.spatial_coherence_file
        if inps.geomap_file:
            seedCmd = seedCmd+' --trans '+inps.geomap_file
        print seedCmd
        os.system(seedCmd)
        inps.ifgram_file = outName


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    print '\n**********  Unwrapping Error Correction  **************'
    if 'pysar.unwrapError' in template.keys() and template['pysar.unwrapError'].lower() in ['y','yes']:
        outName = os.path.splitext(inps.ifgram_file)[0]+'_unwCor.h5'
        if check_isfile(outName):
            print outName+' exists, no need to re-process.'
        else:
            print 'This might take a while depending on the size of your data set!'
            unwCmd='unwrap_error.py -f '+inps.ifgram_file+' -m '+inps.mask_file
            print unwCmd
            os.system(unwCmd)
        inps.ifgram_file = outName
    else:  print 'No unwrapping error correction.'


    #########################################
    # Inversion of Interferograms
    ########################################
    print '\n**********  Network Inversion to Time Series  ********************'
    inps.timeseries_file = 'timeseries.h5'
    if check_isfile(inps.timeseries_file):
        print inps.timeseries_file+' already exists, inversion is not needed.'
    else:
        invertCmd = 'igram_inversion.py '+inps.ifgram_file
        print invertCmd
        os.system(invertCmd)

    ## Check DEM file for tropospheric delay setting
    ## DEM is needed with same coord (radar/geo) as timeseries file
    atr = readfile.read_attribute(inps.timeseries_file)
    if 'X_FIRST' in atr.keys():
        demFile = inps.dem_geo_file
    else:
        demFile = inps.dem_radar_file

    if not demFile or not check_isfile(demFile):
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print 'ERROR:'
        print '    DEM file was not found!'
        if 'pysar.troposphericDelay.method' in template.keys():
            print '    Continue without tropospheric correction ...'
            template.pop('pysar.troposphericDelay.method', None)
        if template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']:
            template.pop('pysar.deramp', None)
        print '++++++++++++++++++++++++++++++++++++++++++++++'


    ##############################################
    # Temporal Coherence: 
    #   A parameter to evaluate the consistency 
    #   of timeseries with the interferograms
    ##############################################
    print '\n********** Temporal Coherence file  *********'
    inps.temp_coherence_file = 'temporal_coherence.h5'
    if check_isfile(inps.temp_coherence_file):
        print inps.temp_coherence_file+' already exists.'
    else:
        tempCohCmd = 'temporal_coherence.py '+inps.ifgram_file+' '+inps.timeseries_file
        print tempCohCmd
        os.system(tempCohCmd)

    print '\nUpdate Mask based on Temporal Coherence ...'
    outName = 'Mask_tempCoh.h5'
    maskCmd = 'generate_mask.py -f '+inps.temp_coherence_file+' -m 0.7 -o '+outName
    print maskCmd
    os.system(maskCmd)
    inps.mask_file = outName


    ##############################################
    # Incident Angle
    ##############################################
    print '\n********** Incident Angle file  *************'
    inps.inc_angle_file = 'incidence_angle.h5'
    if check_isfile(inps.inc_angle_file):
        print inps.inc_angle_file+' already exists, no need to re-generate.'
    else:
        incAngleCmd = 'incidence_angle.py '+inps.timeseries_file
        print incAngleCmd
        os.system(incAngleCmd)


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   when Satellite is Envisat and
    #   Coordinate system is radar
    ############################################## 
    print '\n**********  Local Oscillator Drift correction for Envisat  ********'
    if atr['PLATFORM'].lower() in ['envisat','env']:
        outName = os.path.splitext(inps.timeseries_file)[0]+'_LODcor.h5'
        if check_isfile(outName):
            print inps.timeseries_file+' already exists.'
            inps.timeseries_file = outName
        else:
            if not 'X_FIRST' in atr.keys():
                LODcmd = 'lod.py '+inps.timeseries_file
                print LODcmd
                os.system(LODcmd)
                inps.timeseries_file = outName
            else:
                print 'WARNING: file is in geo coord, LOD correction is only supported for files in radar coord.'
                print 'Skip LOD correction.'
    else:
        print '\nLOD correction is not needed for '+atr['PLATFORM']+' data.'


    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print '\n**********  Tropospheric Delay Correction  ******************'
    if ('pysar.troposphericDelay.method' in template.keys() and 
        template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']):
        print '''
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        WARNING:
            Orbital error correction was BaseTropCor.
            Tropospheric correction was already applied simultaneous with baseline error correction.
            Tropospheric correction can not be applied again.
            To apply the tropospheric correction separated from baseline error correction, \
               choose other existing options for orbital error correction.
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        '''
        template.pop('pysar.troposphericDelay.method', None)
    
    if ('pysar.troposphericDelay.method' in template.keys() and 
        template['pysar.troposphericDelay.method'] not in ['no', 'auto']):
        trop_method = template['pysar.troposphericDelay.method']
        # Height-Correlation
        if trop_method in ['height_correlation']:
            print 'tropospheric delay correction with height-correlation approach'
            outName = os.path.splitext(inps.timeseries_file)[0]+'_tropHgt.h5'
            if check_isfile(outName):
                inps.timeseries_file = outName
                print inps.timeseries_file+' already exists.'
            else:
                try:
                    poly_order = template['pysar.troposphericDelay.polyOrder']
                except:
                    poly_order = '1'
                    print 'Deafult polynomial order for troposphreic correction = 1'
                cmdTrop = 'tropcor_phase_elevation.py'+' -f '+inps.timeseries_file+' -d '+\
                          demFile+' -p '+poly_order+' -m '+inps.mask_file
                print cmdTrop
                os.system(cmdTrop)
                inps.timeseries_file = outName
        
        # PyAPS
        elif trop_method == 'pyaps':
            print 'Atmospheric correction using Weather Re-analysis dataset (using PyAPS software)'
            model = template['pysar.troposphericDelay.weatherModel']
            print 'Weather Re-analysis dataset: '+model
            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+model+'.h5'
            if check_isfile(outName):
                inps.timeseries_file = outName
                print inps.timeseries_file+' already exists.'
            else:
                #acquisition_time = template['pysar.acquisitionTime']
                #print 'acquisition time: '+acquisition_time
                cmdTrop = 'tropcor_pyaps.py '+inps.timeseries_file+' -d '+demFile+' -s '+model
                print cmdTrop
                os.system(cmdTrop)
                inps.timeseries_file = outName
        else:
            print 'ERROR: Unrecognized atmospheric correction method: '+template['pysar.troposphericDelay.method']
    else:
        print 'No atmospheric delay correction.'


    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print '\n**********  Topographic Residual (DEM error) correction  *******'
    outName = os.path.splitext(inps.timeseries_file)[0]+'_demCor.h5'
    if 'pysar.topoError' in template.keys() and template['pysar.topoError'].lower() in ['y','yes','auto']:
        if check_isfile(outName):
            inps.timeseries_file = outName
            print inps.timeseries_file+' already exists.'
        else:
            print 'Correcting topographic residuals using method from Fattahi and Amelung, 2013, TGRS ...'
            topoCmd = 'dem_error.py '+inps.timeseries_file
            print topoCmd
            os.system(topoCmd)
            inps.timeseries_file = outName
    else:
        print 'No correction for topographic residuals.'


    ##############################################
    # Phase Ramp Correction (Optional)
    ##############################################
    print '\n**********  Ramp Removal  ***********************'
    if 'pysar.deramp' in template.keys() and template['pysar.deramp'] not in ['no','auto']:
        deramp_method = template['pysar.deramp']
        print 'Phase Ramp Removal method : '+deramp_method
        
        if deramp_method in ['plane', 'quadratic', 'plane_range', 'quadratic_range',\
                             'plane_azimuth', 'quadratic_azimuth']:
            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+deramp_method+'.h5'
            if check_isfile(outName):
                inps.timeseries_file = outName
                print inps.timeseries_file+' already exists.'
            else:
                derampCmd = 'remove_plane.py '+inps.timeseries_file+' -m '+inps.mask_file+' -s '+deramp_method
                print derampCmd
                os.system(derampCmd)
            inps.timeseries_file = outName

        elif deramp_method in ['baseline_cor','baselinecor']:
            outName = os.path.splitext(inps.timeseries_file)[0]+'_baselineCor.h5'
            if check_isfile(outName):
                inps.timeseries_file = outName
                print inps.timeseries_file+' already exists.'
            else:
                if not 'X_FIRST' in atr.keys():
                    derampCmd = 'baseline_error.py '+inps.timeseries_file+' '+inps.mask_file
                    print derampCmd
                    os.system(derampCmd)
                    inps.timeseries_file = outName
                else:
                    print 'WARNING!'
                    print 'Skipping correction.'
                    print 'baselineCor method can only be applied in radar coordinate'

        elif deramp_method in ['base_trop_cor','basetropcor','baselinetropcor']:
            outName = os.path.splitext(inps.timeseries_file)[0]+'_baseTropCor.h5'
            if check_isfile(outName):
                inps.timeseries_file = outName
                print inps.timeseries_file+' already exists.'
            else:
                if not 'X_FIRST' in atr.keys():
                    print 'Joint estimation of Baseline error and tropospheric delay'+\
                          ' [height-correlation approach]'
                    try:
                        poly_order = template['pysar.troposphericDelay.polyOrder']
                    except:
                        poly_order = '1'
                        print 'Deafult polynomial order for troposphreic correction = 1'
                    derampCmd = 'baseline_trop.py '+inps.timeseries_file+' '+inps.dem_radar_file+' '+\
                                poly_order+' range_and_azimuth'
                    print derampCmd
                    os.system(derampCmd)
                    inps.timeseries_file = outName
                else:
                    print 'WARNING!'
                    print 'Skipping correction.'
                    print 'baselineCor method can only be applied in radar coordinate'
        else:
            print 'WARNING: Unrecognized phase ramp method: '+template['pysar.deramp']
    else:
        print 'No phaes ramp removal.'


    #############################################
    # Velocity and rmse maps
    #############################################
    print '\n**********  Velocity estimation  **********************'
    velCmd = 'timeseries2velocity.py '+inps.timeseries_file
    print velCmd
    os.system(velCmd)
    inps.velocity_file = 'velocity.h5'


    ############################################
    # Post-processing
    # Geocodeing, masking and save to KML 
    ############################################
    print '\n**********  Post-processing  ********************************'
    if 'pysar.geocode' in template.keys() and template['pysar.geocode'].lower() in ['y','yes','auto']: 
        print '\ngeocoding ...\n'
        inps.geo_velocity_file       = check_geocode_file(inps.geomap_file, inps.velocity_file)
        inps.geo_temp_coherence_file = check_geocode_file(inps.geomap_file, inps.temp_coherence_file)
        inps.goe_timeseries_file     = check_geocode_file(inps.geomap_file, inps.timeseries_file)
        
        if inps.geo_velocity_file and inps.geo_temp_coherence_file:
            print 'masking geocoded velocity file: '+inps.geo_velocity_file+' ...'
            maskCmd = 'mask.py '+inps.geo_velocity_file+' -m '+inps.geo_temp_coherence_file+' -t 0.7'
            print maskCmd
            os.system(maskCmd)
            try: inps.geo_velocity_file = glob.glob(os.path.splitext(inps.geo_velocity_file)[0]+'_masked.h5')[0]
            except: pass
        
        if inps.geo_velocity_file:
            print 'creating Google Earth KMZ file for geocoded velocity file: '+inps.geo_velocity_file+' ...'
            kmlCmd = 'save_kml.py '+inps.geo_velocity_file
            print kmlCmd
            os.system(kmlCmd)
        else:
            print 'No geocoded velocity file found, skip creating KML file.'

    print 'masking velocity file: '+inps.velocity_file
    maskCmd = 'mask.py '+inps.velocity_file+' -m '+inps.mask_file
    print maskCmd
    os.system(maskCmd)
    inps.velocity_file = os.path.splitext(inps.velocity_file)[0]+'_masked.h5'


    #############################################
    #                PySAR v1.0                 #
    #############################################
    s = time.time()-start;  m, s = divmod(s, 60);  h, m = divmod(m, 60)
    print '\nTime used: %02d hours %02d mins %02d secs' % (h, m, s)
    print '\n###############################################'
    print 'End of PySAR processing!'
    print '################################################\n'



###########################################################################################
if __name__ == '__main__':
    main(sys.argv[:])

