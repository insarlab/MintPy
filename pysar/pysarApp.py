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
import warnings
import shutil

import h5py
import numpy as np

import pysar
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar.subset as subset
import pysar.load_data as load
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
    #inps.mask_file      = check_subset_file(inps.mask_file, vars(inps))
    inps.dem_radar_file = check_subset_file(inps.dem_radar_file, vars(inps))
    inps.coherence_file = check_subset_file(inps.coherence_file, vars(inps))
    #inps.spatial_coh_file = check_subset_file(inps.spatial_coh_file, vars(inps))

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
            pix_box = subset.bbox_geo2radar(geo_box, atr, geomap_file_orig)
        else:
            print 'use subset input in y/x'
            print 'calculate corresponding bounding box in geo coordinate.'
            geo_box = subset.bbox_radar2geo(pix_box, atr, geomap_file_orig)
        
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

pysar.topoError = yes        #[no], auto for yes
pysar.deramp    = plane      #[plane, quadratic, baseline_cor, base_trop_cor], auto for no
pysar.geocode   = yes        #[no], auto for yes

pysar.save.kml     = yes     #[no], auto for yes
pysar.save.geotiff = no      #[yes], auto for no, not implemented yet
pysar.save.unavco  = no      #[yes], auto for no
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
    inps.template_file = os.path.abspath(inps.template_file)
    inps.project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print 'Project name: '+inps.project_name

    # Work directory
    if not inps.work_dir:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ:
            inps.work_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+"/PYSAR"
            print 'Use file/dir structure in University of Miami.'+\
                  '(To turn it off, change miami_path value to False in pysar/__init__.py)'
        else:
            inps.work_dir = os.getcwd()
    inps.work_dir = os.path.abspath(inps.work_dir)

    if not os.path.isdir(inps.work_dir):
        os.makedirs(inps.work_dir)
    os.chdir(inps.work_dir)
    print "Go to work directory: "+inps.work_dir

    # Read template
    template = readfile.read_template(inps.template_file)
    for key in template.keys():
        if template[key].lower() == 'default':        template[key] = 'auto'
        if template[key].lower() in ['off','false']:  template[key] = 'no'
        if template[key].lower() in ['on','true']:    template[key] = 'yes'
    if 'pysar.deramp' in template.keys():
        template['pysar.deramp'] = template['pysar.deramp'].lower().replace('-','_')
    if 'pysar.troposphericDelay.method' in template.keys():
        template['pysar.troposphericDelay.method'] = template['pysar.troposphericDelay.method'].lower().replace('-','_')

    if not os.path.isfile(os.path.basename(inps.template_file)):
        shutil.copy2(inps.template_file, inps.work_dir)

    # Copy UNAVCO attribute txt file
    file_list = ['unavco_attributes.txt']
    try: inps.unavco_atr_file = ut.get_file_list(file_list)[0]
    except:
        try:
            process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
            shutil.copy2(process_dir+'/'+file_list[0], inps.work_dir)
            inps.unavco_atr_file = inps.work_dir+'/'+file_list[0]
        except:
            inps.unavco_atr_file = None
            print 'No UNAVCO attributes file found in PROCESS directory, skip copy'


    #########################################
    # Loading Data
    #########################################
    print '\n*************** Load Data ****************'
    #inps = load.load_data_from_template(inps.template_file, inps)
    loadCmd = 'load_data.py --template '+inps.template_file+' --dir '+inps.work_dir
    print loadCmd
    os.system(loadCmd)
    os.chdir(inps.work_dir)

    print '--------------------------------------------'
    ## 1. Check loading result
    # Required files
    try:
        inps.ifgram_file = ut.get_file_list([inps.work_dir+'/unwrapIfgram.h5', inps.work_dir+'/LoadedData.h5'], abspath=True)[0]
        atr = readfile.read_attribute(inps.ifgram_file)
        print 'Unwrapped interferograms: '+inps.ifgram_file
    except RuntimeError:
        print '\nNo interferograms file found!\n'

    # Recommended files (None if not found)
    # Spatial coherence for each interferogram
    try:
        inps.coherence_file = ut.get_file_list([inps.work_dir+'/coherence.h5', inps.work_dir+'/Coherence.h5'], abspath=True)[0]
        print 'Coherences: '+inps.coherence_file
    except:
        inps.coherence_file = None
        warnings.warn('No coherences file found. Cannot use coherence-based network modification without it.')

    # DEM in radar coord
    file_list = ['demRadar.h5','radar*.hgt']
    try: file_list.append(os.path.basename(template['pysar.dem.geoCoord']))
    except: pass
    try:
        inps.dem_radar_file = ut.get_file_list(file_list, abspath=True)[0]
        print 'DEM in radar coord: '+str(inps.dem_radar_file)
    except:
        inps.dem_radar_file = None
        if not 'Y_FIRST' in atr.keys():
            warnings.warn('No radar coord DEM found.')

    # DEM in geo coord
    file_list = ['demGeo.h5','*.dem']
    try: file_list.append(os.path.basename(template['pysar.dem.geoCoord']))
    except: pass
    try:
        inps.dem_geo_file = ut.get_file_list(file_list, abspath=True)[0]
        print 'DEM in geo   coord: '+inps.dem_geo_file
    except:
        inps.dem_geo_file = None
        warnings.warn('No geo coord DEM found.')

    # Transform file for geocoding
    file_list = ['geomap*.trans']
    try: file_list.append(os.path.basename(template['pysar.geomap']))
    except: pass
    try:
        inps.geomap_file = ut.get_file_list(file_list, abspath=True)[0]
        print 'Transform file: '+str(inps.geomap_file)
    except:
        inps.geomap_file = None
        if not 'Y_FIRST' in atr.keys():
            warnings.warn('No transform file found! Can not geocoding without it!')


    #########################################
    # Check the subset (Optional)
    #########################################
    print '\n*************** Subset ****************'
    if inps.geomap_file:
        outName = os.path.splitext(inps.geomap_file)[0]+'_tight'+os.path.splitext(inps.geomap_file)[1]
        # Get bounding box of non-zero area in geomap*.trans file
        trans_rg, trans_atr = readfile.read(inps.geomap_file, (), 'range')
        idx_row, idx_col = np.nonzero(trans_rg)
        pix_box = (np.min(idx_col)-10, np.min(idx_row)-10, np.max(idx_col)+10, np.max(idx_row)+10)
        # Subset geomap_file only if it could save > 20% percent of area
        if abs((pix_box[2]-pix_box[0])*(pix_box[3]-pix_box[1])) < 0.8*(trans_rg.shape[0]*trans_rg.shape[1]):
            print "Get tight subset of geomap*.trans file and/or DEM file in geo coord"
            print '--------------------------------------------'
            inps = subset.subset_box2inps(inps, pix_box, None)
            inps.fill_value = 0.0
            inps.geomap_file = check_subset_file(inps.geomap_file, vars(inps), outName)
            
            # Subset DEM in geo coord
            outName = os.path.splitext(inps.dem_geo_file)[0]+'_tight'+os.path.splitext(inps.dem_geo_file)[1]
            geomap_atr = readfile.read_attribute(inps.geomap_file)
            pix_box, geo_box = subset.get_coverage_box(geomap_atr)
            inps = subset.subset_box2inps(inps, pix_box, geo_box)
            inps.dem_geo_file = check_subset_file(inps.dem_geo_file, vars(inps), outName, overwrite=True)

    # Subset based on input template
    if [key for key in template.keys()\
        if ('pysar.subset' in key and template[key].lower() not in ['auto','no',''])]:
        # Read subset option from template file, and return None if lalo/yx is not specified.
        pix_box, geo_box = subset.read_subset_template2box(inps.template_file)
        inps = create_subset_dataset(inps, pix_box, geo_box)
    else:
        print '\nNo Subset selected. Processing the whole area.'


    #########################################
    # Generating Aux files
    #########################################
    print '\n*************** Generate Initla Mask ****************'
    # Initial mask
    inps.mask_file = 'mask.h5'
    if ut.update_file(inps.mask_file, inps.ifgram_file):
        print 'creating mask file using non-zero pixels from file: '+inps.ifgram_file
        inps.mask_file = ut.nonzero_mask(inps.ifgram_file, inps.mask_file)

    # Average spatial coherence
    if inps.coherence_file:
        inps.spatial_coh_file = 'averageSpatialCoherence.h5'
        if ut.update_file(inps.spatial_coh_file, inps.coherence_file):
            inps.spatial_coh_file = ut.temporal_average(inps.coherence_file, inps.spatial_coh_file)
    else:
        inps.spatial_coh_file = None


    #########################################
    # Network Modification (Optional)
    #########################################
    if [key for key in template.keys() if 'pysar.network' in key]:
        print '\n*************** Modify Network ****************'
        #outName = 'Modified_'+os.path.basename(inps.ifgram_file)
        #if ut.update_file(outName, inps.ifgram_file):
        networkCmd = 'modify_network.py --template '+inps.template_file+' --mask '+inps.mask_file+\
                     ' --plot --mark-attribute '+inps.ifgram_file
        if inps.coherence_file:
            networkCmd += ' '+inps.coherence_file
        print networkCmd
        os.system(networkCmd)
        #if not ut.update_file(outName):
        #    inps.ifgram_file = outName
        
        outName = 'Modified_'+os.path.basename(inps.mask_file)
        if not ut.update_file(outName):
            inps.mask_file = outName
        
        if inps.spatial_coh_file:
            outName = 'Modified_'+os.path.basename(inps.spatial_coh_file)
            if not ut.update_file(outName):
                inps.spatial_coh_file = outName
        #if inps.coherence_file:
        #    outName = 'Modified_'+os.path.basename(inps.coherence_file)
        #    if not ut.update_file(outName):
        #        inps.coherence_file = outName


    #########################################
    # Referencing Interferograms in Space
    #########################################
    print '\n**********  Reference in space  ***************'
    atr = readfile.read_attribute(inps.ifgram_file)
    try:
        ref_x = int(atr['ref_x'])
        ref_y = int(atr['ref_y'])
        print 'find reference pixel in y/x: [%d, %d], skip updating.'%(ref_y, ref_x)
    except:
        print 'call seed_data.py to find reference pixel in space'
        seedCmd = 'seed_data.py '+inps.ifgram_file+' -t '+inps.template_file+' -m '+inps.mask_file+' --mark-attribute'
        if inps.spatial_coh_file:
            seedCmd += ' -c '+inps.spatial_coh_file
        if inps.geomap_file:
            seedCmd += ' --trans '+inps.geomap_file
        print seedCmd
        os.system(seedCmd)


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms
    ############################################
    print '\n**********  Unwrapping Error Correction  **************'
    if 'pysar.unwrapError' in template.keys() and template['pysar.unwrapError'].lower() in ['y','yes']:
        outName = os.path.splitext(inps.ifgram_file)[0]+'_unwCor.h5'
        if ut.update_file(outName, inps.ifgram_file):
            print 'This might take a while depending on the size of your data set!'
            unwCmd='unwrap_error.py -f '+inps.ifgram_file+' -m '+inps.mask_file
            print unwCmd
            os.system(unwCmd)
        inps.ifgram_file = outName
    else:
        print 'No unwrapping error correction.'


    #########################################
    # Inversion of Interferograms
    ########################################
    print '\n**********  Network Inversion to Time Series  ********************'
    inps.timeseries_file = 'timeseries.h5'
    if ut.update_file(inps.timeseries_file, inps.ifgram_file):
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

    if ut.update_file(demFile):
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
    inps.temp_coh_file = 'temporalCoherence.h5'
    if ut.update_file(inps.temp_coh_file, inps.timeseries_file):
        tempCohCmd = 'temporal_coherence.py '+inps.ifgram_file+' '+inps.timeseries_file+' '+inps.temp_coh_file
        print tempCohCmd
        os.system(tempCohCmd)

    print '\n--------------------------------------------'
    print 'Update Mask based on Temporal Coherence ...'
    outName = 'maskTempCoh.h5'
    maskCmd = 'generate_mask.py -f '+inps.temp_coh_file+' -m 0.7 -o '+outName
    print maskCmd
    os.system(maskCmd)
    inps.mask_file = outName


    ###############################################
    ## Incident Angle
    ###############################################
    #print '\n********** Incident Angle file  *************'
    #inps.inc_angle_file = 'incidenceAngle.h5'
    #if ut.update_file(inps.inc_angle_file, inps.timeseries_file):
    #    incAngleCmd = 'incidence_angle.py '+inps.timeseries_file+' '+inps.inc_angle_file
    #    print incAngleCmd
    #    os.system(incAngleCmd)


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   when Satellite is Envisat and
    #   Coordinate system is radar
    ############################################## 
    print '\n**********  Local Oscillator Drift correction for Envisat  ********'
    sar_mission = atr['PLATFORM'].lower()
    if sar_mission.startswith('env') and 'Y_FIRST' not in atr.keys():
        outName = os.path.splitext(inps.timeseries_file)[0]+'_LODcor.h5'
        if ut.update_file(outName, inps.timeseries_file):
            LODcmd = 'lod.py '+inps.timeseries_file
            print LODcmd
            os.system(LODcmd)
        inps.timeseries_file = outName
    else:
        if not sar_mission.startswith('env'):
            print '\nLOD correction is not needed for '+sar_mission+' data.'
        else:
            warnings.warn('Can not apply LOD correction for file in radar coord.')

    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print '\n**********  Tropospheric Delay Correction  ******************'
    if ('pysar.troposphericDelay.method' in template.keys() and 
        template['pysar.deramp'] in ['base_trop_cor','basetropcor','baselinetropcor']):
        message='''
        Orbital error correction was BaseTropCor.
        Tropospheric correction was already applied simultaneous with baseline error correction.
        Tropospheric correction can not be applied again.
        To apply the tropospheric correction separated from baseline error correction, \
           choose other existing options for orbital error correction.
        '''
        warnings.warn(message)
        template.pop('pysar.troposphericDelay.method', None)
    
    if ('pysar.troposphericDelay.method' in template.keys() and 
        template['pysar.troposphericDelay.method'] not in ['no', 'auto']):
        trop_method = template['pysar.troposphericDelay.method']
        # Height-Correlation
        if trop_method in ['height_correlation']:
            print 'tropospheric delay correction with height-correlation approach'
            outName = os.path.splitext(inps.timeseries_file)[0]+'_tropHgt.h5'
            if ut.update_file(outName, inps.timeseries_file):
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
            try:    model = template['pysar.troposphericDelay.weatherModel']
            except: model = 'ECMWF'
            print 'Weather Re-analysis dataset: '+model
            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+model+'.h5'
            if ut.update_file(outName, inps.timeseries_file):
                #acquisition_time = template['pysar.acquisitionTime']
                cmdTrop = 'tropcor_pyaps.py '+inps.timeseries_file+' -d '+demFile+' -s '+model+\
                          ' --weather-dir '+inps.work_dir+'/../WEATHER'
                print cmdTrop
                os.system(cmdTrop)
            inps.timeseries_file = outName
        else:
            sys.exit('Unrecognized atmospheric correction method: '+template['pysar.troposphericDelay.method'])
    else:
        print 'No atmospheric delay correction.'


    ##############################################
    # Topographic (DEM) Residuals Correction (Optional)
    ##############################################
    print '\n**********  Topographic Residual (DEM error) correction  *******'
    outName = os.path.splitext(inps.timeseries_file)[0]+'_demErr.h5'
    if 'pysar.topoError' in template.keys() and template['pysar.topoError'].lower() in ['y','yes','auto']:
        if ut.update_file(outName, inps.timeseries_file):
            print 'Correcting topographic residuals using method from Fattahi and Amelung, 2013, TGRS ...'
            topoCmd = 'dem_error.py '+inps.timeseries_file+' -o '+outName
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
            if ut.update_file(outName, inps.timeseries_file):
                derampCmd = 'remove_plane.py '+inps.timeseries_file+' -m '+inps.mask_file+' -s '+deramp_method
                print derampCmd
                os.system(derampCmd)
            inps.timeseries_file = outName

        elif deramp_method in ['baseline_cor','baselinecor']:
            if not 'X_FIRST' in atr.keys():
                outName = os.path.splitext(inps.timeseries_file)[0]+'_baselineCor.h5'
                if ut.update_file(outName, inps.timeseries_file):
                    derampCmd = 'baseline_error.py '+inps.timeseries_file+' '+inps.mask_file
                    print derampCmd
                    os.system(derampCmd)
                inps.timeseries_file = outName
            else:
                warnings.warn('BaselineCor method can only be applied in radar coordinate, skipping correction')

        elif deramp_method in ['base_trop_cor','basetropcor','baselinetropcor']:
            if not 'X_FIRST' in atr.keys():
                print 'Joint estimation of Baseline error and tropospheric delay [height-correlation approach]'
                outName = os.path.splitext(inps.timeseries_file)[0]+'_baseTropCor.h5'
                if ut.update_file(outName, inps.timeseries_file):
                    try:    poly_order = template['pysar.troposphericDelay.polyOrder']
                    except: poly_order = '1'
                    derampCmd = 'baseline_trop.py '+inps.timeseries_file+' '+inps.dem_radar_file+' '+\
                                poly_order+' range_and_azimuth'
                    print derampCmd
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
    inps.velocity_file = 'velocity.h5'
    if ut.update_file(inps.velocity_file, inps.timeseries_file):
        velCmd = 'timeseries2velocity.py '+inps.timeseries_file
        print velCmd
        os.system(velCmd)


    ############################################
    # Post-processing
    # Geocodeing, masking and save to KML 
    ############################################
    print '\n**********  Post-processing  ********************************'
    if 'Y_FIRST' in atr.keys() and 'pysar.geocode' in template.keys() and template['pysar.geocode'].lower() in ['y','yes','auto']:
        template['pysar.geocode'] = 'no'
        print 'dataset is in geo coordinate, no need for geocoding.'

    if 'pysar.geocode' in template.keys() and template['pysar.geocode'].lower() in ['y','yes','auto']: 
        print '\ngeocoding ...\n'
        inps.geo_velocity_file   = check_geocode_file(inps.geomap_file, inps.velocity_file)
        inps.geo_temp_coh_file   = check_geocode_file(inps.geomap_file, inps.temp_coh_file)
        inps.goe_timeseries_file = check_geocode_file(inps.geomap_file, inps.timeseries_file)
        
        if inps.geo_velocity_file and inps.geo_temp_coh_file:
            print 'masking geocoded velocity file: '+inps.geo_velocity_file+' ...'
            maskCmd = 'mask.py '+inps.geo_velocity_file+' -m '+inps.geo_temp_coh_file+' -t 0.7'
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
    outName = os.path.splitext(inps.velocity_file)[0]+'_masked.h5'
    if ut.update_file(outName, inps.velocity_file):
        maskCmd = 'mask.py '+inps.velocity_file+' -m '+inps.mask_file
        print maskCmd
        os.system(maskCmd)
    inps.velocity_file = outName


    #############################################
    # Save to UNAVCO InSAR Archive format
    #############################################
    if 'pysar.save.unavco' in template.keys() and template['pysar.save.unavco'].lower() in ['y','yes']:
        print '\n*********  Output to UNAVCO InSAR Archive Format  ***********'
        if not inps.geomap_file or not inps.dem_geo_file:
            warnings.warn('No geomap*.tran file or DEM in geo coord found! Skip saving.')
        else:
            inps.geo_timeseries_file = check_geocode_file(inps.geomap_file, inps.timeseries_file)
            # Add UNAVCO attributes
            if inps.unavco_atr_file:
                atrCmd = 'add_attribute.py '+inps.geo_timeseries_file+' '+inps.unavco_atr_file
                print atrCmd
                os.system(atrCmd)

            inps.unavco_file = unavco.get_unavco_filename(inps.geo_timeseries_file)
            if ut.update_file(inps.unavco_file, inps.geo_timeseries_file):
                # Temporal Coherence and Mask
                inps.geo_temp_coh_file = check_geocode_file(inps.geomap_file, inps.temp_coh_file)
                inps.geo_mask_file = 'geo_maskTempCoh.h5'
                # Mask file in geo coord
                if inps.geo_temp_coh_file and ut.update_file(inps.geo_mask_file, inps.geo_temp_coh_file):
                    maskCmd = 'generate_mask.py -f '+inps.geo_temp_coh_file+' -m 0.7 -o '+inps.geo_mask_file
                    print maskCmd
                    os.system(maskCmd)
                # Incidence Angle
                inps.inc_angle_file = 'incidenceAngle.h5'
                if ut.update_file(inps.inc_angle_file, inps.timeseries_file):
                    incAngleCmd = 'incidence_angle.py '+inps.timeseries_file+' '+inps.inc_angle_file
                    print incAngleCmd
                    os.system(incAngleCmd)
                inps.geo_inc_angle_file = check_geocode_file(inps.geomap_file, inps.inc_angle_file)

                # Save to UNAVCO format
                unavcoCmd = 'save_unavco.py '+inps.geo_timeseries_file+' -d '+inps.dem_geo_file+\
                            ' -i '+inps.geo_inc_angle_file+' -c '+inps.geo_temp_coh_file
                print unavcoCmd
                os.system(unavcoCmd)


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

