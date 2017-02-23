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


##################################### Sub Functions #######################################

###########################  Usage Function  ###########################
#generate_from: http://patorjk.com/software/taag/
PYSAR_LOGO='''
_________________________________________________
       ____             __     __     ____  
       /    )         /    )   / |    /    )
------/____/----------\-------/__|---/___ /------
     /        /   /    \     /   |  /    |  
____/________(___/_(____/___/____|_/_____|_______
                /                           
            (_ /                            

 A Python package for InSAR time series analysis.
               PySAR v1.0, Jan 2017
 Geodesy Lab, University of Miami, Maimi FL, USA
_________________________________________________
'''

PYSAR_LOGO_OLD='''
*******************************************************
*******************************************************
*******************************************************
*******************************************************
*********   OOOOO      OOOOO     O     OOOO   *********
*********   O   O O O  O        O O    O   O  *********
*********   OOOOO OOO  OOOOO   OOOOO   OOOO   *********  
*********   O       O      O  O     O  O  O   *********
*********   O     OOO  OOOOO O       O O   O  *********
*********                                     *********
*******************************************************
*******************************************************
*******************************************************


*******************************************************
A Python Module for InSAR time-series analysis.
PySAR v1.0 July 2013, Geodesy Lab, RSMAS, University of Miami
'''

TEMPLATE='''template:
  pysar.unwrapFiles    = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*.unw
  pysar.corFiles       = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.cor
  pysar.wrapFiles      = /SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.int     #optional
  pysar.geomap         = /SanAndreasT356EnvD/PROCESS/GEO/*050102-070809*/geomap*.trans
  pysar.dem.radarCoord = /SanAndreasT356EnvD/PROCESS/DONE/*050102-070809*/radar*.hgt
  pysar.dem.geoCoord   = /SanAndreasT356EnvD/DEM/srtm1_30m.dem                   #optional
  
  pysar.network.dropIfgramIndex = 7:9 15 25 26        #start from 1
  pysar.network.dropDate        = 20080520 20090816
  pysar.network.maxTempBaseline = 720
  pysar.network.maxPerpBaseline = 2000
  pysar.network.reference       = Modified_unwrapIfgram.h5
  pysar.network.reference       = Paris.list
  pysar.network.coherenceBase   = yes    #search and use coherence from input
  
  pysar.subset.yx          = 1800:2000,700:800        #optional
  pysar.subset.lalo        = 31.5:32.5,130.5:131.0    #optional, priority: lalo > yx
  
  pysar.reference.yx       = 257 , 151                #optional
  pysar.reference.lalo     = 31.8, 130.8              #[not implemented yet]
  pysar.reference.date     = 20090120
   
  pysar.troposphericDelay.method        = pyaps   #['height_correlation'] 
  pysar.troposphericDelay.polyOrder     = 1       #for 'height_correlation' method
  pysar.troposphericDelay.weatherModel  = ECMWF   #['ERA','MERRA', 'NARR'], for 'pyaps' method
  
  pysar.topoError = yes               #['no'], optional
  
  pysar.deramp    = plane             #[    'plane',     'plane_range',     'plane_azimuth']
                                      #['quadratic', 'quadratic_range', 'quadratic_azimuth']
                                      #['baseline_cor','base_trop_cor']
  
  pysar.geocode   = yes               #['no'], optional
'''

EXAMPLE='''example:
  pysarApp.py  SanAndreasT356EnvD.template
  pysarApp.py  SanAndreasT356EnvD.template  --dir ~/insarlab/SanAndreasT356EnvD/TSSAR
'''

UM_FILE_STRUCT='''
    scratch/                 # $SCRATCHDIR defined in environmental variable
        SanAndreasT356EnvD/  # my_projectName, same as the basename of template file
            DEM/             # DEM file(s) (for topographic phase and geocode)
            DOWNLOAD/        # (optional) Data downloaded from agencies
            PROCESS/         # Interferograms processed by ROI_PAC, Gamma, ISCE, ... 
            RAW/             # (optional) Raw SAR data untared from DOWNLOAD directory
            SLC/             # (optional) SLC SAR data after focusing from RAW directory
            TSSAR/           # PySAR work directory for time series analysis
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description=PYSAR_LOGO,
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)
    
    parser.add_argument('template_file', help='template for PySAR setting')
    parser.add_argument('--dir', dest='work_dir',\
                        help='specify working directory; default is:\n'+\
                             '$SCRATCHDIR/my_projectName/TSSAR\n'+\
                             '    If using University of Miami file/dir structure\n'+\
                             '    To turn it off, change miami_path value to False in pysar/__init__.py\n'
                             './ - current directory\n'+\
                             '    To set this as your default permanetly,'+\
                             ' change miami_path value to False in pysar/__init__.py')
    parser.add_argument('-v','--version', action='version', version='%(prog)s 2.0')

    inps = parser.parse_args()
    return inps


####################################  Main Function  ######################################
def main(argv):
    start = time.time()
    inps = cmdLineParse()


    #########################################
    # Initiation
    #########################################
    print PYSAR_LOGO
    inps.project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print 'Project name: '+inps.project_name
    
    # work directory
    if not inps.work_dir:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ:
            inps.work_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+"/TSSAR"
            print 'Use file/dir structure in University of Miami.'+\
                  '(To turn it off, change miami_path value to False in pysar/__init__.py)'
        else:
            inps.work_dir = os.getcwd()
    if not os.path.isdir(inps.work_dir):   os.mkdir(inps.work_dir)
    os.chdir(inps.work_dir)
    print "Work directory: "+inps.work_dir+'\nGo to work directory'
    
    # Read template
    template = readfile.read_template(inps.template_file)


    #########################################
    # Loading Data
    #########################################
    loadCmd = 'load_data.py '+inps.template_file
    print loadCmd
    os.system(loadCmd)

    ##### Initial files name
    # HDF5 files
    inps.ifgram_file = 'unwrapIfgram.h5'
    inps.coherence_file = 'coherence.h5'
    inps.mask_file = 'Mask.h5'
    inps.spatial_coherence_file = 'average_spatial_coherence.h5'
    
    print '--------------------------------------------'
    if os.path.isfile(inps.ifgram_file):  print 'Unwrapped interferograms: '+inps.ifgram_file
    else:  sys.exit('\nERROR: No interferograms file found!\n')

    if os.path.isfile(inps.coherence_file):  print 'Coherences: '+inps.coherence_file
    else:
        print '\nWARNING: No coherences file found. '+\
              'Cannot use coherence-based network modification without it.\n'

    if not os.path.isfile(inps.mask_file):
        print 'No mask file found. Creating one using non-zero pixels in file: '+inps.ifgram_file
        inps.mask_file = ut.nonzero_mask(inps.ifgram_file, inps.mask_file)
    print 'Mask: '+inps.mask_file

    # DEM in geo coord
    try:    inps.dem_geo_file = os.path.basename(template['pysar.dem.geoCoord'])
    except: inps.dem_geo_file = '*.dem'
    inps.dem_geo_file = glob.glob('./'+inps.dem_geo_file)[0]
    if inps.dem_geo_file:    print 'DEM in geo coord: '+str(inps.dem_geo_file)
    else:  print '\nWARNING: No geo coord DEM found.\n'
    
    # DEM in radar coord
    try:    inps.dem_radar_file = os.path.basename(template['pysar.dem.radarCoord'])
    except: inps.dem_radar_file = 'radar*.hgt'
    inps.dem_radar_file = glob.glob('./'+inps.dem_radar_file)[0]
    if inps.dem_radar_file:  print 'DEM in radar coord: '+str(inps.dem_radar_file)
    else:
        print '\nWARNING: No radar coord DEM found! '+\
              'Cannot use tropospheric delay correction without it.\n'
    
    # Transform file for geocoding
    try:    inps.geomap_file = os.path.basename(template['pysar.geomap'])
    except: inps.geomap_file = 'geomap*.trans'
    inps.geomap_file = glob.glob('./'+inps.geomap_file)[0]
    if inps.geomap_file:     print 'Transform file: '+str(inps.geomap_file)
    else:  print '\nWARNING: No transform file found! Cannot geocoding without it.\n'


    #########################################
    # Network Modification (Optional)
    #########################################
    if [key for key in template.keys() if 'pysar.network' in key]:
        outName = 'Modified_'+inps.ifgram_file
        if os.path.isfile(outName):
            print '\n'+outName+' already existed, no need to re-modify network.\n'
            inps.ifgram_file = outName
        else:
            networkCmd = 'modify_network.py '+inps.ifgram_file+' '+inps.coherence_file+\
                         ' --template '+inps.template_file+' --mask '+inps.mask_file+' --plot'
            print networkCmd
            os.system(networkCmd)

    if os.path.isfile('Modified_'+inps.mask_file):
        inps.mask_file = 'Modified_'+inps.mask_file
    if os.path.isfile('Modified_'+inps.coherence_file):
        inps.coherence_file = 'Modified_'+inps.coherence_file
    if os.path.isfile('Modified_'+inps.spatial_coherence_file):
        inps.spatial_coherence_file = 'Modified_'+inps.spatial_coherence_file


    #########################################
    # Check the subset (Optional)
    #########################################
    # subset geomap*.trans file with --footprint option
    if inps.geomap_file:
        outFile = os.path.splitext(inps.geomap_file)[0]+'_tight'+os.path.splitext(inps.geomap_file)[1]
        subsetCmd = 'subset.py '+inps.geomap_file+' --footprint '+' -o '+outFile
        print subsetCmd
        os.system(subsetCmd)
        inps.geomap_file = outFile
    
    # Subset based on input template
    if [key for key in template.keys() if 'pysar.subset' in key]:
        pix_box, geo_box = subset.read_subset_template2box(inps.template_file)
        if geo_box:
            pix_box = None
            print '--------------------------------------------'
            print 'subseting data in UpperLeft/LowerRight lon/lat: '+str(geo_box)
            inps = subset.subset_box2inps(inps, pix_box, geo_box)
            inps.dem_geo_file = subset.subset_file(inps.dem_geo_file, vars(inps))
            inps.geomap_file  = subset.subset_file(inps.geomap_file, vars(inps))
            
            print '--------------------------------------------'
            print 'calculate bounding box in radar coordinate.'
            lat = np.array([geo_box[3],geo_box[3],geo_box[1],geo_box[1]])
            lon = np.array([geo_box[0],geo_box[2],geo_box[0],geo_box[2]])
            x,y,x_res,y_res = ut.glob2radar(lat, lon, geomapFile=inps.geomap_file)
            buf = 10*(np.max([x_res, y_res]))
            pix_box = (np.min(x)-buf, np.min(y)-buf, np.max(x)+buf, np.max(y)+buf)
            geo_box = None
            print 'subset data in UpperLeft/LowerRight x/y: '+str(pix_box)
            inps = subset.subset_box2inps(inps, pix_box, geo_box)
            inps.ifgram_file    = subset.subset_file(inps.ifgram_file, vars(inps))
            inps.coherence_file = subset.subset_file(inps.coherence_file, vars(inps))
            inps.mask_file      = subset.subset_file(inps.mask_file, vars(inps))

        elif pix_box:
            geo_box = None
            print '--------------------------------------------'
            print 'subset data in UpperLeft/LowerRight x/y: '+str(pix_box)
            inps = subset.subset_box2inps(inps, pix_box, geo_box)
            inps.ifgram_file    = subset.subset_file(inps.ifgram_file, vars(inps))
            inps.coherence_file = subset.subset_file(inps.coherence_file, vars(inps))
            inps.mask_file      = subset.subset_file(inps.mask_file, vars(inps))

            print '--------------------------------------------'
            print 'calculating bounding box in geo coordinate.'
            x = np.array([pix_box[0],pix_box[2],pix_box[0],pix_box[2]])
            y = np.array([pix_box[1],pix_box[1],pix_box[3],pix_box[3]])
            lat,lon,lat_res,lon_res = ut.radar2glob(x,y,inps.ifgram_file,1)
            buf = 10*(np.max([lat_res,lon_res]))
            geo_box = (np.min(lon)-buf, np.max(lat)+buf, np.max(lon)+buf, np.min(lat)-buf)
            pix_box = None
            print 'subset data in UpperLeft/LowerRight lon/lat: '+str(geo_box)
            inps = subset.subset_box2inps(inps, pix_box, geo_box)
            inps.dem_geo_file = subset.subset_file(inps.dem_geo_file, vars(inps))
            inps.geomap_file  = subset.subset_file(inps.geomap_file, vars(inps))  
    else:
        print 'No Subset selected. Processing the whole area.'


    #########################################
    # Referencing Interferograms in Space
    #########################################
    #print '\n**********  Referencing Interferograms  ***************'
    if os.path.isfile('Modified_Seeded_'+inps.ifgram_file):
        inps.ifgram_file = 'Modified_Seeded_'+inps.ifgram_file
        print inps.ifgram_file + ' already exists.'
    elif os.path.isfile('Seeded_'+inps.ifgram_file):
        inps.ifgram_file = 'Seeded_'+inps.ifgram_file
        print inps.ifgram_file + ' already exists.'
    else:
        print 'referncing all interferograms to the same pixel.'
        seedCmd = 'seed_data.py '+inps.ifgram_file+' -t '+inps.template_file+' -m '+inps.mask_file+\
                  ' -c '+inps.spatial_coherence_file
        print seedCmd
        os.system(seedCmd)
        inps.ifgram_file = 'Seeded_'+inps.ifgram_file


    ############################################
    # Unwrapping Error Correction (Optional)
    #    based on the consistency of triplets
    #    of interferograms  
    ############################################
    print '\n**********  Unwrapping Error Correction  **************'
    outName = os.path.splitext(inps.ifgram_file)[0]+'_unwCor.h5'
    try:
        template['pysar.unwrapError']
        if template['pysar.unwrapError'].lower() in ['y','yes']:
            if os.path.isfile(outName):
                inps.ifgram_file = outName
                print inps.ifgram_file+' exists.'
            else:
                print 'This might take a while depending on the size of your data set!'
                unwCmd='unwrap_error.py -f '+inps.ifgram_file+' -m '+inps.mask_file
                os.system(unwCmd)
                inps.ifgram_file = outName
        else:  print 'No unwrapping error correction.'
    except:  print 'No unwrapping error correction.'


    #########################################
    # Inversion of Interferograms
    ########################################
    #print '\n**********  Time Series Inversion  ********************'
    if os.path.isfile('timeseries.h5'):
        print 'timeseries.h5 already exists, inversion is not needed.'
    else:
        invertCmd = 'igram_inversion.py '+inps.ifgram_file
        print invertCmd
        os.system(invertCmd)
    inps.timeseries_file = 'timeseries.h5'
    atr = readfile.read_attribute(inps.timeseries_file)


    ##############################################
    # Temporal Coherence: 
    #   A parameter to evaluate the consistency 
    #   of timeseries with the interferograms
    ##############################################
    #print '\n**********  Generate Temporal Coherence file  *********'
    if os.path.isfile('temporal_coherence.h5'):
        print 'temporal_coherence.h5 already exists.'
    else:
        tempCohCmd='temporal_coherence.py '+inps.ifgram_file+' '+inps.timeseries_file
        print tempCohCmd
        os.system(tempCohCmd)
    inps.temp_coherence_file = 'temporal_coherence.h5'


    ##############################################
    # Update Mask based on temporal coherence (Optional)
    ##############################################
    print '\n**********  Update Mask based on Temporal Coherence  **'
    outName = 'Mask_tempCoh.h5'
    maskCmd='generate_mask.py -f '+inps.temp_coherence_file+' -m 0.7 -o '+outName
    print maskCmd
    os.system(maskCmd)
    inps.mask_file = outName


    ##############################################
    # Incident Angle
    ##############################################
    print '\n**********  Generate Incident Angle file  *************'
    incAngleCmd = 'incidence_angle.py '+inps.timeseries_file
    print incAngleCmd
    os.system(incAngleCmd)
    inps.inc_angle_file = 'incidence_angle.h5'

    ### geoCoord
    #rdr_or_geo = ut.radar_or_geo(inps.timeseries_file)


    ##############################################
    # LOD (Local Oscillator Drift) Correction
    #   when Satellite is Envisat and
    #   Coordinate system is radar
    ############################################## 
    print '\n**********  Local Oscillator Drift correction  ********'
    outName = os.path.splitext(inps.timeseries_file)[0]+'_LODcor.h5'
    if os.path.isfile(outName):
        print inps.timeseries_file+' already exists.'
        inps.timeseries_file = outName
    else:
        platform = atr['PLATFORM']
        if platform in ['envisat','env'] and 'X_FIRST' in atr.keys():
            LODcmd='lod.py '+inps.timeseries_file
            print LODcmd
            os.system(LODcmd)
            inps.timeseries_file = outName
        else:
            print 'No need of LOD correction for '+platform+' data.'
            print 'Only Envisat data in radar coord should be corrected for LOD.'


    ##############################################
    # Tropospheric Delay Correction (Optional)
    ##############################################
    print '\n**********  Tropospheric Correction  ******************'
    # Check conflicts
    if 'pysar.troposphericDelay.method' in template.keys():
        deramp_method = template['pysar.deramp'].lower().replace('-','_')
        # 1. Conflict with Base-Trop ramp removal
        if deramp_method in ['base_trop_cor']:
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

        # 2. If DEM is not existed
        if 'X_FIRST' in atr.keys():
            demFile = inps.dem_geo_file
        else:
            demFile = inps.dem_radar_file
        if not os.path.isfile(demFile):
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print 'ERROR:'
            print '    DEM file was not found!'
            print '    Continue without tropospheric correction ...'
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            template.pop('pysar.troposphericDelay.method', None)
    
    if 'pysar.troposphericDelay.method' in template.keys():
        trop_method = template['pysar.troposphericDelay.method'].lower().replace('-','_')
        # Height-Correlation
        if trop_method in ['height_correlation']:
            print 'tropospheric delay correction with height-correlation approach'
            outName = os.path.splitext(inps.timeseries_file)[0]+'_tropHgt.h5'
            if os.path.isfile(outName):
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
            if os.path.isfile(outName):
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
    print '\n**********  Topographic (DEM) Error correction  *******'
    outName = os.path.splitext(inps.timeseries_file)[0]+'_demCor.h5'
    if 'pysar.topoError' in template.keys() and template['pysar.topoError'].lower() in ['yes','y','on']:
        if os.path.isfile(outName):
            print inps.timeseries_file+' already exists.'
            inps.timeseries_file = outName
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
    if 'pysar.deramp' in template.keys():
        deramp_method = template['pysar.deramp'].lower().replace('-','_')
        print 'Phase Ramp Removal method : '+template['pysar.deramp']

        if deramp_method in ['plane', 'quadratic', 'plane_range', 'quadratic_range',\
                             'plane_azimuth', 'quadratic_azimuth']:
            outName = os.path.splitext(inps.timeseries_file)[0]+'_'+deramp_method+'.h5'
            if os.path.isfile(outName):
                print inps.timeseries_file+' already exists.'
                inps.timeseries_file = outName
            else:
                derampCmd = 'remove_plane.py '+inps.timeseries_file+' -m '+inps.mask_file+' -s '+deramp_method
                print derampCmd
                os.system(derampCmd)
            inps.timeseries_file = outName

        elif deramp_method in ['baseline_cor','baselinecor']:
            outName = os.path.splitext(inps.timeseries_file)[0]+'_baselineCor.h5'
            if os.path.isfile(outName):
                print inps.timeseries_file+' already exists.'
                inps.timeseries_file = outName
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
            if os.path.isfile(outName):
                print inps.timeseries_file+' already exists.'
                inps.timeseries_file = outName
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
    velCmd = 'timeseries2velocity.py -f '+inps.timeseries_file+' -t '+inps.template_file
    print velCmd
    os.system(velCmd)
    inps.velocity_file = 'velocity.h5'


    ############################################
    # Geocoding (Optional)
    ############################################
    print '\n**********  Geocoding  ********************************'
    if 'pysar.geocode' in template.keys() and template['pysar.geocode'].lower() in ['y','yes','on']: 
        for File in [inps.velocity_file, inps.temp_coherence_file, inps.timeseries_file]:
            outName = 'geo_'+os.path.basename(File)
            if os.path.isfile(outName):
                print outName+' already existed.'
            else:
                geocodeCmd = 'geocode.py '+inps.geomap_file+' '+File
                print geocodeCmd
                try: os.system(geocodeCmd)
                except: pass
    else:
        print 'No geocoding applied'
    inps.geo_velocity_file = 'geo_'+os.path.basename(inps.velocity_file)
    inps.geo_temp_coherence_file = 'geo_'+os.path.basename(inps.temp_coherence_file)

    #############################################
    # Masking (Optional)
    #############################################
    print '\n**********  Masking Velocity  *************************'
    maskCmd = 'mask.py -f '+inps.velocity_file+' -m '+inps.mask_file
    print maskCmd
    os.system(maskCmd)
    inps.velocity_file = os.path.splitext(inps.velocity_file)[0]+'_masked.h5'
    
    if os.path.isfile(inps.geo_velocity_file):
        maskCmd = 'mask.py -f '+inps.geo_velocity_file+' -m '+inps.geo_temp_coherence_file+' -t 0.7'
        print maskCmd
        os.system(maskCmd)
        inps.geo_velocity_file = os.path.splitext(inps.geo_velocity_file)[0]+'_masked.h5'
    

    #############################################
    # Google Earth KML file (Optional)
    #############################################
    print '\n*********  Creating Google Earth KML file  ************'
    if os.path.isfile(inps.geo_velocity_file):
        kmlCmd = 'save_kml.py '+inps.geo_velocity_file
        print kmlCmd
        os.system(kmlCmd)
    else:
        print 'No geocoded velocity file found, skip creating KML file.'

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

