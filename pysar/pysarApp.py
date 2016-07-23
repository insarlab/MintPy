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

import os
import sys
import glob
import time

import h5py

import pysar._pysar_utilities as ut
import pysar._readfile as readfile


##################################### Sub Functions #######################################
#############  find accurate file name for template input  #############
def find_filename(template, option, workDir='.'):
  try:
     filename=template[option]
     filename=os.path.basename(filename).split('/')[0]
     filename=glob.glob(workDir+'/'+filename)[0]
     filename=os.path.basename(filename).split('/')[0]
  except:
     if   option == 'pysar.dem.geoCoord':   filename = '*.dem'
     elif option == 'pysar.dem.radarCoord': filename = 'radar*.hgt'
     elif option == 'pysar.geomap':         filename = 'geomap*.trans'
     elif option == 'pysar.mask.file':         
         filename = 'Mask.h5'
         if os.path.isfile('Modified_'+filename):
             filename  = 'Modified_'+filename
     else: print 'Error: Unrecognized option: '+option; sys.exit(1)
     filename=os.path.basename(filename).split('/')[0]
     filename=glob.glob(workDir+'/'+filename)[0]
     filename=os.path.basename(filename).split('/')[0]
  return filename

#############  update the subset of input file  ########################
def check_subset(inName, subset, option='yx', workDir='.'):
  outName='subset_'+inName
  if   os.path.isfile(workDir+'/'+outName):
    print outName + ' already exists.'
  elif os.path.isfile(workDir+'/'+inName):
    if   option == 'yx':
       subsetCmd='subset.py -f '+inName+' -y '+subset[0]+' -x '+subset[1]+' -o '+outName
    elif option == 'lalo':
       subsetCmd='subset.py -f '+inName+' -l '+subset[0]+' -L '+subset[1]+' -o '+outName
    else: print 'Unrecognized option: '+option; return
    print subsetCmd
    os.system(subsetCmd)
  else:
    print 'WARNING: No '+inName+' found.'
    outName = inName
  return outName


#############  update geocoding of input file  #########################
def check_geocode(inName, geomapFile, workDir='.'):
  outName = 'geo_'+inName
  if   os.path.isfile(workDir+'/'+outName):
    print outName+' already exists.'
  elif os.path.isfile(workDir+'/'+inName):
    geoCmd = 'geocode.py '+inName+' '+geomapFile
    print geoCmd
    os.system(geoCmd)
  else:
    outName = inName
    print 'WARNING: No '+inName+' found.'
  return outName


#############  update masking of input file  ###########################
def check_mask(inName, maskFile, workDir='.'):
  outName = inName.split('.')[0]+'_masked.h5'
  if   os.path.isfile(workDir+'/'+outName):
    print outName+' already exists.'
  elif os.path.isfile(workDir+'/'+inName):
    maskCmd = 'masking.py -f '+inName+' -m '+maskFile
    print maskCmd
    os.system(maskCmd)
  else:
    outName = inName
    print 'WARNING: No '+inName+' found.'
  return outName


###########################  Usage Function  ###########################
def Usage():

  print '''
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
PySAR v1.0 July 2013, InSAR Lab, RSMAS, University of Miami

Usage:       
      pysarApp.py TEMPLATEFILE [workDir]

Example:
      pysarApp.py /nethome/hfattahi/SanAndreasT356EnvD.template
      pysarApp.py $TE/SanAndreasT356EnvD.template
      pysarApp.py $TE/SanAndreasT356EnvD.template $SC/TSSAR/SanAndreasT356EnvD

Re-run tips:
      Delete timeseries.h5 to update inversion of time series
      Delete Seeded_*.h5   to update reference point
      Delete all subset_*  to update subset range


*******************************************************
Template file options:
####################################################################

pysar.inputFiles         = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*.unw
pysar.corFiles           = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*rlks.cor         #optional
pysar.wrappedFiles       = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*rlks.int         #optional

pysar.geomap             = $PROCESSDIR/SanAndreasT356EnvD/GEO/geomap_050102-070809/geomap_8rlks.trans
pysar.dem.radarCoord     = $PROCESSDIR/SanAndreasT356EnvD/DONE/IFG_20050102_20070809/radar_8lks.hgt
pysar.dem.geoCoord       = $DEMDIR/SanAndreasT356EnvD/srtm30m.dem                                       #optional
pysar.dem.source         = GSI_10m_DEHM             # DEM source ['SRTM1','SRTM3','ASTER']

pysar.drop.ifgIndex       = 7:9,15,25,26             #optional
pysar.drop.date          =                          #drop date list for time series display and velocity estimation, not in pysarApp.py
pysar.drop.pair          =                          #[not implemented yet]

pysar.subset             = no                       #subset before processing, no by default
pysar.subset.yx          = 1800:2000,700:800        #optional
pysar.subset.lalo        = 31.5:32.5,130.5:131.0    #optional

pysar.seed.yx            = 257 , 151                #optional
pysar.seed.lalo          = 31.8, 130.8              #[not implemented yet]

pysar.unwrapError                     = yes               #['no' (default)], optional
pysar.unwrapError.method              = bonding_point     #['triangular_phase_closure' (default),'bonding_point']
pysar.unwrapError.yx                  = 1160,300,1240,300 #[x_ref1,y_ref1,x1,y1,x_ref2,y_ref2,x2,y2 ... ]

pysar.troposphericDelay               = yes               #['no'], optional
pysar.troposphericDelay.method        = pyaps             #['height-correlation'] 
pysar.troposphericDelay.polyOrder     = 1                 #for 'height-correlation' method
pysar.troposphericDelay.weatherModel  = ECMWF             #['MERRA', 'NARR'], for 'pyaps' mehod
pysar.acquisitionTime                 = 00:00             #['06:00', '12:00', '18:00'], for 'pyaps' method

pysar.topoError          = yes               #['no'], optional

pysar.orbitError         = yes               #['no'], optional
pysar.orbitError.method  = plane             #[    'plane',     'plane_range',     'plane_azimuth']
                                             #['quadratic', 'quadratic_range', 'quadratic_azimuth']
                                             #['baselineCor','BaseTropCor']

pysar.mask               = yes               #['no']
pysar.mask.file          = Modified_Mask.h5  #specify mask rather than the default one (Mask.h5 or Modified_Mask.h5)
pysar.mask.threshold     = 0.7               #threshold for generating mask from temporal coherence [0-1] ; commented - don't update mask
pysar.geocode            = yes               #['no'], optional
pysar.geocode.list       = velocity.h5, timeseries.h5, temporal_coherence.h5  #not implemented yet
pysar.kml                = yes               #['no'], optional

pysar.view.row           = 5
pysar.view.column        = 20
pysar.view.min           = -7
pysar.view.max           = 7
pysar.view.list          = geo_velocity_masked.h5,temporal_coherence.h5       #not implemented yet

####################################################################
  '''


####################################  Main Function  ######################################
def main(argv):
  start = time.time()

  try:     templateFile = argv[1]
  except:  Usage(); sys.exit(1)

  ###########  Path  ############
  projectName = os.path.basename(templateFile).partition('.')[0]
  try:     tssarProjectDir = argv[2]
  except:
     if os.getenv('PARENTDIR'):  tssarProjectDir = os.getenv('SCRATCHDIR')+'/'+projectName+"/TSSAR"
     else:                       tssarProjectDir = os.getenv('TSSARDIR')+'/'+projectName
  print '\n********************************************************'
  print   '*********************    PySAR    **********************'
  print '********************************************************\n'
  print "TSSAR directory: " + os.path.abspath(tssarProjectDir)
  if not os.path.isdir(tssarProjectDir): os.mkdir(tssarProjectDir)
  os.chdir(tssarProjectDir)

  ##########  Initial File Name  ########
  template = readfile.read_template(templateFile)
  try:    template['pysar.subset']
  except: template['pysar.subset'] = 'no'

  igramFile = 'LoadedData.h5'
  if os.path.isfile('Modified_'+igramFile):  igramFile = 'Modified_'+igramFile
  corFile   = 'Coherence.h5'
  if os.path.isfile('Modified_'+corFile):    corFile   = 'Modified_'+corFile

#########################################
# Loading Data
#########################################
  print '\n**********  Loading Data  *****************************'
  if len(glob.glob(igramFile)) > 0:
    print igramFile + ' already exists.'
  else:
    loadCmd='load_data.py ' + templateFile 
    print loadCmd
    os.system(loadCmd)
    #copyDemCmd='copy_dem_trans.py ' + templateFile
    #print copyDemCmd
    #os.system(copyDemCmd)

  #if not os.path.isfile(igramFile): sys.exit('\nERROR: No interferograms file found!\n')

  ##########  Initial File Name - 2  ####
  try:  demGeoFile = find_filename(template, 'pysar.dem.geoCoord');   print 'DEM in geo   coordinate: '+demGeoFile
  except: print '\nWARNING:\n\tNo geo coded DEM found!\n\tMight be a problem in tropospheric delay / orbital error correction!\n'
  try:  demRdrFile = find_filename(template, 'pysar.dem.radarCoord'); print 'DEM in radar coordinate: '+demRdrFile
  except: print '\nWARNING:\n\tNo radar coded DEM found!\n\tWill be a problem in tropospheric delay / orbital error correction!\n'
  try:  geomapFile = find_filename(template, 'pysar.geomap');         print 'geomap file: '+geomapFile
  except: print '\nWARNING:\n\tNo geomap*.trans file found!\n\tWill be a problem in geocoding!\n'
  try:    maskFile = find_filename(template, 'pysar.mask.file');      print 'mask   file: '+maskFile
  except: print '\nWARNING:\n\tNo mask file found!\n\tMight be a problem in the future!\n'

#########################################
# Check the subset (Optional)
#########################################
  print '\n**********  Subseting  ********************************'

  if template['pysar.subset'] == 'yes' and  'pysar.subset.yx' in template.keys():
    print 'subseting data with y/x input...'
    subset = template['pysar.subset.yx'].split(',')

    print 'calculating bounding box in geo coordinate.'
    import numpy as np
    ysub = [float(i) for i in subset[0].split(':')];  ysub.sort()
    xsub = [float(i) for i in subset[1].split(':')];  xsub.sort()
    x = np.array([xsub[0],xsub[1],xsub[0],xsub[1]])
    y = np.array([ysub[0],ysub[0],ysub[1],ysub[1]])
    lat,lon,lat_res,lon_res = ut.radar2glob(x,y,igramFile,1)
    buf = 10*(np.max([lat_res,lon_res]))
    latsub = str(np.min(lat)-buf)+':'+str(np.max(lat)+buf)
    lonsub = str(np.min(lon)-buf)+':'+str(np.max(lon)+buf)
    print '    subset in y - '+subset[0]+'\n    subset in x - '+subset[1]
    print '    subset in lat - '+latsub +'\n    subset in lon - '+lonsub

    ## subset radar coded file
    igramFile  = check_subset(igramFile, subset, option='yx')
    maskFile   = check_subset(maskFile,  subset, option='yx')
    try:    demRdrFile = check_subset(demRdrFile,subset, option='yx')
    except: pass
    ## subset geo coded file
    try:    demGeoFile = check_subset(demGeoFile,[latsub,lonsub], option='lalo')
    except: pass
    try:    geomapFile = check_subset(geomapFile,[latsub,lonsub], option='lalo')
    except: pass

  elif template['pysar.subset'] == 'yes' and 'pysar.subset.lalo' in template.keys():
    print 'subseting data with lat/lon input...'
    subset= template['pysar.subset.lalo'].split(',')

    print 'calculate bounding box in radar coordinate.'
    import numpy as np
    latsub = [float(i) for i in subset[0].split(':')];  latsub.sort()
    lonsub = [float(i) for i in subset[1].split(':')];  lonsub.sort()
    lat = np.array([latsub[0],latsub[0],latsub[1],latsub[1]])
    lon = np.array([lonsub[0],lonsub[1],lonsub[0],lonsub[1]])
    x,y,x_res,y_res = ut.glob2radar(lat,lon)
    buf = 10*(np.max([x_res,y_res]))
    xsub = str(np.min(x)-buf)+':'+str(np.max(x)+buf)
    ysub = str(np.min(y)-buf)+':'+str(np.max(y)+buf)
    print '    subset in lat - '+subset[0]+'\n    subset in lon - '+subset[1]
    print '    subset in y - '  +ysub     +'\n    subset in x - '+xsub

    ## subset geo coded files
    try:    demGeoFile = check_subset(demGeoFile, subset, option='lalo')
    except: pass
    try:    geomapFile = check_subset(geomapFile, subset, option='lalo')
    except: pass
    ## subset radar coded files
    igramFile  = check_subset(igramFile, [ysub,xsub], option='yx')
    maskFile   = check_subset(maskFile,  [ysub,xsub], option='yx')
    try:    demRdrFile = check_subset(demRdrFile,[ysub,xsub], option='yx')
    except: pass

  else: print 'No Subset selected. Processing the whole area'


  #try:
  #  subset= template['pysar.subset.yx'].split(',')
  #  print 'subset in y - '+str(subset[0])
  #  print 'subset in x - '+str(subset[1])
  #  igramFile  = check_subset(igramFile, subset)
  #  corFile    = check_subset(corFile,   subset)
  #  maskFile   = check_subset(maskFile,  subset)
  #  demRdrFile = check_subset(demRdrFile,subset)
  #  demGeoFile = check_subset(demGeoFile,subset)
  #  #geomapFile = check_subset(geomapFile,subset)
  #except:  print   'No Subset selected. Processing the whole area'


#########################################
# Referencing Interferograms
#########################################
  print '\n**********  Referencing Interferograms  ***************'
  if os.path.isfile('Seeded_'+igramFile):
      igramFile = 'Seeded_'+igramFile
      print igramFile + ' already exists.'
  elif os.path.isfile('Modified_Seeded_'+igramFile):
      igramFile = 'Modified_Seeded_'+igramFile
      print igramFile + ' already exists.'
  else:
      print 'referncing all interferograms to the same pixel.'
      seedCmd = 'seed_data.py -f '+igramFile+' -t '+templateFile+' -M '+maskFile
      igramFile = 'Seeded_'+igramFile
      print seedCmd  
      os.system(seedCmd)

  #if os.path.isfile('Modified_'+igramFile):  igramFile = 'Modified_'+igramFile

############################################
# Unwrapping Error Correction (Optional)
#    based on the consistency of triplets
#    of interferograms  
############################################
  print '\n**********  Unwrapping Error Correction  **************'
  outName = igramFile.split('.')[0]+'_unwCor.h5'
  try:
      template['pysar.unwrapError']
      if template['pysar.unwrapError'] in ('Y','y','yes','Yes','YES'):
          if os.path.isfile(outName):
              igramFile = outName
              print igramFile+' exists.'
          else:
              print 'This might take a while depending on the size of your data set!'
              unwCmd='unwrap_error.py -f '+igramFile+' -m '+maskFile
              os.system(unwCmd)
              igramFile = outName
      else:  print 'No unwrapping error correction.'
  except:  print 'No unwrapping error correction.'


#########################################
# Inversion of Interferograms
########################################
  print '\n**********  Time Series Inversion  ********************'
  if os.path.isfile('timeseries.h5'):
     print 'timeseries.h5 already exists, inversion is not needed.'
  else:
     invertCmd = 'igram_inversion.py '+igramFile
     print invertCmd
     os.system(invertCmd)
  timeseriesFile = 'timeseries.h5'


##############################################
# Temporal Coherence: 
#   A parameter to evaluate the consistency 
#   of timeseries with the interferograms
##############################################
  print '\n**********  Generate Temporal Coherence file  *********'
  if os.path.isfile('temporal_coherence.h5'):
     print 'temporal_coherence.h5 already exists.'
  else:
     tempcohCmd='temporal_coherence.py '+igramFile+' '+timeseriesFile
     print tempcohCmd
     os.system(tempcohCmd)
  tempCohFile = 'temporal_coherence.h5'


##############################################
# Update Mask based on temporal coherence (Optional)
##############################################
  print '\n**********  Update Mask based on Temporal Coherence  **'
  outName = 'Mask_tempCoh.h5'
  try:
     template['pysar.mask.threshold']
     try:    cohT = template['pysar.mask.threshold']
     except: cohT = '0.7'
     maskCmd='generate_mask.py -f '+tempCohFile+' -m '+ cohT +' -M 1.0 -o '+outName
     print maskCmd
     os.system(maskCmd)
     maskFile = outName
  except:   print 'No mask update from temporal coherence'


##############################################
# Incident Angle
##############################################
  print '\n**********  Generate Incident Angle file  *************'
  if os.path.isfile('incidence_angle.h5'):
     print 'incidence_angle.h5 already exists.'
  else:
     inciCmd = 'incidence_angle.py -f '+timeseriesFile
     print inciCmd
     os.system(inciCmd)
  incAngleFile = 'incidence_angle.h5'

  ### geoCoord
  rdr_or_geo = ut.radar_or_geo(timeseriesFile)


##############################################
# LOD (Local Oscillator Drift) Correction
#   when Satellite is Envisat and
#   Coordinate system is radar
##############################################
  print '\n**********  Local Oscillator Drift correction  ********'
  outName = timeseriesFile.split('.')[0]+'_LODcor.h5'
  if os.path.isfile(outName):
     timeseriesFile = outName;     print timeseriesFile+' already exists.'
  else:
     h5file   = h5py.File(timeseriesFile,'r')
     platform = h5file['timeseries'].attrs['PLATFORM']
     if platform == 'ENVISAT':
        if rdr_or_geo == 'radar':
           LODcmd='lod.py '+timeseriesFile
           print LODcmd
           os.system(LODcmd)
           timeseriesFile = outName
        else: print 'Cannot correct LOD for geocoded data.'
     else: print 'No need of LOD correction for '+platform
     h5file.close()


##############################################
# Tropospheric Delay Correction (Optional)
##############################################
  print '\n**********  Tropospheric Correction  ******************'
  try:
     if (template['pysar.troposphericDelay'] in ('Y','y','yes','Yes','YES')) and\
         template['pysar.orbitError.method'] in ['BaseTropCor','basetropcor','base_trop_cor']:
        print '''
   +++++++++++++++++++++++++++++++++++++++++++++++++++
   WARNING:
       Orbital error correction was BaseTropCor.
       Tropospheric correction was already applied simultaneous with baseline error correction.
       Tropospheric correction can not be applied again.
       To apply the tropospheric correction separated from baseline error correction, \
          choose other existing options for orbital error correction.
    +++++++++++++++++++++++++++++++++++++++++++++++++++      
        '''
        template['pysar.troposphericDelay']='no'
  except:  print 'Checking the tropospheric delay correction ...'

  if template['pysar.troposphericDelay'] in ('Y','y','yes','Yes','YES'):     
     if   rdr_or_geo == 'radar':  demFile = demRdrFile
     elif rdr_or_geo == 'geo':    demFile = demGeoFile

     if not os.path.isfile(demFile):
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print 'ERROR:'
        print '    DEM file was not found!'
        print '    Continue without tropospheric correction ...'
        print '++++++++++++++++++++++++++++++++++++++++++++++'
     else:
        if template['pysar.troposphericDelay.method'] in ['height-correlation','height_correlation',\
                                                          'Height-Correlation','Height_Correlation']:
           print 'tropospheric delay correction with height-correlation approach'
           outName = timeseriesFile.split('.')[0]+'_tropCor.h5'
           if os.path.isfile(outName):
              timeseriesFile = outName
              print timeseriesFile+' already exists.'
           else:
              try:     poly_order = template['pysar.troposphericDelay.polyOrder']
              except:  poly_order = '1';  print 'Deafult polynomial order for troposphreic correction = 1'
              cmdTrop = 'tropcor_phase_elevation.py'+' -f '+timeseriesFile+' -d '+demFile+' -p '+str(poly_order)+' -m '+maskFile
              print cmdTrop
              os.system(cmdTrop)
              timeseriesFile = outName

        elif template['pysar.troposphericDelay.method'] in ['pyaps','PyAPS','PYAPS']:
           print 'Atmospheric correction using Numerical Weather Models (using PyAPS software)'
           print 'reading DEM, source of NWM and acquisition time from template file'
           source_of_NWM = template['pysar.troposphericDelay.weatherModel']
           print 'Numerical Weather Model: '+source_of_NWM
           outName = timeseriesFile.split('.')[0]+'_'+source_of_NWM+'.h5'
           if os.path.isfile(outName):
              timeseriesFile = outName
              print timeseriesFile+' already exists.'
           else:
              acquisition_time = template['pysar.acquisitionTime']
              print 'acquisition time: '+acquisition_time
              cmdTrop = 'tropcor_pyaps.py -f '+timeseriesFile+' -d '+demFile+' -s '+source_of_NWM+\
                                        ' -h '+acquisition_time+' -i '+incAngleFile
              print cmdTrop
              os.system(cmdTrop)
              #subprocess.Popen(cmdTrop).wait()
              timeseriesFile = outName
        else:  print 'ERROR: Unrecognized atmospheric correction method: '+template['pysar.trop.method']
  else:  print 'No atmospheric delay correction.'


##############################################
# Topographic (DEM) Residuals Correction (Optional)
##############################################
  print '\n**********  Topographic (DEM) Error correction  *******'
  outName = timeseriesFile.split('.')[0]+'_demCor.h5'
  try:
     template['pysar.topoError']
     if template['pysar.topoError'] in ('Y','yes','Yes','YES','y'):
        if os.path.isfile(outName):
           timeseriesFile = outName
           print timeseriesFile+' already exists.'
        else:
           print 'Correcting topographic residuals'
           topoCmd='dem_error.py '+ timeseriesFile
           print topoCmd
           os.system(topoCmd)
           timeseriesFile = outName
     else:  print 'No correction for topographic residuals.'
  except:   print 'No correction for topographic residuals.'


##############################################
# Orbit Correction (Optional)
##############################################
  print '\n**********  Orbital Correction  ***********************'
  try:
     template['pysar.orbitError']
     if template['pysar.orbitError'] in ('Y','yes','Yes','YES','y'):
        try:
           orbit_error_method=template['pysar.orbitError.method'].lower()
           print 'orbit error correction method : '+orbit_error_method
           outName = timeseriesFile.split('.')[0]+'_'+orbit_error_method+'.h5'
           if os.path.isfile(outName):
              timeseriesFile = outName
              print timeseriesFile+' already exists.'

           else:
              if orbit_error_method in [    'plane',     'plane_range',     'plane_azimuth',\
                                        'quadratic', 'quadratic_range', 'quadratic_azimuth']:
                 orbitCmd='remove_plane.py -f '+timeseriesFile+' -t '+templateFile
                 print orbitCmd
                 os.system(orbitCmd)
                 timeseriesFile = outName

              elif orbit_error_method in ['baselineCor','BaselineCor']:
                 try:
                    h5file=h5py.File(timeseriesFile,'r');   daz=float(h5file['timeseries'].attrs['AZIMUTH_PIXEL_SIZE']); h5file.close();
                    orbitCmd='baseline_error.py ' +timeseriesFile #+ ' Mask.h5'
                    print orbitCmd
                    os.system(orbitCmd)
                    timeseriesFile = outName
                 except:
                    print 'WARNING!'
                    print 'Skipping orbital error correction.'
                    print 'baselineCor method can only be applied in radar coordinate'

              elif orbit_error_method in ['BaseTropCor','basetropcor','base_trop_cor']:
                 try:
                    h5file=h5py.File(timeseriesFile,'r');   daz=float(h5file['timeseries'].attrs['AZIMUTH_PIXEL_SIZE']); h5file.close();
                    print 'Joint estimation of Baseline error and tropospheric delay [height-correlation approach]'
                    if   rdr_or_geo == 'radar':  demFile = demRdrFile
                    elif rdr_or_geo == 'geo':    demFile = demGeoFile
                    try:     poly_order = template['pysar.trop.poly_order']
                    except:  poly_order = 1;  print 'Deafult polynomial order for troposphreic correction = 1'
                    orbitCmd='baseline_trop.py '+timeseriesFile+' '+ demFile +' '+ str(poly_order) +' range_and_azimuth'
                    print orbitCmd
                    os.system(orbitCmd)
                    timeseriesFile = outName
                 except:
                    print 'WARNING!'
                    print 'Skipping orbital error correction.'
                    print 'baselineCor method can only be applied in radar coordinate'

              else:
                 print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                 print 'WARNING!'
                 print 'Orbital error correction method was not recognized!'
                 print 'Possible options are:'
                 print '    quadratic, plane, quardatic_range, quadratic_azimiuth'
                 print '    plane_range, plane_azimuth,baselineCor,BaseTropCor'
                 print 'Continue without orbital errors correction...'
                 print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        except:  print 'No orbital errors correction.'
     else:       print 'No orbital errors correction.'
  except:        print 'No orbital errors correction.'


#############################################
# Velocity and rmse maps
#############################################
  print '\n**********  Velocity estimation  **********************'
  if os.path.isfile('velocity.h5'):
    print 'velocity.h5 file already exists.'
  else:
    velCmd='timeseries2velocity.py -f '+timeseriesFile+' -t '+templateFile
    print velCmd
    os.system(velCmd)
  velocityFile = 'velocity.h5'


############################################
# Geocoding (Optional)
############################################
  print '\n**********  Geocoding  ********************************'
  try:
     template['pysar.geocode']
     if template['pysar.geocode'] in ('Y','y','yes','Yes','YES'):
        geoTsFile       = check_geocode(timeseriesFile,geomapFile)
        geoVelocityFile = check_geocode(velocityFile,  geomapFile)
        geoMaskFile     = check_geocode(maskFile,      geomapFile)
     else:  print 'No geocoding applied'
  except:   print 'No geocoding applied'


#############################################
# Masking (Optional)
#############################################
  print '\n**********  Masking Velocity  *************************'
  try:
     template['pysar.mask']
     if template['pysar.mask'] in ['Y','y','yes','Yes','YES']:
        velocityFile  = check_mask(velocityFile, maskFile)
        try:    geoVelocityMaskedFile = check_mask(geoVelocityFile, geoMaskFile)
        except: pass
     else:  print 'No masking applied.'
  except:   print 'No masking applied.'

#############################################
# Google Earth KML file (Optional)
#############################################
  print '\n*********  Creating Google Earth KML file  ************'
  try:
     template['pysar.kml']
     if template['pysar.kml'] in ['Y','y','yes','Yes','YES']:
        if os.path.isfile(geoVelocityFile):
           kmlCmd = 'save_kml.py '+geoVelocityFile; print kmlCmd
           os.system(kmlCmd)
        if os.path.isfile(geoVelocityMaskedFile):
           kmlCmd = 'save_kml.py '+geoVelocityMaskedFile; print kmlCmd
           os.system(kmlCmd)
     else:  print 'No KML file creation applied.'
  except:   print 'No KML file creation applied.'


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

