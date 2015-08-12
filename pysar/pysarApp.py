#! /usr/bin/env python
###############################################################################
# 
# Project: PySAR 
# Purpose: Python Module for InSAR Time-series Analysis
# Author: Heresh Fattahi
# Created: July 2013
# Modified: Yunjun Zhang, Feb 2015
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


import os
import sys
import glob
import time
import _readfile as readfile
import h5py
import subprocess
from pysar._pysar_utilities import check_variable_name


def radar_Or_geo(igramFile):

  h5file=h5py.File(igramFile,'r')
  
  igramList=h5file['interferograms'].keys()

  if 'X_FIRST' in h5file['interferograms'][igramList[0]].attrs.keys():
     rdr_geo='geo'
  else:
     rdr_geo='radar'
  h5file.close()
  return rdr_geo


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

usage:
        
      pysarApp.py TEMPLATEFILE

example:

 pysarApp.py /nethome/hfattahi/SanAndreasT356EnvD.template
 pysarApp.py $TE/SanAndreasT356EnvD.template

*******************************************************
Template file options:

pysar.inputdata=/scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*c10.unw
pysar.CorFiles = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.cor
pysar.wraped = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.int
pysar.geomap = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/GEO/geomap_12/geomap_8rlks.trans
pysar.dem = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG_20050102_20070809/radar_8lks.hgt

pysar.subset.yx = 1800:2000,700:800
pysar.seed.ll=31.5, 67  or  pysar.seed.yx=257 , 151

pysar.unwrap_error = yes [no]

pysar.tropospheric_delay = yes ['no']
pysar.tropospheric_delay.method = pyaps ['height-correlation'] 
pysar.Numerical_Weather_Model = ECMWF ['MERRA', 'NARR']
pysar.acquisition_time = 00:00 ['06:00', '12:00', '18:00']

pysar.topo_error = yes [no]

pysar.orbit_error = yes [np]
pysar.orbit_error.method = plane  ['quadratic', 'plane', 'quardatic_range', 'quadratic_azimiuth', 'plane_range', 'plane_azimuth','baselineCor','BaseTropCor']

pysar.mask=yes
pysar.mask.threshold = 0.7
pysar.geocode = yes
*******************************************************
  '''
#########################################
def main(argv):
  try:
    templateFile = argv[1]
  except:
    Usage();sys.exit(1)

  projectName = os.path.basename(templateFile.partition('.')[0])

  try:
     tssarProjectDir = os.getenv('TSSARDIR') +'/'+projectName
  except:
     tssarProjectDir = os.getenv('SCRATCHDIR') + '/' + projectName + "/TSSAR"     # FA 7/2015: adopted for new directory structure
 
  print "QQ " + tssarProjectDir
  if not os.path.isdir(tssarProjectDir): os.mkdir(tssarProjectDir)
  os.chdir(tssarProjectDir)

  igramFile = 'LoadedData.h5'
  Modified_igramFile = 'Modified_LoadedData.h5'
  if os.path.isfile(Modified_igramFile):
    print Modified_igramFile + ' already exists.'
    igramFile=Modified_igramFile

  template = readfile.read_template(templateFile)
  Rlooks = template['Rlooks_unw']

#########################################
# Loading interferograms
#########################################
  print '******************************************'
  print''
  if os.path.isfile(igramFile):
    print igramFile + ' already exists.'
  else:
    loadCmd='load_data.py ' + templateFile 
    print loadCmd
    os.system(loadCmd)

#    copyDemCmd='copy_dem_trans.py ' + templateFile
#    print copyDemCmd
#    os.system(copyDemCmd)
  print''
  print '******************************************'


#########################################
# Check the subset
#########################################
  try:
    subset= template['pysar.subset.yx'].split(',')
    print subset
    print subset[0]
    subsetOutName='subset_'+igramFile
    subsetCmd='subset.py -f '+ igramFile + ' -y '+subset[0]+' -x '+subset[1] + ' -o ' + subsetOutName
    print '*****************************************'
    print 'Subset the area ...'
    print subsetCmd
    os.system(subsetCmd)
    igramFile=subsetOutName
    print '*****************************************'
  except:
    print '*****************************************'
    print 'No Subset selected. Processing the whole area'
    print '*****************************************'

#########################################
#Referencing all interferograms to the same pixel 
#########################################
  rdr_or_geo=radar_Or_geo(igramFile)  
  print '******************************************'
  print''
  if os.path.isfile('Seeded_'+igramFile):
    igramFile = 'Seeded_'+igramFile
    print igramFile + ' already exists.'
  else:
    print 'referncing all interferograms to the same pixel.'
    if 'pysar.seed.ll' in template.keys():
       'Checking the lat/lon refernce point' 
       lat= template['pysar.seed.ll'].split(',')[0]
       lon= template['pysar.seed.ll'].split(',')[1]
       seedCmd= 'SeedData.py -f ' + igramFile + ' -l ' +lat+ ' -L '+lon
    elif 'pysar.seed.yx' in template.keys():
       'Checking y/x reference point'
       y= template['pysar.seed.yx'].split(',')[0]
       x= template['pysar.seed.yx'].split(',')[1]
       seedCmd= 'seed_data.py -f ' + igramFile + ' -y ' +y+ ' -x '+x
    else: 
       seedCmd= 'seed_data.py -f ' + igramFile

    igramFile = 'Seeded_'+igramFile
    print seedCmd  
    os.system(seedCmd)
  print''
  print '******************************************'

############################################
#unwrapping error correction based on the 
# consistency of triplets of interferograms  
############################################
  print '******************************************'
  print''
  try:
       template['pysar.unwrap_error']
       if template['pysar.unwrap_error'] in ('y','yes','Yes','YES'):
         print 'unwrapping error correction might take a while depending on the size of your data set! '
         unwCmd='unwrap_error.py '+igramFile
         os.system(unwCmd)
         igramFile=igramFile.split('.')[0]+'_unwCor.h5'      
       else:
         print 'No unwrapping error correction.'
  except:
       print 'No unwrapping error correction.'
  print''
  print '******************************************'

#########################################
# inversion of interferograms
########################################
  print '******************************************'
  print''
  if os.path.isfile(igramFile.split('.')[0]+'_unwCor.h5'):
     igramFile = igramFile.split('.')[0]+'_unwCor.h5'
     print igramFile + ' exists.'

  if os.path.isfile('timeseries.h5'):
     print 'timeseries.h5 already exists, inversion is not needed.'
  else:
     invertCmd = 'igram_inversion.py '+ igramFile 
     print invertCmd
     os.system(invertCmd)
  timeseriesFile='timeseries.h5'
  print''
  print '******************************************'

##############################################
#temporal coherence: 
#A parameter to evaluate the consistency of 
# timeseries with the interferograms
##############################################
  print '******************************************'
  print''
#  if os.path.isfile('temporal_coherence.h5'):
#     print 'temporal_coherence.h5 already exists.'
#  else:
#     tempcohCmd='temporal_coherence.py '+igramFile+' '+timeseriesFile
#     print tempcohCmd
#     os.system(tempcohCmd)
  tempcohCmd='temporal_coherence.py '+igramFile+' '+timeseriesFile
  print tempcohCmd
  os.system(tempcohCmd)
  print''
  print '******************************************'


##############################################
#update Mask based on temporal coherence
# add by Yunjun Feb 15, 2015
##############################################
  print '******************************************'
  print''
  try:
       template['pysar.mask']
       if template['pysar.mask'] in ('yes','Yes','YES','y'):
         print 'Updating mask according to temporal coherence'
         cohT=template['pysar.mask.threshold']
         maskCmd='generate_mask.py -f temporal_coherence.h5 -m '+ cohT +' -M 1.0 -o Mask.h5'
         print maskCmd
         os.system(maskCmd)
       else:
         print 'No update for mask.'
  except:
       print 'No update for mask.'
  print''
  print '******************************************'


##############################################
# Generate incident angle
# add by Yunjun Feb 15, 2015 
##############################################
  print '******************************************'
  print''
  inciCmd='incidence_angle.py -f timeseries.h5'
  print inciCmd
  os.system(inciCmd)
  print''
  print '******************************************'

##############################################
#If Satellite is Envisat and if Coordinate
#system is radar then LOD correction
##############################################
  print '******************************************'
  print''
  h5file=h5py.File(timeseriesFile,'r')
  if rdr_or_geo =='radar':
    if h5file['timeseries'].attrs['PLATFORM']=='ENVISAT':
      LODcmd='lod.py '+timeseriesFile
      print LODcmd
      os.system(LODcmd)
      timeseriesFile=timeseriesFile.split('.')[0]+'_LODcor.h5'
  print''
  print '******************************************'


##############################################
# Tropospheric Correction
##############################################
  print '******************************************'
  print''

  try:
     if (template['pysar.tropospheric_delay'] in ('y','yes','Yes','YES')) and template['pysar.orbit_error.method']=='BaseTropCor':
         print '''
   +++++++++++++++++++++++++++++++++++++++++++++++++++
   WARNING:
       Orbital error correction was BaseTropCor.
       Tropospheric correction was already applied simultaneous with baseline error correction.
       Tropospheric correction can not be applied again.
       To apply the tropospheric correction separate from baseline error correction, chhose other existing options for orbital error correction.
    +++++++++++++++++++++++++++++++++++++++++++++++++++      
         '''
         template['pysar.tropospheric_delay']='no'
  except:
     print 'Checking the tropospheric delay correction ...'

  if template['pysar.tropospheric_delay'] in ('y','yes','Yes','YES'):     
    # demFile='radar_'+Rlooks+'rlks.hgt'
     demFile=template['pysar.dem']
     demFile=check_variable_name(demFile)
#     print 'DEM file: '+demFile
     if not os.path.isfile(demFile):
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print 'Error:'
        print 'DEM (radar_*rlks.hgt file) was not found!'
        print 'Continue without tropospheric correction ...'
        print '++++++++++++++++++++++++++++++++++++++++++++++'
     else:       
        if template['pysar.tropospheric_delay.method'] in ['height-correlation','height_correlation','Height-Correlation','Height_Correlation']:
           print 'tropospheric delay correction with height-correlation approach'
           try:
              polyOrder=template['pysar.trop.polyOrder']
           except:
              print 'Deafult polynomial order for troposphreic correction = 1'
              polyOrder='1'
           cmdTrop='tropospheric_correction.py'+ ' -f '+ timeseriesFile + ' -d '+ demfile + ' -p '+ polyOrder
           os.system(cmdTrop)
           timeseriesFile=timeseriesFile.split('.')[0]+'_tropCor.h5'
          
        elif template['pysar.tropospheric_delay.method']=='pyaps':
           print 'Atmospheric correction using Numerical Weather Models (using PyAPS software)'
           print 'reading DEM, source of NWM and acquisition time from template file'
           source_of_NWM=template['pysar.Numerical_Weather_Model']
           print 'Numerical Weather Model: '+source_of_NWM
           acquisition_time=template['pysar.acquisition_time']
           print 'acquisition time: '+acquisition_time
#           cmdTrop = ["tropcor_pyaps.py -f ",timeseriesFile," -d ",demFile," -s ",source_of_NWM," -h ",acquisition_time," -i incidence_angle.h5"]
           cmdTrop = 'tropcor_pyaps.py -f '+timeseriesFile+ ' -d '+ demFile +' -s ' + source_of_NWM + ' -h '+ acquisition_time + ' -i incidence_angle.h5'
           print cmdTrop
           os.system(cmdTrop)
#           subprocess.Popen(cmdTrop).wait()
           timeseriesFile=timeseriesFile.split('.')[0]+'_'+source_of_NWM+'.h5'
        else:
           print 'Atmospheric correction method not recognized.'
  else:
     print 'No atmospheric delay correction.'
  print''
  print '******************************************'

   
##############################################
#topographic residuals
##############################################
  print '******************************************'
  print''
  try:
       template['pysar.topo_error']
       if template['pysar.topo_error'] in ('yes','Yes','YES','y'):
         print 'Correcting topographic residuals'
         topoCmd='dem_error.py '+ timeseriesFile +' '+ igramFile
         print topoCmd
         os.system(topoCmd)
         timeseriesFile=timeseriesFile.split('.')[0]+'_demCor.h5'
       else:
         print 'No correction for topographic residuals.'
  except:
       print 'No correction for topographic residuals.'   
  print''
  print '******************************************'   

##############################################
#Orbit correction  
##############################################
  print '******************************************'
  print''
  try:
       template['pysar.orbit_error']
       if template['pysar.orbit_error'] in ('yes','Yes','YES','y'):
          try:
             orbit_error_method=template['pysar.orbit_error.method']
             print 'orbit error correction method : '+orbit_error_method
             if orbit_error_method in ['quadratic', 'plane', 'quardatic_range', 'quadratic_azimiuth', 'plane_range', 'plane_azimuth']: 
                orbitCmd='remove_plane.py '+timeseriesFile+' '+template['pysar.orbit_error.method'] #+ ' Mask.h5'
                timeseriesFile=timeseriesFile.split('.')[0]+'_'+template['pysar.orbit_error.method']+'.h5'
                print orbitCmd
                os.system(orbitCmd)
             elif orbit_error_method == 'baselineCor':
                orbitCmd='baseline_error.py ' +timeseriesFile #+ ' Mask.h5'
                print orbitCmd
                
                try:
                      h5file=h5py.File(timeseriesFile,'r')
                      daz=float(h5file['timeseries'].attrs['AZIMUTH_PIXEL_SIZE'])
                      os.system(orbitCmd)
                      timeseriesFile=timeseriesFile.split('.')[0]+'_'+template['pysar.orbit_error.method']+'.h5'
                except:
                      print 'WARNING!'
                      print 'Skipping orbital error correction.'
                      print 'baselineCor method can only be applied in radar coordinate'
                        
             elif orbit_error_method =='BaseTropCor':
                demfile=template['pysar.dem']
                demfile=check_variable_name(demfile)
                try:
                    polyOrder=template['pysar.trop.polyOrder']
                except:
                    print 'Deafult polynomial order for troposphreic correction = 1'
                    polyOrder=1  
                try:
                   h5file=h5py.File(timeseriesFile,'r')
                   daz=float(h5file['timeseries'].attrs['AZIMUTH_PIXEL_SIZE']) 
                   orbitCmd='baseline_trop.py '+timeseriesFile+' '+ demfile +' '+ polyOrder +'range_and_azimuth'
                   print 'Joint estimation of Baseline error and tropospheric delay [height-correlation approach]'
                   print orbitCmd
                   os.system(orbitCmd)
                   timeseriesFile=timeseriesFile.split('.')[0]+'_'+template['pysar.orbit_error.method']+'.h5'
                except:
                   print 'WARNING!'
                   print 'Skipping orbital error correction.'
                   print 'baselineCor method can only be applied in radar coordinate'
             else:
                print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                print 'WARNING!'
                print 'Orbital error correction method was not recognized!'
                print 'Possible options are:'
                print 'quadratic, plane, quardatic_range, quadratic_azimiuth, plane_range, plane_azimuth,baselineCor,BaseTropCor'
                print 'Continue without orbital errors correction...'
                print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    
          except:
             print 'No orbital errors correction.'
       else:
          print 'No orbital errors correction.'
  except:
       print 'No orbital errors correction.'
  print''
  print '******************************************'

#############################################
#Velocity and rmse maps
#############################################

  print '******************************************'
  print''
  velCmd='timeseries2velocity.py '+timeseriesFile
  print velCmd
  os.system(velCmd)
  print''
  print '******************************************'

#############################################
#Masking the velocity based on the temporal 
#coherence or rmse if it's specified
#############################################
  print '******************************************'
  print''
  try:
       template['pysar.mask']
       if template['pysar.mask'] in ('yes','Yes','YES','y'):
          try: 
             template['pysar.mask.threshold']
             maskCmd='masking.py -f velocity.h5 -m temporal_coherence.h5 -t '+template['pysar.mask.threshold']
             print 'Masking the velocity file using the temporal coherence with the threshold of '+template['pysar.mask.threshold']
          except:
             maskCmd='Masking.py -f velocity.h5 -m temporal_coherence.h5 -t 0.7'
             print 'Masking the velocity file using the temporal coherence with the threshold of 0.7'
          
          os.system(maskCmd)
#          rmCmd='rm velocity.h5'
#          os.system(rmCmd)
#          mvCmd='mv velocity_masked.h5 velocity.h5'
#          os.system(mvCmd)
       else:
          print 'No masking applied'
  except:
       print 'No masking applied'
  print''
  print '******************************************'

############################################
#Geocoding
############################################
  print '******************************************'
  print''
  try:
       template['pysar.geocode']
       if template['pysar.geocode'] in ('y','yes','Yes','YES'):
          geomapFile='geomap_'+Rlooks+'rlks.trans'
#          geoCmd = 'geocode.py '+timeseriesFile+' '+geomapFile
#          print geoCmd
#          os.system(geoCmd)
          geoCmd = 'geocode.py velocity.h5 '+geomapFile
          print geoCmd
          os.system(geoCmd)
          geoCmd = 'geocode.py Mask.h5 '+geomapFile
          print geoCmd
          os.system(geoCmd)
          
#          maskCmd = 'Masking.py -f geo_'+timeseriesFile+' -m geo_Mask.h5'
#          print maskCmd
#          os.system(maskCmd)
          maskCmd = 'masking.py -f geo_velocity.h5 -m geo_Mask.h5'
          print maskCmd
          os.system(maskCmd)
       else:
          print 'No geocoding applied'
  except: 
       print 'No geocoding applied'
  print''
  print '******************************************'


#############################################
#                   PySAR v1.0              #
#############################################
  print''
  print '###############################################'
  print ''
  print 'End of PySAR processing.'
  print ''
  print '################################################'


if __name__ == '__main__':
  main(sys.argv[:])

