#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Reference:
# Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer
# (2011), Systematic InSAR tropospheric phase delay corrections
# from global meteorological reanalysis data, Geophys. Res. Lett.,
# 38, L17311, doi:10.1029/2011GL048757


import sys
import getopt
import os

try: import pyaps as pa
except: print 'Cannot import pyaps into Python!'; sys.exit(1)
import h5py
from numpy import cos,zeros,pi

import pysar._pysar_utilities as ut


###############################################################
def Usage():
  print '''
  ##############################################################################
  Tropospheric correction using weather models. 
    PyAPS is used to download and calculate the delay for each time-series epoch.
  
  Usage:
      tropcor_pyaps.py -f timeseries.h5 -d demfile.hgt -s source_of_atmospheric_data -h acquisition_time -D Delay_Type -i incidence_angle
        

      -f: timeseries HDF5 file, i.e. timeseries.h5, timeseries_LODcor.h5
      -s: source of the atmospheric data: ECMWF, NARR
      -D: Delay_Type: Dry, Wet , comb [Deafult is comb which calculates both wet and dry delays]   
      -i: incidence_angle  : can be a file containing all incidence angles or can be one average value. 
                                If it's not introduced, average look angle is used.
      -h: time of data (ECMWF takes hh:mm, NARR takes hh only)

  Example:
      
      tropcor_pyaps.py -f timeseries.h5        -d radar_8rlks.hgt -s ECMWF -h 18:00 -i incidence_angle.h5
      tropcor_pyaps.py -f timeseries.h5        -d radar_8rlks.hgt -s NARR  -h 18    -i incidence_angle.h5
      tropcor_pyaps.py -f timeseries.h5        -d radar_8rlks.hgt -s ECMWF -h 18:00 -D Dry -i 23
      tropcor_pyaps.py -f timeseries_LODcor.h5 -d radar_8rlks.hgt -s ECMWF -h 18:00

  ##############################################################################
  '''

###############################################################
def main(argv):

    DelayType='comb'

    try:
      opts, args = getopt.getopt(argv,"f:d:s:h:D:i:")
    except getopt.GetoptError:
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if   opt == '-f':        timeSeriesFile = arg
      elif opt == '-d':        demFile   = arg
      elif opt == '-s':        atmSource = arg.upper()
      elif opt == '-h':        hr        = arg
      elif opt == '-D':        DelayType = arg
      elif opt == '-i':        inc_angle = arg
     
    try:
      timeSeriesFile
      demFile
    except:
      Usage() ; sys.exit(1)

    demFile  = ut.check_variable_name(demFile)
    demCoord = ut.radar_or_geo(demFile)

    h5timeseries = h5py.File(timeSeriesFile)
    yref=h5timeseries['timeseries'].attrs['ref_y']
    xref=h5timeseries['timeseries'].attrs['ref_x']

    ###############################################################
    #incidence angle to map the zenith delay to the slant delay
    try:
       inc_angle
    except:
       print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       print 'WARNING:'
       print 'incedence angle is not specified >>>> Average look angle is used ... '
       print 'For more precise results use input option -i to introduce the incidence angle file or average incodence angle'
       print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       input_inc_angle='None'

       wavelength=float(h5timeseries['timeseries'].attrs['WAVELENGTH'])
       inc_angle1=float(h5timeseries['timeseries'].attrs['LOOK_REF1'])
       inc_angle2=float(h5timeseries['timeseries'].attrs['LOOK_REF2'])
       inc_angle=(inc_angle1+inc_angle2)/2.0
       print '*******************************************************************************'
       print 'Near Look Angle: ' + str(inc_angle1)
       print 'Far  Look Angle:' + str(inc_angle2)
       print 'Average Look Angle (used in pyaps to calculate delay): ' + str(inc_angle)
       print 'Acquisition time is : '+ hr
       print '*******************************************************************************'
       inc_angle=str(inc_angle)

    if os.path.isfile(inc_angle):
       incidenceFile=inc_angle       
       h5incidence=h5py.File(incidenceFile,'r')
       iset=h5incidence['mask'].get('mask')
       inc_angle=iset[0:iset.shape[0],0:iset.shape[1]]
       h5incidence.close()
    else:
       inc_angle=float(inc_angle)
       print 'incidence angle = '+ str(inc_angle)

    inc_angle=inc_angle*pi/180.0
    ################################################################
    dateList = h5timeseries['timeseries'].keys()

    if atmSource == 'ECMWF':
       gribSource='ECMWF'
       if not os.path.isdir('ECMWF'):
          print 'making directory: ECMWF'
          os.mkdir('ECMWF')

       ecmwf_file=[]
       for d in dateList:
          ecm='./ECMWF/ERA-Int_'+d+'_'+hr+'.grb'
          ecmwf_file.append(ecm)
          print [d]
          if not os.path.isfile(ecm):
             pa.ECMWFdload([d],hr,'./ECMWF/')
          else:
             print ecm + ' already exists.'

    elif atmSource == 'NARR':
       gribSource='NARR'
       if not os.path.isdir('NARR'):
          print 'making directory: NARR'
          os.mkdir('NARR')

       ecmwf_file=[]
       for d in dateList:
          ecm='./NARR/narr-a_221_'+d+'_'+hr+'00_000.grb'
          ecmwf_file.append(ecm)
          print [d]
          if not os.path.isfile(ecm):
             pa.NARRdload([d],hr,'./NARR/')
          else:
             print ecm + ' already exists.'

    elif atmSource == 'ERA':
       gribSource='ERA'
       if not os.path.isdir('ERA'):
          print 'making directory: ERA'
          os.mkdir('ERA')

       ecmwf_file=[]
       for d in dateList:
          ecm='./ERA/ERA_'+d+'_'+hr+'.grb'
          ecmwf_file.append(ecm)
          print [d]
          if not os.path.isfile(ecm):
             pa.ERAdload([d],hr,'./ERA/')
          else:
             print ecm + ' already exists.'

    elif atmSource == 'MERRA':
       gribSource='MERRA'
       if not os.path.isdir('MERRA'):
          print 'making directory: MERRA'
          os.mkdir('MERRA')

       ecmwf_file=[]
       for d in dateList:
          ecm='./MERRA/merra-'+d+'-'+hr+'.hdf'
          ecmwf_file.append(ecm)
          print [d]
          if not os.path.isfile(ecm):
             pa.MERRAdload([d],hr,'./MERRA/')
          else:
             print ecm + ' already exists.'

    else:
       Usage();sys.exit(1)

    print '*******************************************************************************'
    print 'Calcualting delay for each epoch.'
    h5phsName=atmSource + '.h5'
    h5phs=h5py.File(h5phsName,'w')
    outName=timeSeriesFile.replace('.h5','_') + atmSource + '.h5'	#Yunjun, Feb 15, 2015
    #outName=timeSeriesFile.replace('.h5','')+'_tropCorPyAPS.h5' 
    h5apsCor=h5py.File(outName,'w')
    group=h5apsCor.create_group('timeseries')
    group_phs=h5phs.create_group('timeseries')

    #if 'X_FIRST' in  h5timeseries['timeseries'].attrs.keys():
    #   demCoord='geo'
    #   print 'The coordinate system is : Geo'
    #else:
    #   demCoord='radar'  
    #   print 'The coordinate system is : radar'      

    print ecmwf_file[0]
    if demCoord=='radar':  aps1 = pa.PyAPS_rdr(str(ecmwf_file[0]),demFile,grib=gribSource,verb=True,Del=DelayType)
    else:                  aps1 = pa.PyAPS_geo(str(ecmwf_file[0]),demFile,grib=gribSource,verb=True,Del=DelayType)

    phs1 = zeros((aps1.ny,aps1.nx))
    aps1.getdelay(phs1)     
    phs1=(phs1 - phs1[yref,xref])/cos(inc_angle)   
    dset = group_phs.create_dataset(dateList[0], data= phs1- phs1, compression='gzip')
    dset = group.create_dataset(dateList[0], data= phs1- phs1, compression='gzip')
    
    for i in range(1,len(ecmwf_file)):
       ecm=ecmwf_file[i] 
       print ecm
       if demCoord=='radar':  aps = pa.PyAPS_rdr(str(ecm),demFile,grib=gribSource,verb=True,Del=DelayType)   
       else:                  aps = pa.PyAPS_geo(str(ecm),demFile,grib=gribSource,verb=True,Del=DelayType)
       phs = zeros((aps.ny,aps.nx))
       aps.getdelay(phs)
       phs=(phs - phs[yref,xref])/cos(inc_angle)
       phs=phs-phs1
       dset  = group_phs.create_dataset(dateList[i], data= phs, compression='gzip')
       dset1 = h5timeseries['timeseries'].get(dateList[i])
       data1 = dset1[0:dset1.shape[0],0:dset1.shape[1]]
       dset  = group.create_dataset(dateList[i], data= data1+phs, compression='gzip')

    for key,value in h5timeseries['timeseries'].attrs.iteritems():
      group.attrs[key] = value
      group_phs.attrs[key] = value

    try:
        dset1 = h5timeseries['mask'].get('mask')[:]
        group = h5apsCor.create_group('mask')
        dset  = group.create_dataset('mask', data=Mask, compression='gzip')
    except: pass


###############################################################
if __name__ == '__main__':

    main(sys.argv[1:])


