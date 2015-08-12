#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import _gmt
import h5py
from numpy import linspace,meshgrid,flipud
import sys
import os


def get_geo_lat_lon(h5Geo):

   k=h5Geo.keys()
   if 'timeseries' in k:
      k[0]='timeseries'  
 
   X_FIRST=float(h5Geo[k[0]].attrs['X_FIRST'])
   Y_FIRST=float(h5Geo[k[0]].attrs['Y_FIRST'])
   X_STEP=float(h5Geo[k[0]].attrs['X_STEP'])
   Y_STEP=float(h5Geo[k[0]].attrs['Y_STEP'])
   W=float(h5Geo[k[0]].attrs['WIDTH'])
   L=float(h5Geo[k[0]].attrs['FILE_LENGTH'])
   Y_END=Y_FIRST+(L-1)*Y_STEP
   X_END=X_FIRST+(W-1)*X_STEP
   X=linspace(X_FIRST,X_END,W)
   Y=linspace(Y_FIRST,Y_END,L)
   XI,YI=meshgrid(X,Y)

   return Y,X


#file='geo_velocity_LOD_DEM_MERISwet_masked_masked.h5'
try:

  file=sys.argv[1]

except:
  print '''
  ******************************************************
  Exporting geocoded pysar velocity file to GMT grd file
  It also exports an epoch of timeseries to the grd file 

  Example:

      save_gmt.py  geocoded_velocity.h5

      save_gmt.py  geocoded_timeseries.h5 20071031

      save_gmt.py  geocoded_timeseries.h5
  *****************************************************
'''
  sys.exit(1)


print 'Input file: '+file

h5=h5py.File(file,'r')
k=h5.keys()

lats,lons=get_geo_lat_lon(h5)

if 'timeseries' in k:
   try:
     d=sys.argv[2]
   except:
     print 'No input date... continue to convert the last date of timeseries'
     dateList=h5['timeseries'].keys()
     d=dateList[-1]
   print 'reading '+d + ' ... '
   dset=h5['timeseries'].get(d)

elif 'interferograms' in k:

   print 'interferograms' 

else:
   dset=h5[k[0]].get(k[0])


z=dset[0:dset.shape[0],0:dset.shape[1]]


#fname='velocity.grd'
outName=os.path.basename(file).replace('.h5','')
fname = outName + '.grd'
print 'Output (GMT grd file): '+ fname
_gmt.write_gmt_simple(lons, flipud(lats), flipud(z), fname, title='default', name='velocity', scale=1.0, offset=0, units='meters')

#_gmt.write_gmt_simple(lons, lats, res, fname, title=title, name=name, units=units)
