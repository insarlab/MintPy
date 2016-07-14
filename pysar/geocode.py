#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jun 2015: Finish 'interferograms' option
# Yunjun, Jul 2015: Add 'coherence'/'wrapped' option
# Emre,   Sep 2015: Add 'dem' option for DEM_error.h5
# Yunjun, Oct 2015: Add geocode_one()
#                   Merge 'interferograms','coherence','wrapped' into one
#                   Add support for subsetted radar coded files


import os
import sys

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile


######################  Geocode one data  ########################
def geocode_one(data,geomapFile,outname):

   print 'writing to roi_pac unw file format'
   writefile.write_float32(data,outname)
   f = open(outname+'.rsc','w')
   f.write('FILE_LENGTH       '+str(data.shape[0])+'\n')
   f.write('WIDTH             '+str(data.shape[1])+'\n')
   f.close()

   geoCmd='geocode.pl '+geomapFile+' '+outname+' geo_'+outname
   print geoCmd
   os.system(geoCmd)

   print 'reading geocoded file...'
   amp,unw,unwrsc = readfile.read_float32('geo_'+outname)

   rmCmd='rm '+outname;                 os.system(rmCmd);       print rmCmd
   rmCmd='rm '+outname+'.rsc';          os.system(rmCmd);       print rmCmd
   rmCmd='rm geo_'+outname;             os.system(rmCmd);       print rmCmd
   rmCmd='rm geo_'+outname+'.rsc';      os.system(rmCmd);       print rmCmd

   return amp, unw, unwrsc


def Usage():
   print '''
*******************************************************************
*******************************************************************
    
    uses roi_pac geocoding functions to geocode the PySAR products

    Usage:
         geocode.py File geocoding_lookupfile

    File: PySAR hdf5 file [velocity, temporal_coherence, mask, rmse, timeseries
                           unwrapped/wrapped interferograms, coherence]
    geocoding_lookupfile: geocoding look-up table of the master interferogram   

    Example:
         geocode.py velocity.h5   geomap_8rlks.trans
         geocode.py timeseries.h5 geomap_8rlks.trans
         geocode.py LoadedData.h5 geomap_8rlks.trans
         geocode.py Coherence.h5  geomap_8rlks.trans
         geocode.py Wrapped.h5    geomap_8rlks.trans
         geocode.py DEM_error.h5  geomap_8rlks.trans
*******************************************************************
*******************************************************************
   '''


def main(argv):

   try:
       file=argv[0]
       geomap=argv[1]
   except:
       Usage();sys.exit(1)

######################################################################################
   
   fileName=os.path.basename(file).split('.')[0]
   h5file=h5py.File(file,'r')
   k=h5file.keys()

   if 'interferograms' in k: k[0] = 'interferograms'
   elif 'coherence'    in k: k[0] = 'coherence'
   elif 'timeseries'   in k: k[0] = 'timeseries'
   if k[0] in ('interferograms','coherence','wrapped'):
       atr  = h5file[k[0]][h5file[k[0]].keys()[0]].attrs
   elif k[0] in ('dem','velocity','mask','temporal_coherence','rmse','timeseries'):
       atr  = h5file[k[0]].attrs

   #### Subsetted radar coded file
   try:
       x0 = float(atr['subset_x0'])
       y0 = float(atr['subset_y0'])
       print '\nSubsetted radar coded file:\n    creating temporary geomap file for it...'
       rg,az,rsc = readfile.read_float32(geomap)
       rg = rg - x0
       az = az - y0
       geomap = 'temp_'+geomap
       print '    writing '+geomap+'\n'
       writefile.write_float32(rg,az,geomap)
       fg = open(geomap+'.rsc','w')
       for kg in rsc.keys():    fg.write(kg+'    '+rsc[kg]+'\n')
       fg.close()
   except: pass


######################################################################################

   if k[0] in ('velocity','temporal_coherence','mask','rmse','dem'):

       dset = h5file[k[0]].get(k[0])
       data = dset[0:dset.shape[0],0:dset.shape[1]]
       outname=fileName+'.unw'

       amp,unw,unwrsc = geocode_one(data,geomap,outname)

       f = h5py.File('geo_'+file,'w')
       group=f.create_group(k[0])
       dset = group.create_dataset(k[0], data=unw, compression='gzip')
       for key,value in h5file[k[0]].attrs.iteritems():   group.attrs[key]=value
       for key,value in unwrsc.iteritems():               group.attrs[key] = value


######################################################################################

   elif  'timeseries' in k:
       print 'geocoding timeseries:'
       outname='epoch_temp.unw'

       f = h5py.File('geo_'+file,'w')
       group = f.create_group('timeseries')
       epochList=h5file['timeseries'].keys()
       for epoch in epochList:
           print 'geocoding '+epoch
           d = h5file['timeseries'].get(epoch)
           data = d[0:d.shape[0],0:d.shape[1]]

           amp,unw,unwrsc = geocode_one(data,geomap,outname)
           dset = group.create_dataset(epoch, data=unw, compression='gzip')

       for key,value in unwrsc.iteritems():              group.attrs[key] = value
       for key,value in h5file[k[0]].attrs.iteritems():  group.attrs[key] = value
       group.attrs['WIDTH']       = unwrsc['WIDTH'] 
       group.attrs['FILE_LENGTH'] = unwrsc['FILE_LENGTH']
    

######################################################################################

   elif k[0] in ['interferograms','coherence','wrapped']:
       print 'geocoding '+k[0]+' ...'
       if   k[0] == 'interferograms': outname = k[0]+'_temp.unw'
       elif k[0] == 'coherence':      outname = k[0]+'_temp.cor'
       else:                          outname = k[0]+'_temp.int'

       f = h5py.File('geo_'+file,'w')
       gg = f.create_group('interferograms')
       igramList=h5file[k[0]].keys()
       for igram in igramList:
           print 'geocoding '+igram
           d = h5file[k[0]][igram].get(igram)
           data = d[0:d.shape[0],0:d.shape[1]]

           amp,unw,unwrsc = geocode_one(data,geomap,outname)

           group = gg.create_group('geo_'+igram)
           dset = group.create_dataset('geo_'+igram, data=unw, compression='gzip')
           for key,value in unwrsc.iteritems():                         group.attrs[key] = value
           for key,value in h5file[k[0]][igram].attrs.iteritems():      group.attrs[key] = value
           group.attrs['WIDTH'] = unwrsc['WIDTH']
           group.attrs['FILE_LENGTH'] = unwrsc['FILE_LENGTH']

       #######################  support of old format  #######################
       ### mask
       try:
           d = h5file['mask'].get('mask')
           data = d[0:d.shape[0],0:d.shape[1]]
           amp,unw,unwrsc = geocode_one(data,geomap,'mask_'+outname)
           gm = f.create_group('mask')
           dset = gm.create_dataset('mask', data=unw, compression='gzip')
       except:  print 'No group for mask found in the file.'
       ### meanCoherence
       try:
           d = h5file['meanCoherence'].get('meanCoherence')
           data = d[0:d.shape[0],0:d.shape[1]]
           amp,unw,unwrsc = geocode_one(data,geomap,'meanCoherence_'+outname)
           gm = f.create_group('meanCoherence')
           dset = gm.create_dataset('meanCoherence', data=unw, compression='gzip')
       except:  print 'No group for meanCoherence found in the file'


######################################################################################
   try:
       atr['subset_x0']
       rmCmd='rm '+geomap;            os.system(rmCmd);       print rmCmd
       rmCmd='rm '+geomap+'.rsc';     os.system(rmCmd);       print rmCmd
   except: pass

   f.close()
   h5file.close()

######################################################################################

if __name__ == '__main__':

  main(sys.argv[1:])


