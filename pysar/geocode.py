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
# Yunjun, Jun 2016: Add geocode_attributes(), use read() and write() for file IO


import os
import sys

import h5py

import pysar._readfile  as readfile
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


######################################################################################
def geocode_attributes(atr_rdr,atr_geo):
   atr = dict()
   for key, value in atr_geo.iteritems():  atr[key] = str(value)
   for key, value in atr_rdr.iteritems():  atr[key] = value
   atr['WIDTH']       = atr_geo['WIDTH']
   atr['FILE_LENGTH'] = atr_geo['FILE_LENGTH']

   return atr


######################################################################################
def Usage():
   print '''
*******************************************************************
*******************************************************************
    
    uses roi_pac geocoding functions to geocode the PySAR products

    Usage:
         geocode.py File geocoding_lookupfile

    File: PySAR hdf5 file, including subseted.
                [velocity, temporal_coherence, mask, rmse, timeseries
                 unwrapped/wrapped interferograms, coherence]
    geocoding_lookupfile: geocoding look-up table of the master interferogram   

    Example:
         geocode.py velocity.h5        geomap_8rlks.trans
         geocode.py subset_velocity.h5 geomap_8rlks.trans
         geocode.py timeseries.h5      geomap_8rlks.trans
         geocode.py LoadedData.h5      geomap_8rlks.trans
         geocode.py Coherence.h5       geomap_8rlks.trans
         geocode.py Wrapped.h5         geomap_8rlks.trans
         geocode.py DEM_error.h5       geomap_8rlks.trans

*******************************************************************
*******************************************************************
   '''


######################################################################################
def main(argv):

   try:
       file=argv[0]
       geomap=argv[1]
   except:
       Usage();sys.exit(1)

   ######################################################################################
   fileName=os.path.basename(file).split('.')[0]
   h5file=h5py.File(file,'r')
   atr = readfile.read_attributes(file)
   k = atr['FILE_TYPE']
   print 'geocoding '+k

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
   if k in ['timeseries']:
       outname='epoch_temp.unw'

       f = h5py.File('geo_'+file,'w')
       group = f.create_group('timeseries')
       epochList=h5file['timeseries'].keys()
       for epoch in epochList:
           print 'geocoding '+epoch
           data = h5file['timeseries'].get(epoch)[:]

           amp,unw,unwrsc = geocode_one(data,geomap,outname)
           dset = group.create_dataset(epoch, data=unw, compression='gzip')

       atr = geocode_attributes(atr,unwrsc)
       for key,value in atr.iteritems():
           group.attrs[key] = value

   ######################################################################################
   elif k in ['interferograms','coherence','wrapped']:
       if   k == 'interferograms': outname = k[0]+'_temp.unw'
       elif k == 'coherence'     : outname = k[0]+'_temp.cor'
       else:                       outname = k[0]+'_temp.int'

       f = h5py.File('geo_'+file,'w')
       gg = f.create_group('interferograms')
       igramList=h5file[k].keys()
       for igram in igramList:
           print 'geocoding '+igram
           data = h5file[k][igram].get(igram)[:]

           amp,unw,unwrsc = geocode_one(data,geomap,outname)

           group = gg.create_group('geo_'+igram)
           dset = group.create_dataset('geo_'+igram, data=unw, compression='gzip')

           atr = geocode_attributes(h5file[k][igram].attrs, unwrsc)
           for key,value in atr.iteritems():
               group.attrs[key] = value

       #######################  support of old format  #######################
       ### mask
       try:
           data = h5file['mask'].get('mask')[:]
           amp,unw,unwrsc = geocode_one(data,geomap,'mask_'+outname)
           gm = f.create_group('mask')
           dset = gm.create_dataset('mask', data=unw, compression='gzip')
       except:  print 'No group for mask found in the file.'
       ### meanCoherence
       try:
           data = h5file['meanCoherence'].get('meanCoherence')[:]
           amp,unw,unwrsc = geocode_one(data,geomap,'meanCoherence_'+outname)
           gm = f.create_group('meanCoherence')
           dset = gm.create_dataset('meanCoherence', data=unw, compression='gzip')
       except:  print 'No group for meanCoherence found in the file'

   ######################################################################################
   else:
       data,atr = readfile.read(file)
       outname=fileName+'.unw'

       amp,unw,unwrsc = geocode_one(data,geomap,outname)
       atr = geocode_attributes(atr,unwrsc)

       writefile.write(unw,atr,'geo_'+file)


   ######################################################################################
   try:
       atr['subset_x0']
       rmCmd='rm '+geomap;            os.system(rmCmd);       print rmCmd
       rmCmd='rm '+geomap+'.rsc';     os.system(rmCmd);       print rmCmd
   except: pass

   try:
       f.close()
       h5file.close()
   except: pass


######################################################################################
if __name__ == '__main__':
  main(sys.argv[1:])


