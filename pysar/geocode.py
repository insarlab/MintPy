#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Finish 'interferograms' option, Yunjun, Jun 2015
# Add 'coherence'/'wrapped' option, Yunjun, Jul 2015
#

import h5py
import _writefile as writefile
import sys
import os
import _readfile as readfile
import numpy as np


def Usage():


    print '''
*******************************************************************
*******************************************************************
    
    uses roi_pac geocoding functions to geocode the PySAR products

    usage:

         geocode.py File geocoding_lookupfile

    File: hdf5 file [velocity, temporal_coherence, mask, rmse, timeseries, unwrapped/wrapped interferograms, coherence]
    geocoding_lookupfile: geocoding look-up table of the master interferogram   

    Example:
         geocode.py velocity.h5   geomap_8rlks.trans
         geocode.py timeseries.h5 geomap_8rlks.trans
         geocode.py LoadedData.h5 geomap_8rlks.trans
         geocode.py Coherence.h5  geomap_8rlks.trans
         geocode.py Wrapped.h5    geomap_8rlks.trans

*******************************************************************
*******************************************************************
    '''

def main(argv):

   try:
      file=argv[0]
      geomap=argv[1]
   except:
      Usage();sys.exit(1)

   

   fileName=os.path.basename(file).split('.')[0]
   h5file=h5py.File(file,'r')
   k=h5file.keys()
   if k[0] in ('velocity','temporal_coherence','mask','rmse') and 'timeseries' not in k:
      
       
       dset = h5file[k[0]].get(k[0])
       data = dset[0:dset.shape[0],0:dset.shape[1]]
       outname=fileName+'.unw'
       print 'writing to roi_pac unw file format'
       writefile.write_float32(data,outname)     
       f = open(outname+'.rsc','w')
       f.write('FILE_LENGTH       '+str(data.shape[0])+'\n')
       f.write('WIDTH             '+str(data.shape[1])+'\n')
       f.close()

       geoCmd='geocode.pl '+geomap+' '+outname+' geo_'+outname
       print geoCmd
       os.system(geoCmd)

       print 'reading geocoded file and write it to h5 format'
       amp,unw,unwrsc = readfile.read_float32('geo_'+outname)

       rmCmd='rm '+outname;                 os.system(rmCmd);       print rmCmd
       rmCmd='rm '+outname+'.rsc';          os.system(rmCmd);       print rmCmd
       rmCmd='rm geo_'+outname;             os.system(rmCmd);       print rmCmd
       rmCmd='rm geo_'+outname+'.rsc';      os.system(rmCmd);       print rmCmd

       f = h5py.File('geo_'+file,'w')
       group=f.create_group(k[0])
       dset = group.create_dataset(k[0], data=unw, compression='gzip')
       for key , value in h5file[k[0]].attrs.iteritems():
           group.attrs[key]=value
       
       for key,value in unwrsc.iteritems():
          group.attrs[key] = value

       f.close()
       h5file.close()

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
           writefile.write_float32(data,outname)

           f = open(outname+'.rsc','w')
           f.write('FILE_LENGTH       '+str(data.shape[0])+'\n')
           f.write('WIDTH             '+str(data.shape[1])+'\n')
           f.close()

           geoCmd='geocode.pl '+geomap+' '+outname+' geo_'+outname
           print geoCmd;	os.system(geoCmd)

           print 'reading geocoded file and add it to '+'geo_'+file
           amp,unw,unwrsc = readfile.read_float32('geo_'+outname)

           dset = group.create_dataset(epoch, data=unw, compression='gzip')

           rmCmd='rm '+outname;                 os.system(rmCmd);       print rmCmd
           rmCmd='rm '+outname+'.rsc';          os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname;             os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname+'.rsc';      os.system(rmCmd);       print rmCmd
        
       for key,value in unwrsc.iteritems():
          group.attrs[key] = value

       for key,value in h5file['timeseries'].attrs.iteritems():
          group.attrs[key] = value
       
       group.attrs['WIDTH'] =unwrsc['WIDTH'] 
       group.attrs['FILE_LENGTH'] =unwrsc['FILE_LENGTH']
    
       f.close()
       h5file.close()
       
   elif k[0] == 'interferograms':
       print 'geocoding interferograms:'
       outname='igram_temp.unw'

       f = h5py.File('geo_'+file,'w')
       gg = f.create_group('interferograms')

       igramList=h5file[k[0]].keys()
       for igram in igramList:
           print 'geocoding '+igram
           group = gg.create_group('geo_'+igram)

           d = h5file['interferograms'][igram].get(igram)
           data = d[0:d.shape[0],0:d.shape[1]]
           writefile.write_float32(data,outname)

           f_temp = open(outname+'.rsc','w')
           f_temp.write('FILE_LENGTH	'+str(data.shape[0])+'\n')
           f_temp.write('WIDTH	'+str(data.shape[1])+'\n')
           f_temp.close()

           geoCmd='geocode.pl '+geomap+' '+outname+' geo_'+outname
           print geoCmd;	os.system(geoCmd)

           print 'reading geocoded file and add it to '+'geo_'+file
           amp,unw,unwrsc = readfile.read_float32('geo_'+outname)
           if igram==igramList[0]:
               MaskZero=np.ones([unw.shape[0],unw.shape[1]])
           MaskZero=amp*MaskZero

           dset = group.create_dataset('geo_'+igram, data=unw, compression='gzip')

           rmCmd='rm '+outname;                 os.system(rmCmd);       print rmCmd
           rmCmd='rm '+outname+'.rsc';          os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname;             os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname+'.rsc';      os.system(rmCmd);       print rmCmd

           for key,value in unwrsc.iteritems():
               group.attrs[key] = value

           for key,value in h5file['interferograms'][igram].attrs.iteritems():
               group.attrs[key] = value
       
#           for key,value in unwrsc.iteritems():
#               group.attrs[key] = value
   
           group.attrs['WIDTH'] = unwrsc['WIDTH']
           group.attrs['FILE_LENGTH'] = unwrsc['FILE_LENGTH']

       Mask=np.ones(MaskZero.shape)
       Mask[MaskZero==0]=0
       gm = f.create_group('mask')
       dset = gm.create_dataset('mask', data=Mask, compression='gzip')

       f.close()
       h5file.close()


   elif k[0] == 'coherence':
       print 'geocoding coherence:'
       outname='cor_temp.unw'

       f = h5py.File('geo_'+file,'w')
       gg = f.create_group(k[0])
       corList=h5file[k[0]].keys()
       for cor in corList:
           print 'geocoding '+cor
           group = gg.create_group('geo_'+cor)
           d = h5file[k[0]][cor].get(cor)
           data = d[0:d.shape[0],0:d.shape[1]]
           writefile.write_float32(data,outname)

           f_temp = open(outname+'.rsc','w')
           f_temp.write('FILE_LENGTH    '+str(data.shape[0])+'\n')
           f_temp.write('WIDTH  '+str(data.shape[1])+'\n')
           f_temp.close()
           geoCmd='geocode.pl '+geomap+' '+outname+' geo_'+outname
           print geoCmd;	os.system(geoCmd)

           print 'reading geocoded file and add it to '+'geo_'+file
           amp,unw,unwrsc = readfile.read_float32('geo_'+outname)
           dset = group.create_dataset('geo_'+cor, data=unw, compression='gzip')
           
           rmCmd='rm '+outname;                 os.system(rmCmd);       print rmCmd
           rmCmd='rm '+outname+'.rsc';          os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname;             os.system(rmCmd);       print rmCmd
           rmCmd='rm geo_'+outname+'.rsc';      os.system(rmCmd);       print rmCmd

           for key,value in unwrsc.iteritems():
               group.attrs[key] = value
           for key,value in h5file[k[0]][cor].attrs.iteritems():
               group.attrs[key] = value
           group.attrs['WIDTH'] = unwrsc['WIDTH']
           group.attrs['FILE_LENGTH'] = unwrsc['FILE_LENGTH']

       f.close()
       h5file.close()

   elif k[0] == 'wrapped':
       print 'geocoding wrapped interferograms:'
       outname='wrap_temp.int'

       f = h5py.File('geo_'+file,'w')
       gg = f.create_group(k[0])
       wrapList=h5file[k[0]].keys()
       for wrap in wrapList:
           print 'geocoding '+wrap
           group = gg.create_group('geo_'+wrap)
           d = h5file[k[0]][wrap].get(wrap)
           data = d[0:d.shape[0],0:d.shape[1]]
           writefile.write_complex64(data,outname)

           f_temp = open(outname+'.rsc','w')
           f_temp.write('FILE_LENGTH    '+str(data.shape[0])+'\n')
           f_temp.write('WIDTH  '+str(data.shape[1])+'\n')
           f_temp.close()
           geoCmd='geocode.pl '+geomap+' '+outname+' geo_'+outname
           print geoCmd;	os.system(geoCmd)

           print 'reading geocoded file and add it to '+'geo_'+file
           amp,unw,unwrsc = readfile.read_complex64('geo_'+outname)
           dset = group.create_dataset('geo_'+wrap, data=unw, compression='gzip')

           rmCmd='rm '+outname;        		os.system(rmCmd);	print rmCmd
           rmCmd='rm '+outname+'.rsc'; 		os.system(rmCmd);	print rmCmd
           rmCmd='rm geo_'+outname;		os.system(rmCmd);	print rmCmd
           rmCmd='rm geo_'+outname+'.rsc';	os.system(rmCmd);	print rmCmd

           for key,value in unwrsc.iteritems():
               group.attrs[key] = value
           for key,value in h5file[k[0]][wrap].attrs.iteritems():
               group.attrs[key] = value
           group.attrs['WIDTH'] = unwrsc['WIDTH']
           group.attrs['FILE_LENGTH'] = unwrsc['FILE_LENGTH']

       f.close()
       h5file.close()


if __name__ == '__main__':

  main(sys.argv[1:])


