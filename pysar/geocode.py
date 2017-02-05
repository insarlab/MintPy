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
# Yunjun, Oct 2015: Add geocode_data()
#                   Merge 'interferograms','coherence','wrapped' into one
#                   Add support for subsetted radar coded files
# Yunjun, Jun 2016: Add geocode_attribute(), use read() and write() for file IO
# Yunjun, Jan 2017: add geocode_file_roipac(), parallel and cmdLineParse()


import os
import sys
import argparse

import h5py
from joblib import Parallel, delayed
import multiprocessing

import pysar._readfile  as readfile
import pysar._writefile as writefile
from pysar._pysar_utilities import get_file_list


def geomap4subset_radar_file(radar_atr, geomap_file):
    ''' Add offset value to geomap file if input radar file has been subsetted.'''
    if 'subset_x0' in radar_atr.keys():
        x0 = float(radar_atr['subset_x0'])
        y0 = float(radar_atr['subset_y0'])
        print '\nInput radar coord file has been subsetted.\n    creating temporary geomap file for it...'

        rg,az,rsc = readfile.read_float32(geomap_file)
        rg = rg - x0
        az = az - y0

        geomap_file = 'temp_'+geomap_file
        print '    writing >>> '+geomap_file+'\n'
        writefile.write_float32(rg, az, geomap_file)
        f = open(geomap+'.rsc','w')
        for key in rsc.keys():
            f.write(kg+'    '+rsc[key]+'\n')
        f.close()
    return geomap_file


######################  Geocode one data  ########################
def geocode_data_roipac(data,geomapFile,outname):

    print 'writing to roi_pac unw file format'
    writefile.write_float32(data, outname)
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
def geocode_attribute(atr_rdr,atr_geo):
    atr = dict()
    for key, value in atr_geo.iteritems():  atr[key] = str(value)
    for key, value in atr_rdr.iteritems():  atr[key] = value
    atr['WIDTH']       = atr_geo['WIDTH']
    atr['FILE_LENGTH'] = atr_geo['FILE_LENGTH']

    return atr


def geocode_file_roipac(infile, geomap_file, outfile=None):
    '''Geocode one file'''
    # Input file info
    atr = readfile.read_attribute(infile)
    k = atr['FILE_TYPE']
    print 'geocoding '+k+' file '+infile+' ...'

    # roipac outfile name info - intermediate product
    ext = os.path.splitext(infile)[1]
    infile_base = os.path.basename(infile).split(ext)[0]
    if k == 'coherence':  roipac_ext = '.cor'
    elif k == 'wrapped':  roipac_ext = '.int'
    else:  roipac_ext = '.unw'
     
    # temporary geomap file - needed for parallel processing
    geomap_file_orig = geomap_file
    geomap_file = geomap_file_orig.split('.trans')[0]+'4'+infile_base+'.trans'
    cpCmd='cp '+geomap_file_orig+' '+geomap_file;              os.system(cpCmd);  print cpCmd
    cpCmd='cp '+geomap_file_orig+'.rsc '+geomap_file+'.rsc';   os.system(cpCmd);  print cpCmd
    
    # Output file name
    if not outfile:
        outfile = 'geo_'+infile
    print 'writing >>> '+outfile

    # Multi-dataset file
    if k in ['timeseries','interferograms','coherence','wrapped']:
        h5 = h5py.File(infile, 'r')
        epochList = sorted(h5[k].keys())
        print 'number of epochs: '+str(len(epochList))
        
        h5out = h5py.File(outfile, 'w')
        group = h5out.create_group(k)
        
        if k in ['interferograms','coherence','wrapped']:
            for epoch in epochList:
                print epoch
                data = h5[k][epoch].get(epoch)[:]
                atr = h5[k][epoch].attrs
                
                roipac_outname = infile_base+'_'+epoch+roipac_ext
                geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, geomap_file, roipac_outname)
                geo_atr = geocode_attribute(atr, geo_rsc)
                
                gg = group.create_group('geo_'+epoch)
                dset = gg.create_dataset('geo_'+epoch, data=geo_data, compression='gzip')
                for key, value in geo_atr.iteritems():
                    gg.attrs[key] = value

        elif k in ['timeseries']:
            for epoch in epochList:
                print epoch
                data = h5[k].get(epoch)[:]
                
                roipac_outname = infile_base+'_'+epoch+roipac_ext
                geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, geomap_file, roipac_outname)
                
                dset = group.create_dataset(epoch, data=geo_data, compression='gzip')
            geo_atr = geocode_attribute(atr, geo_rsc)
            for key, value in geo_atr.iteritems():
                group.attrs[key] = value
                
        h5.close()
        h5out.close()
                
    # Single-dataset file
    else:
        data, atr = readfile.read(infile)
        roipac_outname = infile_base+roipac_ext
        
        geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, geomap_file, roipac_outname)
        geo_atr = geocode_attribute(atr, geo_rsc)
        
        writefile.write(geo_data, geo_atr, outfile)

    # delete temporary geomap file
    rmCmd='rm '+geomap_file;         os.system(rmCmd);   print rmCmd
    rmCmd='rm '+geomap_file+'.rsc';  os.system(rmCmd);   print rmCmd

    return outfile


######################################################################################
EXAMPLE='''example:
  geocode.py  geomap_8rlks.trans  velocity.py
  geocode.py  geomap_8rlks.trans  *velocity*h5
  geocode.py  geomap_8rlks.trans  timeseries_ECMWF_demCor.h5 velocity_ex.h5
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Geocode PySAR products using roi_pac geocoding function',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('lookup_file', \
                        help='geocoding look-up file.\n'
                             'i.e. geomap_*rlks.trans for roi_pac product')
    parser.add_argument('file', nargs='+', help='File(s) to be geocoded')
    parser.add_argument('-o','--outfile', help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    
    inps = parser.parse_args()
    return inps


######################################################################################
def main(argv):
    inps = cmdLineParse()
    inps.file = get_file_list(inps.file)
    print '\n***************** Geocoding *******************'
    if not inps.lookup_file.endswith('.trans'):
        print 'ERROR: Input lookup file is not .trans file: '+inps.lookup_file+'\n'
        sys.exit(1)
    print 'number of file to mask: '+str(len(inps.file))
    print inps.file
    
    # check outfile and parallel option
    if len(inps.file) > 1:
        inps.outfile = None
    elif len(inps.file) == 1 and inps.parallel:
        inps.parallel =  False
        print 'parallel processing is diabled for one input file'

    # Check geomap file for previously subsetted radar coord file
    atr = readfile.read_attribute(inps.file[0])
    if 'subset_x0' in atr.keys():
        inps.lookup_file = geomap4subset_radar_file(atr, inps.lookup_file)

    # Geocode files(s)
    print '----------------------------------------------------'
    if inps.parallel:
        num_cores = multiprocessing.cpu_count()
        print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(geocode_file_roipac)(file, inps.lookup_file) for file in inps.file)
    else:
        geocode_file_roipac(inps.file[0], inps.lookup_file, inps.outfile)

    # clean temporary geomap file for previously subsetted radar coord file
    if 'subset_x0' in atr.keys():
        rmCmd='rm '+geomap;            os.system(rmCmd);       print rmCmd
        rmCmd='rm '+geomap+'.rsc';     os.system(rmCmd);       print rmCmd

    print 'Done.'

######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


