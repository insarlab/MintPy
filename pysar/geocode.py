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
import numpy as np
#from joblib import Parallel, delayed
#import multiprocessing

import pysar._readfile  as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


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
        writefile.write_roipac_rsc(rsc, geomap_file+'.rsc')

    return geomap_file


######################  Geocode one data  ########################
def geocode_data_roipac(data, atr, geomapFile, roipac_name):
    '''Geocode input data with attribute dict using geomapFile
    Inputs:
        data : 2D np.array
        atr  : dict, attributes
        geomapFile  : string, path of geomap*.trans file
        roipac_name : string, path of roipac file
    '''

    print 'writing to roi_pac unw file format'
    writefile.write_float32(data, roipac_name)
    writefile.write_roipac_rsc(atr, roipac_name+'.rsc')
 
    geoCmd = 'geocode.pl '+geomapFile+' '+roipac_name+' geo_'+roipac_name
    os.system(geoCmd)
    print geoCmd
 
    print 'reading geocoded file...'
    amp,unw,unwrsc = readfile.read_float32('geo_'+roipac_name)
 
    rmCmd = 'rm '+roipac_name+' '+roipac_name+'.rsc';           os.system(rmCmd);       print rmCmd
    rmCmd = 'rm geo_'+roipac_name+' geo_'+roipac_name+'.rsc';   os.system(rmCmd);       print rmCmd
 
    return amp, unw, unwrsc


######################################################################################
def geocode_attribute(atr_rdr, atr_geo, transFile=None):
    '''Update attributes after geocoding'''
    atr = dict()
    for key, value in atr_geo.iteritems():  atr[key] = str(value)
    for key, value in atr_rdr.iteritems():  atr[key] = str(value)
    atr['WIDTH']       = atr_geo['WIDTH']
    atr['FILE_LENGTH'] = atr_geo['FILE_LENGTH']
    atr['YMIN'] = str(0)
    atr['YMAX'] = str(int(atr['FILE_LENGTH'])-1)
    atr['XMIN'] = str(0)
    atr['XMAX'] = str(int(atr['WIDTH'])-1)

    # Reference point from y/x to lat/lon
    if transFile and ('ref_x' and 'ref_y' in atr_rdr.keys()):
        ref_x = np.array(int(atr_rdr['ref_x']))
        ref_y = np.array(int(atr_rdr['ref_y']))
        ref_lat, ref_lon = ut.radar2glob(ref_y, ref_x, transFile, atr_rdr)[0:2]
        atr['ref_lat'] = ref_lat
        atr['ref_lon'] = ref_lon
        print 'update reference point info in lat/lon'
    return atr


def geocode_file_roipac(infile, geomap_file, outfile=None):
    '''Geocode one file'''
    # Input file info
    atr = readfile.read_attribute(infile)
    k = atr['FILE_TYPE']
    print 'geocoding '+k+' file: '+infile+' ...'

    # roipac outfile name info - intermediate product
    infile_base, ext = os.path.splitext(infile)
    infile_mark = infile_base+'_'+ext.split('.')[1]
    if k in ['coherence','temporal_coherence','.cor']:
        roipac_ext = '.cor'
    elif k in ['wrapped','.int']:
        roipac_ext = '.int'
    else:
        roipac_ext = '.unw'
     
    # temporary geomap file - needed for parallel processing
    geomap_file2 = geomap_file.split('.trans')[0]+'4'+infile_mark+'.trans'
    cpCmd = 'cp '+geomap_file+' '+geomap_file2;              os.system(cpCmd);  print cpCmd
    cpCmd = 'cp '+geomap_file+'.rsc '+geomap_file2+'.rsc';   os.system(cpCmd);  print cpCmd
    
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
                
                roipac_name = infile_mark+'_'+epoch+roipac_ext
                geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, atr, geomap_file2, roipac_name)
                geo_atr = geocode_attribute(atr, geo_rsc, geomap_file2)
                
                gg = group.create_group('geo_'+epoch)
                dset = gg.create_dataset('geo_'+epoch, data=geo_data, compression='gzip')
                for key, value in geo_atr.iteritems():
                    gg.attrs[key] = value

        elif k in ['timeseries']:
            for epoch in epochList:
                print epoch
                data = h5[k].get(epoch)[:]
                
                roipac_name = infile_mark+'_'+epoch+roipac_ext
                geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, atr, geomap_file2, roipac_name)
                
                dset = group.create_dataset(epoch, data=geo_data, compression='gzip')
            geo_atr = geocode_attribute(atr, geo_rsc, geomap_file2)
            for key, value in geo_atr.iteritems():
                group.attrs[key] = value
                
        h5.close()
        h5out.close()

    # Single-dataset file
    elif atr['PROCESSOR'] == 'roipac':
        rmCmd = 'rm geo_'+infile+' geo_'+infile+'.rsc'
        os.system(rmCmd)
        print rmCmd

        geoCmd = 'geocode.pl '+geomap_file2+' '+infile+' geo_'+infile
        os.system(geoCmd)
        print geoCmd

        atr_rdr = readfile.read_roipac_rsc(infile+'.rsc')
        atr_geo = readfile.read_roipac_rsc('geo_'+infile+'.rsc')
        atr_geo = geocode_attribute(atr_rdr, atr_geo)
        writefile.write_roipac_rsc(atr_geo, 'geo_'+infile+'.rsc')

    else:
        rmCmd = 'rm '+outfile+' '+outfile+'.rsc';    os.system(rmCmd);    print rmCmd
        data, atr = readfile.read(infile)

        roipac_name = infile_mark+roipac_ext
        geo_amp, geo_data, geo_rsc = geocode_data_roipac(data, atr, geomap_file2, roipac_name)
        geo_atr = geocode_attribute(atr, geo_rsc, geomap_file2)
        
        writefile.write(geo_data, geo_atr, outfile)

    # delete temporary geomap file
    rmCmd='rm '+geomap_file2+' '+geomap_file2+'.rsc'
    print rmCmd
    os.system(rmCmd)

    return outfile


######################################################################################
EXAMPLE='''example:
  geocode.py  geomap_8rlks.trans  velocity.py
  geocode.py  geomap_8rlks.trans  *velocity*h5
  geocode.py  geomap_8rlks.trans  timeseries_ECMWF_demCor.h5 velocity_ex.h5
  geocode.py  geomap_8rlks.trans  100901-*.cor
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
    parser.add_argument('--parallel',dest='parallel',action='store_true',\
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    
    inps = parser.parse_args()
    return inps


######################################################################################
def main(argv):
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)
    print 'number of file to mask: '+str(len(inps.file))
    print inps.file

    if not ut.which('geocode.pl'):
        sys.exit("\nERROR: Can not find geocode.pl, it's needed for geocoding.\n")
    
    #print '\n***************** Geocoding *******************'
    if not inps.lookup_file.endswith('.trans'):
        print 'ERROR: Input lookup file is not .trans file: '+inps.lookup_file+'\n'
        sys.exit(1)

    # Check geomap file for previously subsetted radar coord file
    atr = readfile.read_attribute(inps.file[0])
    if 'subset_x0' in atr.keys():
        inps.lookup_file = geomap4subset_radar_file(atr, inps.lookup_file)

    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))

    # Geocoding
    if len(inps.file) == 1:
        geocode_file_roipac(inps.file[0], inps.lookup_file, inps.outfile)

    elif inps.parallel:
        #num_cores = min(multiprocessing.cpu_count(), len(inps.file), pysar.parallel_num)
        #print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(geocode_file_roipac)(file, inps.lookup_file) for file in inps.file)
    else:
        for File in inps.file:
            print '----------------------------------------------------'
            geocode_file_roipac(File, inps.lookup_file)

    # clean temporary geomap file for previously subsetted radar coord file
    if 'subset_x0' in atr.keys():
        rmCmd='rm '+inps.lookup_file+' '+inps.lookup_file+'.rsc'
        os.system(rmCmd)
        print rmCmd

    print 'Done.'
    return

######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


