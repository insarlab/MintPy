#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Dec 2015: add support for ROI_PAC product

import _gmt
import h5py
from numpy import linspace,meshgrid,flipud
import sys
import os


def get_geo_lat_lon(atr):
    X_FIRST=float(atr['X_FIRST'])
    Y_FIRST=float(atr['Y_FIRST'])
    X_STEP =float(atr['X_STEP'])
    Y_STEP =float(atr['Y_STEP'])
    W      =float(atr['WIDTH'])
    L      =float(atr['FILE_LENGTH'])
    Y_END = Y_FIRST+(L-1)*Y_STEP
    X_END = X_FIRST+(W-1)*X_STEP
    X     = linspace(X_FIRST,X_END,W)
    Y     = linspace(Y_FIRST,Y_END,L)
    XI,YI = meshgrid(X,Y)
 
    return Y,X

def Usage():
    print '''
   ******************************************************
   Exporting geocoded pysar velocity file to GMT grd file
   It also exports an epoch of timeseries to the grd file 

   Example:

       save_gmt.py  geo_velocity.h5
       save_gmt.py  geo_timeseries.h5 20071031
       save_gmt.py  geo_timeseries.h5
       save_gmt.py  geo_filt_100608-101024-sim_HDR_16rlks_c10.unw
       save_gmt.py  gsi10m.dem
   *****************************************************
    '''


###############################  Main Function  ####################################
def main(argv):
    try: file = argv[0]
    except: Usage(); sys.exit(1)
  
    print '\n*************** Save to GRD file for GMT ****************'
    print 'Input file: '+file
    ext = os.path.splitext(file)[1]
  
    ########## PySAR HDF5 Files ################
    if ext == '.h5':
        h5=h5py.File(file,'r')
        k=h5.keys()
        if 'interferograms' in k: k[0] = 'interferograms'
        elif 'coherence'    in k: k[0] = 'coherence'
        elif 'timeseries'   in k: k[0] = 'timeseries'
        if k[0] in ('interferograms','coherence','wrapped'):
            atr  = h5[k[0]][h5[k[0]].keys()[0]].attrs
        elif k[0] in ('dem','velocity','mask','temporal_coherence','rmse','timeseries'):
            atr  = h5[k[0]].attrs
   
        if 'timeseries' in k:
            try:
                d=sys.argv[2]
            except:
                print 'No input date... continue to convert the last date of timeseries'
                dateList=h5['timeseries'].keys()
                d=dateList[-1]
            print 'reading '+d + ' ... '
            dset=h5['timeseries'].get(d)
   
        elif k[0] in ['interferograms','coherence','wrapped']:
           print 'interferograms is not supported currently.'; sys.exit(1)
   
        else:   dset=h5[k[0]].get(k[0])
   
        z=dset[0:dset.shape[0],0:dset.shape[1]]
  
    ########## ROI_PAC Files ##################
    elif ext in ['.unw','.cor','.hgt','.dem']:
        import pysar._readfile as readfile
        if ext == '.dem' :                    z,atr = readfile.read_dem(file)
        if ext in ['.unw','.cor','.hgt']:   a,z,atr = readfile.read_float32(file)
  
    else: print 'Unsupported file extension: '+ext;  sys.exit(1)
  
    lats,lons=get_geo_lat_lon(atr)
  
    #fname='velocity.grd'
    outName=os.path.basename(file).replace(ext,'')+'.grd';  print 'Output (GMT grd file): '+outName
    _gmt.write_gmt_simple(lons, flipud(lats), flipud(z), outName, title='default', name='velocity',\
                          scale=1.0, offset=0, units='meters')
    #_gmt.write_gmt_simple(lons, lats, res, fname, title=title, name=name, units=units)


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


