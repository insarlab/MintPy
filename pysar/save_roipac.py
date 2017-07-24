#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: update DATE12 for timeseries option
# Yunjun, Oct 2015: add coherence/wrapped option
#                   add two dates option for timeseries


import sys
import os

from numpy import pi
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime


##############################################################################
def usage():
    print '''usage: save_roipac.py  file  [date_info]

Convert PySAR hdf5 file to ROI_PAC format 

argument:
  file : file to be converted.
         for velocity  : the ouput will be a one year interferogram.
         for timeseries: if date is not specified, the last date will be used
                         if two dates are specified, the earlier date will be
                             used as the reference date.

example:
  save_roipac.py  velocity.h5
  save_roipac.py  timeseries.h5    20050601
  save_roipac.py  timeseries.h5    040728    050601
  save_roipac.py  unwrapIfgram.h5  filt_091225-100723-sim_HDR_8rlks_c10.unw
  save_roipac.py  unwrapIfgram.h5  091225-100723
  save_roipac.py  temporal_coherence.h5
    '''
    return


##############################################################################
def main(argv):
    try:
        File=argv[0]
    except:
        usage();sys.exit(1)
  
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    atr['PROCESSOR'] = 'roipac'
  
    h5file = h5py.File(File,'r')
  
    if k == 'velocity':
        dset = h5file['velocity'].get('velocity')
        data = dset[0:dset.shape[0],0:dset.shape[1]]
        print "converting velocity to a 1 year interferogram."
        wvl=float(h5file[k].attrs['WAVELENGTH'])
        data=(-4*pi/wvl)*data

        outname=File.split('.')[0]+'.unw'
        print 'writing >>> '+outname
        writefile.write(data,atr,outname)

    elif k == 'timeseries':
        dateList=h5file['timeseries'].keys() 
        ## Input
        if   len(sys.argv)==2:
            print 'No input date specified >>> continue with the last date'
            dateList=h5file['timeseries'].keys()
            d=dateList[-1]
        elif len(sys.argv)==3:
            d=sys.argv[2]
        elif len(sys.argv)==4:
            ds=sorted(sys.argv[2:4])
            d_ref = ds[0]
            d     = ds[1]
        else: usage(); sys.exit(1)
        d = ptime.yyyymmdd(d)
        try: d_ref = ptime.yyyymmdd(d_ref)
        except: pass

        ## Data
        print 'reading '+d+' ... '
        data = h5file['timeseries'].get(d)[:]
        try:
            print 'reading '+d_ref+' ... '
            data_ref = h5file['timeseries'].get(d_ref)[:]
            data = data - data_ref
        except: pass
        wvl=float(atr['WAVELENGTH'])
        data *= -4*pi/wvl

        ## outName
        try:      master_d = d_ref
        except:
            try:    master_d = atr['ref_date']
            except: master_d = dateList[0]
        if len(master_d)==8:  master_d=master_d[2:8]
        if len(d)==8:         d=d[2:8]
        outname = master_d+'_'+d+'.unw'

        ## Attributes
        atr['FILE_TYPE']             = '.unw'
        atr['P_BASELINE_TIMESERIES'] = '0.0'
        atr['UNIT']                  = 'radian'
        atr['DATE']                  = master_d
        atr['DATE12']                = master_d+'-'+d

        ## Writing
        print 'writing >>> '+outname
        writefile.write(data,atr,outname)

    elif k in ['interferograms','coherence','wrapped']:
        ## Check input
        igramList=h5file[k].keys()
        try:
            d = sys.argv[2]
            for i in range(len(igramList)):
                if d in igramList[i]:
                    igram = igramList[i]
        except:
            igram = igramList[-1]
            print 'No input date specified >>> continue with the last date'

        ## Read and Write
        print 'reading '+igram+' ... '
        data = h5file[k][igram].get(igram)[:]
        atr = h5file[k][igram].attrs
        atr['PROCESSOR'] = 'roipac'
        outname = igram

        print 'writing >>> '+ outname
        writefile.write(data, atr, outname)  

    else:
        dset = h5file[k].get(k)
        data = dset[0:dset.shape[0],0:dset.shape[1]]
        if k in ['temporal_coherence']:
            outname=File.split('.')[0]+'.cor'
        elif k in ['dem','.hgt','.dem']:
            atr['FILE_TYPE'] = '.dem'
            outname=os.path.splitext(File)[0]+'.dem'
        else:
            outname=File.split('.')[0]+'.unw'
        print 'writing >>> '+ outname
        writefile.write(data,atr,outname)


    h5file.close()
    return


##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
