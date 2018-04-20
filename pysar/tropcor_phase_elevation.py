#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys
import time, datetime
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pysar.utils import readfile, datetime as ptime, utils as ut


############################################################################
EXAMPLE='''example:
  tropcor_phase_elevation.py  timeseries_demErr.h5
  tropcor_phase_elevation.py  timeseries_demErr.h5      -d demRadar.h5  -m maskTempCoh.h5      -p 1
  tropcor_phase_elevation.py  geo_timeseries_demErr.h5  -d demGeo.h5    -m geo_maskTempCoh.h5  -p 1
'''

REFERENCE='''reference:
  Doin, M. P., C. Lasserre, G. Peltzer, O. Cavalie, and C. Doubre (2009), Corrections of stratified 
  tropospheric delays in SAR interferometry: Validation with global atmospheric models, J App. Geophy.,
  69(1), 35-50, doi:http://dx.doi.org/10.1016/j.jappgeo.2009.03.010.
'''

def createParser():
    parser = argparse.ArgumentParser(description='Stratified tropospheric delay correction using height-correlation approach',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='time-series file to be corrected')
    parser.add_argument('-d','--dem', dest='dem_file', help='DEM file used for correlation calculation.')
    parser.add_argument('-t','--threshold', type=float,\
                        help='correlation threshold to apply phase correction.\n'+\
                             'if not set, all dates will be corrected.')
    parser.add_argument('-m','--mask', dest='mask_file', help='mask file for pixels used for correlation calculation')
    parser.add_argument('--poly-order','-p', dest='poly_order', type=int, default=1, choices=[1,2,3],\
                        help='polynomial order of phase-height correlation. Default: 1')
    parser.add_argument('-o','--outfile', help='output corrected timeseries file name')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    inps = parser.parse_args(args=iargs)

    if inps.threshold and (not 0.0 <= inps.threshold <= 1.0):
        raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % inps.threshold)
    return inps


############################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)

    ##### Check default input arguments
    # default output filename
    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.timeseries_file)[0]+'_tropHgt.h5'

    # Basic info
    atr = readfile.read_attribute(inps.timeseries_file)
    k = atr['FILE_TYPE']
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    pix_num = length*width

    # default DEM file
    if not inps.dem_file:
        if 'X_FIRST' in atr.keys():
            inps.dem_file = ['demGeo_tight.h5', 'demGeo.h5']
        else:
            inps.dem_file = ['demRadar.h5']
    try:
        inps.dem_file = ut.get_file_list(inps.dem_file)[0]
    except:
        inps.dem_file = None
        sys.exit('ERROR: No DEM file found!')

    # default Mask file
    if not inps.mask_file:
        if 'X_FIRST' in atr.keys():
            inps.mask_file = 'geo_maskTempCoh.h5'
        else:
            inps.mask_file = 'maskTempCoh.h5'
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None
            sys.exit('ERROR: No mask file found!')

    ##### Read Mask
    print('reading mask from file: '+inps.mask_file)
    mask = readfile.read(inps.mask_file, datasetName='mask')[0].flatten(1)
    ndx = mask != 0
    msk_num = np.sum(ndx)
    print('total            pixel number: %d' % pix_num)
    print('estimating using pixel number: %d' % msk_num)

    ##### Read DEM
    print('read DEM from file: '+inps.dem_file)
    dem = readfile.read(inps.dem_file, datasetName='height')[0]

    ref_y = int(atr['REF_Y'])
    ref_x = int(atr['REF_X'])
    dem -= dem[ref_y,ref_x]

    print('considering the incidence angle of each pixel ...')
    inc_angle = ut.incidence_angle(atr, dimension=2)
    dem *= 1.0/np.cos(inc_angle*np.pi/180.0)

    ##### Design matrix for elevation v.s. phase
    dem = dem.flatten(1)
    if inps.poly_order == 1:
        A = np.vstack((dem[ndx], np.ones(msk_num))).T
        B = np.vstack((dem,      np.ones(pix_num))).T
    elif inps.poly_order == 2: 
        A = np.vstack((dem[ndx]**2, dem[ndx], np.ones(msk_num))).T
        B = np.vstack((dem**2,      dem,      np.ones(pix_num))).T  
    elif inps.poly_order == 3:
        A = np.vstack((dem[ndx]**3, dem[ndx]**2, dem[ndx], np.ones(msk_num))).T
        B = np.vstack((dem**3,      dem**2,      dem,      np.ones(pix_num))).T
    print('polynomial order: %d' % inps.poly_order)

    A_inv = np.linalg.pinv(A)

    ##### Calculate correlation coefficient
    print('Estimating the tropospheric effect between the differences of the subsequent epochs and DEM')

    h5 = h5py.File(inps.timeseries_file)
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)
    print('number of acquisitions: '+str(date_num))
    try:    ref_date = atr['REF_DATE']
    except: ref_date = date_list[0]

    print('----------------------------------------------------------')
    print('correlation of DEM with each time-series epoch:')
    corr_array = np.zeros(date_num)
    par_dict = {}
    for i in range(date_num):
        date = date_list[i]
        if date == ref_date:
            cc = 0.0
            par = np.zeros(inps.poly_order+1)
        else:
            data = h5[k].get(date)[:].flatten(1)

            C = np.zeros((2, msk_num))
            C[0,:] = dem[ndx]
            C[1,:] = data[ndx]
            cc = np.corrcoef(C)[0,1]

            corr_array[i] = cc
            if inps.threshold and np.abs(cc) < inps.threshold:
                par = np.zeros(inps.poly_order+1)
            else:
                par = np.dot(A_inv, data[ndx])
        print('%s: %.2f' % (date, cc))
        par_dict[date] = par

    average_phase_height_corr = np.nansum(np.abs(corr_array))/(date_num-1)
    print('----------------------------------------------------------')
    print('Average Correlation of DEM with time-series epochs: %.2f' % average_phase_height_corr)

    # Correlation of DEM with Difference of subsequent epochs (Not used for now)
    corr_diff_dict = {}
    par_diff_dict = {}
    for i in range(date_num-1):
        date1 = date_list[i]
        date2 = date_list[i+1]
        date12 = date1+'-'+date2

        data1 = h5[k].get(date1)[:].flatten(1)
        data2 = h5[k].get(date2)[:].flatten(1)
        data_diff = data2-data1

        C_diff = np.zeros((2, msk_num))
        C_diff[0,:] = dem[ndx]
        C_diff[1,:] = data_diff[ndx]
        cc_diff = np.corrcoef(C_diff)[0,1]

        corr_diff_dict[date12] = cc_diff
        par = np.dot(A_inv, data_diff[ndx])
        par_diff_dict[date12] = par


    ##### Correct and write time-series file
    print('----------------------------------------------------------')
    print('removing the stratified tropospheric delay from each epoch')
    print('writing >>> '+inps.outfile)
    h5out = h5py.File(inps.outfile,'w')
    group = h5out.create_group(k)

    prog_bar = ptime.progressBar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        data = h5[k].get(date)[:]

        if date != ref_date:
            par = par_dict[date]
            trop_delay = np.reshape(np.dot(B, par), [width, length]).T
            trop_delay -= trop_delay[ref_y, ref_x]
            data -= trop_delay

        dset = group.create_dataset(date, data=data)
        prog_bar.update(i+1, suffix=date)

    for key,value in iter(atr.items()):
        group.attrs[key] = value

    prog_bar.close()
    h5out.close()
    h5.close()

    print('Done.')
    return inps.outfile


############################################################################
if __name__ == '__main__':
    main()


