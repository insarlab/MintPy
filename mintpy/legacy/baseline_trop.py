#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import os
import sys

import h5py
import matplotlib
import numpy as np

from mintpy.utils import readfile


####################################################################################
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'


####################################################################################
def usage():
    print('''usage: baseline_trop.py timeseries_file dem_file poly_order {range_and_azimuth|range|azimuth} mask_file

Simultaneously correcting the baseline error and stratified tropospheric delay correlated with DEM.

reference:
  Jo, M.-J., J.-S. Won, S.-W. Kim, and H.-S. Jung (2010), A time-series SAR observation of surface
  deformation at the southern end of the San Andreas Fault Zone, Geos J., 14(3), 277-287..

example:
  baseline_trop.py  timeseries.h5 radar.hgt 1
  baseline_trop.py  timeseries.h5 radar.hgt 1 range
  baseline_trop.py  timeseries.h5 radar.hgt 1 range_and_azimuth mask.h5
    ''')
    return


####################################################################################
def main(argv):

    try:
        File = argv[0]
        demFile = argv[1]
        p = int(argv[2])
    except:
        usage()
        sys.exit(1)

    try:
        baseline_error = argv[3]
    except:
        baseline_error = 'range_and_azimuth'
    print(baseline_error)
    ##################################
    h5file = h5py.File(File)
    dateList = list(h5file['timeseries'].keys())
    ##################################

    try:
        maskFile = argv[4]
    except:
        if os.path.isfile('Modified_Mask.h5'):
            maskFile = 'Modified_Mask.h5'
        elif os.path.isfile('Mask.h5'):
            maskFile = 'Mask.h5'
        else:
            print('No mask found!')
            sys.exit(1)
    try:
        Mask, Matr = readfile.read(maskFile, datasetName='mask')
        print('mask: '+maskFile)
    except:
        print('Can not open mask file: '+maskFile)
        sys.exit(1)

    ##################################
    Mask = Mask.flatten(1)
    ndx = Mask != 0
    ##################################
    # h5file = h5py.File(File)
    # dateList = h5file['timeseries'].keys()
    ##################################
    nt = float(h5file['timeseries'].attrs['LOOK_REF1'])
    ft = float(h5file['timeseries'].attrs['LOOK_REF2'])
    sy, sx = np.shape(Mask)
    npixel = sx*sy
    lookangle = np.tile(np.linspace(nt, ft, sx), [sy, 1])
    lookangle = lookangle.flatten(1)*np.pi/180.0
    Fh = -np.sin(lookangle)
    Fv = -np.cos(lookangle)

    print('Looking for azimuth pixel size')
    try:
        daz = float(h5file['timeseries'].attrs['AZIMUTH_PIXEL_SIZE'])
    except:
        print('''
        ERROR!
        The attribute AZIMUTH_PIXEL_SIZE was not found!
        Possible cause of error: Geo coordinate.
        This function works only in radar coordinate system.
        ''')
        sys.exit(1)

    lines = np.tile(np.arange(0, sy, 1), [1, sx])
    lines = lines.flatten(1)
    rs = lines*daz

    if baseline_error == 'range_and_azimuth':
        A = np.zeros([npixel, 4])

        A[:, 0] = Fh
        A[:, 1] = Fh*rs
        A[:, 2] = Fv
        A[:, 3] = Fv*rs
        num_base_par = 4
    elif baseline_error == 'range':
        A = np.zeros([npixel, 2])

        A[:, 0] = Fh
        A[:, 1] = Fv
        num_base_par = 2

    ###########################################
    yref = int(h5file['timeseries'].attrs['REF_Y'])
    xref = int(h5file['timeseries'].attrs['REF_X'])
    ###########################################
    if os.path.basename(demFile).split('.')[1] == 'hgt':
        amp, dem, demRsc = readfile.read_float32(demFile)
    elif os.path.basename(demFile).split('.')[1] == 'dem':
        dem, demRsc = readfile.read_real_int16(demFile)

    dem = dem-dem[yref][xref]
    dem = dem.flatten(1)
    ###################################################
    if p == 1:
        # A=np.vstack((dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem, np.ones(len(dem)))).T
    elif p == 2:
        # A=np.vstack((dem[ndx]**2,dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem**2, dem, np.ones(len(dem)))).T
    elif p == 3:
        #  A = np.vstack((dem[ndx]**3,dem[ndx]**2,dem[ndx],np.ones(len(dem[ndx])))).T
        B = np.vstack((dem**3, dem**2, dem, np.ones(len(dem)))).T
    print(np.shape(A))

    ###################################################

    Bh = []
    Bv = []
    Bhrate = []
    Bvrate = []
    Be = np.zeros([len(dateList), num_base_par+p+1])
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    for i in range(1, len(dateList)):
        dset = h5file['timeseries'].get(dateList[i])
        data = dset[0:dset.shape[0], 0:dset.shape[1]]
        L = data.flatten(1)
        M = np.hstack((A, B))
        Berror = np.dot(np.linalg.pinv(M[ndx]), L[ndx])
        Bh.append(Berror[0])
        Bhrate.append(Berror[1])
        Bv.append(Berror[2])
        Bvrate.append(Berror[3])
        Be[i, :] = Berror
        print(Berror)
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('baseline error           mean                          std')
    print('       bh     :  ' + str(np.mean(Bh)) + '     ,  '+str(np.std(Bh)))
    print('     bh rate  :  ' + str(np.mean(Bhrate)) +
          '     ,  '+str(np.std(Bhrate)))
    print('       bv     :  ' + str(np.mean(Bv)) + '     ,  '+str(np.std(Bv)))
    print('     bv rate  :  ' + str(np.mean(Bvrate)) +
          '     ,  '+str(np.std(Bvrate)))
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    orbEffect = np.zeros([len(dateList), sy, sx])
    for i in range(1, len(dateList)):
        effect = np.dot(M, Be[i, :])
        effect = np.reshape(effect, [sx, sy]).T
        # orbEffect[i,:,:]=orbEffect[i-1,:,:]+effect
        # orbEffect[i,:,:]=orbEffect[i,:,:]-orbEffect[i,yref,xref]
        orbEffect[i, :, :] = effect - effect[yref, xref]
        del effect

    print('Correctiing the time series ')
    outName = File.replace('.h5', '')+'_baseTropCor.h5'
    h5orbCor = h5py.File(outName, 'w')
    group = h5orbCor.create_group('timeseries')
    for i in range(len(dateList)):
        dset1 = h5file['timeseries'].get(dateList[i])
        data = dset1[0:dset1.shape[0], 0:dset1.shape[1]] - orbEffect[i, :, :]
        dset = group.create_dataset(dateList[i], data=data)

    for key, value in h5file['timeseries'].attrs.items():
        group.attrs[key] = value

    dset1 = h5file['mask'].get('mask')
    group = h5orbCor.create_group('mask')
    dset = group.create_dataset('mask', data=dset1)

    h5file.close()
    h5orbCor.close()
    return


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
