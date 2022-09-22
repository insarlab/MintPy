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


############################################################
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'


############################################################
def usage():
    print('''usage: baseline_error.py [-h] timeseries_file [mask_file] [exclude_date]

Estimating the errors in baseline components (bh,bv,dbh,dbv), and correcting the time-series.

arguments:
  timeseries_file : The timeseries in HDF5 format.
  mask_file       : a mask file to mask out points with high deformation or located in noisy areas

reference:
  Gourmelen, N., F. Amelung, and R. Lanari (2010), Interferometric synthetic aperture radar-GPS
  integration: Interseismic strain accumulation across the Hunter Mountain fault in the eastern
  California shear zone, JGR, 115(B9).

example:
  baseline_error.py  timeseries.h5 Mask.h5
    ''')
    return


############################################################
def main(argv):
    try:
        if argv[0] in ['-h', '--help']:
            usage()
            sys.exit(1)
        else:
            File = argv[0]
    except:
        usage()
        sys.exit(1)

    try:
        maskFile = argv[1]
    except:
        pass

    ##################################
    h5file = h5py.File(File)
    dateList = list(h5file['timeseries'].keys())
    ##################################

    # Read Mask File
    # Priority:
    # Input mask file > mintpy.mask.file > existed Modified_Mask.h5 > existed Mask.h5
    try:
        maskFile
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

    sy, sx = np.shape(Mask)
    npixel = sx*sy

    ##################################
    Mask = Mask.flatten(1)
    ndx = Mask != 0
    ##################################
    nt = float(h5file['timeseries'].attrs['LOOK_REF1'])
    ft = float(h5file['timeseries'].attrs['LOOK_REF2'])
    lookangle = np.tile(np.linspace(nt, ft, sx), [sy, 1])
    lookangle = lookangle.flatten(1)*np.pi/180.0
    Fh = -np.sin(lookangle)
    Fv = -np.cos(lookangle)

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

    A = np.zeros([npixel, 4])

    A[:, 0] = Fh
    A[:, 1] = Fh*rs
    A[:, 2] = Fv
    A[:, 3] = Fv*rs

    Bh = []
    Bv = []
    Bhrate = []
    Bvrate = []
    Be = np.zeros([len(dateList), 4])
    try:
        excludedDates = argv[2]
    except:
        excludedDates = []

    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    for i in range(1, len(dateList)):
        if not dateList[i] in excludedDates:
            dset = h5file['timeseries'].get(dateList[i])
            data = dset[0:dset.shape[0], 0:dset.shape[1]]
            L = data.flatten(1)
            Berror = np.dot(np.linalg.pinv(A[ndx]), L[ndx])
            Bh.append(Berror[0])
            Bhrate.append(Berror[1])
            Bv.append(Berror[2])
            Bvrate.append(Berror[3])
            Be[i, :] = Berror
        else:
            print(str(dateList[i]) +
                  ' is not considered for Baseline Error estimation')

    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('baseline error           mean                          std')
    print('       bh     :  ' + str(np.mean(Bh)) + '     ,  '+str(np.std(Bh)))
    print('     bh rate  :  ' + str(np.mean(Bhrate)) +
          '     ,  '+str(np.std(Bhrate)))
    print('       bv     :  ' + str(np.mean(Bv)) + '     ,  '+str(np.std(Bv)))
    print('     bv rate  :  ' + str(np.mean(Bvrate)) +
          '     ,  '+str(np.std(Bvrate)))
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('bh error of each epoch:')
    print(Bh)
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('bv error of each epoch:')
    print(Bv)
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('Estimating Baseline error from each differences ...')

    Bedif = np.zeros([len(dateList), 4])
    for i in range(1, len(dateList)):
        dset1 = h5file['timeseries'].get(dateList[i-1])
        data1 = dset1[0:dset1.shape[0], 0:dset1.shape[1]]
        dset2 = h5file['timeseries'].get(dateList[i])
        data2 = dset2[0:dset2.shape[0], 0:dset2.shape[1]]
        data = data2-data1
        L = data.flatten(1)
        Berrord = np.dot(np.linalg.pinv(A[ndx]), L[ndx])
        Bedif[i, :] = Berrord

    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    yref = int(h5file['timeseries'].attrs['REF_Y'])
    xref = int(h5file['timeseries'].attrs['REF_X'])

    orbEffect = np.zeros([len(dateList), sy, sx])
    for i in range(1, len(dateList)):
        effect = np.dot(A, Be[i, :])
        effect = np.reshape(effect, [sx, sy]).T
        # orbEffect[i,:,:]=orbEffect[i-1,:,:]+effect
        # orbEffect[i,:,:]=orbEffect[i,:,:]-orbEffect[i,yref,xref]
        orbEffect[i, :, :] = effect - effect[yref, xref]
        del effect

    print('Correctiing the time series ')
    outName = File.replace('.h5', '')+'_baselineCor.h5'
    h5orbCor = h5py.File(outName, 'w')
    group = h5orbCor.create_group('timeseries')
    for i in range(len(dateList)):
        dset1 = h5file['timeseries'].get(dateList[i])
        data = dset1[0:dset1.shape[0], 0:dset1.shape[1]] - orbEffect[i, :, :]
        dset = group.create_dataset(dateList[i], data=data)

    for key, value in h5file['timeseries'].attrs.items():
        group.attrs[key] = value

    try:
        dset1 = h5file['mask'].get('mask')
        group = h5orbCor.create_group('mask')
        dset = group.create_dataset('mask', data=dset1)
    except:
        pass

    h5file.close()
    h5orbCor.close()


if __name__ == '__main__':
    main(sys.argv[1:])
