#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Reference:
# Fattahi, H. and F. Amelung, (2013), DEM-error correction in 
# InSAR time-series analysis, IEEE TGRS, vol. no.99,
# doi: 10.1109/TGRS.2012.2227761.
#
# Yunjun, Jun 2016: Add phase velocity approach from the paper.
#                   Use different range and look angle for each column
# Yunjun, Apr 2017: use variable P_BASELINE(_TOP/BOTTOM)_TIMESERIES
#                   support geocoded file


import os
import sys
import argparse

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar._writefile as writefile


######################################
EXAMPLE='''example:
  dem_error.py  timeseries_ECMWF.h5
  dem_error.py  timeseries_ECMWF.h5  --phase-velocity
  dem_error.py  timeseries_ECMWF.h5  -d dem_radar.h5
  dem_error.py  geo_timeseries.h5    -i geo_incidence_angle.h5  -r geo_range.h5
'''

REFERENCE='''reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
  IEEE TGRS, 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='DEM Error Correction.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file to be corrrected')
    parser.add_argument('-o','--outfile', help='Output file name for corrected time series')
    
    parser.add_argument('-i', dest='incidence_angle', help='Incidence angle value/file in degrees')
    parser.add_argument('-r','--range-distance', dest='range_dis', help='Range distance value/file')
    #parser.add_argument('-d','--dem', dest='dem_file', help='DEM file to be updated')
    parser.add_argument('--interferograms', dest='ifgram_file',\
                        help='Interferograms file to get perpendicular baseline time series\n'+\
                             'for old PySAR time series file only, will be removed in the future.')
    parser.add_argument('--phase-velocity', dest='phase_velocity', action='store_true',\
                        help='Use phase velocity instead of phase for inversion constrain.')
    parser.add_argument('--no-update-timeseries', dest='update_timeseries', action='store_false',\
                        help='Do not update timeseries; if specified, only DEM error will be calculated.')
    parser.add_argument('--poly-order', dest='poly_order', type=int, default=2, choices=[1,2,3],\
                        help='polynomial order number of temporal deformation model, default = 2')

    inps = parser.parse_args()
    return inps


######################################
def main(argv):
    inps = cmdLineParse()
    suffix = '_demErr'
    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.timeseries_file)[0]+suffix+os.path.splitext(inps.timeseries_file)[1]

    # Read Time Series
    print "loading time series: " + inps.timeseries_file
    atr = readfile.read_attribute(inps.timeseries_file)
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    h5 = h5py.File(inps.timeseries_file)
    date_list = sorted(h5['timeseries'].keys())
    date_num = len(date_list)
    print 'number of acquisitions: '+str(date_num)

    timeseries = np.zeros((len(date_list),length*width),np.float32)
    prog_bar = ut.progress_bar(maxValue=date_num, prefix='loading:')
    for i in range(date_num):
        date = date_list[i]
        d = h5['timeseries'].get(date)[:]
        timeseries[i][:] = d.flatten('F')
        prog_bar.update(i+1, suffix=date)
    del d
    h5.close()
    prog_bar.close()

    # Perpendicular Baseline
    print 'read perpendicular baseline'
    try:
        inps.pbase = ut.perp_baseline_timeseries(atr, dimension=1)
        if inps.pbase.shape[1] > 1:
            print '\tconsider P_BASELINE variation in azimuth direction'
        else:
            pbase = inps.pbase
    except:
        print '\tCannot find P_BASELINE_TIMESERIES from timeseries file.'
        print '\tTrying to calculate it from interferograms file'
        if inps.ifgram_file:
            inps.pbase = np.array(ut.perp_baseline_ifgram2timeseries(inps.ifgram_file)[0]).reshape(date_num,1)
        else:
            message = 'No interferogram file input!\n'+\
                      'Can not correct for DEM residula without perpendicular base info!'
            raise Exception(message)

    # Temporal Baseline
    print 'read temporal baseline'
    inps.tbase = np.array(ptime.date_list2tbase(date_list)[0]).reshape(date_num,1)

    # Incidence angle (look angle in the paper)
    if inps.incidence_angle:
        if os.path.isfile(inps.incidence_angle):
            print 'reading incidence angle from file: '+inps.incidence_angle
            inps.incidence_angle = readfile.read(inps.incidence_angle)[0]
        else:
            try:
                inps.incidence_angle = np.array(float(inps.incidence_angle))
                print 'use input incidence angle : '+str(inps.incidence_angle)
            except:
                raise ValueError('Can not read input incidence angle: '+str(inps.incidence_angle))
    else:
        print 'calculate incidence angle using attributes of time series file'
        if inps.pbase.shape[1] > 1:
            inps.incidence_angle = ut.incidence_angle(atr, dimension=2)
        else:
            inps.incidence_angle = ut.incidence_angle(atr, dimension=1)
    inps.incidence_angle *= np.pi/180.0

    # Range distance
    if inps.range_dis:
        if os.path.isfile(inps.range_dis):
            print 'reading range distance from file: '+inps.range_dis
            inps.range_dis = readfile.read(inps.range_dis)[0]
        else:
            try:
                inps.range_dis = np.array(float(inps.range_dis))
                print 'use input range distance : '+str(inps.range_dis)
            except:
                raise ValueError('Can not read input incidence angle: '+str(inps.range_dis))
    else:
        print 'calculate range distance using attributes from time series file'
        if inps.pbase.shape[1] > 1:
            inps.range_dis = ut.range_distance(atr, dimension=2)
        else:
            inps.range_dis = ut.range_distance(atr, dimension=1)


    # Design matrix - temporal deformation model using tbase
    if inps.phase_velocity:
        print 'using phase velocity history'
        A1 = np.ones((date_num-1, 1))
        A2 = (inps.tbase[1:date_num] + inps.tbase[0:date_num-1]) / 2.0
        A3 = (inps.tbase[1:date_num]**3 - inps.tbase[0:date_num-1]**3) / np.diff(inps.tbase, axis=0) / 6.0
        #A3 = (inps.tbase[1:date_num]**2 + inps.tbase[1:date_num]*inps.tbase[0:date_num-1] +\
        #      inps.tbase[0:date_num-1]**2) / 6.0
    else:
        print 'using phase history'
        A1 = np.hstack((np.ones((date_num, 1)), inps.tbase))
        A2 = inps.tbase**2 / 2.0
        A3 = inps.tbase**3 / 6.0
    # Polynomial order of model
    print "temporal deformation model's polynomial order = "+str(inps.poly_order)
    if   inps.poly_order == 1:  A_def = A1
    elif inps.poly_order == 2:  A_def = np.hstack((A1,A2))
    elif inps.poly_order == 3:  A_def = np.hstack((A1,A2,A3))
    # Heresh's original code for phase history approach
    #A_def = np.hstack((A2,A1,np.ones((date_num,1))))


    ##---------------------------------------- Loop for L2-norm inversion  -----------------------------------##
    delta_z_mat = np.zeros([length, width])
    #delta_a_mat = np.zeros([length, width])
    if inps.incidence_angle.ndim == 2 and inps.range_dis.ndim == 2:
        print 'inversing using L2-norm minimization (unweighted least squares)'\
              ' pixel by pixel: %d loops in total' % (length*width)
        prog_bar = ut.progress_bar(maxValue=length*width, prefix='calculating:')
        for i in range(length*width):
            row = i%length
            col = i/length
            range_dis = inps.range_dis[row, col]
            inc_angle = inps.incidence_angle[row, col]
            # Consider P_BASELINE variation within one interferogram
            if inps.pbase.shape[1] > 1:
                pbase = inps.pbase[:,row].reshape(date_num, 1)

            # Design matrix - DEM error using pbase, range distance and incidence angle
            A_delta_z = pbase / (range_dis * np.sin(inc_angle))
            if inps.phase_velocity:
                pbase_v = np.diff(pbase, axis=0) / np.diff(inps.tbase, axis=0)
                A_delta_z_v = pbase_v / (range_dis * np.sin(inc_angle))
                A = np.hstack((A_delta_z_v, A_def))
            else:
                A = np.hstack((A_delta_z, A_def))

            # L-2 norm inversion
            A_inv = np.linalg.pinv(A)

            # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
            ts_dis = timeseries[:,i]
            if inps.phase_velocity:
                ts_dis = np.diff(ts_dis, axis=0) / np.diff(inps.tbase, axis=0)
            X = np.dot(A_inv, ts_dis)

            # Update DEM error / timeseries matrix
            delta_z = X[0]
            delta_z_mat[row, col] = delta_z
            if inps.update_timeseries:
                timeseries[:,i] -= np.dot(A_delta_z, delta_z).flatten()
            prog_bar.update(i+1, every=length*width/100)
        prog_bar.close()


    elif inps.incidence_angle.ndim == 1 and inps.range_dis.ndim == 1:
        print 'inversing using L2-norm minimization (unweighted least squares)'\
              ' column by column: %d loops in total' % (width)
        prog_bar = ut.progress_bar(maxValue=width, prefix='calculating:')
        for i in range(width):
            range_dis = inps.range_dis[i]
            inc_angle = inps.incidence_angle[i]

            # Design matrix - DEM error using pbase, range distance and incidence angle
            A_delta_z = pbase / (range_dis * np.sin(inc_angle))
            if inps.phase_velocity:
                pbase_v = np.diff(pbase, axis=0) / np.diff(inps.tbase, axis=0)
                A_delta_z_v = pbase_v / (range_dis * np.sin(inc_angle))
                A = np.hstack((A_delta_z_v, A_def))
            else:
                A = np.hstack((A_delta_z, A_def))

            # L-2 norm inversion
            A_inv = np.linalg.pinv(A)

            # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
            ts_dis = timeseries[:,i*length:(i+1)*length]
            if inps.phase_velocity:
                ts_dis = np.diff(ts_dis, axis=0) / np.diff(inps.tbase, axis=0)
            X = np.dot(A_inv, ts_dis)

            # Update DEM error / timeseries matrix
            delta_z = X[0].reshape((1,length))
            delta_z_mat[:, i] = delta_z
            if inps.update_timeseries:
                timeseries[:,i*length:(i+1)*length] -= np.dot(A_delta_z, delta_z)
            prog_bar.update(i+1, every=width/100)
        prog_bar.close()


    elif inps.incidence_angle.ndim == 0 and inps.range_dis.ndim == 0:
        print 'inversing using L2-norm minimization (unweighted least squares) for the whole area'
        
        # Design matrix - DEM error using pbase, range distance and incidence angle
        A_delta_z = pbase / (inps.range_dis * np.sin(inps.incidence_angle))
        if inps.phase_velocity:
            pbase_v = np.diff(pbase, axis=0) / np.diff(inps.tbase, axis=0)
            A_delta_z_v = pbase_v / (inps.range_dis * np.sin(inps.incidence_angle))
            A = np.hstack((A_delta_z_v, A_def))
        else:
            A = np.hstack((A_delta_z, A_def))

        # L-2 norm inversion
        A_inv = np.linalg.pinv(A)

        # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
        if inps.phase_velocity:
            timeseries = np.diff(timeseries, axis=0) / np.diff(inps.tbase, axis=0)
        X = np.dot(A_inv, timeseries)

        # Update DEM error / timeseries matrix
        delta_z_mat = X[0].reshape((1, length*width))
        if inps.update_timeseries:
            timeseries -= np.dot(A_delta_z, delta_z_mat)
        delta_z_mat = np.reshape(delta_z_mat, [length, width], order='F')

    else:
        print 'ERROR: Script only support same dimension for both incidence angle and range distance matrix.'
        print 'dimension of incidence angle: '+str(inps.incidence_angle.ndim)
        print 'dimension of range distance: '+str(inps.range_dis.ndim)
        sys.exit(1)


    ##------------------------------------------------ Output  --------------------------------------------##
    # DEM error file
    if 'Y_FIRST' in atr.keys():
        dem_error_file = 'demGeo_error.h5'
    else:
        dem_error_file = 'demRadar_error.h5'
    #if inps.phase_velocity:  suffix = '_pha_poly'+str(inps.poly_order)
    #else:                    suffix = '_vel_poly'+str(inps.poly_order)
    #dem_error_file = os.path.splitext(dem_error_file)[0]+suffix+os.path.splitext(dem_error_file)[1]
    print 'writing >>> '+dem_error_file
    atr_dem_error = atr.copy()
    atr_dem_error['FILE_TYPE'] = 'dem'
    atr_dem_error['UNIT'] = 'm'
    writefile.write(delta_z_mat, atr_dem_error, dem_error_file)
    
    ## Corrected DEM file
    #if inps.dem_file:
    #    inps.dem_outfile = os.path.splitext(inps.dem_file)[0]+suffix+os.path.splitext(inps.dem_file)[1]
    #    print '--------------------------------------'
    #    print 'writing >>> '+inps.dem_outfile
    #    dem, atr_dem = readfile.read(inps.dem_file)
    #    writefile.write(dem+delta_z_mat, atr_dem, inps.dem_outfile)
    
    #outfile = 'delta_acc.h5'
    #print 'writing >>> '+outfile
    #atr_dem_error = atr.copy()
    #atr_dem_error['FILE_TYPE'] = 'velocity'
    #atr_dem_error['UNIT'] = 'm/s'
    #writefile.write(delta_a_mat, atr_dem_error, outfile)
    #print '**************************************'

    # Corrected Time Series
    if inps.update_timeseries:
        print 'writing >>> '+inps.outfile
        print 'number of dates: '+str(len(date_list))
        h5out = h5py.File(inps.outfile,'w')
        group = h5out.create_group('timeseries')
        prog_bar = ut.progress_bar(maxValue=date_num, prefix='writing:')
        for i in range(date_num):
            date = date_list[i]
            d = np.reshape(timeseries[i][:], [length,width], order='F')
            dset = group.create_dataset(date, data=d, compression='gzip')
            prog_bar.update(i+1, suffix=date)
        prog_bar.close()
        for key,value in atr.iteritems():
            group.attrs[key] = value
        h5out.close()
    
    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  



