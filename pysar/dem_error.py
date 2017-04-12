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
#


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
    for i in range(date_num):
        date = date_list[i]
        ut.print_progress(i+1, date_num, prefix='loading:', suffix=date)
        d = h5['timeseries'].get(date)[:]
        timeseries[i][:] = d.flatten('F')
    del d
    h5.close()
    print '-------------------------------------------'

    # Temporal Baseline
    print 'read temporal baseline'
    tbase = ptime.date_list2tbase(date_list)[0]
    tbase = np.array(tbase).reshape(date_num,1)
    # Perpendicular Baseline
    try:
        pbase = [float(i) for i in atr['P_BASELINE_TIMESERIES'].split()]
        pbase = np.array(pbase).reshape(date_num,1)
    except:
        print 'Cannot find P_BASELINE_TIMESERIES from timeseries file.'
        print 'Trying to calculate it from interferograms file'
        if inps.ifgram_file:
            pbase = ut.Baseline_timeseries(inps.ifgram_file)
            pbase = np.array(pbase).reshape(date_num,1)
        else:
            message = 'No interferogram file input!\n'+\
                      'Can not correct for DEM residula without perpendicular base info!'
            raise Exception(message)
    # perpendicular baseline velocity
    pbase_v = np.diff(pbase, axis=0) / np.diff(tbase, axis=0)

    # Incidence angle (look angle in the paper)
    if inps.incidence_angle:
        if os.path.isfile(inps.incidence_angle):
            print 'reading incidence angle from file: '+inps.incidence_angle
            inps.incidence_angle = readfile.read(inps.incidence_angle)[0]
        else:
            try: inps.incidence_angle = float(inps.incidence_angle)
            except: raise ValueError('Can not read input incidence angle: '+str(inps.incidence_angle))
    else:
        print 'calculate incidence angle using attributes of time series file'
        inps.incidence_angle = ut.incidence_angle(atr, dimension=1)
    inps.incidence_angle *= np.pi/180.0
    #inps.incidence_angle = inps.incidence_angle.flatten('F')

    # Range distance
    if inps.range_dis:
        if os.path.isfile(inps.range_dis):
            print 'reading range distance from file: '+inps.range_dis
            inps.range_dis = readfile.read(inps.range_dis)[0]
        else:
            try: inps.range_dis = float(inps.range_dis)
            except: raise ValueError('Can not read input incidence angle: '+str(inps.range_dis))
    else:
        print 'calculate range distance using attributes from time series file'
        inps.range_dis = ut.range_distance(atr, dimension=1)
    #inps.range_dis = inps.range_dis.flatten('F')
    
    # Design matrix - temporal deformation model using tbase
    if inps.phase_velocity:
        print 'using phase velocity history'
        A1 = np.ones((date_num-1, 1))
        A2 = (tbase[1:date_num] + tbase[0:date_num-1]) / 2.0
        A3 = (tbase[1:date_num]**2 + tbase[1:date_num]*tbase[0:date_num-1] + tbase[0:date_num-1]**2) / 6.0
    else:
        print 'using phase history'
        A1 = np.hstack((np.ones((date_num, 1)), tbase))
        A2 = tbase**2 / 2.0
        A3 = tbase**3 / 6.0
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
        print 'inversing using L2-norm minimization (unweighted least squares) pixel by pixel'
        for i in range(length*width):
            row = i%length
            col = i/length
            range_dis = inps.range_dis[row, col]
            inc_angle = inps.incidence_angle[row, col]

            # Design matrix - DEM error using pbase, range distance and incidence angle
            A_delta_z = pbase / (range_dis * np.sin(inc_angle))
            if inps.phase_velocity:
                A_delta_z_v = pbase_v / (range_dis * np.sin(inc_angle))
                A = np.hstack((A_delta_z_v, A_def))
            else:
                A = np.hstack((A_delta_z, A_def))

            # L-2 norm inversion
            A_inv = np.linalg.pinv(A)

            # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
            ts_dis = timeseries[:,i]
            if inps.phase_velocity:
                ts_dis = np.diff(ts_dis, axis=0) / np.diff(tbase, axis=0)
            X = np.dot(A_inv, ts_dis)

            # Update DEM error / timeseries matrix
            delta_z = X[0]
            delta_z_mat[row, col] = delta_z
            if inps.update_timeseries:
                timeseries[:,i] -= np.dot(A_delta_z, delta_z).flatten()
            ut.print_progress(i+1, length*width)


    elif inps.incidence_angle.ndim == 1 and inps.range_dis.ndim == 1:
        print 'inversing using L2-norm minimization (unweighted least squares) column by column'
        for i in range(width):
            range_dis = inps.range_dis[i]
            inc_angle = inps.incidence_angle[i]

            # Design matrix - DEM error using pbase, range distance and incidence angle
            A_delta_z = pbase / (range_dis * np.sin(inc_angle))
            if inps.phase_velocity:
                A_delta_z_v = pbase_v / (range_dis * np.sin(inc_angle))
                A = np.hstack((A_delta_z_v, A_def))
            else:
                A = np.hstack((A_delta_z, A_def))

            # L-2 norm inversion
            A_inv = np.linalg.pinv(A)

            # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
            ts_dis = timeseries[:,i*length:(i+1)*length]
            if inps.phase_velocity:
                ts_dis = np.diff(ts_dis, axis=0) / np.diff(tbase, axis=0)
            X = np.dot(A_inv, ts_dis)

            # Update DEM error / timeseries matrix
            delta_z = X[0].reshape((1,length))
            delta_z_mat[:, i] = delta_z
            if inps.update_timeseries:
                timeseries[:,i*length:(i+1)*length] -= np.dot(A_delta_z, delta_z)
            ut.print_progress(i+1, width)


    elif inps.incidence_angle.ndim == 0 and inps.range_dis.ndim == 0:
        print 'inversing using L2-norm minimization (unweighted least squares) for the whole area'
        
        # Design matrix - DEM error using pbase, range distance and incidence angle
        A_delta_z = pbase / (inps.range_dis * np.sin(inps.incidence_angle))
        if inps.phase_velocity:
            A_delta_z_v = pbase_v / (inps.range_dis * np.sin(inps.incidence_angle))
            A = np.hstack((A_delta_z_v, A_def))
        else:
            A = np.hstack((A_delta_z, A_def))

        # L-2 norm inversion
        A_inv = np.linalg.pinv(A)

        # Get unknown parameters X = [delta_z, vel, acc, delta_acc, ...]
        if inps.phase_velocity:
            timeseries = np.diff(timeseries, axis=0) / np.diff(tbase, axis=0)
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

    
    #C1_v = pbase_v / (center_range * np.sin(center_look_angle))
    #C1   = pbase   / (center_range * np.sin(center_look_angle))
    #timeseries_v = (timeseries[1:date_num,:] - timeseries[0:date_num-1,:]) / (tbase[1:date_num] - tbase[0:date_num-1])    
    ###### Inversion column by column
    #print 'inversing using L2-norm minimization (unweighted least squares)...'
    #dz = np.zeros([1,length*width])
    #
    #for i in range(width):
    #    ## Design Matrix Inversion
    #    C1_v = pbase_v / (range_x[i] * np.sin(look_angle_x[i]))
    #    C1   = pbase   / (range_x[i] * np.sin(look_angle_x[i]))
    #    if inps.phase_velocity:  C = np.hstack((M,C1_v))
    #    else:                    C = np.hstack((M,C1))
    #
    #    #print '    rank of the design matrix : '+str(np.linalg.matrix_rank(C))
    #    #if np.linalg.matrix_rank(C) == 4:  print '    design matrix has full rank'
    #    Cinv = np.linalg.pinv(C)
    #
    #    ## (Phase) Velocity History
    #    ts_x  = timeseries[:,i*length:(i+1)*length]
    #    ts_xv = (ts_x[1:date_num,:] - ts_x[0:date_num-1,:]) / (tbase[1:date_num] - tbase[0:date_num-1])
    #
    #    ## DEM error
    #    if inps.phase_velocity:    par  = np.dot(Cinv,ts_xv)
    #    else:                      par  = np.dot(Cinv,ts_x)
    #    dz_x = par[3].reshape((1,length))
    #
    #    ## Update DEM error matrix and timeseries matrix
    #    dz[0][i*length:(i+1)*length]         = dz_x
    #    timeseries[:,i*length:(i+1)*length] -= np.dot(C1,dz_x)
    #
    #    ut.print_progress(i+1,width)
    #
    ##dz[0][:] = par[3][:]
    #dz = np.reshape(dz,[length,width],order='F')


    ########## Output
    # DEM error file
    if 'Y_FIRST' in atr.keys():
        dem_error_file = 'demGeo_error.h5'
    else:
        dem_error_file = 'demRadar_error.h5'
    #if inps.phase_velocity:  suffix = '_pha_poly'+str(inps.poly_order)
    #else:                    suffix = '_vel_poly'+str(inps.poly_order)
    #dem_error_file = os.path.splitext(dem_error_file)[0]+suffix+os.path.splitext(dem_error_file)[1]
    print '--------------------------------------'
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
        print '--------------------------------------'
        print 'writing >>> '+inps.outfile
        print 'number of dates: '+str(len(date_list))
        h5out = h5py.File(inps.outfile,'w')
        group = h5out.create_group('timeseries')
        for i in range(date_num):
            date = date_list[i]
            ut.print_progress(i+1, date_num, prefix='writing:', suffix=date)
            d = np.reshape(timeseries[i][:], [length,width], order='F')
            dset = group.create_dataset(date, data=d, compression='gzip')
        for key,value in atr.iteritems():
            group.attrs[key] = value
        h5out.close()
    
    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  



