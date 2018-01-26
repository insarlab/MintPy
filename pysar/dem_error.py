#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse

import h5py
import numpy as np
from scipy.special import gamma

import _datetime as ptime
import _pysar_utilities as ut
import _readfile as readfile
import _writefile as writefile


def topographic_residual_inversion(ts0, A0, inps):
    '''
    Inputs:
        ts0  - 2D np.array in size of (date_num, pixel_num), original time series displacement
        A0   - 2D np.array in size of (date_num, model_num), design matrix in [A_deltaZ, A_def]
        inps - Namespace with the following settings:
               tbase     - 2D np.array in size of (date_num, 1), temporal baseline
               date_flag - 1D np.array in bool data type, mark the date used in the estimation
               phase_velocity - bool, use phase history or phase velocity for minimization
    Outputs:
        deltaZ - 2D np.array in size of (1,        pixel_num), estimated DEM residual
        tsCor  - 2D np.array in size of (date_num, pixel_num), corrected timeseries = tsOrig - topoRes
        tsRes  - 2D np.array in size of (date_num, pixel_num), resudal   timeseries = tsOrig - topoRes - defModel
        stepEst- 2D np.array in size of (step_num, pixel_num), estimated step deformation
    Example:
        deltaZ, tsCor, tsRes = topographic_residual_inversion(ts, A, inps)
    '''
    if len(ts0.shape) == 1:
        ts0 = ts0.reshape(-1,1)
    date_num, pixel_num = ts0.shape

    ##Prepare Design matrix A and observations ts for inversion
    A = A0[inps.date_flag,:]
    ts = ts0[inps.date_flag,:]
    if inps.phase_velocity:
        ts = np.diff(ts, axis=0) / np.diff(inps.tbase, axis=0)
        A = np.diff(A, axis=0) / np.diff(inps.tbase, axis=0)

    ##Inverse using L-2 norm to get unknown parameters X = [deltaZ, vel, acc, deltaAcc, ...]
    X = np.linalg.pinv(A).dot(ts)   #equivalent to X = np.linalg.inv(A.T.dot(A)).dot(A.T).dot(ts)

    ##Prepare Outputs
    deltaZ = X[0,:]
    tsCor = ts0 - np.dot(A0[:,0].reshape(-1,1), deltaZ.reshape(1,-1))
    tsRes = ts0 - np.dot(A0, X)

    stepEst = None
    if inps.step_num > 0:
        s0 = inps.poly_order+2
        s1 = inps.poly_order+2+inps.step_num
        stepEst = X[s0:s1,:].reshape(inps.step_num,-1)

    return deltaZ, tsCor, tsRes, stepEst


def read_template2inps(template_file, inps=None):
    '''Read input template file into inps.ex_date'''
    if not inps:
        inps = cmdLineParse()
    template = readfile.read_template(template_file)
    key_list = list(template.keys())

    # Read template option
    prefix = 'pysar.topographicResidual.'

    key = prefix+'polyOrder'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.poly_order = 2
        else:
            inps.poly_order = int(value)

    key = prefix+'excludeDate'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            value = value.replace(',',' ').split()
            value = ptime.yyyymmdd(value)
            inps.ex_date += value

    key = prefix+'stepFuncDate'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            value = value.replace(',',' ').split()
            value = ptime.yyyymmdd(value)
            inps.step_date += value

    return inps


def check_exclude_date(exDateIn, dateList):
    '''Read exclude dates info
    Inputs:
        exDateIn - list of string, date in YYMMDD or YYYYMMDD format,
                   or text file with date in it
        dateList - list of string, date in YYYYMMDD format
    Output:
        exDateOut - list of string, date in YYYYMMDD format
    '''
    if not exDateIn:
        return []

    exDateOut = []
    for exDate in exDateIn:
        if os.path.isfile(exDate):
            exDate = ptime.read_date_list(exDate)
        else:
            exDate = [ptime.yyyymmdd(exDate)]
        exDateOut += exDate
    exDateOut = sorted(list(set(exDateOut).intersection(dateList)))
    print('Exclude date for DEM error estimation:')
    print exDateOut
    return exDateOut


######################################
TEMPLATE='''
## 8. Topographic (DEM) Residual Correction (Fattahi and Amelung, 2013, IEEE-TGRS)
## Specify stepFuncDate option if you know there are sudden displacement jump in your area,
## i.e. volcanic eruption, or earthquake, and check timeseriesStepModel.h5 afterward for their estimation.
pysar.topographicResidual              = auto  #[yes / no], auto for yes
pysar.topographicResidual.polyOrder    = auto  #[1-inf], auto for 2, polynomial order of temporal deformation model
pysar.topographicResidual.stepFuncDate = auto  #[20080529,20100611 / no], auto for no, date of step jump
pysar.topographicResidual.excludeDate  = auto  #[20070321 / txtFile / no], auto for no, date exlcuded for error estimation
'''

EXAMPLE='''example:
  dem_error.py  timeseries_ECMWF.h5
  dem_error.py  timeseries_ECMWF.h5  --phase-velocity
  dem_error.py  timeseries_ECMWF.h5  -d dem_radar.h5
  dem_error.py  geo_timeseries.h5    -i geo_incidence_angle.h5  -r geo_range.h5

  dem_error.py  timeseries_ECMWF.h5 --template pysarApp_template.txt
'''

REFERENCE='''reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
  IEEE TGRS, 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Topographic (DEM) Residual Correction',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file to be corrrected')
    parser.add_argument('-o','--outfile', help='Output file name for corrected time series')
    parser.add_argument('--exclude','--ex', dest='ex_date', nargs='*', default=[],\
                        help='Exclude date(s) for DEM error estimation.\n'+\
                             'All dates will be corrected for DEM residual phase still.')
    parser.add_argument('--template','-t', dest='template_file',\
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('--step-date', dest='step_date', nargs='*', default=[],\
                        help='Date of step jump for temporal deformation model, i.e. date of earthquake/volcanic eruption')

    parser.add_argument('-i', dest='inc_angle_file', default=['geometry*.h5','*incidenceAngle.h5'],\
                        help='Incidence angle file in degrees')
    parser.add_argument('-r','--range-distance', dest='range_dist_file', default=['geometry*.h5','*angeDistance.h5'],\
                        help='Range distance file with value in meter')
    parser.add_argument('--phase-velocity', dest='phase_velocity', action='store_true',\
                        help='Use phase velocity instead of phase for inversion constrain.')
    #parser.add_argument('--no-update-timeseries', dest='update_timeseries', action='store_false',\
    #                    help='Do not update timeseries; if specified, only DEM error will be calculated.')
    parser.add_argument('--poly-order','-p', dest='poly_order', type=int, default=2,\
                        help='polynomial order number of temporal deformation model, default = 2')

    inps = parser.parse_args()
    if inps.poly_order < 1:
        raise argparse.ArgumentTypeError("Minimum polynomial order is 1")
    return inps


######################################
def main(argv):
    inps = cmdLineParse()
    suffix = '_demErr'
    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.timeseries_file)[0]+suffix+os.path.splitext(inps.timeseries_file)[1]
    if inps.template_file:
        print('read option from template file: '+inps.template_file)
        inps = read_template2inps(inps.template_file, inps)

    atr = readfile.read_attribute(inps.timeseries_file)
    coordType = 'radar'
    if 'Y_FIRST' in atr.keys():
        coordType = 'geo'

    # 1. Incidence angle
    try:
        inps.inc_angle_file = ut.get_file_list(inps.inc_angle_file, coord=coordType)[0]
    except ValueError:
        print 'No incidence angle file found!\nRun incidence_angle.py to generate it.'
    print 'read incidence angle from file: '+str(inps.inc_angle_file)
    inps.inc_angle = readfile.read(inps.inc_angle_file, epoch='incidenceAngle')[0].flatten()
    inps.inc_angle *= np.pi/180.0

    # 2. Slant Range distance
    try:
        inps.range_dist_file = ut.get_file_list(inps.range_dist_file, coord=coordType)[0]
    except ValueError:
        print 'No range distance file found!\nRun range_distance.py to generate it.'
    print 'read slant range distance from file: '+str(inps.range_dist_file)
    inps.range_dist = readfile.read(inps.range_dist_file, epoch='slantRangeDistance')[0].flatten()

    # 3. Perp Baseline - 1D in time, 0D/1D in space (azimuth)
    print 'read perpendicular baseline'
    try:
        inps.pbase = ut.perp_baseline_timeseries(atr, dimension=1)
        if inps.pbase.shape[1] > 1:
            print 'consider perp baseline variance in azimuth direction'
    except valueError:
        print 'No P_BASELINE_TIMESERIES found in timeseries file.\n'+\
              'Can not correct for DEM residula without it!'

    # 4. Time Series - 1D in time, 1D in space (flattened)
    print "read time series file: " + inps.timeseries_file
    h5 = h5py.File(inps.timeseries_file)
    date_list = sorted(h5['timeseries'].keys())
    date_num = len(date_list)

    inps.tbase = np.array(ptime.date_list2tbase(date_list)[0]).reshape(-1,1)

    #Mark dates used in the estimation
    inps.ex_date = check_exclude_date(inps.ex_date, date_list)
    inps.date_flag = np.array([i not in inps.ex_date for i in date_list], dtype=np.bool_)
    if inps.poly_order > np.sum(inps.date_flag):
        raise ValueError("ERROR: input polynomial order=%d is larger than number of acquisition=%d used in estimation!" %\
                         (inps.poly_order, np.sum(inps.date_flag)))

    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])
    pixel_num = length*width
    timeseries = np.zeros((date_num, pixel_num),np.float32)
    for i in range(date_num):
        timeseries[i] = h5['timeseries'].get(date_list[i])[:].flatten()
        sys.stdout.write('\rreading acquisition %3d/%3d ...' % (i+1, date_num))
        sys.stdout.flush()
    h5.close()
    print('')


    ##### Design matrix - temporal deformation model
    print('-------------------------------------------------')
    print('Correct topographic phase residual using Fattahi and Amelung (2013, IEEE-TGRS)')
    msg = 'minimum-norm constrain on: phase'
    if inps.phase_velocity:
        msg += ' velocity'
    print msg

    # Heresh's original code for phase history approach
    #A1 = np.hstack((np.ones((date_num, 1)), inps.tbase))
    #A2 = inps.tbase**2 / 2.0
    #A_def = np.hstack((A2,A1,np.ones((date_num,1))))

    # 1. Polynomial - 2D matrix in size of (date_num, polyOrder+1)
    print "temporal deformation model: polynomial order = "+str(inps.poly_order)
    A_def = np.ones((date_num, 1), np.float32)
    for i in range(inps.poly_order):
        Ai = inps.tbase**(i+1) / gamma(i+2)
        Ai = np.array(Ai, np.float32).reshape(-1,1)
        A_def = np.hstack((A_def, Ai))

    # 2. Step function - 2D matrix in size of (date_num, stepNum)
    if inps.step_date:
        print "temporal deformation model: step functions at "+str(inps.step_date)
        yySteps = ptime.yyyymmdd2years(inps.step_date)
        yyList = np.array(ptime.yyyymmdd2years(date_list)).reshape(-1,1)
        for yyStep in yySteps:
            Ai = yyList > yyStep
            Ai = np.array(Ai, np.float32).reshape(-1,1)
            A_def = np.hstack((A_def, Ai))
    inps.step_num = len(inps.step_date)

    print '-------------------------------------------------'



    ##---------------------------------------- Loop for L2-norm inversion  -----------------------------------##
    ## Output estimated steps 
    print 'ordinal least squares (OLS) inversion using L2-norm minimization'
    timeseriesCor = np.zeros((date_num, pixel_num), dtype=np.float32)
    timeseriesRes = np.zeros((date_num, pixel_num), dtype=np.float32)
    topoRes = np.zeros(pixel_num, dtype=np.float32)
    constC  = np.zeros(pixel_num, dtype=np.float32)
    if inps.step_num > 0:
        stepModel = np.zeros((inps.step_num, pixel_num), dtype=np.float32)

    print 'skip pixels with zero/nan value in geometry files - incidence angle and range distance'
    mask = np.multiply(~np.isnan(inps.inc_angle), ~np.isnan(inps.range_dist))
    mask[inps.inc_angle == 0.] = 0
    mask[inps.range_dist == 0.] = 0
    pixel_num2inv = np.sum(mask)
    pixel_idx2inv = np.where(mask)[0]
    print 'number of pixels in the file: %d' % (pixel_num)
    print 'number of pixels to  inverse: %d' % (pixel_num2inv)

    if inps.pbase.shape[1] == 1:
        pbase = inps.pbase
    prog_bar = ptime.progress_bar(maxValue=pixel_num)
    for i in range(pixel_num2inv):
        prog_bar.update(i+1, every=1000, suffix='%s/%s pixels'%(str(i+1), str(pixel_num2inv)))
        idx = pixel_idx2inv[i]

        r = inps.range_dist[idx]
        inc_angle = inps.inc_angle[idx]
        if inps.pbase.shape[1] > 1:
            pbase = inps.pbase[:, int(idx/width)].reshape(-1,1)
        A_deltaZ = pbase / (r * np.sin(inc_angle))

        A = np.hstack((A_deltaZ, A_def))
        ts = timeseries[:,idx].reshape(date_num,-1)
        deltaZ, tsCor, tsRes, stepEst = topographic_residual_inversion(ts, A, inps)
        topoRes[idx:idx+1] = deltaZ
        timeseriesCor[:,idx:idx+1] = tsCor
        timeseriesRes[:,idx:idx+1] = tsRes
        if inps.step_num > 0:
            stepModel[:,idx:idx+1] = stepEst
    prog_bar.close()


    ##------------------------------------------------ Output  --------------------------------------------##
    # 1. DEM error file
    if 'Y_FIRST' in atr.keys():
        deltaZFile = 'demGeo_error.h5'
    else:
        deltaZFile = 'demRadar_error.h5'
    print 'writing >>> '+deltaZFile
    atrDeltaZ = atr.copy()
    atrDeltaZ['FILE_TYPE'] = 'dem'
    atrDeltaZ['UNIT'] = 'm'
    writefile.write(topoRes.reshape(length, width), atrDeltaZ, deltaZFile)

    # 2. Topo Residual Corrected Time Series
    print 'writing >>> '+inps.outfile
    h5 = h5py.File(inps.outfile,'w')
    group = h5.create_group('timeseries')
    for i in range(date_num):
        sys.stdout.write('\rwriting acquisition %3d/%3d ...' % (i+1, date_num))
        sys.stdout.flush()
        dset = group.create_dataset(date_list[i], data=timeseriesCor[i].reshape(length, width), compression='gzip')
    print ''
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5.close()

    # 3. Inversion residual Time Series
    tsResFile = os.path.join(os.path.dirname(inps.outfile), 'timeseriesResidual.h5')
    print 'writing >>> '+os.path.basename(tsResFile)
    h5 = h5py.File(tsResFile,'w')
    group = h5.create_group('timeseries')
    for i in range(date_num):
        sys.stdout.write('\rwriting acquisition %3d/%3d ...' % (i+1, date_num))
        sys.stdout.flush()
        dset = group.create_dataset(date_list[i], data=timeseriesRes[i].reshape(length, width), compression='gzip')
    print ''

    # Attribute
    for key,value in atr.items():
        group.attrs[key] = value
    h5.close()

    # 4. Step temporal Model estimation
    if inps.step_num > 0:
        stepFile = os.path.join(os.path.dirname(inps.outfile), 'timeseriesStepModel.h5')
        print 'writing >>> '+os.path.basename(stepFile)
        h5 = h5py.File(stepFile,'w')
        group = h5.create_group('timeseries')
        for i in range(inps.step_num):
            sys.stdout.write('\rwriting acquisition %3d/%3d ...' % (i+1, inps.step_num))
            sys.stdout.flush()
            dset = group.create_dataset(inps.step_date[i], data=stepModel[i].reshape(length, width), compression='gzip')
        print ''
        # Attribute
        for key,value in atr.iteritems():
            group.attrs[key] = value
        group.attrs.pop('ref_date')
        h5.close()

    print 'Done.'
    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  



