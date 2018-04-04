#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys
import argparse
import h5py
import numpy as np
from scipy.special import gamma
from pysar.utils import datetime as ptime, readfile, writefile, utils as ut
from pysar.objects import timeseries, geometry


############################################################################
TEMPLATE='''
## 8. Topographic Residual (DEM Error) Correction (Fattahi and Amelung, 2013, IEEE-TGRS)
## Specify stepFuncDate option if you know there are sudden displacement jump in your area,
## i.e. volcanic eruption, or earthquake, and check timeseriesStepModel.h5 afterward for their estimation.
pysar.topographicResidual               = auto  #[yes / no], auto for yes
pysar.topographicResidual.polyOrder     = auto  #[1-inf], auto for 2, poly order of temporal deformation model
pysar.topographicResidual.stepFuncDate  = auto  #[20080529,20100611 / no], auto for no, date of step jump
pysar.topographicResidual.excludeDate   = auto  #[20070321 / txtFile / no], auto for no, date exlcuded for error estimation
pysar.topographicResidual.phaseVelocity = auto  #[yes / no], auto for no - phase, use phase velocity for error estimation
'''

EXAMPLE='''example:
  dem_error.py  timeseries_ECMWF.h5 -g geometryRadar.h5 -t pysarApp_template.txt
  dem_error.py  timeseries_ECMWF.h5
'''

REFERENCE='''reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
  IEEE TGRS, 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
'''

def createParser():
    parser = argparse.ArgumentParser(description='DEM Error (Topographic Residual) Correction',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file to be corrrected')
    parser.add_argument('-g','--geometry', dest='geom_file',\
                        help='geometry file including datasets of incidence angle, slant range distance\n'+\
                             'and/or 3D perpendicular baseline')
    parser.add_argument('-o','--outfile', help='Output file name for corrected time series')

    parser.add_argument('--ex','--exclude', dest='ex_date', nargs='*', default=[],\
                        help='Exclude date(s) for DEM error estimation.\n'+\
                             'All dates will be corrected for DEM residual phase still.')
    parser.add_argument('-t','--template', dest='template_file',\
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-s','--step-date', dest='step_date', nargs='*', default=[],\
                        help='Date of step jump for temporal deformation model, i.e. date of earthquake/volcanic eruption')
    parser.add_argument('--phase-velocity', dest='phase_velocity', action='store_true',\
                        help='Use phase velocity instead of phase for inversion constrain.')
    parser.add_argument('-p','--poly-order', dest='poly_order', type=int, default=2,\
                        help='polynomial order number of temporal deformation model, default = 2')
    return parser


def cmdLineParse(iargs=None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    if inps.poly_order < 1:
        raise argparse.ArgumentTypeError("Minimum polynomial order is 1")
    return inps


############################################################################
def read_template2inps(template_file, inps=None):
    '''Read input template file into inps.ex_date'''
    if not inps:
        inps = cmdLineParse()
    template = readfile.read_template(template_file)

    # Read template option
    prefix = 'pysar.topographicResidual.'

    key = prefix+'polyOrder'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.poly_order = 2
        else:
            inps.poly_order = int(value)

    key = prefix+'excludeDate'
    if key in template.keys():
        value = template[key]
        if value not in ['auto','no']:
            value = value.replace(',',' ').split()
            value = ptime.yyyymmdd(value)
            inps.ex_date += value

    key = prefix+'stepFuncDate'
    if key in template.keys():
        value = template[key]
        if value not in ['auto','no']:
            value = value.replace(',',' ').split()
            value = ptime.yyyymmdd(value)
            inps.step_date += value

    key = prefix='phaseVelocity'
    if key in template.keys():
        value = template[key]
        if value.lower() not in ['auto','no']:
            inps.phase_velocity = True

    return inps


def read_exclude_date(inps, dateList):
    '''Read exclude dates info
    Inputs:
        exDateIn - list of string, date in YYMMDD or YYYYMMDD format,
                   or text file with date in it
        dateList - list of string, date in YYYYMMDD format
    Output:
        exDateOut - list of string, date in YYYYMMDD format
    '''
    ##Read exclude date input
    if inps.ex_date:
        tempList = []
        for d in inps.ex_date:
            if os.path.isfile(d):
                tempList += ptime.read_date_list(d)
            else:
                tempList.append(d)
        inps.ex_date = sorted(ptime.yyyymmdd(tempList))
        print('exclude the following dates for DEM error estimation: ({})\n{}'.format(len(inps.ex_date), inps.ex_date))
    else:
        inps.ex_date = []

    ##convert to mark array
    inps.dropDate = np.array([i not in inps.ex_date for i in dateList], dtype=np.bool_)
    if inps.poly_order > np.sum(inps.dropDate):
        raise ValueError("ERROR: input poly order {} > number of acquisition {}! Reduce it!".format(inps.poly_order,\
                                                                                                    np.sum(inps.dropDate)))
    return inps.dropDate


def design_matrix4deformation(inps):
    ## Date Info
    tsobj = timeseries(inps.timeseries_file)
    tsobj.open()

    ## Design matrix - temporal deformation model
    print('-'*50+'\ncorrect topographic phase residual (DEM error) using Fattahi and Amelung (2013, IEEE-TGRS)')
    msg = 'ordinal least squares (OLS) inversion with L2-norm minimization on: phase'
    if inps.phase_velocity:
        msg += ' velocity'
    if inps.rangeDist.size != 1:
        msg += ' (pixel-wisely)'
    print(msg)

    # 1. Polynomial - 2D matrix in size of (numDate, polyOrder+1)
    print("temporal deformation model: polynomial order = "+str(inps.poly_order))
    A_def = np.ones((tsobj.numDate, 1), np.float32)
    for i in range(inps.poly_order):
        Ai = np.array(tsobj.btemp**(i+1) / gamma(i+2), np.float32).reshape(-1,1)
        A_def = np.hstack((A_def, Ai))

    # 2. Step function - 2D matrix in size of (numDate, len(inps.step_date))
    if inps.step_date:
        print("temporal deformation model: step functions at "+str(inps.step_date))
        t_steps = ptime.yyyymmdd2years(inps.step_date)
        t = np.array(ptime.yyyymmdd2years(tsobj.dateList))
        for t_step in t_steps:
            Ai = np.array(t > t_step, np.float32).reshape(-1,1)
            A_def = np.hstack((A_def, Ai))
    print('-'*50)
    return A_def


def read_geometry(inps):
    tsobj = timeseries(inps.timeseries_file)
    tsobj.open(printMsg=False)
    ## 2D / 3D geometry
    if inps.geom_file:
        geomobj = geometry(inps.geom_file)
        geomobj.open()
        print('read 2D incidenceAngle,slantRangeDistance from {} file: {}'.format(geomobj.key, geomobj.file))
        inps.incAngle = geomobj.read(datasetName='incidenceAngle').flatten() * np.pi/180.
        inps.rangeDist = geomobj.read(datasetName='slantRangeDistance').flatten()
        if 'bperp' in geomobj.f[geomobj.key].keys():
            print('reading 3D bperp from {} file: {} ...'.format(geomobj.key, geomobj.file))
            inps.bperp = geomobj.read(datasetName='bperp').reshape((geomobj.numDate, -1))
            inps.bperp -= inps.bperp[tsobj.refIndex]
        else:
            print('read mean bperp from {} file'.format(tsobj.key))
            inps.bperp = tsobj.bperp.reshape((-1,1))
    ## 0D geometry
    else:
        print('read mean incidenceAngle,slantRangeDistance,bperp value from {} file'.format(tsobj.key))
        inps.incAngle = ut.incidence_angle(tsobj.metadata, dimension=0)
        inps.rangeDist = ut.range_distance(tsobj.metadata, dimension=0)
        inps.bperp = tsobj.bperp.reshape((-1,1))
    inps.sinIncAngle = np.sin(inps.incAngle)
    return inps


def estimate_dem_error_data(ts0, A0, inps):
    '''
    Parameters: ts0 : 2D np.array in size of (numDate, numPixel)
                    original time series displacement
                A0 : 2D np.array in size of (numDate, model_num)
                    design matrix in [A_geom, A_def]
                inps : Namespace with the following settings:
                    btemp : 2D np.array in size of (numDate, 1), temporal baseline
                    dropDate : 1D np.array in bool data type, mark the date used in the estimation
                    phase_velocity : bool, use phase history or phase velocity for minimization
    Returns:    deltaZ : 2D np.array in size of (1,       numPixel)
                    estimated DEM residual
                tsCor : 2D np.array in size of (numDate, numPixel)
                    corrected timeseries = tsOrig - deltaZphase
                tsRes : 2D np.array in size of (numDate, numPixel)
                    residual timeseries = tsOrig - deltaZphase - defModel
                stepEst : 2D np.array in size of (numStep, numPixel)
                    estimated step deformation
    Example:
        deltaZ, tsCor, tsRes = estimate_dem_error_data(ts, A, inps)
    '''
    if len(ts0.shape) == 1:
        ts0 = ts0.reshape(-1,1)
    numDate, numPixel = ts0.shape

    ##Prepare Design matrix A and observations ts for inversion
    A = A0[inps.dropDate,:]
    ts = ts0[inps.dropDate,:]
    btemp = inps.btemp[inps.dropDate,:]
    if inps.phase_velocity:
        ts = np.diff(ts, axis=0) / np.diff(btemp, axis=0)
        A = np.diff(A, axis=0) / np.diff(btemp, axis=0)

    ##Inverse using L-2 norm to get unknown parameters X = [deltaZ, vel, acc, deltaAcc, ...]
    X = np.linalg.pinv(A).dot(ts)   #equivalent to X = np.linalg.inv(A.T.dot(A)).dot(A.T).dot(ts)

    ##Prepare Outputs
    deltaZ = X[0,:]
    tsCor = ts0 - np.dot(A0[:,0].reshape(-1,1), deltaZ.reshape(1,-1))
    tsRes = ts0 - np.dot(A0, X)
    stepEst = None
    if inps.numStep > 0:
        s0 = inps.poly_order+2
        s1 = inps.poly_order+2+inps.numStep
        stepEst = X[s0:s1,:].reshape(inps.numStep,-1)
    return deltaZ, tsCor, tsRes, stepEst


def estimate_dem_error(inps, A_def):
    ##### Read Date Info
    tsobj = timeseries(inps.timeseries_file)
    tsobj.open()
    inps.btemp = tsobj.btemp.reshape(-1,1)
    inps.dropDate = read_exclude_date(inps, tsobj.dateList)
    inps.numStep = len(inps.step_date)

    ##### Read time-series Data
    print('reading time-series data ...')
    tsData = tsobj.read().reshape((tsobj.numDate, -1))

    ##---------------------------------------- Loop for L2-norm inversion  -----------------------------------##
    if inps.rangeDist.size == 1:
        A_geom = inps.bperp / (inps.rangeDist * inps.sinIncAngle)
        A = np.hstack((A_geom, A_def))
        deltaZ, tsCor, tsRes, stepEst = estimate_dem_error_data(tsData, A, inps)
    else:
        tsCor = np.zeros((tsobj.numDate, tsobj.numPixel), dtype=np.float32)
        tsRes = np.zeros((tsobj.numDate, tsobj.numPixel), dtype=np.float32)
        deltaZ = np.zeros(tsobj.numPixel, dtype=np.float32)
        #constC = np.zeros(tsobj.numPixel, dtype=np.float32)
        if inps.numStep > 0:
            stepModel = np.zeros((inps.numStep, tsobj.numPixel), dtype=np.float32)

        print('skip pixels with zero/nan value in geometry: incidence angle or range distance')
        mask = np.ones_like(inps.sinIncAngle, dtype=np.bool_)
        mask[inps.sinIncAngle == 0.] = 0
        mask[inps.rangeDist == 0.] = 0
        numPixel2inv = np.sum(mask)
        idxPixel2inv = np.where(mask)[0]
        print('number of pixels in   file: {}'.format(tsobj.numPixel))
        print('number of pixels to invert: {}'.format(numPixel2inv))

        progBar = ptime.progress_bar(maxValue=tsobj.numPixel)
        for i in range(numPixel2inv):
            progBar.update(i+1, every=1000, suffix='{}/{} pixels'.format(i+1,numPixel2inv))
            idx = idxPixel2inv[i]

            if inps.bperp.shape[1] == 1:
                bperp = inps.bperp
            else:
                bperp = inps.bperp[:,idx].reshape(-1,1)
            A_geom = bperp / (inps.rangeDist[idx] * inps.sinIncAngle[idx])
            A = np.hstack((A_geom, A_def))

            #ts = tsData[:,idx].reshape(tsobj.numDate,-1)
            deltaZ_i, tsCor_i, tsRes_i, stepEst_i = estimate_dem_error_data(tsData[:,idx], A, inps)
            deltaZ[idx:idx+1] = deltaZ_i
            tsCor[:,idx:idx+1] = tsCor_i
            tsRes[:,idx:idx+1] = tsRes_i
            if inps.numStep > 0:
                stepModel[:,idx:idx+1] = stepEst_i
        progBar.close()

    ##------------------------------------------------ Output  --------------------------------------------##
    print('-'*50)
    ##prepare for output
    tsCor = tsCor.reshape((tsobj.numDate, tsobj.length, tsobj.width))
    tsRes = tsRes.reshape((tsobj.numDate, tsobj.length, tsobj.width))
    deltaZ = deltaZ.reshape((tsobj.length, tsobj.width))
    if inps.numStep > 0:
        stepModel = stepModel.reshape((inps.numStep, tsobj.length, tsobj.width))
    atr = dict(tsobj.metadata)

    # 1. Estimated DEM error
    outfile = 'demErr.h5'
    print('writing >>> '+outfile)
    atr['FILE_TYPE'] = 'dem'
    atr['UNIT'] = 'm'
    writefile.write(deltaZ, atr, outfile)

    # 2. Time-series corrected for DEM error
    if not inps.outfile:
        inps.outfile = '{base}_demErr{ext}'.format(base=os.path.splitext(inps.timeseries_file)[0],
                                                   ext=os.path.splitext(inps.timeseries_file)[1])
    tsCorobj = timeseries(inps.outfile)
    tsCorobj.write2hdf5(data=tsCor, refFile=tsobj.file)

    # 3. Time-series of inversion residual
    tsResobj = timeseries(os.path.join(os.path.dirname(inps.outfile), 'timeseriesResidual.h5'))
    tsResobj.write2hdf5(data=tsRes, refFile=tsobj.file)

    # 4. Time-series of estimated Step Model
    if inps.numStep > 0:
        atr['FILE_TYPE'] = 'timeseries'
        atr.pop('REF_DATE')
        stepObj = timeseries(os.path.join(os.path.dirname(inps.outfile), 'timeseriesStepModel.h5'))
        stepObj.write2hdf5(data=stepModel, metadata=atr, dates=inps.step_date)
    return inps


############################################################################
def main(iargs=None):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    inps = read_geometry(inps)
    A_def = design_matrix4deformation(inps)
    inps = estimate_dem_error(inps, A_def)

    print('Done.')
    return

################################################################################
if __name__ == '__main__':
    main()  

