#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import time
import argparse
import numpy as np
from scipy.special import gamma
from pysar.utils import ptime, readfile, writefile, utils as ut
from pysar.objects import timeseries, geometry


############################################################################
TEMPLATE = """
## 8. Topographic Residual (DEM Error) Correction (Fattahi and Amelung, 2013, IEEE-TGRS)
## Specify stepFuncDate option if you know there are sudden displacement jump in your area,
## i.e. volcanic eruption, or earthquake, and check timeseriesStepModel.h5 afterward for their estimation.
pysar.topographicResidual               = auto  #[yes / no], auto for yes
pysar.topographicResidual.polyOrder     = auto  #[1-inf], auto for 2, poly order of temporal deformation model
pysar.topographicResidual.stepFuncDate  = auto  #[20080529,20100611 / no], auto for no, date of step jump
pysar.topographicResidual.excludeDate   = auto  #[20070321 / txtFile / no], auto for no, date exlcuded for error estimation
pysar.topographicResidual.phaseVelocity = auto  #[yes / no], auto for no - phase, use phase velocity for error estimation
"""

EXAMPLE = """example:
  # correct DEM error with pixel-wise geometry parameters
  dem_error.py  timeseries_ECMWF.h5 -g INPUTS/geometryRadar.h5 -t pysarApp_template.txt

  # correct DEM error with mean geometry parameters
  dem_error.py  timeseries_ECMWF.h5

  # get time-series of estimated deformation model
  diff.py timeseries_ECMWF_demErr.h5 timeseriesResidual.h5 -o timeseriesDefModel.h5
"""

REFERENCE = """reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
  IEEE TGRS, 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='DEM Error (Topographic Residual) Correction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file',
                        help='Timeseries file to be corrrected')
    parser.add_argument('-g', '--geometry', dest='geom_file',
                        help='geometry file including datasets:\n'+
                             'incidence angle\n'+
                             'slant range distance\n' +
                             'and/or 3D perpendicular baseline')
    parser.add_argument('-o', '--outfile',
                        help='Output file name for corrected time series')

    parser.add_argument('--ex', '--exclude', dest='ex_date', nargs='*', default=[],
                        help='Exclude date(s) for DEM error estimation.\n' +
                             'All dates will be corrected for DEM residual phase still.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-s', '--step-date', dest='step_date', nargs='*', default=[],
                        help='Date of step jump for temporal deformation model,'+
                             ' i.e. date of earthquake/volcanic eruption')
    parser.add_argument('--phase-velocity', dest='min_phase_velocity', action='store_true',
                        help='Use phase velocity instead of phase for inversion constrain.')
    parser.add_argument('-p', '--poly-order', dest='poly_order', type=int, default=2,
                        help='polynomial order number of temporal deformation model, default = 2')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if inps.poly_order < 1:
        raise argparse.ArgumentTypeError("Minimum polynomial order is 1")
    return inps


############################################################################
def read_template2inps(template_file, inps=None):
    """Read input template file into inps.ex_date"""
    if not inps:
        inps = cmd_line_parse()
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
        if value not in ['auto', 'no']:
            value = value.replace(',', ' ').split()
            value = ptime.yyyymmdd(value)
            inps.ex_date += value

    key = prefix+'stepFuncDate'
    if key in template.keys():
        value = template[key]
        if value not in ['auto', 'no']:
            value = value.replace(',', ' ').split()
            value = ptime.yyyymmdd(value)
            inps.step_date += value

    key = prefix+'phaseVelocity'
    if key in template.keys():
        value = template[key]
        if value.lower() not in ['auto', 'no']:
            inps.min_phase_velocity = True

    return inps


def read_exclude_date(ex_date_list, date_list_all):
    """Read exclude dates info
    Parameters: ex_date_list  : list of string, date in YYMMDD or YYYYMMDD format,
                                or text file with date in it
                date_list_all : list of string, date in YYYYMMDD format
    Returns:    drop_date     : 1D array of bool in size of (num_date,)
    """
    # Read exclude date input
    if ex_date_list:
        tempList = []
        for d in ex_date_list:
            if os.path.isfile(d):
                tempList += ptime.read_date_list(d)
            else:
                tempList.append(d)
        ex_date_list = sorted(ptime.yyyymmdd(tempList))
        print(('exclude the following dates for DEM error estimation:'
               ' ({})\n{}').format(len(ex_date_list), ex_date_list))
    else:
        ex_date_list = []

    # convert to mark array
    drop_date = np.array([i not in ex_date_list for i in date_list_all],
                         dtype=np.bool_)
    return drop_date


def design_matrix4deformation(inps):
    # Date Info
    ts_obj = timeseries(inps.timeseries_file)
    ts_obj.open()

    # Design matrix - temporal deformation model
    print('-'*50)
    print('correct topographic phase residual (DEM error) using Fattahi and Amelung (2013, IEEE-TGRS)')
    msg = 'ordinal least squares (OLS) inversion with L2-norm minimization on: phase'
    if inps.min_phase_velocity:
        msg += ' velocity'
    if inps.rangeDist.size != 1:
        msg += ' (pixel-wisely)'
    print(msg)

    tbase = np.array(ts_obj.tbase, np.float32) / 365.25

    # 1. Polynomial - 2D matrix in size of (numDate, polyOrder+1)
    print("temporal deformation model: polynomial order = "+str(inps.poly_order))
    A_def = np.ones((ts_obj.numDate, 1), np.float32)
    for i in range(inps.poly_order):
        Ai = np.array(tbase**(i+1) / gamma(i+2), np.float32).reshape(-1, 1)
        A_def = np.hstack((A_def, Ai))

    # 2. Step function - 2D matrix in size of (numDate, len(inps.step_date))
    if inps.step_date:
        print("temporal deformation model: step functions at "+str(inps.step_date))
        t_steps = ptime.yyyymmdd2years(inps.step_date)
        t = np.array(ptime.yyyymmdd2years(ts_obj.dateList))
        for t_step in t_steps:
            Ai = np.array(t > t_step, np.float32).reshape(-1, 1)
            A_def = np.hstack((A_def, Ai))
    print('-'*50)
    return A_def


def read_geometry(inps):
    ts_obj = timeseries(inps.timeseries_file)
    ts_obj.open(print_msg=False)
    # 2D / 3D geometry
    if inps.geom_file:
        geom_obj = geometry(inps.geom_file)
        geom_obj.open()
        print(('read 2D incidenceAngle,slantRangeDistance from {} file:'
               ' {}').format(geom_obj.name, os.path.basename(geom_obj.file)))
        inps.incAngle = geom_obj.read(datasetName='incidenceAngle', print_msg=False).flatten()
        inps.rangeDist = geom_obj.read(datasetName='slantRangeDistance', print_msg=False).flatten()
        if 'bperp' in geom_obj.datasetNames:
            print('read 3D bperp from {} file: {} ...'.format(geom_obj.name, os.path.basename(geom_obj.file)))
            inps.pbase = geom_obj.read(datasetName='bperp', print_msg=False).reshape((geom_obj.numDate, -1))
            inps.pbase -= inps.pbase[ts_obj.refIndex]
        else:
            print('read mean bperp from {} file'.format(ts_obj.name))
            inps.pbase = ts_obj.pbase.reshape((-1, 1))

    # 0D geometry
    else:
        print('read mean incidenceAngle,slantRangeDistance,bperp value from {} file'.format(ts_obj.name))
        inps.incAngle = ut.incidence_angle(ts_obj.metadata, dimension=0)
        inps.rangeDist = ut.range_distance(ts_obj.metadata, dimension=0)
        inps.pbase = ts_obj.pbase.reshape((-1, 1))

    inps.sinIncAngle = np.sin(inps.incAngle * np.pi / 180.)
    return inps


def estimate_dem_error(ts0, A0, tbase, drop_date=None, min_phase_velocity=False, num_step=0):
    """Estimate DEM error with least square optimization.
    Parameters: ts0 : 2D np.array in size of (numDate, numPixel), original time series displacement
                A0  : 2D np.array in size of (numDate, model_num), design matrix in [A_geom, A_def]
                tbase : 2D np.array in size of (numDate, 1), temporal baseline
                drop_date : 1D np.array in bool data type, mark the date used in the estimation
                min_phase_velocity : bool, use phase history or phase velocity for minimization
    Returns:    delta_z: 2D np.array in size of (1,       numPixel) estimated DEM residual
                ts_cor : 2D np.array in size of (numDate, numPixel),
                            corrected timeseries = tsOrig - delta_z_phase
                ts_res : 2D np.array in size of (numDate, numPixel),
                            residual timeseries = tsOrig - delta_z_phase - defModel
    Example:    delta_z, ts_cor, ts_res = estimate_dem_error(ts, A, tbase, drop_date)
    """
    if len(ts0.shape) == 1:
        ts0 = ts0.reshape(-1, 1)
    if drop_date is None:
        drop_date = np.ones(ts0.shape[0], np.bool_)

    # Prepare Design matrix A and observations ts for inversion
    A = A0[drop_date, :]
    ts = ts0[drop_date, :]
    if min_phase_velocity:
        tbase = tbase[drop_date, :]
        ts = np.diff(ts, axis=0) / np.diff(tbase, axis=0)
        A = np.diff(A, axis=0) / np.diff(tbase, axis=0)

    # Inverse using L-2 norm to get unknown parameters X
    # X = [delta_z, constC, vel, acc, deltaAcc, ..., step1, step2, ...]
    # equivalent to X = np.linalg.inv(A.T.dot(A)).dot(A.T).dot(ts)
    X = np.linalg.pinv(A).dot(ts)

    # Prepare Outputs
    delta_z = X[0, :]
    ts_cor = ts0 - np.dot(A0[:, 0].reshape(-1, 1), delta_z.reshape(1, -1))
    ts_res = ts0 - np.dot(A0, X)

    step_def = None
    if num_step > 0:
        step_def = X[-1*num_step:, :].reshape(num_step, -1)

    # for debug
    debug_mode = False
    if debug_mode:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8))
        ts_all = np.hstack((ts0, ts_res, ts_cor))
        ymin = np.min(ts_all)
        ymax = np.max(ts_all)
        ax1.plot(ts0, '.');           ax1.set_ylim((ymin, ymax)); ax1.set_title('Original  Timeseries')
        ax2.plot(ts_cor, '.');        ax2.set_ylim((ymin, ymax)); ax2.set_title('Corrected Timeseries')
        ax3.plot(ts_res, '.');        ax3.set_ylim((ymin, ymax)); ax3.set_title('Fitting Residual')
        ax4.plot(ts_cor-ts_res, '.'); ax4.set_ylim((ymin, ymax)); ax4.set_title('Fitted Deformation Model')
        plt.show()

    return delta_z, ts_cor, ts_res, step_def


def correct_dem_error(inps, A_def):
    """Correct DEM error of input timeseries file"""
    # Read Date Info
    ts_obj = timeseries(inps.timeseries_file)
    ts_obj.open()
    num_date = ts_obj.numDate
    num_pixel = ts_obj.numPixel
    tbase = np.array(ts_obj.tbase, np.float32) / 365.25

    num_step = len(inps.step_date)
    drop_date = read_exclude_date(inps.ex_date, ts_obj.dateList)
    if inps.poly_order > np.sum(drop_date):
        raise ValueError(("ERROR: input poly order {} > number of acquisition {}!"
                          " Reduce it!").format(inps.poly_order, np.sum(drop_date)))

    # Read time-series Data
    print('reading time-series data ...')
    ts_data = ts_obj.read().reshape((num_date, -1))

    ##-------------------------------- Loop for L2-norm inversion  --------------------------------##
    print('inverting DEM error ...')
    if inps.rangeDist.size == 1:
        A_geom = inps.pbase / (inps.rangeDist * inps.sinIncAngle)
        A = np.hstack((A_geom, A_def))
        (delta_z,
         ts_cor,
         ts_res,
         step_model) = estimate_dem_error(ts_data,
                                          A,
                                          tbase=tbase,
                                          drop_date=drop_date,
                                          min_phase_velocity=inps.min_phase_velocity,
                                          num_step=num_step)

    else:
        ts_cor = np.zeros((num_date, num_pixel), dtype=np.float32)
        ts_res = np.zeros((num_date, num_pixel), dtype=np.float32)
        delta_z = np.zeros(num_pixel, dtype=np.float32)
        #constC = np.zeros(num_pixel, dtype=np.float32)
        if num_step > 0:
            step_model = np.zeros((num_step, num_pixel), dtype=np.float32)

        # mask
        print('skip pixels with zero/nan value in geometry: incidence angle or range distance')
        mask = np.multiply(inps.sinIncAngle != 0., inps.rangeDist != 0.)
        print('skip pixels with zero value in all acquisitions')
        ts_mean = np.nanmean(ts_data, axis=0)
        mask *= ts_mean != 0.
        del ts_mean

        num_pixel2inv = np.sum(mask)
        idx_pixel2inv = np.where(mask)[0]
        print(('number of pixels to invert: {} out of {}'
               ' ({:.1f}%)').format(num_pixel2inv,
                                    num_pixel,
                                    num_pixel2inv/num_pixel*100))

        # update data matrix to save memory and IO
        ts_data = ts_data[:, mask]
        inps.rangeDist = inps.rangeDist[mask]
        inps.sinIncAngle = inps.sinIncAngle[mask]
        if inps.pbase.shape[1] != 1:
            inps.pbase = inps.pbase[:, mask]

        # loop pixel by pixel
        prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
        for i in range(num_pixel2inv):
            prog_bar.update(i+1, every=1000, suffix='{}/{}'.format(i+1, num_pixel2inv))
            idx = idx_pixel2inv[i]

            # design matrix
            if inps.pbase.shape[1] == 1:
                pbase = inps.pbase
            else:
                pbase = inps.pbase[:, i].reshape(-1, 1)
            A_geom = pbase / (inps.rangeDist[i] * inps.sinIncAngle[i])
            A = np.hstack((A_geom, A_def))

            (delta_z_i,
             ts_cor_i,
             ts_res_i,
             step_model_i) = estimate_dem_error(ts_data[:, i],
                                                A,
                                                tbase=tbase,
                                                drop_date=drop_date,
                                                min_phase_velocity=inps.min_phase_velocity,
                                                num_step=num_step)
            delta_z[idx:idx+1] = delta_z_i
            ts_cor[:, idx:idx+1] = ts_cor_i
            ts_res[:, idx:idx+1] = ts_res_i
            if num_step > 0:
                step_model[:, idx:idx+1] = step_model_i
        prog_bar.close()
    del ts_data

    ##---------------------------------------- Output  -----------------------------------------##
    # prepare for output
    ts_cor = ts_cor.reshape((num_date, ts_obj.length, ts_obj.width))
    ts_res = ts_res.reshape((num_date, ts_obj.length, ts_obj.width))
    delta_z = delta_z.reshape((ts_obj.length, ts_obj.width))
    if num_step > 0:
        step_model = step_model.reshape((num_step, ts_obj.length, ts_obj.width))
    atr = dict(ts_obj.metadata)

    # 1. Estimated DEM error
    outfile = 'demErr.h5'
    atr['FILE_TYPE'] = 'dem'
    atr['UNIT'] = 'm'
    writefile.write(delta_z, out_file=outfile, metadata=atr)

    # 2. Time-series corrected for DEM error
    if not inps.outfile:
        inps.outfile = '{}_demErr.h5'.format(os.path.splitext(inps.timeseries_file)[0])
    ts_cor_obj = timeseries(inps.outfile)
    ts_cor_obj.write2hdf5(data=ts_cor, refFile=ts_obj.file)

    # 3. Time-series of inversion residual
    ts_res_obj = timeseries(os.path.join(os.path.dirname(inps.outfile), 'timeseriesResidual.h5'))
    ts_res_obj.write2hdf5(data=ts_res, refFile=ts_obj.file)

    # 4. Time-series of estimated Step Model
    if num_step > 0:
        atr['FILE_TYPE'] = 'timeseries'
        atr.pop('REF_DATE')
        step_obj = timeseries(os.path.join(os.path.dirname(inps.outfile), 'timeseriesStepModel.h5'))
        step_obj.write2hdf5(data=step_model, metadata=atr, dates=inps.step_date)

    ## 5. Time-series of estimated Deformation Model = poly model + step model
    #ts_def_obj = timeseries(os.path.join(os.path.dirname(inps.outfile), 'timeseriesDefModel.h5'))
    #ts_def_obj.write2hdf5(data=ts_cor - ts_res, refFile=ts_obj.file)
    #del ts_cor, ts_res

    return inps


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    start_time = time.time()
    inps = read_geometry(inps)
    A_def = design_matrix4deformation(inps)
    inps = correct_dem_error(inps, A_def)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return


################################################################################
if __name__ == '__main__':
    main()
