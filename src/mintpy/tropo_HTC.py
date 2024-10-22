############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import cv2
import numpy as np
from scipy.interpolate import interp2d

from mintpy.mask import mask_matrix
from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile

############################################################################

def read_topographic_data(geom_file, meta):
    print('read height & incidenceAngle from file: '+geom_file)
    dem = readfile.read(geom_file, datasetName='height', print_msg=False)[0]
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle', print_msg=False)[0]
    dem *= 1.0/np.cos(inc_angle*np.pi/180.0)

    ref_y = int(meta['REF_Y'])
    ref_x = int(meta['REF_X'])
    dem -= dem[ref_y, ref_x]

    # Design matrix for elevation v.s. phase
    # dem = dem.flatten()
    return dem

def read_velocity_data(velo_file, meta):
    print('read velocity from file: '+velo_file)
    velocity = readfile.read(velo_file, datasetName='velocity', print_msg=False)[0] #TODO

    return velocity


def estimate_local_slope(dem, ts_data, inps, n_ref):
    """Estimate local slope based on texture correlation for each acquisition of timeseries
    Parameters: dem     : 2D array in size of (          length, width)
                ts_data : 3D array in size of (num_date, length, width)
                inps    : Namespace
                n_ref   : Index of reference SAR image in timeseries
    Returns:    X       : 2D array in size of (poly_num+1, num_date)  
    """#TODO

    """Filtering parameters for obtaining high-frequency texture correlation"""
    w1 = 9
    w2 = 13
    res_step = 0.5
    Iteration = 10
    rg = 40
    step = 0.0001

    velocity = readfile.read(inps.velo_file, datasetName='velocity')[0]

    print('reading mask from file: '+inps.mask_file)
    mask_coh = readfile.read(inps.mask_file, datasetName='temporalCoherence')[0]
    maskdef = 0 #TODO
    mask_def = create_deformation_mask(dem, velocity, maskdef)
    mask = mask_def * mask_coh


    Na, Nr, N = ts_data.shape
    #Na, Nr = dem.shape
    W = inps.windowsize
    w = (W-1)/2
    overlap = round(inps.overlapratio*(2*w+1))
    Na_C = np.arange(w + 1, Na - w - overlap + 1, 2 * w - overlap)
    Nr_C = np.arange(w + 1, Nr - w - overlap + 1, 2 * w - overlap)
    
    k_LLF = np.zeros(len(Na_C), len(Nr_C), N)
    d_LLF = np.zeros(len(Na_C), len(Nr_C), N)
    k_HTC = np.zeros(len(Na_C), len(Nr_C), N)

    for na in range(0, len(Na_C)):
        for nr in range(0, len(Nr_C)):
            # patch
            nac = Na_C[na]
            nrc = Nr_C[nr]
            u = nac - w
            U = np.where(u < 1, 1, u)
            d = nac + w
            D = np.where(d > Na, Na, d)
            l = nrc - w
            L = np.where(l < 1, 1, l)
            r = nrc + w
            R = np.where(r > Nr, Nr, r)

            if na == len(Na_C):
                D = Na
            if nr == len(Nr_C):
                R = Nr
            tmp = np.full((Na, Nr), np.nan)
            tmp[U:D, L:R] = 1
            mask_process = mask
            mask_process[np.isnan(tmp)] = np.nan

            # solve
            result_compare = np.zeros(4, N)
            for n in range(0, N):
                if n == n_ref:
                    continue

                # ----------- initial value ----------- #
                # mask
                mask_std = np.ones(Na, Nr) * mask_process

                # iterative linear fitting
                # res+step = 0.5 # step [rad]
                # Iteration = 10
                for iter in range(1, Iteration + 1):
                    # linear fitting
                    phase_tmp = ts_data[:, :, n] * mask_process
                    dem_tmp = dem * mask_std
                    coe = np.polyfit(dem[~np.isnan(mask_std)], phase_tmp[~np.isnan(mask_std)], 1)
                    cor_tmp = phase_tmp - coe[0]*dem_tmp - coe[1]

                    # result recording
                    if iter == 1:
                        result_compare[0, n] = coe[0]
                        result_compare[3, n] = coe[1]
                    
                    # mask uploading
                    max_tmp = np.max(np.max(np.abs(cor_tmp))) #TODO
                    max_tmp = np.where(max_tmp > res_step, max_tmp - res_step, max_tmp)
                    mask_std[np.abs(cor_tmp) > max_tmp] = np.nan

                    if sum(np.sum(~np.isnan(mask_std))) < sum(np.sum(~np.isnan(mask_process))) / 10: #TODO
                        break

                # ----------- texture correlation ----------- #        
                mask_tmp = mask[U:D, L:R]
                
                A = dem[U:D, L:R]
                A_LP = cv2.GaussianBlur(A, (w2, w2), sigmaX=w1, sigmaY=w1)
                A = A - A_LP
                A_line = A[~np.isnan(mask_tmp)]
                A_line = A_line / np.linalg.norm(A_line)

                # range = 40
                # step = 0.0001
                # left
                k_left = -rg * step + coe[0]
                phase_ts_Scor_tmp = ts_data[:, :, n] - k_left * dem
                C = phase_ts_Scor_tmp[U:D, L:R]
                C[np.isnan(C)] = 0
                C_LP = cv2.GaussianBlur(A, (w2, w2), sigmaX=w1, sigmaY=w1)
                C = C - C_LP
                C_line = C[~np.isnan(mask_tmp)]
                C_line = C_line / np.linalg.norm(C_line)
                conv_AC = np.sum(A_line * C_line)
                record_left = np.abs(conv_AC)
                # right
                k_right = rg * step + coe[0]
                phase_ts_Scor_tmp = ts_data[:, :, n] - k_right * dem
                C = phase_ts_Scor_tmp[U:D, L:R]
                C[np.isnan(C)] = 0
                C_LP = cv2.GaussianBlur(A, (w2, w2), sigmaX=w1, sigmaY=w1)
                C = C - C_LP
                C_line = C[~np.isnan(mask_tmp)]
                C_line = C_line / np.linalg.norm(C_line)
                conv_AC = np.sum(A_line * C_line)
                record_right = np.abs(conv_AC)

                if record_left >= record_right:
                    k = np.arange(0, rg + 1) * step + coe[0]
                    record = np.zeros(len(k))
                else:
                    k = np.arange(-rg, 1) * step + coe[0]
                    record = np.zeros(len(k))

                for i in range(0, len(k)):
                    phase_ts_Scor_tmp = ts_data[:, :, n] - k[i] * dem
                    C = phase_ts_Scor_tmp[U:D, L:R]
                    C[np.isnan(C)] = 0
                    C_LP = cv2.GaussianBlur(A, (w2, w2), sigmaX=w1, sigmaY=w1)
                    C = C - C_LP
                    C_line = C[~np.isnan(mask_tmp)]
                    C_line = C_line / np.linalg.norm(C_line)

                    conv_AC = np.sum(A_line * C_line)
                    record[i] = np.abs(conv_AC)

                # slope
                index = np.where(record == np.min(record))[0]
                result_compare[1, n] = coe[0]
                result_compare[2, n] = k[index]
            
            # result recording
            k_LLF[na, nr, :] = result_compare[0, :]
            d_LLF[na, nr, :] = result_compare[3, :]
            k_HTC[na, nr, :] = result_compare[2, :]
    return [k_LLF, d_LLF, k_HTC]

def slope_interpolation(ts_data, inps, k_HTC):
    """Filtering parameters for obtaining high-frequency texture correlation"""
    sigma_slope = 7
    w_slope = 7   

    Na, Nr, N = ts_data.shape
    #Na, Nr = dem.shape
    W = inps.windowsize
    w = (W-1)/2
    overlap = round(inps.overlapratio*(2*w+1))
    Na_C = np.arange(w + 1, Na - w - overlap + 1, 2 * w - overlap)
    Nr_C = np.arange(w + 1, Nr - w - overlap + 1, 2 * w - overlap)
    # grid
    Y = Na_C.reshape(-1, 1) * np.ones((1, len(Nr_C)))
    X = np.ones((len(Na_C), 1)) * Nr_C.reshape(-1, 1)

    Xq, Yq = np.meshgrid(np.arange(1, Nr+1), np.arange(1, Na+1))

    # K_HTC
    k_HTC_interp = np.zeros(Na, Nr, N)
    for n in range(0, N):
        k_HTC_filt = cv2.GaussianBlur(k_HTC[:, :, n], (w_slope, w_slope), sigmaX=sigma_slope, sigmaY=sigma_slope)
        f_interp = interp2d(X, Y, k_HTC_filt, kind='spline')
        k_HTC_interp[:, :, n] = f_interp(Xq, Yq).reshape(k_HTC_interp[:, :, n].shape) #TODO 
    
    return k_HTC_interp

def intercept_filtering(dem, ts_data, inps, k_HTC_interp, meta):
    """Filtering parameters for obtaining high-frequency texture correlation"""
    sigma_intercept = 251 
    w_intercept = 251 

    ref_y = int(meta['REF_Y'])
    ref_x = int(meta['REF_X'])

    phase_ts_HTC_low = ts_data
    intercept = np.zeros(phase_ts_HTC_low.shape)
    for n in range(0, N):
        tmp = ts_data[:, :, n] - k_HTC_interp[:, :, n] * dem
        tmp_filt = cv2.GaussianBlur(tmp, (w_intercept, w_intercept), sigmaX=sigma_intercept, sigmaY=sigma_intercept)
        tmp = tmp - tmp_filt
        phase_ts_HTC_low[:, :, n] = tmp
        intercept[:, :, n] = tmp_filt
    
    phase_ts_HTC_low = phase_ts_HTC_low - phase_ts_HTC_low[ref_y, ref_x, :] #TODO

    return phase_ts_HTC_low
    
    

def create_deformation_mask(dem, rate, maskdef):
    if maskdef == 1:
        mask_def = rate
        mask_def[np.abs(rate)>1.5] = np.nan
        mask_def[~np.isnan(mask_def)] = 1
    else:
        mask_def = np.ones(dem.shape) 
    return mask_def

############################################################################

def run_tropo_HTC(inps):

    # read time-series data
    ts_obj = timeseries(inps.timeseries_file)
    ts_obj.open()
    ts_data = ts_obj.read()
    inps.date_list = list(ts_obj.dateList)
    n_ref = inps.date_list.index(ts_obj.metadata['REF_DATE'])

    # read topographic data (DEM)
    dem = read_topographic_data(inps.geom_file, ts_obj.metadata)
    #TODO
    # slope estimation
    k_HTC = estimate_local_slope(dem, ts_data, inps, n_ref)[2]
    # slope interpolation
    k_HTC_interp = slope_interpolation(ts_data, inps, k_HTC)
    # intercept filtering
    ts_HTC_low = intercept_filtering(dem, ts_data, inps, k_HTC_interp, ts_obj.metadata)

    # write corrected time-series file
    meta = dict(ts_obj.metadata)
    meta['mintpy.troposphericDelay.polyOrder'] = str(inps.poly_order)
    if not inps.outfile:
        fbase = os.path.splitext(inps.timeseries_file)[0]
        inps.outfile - f'{fbase}_tropHTC.h5'
    
    writefile.write(
        ts_data,
        out_file=inps.outfile,
        metadata=meta,
        ref_file=inps.timeseries_file
    )

    return

