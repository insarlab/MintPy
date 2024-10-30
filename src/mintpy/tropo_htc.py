############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import numpy as np
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile

############################################################################

def read_topographic_data(geom_file, meta):
    print('read height & incidenceAngle from file: '+geom_file)
    dem = readfile.read(geom_file, datasetName='height', print_msg=False)[0]
    Na, Nr = dem.shape
    dataLine = dem.flatten()
    
    # 找到非NaN值的索引
    index = np.where(~np.isnan(dataLine))[0]
    x = np.ceil(index / Na).astype(int)
    y = np.mod(index, Na).astype(int)
    
    # 创建插值网格
    xi = np.arange(0, Nr)
    yi = np.arange(0, Na)
    
    # 使用三次样条插值
    dataInp = griddata((x, y), dataLine[index], (xi[None, :], yi[:, None]), method='cubic')
    dem[np.isnan(dem)] = dataInp[np.isnan(dem)]
    return dem

def read_velocity_data(velo_file, meta):
    print('read velocity from file: '+velo_file)
    velocity = readfile.read(velo_file, datasetName='velocity', print_msg=False)[0] #TODO

    return velocity


def estimate_local_slope(dem, ts_data, inps, n_ref, meta):
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
    lamda = 0.3284035 # TODO
    # w2 = 2 * int(truncate * sigma + 0.5) + 1
    truncate = ((w2 - 1)/2 - 0.5)/w1
    
    # TODO
    #ref_y = int(meta['REF_Y'])
    #ref_x = int(meta['REF_X'])
    ref_y = 282
    ref_x = 204

    velocity = readfile.read(inps.velo_file, datasetName='velocity')[0]

    print('reading mask from file: '+inps.mask_file)
    mask_coh = readfile.read(inps.mask_file, datasetName='temporalCoherence')[0]
    mask_coh = np.where(mask_coh < 0.87, np.nan, 1) #TODO
    
    maskdef = 0 #TODO
    mask_def = create_deformation_mask(dem, velocity, maskdef)
    mask = mask_def * mask_coh

    N, Na, Nr = ts_data.shape
    #Na, Nr = dem.shape
    W = inps.windowsize
    w = (W-1)/2
    overlap = round(inps.overlapratio*(2*w+1))
    Na_C = np.arange(w + 1, Na + 1, 2 * w - overlap)
    Nr_C = np.arange(w + 1, Nr + 1, 2 * w - overlap)
    
    k_LLF = np.zeros((N, len(Na_C), len(Nr_C)))
    d_LLF = np.zeros((N, len(Na_C), len(Nr_C)))
    k_htc = np.zeros((N, len(Na_C), len(Nr_C)))

    
    ts_data = 4 * np.pi / lamda * ts_data[:N, :, :]
    reference_value = ts_data[:, ref_y-1, ref_x-1]
    ts_data = ts_data - reference_value[:, np.newaxis, np.newaxis]
    ts_data[np.isnan(ts_data)] = 0

    for na, nac in enumerate(Na_C):
        for nr, nrc in enumerate(Nr_C):
            # patch
            u = nac - w
            U = np.where(u < 1, 1, u)
            d = nac + w
            D = np.where(d > Na, Na, d)
            l = nrc - w
            L = np.where(l < 1, 1, l)
            r = nrc + w
            R = np.where(r > Nr, Nr, r)

            if na == len(Na_C) - 1:
                D = Na
            if nr == len(Nr_C) - 1:
                R = Nr
            U = int(U)
            D = int(D)
            L = int(L)
            R = int(R)    
            tmp = np.full((Na, Nr), np.nan)
            tmp[U-1:D, L-1:R] = 1
            mask_process = mask.copy()
            mask_process[np.isnan(tmp)] = np.nan

            # solve
            result_compare = np.zeros((N, 4))
            for n in range(N):
                if n == n_ref:
                    continue

                # ----------- initial value ----------- #
                # mask
                mask_std = np.ones((Na, Nr)) * mask_process

                # iterative linear fitting
                # res+step = 0.5 # step [rad]
                # Iteration = 10
                for iter in range(Iteration):
                    # linear fitting
                    phase_tmp = ts_data[n, :, :] * mask_std
                    dem_tmp = dem * mask_std
                    valid_idx = ~np.isnan(phase_tmp) & ~np.isnan(dem_tmp)
                    if np.nansum(valid_idx) == 0:
                        continue
                    coe = np.polyfit(dem[valid_idx], phase_tmp[valid_idx], 1)
                    cor_tmp = phase_tmp - coe[0]*dem_tmp - coe[1]

                    # result recording
                    if iter == 0:
                        result_compare[n, 0] = coe[0]
                        result_compare[n, 3] = coe[1]
                    
                    # mask uploading
                    max_tmp = np.nanmax(np.abs(cor_tmp)) #TODO
                    max_tmp = np.where(max_tmp > res_step, max_tmp - res_step, max_tmp)
                    mask_std[np.abs(cor_tmp) > max_tmp] = np.nan

                    if np.nansum(~np.isnan(mask_std)) < np.nansum(~np.isnan(mask_process)) / 10: #TODO
                        break

                # ----------- texture correlation ----------- #        
                mask_tmp = mask[U-1:D, L-1:R]
                
                A = dem[U-1:D, L-1:R]
                A_LP = scipy.ndimage.gaussian_filter(A, sigma=w1, truncate=truncate, mode='nearest')
                A = A - A_LP
                A_line = A[~np.isnan(mask_tmp)]
                A_line = A_line / np.linalg.norm(A_line)
                
                # range = 40
                # step = 0.0001
                # left
                k_left = -rg * step + coe[0]
                phase_ts_Scor_tmp = ts_data[n, :, :] - k_left * dem
                C = phase_ts_Scor_tmp[U-1:D, L-1:R]
                C[np.isnan(C)] = 0
                C_LP = scipy.ndimage.gaussian_filter(C, sigma=w1, truncate=truncate, mode='nearest')
                C = C - C_LP
                C_line = C[~np.isnan(mask_tmp)]
                C_line = C_line / np.linalg.norm(C_line)
                conv_AC = np.nansum(A_line * C_line)
                record_left = np.abs(conv_AC)
                # right
                k_right = rg * step + coe[0]
                phase_ts_Scor_tmp = ts_data[n, :, :] - k_right * dem
                C = phase_ts_Scor_tmp[U-1:D, L-1:R]
                C[np.isnan(C)] = 0
                C_LP = scipy.ndimage.gaussian_filter(C, sigma=w1, truncate=truncate, mode='nearest')
                C = C - C_LP
                C_line = C[~np.isnan(mask_tmp)]
                C_line = C_line / np.linalg.norm(C_line)
                conv_AC = np.nansum(A_line * C_line)
                record_right = np.abs(conv_AC)

                if record_left >= record_right:
                    k = np.arange(0, rg + 1) * step + coe[0]
                    record = np.zeros(len(k))
                else:
                    k = np.arange(-rg, 1) * step + coe[0]
                    record = np.zeros(len(k))

                for i in range(0, len(k)):
                    phase_ts_Scor_tmp = ts_data[n, :, :] - k[i] * dem
                    C = phase_ts_Scor_tmp[U-1:D, L-1:R]
                    C[np.isnan(C)] = 0
                    C_LP = scipy.ndimage.gaussian_filter(C, sigma=w1, truncate=truncate, mode='nearest')
                    C = C - C_LP
                    C_line = C[~np.isnan(mask_tmp)]
                    C_line = C_line / np.linalg.norm(C_line)

                    conv_AC = np.nansum(A_line * C_line)
                    record[i] = np.abs(conv_AC)

                # slope
                #index = np.where(record == np.min(record))[0]
                index = np.argmin(record)
                result_compare[n, 1] = coe[0]
                result_compare[n, 2] = k[index]
            
            # result recording
            k_LLF[:, na, nr] = result_compare[:, 0]
            d_LLF[:, na, nr] = result_compare[:, 3]
            k_htc[:, na, nr] = result_compare[:, 2]
    return [k_LLF, d_LLF, k_htc]

def slope_interpolation(ts_data, inps, k_htc):
    """Filtering parameters for obtaining high-frequency texture correlation"""
    sigma_slope = 7
    w_slope = 7   
    truncate_slope = ((w_slope - 1)/2 - 0.5)/sigma_slope
    N, Na, Nr = ts_data.shape
    #Na, Nr = dem.shape
    W = inps.windowsize
    w = (W-1)/2
    overlap = round(inps.overlapratio*(2*w+1))
    Na_C = np.arange(w + 1, Na + 1, 2 * w - overlap)
    Nr_C = np.arange(w + 1, Na + 1, 2 * w - overlap)
    # grid
    Y = Na_C.reshape(-1, 1) * np.ones((1, len(Nr_C)))
    X = np.ones((len(Na_C), 1)) * Nr_C

    Xq, Yq = np.meshgrid(np.arange(0, Nr), np.arange(0, Na))
    
    # K_htc
    k_htc_interp = np.zeros((N, Na, Nr))
    for n in range(0, N):
        k_htc_filt = scipy.ndimage.gaussian_filter(k_htc[n, :, :], sigma=sigma_slope, truncate=truncate_slope, mode='nearest')
        f = interp2d(Y, X, k_htc_filt, kind='cubic')
        k_htc_interp[n, :, :] = f(Yq, Xq)
        #f_interp = interp2d(X.flatten(), Y.flatten(), k_htc_filt, kind='spline')
        #k_htc_interp[n, :, :] = f_interp(Xq, Yq).reshape(Xq.shape) #TODO 
    
    return k_htc_interp

def intercept_filtering(dem, ts_data, inps, k_htc_interp, meta):
    """Filtering parameters for obtaining high-frequency texture correlation"""
    sigma_intercept = 251 
    w_intercept = 251 
    truncate_intercept = ((w_intercept - 1)/2 - 0.5)/sigma_intercept
    # TODO
    #ref_y = int(meta['REF_Y'])
    #ref_x = int(meta['REF_X'])
    ref_y = 282
    ref_x = 204
    N, Na, Nr = ts_data.shape
    lamda = 0.3284035

    ts_data = 4 * np.pi / lamda * ts_data[:N, :, :]
    reference_value = ts_data[:, ref_y-1, ref_x-1]
    ts_data = ts_data - reference_value[:, np.newaxis, np.newaxis]
    ts_data[np.isnan(ts_data)] = 0

    phase_ts_htc_low = ts_data
    
    intercept = np.zeros(phase_ts_htc_low.shape)
    for n in range(0, N):
        tmp = ts_data[n, :, :] - k_htc_interp[n, :, :] * dem
        tmp_filt = scipy.ndimage.gaussian_filter(tmp, sigma=w_intercept, truncate=truncate_intercept, mode='nearest')
        tmp = tmp - tmp_filt
        phase_ts_htc_low[n, :, :] = tmp
        intercept[n, :, :] = tmp_filt
    reference_phase = phase_ts_htc_low[:, ref_y - 1, ref_x - 1]   
    phase_ts_htc_low = phase_ts_htc_low - reference_phase[:, np.newaxis, np.newaxis] #TODO

    return phase_ts_htc_low
    
    

def create_deformation_mask(dem, rate, maskdef):
    if maskdef == 1:
        mask_def = rate
        mask_def[np.abs(rate)>1.5] = np.nan
        mask_def[~np.isnan(mask_def)] = 1
    else:
        mask_def = np.ones(dem.shape) 
    return mask_def

############################################################################

def run_tropo_htc(inps):

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
    k_htc = estimate_local_slope(dem, ts_data, inps, n_ref, ts_obj.metadata)[2]
    # slope interpolation
    k_htc_interp = slope_interpolation(ts_data, inps, k_htc)
    # intercept filtering
    ts_htc_low = intercept_filtering(dem, ts_data, inps, k_htc_interp, ts_obj.metadata)
    lamda = 0.3284035 #TODO
    ts_htc_data = lamda / 4 /np.pi * ts_htc_low
    # write corrected time-series file
    meta = dict(ts_obj.metadata)
    if not inps.outfile:
        fbase = os.path.splitext(inps.timeseries_file)[0]
        inps.outfile = f'{fbase}_trophtc.h5'
    
    writefile.write(
        ts_htc_data,
        out_file=inps.outfile,
        metadata=meta,
        ref_file=inps.timeseries_file
    )

    return

