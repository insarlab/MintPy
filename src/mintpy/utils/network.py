"""Utilities for interferogram network selection."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################
# Recommend import:
#   from mintpy.utils import network as pnet


import itertools
import os
import sys

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from scipy import sparse

from mintpy.objects import ifgramStack, sensor
from mintpy.utils import ptime, readfile

SPEED_OF_LIGHT = 299792458   # m/s, speed of light


##################################################################
BASELINE_LIST_FILE = """
# Date  Bperp    dop0/PRF  dop1/PRF   dop2/PRF   PRF    slcDir
070106     0.0   0.03      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070106/
070709  2631.9   0.07      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070709/
070824  2787.3   0.07      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070824/
...
"""

IFGRAM_LIST_FILE = """
060713-070113
060828-070113
060828-070831
...
"""


##################################################################
def read_pairs_list(date12ListFile, dateList=[]):
    """Read Pairs List file like below:
    070311-070426
    070311-070611
    ...
    """
    # Read date12 list file
    date12List = sorted(list(np.loadtxt(date12ListFile, dtype=bytes).astype(str)))

    # Get dateList from date12List
    if not dateList:
        dateList = []
        for date12 in date12List:
            dates = date12.split('-')
            if not dates[0] in dateList:
                dateList.append(dates[0])
            if not dates[1] in dateList:
                dateList.append(dates[1])
        dateList.sort()
    date6List = ptime.yymmdd(dateList)

    # Get pair index
    pairs_idx = []
    for date12 in date12List:
        dates = date12.split('-')
        pair_idx = [date6List.index(dates[0]), date6List.index(dates[1])]
        pairs_idx.append(pair_idx)

    return pairs_idx


def write_pairs_list(pairs, dateList, outName):
    """Write pairs list file."""
    dateList6 = ptime.yymmdd(dateList)
    fl = open(outName, 'w')
    for idx in pairs:
        date12 = dateList6[idx[0]]+'-'+dateList6[idx[1]]+'\n'
        fl.write(date12)
    fl.close()
    return 1


def read_baseline_file(baselineFile, exDateList=[]):
    """Read bl_list.txt without dates listed in exDateList
    # Date  Bperp    dop0/PRF  dop1/PRF   dop2/PRF      PRF    slcDir
    070106     0.0   0.03      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070106/
    070709  2631.9   0.07      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070709/
    070824  2787.3   0.07      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070824/
    ...

    Examples:
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile)
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile,['080520','100726'])
        date8List, perpBaseList = read_baseline_file(baselineFile)[0:2]
    """
    exDateList = ptime.yymmdd(exDateList)
    if not exDateList:
        exDateList = []

    # Read baseline file into lines
    fb = open(baselineFile)
    lines = []
    for line in fb:
        line = str.replace(line, '\n', '').strip()
        lines.append(line)
    fb.close()

    # Read each line and put the values into arrays
    date6List = []
    perpBaseList = []
    dopplerList = []
    slcDirList = []
    for line in lines:
        c = line.split()    # splits on white space
        date = c[0]
        if date not in exDateList:
            date6List.append(date)
            perpBaseList.append(float(c[1]))
            try:
                dop = np.array([float(c[2]), float(c[3]), float(c[4])])
                prf = float(c[5])
                dop *= prf
                dopplerList.append(dop)
            except:
                pass
            try:
                slcDirList.append(c[6])
            except:
                pass

    date8List = ptime.yyyymmdd(date6List)
    return date8List, perpBaseList, dopplerList, slcDirList


def date12_list2index(date12_list, date_list=[]):
    """Convert list of date12 string into list of index"""
    # Get dateList from date12List
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = list(set(m_dates + s_dates))
    date6_list = ptime.yymmdd(sorted(ptime.yyyymmdd(date_list)))

    # Get pair index
    pairs_idx = []
    for date12 in date12_list:
        dates = date12.split('-')
        pair_idx = [date6_list.index(dates[0]), date6_list.index(dates[1])]
        pairs_idx.append(pair_idx)

    return pairs_idx


def get_date12_list(fname, dropIfgram=False):
    """Read date12 info from input file: Pairs.list or multi-group hdf5 file
    Parameters: fname       - string, path/name of input multi-group hdf5 file or text file
                dropIfgram  - bool, check the "dropIfgram" dataset in ifgramStack hdf5 file
    Returns:    date12_list - list of string in YYYYMMDD_YYYYMMDD format
    Example:
        date12List = get_date12_list('ifgramStack.h5')
        date12List = get_date12_list('ifgramStack.h5', dropIfgram=True)
        date12List = get_date12_list('coherenceSpatialAvg.txt')
    """
    date12_list = []
    ext = os.path.splitext(fname)[1].lower()
    if ext == '.h5':
        ftype = readfile.read_attribute(fname)['FILE_TYPE']
        if ftype == 'ifgramStack':
            date12_list = ifgramStack(fname).get_date12_list(dropIfgram=dropIfgram)
        else:
            return None
    else:
        date12_list = np.loadtxt(fname, dtype=bytes, usecols=0).astype(str).tolist()
        # for txt file with only one interferogram
        if isinstance(date12_list, str):
            date12_list = [date12_list]

    date12_list = sorted(date12_list)
    date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def critical_perp_baseline(sensor_name, inc_angle, print_msg=False):
    """Critical Perpendicular Baseline for each satellite"""
    # Jers: 5.712e3 m (near_range=688849.0551m)
    # Alos: 6.331e3 m (near_range=842663.2917m)
    # Tsx : 8.053e3 m (near_range=634509.1271m)
    sensor_dict = sensor.SENSOR_DICT[sensor_name.lower()]
    wvl = SPEED_OF_LIGHT / sensor_dict['carrier_frequency']
    near_range = 688849  # Yunjun 5/2016, case for Jers, need a automatic way to get this number
    rg_bandwidth = sensor_dict['chirp_bandwidth']
    bperp_c = wvl * (rg_bandwidth / SPEED_OF_LIGHT) * near_range * np.tan(inc_angle * np.pi / 180.0)
    if print_msg:
        print(f'Critical Perpendicular Baseline: {bperp_c} m')
    return bperp_c


def calculate_doppler_overlap(dop_a, dop_b, bandwidth_az):
    """Calculate Overlap Percentage of Doppler frequency in azimuth direction
    Inputs:
        dop_a/b      : np.array of 3 floats, doppler frequency
        bandwidth_az : float, azimuth bandwidth
    Output:
        dop_overlap  : float, doppler frequency overlap between a & b.
    """
    # Calculate mean Doppler difference between a and b
    no_of_rangepix = 5000
    ddiff = np.zeros(10)
    for i in range(10):
        rangepix = (i-1)*no_of_rangepix/10 + 1
        da = dop_a[0]+(rangepix-1)*dop_a[1]+(rangepix-1)**2*dop_a[2]
        db = dop_b[0]+(rangepix-1)*dop_b[1]+(rangepix-1)**2*dop_b[2]
        ddiff[i] = abs(da - db)
    ddiff_mean = np.mean(ddiff)

    #dopOverlap_prf = bandwidth_az - ddiff
    #dopOverlap_percent = np.mean(dopOverlap_prf / bandwidth_az * 100)
    ##overlapLessthanZero = np.nonzero(dopOverlap_percent < 0)
    #dopOverlap_percent = np.mean(dopOverlap_percent)
    # return dopOverlap_percent

    # Doppler overlap
    dop_overlap = (bandwidth_az - ddiff_mean) / bandwidth_az
    return dop_overlap


def simulate_coherence_v2(date12_list, decor_time=200.0, coh_resid=0.2, inc_angle=40, sensor_name='Sen',
                          display=False):
    """Simulate coherence version 2 (without using bl_list.txt file).
    Parameters: date12_list - list of string in YYYYMMDD_YYYYMMDD format, indicating pairs configuration
                decor_time  - float, decorrelation rate in days, time for coherence to drop to 1/e of its initial value
                coh_resid   - float, long-term coherence, minimum attainable coherence value
                inc_angle   - float, incidence angle in degrees
                sensor_name - string, SAR sensor name
                display     - bool, display result as matrix or not
    Returns:    coh         - 2D np.array in size of (ifgram_num)
    """
    num_pair = len(date12_list)
    date1s = [x.split('_')[0] for x in date12_list]
    date2s = [x.split('_')[1] for x in date12_list]
    date_list = sorted(list(set(date1s + date2s)))
    tbase_list = ptime.date_list2tbase(date_list)[0]

    SNR = 22  # NESZ = -22 dB from Table 1 in https://sentinels.copernicus.eu/web/sentinel/
    coh_thermal = 1. / (1. + 1./SNR)

    # bperp
    rng = np.random.default_rng(2)
    pbase_list = rng.normal(0, 50, num_pair).tolist()
    pbase_c = critical_perp_baseline(sensor_name, inc_angle)

    coh = np.zeros(num_pair, dtype=np.float32)
    for i in range(num_pair):
        date1, date2 = date12_list[i].split('_')
        ind1, ind2 = date_list.index(date1), date_list.index(date2)
        tbase = tbase_list[ind2] - tbase_list[ind1]
        pbase = pbase_list[ind2] - pbase_list[ind1]

        coh_geom = (pbase_c - abs(pbase)) / pbase_c
        coh_temp = np.multiply((coh_thermal - coh_resid), np.exp(-1*abs(tbase)/decor_time)) + coh_resid
        coh[i] = coh_geom * coh_temp

    if display:
        print('critical perp baseline: %.f m' % pbase_c)
        coh_mat = coherence_matrix(date12_list, coh)
        plt.figure()
        plt.imshow(coh_mat, vmin=0.0, vmax=1.0, cmap='jet')
        plt.xlabel('Image number')
        plt.ylabel('Image number')
        cbar = plt.colorbar()
        cbar.set_label('Coherence')
        plt.title('Coherence matrix')
        plt.show()

    return coh


def simulate_coherence(date12_list, baseline_file='bl_list.txt', sensor_name='Env', inc_angle=22.8,
                       decor_time=200.0, coh_resid=0.2, display=False):
    """Simulate coherence for a given set of interferograms
    Inputs:
        date12_list  - list of string in YYMMDD-YYMMDD format, indicating pairs configuration
        baseline_file - string, path of baseline list text file
        sensor_name     - string, SAR sensor name
        inc_angle  - float, incidence angle
        decor_time - float / 2D np.array in size of (1, pixel_num)
                     decorrelation rate in days, time for coherence to drop to 1/e of its initial value
        coh_resid  - float / 2D np.array in size of (1, pixel_num)
                     long-term coherence, minimum attainable coherence value
        display    - bool, display result as matrix or not
    Output:
        cohs       - 2D np.array in size of (ifgram_num, pixel_num)
    Example:
        date12_list = pnet.get_date12_list('ifgram_list.txt')
        cohs = simulate_coherences(date12_list, 'bl_list.txt', sensor_name='Tsx')

    References:
        Guarnieri, A. M. (2013), Introduction to RADAR, Politecnico di Milano Dipartimento di Elettronica
            e Informazione, Milano.
        Zebker, H. A., & Villasenor, J. (1992). Decorrelation in interferometric radar echoes.
            IEEE-TGRS, 30(5), 950-959.
        Hanssen, R. F. (2001). Radar interferometry: data interpretation and error analysis
            (Vol. 2). Dordrecht, Netherlands: Kluwer Academic Pub.
        Morishita, Y., & Hanssen, R. F. (2015). Temporal decorrelation in L-, C-, and X-band satellite
            radar interferometry for pasture on drained peat soils. IEEE-TGRS, 53(2), 1096-1104.
        Parizzi, A., Cong, X., & Eineder, M. (2009). First Results from Multifrequency Interferometry.
            A comparison of different decorrelation time constants at L, C, and X Band. ESA Scientific
            Publications(SP-677), 1-5.
    """
    date_list, pbase_list, dop_list = read_baseline_file(baseline_file)[0:3]
    tbase_list = ptime.date_list2tbase(date_list)[0]

    # Thermal decorrelation (Zebker and Villasenor, 1992, Eq.4)
    SNR = 19.5  # hardwired for Envisat (Guarnieri, 2013)
    coh_thermal = 1. / (1. + 1./SNR)

    pbase_c = critical_perp_baseline(sensor_name, inc_angle)
    bandwidth_az = sensor.SENSOR_DICT[sensor_name.lower()]['doppler_bandwidth']

    date12_list = ptime.yyyymmdd_date12(date12_list)
    ifgram_num = len(date12_list)

    if isinstance(decor_time, (int, np.int16, np.int32, float, np.float32, np.float64)):
        pixel_num = 1
        decor_time = float(decor_time)
    else:
        pixel_num = decor_time.shape[1]
    if decor_time == 0.:
        decor_time = 0.01
    cohs = np.zeros((ifgram_num, pixel_num), np.float32)
    for i in range(ifgram_num):
        if display:
            sys.stdout.write('\rinterferogram = %4d/%4d' % (i, ifgram_num))
            sys.stdout.flush()
        m_date, s_date = date12_list[i].split('_')
        m_idx = date_list.index(m_date)
        s_idx = date_list.index(s_date)

        pbase = pbase_list[s_idx] - pbase_list[m_idx]
        tbase = tbase_list[s_idx] - tbase_list[m_idx]

        # Geometric decorrelation (Hanssen, 2001, Eq. 4.4.12)
        coh_geom = (pbase_c - abs(pbase)) / pbase_c
        if coh_geom < 0.:
            coh_geom = 0.

        # Doppler centroid decorrelation (Hanssen, 2001, Eq. 4.4.13)
        if not dop_list:
            coh_dc = 1.
        else:
            coh_dc = calculate_doppler_overlap(dop_list[m_idx],
                                               dop_list[s_idx],
                                               bandwidth_az)
            if coh_dc < 0.:
                coh_dc = 0.

        # Option 1: Temporal decorrelation - exponential delay model (Parizzi et al., 2009; Morishita and Hanssen, 2015)
        coh_temp = np.multiply((coh_thermal - coh_resid), np.exp(-1*abs(tbase)/decor_time)) + coh_resid

        coh = coh_geom * coh_dc * coh_temp
        cohs[i, :] = coh
    #epsilon = 1e-3
    #cohs[cohs < epsilon] = epsilon
    if display:
        print('')

    if display:
        print('critical perp baseline: %.f m' % pbase_c)
        cohs_mat = coherence_matrix(date12_list, cohs)
        plt.figure()
        plt.imshow(cohs_mat, vmin=0.0, vmax=1.0, cmap='jet')
        plt.xlabel('Image number')
        plt.ylabel('Image number')
        cbar = plt.colorbar()
        cbar.set_label('Coherence')
        plt.title('Coherence matrix')
        plt.show()
    return cohs


##################################################################
def threshold_doppler_overlap(date12_list, date_list, dop_list, bandwidth_az, dop_overlap_min=0.15):
    """Remove pairs/interoferogram with doppler overlap larger than critical value
    Inputs:
        date12_list : list of string, for date12 in YYMMDD-YYMMDD format
        date_list   : list of string, for date in YYMMDD/YYYYMMDD format, optional
        dop_list    : list of list of 3 float, for centroid Doppler frequency
        bandwidth_az    : float, bandwidth in azimuth direction
        dop_overlap_min : float, minimum overlap of azimuth Doppler frequency
    Outputs:
        date12_list : list of string, for date12 in YYMMDD-YYMMDD format
    """
    if not date12_list:
        return []
    # Get date6_list
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))
    date6_list = ptime.yymmdd(date_list)

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        dop_overlap = calculate_doppler_overlap(dop_list[idx1],
                                                dop_list[idx2],
                                                bandwidth_az)
        if dop_overlap >= dop_overlap_min:
            date12_list_out.append(date12)
    return date12_list_out


def threshold_perp_baseline(date12_list, date_list, pbase_list, pbase_max, pbase_min=0.0):
    """Remove pairs/interoferogram out of [pbase_min, pbase_max]
    Inputs:
        date12_list : list of string for date12 in YYMMDD-YYMMDD format
        date_list   : list of string for date in YYMMDD/YYYYMMDD format, optional
        pbase_list  : list of float for perpendicular spatial baseline
        pbase_max   : float, maximum perpendicular baseline
        pbase_min   : float, minimum perpendicular baseline
    Output:
        date12_list_out : list of string for date12 in YYMMDD-YYMMDD format
    Example:
        date12_list = threshold_perp_baseline(date12_list, date_list, pbase_list, 500)
    """
    if not date12_list:
        return []

    # Get date_list
    delimiter = [i for i in ['-', '_'] if i in date12_list[0]][0]
    if not date_list:
        m_dates = [date12.split(delimiter)[0] for date12 in date12_list]
        s_dates = [date12.split(delimiter)[1] for date12 in date12_list]
        date_list = sorted(list(set(m_dates + s_dates)))

    if not len(date_list) == len(pbase_list):
        print('ERROR: number of existing dates is not equal to number of perp baseline!')
        print('date list is needed for threshold filtering!')
        print('skip filtering.')
        return date12_list

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split(delimiter)
        pbase1 = pbase_list[date_list.index(date1)]
        pbase2 = pbase_list[date_list.index(date2)]
        pbase = abs(pbase2 - pbase1)

        if pbase_min <= pbase <= pbase_max:
            date12_list_out.append(date12)

    return date12_list_out


def threshold_temporal_baseline(date12_list, btemp_max, keep_seasonal=True, btemp_min=0.0):
    """Remove pairs/interferograms out of min/max/seasonal temporal baseline limits
    Inputs:
        date12_list : list of string for date12 in YYMMDD-YYMMDD or YYYYMMDD_YYYYMMDD format
        btemp_max   : float, maximum temporal baseline
        btemp_min   : float, minimum temporal baseline
        keep_seasonal : keep interferograms with seasonal temporal baseline
    Output:
        date12_list_out : list of string for date12
    Example:
        date12_list = threshold_temporal_baseline(date12_list, 200)
        date12_list = threshold_temporal_baseline(date12_list, 200, False)
    """
    if not date12_list:
        return []

    # Get date list and tbase list
    delimiter = [i for i in ['-', '_'] if i in date12_list[0]][0]
    m_dates = [date12.split(delimiter)[0] for date12 in date12_list]
    s_dates = [date12.split(delimiter)[1] for date12 in date12_list]
    date_list = sorted(list(set(m_dates + s_dates)))
    tbase_list = ptime.date_list2tbase(date_list)[0]

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split(delimiter)
        tbase1 = tbase_list[date_list.index(date1)]
        tbase2 = tbase_list[date_list.index(date2)]
        tbase = int(abs(tbase2 - tbase1))

        if btemp_min <= tbase <= btemp_max:
            date12_list_out.append(date12)

        elif keep_seasonal and tbase/30 in [11, 12]:
            date12_list_out.append(date12)

    return date12_list_out


def coherence_matrix(date12_list, coh_list, diag_value=np.nan, fill_triangle='both', date_list=None):
    """Return coherence matrix based on input date12 list and its coherence
    Inputs:
        date12_list - list of string in YYMMDD-YYMMDD format
        coh_list    - list of float, average coherence for each interferograms
        diag_value  - number, value to be filled in the diagonal
        fill_triangle - str, 'both', 'upper', 'lower'
    Output:
        coh_matrix  - 2D np.array with dimension length = date num
                      np.nan value for interferograms non-existed.
                      1.0 for diagonal elements
    """
    # Get date list
    date12_list = ptime.yymmdd_date12(date12_list)
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = sorted(list(set(m_dates + s_dates)))
    date_list = ptime.yymmdd(date_list)
    date_num = len(date_list)

    coh_mat = np.zeros([date_num, date_num])
    coh_mat[:] = np.nan
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date_list.index(date1)
        idx2 = date_list.index(date2)
        coh = coh_list[date12_list.index(date12)]
        if fill_triangle in ['upper', 'both']:
            coh_mat[idx1, idx2] = coh  # symmetric
        if fill_triangle in ['lower', 'both']:
            coh_mat[idx2, idx1] = coh

    if diag_value is not np.nan:
        for i in range(date_num):    # diagonal value
            coh_mat[i, i] = diag_value
    return coh_mat


def threshold_coherence_based_mst(date12_list, coh_list):
    """Return a minimum spanning tree of network based on the coherence inverse.
    Inputs:
        date12_list - list of string in YYMMDD-YYMMDD format
        coh_list    - list of float, average coherence for each interferogram
    Output:
        mst_date12_list - list of string in YYMMDD-YYMMDD format, for MST network of interferograms
    """
    # coh_list --> coh_mat --> weight_mat
    coh_mat = coherence_matrix(date12_list, coh_list)
    mask = ~np.isnan(coh_mat)
    wei_mat = np.zeros(coh_mat.shape)
    wei_mat[:] = np.inf
    wei_mat[mask] = 1/coh_mat[mask]

    # MST path based on weight matrix
    wei_mat_csr = sparse.csr_matrix(wei_mat)
    mst_mat_csr = sparse.csgraph.minimum_spanning_tree(wei_mat_csr)

    # Get date6_list
    date12_list = ptime.yymmdd_date12(date12_list)
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date6_list = ptime.yymmdd(sorted(ptime.yyyymmdd(list(set(m_dates + s_dates)))))

    # Convert MST index matrix into date12 list
    [s_idx_list, m_idx_list] = [date_idx_array.tolist()
                                for date_idx_array in sparse.find(mst_mat_csr)[0:2]]
    mst_date12_list = []
    for m_idx, s_idx in zip(m_idx_list, s_idx_list):
        idx = sorted([m_idx, s_idx])
        date12 = date6_list[idx[0]]+'-'+date6_list[idx[1]]
        mst_date12_list.append(date12)
    return mst_date12_list


def pair_sort(pairs):
    for i, pair in enumerate(pairs):
        if pair[0] > pair[1]:
            pairs[i][0] = pair[1]
            pairs[i][1] = pair[0]
    pairs = sorted(pairs)
    return pairs


def pair_merge(pairs1, pairs2):
    pairs = pairs1
    for pair in pairs2:
        if pair not in pairs:
            pairs.append(pair)

    pairs = sorted(pairs)
    return pairs


def select_pairs_all(date_list, date_format='YYMMDD'):
    """Select All Possible Pairs/Interferograms
    Input : date_list   - list of date in YYMMDD/YYYYMMDD format
    Output: date12_list - list date12 in YYMMDD-YYMMDD format
    Reference:
        Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
        based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
    """
    date8_list = sorted(ptime.yyyymmdd(date_list))
    date6_list = ptime.yymmdd(date8_list)
    date12_list = list(itertools.combinations(date6_list, 2))
    date12_list = [date12[0]+'-'+date12[1] for date12 in date12_list]
    if date_format == 'YYYYMMDD':
        date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def select_pairs_sequential(date_list, num_conn=2, date_format=None):
    """Select Pairs in a Sequential way:
        For each acquisition, find its num_connection nearest acquisitions in the past time.

    Parameters: date_list   - list of str for date
                num_conn    - int, number of sequential connections
                date_format - str / None, output date format
    Returns:    date12_list - list of str for date12
    """

    date_list = sorted(date_list)
    date_inds = list(range(len(date_list)))

    # Get pairs index list
    date12_inds = []
    for date_ind in date_inds:
        for i in range(num_conn):
            if date_ind-i-1 >= 0:
                date12_inds.append([date_ind-i-1, date_ind])
    date12_inds = [sorted(i) for i in sorted(date12_inds)]

    # Convert index into date12
    date12_list = [f'{date_list[ind12[0]]}_{date_list[ind12[1]]}'
                  for ind12 in date12_inds]

    # adjust output date format
    if date_format is not None:
        if date_format == 'YYYYMMDD':
            date12_list = ptime.yyyymmdd_date12(date12_list)
        elif date_format == 'YYMMDD':
            date12_list = ptime.yymmdd_date12(date12_list)
        else:
            raise ValueError(f'un-supported date format: {date_format}!')

    return date12_list


def select_pairs_hierarchical(date_list, pbase_list, temp_perp_list, date_format='YYMMDD'):
    """Select Pairs in a hierarchical way using list of temporal and perpendicular baseline thresholds
        For each temporal/perpendicular combination, select all possible pairs; and then merge all combination results
        together for the final output (Zhao, 2015).
    Inputs:
        date_list  : list of date in YYMMDD/YYYYMMDD format
        pbase_list : list of float, perpendicular spatial baseline
        temp_perp_list : list of list of 2 floats, for list of temporal/perp baseline, e.g.
                         [[32.0, 800.0], [48.0, 600.0], [64.0, 200.0]]
    Examples:
        pairs = select_pairs_hierarchical(date_list, pbase_list, [[32.0, 800.0], [48.0, 600.0], [64.0, 200.0]])
    Reference:
        Zhao, W., (2015), Small deformation detected from InSAR time-series and their applications in geophysics, Doctoral
        dissertation, Univ. of Miami, Section 6.3.
    """
    # Get all date12
    date12_list_all = select_pairs_all(date_list)

    # Loop of Threshold
    print('List of temporal and perpendicular spatial baseline thresholds:')
    print(temp_perp_list)
    date12_list = []
    for temp_perp in temp_perp_list:
        tbase_max = temp_perp[0]
        pbase_max = temp_perp[1]
        date12_list_tmp = threshold_temporal_baseline(date12_list_all,
                                                      tbase_max,
                                                      keep_seasonal=False)
        date12_list_tmp = threshold_perp_baseline(date12_list_tmp,
                                                  date_list,
                                                  pbase_list,
                                                  pbase_max)
        date12_list += date12_list_tmp
    date12_list = sorted(list(set(date12_list)))
    if date_format == 'YYYYMMDD':
        date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def select_pairs_delaunay(date_list, pbase_list, norm=True, date_format='YYMMDD'):
    """Select Pairs using Delaunay Triangulation based on temporal/perpendicular baselines
    Inputs:
        date_list  : list of date in YYMMDD/YYYYMMDD format
        pbase_list : list of float, perpendicular spatial baseline
        norm       : normalize temporal baseline to perpendicular baseline
    Key points
        1. Define a ratio between perpendicular and temporal baseline axis units (Pepe and Lanari, 2006, TGRS).
        2. Pairs with too large perpendicular / temporal baseline or Doppler centroid difference should be removed
           after this, using a threshold, to avoid strong decorrelations (Zebker and Villasenor, 1992, TGRS).
    Reference:
        Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
        of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
        Zebker, H. A., and J. Villasenor (1992), Decorrelation in interferometric radar echoes, IEEE TGRS, 30(5), 950-959.
    """
    # Get temporal baseline in days
    date6_list = ptime.yymmdd(date_list)
    date8_list = ptime.yyyymmdd(date_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]

    # Normalization (Pepe and Lanari, 2006, TGRS)
    if norm:
        temp2perp_scale = (max(pbase_list)-min(pbase_list)) / (max(tbase_list)-min(tbase_list))
        tbase_list = [tbase*temp2perp_scale for tbase in tbase_list]

    # Generate Delaunay Triangulation
    date12_idx_list = Triangulation(tbase_list, pbase_list).edges.tolist()
    date12_idx_list = [sorted(idx) for idx in sorted(date12_idx_list)]

    # Convert index into date12
    date12_list = [date6_list[idx[0]]+'-'+date6_list[idx[1]]
                   for idx in date12_idx_list]
    if date_format == 'YYYYMMDD_YYYYMMDD':
        date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def select_pairs_mst(date_list, pbase_list, date_format='YYMMDD'):
    """Select Pairs using Minimum Spanning Tree technique
        Connection Cost is calculated using the baseline distance in perp and scaled temporal baseline (Pepe and Lanari,
        2006, TGRS) plane.
    Inputs:
        date_list  : list of date in YYMMDD/YYYYMMDD format
        pbase_list : list of float, perpendicular spatial baseline
    References:
        Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
        of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
        Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
    """
    # Get temporal baseline in days
    date6_list = ptime.yymmdd(date_list)
    date8_list = ptime.yyyymmdd(date_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    # Normalization (Pepe and Lanari, 2006, TGRS)
    temp2perp_scale = (max(pbase_list)-min(pbase_list)) / (max(tbase_list)-min(tbase_list))
    tbase_list = [tbase*temp2perp_scale for tbase in tbase_list]

    # Get weight matrix
    ttMat1, ttMat2 = np.meshgrid(np.array(tbase_list), np.array(tbase_list))
    ppMat1, ppMat2 = np.meshgrid(np.array(pbase_list), np.array(pbase_list))
    ttMat = np.abs(ttMat1 - ttMat2)  # temporal distance matrix
    ppMat = np.abs(ppMat1 - ppMat2)  # spatial distance matrix

    # 2D distance matrix in temp/perp domain
    weightMat = np.sqrt(np.square(ttMat) + np.square(ppMat))
    weightMat = sparse.csr_matrix(weightMat)  # compress sparse row matrix

    # MST path based on weight matrix
    mstMat = sparse.csgraph.minimum_spanning_tree(weightMat)

    # Convert MST index matrix into date12 list
    [s_idx_list, m_idx_list] = [date_idx_array.tolist()
                                for date_idx_array in sparse.find(mstMat)[0:2]]
    date12_list = []
    for i in range(len(m_idx_list)):
        idx = sorted([m_idx_list[i], s_idx_list[i]])
        date12 = date6_list[idx[0]]+'-'+date6_list[idx[1]]
        date12_list.append(date12)
    if date_format == 'YYYYMMDD':
        date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def select_pairs_star(date_list, m_date=None, pbase_list=[], date_format='YYMMDD'):
    """Select Star-like network/interferograms/pairs, it's a single reference network, similar to PS approach.
    Usage:
        m_date : reference date, choose it based on the following cretiria:
                 1) near the center in temporal and spatial baseline
                 2) prefer winter season than summer season for less temporal decorrelation
    Reference:
        Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
    """
    date8_list = sorted(ptime.yyyymmdd(date_list))
    date6_list = ptime.yymmdd(date8_list)

    # Select reference date if not chosen
    if not m_date:
        m_date = select_reference_date(date8_list, pbase_list)
        print('auto select reference date: '+m_date)

    # Check input reference date
    m_date8 = ptime.yyyymmdd(m_date)
    if m_date8 not in date8_list:
        print('Input reference date does not exist in date list!')
        print(f'Input reference date: {m_date8}')
        print(f'Input date list: {date8_list}')
        m_date8 = None

    # Generate star/ps network
    m_idx = date8_list.index(m_date8)
    date12_idx_list = [sorted([m_idx, s_idx]) for s_idx in range(len(date8_list))
                       if s_idx is not m_idx]
    date12_list = [date6_list[idx[0]]+'-'+date6_list[idx[1]]
                   for idx in date12_idx_list]
    if date_format == 'YYYYMMDD':
        date12_list = ptime.yyyymmdd_date12(date12_list)
    return date12_list


def select_reference_date(date_list, pbase_list=[]):
    """Select super reference date based on input temporal and/or perpendicular baseline info.
    Return reference date in YYYYMMDD format.
    """
    date8_list = ptime.yyyymmdd(date_list)
    if not pbase_list:
        # Choose date in the middle
        m_date8 = date8_list[int(len(date8_list)/2)]
    else:
        # Get temporal baseline list
        tbase_list = ptime.date_list2tbase(date8_list)[0]
        # Normalization (Pepe and Lanari, 2006, TGRS)
        temp2perp_scale = (max(pbase_list)-min(pbase_list)) / (max(tbase_list)-min(tbase_list))
        tbase_list = [tbase*temp2perp_scale for tbase in tbase_list]
        # Get distance matrix
        ttMat1, ttMat2 = np.meshgrid(np.array(tbase_list),
                                     np.array(tbase_list))
        ppMat1, ppMat2 = np.meshgrid(np.array(pbase_list),
                                     np.array(pbase_list))
        ttMat = np.abs(ttMat1 - ttMat2)  # temporal distance matrix
        ppMat = np.abs(ppMat1 - ppMat2)  # spatial distance matrix
        # 2D distance matrix in temp/perp domain
        disMat = np.sqrt(np.square(ttMat) + np.square(ppMat))

        # Choose date minimize the total distance of temp/perp baseline
        disMean = np.mean(disMat, 0)
        m_idx = np.argmin(disMean)
        m_date8 = date8_list[m_idx]
    return m_date8


def select_reference_interferogram(date12_list, date_list, pbase_list, m_date=None):
    """Select reference interferogram based on input temp/perp baseline info
    If m_date is specified, select its closest s_date, which is newer than m_date;
        otherwise, choose the closest pair among all pairs as reference interferogram.
    Example:
        m_date12   = pnet.select_reference_ifgram(date12_list, date_list, pbase_list)
        '080211-080326' = pnet.select_reference_ifgram(date12_list, date_list, pbase_list, m_date='080211')
    """
    pbase_array = np.array(pbase_list, dtype='float64')
    # Get temporal baseline
    date8_list = ptime.yyyymmdd(date_list)
    date6_list = ptime.yymmdd(date8_list)
    tbase_array = np.array(ptime.date_list2tbase(date8_list)[0], dtype='float64')
    # Normalization (Pepe and Lanari, 2006, TGRS)
    temp2perp_scale = (max(pbase_array)-min(pbase_array)) / (max(tbase_array)-min(tbase_array))
    tbase_array *= temp2perp_scale

    # Calculate sqrt of temp/perp baseline for input pairs
    idx1 = np.array([date6_list.index(date12.split('-')[0]) for date12 in date12_list])
    idx2 = np.array([date6_list.index(date12.split('-')[1]) for date12 in date12_list])
    base_distance = np.sqrt((tbase_array[idx2] - tbase_array[idx1])**2 +
                            (pbase_array[idx2] - pbase_array[idx1])**2)

    # Get reference interferogram index
    if not m_date:
        # Choose pair with shortest temp/perp baseline
        m_date12_idx = np.argmin(base_distance)
    else:
        m_date = ptime.yymmdd(m_date)
        # Choose pair contains m_date with shortest temp/perp baseline
        m_date12_idx_array = np.array([date12_list.index(date12) for date12 in date12_list
                                       if m_date+'-' in date12])
        min_base_distance = np.min(base_distance[m_date12_idx_array])
        m_date12_idx = np.where(base_distance == min_base_distance)[0][0]

    m_date12 = date12_list[m_date12_idx]
    return m_date12
