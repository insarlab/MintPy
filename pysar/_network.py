#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2016, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# Recommended Usage:
#   import pysar._network as pnet
#


import os
import sys
import datetime
import itertools

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.sparse as sparse
from scipy.sparse.csgraph import minimum_spanning_tree

import pysar._datetime as ptime
import pysar._readfile as readfile
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


##################################################################
BASELINE_LIST_FILE='''
# Date  Bperp    dop0/PRF  dop1/PRF   dop2/PRF   PRF    slcDir
070106     0.0   0.03      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070106/
070709  2631.9   0.07      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070709/
070824  2787.3   0.07      0.0000000  0.000000   2155.2 /KyushuT422F650AlosA/SLC/070824/
...
'''

IFGRAM_LIST_FILE='''
060713-070113
060828-070113
060828-070831
...
'''


##################################################################
def read_pairs_list(date12ListFile, dateList=[]):
    '''Read Pairs List file like below:
    070311-070426
    070311-070611
    ...
    '''
    # Read date12 list file
    date12List = sorted(list(np.loadtxt(date12ListFile, dtype=str)))

    # Get dateList from date12List
    if not dateList:
        dateList = []
        for date12 in date12List:
            dates = date12.split('-')
            if not dates[0] in dateList: dateList.append(dates[0])
            if not dates[1] in dateList: dateList.append(dates[1])
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
    '''Write pairs list file.'''
    dateList6 = ptime.yymmdd(dateList)
    fl = open(outName,'w')
    for idx in pairs:
        date12 = dateList6[idx[0]]+'-'+dateList6[idx[1]]+'\n'
        fl.write(date12)
    fl.close()
    return 1


def read_igram_pairs(igramFile):
    '''Read pairs index from hdf5 file'''
    ## Read Igram file
    h5file = h5py.File(igramFile,'r')
    k = h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    if k[0] not in  ['interferograms','coherence','wrapped']:
        print 'Only interferograms / coherence / wrapped are supported.';  sys.exit(1)

    dateList  = ptime.ifgram_date_list(igramFile)
    dateList6 = ptime.yymmdd(dateList)

    pairs = []
    igramList=h5file[k[0]].keys()
    for igram in igramList:
        date12 = h5file[k[0]][igram].attrs['DATE12'].split('-')
        pairs.append([dateList6.index(date12[0]),dateList6.index(date12[1])])
    h5file.close()

    pairs = pair_sort(pairs)

    return pairs


def read_baseline_file(baselineFile, exDateList=[]):
    '''Read bl_list.txt without dates listed in exDateList
    # Date  Bperp    dop0/PRF  dop1/PRF   dop2/PRF      PRF    slcDir
    070106     0.0   0.03      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070106/
    070709  2631.9   0.07      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070709/
    070824  2787.3   0.07      0.0000000  0.00000000000 2155.2 /scratch/KyushuT422F650AlosA/SLC/070824/
    ...
    
    Examples:
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile)
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile,['080520','100726'])
        date8List, perpBaseList = read_baseline_file(baselineFile)[0:2]
    '''
    exDateList = ptime.yymmdd(exDateList)
    if not exDateList:  exDateList = []

    ## Read baseline file into lines
    fb = open(baselineFile)
    lines = []
    for line in fb.xreadlines():
        l = str.replace(line,'\n','').strip()
        lines.append(l)
    fb.close()

    ## Read each line and put the values into arrays
    date6List    = []
    perpBaseList = []
    dopplerList  = []
    slcDirList   = []
    for line in lines:
        c = line.split()    # splits on white space
        date = c[0]
        if not date in exDateList:
            date6List.append(date)
            perpBaseList.append(float(c[1]))
            dop = np.array([float(c[2]), float(c[3]), float(c[4])])
            prf = float(c[5])
            dop *= prf
            dopplerList.append(dop)
            slcDirList.append(c[6])

    date8List = ptime.yyyymmdd(date6List)
    return date8List, perpBaseList, dopplerList, slcDirList


def date12_list2index(date12_list, date_list=[]):
    '''Convert list of date12 string into list of index'''
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


def get_date12_list(File, check_drop_ifgram=False):
    '''Read Date12 info from input file: Pairs.list or multi-group hdf5 file
    Inputs:
        File - string, path/name of input multi-group hdf5 file or text file
        check_drop_ifgram - bool, check the "drop_ifgram" attribute or not for multi-group hdf5 file
    Output:
        date12_list - list of string in YYMMDD-YYMMDD format
    Example:
        date12List = get_date12_list('unwrapIfgram.h5')
        date12List = get_date12_list('unwrapIfgram.h5', check_drop_ifgram=True)
        date12List = get_date12_list('Pairs.list')
    '''
    #print 'read pairs info from '+File
    date12_list = []
    ext = os.path.splitext(File)[1].lower()
    if ext == '.h5':
        k = readfile.read_attribute(File)['FILE_TYPE']
        h5 = h5py.File(File, 'r')
        epochList = sorted(h5[k].keys())
        for epoch in epochList:
            if not check_drop_ifgram or h5[k][epoch].attrs['drop_ifgram'] == 'no':
                date12 = h5[k][epoch].attrs['DATE12']
                date12_list.append(date12)
        h5.close()
    else:
        txtContent = np.loadtxt(File, dtype=str)
        if len(txtContent.shape) == 1:
            txtContent = txtContent.reshape(-1,1)
        date12_list = [i for i in txtContent[:,0]]

    date12_list = sorted(date12_list)
    return date12_list


def igram_perp_baseline_list(File):
    '''Get perpendicular baseline list from input multi_group hdf5 file'''
    print 'read perp baseline info from '+File
    k = readfile.read_attribute(File)['FILE_TYPE']
    h5 = h5py.File(File, 'r')
    epochList = sorted(h5[k].keys())
    p_baseline_list = []
    for epoch in epochList:
        p_baseline = (float(h5[k][epoch].attrs['P_BASELINE_BOTTOM_HDR'])+\
                      float(h5[k][epoch].attrs['P_BASELINE_TOP_HDR']))/2
        p_baseline_list.append(p_baseline)
    h5.close()
    return p_baseline_list


##################################################################
def azimuth_bandwidth(sensor):
    '''Find the hardwired azimuth bandwidth in hertz for the given satellite'''
    if    sensor == 'Ers'  :  bandwidth =  1300.0
    elif  sensor == 'Env'  :  bandwidth =  1340.0
    elif  sensor == 'S1'   :  bandwidth =  4000.0   # shong 08/2016 sould be checked
    elif  sensor == 'Rsat' :  bandwidth =   900.0
    elif  sensor == 'Rsat2':  bandwidth =   900.0

    elif  sensor == 'Jers' :  bandwidth =   900.0   # FA 11/2015 just copied, need to research
    elif  sensor == 'Alos' :  bandwidth =  1720.0

    elif  sensor == 'Tsx'  :  bandwidth = 15000.0
    elif  sensor == 'Csk'  :  bandwidth = 15000.0   # FA 9/2015  shoud be checked
    elif  sensor == 'Kmps5':  bandwidth = 15000.0   # shong 08/2016 sould be checked
    else: print 'satellite not found'; bandwidth = None
    return bandwidth


def range_bandwidth(sensor):
    ## Range Bandwidth in Hz for the given satellite
    if    sensor == 'Ers' :  bandwidth = 15.55e6
    elif  sensor == 'Env' :  bandwidth = 16.00e6

    elif  sensor == 'Jers':  bandwidth = 15e6      # Jers only has HH pol
    elif  sensor == 'Alos':  bandwidth = 14e6      # for FBD, 28MHz for FBS

    elif  sensor == 'Tsx' :  bandwidth = 150e6
    return bandwidth


def wavelength(sensor):
    if    sensor == 'Ers'  :  center_frequency = 5.300e9
    elif  sensor == 'Env'  :  center_frequency = 5.331e9
    elif  sensor == 'S1'   :  center_frequency = 5.405e9
    elif  sensor == 'Rsat' :  center_frequency = 5.300e9
    elif  sensor == 'Rsat2':  center_frequency = 5.405e9

    elif  sensor == 'Jers' :  center_frequency = 1.275e9
    elif  sensor == 'Alos' :  center_frequency = 1.270e9
    elif  sensor == 'Alos2':  center_frequency = 1.270e9

    elif  sensor == 'Tsx'  :  center_frequency = 9.65e9
    elif  sensor == 'Csk'  :  center_frequency = 9.60e9
    elif  sensor == 'Kmps5':  center_frequency = 9.66e9

    c = 299792458;   # m/s, speed of light
    wavelength = c / center_frequency
    return wavelength


def incidence_angle(sensor, inc_angle=None):
    if not inc_angle:
        if   sensor == 'Ers' :  inc_angle = 34.3
        elif sensor == 'Env' :  inc_angle = 34.3

        elif sensor == 'Jers':  inc_angle = 35.21
        elif sensor == 'Alos':  inc_angle = 34.3     # degree, for ALOS PALSAR Fine mode

        elif sensor == 'Tsx' :  inc_angle = 39.23    # Yunjun 5/2016, for TaizhouTsx, not sure it's for all cases.
    return inc_angle


def signal2noise_ratio(sensor):
    '''Fine the Signal to Noise Ratio in dB for the given satellite
    Reference:
        ERS - Zebker et al., 1994, TGRS
        Envisat - Guarnieri, A.M., 2013. Introduction to RADAR. POLIMI DEI, Milano.
        JERS - https://directory.eoportal.org/web/eoportal/satellite-missions/j/jers-1
    '''
    if   sensor.startswith('Ers') :  SNR = 11.7
    elif sensor.startswith('Env') :  SNR = 19.5
    elif sensor.startswith('Jers'):  SNR = 14
    else: print 'satellite not found'; SNR = None
    return SNR


def critical_perp_baseline(sensor, inc_angle=None, print_msg=False):
    '''Critical Perpendicular Baseline for each satellite'''
    # Jers: 5.712e3 m (near_range=688849.0551m)
    # Alos: 6.331e3 m (near_range=842663.2917m)
    # Tsx : 8.053e3 m (near_range=634509.1271m)

    c = 299792458;   # m/s, speed of light
    wvl = wavelength(sensor)
    near_range = 688849;        # Yunjun 5/2016, case for Jers, need a automatic way to get this number
    rg_bandwidth = range_bandwidth(sensor)
    inc_angle    = incidence_angle(sensor, inc_angle) / 180 * np.pi
    Bperp_c      = wvl * (rg_bandwidth/c) * near_range * np.tan(inc_angle)
    if print_msg:
        print 'Critical Perpendicular Baseline: '+str(Bperp_c)+' m'
    return Bperp_c


def calculate_doppler_overlap(dop_a, dop_b, bandwidth_az):
    '''Calculate Overlap Percentage of Doppler frequency in azimuth direction
    Inputs:
        dop_a/b      : np.array of 3 floats, doppler frequency
        bandwidth_az : float, azimuth bandwidth
    Output:
        dop_overlap  : float, doppler frequency overlap between a & b.
    '''
    # Calculate mean Doppler difference between a and b
    no_of_rangepix = 5000
    ddiff = np.zeros(10)
    for i in range (10):
        rangepix = (i-1)*no_of_rangepix/10 + 1
        da       = dop_a[0]+(rangepix-1)*dop_a[1]+(rangepix-1)**2*dop_a[2]
        db       = dop_b[0]+(rangepix-1)*dop_b[1]+(rangepix-1)**2*dop_b[2]
        ddiff[i] =  abs(da - db)
    ddiff_mean = np.mean(ddiff)

    #dopOverlap_prf = bandwidth_az - ddiff
    #dopOverlap_percent = np.mean(dopOverlap_prf / bandwidth_az * 100)
    ##overlapLessthanZero = np.nonzero(dopOverlap_percent < 0)
    #dopOverlap_percent = np.mean(dopOverlap_percent)
    #return dopOverlap_percent
    
    # Doppler overlap
    dop_overlap = (bandwidth_az - ddiff_mean) / bandwidth_az
    return dop_overlap


def simulate_coherence(date12_list, baselineFile='bl_list.txt', sensor='Env', inc_angle=22.8,\
                       decor_time=200.0, coh_resid=0.2, display=False):
    '''Simulate coherence for a given set of interferograms
    Inputs:
        date12_list  - list of string in YYMMDD-YYMMDD format, indicating pairs configuration
        baselineFile - string, path of baseline list text file
        sensor     - string, SAR sensor
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
        cohs = simulate_coherences(date12_list, 'bl_list.txt', sensor='Tsx')

    References:
        Zebker, H. A., & Villasenor, J. (1992). Decorrelation in interferometric radar echoes.
            IEEE-TGRS, 30(5), 950-959. 
        Hanssen, R. F. (2001). Radar interferometry: data interpretation and error analysis
            (Vol. 2). Dordrecht, Netherlands: Kluwer Academic Pub.
        Morishita, Y., & Hanssen, R. F. (2015). Temporal decorrelation in L-, C-, and X-band satellite
            radar interferometry for pasture on drained peat soils. IEEE-TGRS, 53(2), 1096-1104. 
        Parizzi, A., Cong, X., & Eineder, M. (2009). First Results from Multifrequency Interferometry.
            A comparison of different decorrelation time constants at L, C, and X Band. ESA Scientific
            Publications(SP-677), 1-5. 
    '''
    date8_list, pbase_list, dop_list = read_baseline_file(baselineFile)[0:3]
    date6_list = ptime.yymmdd(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]

    #Thermal decorrelation (Zebker and Villasenor, 1992, Eq.4)
    SNR = signal2noise_ratio(sensor)
    coh_thermal = 1. / (1. + 1./SNR)

    pbase_c = critical_perp_baseline(sensor, inc_angle)
    bandwidth_az = azimuth_bandwidth(sensor)

    ifgram_num = len(date12_list)
    if isinstance(decor_time, (int, float)):
        pixel_num = 1
        decor_time = float(decor_time)
    else:
        pixel_num = decor_time.shape[1]
    cohs = np.zeros((ifgram_num, pixel_num), np.float32)
    for i in range(ifgram_num):
        if display:
            sys.stdout.write('\rinterferogram = %4d/%4d' % (i, ifgram_num))
            sys.stdout.flush()
        m_date, s_date = date12_list[i].split('-')
        m_idx = date6_list.index(m_date)
        s_idx = date6_list.index(s_date)

        pbase = pbase_list[s_idx] - pbase_list[m_idx]
        tbase = tbase_list[s_idx] - tbase_list[m_idx]
        m_dop = dop_list[m_idx]
        s_dop = dop_list[s_idx]

        #Geometric decorrelation (Hanssen, 2001, Eq. 4.4.12)
        coh_geom = (pbase_c - abs(pbase)) / pbase_c
        if coh_geom < 0.:
            coh_geom = 0.

        #Doppler centroid decorrelation (Hanssen, 2001, Eq. 4.4.13)
        coh_dc = calculate_doppler_overlap(m_dop, s_dop, bandwidth_az)
        if coh_dc < 0.:
            coh_dc = 0.

        #Option 1: Temporal decorrelation - exponential delay model (Parizzi et al., 2009; Morishita and Hanssen, 2015)
        coh_temp = np.multiply((coh_thermal - coh_resid), np.exp(-1*abs(tbase)/decor_time)) + coh_resid

        coh = coh_geom * coh_dc * coh_temp
        cohs[i,:] = coh
    epsilon = 1e-3
    cohs[cohs < epsilon] = epsilon
    if display:
        print ''

    if display:
        print 'critical perp baseline: %.f m' % pbase_c
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
    '''Remove pairs/interoferogram with doppler overlap larger than critical value
    Inputs:
        date12_list : list of string, for date12 in YYMMDD-YYMMDD format
        date_list   : list of string, for date in YYMMDD/YYYYMMDD format, optional
        dop_list    : list of list of 3 float, for centroid Doppler frequency
        bandwidth_az    : float, bandwidth in azimuth direction
        dop_overlap_min : float, minimum overlap of azimuth Doppler frequency
    Outputs:
        date12_list : list of string, for date12 in YYMMDD-YYMMDD format
    '''
    if not date12_list:  return []
    # Get date6_list
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))
        if not len(date_list) == len(pbase_list):
            print 'ERROR: number of existing dates is not equal to number of perp baseline!'
            print 'date list is needed for threshold filtering!'
            print 'skip filtering.'
            return date12_list
    date6_list = ptime.yymmdd(date_list)

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        dop_overlap = calculate_doppler_overlap(dop_list[idx1], dop_list[idx2], bandwidth_az)
        if dop_overlap >= dop_overlap_min:
            date12_list_out.append(date12)
    return date12_list_out


def threshold_perp_baseline(date12_list, date_list, pbase_list, pbase_max, pbase_min=0.0):
    '''Remove pairs/interoferogram out of [pbase_min, pbase_max]
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
    '''
    if not date12_list:  return []
    # Get date6_list
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))
        if not len(date_list) == len(pbase_list):
            print 'ERROR: number of existing dates is not equal to number of perp baseline!'
            print 'date list is needed for threshold filtering!'
            print 'skip filtering.'
            return date12_list
    date6_list = ptime.yymmdd(date_list)

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        pbase = abs(pbase_list[idx1] - pbase_list[idx2])
        if pbase_min <= pbase <= pbase_max:
            date12_list_out.append(date12)
    return date12_list_out


def threshold_temporal_baseline(date12_list, btemp_max, keep_seasonal=True, btemp_min=0.0):
    '''Remove pairs/interferograms out of min/max/seasonal temporal baseline limits
    Inputs:
        date12_list : list of string for date12 in YYMMDD-YYMMDD format
        btemp_max   : float, maximum temporal baseline
        btemp_min   : float, minimum temporal baseline
        keep_seasonal : keep interferograms with seasonal temporal baseline
    Output:
        date12_list_out : list of string for date12 in YYMMDD-YYMMDD format
    Example:
        date12_list = threshold_temporal_baseline(date12_list, 200)
        date12_list = threshold_temporal_baseline(date12_list, 200, False)
    '''
    if not date12_list:  return []
    # Get date list and tbase list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date8_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))
    date6_list = ptime.yymmdd(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]

    # Threshold
    date12_list_out = []
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        tbase = int(abs(tbase_list[idx1] - tbase_list[idx2]))
        if btemp_min <= tbase <= btemp_max:
            date12_list_out.append(date12)
        elif keep_seasonal and tbase/30 in [11,12]:
            date12_list_out.append(date12)    
    return date12_list_out


def coherence_matrix(date12_list, coh_list, diagValue=np.nan):
    '''Return coherence matrix based on input date12 list and its coherence
    Inputs:
        date12_list - list of string in YYMMDD-YYMMDD format
        coh_list    - list of float, average coherence for each interferograms
    Output:
        coh_matrix  - 2D np.array with dimension length = date num
                      np.nan value for interferograms non-existed.
                      1.0 for diagonal elements
    '''
    # Get date list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date6_list = ptime.yymmdd(sorted(ptime.yyyymmdd(list(set(m_dates + s_dates)))))
    date_num = len(date6_list)

    coh_mat = np.zeros([date_num, date_num])
    coh_mat[:] = np.nan
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        coh = coh_list[date12_list.index(date12)]
        coh_mat[idx1, idx2] = coh    #symmetric
        coh_mat[idx2, idx1] = coh

    if diagValue is not np.nan:
        for i in range(date_num):    # diagonal value
            coh_mat[i, i] = diagValue

    return coh_mat


def threshold_coherence_based_mst(date12_list, coh_list):
    '''Return a minimum spanning tree of network based on the coherence inverse.
    Inputs:
        date12_list - list of string in YYMMDD-YYMMDD format
        coh_list    - list of float, average coherence for each interferogram
    Output:
        mst_date12_list - list of string in YYMMDD-YYMMDD format, for MST network of interferograms 
    '''
    # coh_list --> coh_mat --> weight_mat
    coh_mat = coherence_matrix(date12_list, coh_list)
    mask = ~np.isnan(coh_mat)
    wei_mat = np.zeros(coh_mat.shape)
    wei_mat[:] = np.inf
    wei_mat[mask] = 1/coh_mat[mask]

    # MST path based on weight matrix
    wei_mat_csr = sparse.csr_matrix(wei_mat)
    mst_mat_csr = minimum_spanning_tree(wei_mat_csr)

    # Get date6_list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date6_list = ptime.yymmdd(sorted(ptime.yyyymmdd(list(set(m_dates + s_dates)))))

    # Convert MST index matrix into date12 list
    [s_idx_list, m_idx_list] = [date_idx_array.tolist() for date_idx_array in sparse.find(mst_mat_csr)[0:2]]
    mst_date12_list = []
    for i in range(len(m_idx_list)):
        idx = sorted([m_idx_list[i], s_idx_list[i]])
        date12 = date6_list[idx[0]]+'-'+date6_list[idx[1]]
        mst_date12_list.append(date12)
    return mst_date12_list


def pair_sort(pairs):
    for idx in range(len(pairs)):
        if pairs[idx][0] > pairs[idx][1]:
            index1 = pairs[idx][1]
            index2 = pairs[idx][0]
            pairs[idx][0] = index1
            pairs[idx][1] = index2
    pairs=sorted(pairs)
    return pairs


def pair_merge(pairs1,pairs2):
    pairs = pairs1
    for pair in pairs2:
        if not pair in pairs:
            pairs.append(pair)

    pairs=sorted(pairs)
    return pairs


def select_pairs_all(date_list):
    '''Select All Possible Pairs/Interferograms
    Input : date_list   - list of date in YYMMDD/YYYYMMDD format
    Output: date12_list - list date12 in YYMMDD-YYMMDD format
    Reference:
        Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
        based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
    '''
    date8_list = sorted(ptime.yyyymmdd(date_list))
    date6_list = ptime.yymmdd(date8_list)
    date12_list = list(itertools.combinations(date6_list, 2))
    date12_list = [date12[0]+'-'+date12[1] for date12 in date12_list]
    return date12_list


def select_pairs_sequential(date_list, increment_num=2):
    '''Select Pairs in a Sequential way:
        For each acquisition, find its increment_num nearest acquisitions in the past time.
    Inputs:
        date_list  : list of date in YYMMDD/YYYYMMDD format
    Reference:
        Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series, IEEE TGRS, 51(7), 4249-4259.
    '''
    date8_list = sorted(ptime.yyyymmdd(date_list))
    date6_list = ptime.yymmdd(date8_list)
    date_idx_list = list(range(len(date6_list)))
    
    # Get pairs index list
    date12_idx_list = []
    for date_idx in date_idx_list:
        for i in range(increment_num):
            if date_idx-i-1 >= 0:
                date12_idx_list.append([date_idx-i-1, date_idx])
    date12_idx_list = [sorted(idx) for idx in sorted(date12_idx_list)]

    # Convert index into date12
    date12_list = [date6_list[idx[0]]+'-'+date6_list[idx[1]] for idx in date12_idx_list]
    return date12_list


def select_pairs_hierarchical(date_list, pbase_list, temp_perp_list):
    '''Select Pairs in a hierarchical way using list of temporal and perpendicular baseline thresholds
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
    '''
    # Get all date12
    date12_list_all = select_pairs_all(date_list)

    # Loop of Threshold
    print 'List of temporal and perpendicular spatial baseline thresholds:'
    print temp_perp_list
    date12_list = []
    for temp_perp in temp_perp_list:
        tbase_max = temp_perp[0]
        pbase_max = temp_perp[1]
        date12_list_tmp = threshold_temporal_baseline(date12_list_all, tbase_max, keep_seasonal=False)
        date12_list_tmp = threshold_perp_baseline(date12_list_tmp, date_list, pbase_list, pbase_max)
        date12_list += date12_list_tmp
    date12_list = sorted(list(set(date12_list)))
    return date12_list


def select_pairs_delaunay(date_list, pbase_list, norm=True):
    '''Select Pairs using Delaunay Triangulation based on temporal/perpendicular baselines
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
    '''
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
    date12_list = [date6_list[idx[0]]+'-'+date6_list[idx[1]] for idx in date12_idx_list]
    return date12_list


def select_pairs_mst(date_list, pbase_list):
    '''Select Pairs using Minimum Spanning Tree technique
        Connection Cost is calculated using the baseline distance in perp and scaled temporal baseline (Pepe and Lanari,
        2006, TGRS) plane.
    Inputs:
        date_list  : list of date in YYMMDD/YYYYMMDD format
        pbase_list : list of float, perpendicular spatial baseline
    References:
        Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
        of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
        Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
    '''
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

    weightMat = np.sqrt(np.square(ttMat) + np.square(ppMat))  # 2D distance matrix in temp/perp domain
    weightMat = sparse.csr_matrix(weightMat)  # compress sparse row matrix

    # MST path based on weight matrix
    mstMat = sparse.csgraph.minimum_spanning_tree(weightMat)

    # Convert MST index matrix into date12 list
    [s_idx_list, m_idx_list] = [date_idx_array.tolist() for date_idx_array in sparse.find(mstMat)[0:2]]
    date12_list = []
    for i in range(len(m_idx_list)):
        idx = sorted([m_idx_list[i], s_idx_list[i]])
        date12 = date6_list[idx[0]]+'-'+date6_list[idx[1]]
        date12_list.append(date12)
    return date12_list


def select_pairs_star(date_list, m_date=None, pbase_list=[]):
    '''Select Star-like network/interferograms/pairs, it's a single master network, similar to PS approach.
    Usage:
        m_date : master date, choose it based on the following cretiria:
                 1) near the center in temporal and spatial baseline
                 2) prefer winter season than summer season for less temporal decorrelation
    Reference:
        Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
    '''
    date8_list = sorted(ptime.yyyymmdd(date_list))
    date6_list = ptime.yymmdd(date8_list)
    
    # Select master date if not existed
    if not m_date:
        m_date = select_master_date(date8_list, pbase_list)
        print 'auto select master date: '+m_date
    
    # Check input master date
    m_date8 = ptime.yyyymmdd(m_date)
    if m_date8 not in date8_list:
        print 'Input master date is not existed in date list!'
        print 'Input master date: '+str(m_date8)
        print 'Input date list: '+str(date8_list)
        m_date8 = None
    
    # Generate star/ps network
    m_idx = date8_list.index(m_date8)
    date12_idx_list = [sorted([m_idx, s_idx]) for s_idx in range(len(date8_list)) if s_idx is not m_idx]
    date12_list = [date6_list[idx[0]]+'-'+date6_list[idx[1]] for idx in date12_idx_list]
    
    return date12_list


def select_master_date(date_list, pbase_list=[]):
    '''Select super master date based on input temporal and/or perpendicular baseline info.
    Return master date in YYYYMMDD format.
    '''
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
        ttMat1, ttMat2 = np.meshgrid(np.array(tbase_list), np.array(tbase_list))
        ppMat1, ppMat2 = np.meshgrid(np.array(pbase_list), np.array(pbase_list))
        ttMat = np.abs(ttMat1 - ttMat2)  # temporal distance matrix
        ppMat = np.abs(ppMat1 - ppMat2)  # spatial distance matrix
        disMat = np.sqrt(np.square(ttMat) + np.square(ppMat))  # 2D distance matrix in temp/perp domain
        
        # Choose date minimize the total distance of temp/perp baseline
        disMean = np.mean(disMat, 0)
        m_idx = np.argmin(disMean)
        m_date8 = date8_list[m_idx]
    return m_date8


def select_master_interferogram(date12_list, date_list, pbase_list, m_date=None):
    '''Select reference interferogram based on input temp/perp baseline info
    If master_date is specified, select its closest slave_date, which is newer than master_date;
        otherwise, choose the closest pair among all pairs as master interferogram.
    Example:
        master_date12   = pnet.select_master_ifgram(date12_list, date_list, pbase_list)
        '080211-080326' = pnet.select_master_ifgram(date12_list, date_list, pbase_list, m_date='080211')
    '''
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
    base_distance = np.sqrt((tbase_array[idx2] - tbase_array[idx1])**2 + (pbase_array[idx2] - pbase_array[idx1])**2)
    
    # Get master interferogram index
    if not m_date:
        # Choose pair with shortest temp/perp baseline
        m_date12_idx = np.argmin(base_distance)
    else:
        m_date = ptime.yymmdd(m_date)
        # Choose pair contains m_date with shortest temp/perp baseline
        m_date12_idx_array = np.array([date12_list.index(date12) for date12 in date12_list if m_date+'-' in date12])
        min_base_distance = np.min(base_distance[m_date12_idx_array])
        m_date12_idx = np.where(base_distance == min_base_distance)[0][0]
    
    m_date12 = date12_list[m_date12_idx]
    return m_date12


##################################################################
def plot_network(ax, date12_list, date_list, pbase_list, plot_dict={}, date12_list_drop=[], print_msg=True):
    '''Plot Temporal-Perp baseline Network
    Inputs
        ax : matplotlib axes object
        date12_list : list of string for date12 in YYMMDD-YYMMDD format
        date_list   : list of string, for date in YYYYMMDD/YYMMDD format
        pbase_list  : list of float, perp baseline, len=number of acquisition
        plot_dict   : dictionary with the following items:
                      fontsize
                      linewidth
                      markercolor
                      markersize

                      coherence_list : list of float, coherence value of each interferogram, len = number of ifgrams
                      disp_min/max :  float, min/max range of the color display based on coherence_list
                      colormap : string, colormap name
                      coh_thres : float, coherence of where to cut the colormap for display
                      disp_title : bool, show figure title or not, default: True
                      disp_drop: bool, show dropped interferograms or not, default: True
    Output
        ax : matplotlib axes object
    '''
    
    # Figure Setting
    keyList = plot_dict.keys()
    if not 'fontsize'    in keyList:   plot_dict['fontsize']    = 12
    if not 'linewidth'   in keyList:   plot_dict['linewidth']   = 2
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16
    # For colorful display of coherence
    if not 'coherence_list' in keyList:  plot_dict['coherence_list'] = None
    if not 'cbar_label'     in keyList:  plot_dict['cbar_label']     = 'Coherence'
    if not 'disp_min'       in keyList:  plot_dict['disp_min']       = 0.2
    if not 'disp_max'       in keyList:  plot_dict['disp_max']       = 1.0
    if not 'colormap'       in keyList:  plot_dict['colormap']       = 'RdBu'
    if not 'disp_title'     in keyList:  plot_dict['disp_title']     = True
    if not 'coh_thres'      in keyList:  plot_dict['coh_thres']      = None
    if not 'disp_drop'      in keyList:  plot_dict['disp_drop']      = True
    coh_list = plot_dict['coherence_list']
    disp_min = plot_dict['disp_min']
    disp_max = plot_dict['disp_max']
    coh_thres = plot_dict['coh_thres']
    transparency = 0.7

    # Date Convert
    date8_list = ptime.yyyymmdd(sorted(date_list))
    date6_list = ptime.yymmdd(date8_list)
    dates, datevector = ptime.date_list2vector(date8_list)
    tbase_list = ptime.date_list2tbase(date8_list)[0]

    ## maxBperp and maxBtemp
    ifgram_num = len(date12_list)
    pbase12 = np.zeros(ifgram_num)
    tbase12 = np.zeros(ifgram_num)
    for i in range(ifgram_num):
        m_date, s_date = date12_list[i].split('-')
        m_idx = date6_list.index(m_date)
        s_idx = date6_list.index(s_date)
        pbase12[i] = pbase_list[s_idx] - pbase_list[m_idx]
        tbase12[i] = tbase_list[s_idx] - tbase_list[m_idx]
    print 'max perpendicular baseline: %.2f m' % (np.max(np.abs(pbase12)))
    print 'max temporal      baseline: %d days' % (np.max(tbase12))

    ## Keep/Drop - date12
    date12_list_keep = sorted(list(set(date12_list) - set(date12_list_drop)))
    idx_date12_keep = [date12_list.index(i) for i in date12_list_keep]
    idx_date12_drop = [date12_list.index(i) for i in date12_list_drop]
    if not date12_list_drop:
        plot_dict['disp_drop'] = False

    ## Keep/Drop - date
    m_dates = [i.split('-')[0] for i in date12_list_keep]
    s_dates = [i.split('-')[1] for i in date12_list_keep]
    date8_list_keep = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
    date8_list_drop = sorted(list(set(date8_list) - set(date8_list_keep)))
    idx_date_keep = [date8_list.index(i) for i in date8_list_keep]
    idx_date_drop = [date8_list.index(i) for i in date8_list_drop]


    # Ploting
    #ax=fig.add_subplot(111)
    ## Colorbar when conherence is colored
    if coh_list:
        data_min = min(coh_list)
        data_max = max(coh_list)
        # Normalize
        normalization = False
        if normalization:
            coh_list = [(coh-data_min) / (data_min-data_min) for coh in coh_list]
            disp_min = data_min
            disp_max = data_max

        if print_msg:
            print 'showing coherence'
            print 'colormap: '+plot_dict['colormap']
            print 'display range: '+str([disp_min, disp_max])
            print 'data    range: '+str([data_min, data_max])

        # Use lower/upper part of colormap to emphasis dropped interferograms
        if not coh_thres:
            # Find proper cut percentage so that all keep pairs are blue and drop pairs are red
            coh_list_keep = [coh_list[i] for i in idx_date12_keep]
            coh_list_drop = [coh_list[i] for i in idx_date12_drop]
            if coh_list_drop:
                coh_thres = max(coh_list_drop)
            else:
                coh_thres = min(coh_list_keep)

        if coh_thres < disp_min:
            disp_min = 0.0
            if print_msg:
                print 'data range exceed orginal display range, set new display range to: [0.0, %f]' % (disp_max)
        c1_num = np.ceil(200.0 * (coh_thres - disp_min) / (disp_max - disp_min)).astype('int')
        coh_thres = c1_num / 200.0 * (disp_max-disp_min) + disp_min
        cmap = plt.get_cmap(plot_dict['colormap'])
        colors1 = cmap(np.linspace(0.0, 0.3, c1_num))
        colors2 = cmap(np.linspace(0.6, 1.0, 200 - c1_num))
        cmap = colors.LinearSegmentedColormap.from_list('truncate_RdBu', np.vstack((colors1, colors2)))
        if print_msg:
            print 'color jump at '+str(coh_thres)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        norm = mpl.colors.Normalize(vmin=disp_min, vmax=disp_max)
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        cbar.set_label(plot_dict['cbar_label'], fontsize=plot_dict['fontsize'])

        #plot low coherent ifgram first and high coherence ifgram later
        coh_list_keep = [coh_list[date12_list.index(i)] for i in date12_list_keep]
        date12_list_keep = [x for _,x in sorted(zip(coh_list_keep, date12_list_keep))]

    ## Dot - SAR Acquisition
    if idx_date_keep:
        x_list = [dates[i] for i in idx_date_keep]
        y_list = [pbase_list[i] for i in idx_date_keep]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    if idx_date_drop:
        x_list = [dates[i] for i in idx_date_drop]
        y_list = [pbase_list[i] for i in idx_date_drop]
        ax.plot(x_list, y_list, 'ko', alpha=0.7, ms=plot_dict['markersize'], mfc='gray')

    ## Line - Pair/Interferogram        
    # interferograms dropped
    if plot_dict['disp_drop']:
        for date12 in date12_list_drop:
            date1, date2 = date12.split('-')
            idx1 = date6_list.index(date1)
            idx2 = date6_list.index(date2)
            x = np.array([dates[idx1], dates[idx2]])
            y = np.array([pbase_list[idx1], pbase_list[idx2]])
            if coh_list:
                coh = coh_list[date12_list.index(date12)]
                coh_idx = (coh - disp_min) / (disp_max - disp_min)
                ax.plot(x, y, '--', lw=plot_dict['linewidth'], alpha=transparency, c=cmap(coh_idx)) 
            else:
                ax.plot(x, y, '--', lw=plot_dict['linewidth'], alpha=transparency, c='k')

    # interferograms kept
    for date12 in date12_list_keep:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        x = np.array([dates[idx1], dates[idx2]])
        y = np.array([pbase_list[idx1], pbase_list[idx2]])
        if coh_list:
            coh = coh_list[date12_list.index(date12)]
            coh_idx = (coh - disp_min) / (disp_max - disp_min)
            ax.plot(x, y, '-', lw=plot_dict['linewidth'], alpha=transparency, c=cmap(coh_idx)) 
        else:
            ax.plot(x, y, '-', lw=plot_dict['linewidth'], alpha=transparency, c='k')

    if plot_dict['disp_title']:
        ax.set_title('Interferogram Network', fontsize=plot_dict['fontsize'])

    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])[0]
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perp Baseline [m]',fontsize=plot_dict['fontsize'])

    # Legend
    if plot_dict['disp_drop']:
        solid_line = mlines.Line2D([],[],color='k',ls='solid', label='Interferograms')
        dash_line  = mlines.Line2D([],[],color='k',ls='dashed', label='Interferograms dropped')
        ax.legend(handles=[solid_line,dash_line])

    return ax


def plot_perp_baseline_hist(ax, date8_list, pbase_list, plot_dict={}, date8_list_drop=[]):
    ''' Plot Perpendicular Spatial Baseline History
    Inputs
        ax : matplotlib axes object
        date8_list : list of string, date in YYYYMMDD format
        pbase_list : list of float, perp baseline 
        plot_dict : dictionary with the following items:
                    fontsize
                    linewidth
                    markercolor
                    markersize
                    disp_title : bool, show figure title or not, default: True
        date8_list_drop : list of string, date dropped in YYYYMMDD format
                          e.g. ['20080711', '20081011']
    Output:
        ax : matplotlib axes object
    '''
    # Figure Setting
    keyList = plot_dict.keys()
    if not 'fontsize'    in keyList:   plot_dict['fontsize']    = 12
    if not 'linewidth'   in keyList:   plot_dict['linewidth']   = 2
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16
    if not 'disp_title'  in keyList:   plot_dict['disp_title']  = True
    transparency = 0.7

    # Date Convert
    dates, datevector = ptime.date_list2vector(date8_list)

    # Get index of date used and dropped
    #date8_list_drop = ['20080711', '20081011']  # for debug
    idx_keep = range(len(date8_list))
    idx_drop = []
    for i in date8_list_drop:
        idx = date8_list.index(i)
        idx_keep.remove(idx)
        idx_drop.append(idx)

    # Plot
    #ax=fig.add_subplot(111)

    # Plot date used
    if idx_keep:
        x_list = [dates[i] for i in idx_keep]
        y_list = [pbase_list[i] for i in idx_keep]
        ax.plot(x_list, y_list, '-ko', alpha=transparency, lw=plot_dict['linewidth'], \
                ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    
    # Plot date dropped
    if idx_drop:
        x_list = [dates[i] for i in idx_drop]
        y_list = [pbase_list[i] for i in idx_drop]
        ax.plot(x_list, y_list, 'ko', alpha=transparency, ms=plot_dict['markersize'], mfc='gray')

    if plot_dict['disp_title']:
        ax.set_title('Perpendicular Baseline History',fontsize=plot_dict['fontsize'])

    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])[0]
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]',fontsize=plot_dict['fontsize'])

    return ax


def plot_coherence_matrix(ax, date12_list, coherence_list, date12_list_drop=[], plot_dict={}):
    '''Plot Coherence Matrix of input network
    
    if date12_list_drop is not empty, plot KEPT pairs in the upper triangle and
                                           ALL  pairs in the lower triangle.
    '''
    # Figure Setting
    keyList = plot_dict.keys()
    if not 'fontsize'    in keyList:   plot_dict['fontsize']    = 12
    if not 'linewidth'   in keyList:   plot_dict['linewidth']   = 2
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16
    if not 'disp_title'  in keyList:   plot_dict['disp_title']  = True
    if not 'cbar_label'  in keyList:   plot_dict['cbar_label']  = 'Coherence'

    coh_mat = coherence_matrix(date12_list, coherence_list)

    if date12_list_drop:
        # Date Convert
        m_dates = [i.split('-')[0] for i in date12_list]
        s_dates = [i.split('-')[1] for i in date12_list]
        date6_list = ptime.yymmdd(sorted(list(set(m_dates + s_dates))))
        # Set dropped pairs' value to nan, in upper triangle only.
        for date12 in date12_list_drop:
            idx1, idx2 = [date6_list.index(i) for i in date12.split('-')]
            coh_mat[idx1, idx2] = np.nan

    #Show diagonal value as black, to be distinguished from un-selected interferograms
    diag_mat = np.diag(np.ones(coh_mat.shape[0]))
    diag_mat[diag_mat == 0.] = np.nan
    im = ax.imshow(diag_mat, cmap='gray_r', vmin=0.0, vmax=1.0, interpolation='nearest')
    #Show coherence matrix
    im = ax.imshow(coh_mat, cmap='jet', vmin=0.0, vmax=1.0, interpolation='nearest')

    date_num = coh_mat.shape[0]
    if date_num < 30:
        tick_list = range(0,date_num,5)
    else:
        tick_list = range(0,date_num,10)
    ax.get_xaxis().set_ticks(tick_list)
    ax.get_yaxis().set_ticks(tick_list)
    ax.set_xlabel('Image Number', fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Image Number', fontsize=plot_dict['fontsize'])

    if plot_dict['disp_title']:
        ax.set_title('Coherence Matrix')

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="3%")
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(plot_dict['cbar_label'], fontsize=plot_dict['fontsize'])

    # Legend
    if date12_list_drop:
        ax.plot([],[],label='Upper: used ifgrams')
        ax.plot([],[],label='Lower: all ifgrams')
        ax.legend(handlelength=0)

    return ax


def mode (thelist):
    '''Find Mode (most common) item in the list
    Borrowded from pysar._pysar_utilities
    '''
    if not thelist:
        return None
    if len(thelist) == 1:
        return thelist[0]

    counts = {}
    for item in thelist:
        counts[item] = counts.get(item, 0) + 1
    maxcount = 0
    maxitem  = None
    for k, v in counts.items():
        if v > maxcount:
            maxitem  = k
            maxcount = v

    if maxcount == 1:
        print "All values only appear once"
        return None
    elif counts.values().count(maxcount) > 1:
        print "List has multiple modes"
        return None
    else:
        return maxitem


def plot_coherence_history(ax, date12_list, coherence_list, plot_dict={}):
    '''Plot min/max Coherence of all interferograms for each date'''
    # Figure Setting
    keyList = plot_dict.keys()
    if not 'fontsize'    in keyList:   plot_dict['fontsize']    = 12
    if not 'linewidth'   in keyList:   plot_dict['linewidth']   = 2
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16
    if not 'disp_title'  in keyList:   plot_dict['disp_title']  = True

    # Get date list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    date8_list = sorted(ptime.yyyymmdd(list(set(m_dates + s_dates))))

    dates, datevector = ptime.date_list2vector(date8_list)
    bar_width = mode(np.diff(dates).tolist())*3/4
    x_list = [i-bar_width/2 for i in dates]

    coh_mat = coherence_matrix(date12_list, coherence_list)

    ax.bar(x_list, np.nanmax(coh_mat, axis=0), bar_width.days, label='Max Coherence')
    ax.bar(x_list, np.nanmin(coh_mat, axis=0), bar_width.days, label='Min Coherence')

    if plot_dict['disp_title']:
        ax.set_title('Coherence History of All Related Interferograms')

    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])[0]
    ax.set_ylim([0.0,1.0])

    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Coherence',fontsize=plot_dict['fontsize'])
    ax.legend(loc='lower right')

    return ax


def auto_adjust_yaxis(ax, dataList, fontSize=12, ymin=None, ymax=None):
    '''Adjust Y axis
    Input:
        ax       : matplot figure axes object
        dataList : list of float, value in y axis
        fontSize : float, font size
        ymin     : float, lower y axis limit
        ymax     : float, upper y axis limit
    Output:
        ax
    '''
    # Min/Max
    dataRange = max(dataList) - min(dataList)
    if ymin is None:  ymin = min(dataList) - 0.1*dataRange
    if ymax is None:  ymax = max(dataList) + 0.1*dataRange
    ax.set_ylim([ymin, ymax])
    ## Tick/Label setting
    #xticklabels = plt.getp(ax, 'xticklabels')
    #yticklabels = plt.getp(ax, 'yticklabels')
    #plt.setp(yticklabels, 'color', 'k', fontsize=fontSize)
    #plt.setp(xticklabels, 'color', 'k', fontsize=fontSize)

    return ax

