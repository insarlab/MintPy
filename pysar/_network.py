#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar._network as pnet
#


import os
import datetime
import itertools

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.sparse import csr_matrix, find
from scipy.sparse.csgraph import minimum_spanning_tree

import pysar._datetime as ptime
import pysar._readfile as readfile


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

    dateList  = ptime.igram_date_list(igramFile)
    dateList6 = ptime.yymmdd(dateList)

    pairs = []
    igramList=h5file[k[0]].keys()
    for igram in igramList:
        date12 = h5file[k[0]][igram].attrs['DATE12'].split('-')
        pairs.append([dateList6.index(date12[0]),dateList6.index(date12[1])])
    h5file.close()

    pairs = pair_sort(pairs)

    return pairs


def read_baseline_file(baselineFile,exDateList=[]):
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


def get_date12_list(File):
    '''Read Date12 info from input file: Pairs.list or multi-group hdf5 file
    Example:
        date12List = get_date12_list('unwrapIfgram.h5')
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
            date12 = h5[k][epoch].attrs['DATE12']
            date12_list.append(date12)
        h5.close()
    else:
        date12_list = list(np.loadtxt(File, dtype=str))
    
    date12_list = sorted(date12_list)
    return date12_list


def igram_perp_baseline_list(File):
    '''Get perpendicular baseline list from input multi_group hdf5 file'''
    print 'read perp baseline info from '+File
    p_baseline_list = []
    k = readfile.read_attribute(File)['FILE_TYPE']
    h5 = h5py.File(File, 'r')
    epochList = sorted(h5[k].keys())
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
    elif  sensor == 'Jers' :  bandwidth =   900.0   # FA 11/2015 just copied, need to research
    elif  sensor == 'Alos' :  bandwidth =  1720.0
    elif  sensor == 'Tsx'  :  bandwidth = 15000.0
    elif  sensor == 'Csk'  :  bandwidth = 15000.0   # FA 9/2015  shoud be checked
    elif  sensor == 'Rsat' :  bandwidth =   900.0
    elif  sensor == 'Rsat2':  bandwidth =   900.0
    elif  sensor == 'S1'   :  bandwidth =  4000.0   # shong 08/2016 sould be checked
    elif  sensor == 'Kmps5':  bandwidth = 15000.0   # shong 08/2016 sould be checked
    else: print 'satellite not found'; bandwidth = None
    return bandwidth


def range_bandwidth(sensor):
    ## Range Bandwidth in Hz for the given satellite
    if    sensor == 'Alos':  bandwidth = 14e6      # for FBD, 28MHz for FBS
    elif  sensor == 'Ers' :  bandwidth = 15.55e6
    elif  sensor == 'Jers':  bandwidth = 15e6      # Jers only has HH pol
    elif  sensor == 'Tsx' :  bandwidth = 150e6
    return bandwidth


def wavelength(sensor):
    c = 299792458;   # m/s, speed of light
    if    sensor == 'Alos':  center_frequency = 1.270e9
    elif  sensor == 'Ers' :  center_frequency = 5.3e9
    elif  sensor == 'Jers':  center_frequency = 1.275e9  # Hz
    elif  sensor == 'Tsx' :  center_frequency = 9.65e9
    wavelength = c / center_frequency
    return wavelength


def incidence_angle(sensor):
    if   sensor == 'Alos':  inc_angle = 34.3     # degree, for ALOS PALSAR Fine mode
    elif sensor == 'Jers':  inc_angle = 35.21
    elif sensor == 'Tsx' :  inc_angle = 39.23    # Yunjun 5/2016, for TaizhouTsx, not sure it's for all cases.
    return inc_angle


def signal2noise_ratio(sensor):
    '''Fine the Signal to Noise Ratio in dB for the given satellite
    Reference:
        ERS - Zebker et al., 1994, TGRS
        Envisat - Guarnieri, A.M., 2013. Introduction to RADAR. POLIMI DEI, Milano.
        JERS - https://directory.eoportal.org/web/eoportal/satellite-missions/j/jers-1
    '''
    if   sensor == 'Ers' :  SNR = 11.7
    elif sensor == 'Env' :  SNR = 19.5
    elif sensor == 'Jers':  SNR = 14
    else: print 'satellite not found'; SNR = None
    return SNR


def critical_perp_baseline(sensor):
    '''Critical Perpendicular Baseline for each satellite'''
    # Jers: 5.712e3 m (near_range=688849.0551m)
    # Alos: 6.331e3 m (near_range=842663.2917m)
    # Tsx : 8.053e3 m (near_range=634509.1271m)

    c = 299792458;   # m/s, speed of light
    wvl = wavelength(sensor)
    near_range = 688849;        # Yunjun 5/2016, case for Jers, need a automatic way to get this number
    rg_bandwidth = range_bandwidth(sensor)
    inc_angle    = incidence_angle(sensor) / 180 * np.pi;
    Bperp_c      = wvl * (rg_bandwidth/c) * near_range * math.tan(inc_angle)
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
    date12_list = select_pairs_all(date_list)
    
    # Loop of Threshold
    for temp_perp in temp_perp_list:
        tbase_max = temp_perp[0]
        pbase_max = temp_perp[1]
        date12_list = threshold_temporal_baseline(date12_list, tbase_max, keep_seasonal=False)
        date12_list = threshold_perp_baseline(date12_list, date_list, pbase_list, pbase_max)
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
    weightMat = csr_matrix(weightMat)  # compress sparse row matrix

    # MST path based on weight matrix
    mstMat = minimum_spanning_tree(weightMat)

    # Convert MST index matrix into date12 list
    [s_idx_list, m_idx_list] = [date_idx_array.tolist() for date_idx_array in find(mstMat)[0:2]]
    date12_list = [date6_list[m_idx_list[i]]+'-'+date6_list[s_idx_list[i]] for i in range(len(m_idx_list))]
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
    
    # Check input master date
    m_date8 = ptime.yyyymmdd(m_date)
    if m_date8 not in date8_list:
        print 'Input master date is not existed in date list!'
        print 'Input master date: '+m_date8
        print 'Input date list: '+str(date8_list)
        m_date8 = None
    
    # Select master date if not existed
    if not m_date8:
        m_date8 = select_master_date(date8_list, pbase_list)
    
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


def select_master_ifgram(date12_list, date_list, pbase_list, m_date=None):
    '''Select reference interferogram based on input temp/perp baseline info
    If master_date is specified, select its closest slave_date; otherwise, choose the closest pair
    among all pairs as master interferogram.
    Example:
        master_date12 = pnet.select_master_ifgram(date12_list, date_list, pbase_list)
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
        m_date12_idx_array = np.array([date12_list.index(date12) for date12 in date12_list if m_date in date12])
        min_base_distance = np.min(base_distance[m_date12_idx_array])
        m_date12_idx = np.where(base_distance==min_dis)[0][0]
    
    m_date12 = date12_list[m_date12_idx]
    return m_date12


##################################################################
def plot_network(ax, date12_list, date_list, pbase_list, plot_dict={}):
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
    if not 'disp_min'       in keyList:  plot_dict['disp_min']       = 0.4
    if not 'disp_max'       in keyList:  plot_dict['disp_max']       = 1.0
    if not 'colormap'       in keyList:  plot_dict['colormap']       = 'RdBu'

    # Date Convert
    date8_list = ptime.yyyymmdd(date_list)
    date6_list = ptime.yymmdd(date8_list)
    dates, datevector = ptime.date_list2vector(date8_list)

    # Ploting
    #ax=fig.add_subplot(111)
    # Colorbar when conherence is colored
    if plot_dict['coherence_list']:
        data_min = min(plot_dict['coherence_list'])
        data_max = max(plot_dict['coherence_list'])
        # Normalize
        normalization = False
        if normalization:
            plot_dict['coherence_list'] = [(coh-data_min) / (data_min-data_min) for coh in plot_dict['coherence_list']]
            plot_dict['disp_min'] = data_min
            plot_dict['disp_max'] = data_max
        
        print 'showing coherence'
        print 'colormap: '+plot_dict['colormap']
        print 'display range: '+str([plot_dict['disp_min'], plot_dict['disp_max']])
        print 'data    range: '+str([data_min, data_max])
        cmap = plt.get_cmap(plot_dict['colormap'])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        norm = mpl.colors.Normalize(vmin=plot_dict['disp_min'], vmax=plot_dict['disp_max'])
        cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        cbar.set_label('Coherence')

    # Dot - SAR Acquisition
    ax.plot(dates, pbase_list, 'ko', lw=plot_dict['linewidth'], alpha=0.7,\
            ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])

    # Line - Pair/Interferogram
    for date12 in date12_list:
        date1, date2 = date12.split('-')
        idx1 = date6_list.index(date1)
        idx2 = date6_list.index(date2)
        x = np.array([dates[idx1], dates[idx2]])
        y = np.array([pbase_list[idx1], pbase_list[idx2]])
        if plot_dict['coherence_list']:
            coh_idx = (plot_dict['coherence_list'][i] - plot_dict['disp_min']) /\
                      (plot_dict['disp_max'] - plot_dict['disp_min'])
            ax.plot(x, y, lw=plot_dict['linewidth'], alpha=0.7, c=cmap(coh_idx)) 
        else:
            ax.plot(x, y, lw=plot_dict['linewidth'], alpha=0.7, c='k')

    ax.set_title('Interferogram Network', fontsize=plot_dict['fontsize'])
    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]',fontsize=plot_dict['fontsize'])

    return ax


def plot_perp_baseline_hist(ax, date8_list, pbase_list, plot_dict={}):
    ''' Plot Perpendicular Spatial Baseline History
    Inputs
        ax : matplotlib axes object
        date8_list : list of 8-digit string, date 
        pbase_list : list of float, perp baseline 
        plot_dict : dictionary with the following items:
                    fontsize
                    linewidth
                    markercolor
                    markersize
    Output:
        ax : matplotlib axes object
    '''
    # Figure Setting
    keyList = plot_dict.keys()
    if not 'fontsize'    in keyList:   plot_dict['fontsize']    = 12
    if not 'linewidth'   in keyList:   plot_dict['linewidth']   = 2
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'orange'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16

    # Date Convert
    dates, datevector = ptime.date_list2vector(date8_list)

    # Plot
    #ax=fig.add_subplot(111)
    ax.plot(dates, pbase_list, '-ko', lw=plot_dict['linewidth'], alpha=0.7,\
            ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    ax.set_title('Perpendicular Baseline History',fontsize=plot_dict['fontsize'])

    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])
    ax = auto_adjust_yaxis(ax, pbase_list, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]',fontsize=plot_dict['fontsize'])

    return ax


def auto_adjust_yaxis(ax, dataList, fontSize=12):
    '''Adjust Y axis
    Input:
        ax - matplot figure axes object
        dataList : list of float, value in y axis
    '''
    # Min/Max
    dataRange = max(dataList) - min(dataList)
    ax.set_ylim(min(dataList) - 0.1*dataRange,\
                max(dataList) + 0.1*dataRange)
    ## Tick/Label setting
    #xticklabels = plt.getp(ax, 'xticklabels')
    #yticklabels = plt.getp(ax, 'yticklabels')
    #plt.setp(yticklabels, 'color', 'k', fontsize=fontSize)
    #plt.setp(xticklabels, 'color', 'k', fontsize=fontSize)
    
    return ax

