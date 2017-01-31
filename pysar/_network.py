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


################################# Basic File I/O #################################
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
    Examples:
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile)
        date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile,['080520','100726'])
        date8List, perpBaseList = read_baseline_file(baselineFile)[0:2]
    '''

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
    prfList      = []
    slcDirList   = []
    for line in lines:
        c = line.split()    # splits on white space
        date = c[0]
        if not date in exDateList:
            date6List.append(date)
            perpBaseList.append(float(c[1]))
            dopplerList.append(np.array([float(c[2]), float(c[3]), float(c[4])]))
            prfList.append(float(c[5]))
            slcDirList.append(c[6])

    date8List = ptime.yyyymmdd(date6List)

    return date8List, perpBaseList, dopplerList, prfList, slcDirList


def date12_list2index(date12_list, date_list=[]):
    '''Convert list of date12 string into list of index'''
    # Get dateList from date12List
    if not date_list:
        m_dates = [date12.split('-')[0] for date12 in date12_list]
        s_dates = [date12.split('-')[1] for date12 in date12_list]
        date_list = sorted(list(set(m_dates + s_dates)))
    date6_list = ptime.yymmdd(date_list)
    
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


################################# Network Selection #################################
def threshold_perp_baseline(igramIdxList,perpBaseList,perpBaseMax=800,perpBaseMin=0):
    '''Remove pairs/interoferogram out of [perpBaseMin, perpBaseMax]
    Example:
        pairs = threshold_perp_baseline(pairs,perpBaseList,500)
    '''
    igramIdxListOut = []
    for igramIdx in igramIdxList:
        idx1 = igramIdx[0]
        idx2 = igramIdx[1]
        perpBase = abs(perpBaseList[idx1]-perpBaseList[idx2])
        if perpBaseMin <= perpBase <= perpBaseMax:
            igramIdxListOut.append(igramIdx)

    return igramIdxListOut


def threshold_temporal_baseline(igramIdxList,tempBaseList,tempBaseMax=365,seasonal=1,tempBaseMin=0):
    '''Remove pairs/interferograms out of min/max/seasonal temporal baseline limits
    Usage:
        seasonal : keep interferograms with seasonal temporal baseline
                   1 - keep them, by default
                   0 - do not keep them
    Example:
        pairs = threshold_temporal_baseline(pairs,tempBaseList,80)
        pairs = threshold_temporal_baseline(pairs,tempBaseList,80,0)  # disable seasonal checking
    '''
    igramIdxListOut = []
    
    for igramIdx in igramIdxList:
        idx1 = igramIdx[0]
        idx2 = igramIdx[1]
        tempBase = abs(tempBaseList[idx1]-tempBaseList[idx2])
        if seasonal == 1:
            if tempBaseMin <= tempBase <= tempBaseMax or tempBase/30 in [11,12]:
                igramIdxListOut.append(igramIdx)
        else:
            if tempBaseMin <= tempBase <= tempBaseMax:
                igramIdxListOut.append(igramIdx)

    return igramIdxListOut


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


def select_pairs_all(dateList):
    '''Select All Possible Pairs/Interferograms
    Reference:
        Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
        based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
    '''

    ## Get date index
    indexList = []
    for ni in range(len(dateList)):  indexList.append(ni)

    ## Generate All Pairs
    allPairs = list(itertools.combinations(indexList,2))

    ## Convert tuple object to list object
    allPairs2=[]
    for pair in allPairs: allPairs2.append(list(pair))

    return allPairs2


def select_pairs_delaunay(tempBaseList,perpBaseList,normalize=1):
    '''Select Pairs using Delaunay Triangulation based on temporal/perpendicular baselines
    Usage:
        tempBaseList : list of temporal baseline
        perpBaseList : list of perpendicular spatial baseline
        normalize    : normalize temporal baseline to perpendicular baseline
                       1 - enable  normalization, default
                       0 - disable normalization
    Key points
        1. Define a ratio between perpendicular and temporal baseline axis units (Pepe and Lanari, 2006, TGRS).
        2. Pairs with too large perpendicular / temporal baseline or Doppler centroid difference should be removed
           after this, using a threshold, to avoid strong decorrelations (Zebker and Villasenor, 1992, TGRS).
    Reference:
        Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
        of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
        Zebker, H. A., and J. Villasenor (1992), Decorrelation in interferometric radar echoes, IEEE TGRS, 30(5), 950-959.
    '''

    ##### Generate Delaunay Triangulation based on temporal and spatial perpendicular baselines
    if normalize == 0:
        delaunayPairs = Triangulation(tempBase,perpBase).edges.tolist()
    else:
        ##### Ratio between perpendicular and temporal baselines (Pepe and Lanari, 2006, TGRS)
        tempBaseFacList = []
        temp2perp_scale = (max(perpBaseList)-min(perpBaseList)) / (max(tempBaseList)-min(tempBaseList))
        for tempBase in tempBaseList:
            #tempBaseFac = (tempBase - min(tempBaseList)) * temp2perp_scale + min(perpBaseList)
            tempBaseFac = tempBase * temp2perp_scale   # giving same result as the line above
            tempBaseFacList.append(tempBaseFac)
        delaunayPairs = Triangulation(tempBaseFacList,perpBaseList).edges.tolist()

    ## The delaunay pairs do not necessarily have the indexes with lowest 
    ## first, so let's arrange and sort the delaunay pairs
    delaunayPairs = pair_sort(delaunayPairs)

    return delaunayPairs


def select_pairs_sequential(dateList, num_incr=2):
    '''Select Pairs in a Sequential way:
        For each acquisition, find its num_incr nearest acquisitions in the past time.
    Reference:
        Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series, IEEE TGRS, 51(7), 4249-4259.
    '''

    ## Get date index
    indexList = list(range(len(dateList)))

    ## Generate Pairs
    pairs = []
    for idx in indexList:
        for i in range(num_incr):
            if not idx-i-1 < 0:  pairs.append([idx-i-1,idx])
    pairs = pair_sort(pairs)

    return pairs


def select_pairs_hierarchical(tempBaseList,perpBaseList,tempPerpList):
    '''Select Pairs in a hierarchical way using list of temporal and perpendicular baseline thresholds
        For each temporal/perpendicular combination, select all possible pairs; and then merge all combination results
        together for the final output (Zhao, 2015).
    Examples:
        pairs = select_pairs_hierarchical(tempBaseList,perpBaseList,[[32, 800], [48, 600], [64, 200]])
    Reference:
        Zhao, W., (2015), Small deformation detected from InSAR time-series and their applications in geophysics, Doctoral
        dissertation, Univ. of Miami, Section 6.3.
    '''

    ## Get date index
    indexList = list(range(len(tempBaseList)))

    ## Generate All possible pairs
    pairsAll = []
    pairsAll0 = list(itertools.combinations(indexList,2))    ## Generate All Pairs
    for pair in pairsAll0:  pairsAll.append(list(pair))      ## Convert tuple object to list object

    ## Filter Pairs in a Hierarchical way through loops
    pairs = []
    for i in range(len(tempPerpList)):
        tempBaseMax = tempPerpList[i][0]
        perpBaseMax = tempPerpList[i][1]

        pairs_tmp = threshold_temporal_baseline(pairsAll,tempBaseList,tempBaseMax,0)
        pairs_tmp = threshold_perp_baseline(pairs_tmp,perpBaseList,perpBaseMax)

        pairs = pair_merge(pairs,pairs_tmp)

    pairs = pair_sort(pairs)

    return pairs


def select_pairs_mst(tempBaseList,perpBaseList,normalize=1):
    '''Select Pairs using Minimum Spanning Tree technique
        Connection Cost is calculated using the baseline distance in perp and scaled temporal baseline (Pepe and Lanari,
        2006, TGRS) plane.
    References:
        Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
        of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
        Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
    '''

    ##### Ratio between perpendicular and temporal baselines (Pepe and Lanari, 2006, TGRS)
    temp2perp_scale = (max(perpBaseList)-min(perpBaseList)) / (max(tempBaseList)-min(tempBaseList))
    tempBaseFacList = []
    for tempBase in tempBaseList:  tempBaseFacList.append(tempBase * temp2perp_scale )

    ##### Minimum Spanning Tree
    ## Prepare Input graph from temporal/perpendicular baseline
    ttMat1, ttMat2 = np.meshgrid(np.array(tempBaseFacList), np.array(tempBaseFacList))
    ppMat1, ppMat2 = np.meshgrid(np.array(perpBaseList), np.array(perpBaseList))
    ttMat = np.abs(ttMat1-ttMat2)
    ppMat = np.abs(ppMat1-ppMat2)

    #baseMat = ttMat + ppMat
    baseMat = np.sqrt(np.square(ttMat) + np.square(ppMat))    ## baseline distance in perp - scaled_temp plane.
    baseMat = csr_matrix(baseMat)

    ## Find MST
    pairMat = minimum_spanning_tree(baseMat)

    ## Convert to Pairs Index List
    idxArray2,idxArray1 = find(pairMat)[0:2]
    idxList1 = np.ndarray.tolist(idxArray1)
    idxList2 = np.ndarray.tolist(idxArray2)

    pairs = []
    for i in range(len(idxList1)):  pairs.append([idxList1[i],idxList2[i]])
    pairs = pair_sort(pairs)

    return pairs


def select_pairs_star(dateList, m_date):
    '''Select Star-like network/interferograms/pairs, it's a single master network, similar to PS approach.
    Usage:
        m_date : master date, choose it based on the following cretiria:
                 1) near the center in temporal and spatial baseline
                 2) prefer winter season than summer season for less temporal decorrelation
    Reference:
        Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
    '''

    ## Pre-process Inputs
    dateList = ptime.yymmdd(sorted(dateList))
    m_date   = ptime.yymmdd(m_date)

    ## Get date index
    idxList = list(range(len(dateList)))
    m_idx = dateList.index(m_date)

    pairs = []
    for idx in idxList:
        if not idx == m_idx:
            pairs.append([m_idx,idx])
    pairs = pair_sort(pairs)

    return pairs



################################# Plotting #################################
def plot_network(ax, pairs_idx, date8List, bperpList, plot_dict={}):
    '''Plot Temporal-Perp baseline Network
    Inputs
        ax : matplotlib axes object
        pairs_idx : list of list of 2 int, pairs index, len = number of interferograms
        date8List : list of 8-digit string, date, len=number of acquisition
        bperpList : list of float, perp baseline, len=number of acquisition
        plot_dict : dictionary with the following items:
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
    if not 'markercolor' in keyList:   plot_dict['markercolor'] = 'k'
    if not 'markersize'  in keyList:   plot_dict['markersize']  = 16
    # For colorful display of coherence
    if not 'coherence_list' in keyList:  plot_dict['coherence_list'] = None
    if not 'disp_min'       in keyList:  plot_dict['disp_min']       = 0.4
    if not 'disp_max'       in keyList:  plot_dict['disp_max']       = 1.0
    if not 'colormap'       in keyList:  plot_dict['colormap']       = 'RdBu'

    # Date Convert
    dates, datevector = ptime.date_list2vector(date8List)

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
    ax.plot(dates, bperpList, 'o', lw=plot_dict['linewidth'], alpha=0.7,\
            ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])

    # Line - Pair/Interferogram
    for i in range(len(pairs_idx)):
        idx = pairs_idx[i]
        x = np.array([dates[idx[0]], dates[idx[1]]])
        y = np.array([bperpList[idx[0]], bperpList[idx[1]]])
        if plot_dict['coherence_list']:
            coherence = plot_dict['coherence_list'][i]
            ax.plot(x, y, lw=plot_dict['linewidth'], alpha=0.7, c=cmap(coherence)) 
        else:
            ax.plot(x, y, lw=plot_dict['linewidth'], alpha=0.7, c='k')
    
    ax.set_title('Interferogram Network', fontsize=plot_dict['fontsize'])
    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])
    ax = auto_adjust_yaxis(ax, bperpList, plot_dict['fontsize'])
    ax.set_xlabel('Time [years]',fontsize=plot_dict['fontsize'])
    ax.set_ylabel('Perpendicular Baseline [m]',fontsize=plot_dict['fontsize'])

    return ax


def plot_perp_baseline_hist(ax, date8List, bperpList, plot_dict={}):
    ''' Plot Perpendicular Spatial Baseline History
    Inputs
        ax : matplotlib axes object
        date8List : list of 8-digit string, date 
        bperpList : list of float, perp baseline 
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
    dates, datevector = ptime.date_list2vector(date8List)

    # Plot
    #ax=fig.add_subplot(111)
    ax.plot(dates, bperpList, '-ko', lw=plot_dict['linewidth'], alpha=0.7,\
            ms=plot_dict['markersize'], mfc=plot_dict['markercolor'])
    ax.set_title('Perpendicular Baseline History',fontsize=plot_dict['fontsize'])

    # axis format
    ax = ptime.auto_adjust_xaxis_date(ax, datevector, plot_dict['fontsize'])
    ax = auto_adjust_yaxis(ax, bperpList, plot_dict['fontsize'])
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



