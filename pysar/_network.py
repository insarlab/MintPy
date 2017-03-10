#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar._network as pnet
#


import h5py
import numpy as np
import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import itertools

import pysar._datetime as ptime


############################################################
def read_pairs_list(listFile,dateList):
    ## Read Pairs List file like below:
    ## 070311-070426
    ## 070311-070611
    ## ...

    dateList6 = ptime.yymmdd(dateList)
    pairs=[]
    fl = open(listFile,'r')
    lines = []
    lines = fl.read().splitlines()
    for line in lines:
        date12 = line.split('-')
        pairs.append([dateList6.index(date12[0]),dateList6.index(date12[1])])
    fl.close()

    pairs = pair_sort(pairs)

    return pairs


############################################################
def write_pairs_list(pairs,dateList,outName):
    dateList6 = ptime.yymmdd(dateList)
    fl = open(outName,'w')
    for idx in pairs:
        date12 = dateList6[idx[0]]+'-'+dateList6[idx[1]]+'\n'
        fl.write(date12)
    fl.close()
    return 1


############################################################
def read_igram_pairs(igramFile):
    ## Read Igram file
    h5file = h5py.File(igramFile,'r')
    k = h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    if k[0] not in  ['interferograms','coherence','wrapped']:
        print 'Only interferograms / coherence / wrapped are supported.';  sys.exit(1)

    dateList  = ptime.date_list(igramFile)
    dateList6 = ptime.yymmdd(dateList)

    pairs = []
    igramList=h5file[k[0]].keys()
    for igram in igramList:
        date12 = h5file[k[0]][igram].attrs['DATE12'].split('-')
        pairs.append([dateList6.index(date12[0]),dateList6.index(date12[1])])
    h5file.close()

    pairs = pair_sort(pairs)

    return pairs


#############################################################
def read_baseline_file(baselineFile,exDateList=[]):
    ## Read bl_list.txt without dates listed in exDateList
    ##
    ## Examples:
    ##     date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile)
    ##     date8List, perpBaseList, dopList, prfList, slcDirList = read_baseline_file(baselineFile,['080520','100726'])
    ##     date8List, perpBaseList = read_baseline_file(baselineFile)[0:2]

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

############################################################
def threshold_perp_baseline(igramIdxList,perpBaseList,perpBaseMax=800,perpBaseMin=0):
    ## Remove pairs/interoferogram that have perpendicular baseline out of [perpBaseMin, perpBaseMax]
    ##
    ## Example:
    ##     pairs = threshold_perp_baseline(pairs,perpBaseList,500)
    ##

    igramIdxListOut = []
    for igramIdx in igramIdxList:
        idx1 = igramIdx[0]
        idx2 = igramIdx[1]
        perpBase = abs(perpBaseList[idx1]-perpBaseList[idx2])
        if perpBaseMin <= perpBase <= perpBaseMax:
            igramIdxListOut.append(igramIdx)

    return igramIdxListOut


############################################################
def threshold_temporal_baseline(igramIdxList,tempBaseList,tempBaseMax=365,seasonal=1,tempBaseMin=0):
    ## Remove pairs/interferograms out of min/max/seasonal temporal baseline limits
    ##
    ## Usage:
    ##     seasonal : keep interferograms with seasonal temporal baseline
    ##                1 - keep them, by default
    ##                0 - do not keep them
    ##
    ## Example:
    ##     pairs = threshold_temporal_baseline(pairs,tempBaseList,80)
    ##     pairs = threshold_temporal_baseline(pairs,tempBaseList,80,0)  # disable seasonal checking

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


############################################################
def pair_sort(pairs):
    for idx in range(len(pairs)):
        if pairs[idx][0] > pairs[idx][1]:
            index1 = pairs[idx][1]
            index2 = pairs[idx][0]
            pairs[idx][0] = index1
            pairs[idx][1] = index2
    pairs=sorted(pairs)
    return pairs


############################################################
def pair_merge(pairs1,pairs2):
    pairs = pairs1
    for pair in pairs2:
        if not pair in pairs:
            pairs.append(pair)

    pairs=sorted(pairs)
    return pairs

############################################################
def select_pairs_all(dateList):
    ## Select All Possible Pairs/Interferograms
    ##
    ## Reference:
    ##     Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
    ## based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
    ##

    ## Get date index
    indexList = []
    for ni in range(len(dateList)):  indexList.append(ni)

    ## Generate All Pairs
    allPairs = list(itertools.combinations(indexList,2))

    ## Convert tuple object to list object
    allPairs2=[]
    for pair in allPairs: allPairs2.append(list(pair))

    return allPairs2


############################################################
def select_pairs_delaunay(tempBaseList,perpBaseList,normalize=1):
    ## Select Pairs using Delaunay Triangulation based on temporal/perpendicular baselines
    ##
    ## Usage:
    ##     tempBaseList : list of temporal baseline
    ##     perpBaseList : list of perpendicular spatial baseline
    ##     normalize    : normalize temporal baseline to perpendicular baseline
    ##                    1 - enable  normalization, default
    ##                    0 - disable normalization
    ##
    ## Key points
    ##     1. Define a ratio between perpendicular and temporal baseline axis units (Pepe and Lanari, 2006, TGRS).
    ##     2. Pairs with too large perpendicular / temporal baseline or Doppler centroid difference should be removed
    ##        after this, using a threshold, to avoid strong decorrelations (Zebker and Villasenor, 1992, TGRS).
    ##
    ## Reference:
    ##     Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    ## of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
    ##     Zebker, H. A., and J. Villasenor (1992), Decorrelation in interferometric radar echoes, IEEE TGRS, 30(5), 950-959.
    ##


    import matplotlib.delaunay as md

    ##### Generate Delaunay Triangulation based on temporal and spatial perpendicular baselines
    if normalize == 0:
        centers,edges,tri,neighbors = md.delaunay(tempBase,perpBase)
    else:
        ##### Ratio between perpendicular and temporal baselines (Pepe and Lanari, 2006, TGRS)
        tempBaseFacList = []
        temp2perp_scale = (max(perpBaseList)-min(perpBaseList)) / (max(tempBaseList)-min(tempBaseList))
        for tempBase in tempBaseList:
            #tempBaseFac = (tempBase - min(tempBaseList)) * temp2perp_scale + min(perpBaseList)
            tempBaseFac = tempBase * temp2perp_scale   # giving same result as the line above
            tempBaseFacList.append(tempBaseFac)
        centers,edges,tri,neighbors = md.delaunay(tempBaseFacList,perpBaseList)

    ## The delaunay pairs do not necessarily have the indexes with lowest 
    ## first, so let's arrange and sort the delaunay pairs
    delaunayPairs = edges.tolist()
    delaunayPairs = pair_sort(delaunayPairs)

    return delaunayPairs


############################################################
def select_pairs_sequential(dateList, num_incr=2):
    ## Select Pairs in a Sequential way:
    ##     For each acquisition, find its num_incr nearest acquisitions in the past time.
    ##
    ## Reference:
    ##     Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series, IEEE TGRS, 51(7), 4249-4259.
    ##

    ## Get date index
    indexList = list(range(len(dateList)))

    ## Generate Pairs
    pairs = []
    for idx in indexList:
        for i in range(num_incr):
            if not idx-i-1 < 0:  pairs.append([idx-i-1,idx])

    pairs = pair_sort(pairs)

    return pairs


############################################################
def select_pairs_hierarchical(tempBaseList,perpBaseList,tempPerpList):
    ## Select Pairs in a hierarchical way using list of temporal and perpendicular baseline thresholds
    ##     For each temporal/perpendicular combination, select all possible pairs; and then merge all combination results
    ##     together for the final output (Zhao, 2015).
    ##
    ## Examples:
    ##     pairs = select_pairs_hierarchical(tempBaseList,perpBaseList,[[32, 800], [48, 600], [64, 200]])
    ##
    ## Reference:
    ##     Zhao, W., (2015), Small deformation detected from InSAR time-series and their applications in geophysics, Doctoral
    ## dissertation, Univ. of Miami, Section 6.3.
    ##

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


############################################################
def select_pairs_mst(tempBaseList,perpBaseList,normalize=1):
    ## Select Pairs using Minimum Spanning Tree technique
    ##     Connection Cost is calculated using the baseline distance in perp and scaled temporal baseline (Pepe and Lanari,
    ##     2006, TGRS) plane.
    ##
    ## References:
    ##     Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    ## of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
    ##     Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
    ##

    from scipy.sparse import csr_matrix, find
    from scipy.sparse.csgraph import minimum_spanning_tree

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


def select_pairs_star(dateList,m_date):
    ## Select Star-like network/interferograms/pairs, it's a single master network, similar to PS approach.
    ##
    ## Usage:
    ##     m_date : master date, choose it based on the following cretiria:
    ##              1) near the center in temporal and spatial baseline
    ##              2) prefer winter season than summer season for less temporal decorrelation
    ##
    ## Reference:
    ##     Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
    ##

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








