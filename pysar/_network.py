#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar._network as pnet


import h5py
import numpy as np
import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

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
def read_baseline_file(baselineFile):
    ## Read bl_list.txt and put each line into an array
    fb = open(baselineFile)
    lines = []
    for line in fb.xreadlines():
        l = str.replace(line,'\n','').strip()
        lines.append(l)
    fb.close()
  
    ## Read each line and put the values into arrays
    dates   = []
    pbase   = []
    doppler = []
    prf     = []
    SLCdirs = []
    for line in lines:
        c = line.split()                # splits on white space
        dates.append(c[0])
        pbase.append(float(c[1]))
        doppler.append(np.array([float(c[2]), float(c[3]), float(c[4])]))
        prf.append(float(c[5]))
        SLCdirs.append(c[6])
  
    dateList = ptime.yyyymmdd(dates)
  
    return dateList, pbase


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




