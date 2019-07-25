#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2019    Yunmeng Cao                         #
# Author:  Yunmeng Cao                                     #
############################################################

import os
import sys
import glob
import time
import argparse

import h5py
import numpy as np
from mintpy.utils import readfile, writefile, ptime, utils as ut
#########################################################################

INTRODUCTION = '''
#############################################################################
   
   Convert the Geo-coordinates based lookup-table (GAMMA, ROI_PAC) into the Radar-coordinates based lookup-table (ISCE).
'''

EXAMPLE = '''
    Usage:
            generate_rdc_lt.py geometryGeo.h5 
            generate_rdc_lt.py geometryGeo.h5 -w geometryRadar.h5 
            generate_rdc_lt.py geometryGeo.h5 -w geometryRadar.h5 -n 2

##############################################################################
'''
def write_h5(datasetDict, out_file, metadata=None, ref_file=None, compression=None):
    #output = 'variogramStack.h5'
    'lags                  1 x N '
    'semivariance          M x N '
    'sills                 M x 1 '
    'ranges                M x 1 '
    'nuggets               M x 1 '
    
    if os.path.isfile(out_file):
        print('delete exsited file: {}'.format(out_file))
        os.remove(out_file)

    print('create HDF5 file: {} with w mode'.format(out_file))
    dt = h5py.special_dtype(vlen=np.dtype('float64'))

    
    with h5py.File(out_file, 'w') as f:
        for dsName in datasetDict.keys():
            data = datasetDict[dsName]
            ds = f.create_dataset(dsName,
                              data=data,
                              compression=compression)
        
        for key, value in metadata.items():
            f.attrs[key] = str(value)
            #print(key + ': ' +  value)
    print('finished writing to {}'.format(out_file))
        
    return out_file  

def get_idx_cpx(RangeCoord,AzimuthCoord,WIDTH,LENGTH,n_lines):
    
    NN = n_lines # search n lines one time
    x = np.arange(0,WIDTH)+1
    y = np.arange(0,LENGTH)+1
    xv, yv = np.meshgrid(x, y)
    yv = yv.flatten()
    xv = xv.flatten()
    radar_cpx = xv + yv*1j
    #print(len(radar_cpx))

    
    Ns = len(xv)
    short_idx = np.zeros((Ns,),dtype=int)
    
    
    gx = RangeCoord.flatten()
    gy = AzimuthCoord.flatten()
    gy[gy==0] = 10000
    gy[gy<0] = 10000
    geo_cpx = gx + gy*1j
    
    sort_gy_idx = np.argsort(gy)
    sort_gy_array = gy[sort_gy_idx]
    sort_gx_array = gx[sort_gy_idx]
    sort_geo_cpx = sort_gx_array + sort_gy_array*1j
    #print(sort_gy_array[100:200])
    k = round(LENGTH/NN)+1
    
    prog_bar = ptime.progressBar(maxValue=k)
    for i in range(k):
        i = i +1
        a0 = NN*(i-1)
        if NN*i < LENGTH:
            b0 = NN*i
        else:
            b0 = LENGTH
        
        if not a0 > b0:
            idx0 = a0*WIDTH
            idx1 = b0*WIDTH
            radar_cpx0 = radar_cpx[idx0:idx1]
            #print(a0-0.5)
            #print(b0+0.5)
            #idx_geo0 = np.where((a0-0.5)<sort_gy_array & sort_gy_array < (b0+0.5)
            idx_geo = np.where(((a0-0.5)<sort_gy_array)*(sort_gy_array < (b0+0.5)))
            sort_geo_cpx0 = sort_geo_cpx[idx_geo]
            
            #print(len(sort_geo_cpx0))
            Ns0 = len(radar_cpx0)
            #print(Ns0)
            short_idx0 = np.zeros((Ns0,))
            sort_idx0 = sort_gy_idx[idx_geo]
            
            
            for j in range(Ns0):
                id0 = find_nearest_cpx(sort_geo_cpx0, radar_cpx0[j])
                short_idx0[j] = sort_idx0[id0]
            #print(short_idx0)
            short_idx[idx0:idx1] = short_idx0
        prog_bar.update(i+1, every=round(k/100), suffix='{}/{} lines'.format((i+1)*NN, LENGTH))
    prog_bar.close()
    
    return short_idx

def find_nearest_cpx(cpx_long, cpx0):
    cpx_long = np.asarray(cpx_long)
    idx =  (np.abs(cpx_long - cpx0)).argmin()
    return idx 


def get_dataNames(FILE):
    with h5py.File(FILE, 'r') as f:
        dataNames = []
        for k0 in f.keys():
            dataNames.append(k0)
    return dataNames


def find_nearest(array, value):
    array = np.asarray(array,dtype=np.float64)
    idx = (np.abs(array - value)).argmin()
    return idx

def get_idx(long_array,short_array,n=2):
        sort_idx = np.argsort(long_array)
        long_array_sort = long_array[sort_idx]
        #print(long_array_sort[0:100])
        
        Ns = len(short_array)
        k = round(Ns/n)+1
        short_idx = np.zeros((Ns,),dtype=int)
        
        prog_bar = ptime.progressBar(maxValue=n)
        for i in range(n):
            i=i+1
            a0 = k*(i-1)
            if k*(i+1) < Ns:
                b0 = k*(i+1)
            else:
                b0 = Ns
            
            if not a0 > b0:
                idx0 = np.arange(a0,b0)
                #print(idx0)
                short_array0 = short_array[idx0]
                Ns0 = len(short_array0)
            
                idx_value = find_nearest(long_array_sort, short_array0[0])
                long_array0 = long_array_sort[idx_value:(idx_value+2*k)]
                short_idx0 = np.zeros((Ns0,))
                sort_idx0 = sort_idx[idx_value:(idx_value+2*k)]
            
                for j in range(Ns0):
                    id0 = find_nearest(long_array0, short_array0[j])
                    short_idx0[j] = sort_idx0[id0]
            
                short_idx[idx0] = short_idx0
            prog_bar.update(i+1, every=round(n/100), suffix='{}/{} pixels'.format(i+1, n))
        prog_bar.close() 
        return short_idx

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Generate Radar-coordinates based lookup-table.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=INTRODUCTION+'\n'+EXAMPLE)

    parser.add_argument('geometryGeo',help='geometryGeo file which includes geo-coordinates based lookup-table')
    parser.add_argument('-w','--write', dest='write', metavar='FILE', default = 'geometryRadar.h5',
                      help='update geometryRadar.h5 file by adding the radar-coordinates based lookup-table.')
    parser.add_argument('-n','--numb', dest='search_line_numb', type=int, metavar='NUM',default = 2,
                      help='the lines used as the length of the search pool.')
    
    inps = parser.parse_args()

    return inps

################################################################################    
    
    
def main(argv):
    
    inps = cmdLineParse() 
    geom = inps.geometryGeo
    rangeCoord = readfile.read(geom,datasetName = 'rangeCoord')[0]
    azimuthCoord = readfile.read(geom,datasetName = 'azimuthCoord')[0]
    Search_Linear_Number = inps.search_line_numb
    rangeCoord = rangeCoord.astype(np.float64)
    azimuthCoord = azimuthCoord.astype(np.float64)
    #CPX_lt =complex(rangeCoord + '+' + azimuthCoord+'j')
    CPX_lt = rangeCoord  + 1j *azimuthCoord
    
    meta_geo = readfile.read_attribute(geom)
    post_Lat = meta_geo['Y_STEP']
    post_Lon = meta_geo['X_STEP']
    Corner_LAT = meta_geo['Y_FIRST']
    Corner_LON = meta_geo['X_FIRST']
            
    if inps.write:
        meta = readfile.read_attribute(inps.write)
    elif inps.reference:
        meta = readfile.read_attribute(inps.reference)
    else:
        print('write_file or the reference_file should be provided at least one.')
        sys.exit(1)
    
    WIDTH_geo  = int(meta_geo['WIDTH'])
    LENGTH_geo  = int(meta_geo['LENGTH'])
    
    x = np.arange(0,WIDTH_geo)
    y = np.arange(0,LENGTH_geo)
    xv, yv = np.meshgrid(x, y)
    
    LAT = float(Corner_LAT) + yv*float(post_Lat)
    LON = float(Corner_LON) + xv*float(post_Lon)
    LAT = LAT.flatten()
    LON = LON.flatten()
    #print(len(LAT))   
        
    WIDTH  = int(meta['WIDTH'])
    LENGTH  = int(meta['LENGTH'])
        
    Radar = np.arange(0,WIDTH*LENGTH) +1
    Radar = Radar.reshape(LENGTH,WIDTH)
    Radar = Radar.flatten()
    
    #IDX = get_idx(Rdc_IDX,Radar,n=4000)
    IDX = get_idx_cpx(rangeCoord,azimuthCoord,WIDTH,LENGTH,Search_Linear_Number) #every two lines update the search pool.
    
    lat_sar = LAT[IDX]
    lon_sar = LON[IDX]
    
    lat_sar = lat_sar.reshape(LENGTH,WIDTH)
    lon_sar = lon_sar.reshape(LENGTH,WIDTH)
    
    dataNames = get_dataNames(inps.write)
    datasetDict = dict()
    meta = readfile.read_attribute(inps.write)
    
    for k0 in dataNames:
        datasetDict[k0] = readfile.read(inps.write,datasetName = k0)[0]
     
    DEM  = readfile.read(inps.write,datasetName = 'height')[0]
    lat_sar[DEM==0] = 0
    lon_sar[DEM==0] = 0
    datasetDict['latitude'] = lat_sar.astype(np.float32)
    datasetDict['longitude'] = lon_sar.astype(np.float32)
    write_h5(datasetDict, inps.write, metadata=meta, ref_file=None, compression=None)
    
    print('done.')
    
    sys.exit(1)
##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
