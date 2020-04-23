#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, Emre Havazli,	   #
#	  Forrest Williams 2020				   #
############################################################

import os
import re
import glob
import h5py
import time
import argparse
import numpy as np
from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile, ptime, utils
from mintpy.prep_isce import read_baseline_timeseries
try:
    import gdal
except ImportError:
    raise ImportError("Can not import gdal [version>=3.0]!")


####################################################################################
EXAMPLE = """example:
prep_fringe.py -s merged/SLC -i unwrap -f *.unw -m IW*.xml -g merged/geom_master -c tcorr_ds_ps.bin -b baselines -B 4966 1145 5485 1349
"""

def create_parser():
    """Command Line Parser"""
    parser = argparse.ArgumentParser(description="Create MintPy objects from FRInGE output",
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-s", "--slc_dir", dest="slcDir", help="ISCE generated SLC directory")
    parser.add_argument("-i", "--ifg_dir", dest="ifgDir", help="unwrapped interferogram directory")
    parser.add_argument("-f", "--file", dest="file", help="unwrapped interferogram file pattern (*.unw)")
    parser.add_argument("-m", "--meta-file", dest="metafile", help="Metadata file to extract common metada for the stack")
    parser.add_argument("-g", "--geom_dir", dest="geomDir", help="geometry directory")
    parser.add_argument("-c", "--coherence", dest="cohFile", help="coherence file")
    parser.add_argument("-b", "--bperp_dir", dest="baselineDir", help="baseline text file directory")
    parser.add_argument("-B", "--bbox", dest="bbox", type=int, nargs="+", help="pixel bounding box used in FRInGE: xMin yMin xMax yMax")

    return parser

def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if all(not i for i in [inps.ifgDir, inps.geomDir, inps.cohFile, inps.baselineDir, inps.bbox]):
        parser.print_usage()
        raise SystemExit("error: all of the following arguments are required: -i, -f, -g, -c, -p, -b")
    return inps

####################################################################################
def layout_hdf5(fname, dsNameDict, metadata):
    print("-"*50)
    print("create HDF5 file {} with w mode".format(fname))
    h5 = h5py.File(fname, "w")

    # initiate dataset
    for key in dsNameDict.keys():
        compression = None

        # changable dataset shape
        if len(dsNameDict[key][1]) == 3:
            maxShape = (None, dsNameDict[key][1][1], dsNameDict[key][1][2])
        else:
            maxShape = dsNameDict[key][1]

        print("create dataset: {d:<25} of {t:<25} in size of {s}".format(
            d=key,
            t=str(dsNameDict[key][0]),
            s=dsNameDict[key][1]))
        h5.create_dataset(key,
                          shape=dsNameDict[key][1],
                          maxshape=maxShape,
                          dtype=dsNameDict[key][0],
                          chunks=True,
                          compression=compression)

    # write attributes
    for key in metadata.keys():
        h5.attrs[key] = metadata[key]

    h5.close()
    print("close HDF5 file {}".format(fname))
    return

def write_geometry(outFile, geomDir, bbox):
    f = h5py.File(outFile, "a")
    
    hgt = os.path.join(geomDir, "hgt.rdr.full")
    lat = os.path.join(geomDir, "lat.rdr.full")
    lon = os.path.join(geomDir, "lon.rdr.full")
    los = os.path.join(geomDir, "los.rdr.full")
    shdw = os.path.join(geomDir, "shadowMask.rdr.full")
    
    geomDict = {"height":hgt, "latitude":lat, "longitude":lon,
                   "incidenceAngle":los, "azimuthAngle":los,"shadowMask":shdw}
    
    for key in geomDict:
        ds = gdal.Open(geomDict[key], gdal.GA_ReadOnly)
        
        if key == "ShadowMaxk":
            data = np.array(ds.ReadAsArray(), dtype=np.byte)
        else:
            data = np.array(ds.ReadAsArray(), dtype=np.float32)
        
        data[data == ds.GetRasterBand(1).GetNoDataValue()] = np.nan 

        if key == "incidenceAngle":
            data = data[0]
        elif key == "azimuthAngle":
           data = data[1]
        
        f[key][:,:] = data[bbox[1]:bbox[3], bbox[0]:bbox[2]]
    
    meta = dict(f.attrs)
    meta["LENGTH"], meta["WIDTH"] = data.shape
    slantRangeFull = utils.range_distance(meta, dimension=2)
    f["slantRangeDistance"][:,:] = slantRangeFull[bbox[1]:bbox[3], bbox[0]:bbox[2]]
    f.close()
    
    return outFile

def write_timeseries(outFile, ifgs, baselineDir):
    #Prepare timeseries data
    f = h5py.File(outFile, "a")
    
    bperpDict = read_baseline_timeseries(baselineDir, "tops")
    
    dates = list(bperpDict.keys())
    dates.sort()
    f["date"][:,] = [np.string_(x) for x in dates]
    
    bperps = [bperpDict[x][0] for x in dates]
    f["bperp"][:,] = bperps
        
    #Add data to arrays
    for i in range(0, len(ifgs)):
        ifg = gdal.Open(ifgs[i])
        array = np.array(ifg.GetRasterBand(2).ReadAsArray())
        f["timeseries"][i+1] = array
    
    f.close()
    return

def write_temporalcoherence(outFile, cohFile):
    f = h5py.File(outFile, "a")
    
    cohGdal = gdal.Open(cohFile)
    cohArray = np.array(cohGdal.GetRasterBand(1).ReadAsArray(), dtype = "float32")
    cohArray[cohArray < 0] = 0
    f["temporalCoherence"][:,:] = cohArray

    return

####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    
    # find interferograms
    ifgs = glob.glob(os.path.join(inps.ifgDir, inps.file))
    
    # prepare metadata
    script_name = "prep_isce.py"
    cmd = script_name+" -d {} -f *.slc.full -m {} -b {} -g {} --force".format(inps.slcDir, inps.metafile, inps.baselineDir, inps.geomDir)
    os.system(cmd)
    
    metadata = readfile.read_attribute(ifgs[0])
    width = inps.bbox[2] - inps.bbox[0]
    length = inps.bbox[3] - inps.bbox[1]
    metadata["LENGTH"] = length
    metadata["WIDTH"] = width
    numPairs = len(ifgs)
    numDates = numPairs + 1
    metadata["NUMBER_OF_PAIRS"] = numPairs
    metadata["REF_DATE"] = metadata["startUTC"][0:10].replace("-", "") #could also be done using datetime
    

    # prepare output directory
    inputDir = os.path.join("inputs")
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)

    # 1. geometryGeo
    # define dataset structure
    dsNameDict = {
        "height": (np.float32, (length, width)),
        "latitude": (np.float32, (length, width)),
        "longitude": (np.float32, (length, width)),
        "incidenceAngle": (np.float32, (length, width)),
        "azimuthAngle": (np.float32, (length, width)),
        "shadowMask": (np.byte, (length, width)),
        "slantRangeDistance": (np.float32, (length, width))
    }

    # write data to disk
    geom_file = os.path.join(inputDir, "geometryRadar.h5")
    metadata["FILE_TYPE"] = "geometry"
    layout_hdf5(geom_file, dsNameDict, metadata)
    write_geometry(geom_file,inps.geomDir,inps.bbox)

    # 2. timeseries
    # define dataset structure
    dsNameDict = {
        "bperp": (np.float32, (numDates,)),
        "date": (np.dtype("S8"), (numDates,)),
        "timeseries": (np.float32, (numDates, length, width))
    }

    # write data to disk
    timeseries_file = "timeseries.h5"
    metadata["FILE_TYPE"] = "timeseries"
    layout_hdf5(timeseries_file, dsNameDict, metadata)
    write_timeseries(timeseries_file, ifgs, inps.baselineDir)

    # 3. temporalCoherence
    # define dataset structure
    dsNameDict = {
        "temporalCoherence": (np.float32, (length, width))
    }

    # write data to disk
    coherence_file = "temporalCoherence.h5"
    metadata["FILE_TYPE"] = "temporalCoherence"
    layout_hdf5(coherence_file, dsNameDict, metadata)
    write_temporalcoherence(coherence_file, inps.cohFile)
    
    return

####################################################################################
if __name__=="__main__":
    main()