#!/usr/bin/env python3
import os
import gdal
import h5py
import argparse
import numpy as np
from mintpy.utils import ptime

EXAMPLE = """example:
  prep_aria.py -w mintpy_SanFran -s stack/stack/ -i stack/incidenceAngle/20150605_20150512.vrt
"""

def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ARIA processed products for MintPy.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('-w', '--work-dir', dest='workDir', type=str, default="./",
                        help='The working directory for MintPy.'
                        )
    parser.add_argument('-s', '--stack-dir', dest='stackDir', type=str, required=True,
                        help='The directory which contains stack VRT files.'
                        )
    parser.add_argument('-u', '--unwrap-stack-name', dest='unwStack', type=str,
                        default="unwrapStack.vrt",
                        help='Name of the stack VRT file of unwrapped data.\n'+
                              ' default: unwrapStack.vrt')
    parser.add_argument('-c', '--coherence-stack-name', dest='cohStack', type=str,
                        default="cohStack.vrt",
                        help='Name of the stack VRT file of coherence data.\n'+
                              'default: cohStack.vrt')
    parser.add_argument('-l', '--conn-comp-name', dest='connCompStack', type=str,
                        default="connCompStack.vrt",
                        help='Name of the stack VRT file of connected component data.\n' +
                              'default: connCompStack.vrt')
    parser.add_argument('-i', '--incidence-angle', dest='incidenceAngle', type=str,
                        required=True,
                        help='Name of the incidence angle file')
    return parser

def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

def extract_metadata(stack):

    ds = gdal.Open(stack, gdal.GA_ReadOnly)

    meta = {}
    meta["PROCESSOR"] = "isce"
    meta["FILE_LENGTH"] = ds.RasterYSize
    meta["LENGTH"] = ds.RasterYSize
    meta["ORBIT_DIRECTION"] = "DESCENDING"
    meta["PLATFORM"] = "Sen"
    meta["STARTING_RANGE"] = float(ds.GetRasterBand(1).GetMetadata("unwrappedPhase")['startRange'][1:-1])
    meta["WAVELENGTH"] = float(ds.GetRasterBand(1).GetMetadata("unwrappedPhase")['Wavelength (m)'])
    meta["WIDTH"] = ds.RasterXSize
    meta["NUMBER_OF_PAIRS"] = ds.RasterCount

    # ARIA standard products currently don't number of range and
    # azimuth looks. They are however fixed to the following values
    meta['ALOOKS'] = 7
    meta['RLOOKS'] = 19

    # geo transformation
    geoTrans = ds.GetGeoTransform()
    meta["X_FIRST"] =   geoTrans[0]
    meta["Y_FIRST"] =   geoTrans[3]
    meta["X_STEP"]  =   geoTrans[1]
    meta["Y_STEP"]  =   geoTrans[5]
    meta["X_UNIT"]  =   "degrees"
    meta["Y_UNIT"]  =   "degrees"

    # following values probably won't be used anywhere for the geocoded data
    # earth radius
    meta["EARTH_RADIUS"] = 6337286.638938101
    # nominal altitude of Sentinel1 orbit
    meta["HEIGHT"] = 693000.0

    ds = None

    return meta

def layout_hdf5(filename, dsNameDict, metadata):

    h5 = h5py.File(filename, "w")

    for key in dsNameDict.keys():
        print("creating dataset for {0}".format(key))
        h5.create_dataset(key, shape=dsNameDict[key][1], dtype=dsNameDict[key][0])

    for key in metadata.keys():
        h5.attrs[key] = metadata[key]

    h5.close()

    return


def add_unwrapped_phase(h5File, unwStack, cohStack, connCompStack):

    dsUnw = gdal.Open(unwStack, gdal.GA_ReadOnly)
    dsCoh = gdal.Open(cohStack, gdal.GA_ReadOnly)
    dsComp = gdal.Open(connCompStack, gdal.GA_ReadOnly)

    nPairs = dsUnw.RasterCount

    h5 = h5py.File(h5File, "a")

    prog_bar = ptime.progressBar(maxValue=nPairs)
    for ii in range(nPairs):
        bnd = dsUnw.GetRasterBand(ii+1)
        h5["unwrapPhase"][ii,:,:] = -1.0*bnd.ReadAsArray()

        d12 = bnd.GetMetadata("unwrappedPhase")["Dates"]
        h5["date"][ii,0] = d12.split("_")[1].encode("utf-8")
        h5["date"][ii,1] = d12.split("_")[0].encode("utf-8")

        bnd = dsCoh.GetRasterBand(ii+1)
        h5["coherence"][ii,:,:] = bnd.ReadAsArray()

        bnd = dsComp.GetRasterBand(ii+1)
        raster = bnd.ReadAsArray()
        raster[raster<0]=0   # assign pixel with no-data [-1] to zero
        h5["connectComponent"][ii,:,:] = raster

        bperp = float(dsUnw.GetRasterBand(ii+1).GetMetadata("unwrappedPhase")["perpendicularBaseline"])
        h5["bperp"][ii] = -1.0*bperp

        h5["dropIfgram"][ii] = True
        prog_bar.update(ii+1, suffix='{}_{}'.format(d12.split("_")[1], d12.split("_")[0]))

    prog_bar.close()
    h5.close()
    dsUnw = None
    dsCoh = None
    dsComp = None

    return

def add_geometry(h5File, incAngleFile):

    h5 = h5py.File(h5File, 'a')

    startRange = h5.attrs['STARTING_RANGE']

    ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
    data = ds.ReadAsArray()
    h5['incidenceAngle'][:,:] = data
    data[data!=0] = startRange
    h5['slantRangeDistance'][:,:] = data

def main(iargs=None):
    inps = cmd_line_parse(iargs)

    inps.stackDir = os.path.abspath(inps.stackDir)
    cohStack = os.path.join(inps.stackDir, inps.cohStack)
    unwStack = os.path.join(inps.stackDir, inps.unwStack)
    connCompStack = os.path.join(inps.stackDir, inps.connCompStack)

    metadata = extract_metadata(unwStack)

    length = metadata["LENGTH"]
    width = metadata["WIDTH"]
    numPairs = metadata["NUMBER_OF_PAIRS"]

    dsNameDict = {"bperp": (np.float32, (numPairs,)),
                  "coherence": (np.float32, (numPairs, length, width)),
                  "connectComponent": (np.byte, (numPairs, length, width)),
                  "date": (np.dtype('S8'), (numPairs, 2)),
                  "dropIfgram": (np.bool, (numPairs,)),
                  "unwrapPhase": (np.float32, (numPairs, length, width))}

    workDir = os.path.abspath(inps.workDir)
    if not os.path.exists(workDir):
        os.makedirs(workDir)

    inputDir = os.path.join(workDir , "inputs")
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)

    h5Filename = os.path.join(inputDir, "ifgramStack.h5")

    layout_hdf5(h5Filename, dsNameDict, metadata)

    add_unwrapped_phase(h5Filename, unwStack, cohStack, connCompStack)

    # setting up geometry:
    dsNameDict = {"azimuthAngle": (np.float32, (length, width)),
                  "height": (np.float32, (length, width)),
                  "incidenceAngle": (np.float32, (length, width)),
                  "shadowMask": (np.bool, (length, width)),
                  "slantRangeDistance": (np.float32, (length, width))}

    h5Filename = os.path.join(inputDir, "geometryGeo.h5")
    layout_hdf5(h5Filename, dsNameDict, metadata)
    add_geometry(h5Filename, inps.incidenceAngle)

if __name__=="__main__":
    main()
