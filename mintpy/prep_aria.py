#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2019, Heresh Fattahi, Zhang Yunjun          #
# Author: Heresh Fattahi, Zhang Yunjun, Jul 2019           #
############################################################


import os
import time
import argparse
import h5py
import numpy as np
from mintpy.utils import ptime
try:
    import gdal
except ImportError:
    raise ImportError('Can not import gdal [version>=3.0]!')


####################################################################################
EXAMPLE = """example:
  # run ARIA-tools to download and extract inteferograms stack for time-series analysis.
  # reference: https://github.com/aria-tools/ARIA-tools
  conda activate ARIA-tools
  ariaDownload.py -b '37.25 38.1 -122.6 -121.75' --track 42
  ariaTSsetup.py -f 'products/*.nc' -b '37.25 38.1 -122.6 -121.75' --mask Download

  prep_aria.py -w mintpy -s stack/ -d DEM/SRTM_3arcsec.dem -i incidenceAngle/20150605_20150512.vrt 
  prep_aria.py -w mintpy -s stack/ -d DEM/SRTM_3arcsec.dem -i incidenceAngle/20150605_20150512.vrt --water-mask mask/watermask.msk
"""

def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ARIA processed products for MintPy.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('-w','--work-dir', dest='workDir', type=str, default="./",
                        help='The working directory for MintPy.')

    stack = parser.add_argument_group('Interferogram stack')
    stack.add_argument('-s','--stack-dir', dest='stackDir', type=str, required=True,
                       help='The directory which contains stack VRT files.')

    stack.add_argument('-u','--unwrap-stack-name', dest='unwStack', type=str,
                       default="unwrapStack.vrt",
                       help='Name of the stack VRT file of unwrapped data.\n'+
                            'default: unwrapStack.vrt')
    stack.add_argument('-c','--coherence-stack-name', dest='cohStack', type=str,
                       default="cohStack.vrt",
                       help='Name of the stack VRT file of coherence data.\n'+
                            'default: cohStack.vrt')
    stack.add_argument('-l','--conn-comp-name', dest='connCompStack', type=str,
                       default="connCompStack.vrt",
                       help='Name of the stack VRT file of connected component data.\n' +
                            'default: connCompStack.vrt')

    geom = parser.add_argument_group('Geometry')
    geom.add_argument('-d','--dem', dest='dem', type=str, required=True,
                      help='Name of the DEM file')
    geom.add_argument('-i','--incidence-angle', dest='incAngle', type=str, required=True,
                      help='Name of the incidence angle file')
    geom.add_argument('-a','--az-angle','--azimuth-angle', dest='azAngle', type=str,
                      help='Name of the azimuth angle file. [not recommend for now]')
    geom.add_argument('--water-mask', dest='waterMask', type=str,
                      help='Name of the water mask file')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input dir and files
    inps.stackDir = os.path.abspath(inps.stackDir)
    inps.cohStack = os.path.join(inps.stackDir, inps.cohStack)
    inps.unwStack = os.path.join(inps.stackDir, inps.unwStack)
    inps.connCompStack = os.path.join(inps.stackDir, inps.connCompStack)
    for vrtFile in [inps.cohStack, inps.unwStack, inps.connCompStack]:
        if not os.path.isfile(vrtFile):
            raise FileNotFoundError(vrtFile)

    return inps


####################################################################################
def extract_metadata(stack):

    print('extract metadata from {}'.format(stack))
    ds = gdal.Open(stack, gdal.GA_ReadOnly)
    metaUnw = ds.GetRasterBand(1).GetMetadata("unwrappedPhase")

    # copy over all metadata from unwrappedPhase
    meta = {}
    for key, value in metaUnw.items():
        if key not in ["Dates", "perpendicularBaseline"]:
            meta[key] = value

    meta["ANTENNA_SIDE"] = -1
    meta["PROCESSOR"] = "isce"
    meta["FILE_LENGTH"] = ds.RasterYSize
    meta["LENGTH"] = ds.RasterYSize
    meta["ORBIT_DIRECTION"] = meta["orbitDirection"].upper()
    meta["PLATFORM"] = "Sen"
    meta["WAVELENGTH"] = float(meta["Wavelength (m)"])
    meta["WIDTH"] = ds.RasterXSize
    meta["NUMBER_OF_PAIRS"] = ds.RasterCount

    # Note from YZ, 2019-07-25
    # convert isce azimuth angle to roipac orbit heading angle
    # This value is not consistent with band2 of los.rdr from ISCE/topsStack
    # need to check with ARIA-tools team. 
    # use hardwired value for now
    #az_angle = float(meta["azimuthAngle"])
    #head_angle = -1 * (270 + az_angle)
    #head_angle -= np.round(head_angle / 360.) * 360.
    #meta['HEADING'] = head_angle
    if meta["ORBIT_DIRECTION"].startswith("D"):
        meta["HEADING"] = -168
    else:
        meta["HEADING"] = -12

    # ARIA standard products currently don't have number of range and
    # azimuth looks. They are however fixed to the following values
    meta["ALOOKS"] = 7
    meta["RLOOKS"] = 19
    meta["RANGE_PIXEL_SIZE"] = float(meta["slantRangeSpacing"]) * meta["RLOOKS"]
    meta["STARTING_RANGE"] = float(meta["startRange"])

    # geo transformation
    geoTrans = ds.GetGeoTransform()
    lon0 = geoTrans[0]
    lat0 = geoTrans[3]
    lon_step = geoTrans[1]
    lat_step = geoTrans[5]
    lon1 = lon0 + lon_step * meta["WIDTH"]
    lat1 = lat0 + lat_step * meta["LENGTH"]
    meta["X_FIRST"] = lon0
    meta["Y_FIRST"] = lat0
    meta["X_STEP"] = lon_step
    meta["Y_STEP"] = lat_step
    meta["X_UNIT"] = "degrees"
    meta["Y_UNIT"] = "degrees"

    utc = meta["UTCTime (HH:MM:SS.ss)"]
    utc = time.strptime(utc, "%H:%M:%S.%f")
    meta["CENTER_LINE_UTC"] = utc.tm_hour*3600.0 + utc.tm_min*60.0 + utc.tm_sec

    # following values probably won't be used anywhere for the geocoded data
    # earth radius
    meta["EARTH_RADIUS"] = 6337286.638938101
    # nominal altitude of Sentinel1 orbit
    meta["HEIGHT"] = 693000.0

    meta["LON_REF1"] = lon0
    meta["LON_REF2"] = lon1
    meta["LON_REF3"] = lon0
    meta["LON_REF4"] = lon1

    meta["LAT_REF1"] = lat0
    meta["LAT_REF2"] = lat0
    meta["LAT_REF3"] = lat1
    meta["LAT_REF4"] = lat1

    ds = None
    return meta


def layout_hdf5(fname, dsNameDict, metadata):
    print('-'*50)
    print('create HDF5 file {} with w mode'.format(fname))
    h5 = h5py.File(fname, "w")

    # initiate dataset
    for key in dsNameDict.keys():
        # turn ON compression
        if key in ['connectComponent']:
            compression = 'lzf'
        else:
            compression = None

        # changable dataset shape
        if len(dsNameDict[key][1]) == 3:
            maxShape = (None, dsNameDict[key][1][1], dsNameDict[key][1][2])
        else:
            maxShape = dsNameDict[key][1]

        print("create dataset: {}".format(key))
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
    print('close HDF5 file {}'.format(fname))
    return


def write_geometry(h5File, demFile, incAngleFile, azAngleFile=None, waterMaskFile=None):
    print('-'*50)
    print('writing data to HDF5 file {} with a mode ...'.format(h5File))
    h5 = h5py.File(h5File, 'a')

    # height
    ds = gdal.Open(demFile, gdal.GA_ReadOnly)
    data = np.array(ds.ReadAsArray(), dtype=np.float32)
    data[data == ds.GetRasterBand(1).GetNoDataValue()] = np.nan
    h5['height'][:,:] = data

    # slantRangeDistance
    h5['slantRangeDistance'][:,:] = float(h5.attrs['STARTING_RANGE'])

    # incidenceAngle
    ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
    data = ds.ReadAsArray()
    data[data == ds.GetRasterBand(1).GetNoDataValue()] = np.nan
    h5['incidenceAngle'][:,:] = data

    # azimuthAngle
    if azAngleFile is not None:
        ds = gdal.Open(azAngleFile, gdal.GA_ReadOnly)
        data = ds.ReadAsArray()
        data[data == ds.GetRasterBand(1).GetNoDataValue()] = np.nan
        h5['azimuthAngle'][:,:] = data

    # waterMask
    if waterMaskFile is not None:
        ds = gdal.Open(waterMaskFile, gdal.GA_ReadOnly)
        water_mask = ds.ReadAsArray()
        water_mask[water_mask == ds.GetRasterBand(1).GetNoDataValue()] = False

        # assign False to invalid pixels based on incAngle data
        ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
        data = ds.ReadAsArray()
        water_mask[data == ds.GetRasterBand(1).GetNoDataValue()] = False
        h5['waterMask'][:,:] = water_mask

    h5.close()
    print('finished writing to HD5 file: {}'.format(h5File))
    return h5File


def write_ifgram_stack(h5File, unwStack, cohStack, connCompStack):

    print('-'*50)
    print('opening {}, {}, {} with gdal ...'.format(os.path.basename(unwStack),
                                                    os.path.basename(cohStack),
                                                    os.path.basename(connCompStack)))
    dsUnw = gdal.Open(unwStack, gdal.GA_ReadOnly)
    dsCoh = gdal.Open(cohStack, gdal.GA_ReadOnly)
    dsComp = gdal.Open(connCompStack, gdal.GA_ReadOnly)

    # extract NoDataValue (from the 2nd */date2_date1.vrt file for example)
    ds = gdal.Open(dsUnw.GetFileList()[2], gdal.GA_ReadOnly)
    noDataValueUnw = ds.GetRasterBand(1).GetNoDataValue()
    print('grab NoDataValue for unwrapPhase:      {:<5} and convert to 0.'.format(noDataValueUnw))

    ds = gdal.Open(dsCoh.GetFileList()[2], gdal.GA_ReadOnly)
    noDataValueCoh = ds.GetRasterBand(1).GetNoDataValue()
    print('grab NoDataValue for coherence:        {:<5} and convert to 0.'.format(noDataValueCoh))

    ds = gdal.Open(dsComp.GetFileList()[2], gdal.GA_ReadOnly)
    noDataValueComp = ds.GetRasterBand(1).GetNoDataValue()
    print('grab NoDataValue for connectComponent: {:<5} and convert to 0.'.format(noDataValueComp))
    ds = None

    # sort the order of interferograms based on date1_date2 with date1 < date2
    nPairs = dsUnw.RasterCount
    d12BandDict = {}
    for ii in range(nPairs):
        bnd = dsUnw.GetRasterBand(ii+1)
        d12 = bnd.GetMetadata("unwrappedPhase")["Dates"]
        d12 = sorted(d12.split("_"))
        d12 = '{}_{}'.format(d12[0], d12[1])
        d12BandDict[d12] = ii+1
    d12List = sorted(d12BandDict.keys())

    print('writing data to HDF5 file {} with a mode ...'.format(h5File))
    h5 = h5py.File(h5File, "a")

    prog_bar = ptime.progressBar(maxValue=nPairs)
    for ii in range(nPairs):
        d12 = d12List[ii]
        bndIdx = d12BandDict[d12]

        h5["date"][ii,0] = d12.split("_")[0].encode("utf-8")
        h5["date"][ii,1] = d12.split("_")[1].encode("utf-8")
        h5["dropIfgram"][ii] = True

        bnd = dsUnw.GetRasterBand(bndIdx)
        data = bnd.ReadAsArray()
        data[data == noDataValueUnw] = 0      #assign pixel with no-data to 0
        h5["unwrapPhase"][ii,:,:] = -1.0*data #date2_date1 -> date1_date2

        bperp = float(bnd.GetMetadata("unwrappedPhase")["perpendicularBaseline"])
        h5["bperp"][ii] = -1.0*bperp          #date2_date1 -> date1_date2

        bnd = dsCoh.GetRasterBand(bndIdx)
        data = bnd.ReadAsArray()
        data[data == noDataValueCoh] = 0      #assign pixel with no-data to 0
        h5["coherence"][ii,:,:] = data

        bnd = dsComp.GetRasterBand(bndIdx)
        data = bnd.ReadAsArray()
        data[data == noDataValueComp] = 0     #assign pixel with no-data to 0
        h5["connectComponent"][ii,:,:] = data

        prog_bar.update(ii+1, suffix='{}'.format(d12))
    prog_bar.close()

    # add MODIFICATION_TIME metadata to each 3D dataset
    h5["unwrapPhase"].attrs['MODIFICATION_TIME'] = str(time.time())
    h5["coherence"].attrs['MODIFICATION_TIME'] = str(time.time())
    h5["connectComponent"].attrs['MODIFICATION_TIME'] = str(time.time())

    h5.close()
    print('finished writing to HD5 file: {}'.format(h5File))
    dsUnw = None
    dsCoh = None
    dsComp = None
    return h5File


####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # prepare metadata
    metadata = extract_metadata(inps.unwStack)
    length = metadata["LENGTH"]
    width = metadata["WIDTH"]
    numPairs = metadata["NUMBER_OF_PAIRS"]

    # prepare mintpy working directory
    inps.workDir = os.path.abspath(inps.workDir)
    if not os.path.exists(inps.workDir):
        os.makedirs(inps.workDir)

    inputDir = os.path.join(inps.workDir , "inputs")
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)

    # write geometryGeo
    dsNameDict = {
        "height": (np.float32, (length, width)),
        "incidenceAngle": (np.float32, (length, width)),
        "slantRangeDistance": (np.float32, (length, width)),
    }
    if inps.azAngle is not None:
        dsNameDict["azimuthAngle"] = (np.float32, (length, width))
    if inps.waterMask is not None:
        dsNameDict["waterMask"] = (np.bool_, (length, width))

    geom_file = os.path.join(inputDir, "geometryGeo.h5")
    metadata['FILE_TYPE'] = 'geometry'
    layout_hdf5(geom_file, dsNameDict, metadata)
    write_geometry(geom_file,
                   demFile=inps.dem,
                   incAngleFile=inps.incAngle,
                   azAngleFile=inps.azAngle,
                   waterMaskFile=inps.waterMask)

    # write ifgramStack
    dsNameDict = {
        "bperp": (np.float32, (numPairs,)),
        "coherence": (np.float32, (numPairs, length, width)),
        "connectComponent": (np.int16, (numPairs, length, width)),
        "date": (np.dtype('S8'), (numPairs, 2)),
        "dropIfgram": (np.bool_, (numPairs,)),
        "unwrapPhase": (np.float32, (numPairs, length, width)),
    }

    ifgram_file = os.path.join(inputDir, "ifgramStack.h5")
    metadata['FILE_TYPE'] = 'ifgramStack'
    layout_hdf5(ifgram_file, dsNameDict, metadata)
    write_ifgram_stack(ifgram_file,
                       inps.unwStack,
                       inps.cohStack,
                       inps.connCompStack)
    return ifgram_file, geom_file


####################################################################################
if __name__=="__main__":
    main()
