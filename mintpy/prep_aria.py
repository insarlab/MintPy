#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, Emre Havazli, 2019 #
############################################################


import os
import time
import argparse
import h5py
import numpy as np
import glob
from mintpy.objects import ifgramStack, geometry, sensor
from mintpy.utils import ptime, readfile, writefile, utils as ut
try:
    from osgeo import gdal
except ImportError:
    raise ImportError('Can not import gdal [version>=3.0]!')


####################################################################################
EXAMPLE = """example:
  prep_aria.py -t SanFranSenDT42.txt --update
  prep_aria.py -s stack/ -d DEM/SRTM_3arcsec.dem -i incidenceAngle/*.vrt
  prep_aria.py -s stack/ -d DEM/SRTM_3arcsec.dem -i incidenceAngle/*.vrt  -a azimuthAngle/*.vrt --water-mask mask/watermask.msk

  # before above, one should run ARIA-tools to download / extract / prepare inteferograms stack.
  # reference: https://github.com/aria-tools/ARIA-tools
  ariaDownload.py -b '37.25 38.1 -122.6 -121.75' --track 42
  ariaTSsetup.py -f 'products/*.nc' -b '37.25 38.1 -122.6 -121.75' --mask Download
"""

TEMPLATE = """template options:
  ########## 1. load_data
  ## no   - save   0% disk usage, fast [default]
  ## lzf  - save ~57% disk usage, relative slow
  ## gzip - save ~62% disk usage, very slow [not recommend]
  mintpy.load.processor      = aria  #[isce, aria, snap, gamma, roipac], auto for isce
  mintpy.load.updateMode     = auto  #[yes / no], auto for yes, skip re-loading if HDF5 files are complete
  mintpy.load.compression    = auto  #[gzip / lzf / no], auto for no.
  ##---------interferogram datasets:
  mintpy.load.unwFile        = ../stack/unwrapStack.vrt
  mintpy.load.corFile        = ../stack/cohStack.vrt
  mintpy.load.connCompFile   = ../stack/connCompStack.vrt
  ##---------geometry datasets:
  mintpy.load.demFile        = ../DEM/SRTM_3arcsec.dem
  mintpy.load.incAngleFile   = ../incidenceAngle/*.vrt
  mintpy.load.azAngleFile    = ../azimuthAngle/*.vrt
  mintpy.load.waterMaskFile  = ../mask/watermask.msk
"""


def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ARIA processed products for MintPy.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('-t','--template', dest='template_file', type=str,
                        help='template file with the options')
    parser.add_argument('--update', dest='updateMode', action='store_true',
                        help='Enable the update mode: checking dataset already loaded.')
    parser.add_argument('-o', '--output', type=str, nargs=2, dest='outfile',
                        default=['./inputs/ifgramStack.h5',
                                 './inputs/geometryGeo.h5'],
                        help='output HDF5 file')

    # ifgramStack
    stack = parser.add_argument_group('interferogram stack')
    stack.add_argument('-s','--stack-dir', dest='stackDir', type=str,
                       help='The directory which contains stack VRT files.')
    stack.add_argument('-u','--unwrap-stack-name', dest='unwFile', type=str,
                       default="unwrapStack.vrt",
                       help='Name of the stack VRT file of unwrapped data.\n'+
                            'default: %(default)s')
    stack.add_argument('-c','--coherence-stack-name', dest='corFile', type=str,
                       default="cohStack.vrt",
                       help='Name of the stack VRT file of coherence data.\n'+
                            'default: %(default)s')
    stack.add_argument('-l','--conn-comp-name', dest='connCompFile', type=str,
                       default="connCompStack.vrt",
                       help='Name of the stack VRT file of connected component data.\n' +
                            'default: %(default)s')

    # geometryGeo
    geom = parser.add_argument_group('geometry')
    geom.add_argument('-d','--dem', dest='demFile', type=str,
                      help='Name of the DEM file')
    geom.add_argument('-i','--incidence-angle', dest='incAngleFile', type=str,
                      help='Name of the incidence angle file')
    geom.add_argument('-a','--az-angle','--azimuth-angle', dest='azAngleFile', type=str,
                      help='Name of the azimuth angle file.')
    geom.add_argument('--water-mask', dest='waterMaskFile', type=str,
                      help='Name of the water mask file')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # --template
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # --stack-dir
    elif inps.stackDir is not None:
        inps.stackDir = os.path.abspath(inps.stackDir)
        inps.corFile = os.path.join(inps.stackDir, inps.corFile)
        inps.unwFile = os.path.join(inps.stackDir, inps.unwFile)
        inps.connCompFile = os.path.join(inps.stackDir, inps.connCompFile)

    # check datasets
    # 1. translate wildcard path input with search result
    # 2. raise error if required datasets are missing
    iDict = vars(inps)
    ds_keys = [key for key in list(iDict.keys()) if key.endswith('File')]
    required_ds_keys = ['unwFile', 'corFile', 'demFile', 'incAngleFile']

    for key in ds_keys:
        fname = iDict[key]

        # search for wildcard pattern
        if fname:
            fnames = glob.glob(fname)
        else:
            fnames = []

        # user the first element if more than one exist
        if len(fnames) > 0:
            iDict[key] = fnames[0]

        elif key in required_ds_keys:
            # raise exception if any required DS is missing
            parser.print_usage()
            raise SystemExit('ERROR: no file found for {} in input path: "{}"!'.format(key, iDict[key]))

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template file into inps"""
    if not inps:
        inps = cmd_line_parse()
    iDict = vars(inps)

    print('read options from template file: {}'.format(os.path.basename(template_file)))
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    key_prefix = 'mintpy.load.'
    keys = [i for i in list(iDict.keys()) if key_prefix+i in template.keys()]
    for key in keys:
        value = template[key_prefix+key]
        if value:
            iDict[key] = str(value)

    return inps


def run_or_skip(inps, dsNameDict, out_file):
    flag = 'run'

    # check 1 - update mode status
    if not inps.updateMode:
        return flag

    # check 2 - output file existance
    if ut.run_or_skip(out_file, check_readable=True) == 'run':
        return flag

    # check 3 - output dataset info
    in_size = (inps.length, inps.width)

    if 'unwrapPhase' in dsNameDict.keys():
        # compare date12 and size
        ds = gdal.Open(inps.unwFile, gdal.GA_ReadOnly)
        in_date12_list = [ds.GetRasterBand(i+1).GetMetadata("unwrappedPhase")['Dates']
                          for i in range(inps.num_pair)]
        in_date12_list = ['_'.join(d.split('_')[::-1]) for d in in_date12_list]

        out_obj = ifgramStack(out_file)
        out_obj.open(print_msg=False)
        out_size = (out_obj.length, out_obj.width)
        out_date12_list = out_obj.get_date12_list(dropIfgram=False)

        if out_size == in_size and set(in_date12_list).issubset(set(out_date12_list)):
            print(('All date12   exists in file {} with same size as required,'
                   ' no need to re-load.'.format(os.path.basename(out_file))))
            flag = 'skip'

    elif 'height' in dsNameDict.keys():
        # compare dataset names and size
        in_dsNames = list(dsNameDict.keys())

        out_obj = geometry(out_file)
        out_obj.open(print_msg=False)
        out_size = (out_obj.length, out_obj.width)
        out_dsNames = out_obj.datasetNames

        if out_size == in_size and set(in_dsNames).issubset(set(out_dsNames)):
            print(('All datasets exists in file {} with same size as required,'
                   ' no need to re-load.'.format(os.path.basename(out_file))))
            flag = 'skip'

    return flag


####################################################################################
def extract_metadata(stack):

    meta = {}
    ds = gdal.Open(stack, gdal.GA_ReadOnly)
    if not ds:
        raise RuntimeError('Failed to open file {} with GDAL.'.format(stack))

    # read metadata from unwrapStack.vrt file
    print('extract metadata from {}'.format(stack))
    metaUnw = ds.GetRasterBand(1).GetMetadata("unwrappedPhase")

    # copy over all metadata from unwrapStack
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
    meta["STARTING_RANGE"] = float(meta["startRange"])

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

    # number of independent looks
    sen_dict = sensor.SENSOR_DICT['sen']
    rgfact = sen_dict['IW2']['range_resolution'] / sen_dict['range_pixel_size']
    azfact = sen_dict['IW2']['azimuth_resolution'] / sen_dict['azimuth_pixel_size']
    meta['NCORRLOOKS'] = meta['RLOOKS'] * meta['ALOOKS'] / (rgfact * azfact)

    # geo transformation
    geoTrans = ds.GetGeoTransform()
    lon0 = geoTrans[0]
    lat0 = geoTrans[3]
    lon_step = geoTrans[1]
    lat_step = geoTrans[5]
    lon1 = lon0 + lon_step * meta["WIDTH"]
    lat1 = lat0 + lat_step * meta["LENGTH"]
    meta["X_FIRST"] = '{:.9f}'.format(lon0)
    meta["Y_FIRST"] = '{:.9f}'.format(lat0)
    meta["X_STEP"] = '{:.9f}'.format(lon_step)
    meta["Y_STEP"] = '{:.9f}'.format(lat_step)
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


def write_geometry(outfile, demFile, incAngleFile, azAngleFile=None, waterMaskFile=None):
    print('-'*50)
    print('writing data to HDF5 file {} with a mode ...'.format(outfile))
    h5 = h5py.File(outfile, 'a')

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
        # azimuth angle of the line-of-sight vector:
        # ARIA: vector from target to sensor measured from the east  in counterclockwise direction
        # ISCE: vector from sensor to target measured from the north in counterclockwise direction
        # convert ARIA format to ISCE format, which is used in mintpy
        data -= 90
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
    print('finished writing to HD5 file: {}'.format(outfile))
    return outfile


def write_ifgram_stack(outfile, unwStack, cohStack, connCompStack):

    print('-'*50)
    print('opening {}, {}, {} with gdal ...'.format(os.path.basename(unwStack),
                                                    os.path.basename(cohStack),
                                                    os.path.basename(connCompStack)))
    dsUnw = gdal.Open(unwStack, gdal.GA_ReadOnly)
    dsCoh = gdal.Open(cohStack, gdal.GA_ReadOnly)
    dsComp = gdal.Open(connCompStack, gdal.GA_ReadOnly)

    # extract NoDataValue (from the last */date2_date1.vrt file for example)
    ds = gdal.Open(dsUnw.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValueUnw = ds.GetRasterBand(1).GetNoDataValue()
    print('grab NoDataValue for unwrapPhase:      {:<5} and convert to 0.'.format(noDataValueUnw))

    ds = gdal.Open(dsCoh.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValueCoh = ds.GetRasterBand(1).GetNoDataValue()
    print('grab NoDataValue for coherence:        {:<5} and convert to 0.'.format(noDataValueCoh))

    ds = gdal.Open(dsComp.GetFileList()[-1], gdal.GA_ReadOnly)
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

    print('writing data to HDF5 file {} with a mode ...'.format(outfile))
    h5 = h5py.File(outfile, "a")

    prog_bar = ptime.progressBar(maxValue=nPairs)
    for ii in range(nPairs):
        d12 = d12List[ii]
        bndIdx = d12BandDict[d12]
        prog_bar.update(ii+1, suffix='{}'.format(d12))

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

    prog_bar.close()

    # add MODIFICATION_TIME metadata to each 3D dataset
    for dsName in ['unwrapPhase','coherence','connectComponent']:
        h5[dsName].attrs['MODIFICATION_TIME'] = str(time.time())

    h5.close()
    print('finished writing to HD5 file: {}'.format(outfile))
    dsUnw = None
    dsCoh = None
    dsComp = None
    return outfile


####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.updateMode:
        print('update mode: ON')
    else:
        print('update mode: OFF')

    # extract metadata
    meta = extract_metadata(inps.unwFile)
    inps.length = meta["LENGTH"]
    inps.width = meta["WIDTH"]
    inps.num_pair = meta["NUMBER_OF_PAIRS"]

    # prepare output directory
    out_dir = os.path.dirname(inps.outfile[0])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ########## output file 1 - ifgramStack
    # define dataset structure for ifgramStack
    dsNameDict = {
        "date"             : (np.dtype('S8'), (inps.num_pair, 2)),
        "dropIfgram"       : (np.bool_,       (inps.num_pair,)),
        "bperp"            : (np.float32,     (inps.num_pair,)),
        "unwrapPhase"      : (np.float32,     (inps.num_pair, inps.length, inps.width)),
        "coherence"        : (np.float32,     (inps.num_pair, inps.length, inps.width)),
        "connectComponent" : (np.int16,       (inps.num_pair, inps.length, inps.width)),
    }

    if run_or_skip(inps, dsNameDict, out_file=inps.outfile[0]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'ifgramStack'
        writefile.layout_hdf5(inps.outfile[0], dsNameDict, meta)

        # write data to h5 file in disk
        write_ifgram_stack(inps.outfile[0],
                           inps.unwFile,
                           inps.corFile,
                           inps.connCompFile)

    ########## output file 2 - geometryGeo
    # define dataset structure for geometry
    dsNameDict = {
        "height"             : (np.float32, (inps.length, inps.width)),
        "incidenceAngle"     : (np.float32, (inps.length, inps.width)),
        "slantRangeDistance" : (np.float32, (inps.length, inps.width)),
    }
    if inps.azAngleFile is not None:
        dsNameDict["azimuthAngle"] = (np.float32, (inps.length, inps.width))
    if inps.waterMaskFile is not None:
        dsNameDict["waterMask"]    = (np.bool_,   (inps.length, inps.width))

    if run_or_skip(inps, dsNameDict, out_file=inps.outfile[1]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'geometry'
        writefile.layout_hdf5(inps.outfile[1], dsNameDict, meta)

        # write data to disk
        write_geometry(inps.outfile[1],
                       demFile=inps.demFile,
                       incAngleFile=inps.incAngleFile,
                       azAngleFile=inps.azAngleFile,
                       waterMaskFile=inps.waterMaskFile)

    print('-'*50)
    return inps.outfile


####################################################################################
if __name__=="__main__":
    main()
