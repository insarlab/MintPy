"""Utilities wrapped around ISCE."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, Apr 2020           #
############################################################
# 2020-07: Talib Oliver-Cabrera, add UAVSAR support
# 2020-10: Cunren Liang, add alosStack support
# 2022-06: Yujie Zheng, add standard processing from isce2
# Group contents:
#   metadata
#   geometry
#   baseline
#   multilook
#   miscellaneous
# Recommend import:
#   from mintpy.utils import isce_utils


import datetime
import glob
import logging
import os
import shelve
import time

import numpy as np
from scipy import ndimage

from mintpy.constants import EARTH_RADIUS, SPEED_OF_LIGHT
from mintpy.objects import sensor
from mintpy.utils import (
    attribute as attr,
    ptime,
    readfile,
    utils1 as ut,
    writefile,
)

# suppress matplotlib DEBUG message
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)



def get_processor(meta_file):
    """
    Get the name of ISCE processor (imaging mode)
    """
    meta_dir = os.path.dirname(meta_file)
    tops_meta_file = os.path.join(meta_dir, 'IW*.xml')
    alos_meta_file = os.path.join(meta_dir, '*.track.xml')
    stripmap_meta_files = [os.path.join(meta_dir, i) for i in ['data.db', 'data.dat', 'data']]

    processor = None
    if len(glob.glob(tops_meta_file)) > 0:
        # topsStack
        processor = 'tops'

    elif len(glob.glob(alos_meta_file)) > 0:
        # alosStack / alos2App
        processor = 'alosStack'

    elif any(os.path.isfile(i) for i in stripmap_meta_files):
        # stripmapStack
        processor = 'stripmap'

    elif meta_file.endswith('.xml'):
        # stripmapApp
        processor = 'stripmap'

    else:
        raise ValueError(f'Un-recognized ISCE processor for metadata file: {meta_file}')
    return processor


#####################################  metadata  #######################################
def load_product(xml_name):
    """Load the product using Product Manager."""
    import isce
    from iscesys.Component.ProductManager import ProductManager
    pm = ProductManager()
    pm.configure()
    return pm.loadProduct(xml_name)


def extract_isce_metadata(meta_file, geom_dir=None, rsc_file=None, update_mode=True):
    """Extract metadata from ISCE stack products
    Parameters: meta_file : str, path of metadata file, reference/IW1.xml or referenceShelve/data.dat
                geom_dir  : str, path of geometry directory.
                rsc_file  : str, output file name of ROIPAC format rsc file. None for not write to disk.
    Returns:    meta      : dict
                frame     : object, isceobj.Scene.Frame.Frame / isceobj.Scene.Burst.Burst
    """
    # check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=meta_file, readable=False) == 'skip':
        return readfile.read_roipac_rsc(rsc_file), None

    # 1. read/extract metadata from XML / shelve file
    processor = get_processor(meta_file)
    if processor == 'tops':
        print('extract metadata from ISCE/topsStack xml file:', meta_file)
        meta, frame = extract_tops_metadata(meta_file)

    elif processor == 'alosStack':
        print('extract metadata from ISCE/alosStack xml file:', meta_file)
        meta, frame = extract_alosStack_metadata(meta_file, geom_dir)

    else:
        print('extract metadata from ISCE/stripmapStack shelve file:', meta_file)
        meta, frame = extract_stripmap_metadata(meta_file)

    # 2. extract metadata from geometry file
    if geom_dir:
        if processor != 'alosStack':
            meta = extract_geometry_metadata(geom_dir, meta)

    # 3. common metadata
    meta['PROCESSOR'] = 'isce'
    if 'ANTENNA_SIDE' not in meta.keys():
        meta['ANTENNA_SIDE'] = '-1'

    # convert all value to string format
    for key, value in meta.items():
        meta[key] = str(value)

    # write to .rsc file
    meta = readfile.standardize_metadata(meta)
    if rsc_file:
        print('writing ', rsc_file)
        writefile.write_roipac_rsc(meta, rsc_file)
    return meta, frame


def extract_tops_metadata(xml_file):
    """Read metadata from xml file for Sentinel-1/TOPS
    Parameters: xml_file : str, path of the .xml file, i.e. reference/IW1.xml
    Returns:    meta     : dict, metadata
                burst    : isceobj.Sensor.TOPS.BurstSLC.BurstSLC object
    """
    import isce
    import isceobj
    from isceobj.Planet.Planet import Planet

    obj = load_product(xml_file)
    burst = obj.bursts[0]
    burstEnd = obj.bursts[-1]

    meta = {}
    meta['prf']             = burst.prf
    meta['startUTC']        = burst.burstStartUTC
    meta['stopUTC']         = burstEnd.burstStopUTC
    meta['radarWavelength'] = burst.radarWavelength
    meta['startingRange']   = burst.startingRange
    meta['passDirection']   = burst.passDirection
    meta['polarization']    = burst.polarization
    meta['trackNumber']     = burst.trackNumber
    meta['orbitNumber']     = burst.orbitNumber

    try:
        meta['PLATFORM'] = sensor.standardize_sensor_name(obj.spacecraftName)
    except:
        if os.path.basename(xml_file).startswith('IW'):
            meta['PLATFORM'] = 'sen'

    time_seconds = (burst.sensingMid.hour * 3600.0 +
                    burst.sensingMid.minute * 60.0 +
                    burst.sensingMid.second)
    meta['CENTER_LINE_UTC'] = time_seconds

    orbit = burst.orbit
    peg = orbit.interpolateOrbit(burst.sensingMid, method='hermite')

    # Sentinel-1 TOPS pixel spacing
    Vs = np.linalg.norm(peg.getVelocity())   #satellite speed
    meta['satelliteSpeed'] = Vs
    meta['azimuthPixelSize'] = Vs * burst.azimuthTimeInterval
    meta['rangePixelSize'] = burst.rangePixelSize

    # Sentinel-1 TOPS spatial resolution
    iw_str = 'IW2'
    if os.path.basename(xml_file).startswith('IW'):
        iw_str = os.path.splitext(os.path.basename(xml_file))[0]
    meta['azimuthResolution'] = sensor.SENSOR_DICT['sen'][iw_str]['azimuth_resolution']
    meta['rangeResolution']   = sensor.SENSOR_DICT['sen'][iw_str]['range_resolution']

    elp = Planet(pname='Earth').ellipsoid
    llh = elp.xyz_to_llh(peg.getPosition())
    elp.setSCH(llh[0], llh[1], orbit.getENUHeading(burst.sensingMid))
    meta['HEADING'] = orbit.getENUHeading(burst.sensingMid)
    meta['earthRadius'] = elp.pegRadCur
    meta['altitude'] = llh[2]

    # for Sentinel-1
    meta['beam_mode'] = 'IW'
    meta['swathNumber'] = burst.swathNumber
    # 1. multipel subswaths
    xml_files = glob.glob(os.path.join(os.path.dirname(xml_file), 'IW*.xml'))
    if len(xml_files) > 1:
        swath_num = [load_product(fname).bursts[0].swathNumber for fname in xml_files]
        meta['swathNumber'] = ''.join(str(i) for i in sorted(swath_num))

    # 2. calculate ASF frame number for Sentinel-1
    meta['firstFrameNumber'] = int(0.2 * (burst.burstStartUTC   - obj.ascendingNodeTime).total_seconds())
    meta['lastFrameNumber']  = int(0.2 * (burstEnd.burstStopUTC - obj.ascendingNodeTime).total_seconds())
    return meta, burst


def extract_stripmap_metadata(meta_file):
    """Read metadata from shelve file for StripMap stack from ISCE
    Parameters: meta_file : str, path of the shelve file, i.e. referenceShelve/data.dat
    Returns:    meta      : dict, metadata
                frame     : isceobj.Scene.Frame.Frame object
    """
    import isce
    import isceobj
    from isceobj.Planet.Planet import Planet

    if os.path.basename(meta_file).startswith('data'):
        # shelve file from stripmapStack
        # referenceShelve/data     for uavsar
        # referenceShelve/data.dat for all the others
        fbase = os.path.splitext(meta_file)[0]
        with shelve.open(fbase, flag='r') as mdb:
            frame = mdb['frame']

    elif meta_file.endswith(".xml"):   #XML file from stripmapApp
        frame = load_product(meta_file)

    else:
        raise ValueError(f'un-recognized isce/stripmap metadata file: {meta_file}')

    meta = {}
    meta['prf']             = frame.PRF
    meta['startUTC']        = frame.sensingStart
    meta['stopUTC']         = frame.sensingStop
    meta['radarWavelength'] = frame.radarWavelegth
    meta['startingRange']   = frame.startingRange
    meta['trackNumber']     = frame.trackNumber
    meta['orbitNumber']     = frame.orbitNumber
    meta['PLATFORM'] = sensor.standardize_sensor_name(frame.platform.getSpacecraftName())
    meta['polarization'] = str(frame.polarization).replace('/', '')
    if meta['polarization'].startswith("b'"):
        meta['polarization'] = meta['polarization'][2:4]

    time_seconds = (frame.sensingMid.hour * 3600.0 +
                    frame.sensingMid.minute * 60.0 +
                    frame.sensingMid.second)
    meta['CENTER_LINE_UTC'] = time_seconds

    orbit = frame.orbit
    peg = orbit.interpolateOrbit(frame.sensingMid, method='hermite')

    Vs = np.linalg.norm(peg.getVelocity())  #satellite speed
    meta['satelliteSpeed'] = Vs
    meta['azimuthResolution'] = frame.platform.antennaLength / 2.0
    meta['azimuthPixelSize'] = Vs / frame.PRF

    frame.getInstrument()
    rgBandwidth = frame.instrument.pulseLength * frame.instrument.chirpSlope
    meta['rangeResolution'] = abs(SPEED_OF_LIGHT / (2.0 * rgBandwidth))
    meta['rangePixelSize'] = frame.instrument.rangePixelSize

    elp = Planet(pname='Earth').ellipsoid
    llh = elp.xyz_to_llh(peg.getPosition())
    elp.setSCH(llh[0], llh[1], orbit.getENUHeading(frame.sensingMid))
    meta['HEADING'] = orbit.getENUHeading(frame.sensingMid)
    meta['earthRadius'] = elp.pegRadCur
    meta['altitude'] = llh[2]

    # for StripMap
    meta['beam_mode'] = 'SM'
    return meta, frame


def extract_alosStack_metadata(meta_file, geom_dir):
    """Read metadata for ISCE/alosStack from the following files:

    pairs/*-*/
        {date1}.track.xml
        f1_{frame1}/{date1}.frame.xml
        f2_{frame2}/{date1}.frame.xml
        ...

    Parameters: meta_file : str, path of the track xml file, i.e. pairs/*-*/150408/track.xml
                geom_dir  : str, path of the geometry directory, i.e. dates_resampled/150408/insar
    Returns:    meta      : dict, metadata
                track     : isceobj.Sensor.MultiMode.Track.Track object
    """

    import isce
    import isceobj
    from isceobj.Planet.Planet import Planet

    # default geom_dir
    if not geom_dir:
        geom_dir_cand = os.path.join(os.path.dirname(meta_file), 'insar')
        if os.path.isdir(geom_dir_cand):
            geom_dir = geom_dir_cand

    track = load_track(os.path.dirname(meta_file), dateStr=os.path.basename(meta_file).strip('.track.xml'))
    rlooks, alooks, width, length = extract_image_size_alosStack(geom_dir)
    spotlightModes, stripmapModes, scansarNominalModes, scansarWideModes, scansarModes = alos2_acquisition_modes()

    meta = {}
    meta['prf']             = track.prf
    meta['startUTC']        = track.sensingStart + datetime.timedelta(seconds=(alooks-1.0)/2.0*track.azimuthLineInterval)
    meta['stopUTC']         = meta['startUTC'] + datetime.timedelta(seconds=(length-1)*alooks*track.azimuthLineInterval)
    meta['radarWavelength'] = track.radarWavelength
    meta['startingRange']   = track.startingRange + (rlooks-1.0)/2.0*track.rangePixelSize
    meta['passDirection']   = track.passDirection.upper()
    meta['polarization']    = track.frames[0].swaths[0].polarization
    #meta['trackNumber']     = track.trackNumber
    #meta['orbitNumber']     = track.orbitNumber

    meta['PLATFORM'] = sensor.standardize_sensor_name('alos2')

    sensingSec = (meta['stopUTC'] - meta['startUTC']).total_seconds()
    sensingMid = meta['startUTC'] + datetime.timedelta(seconds=sensingSec/2.0)
    time_seconds = (sensingMid.hour * 3600.0 +
                    sensingMid.minute * 60.0 +
                    sensingMid.second)
    meta['CENTER_LINE_UTC'] = time_seconds

    peg = track.orbit.interpolateOrbit(sensingMid, method='hermite')
    Vs = np.linalg.norm(peg.getVelocity())
    meta['satelliteSpeed'] = Vs
    meta['azimuthPixelSize'] = Vs * track.azimuthLineInterval
    meta['rangePixelSize'] = track.rangePixelSize

    azBandwidth = track.prf * 0.8
    if track.operationMode in scansarNominalModes:
        azBandwidth /= 5.0
    if track.operationMode in scansarWideModes:
        azBandwidth /= 7.0
    #use a mean burst synchronizatino here
    if track.operationMode in scansarModes:
        azBandwidth *= 0.85

    meta['azimuthResolution'] = Vs * (1.0/azBandwidth)
    meta['rangeResolution']   = 0.5 * SPEED_OF_LIGHT * (1.0/track.frames[0].swaths[0].rangeBandwidth)

    elp = Planet(pname='Earth').ellipsoid
    llh = elp.xyz_to_llh(peg.getPosition())
    elp.setSCH(llh[0], llh[1], track.orbit.getENUHeading(sensingMid))
    meta['HEADING'] = track.orbit.getENUHeading(sensingMid)
    meta['earthRadius'] = elp.pegRadCur
    meta['altitude'] = llh[2]

    meta['beam_mode'] = track.operationMode
    meta['swathNumber'] = ''.join(str(swath.swathNumber) for swath in track.frames[0].swaths)

    meta['firstFrameNumber'] = track.frames[0].frameNumber
    meta['lastFrameNumber']  = track.frames[-1].frameNumber

    meta['ALOOKS'] = alooks
    meta['RLOOKS'] = rlooks

    # NCORRLOOKS for coherence calibration
    rgfact = float(meta['rangeResolution']) / float(meta['rangePixelSize'])
    azfact = float(meta['azimuthResolution']) / float(meta['azimuthPixelSize'])
    meta['NCORRLOOKS'] = meta['RLOOKS'] * meta['ALOOKS'] / (rgfact * azfact)

    # update pixel_size for multilooked data
    meta['rangePixelSize'] *= meta['RLOOKS']
    meta['azimuthPixelSize'] *= meta['ALOOKS']

    # LAT/LON_REF1/2/3/4
    edge = 3
    lat_file = glob.glob(os.path.join(geom_dir, f'*_{rlooks}rlks_{alooks}alks.lat'))[0]
    img = isceobj.createImage()
    img.load(lat_file+'.xml')
    width = img.width
    length = img.length
    data = np.memmap(lat_file, dtype='float64', mode='r', shape=(length, width))
    meta['LAT_REF1'] = str(data[ 0+edge,  0+edge])
    meta['LAT_REF2'] = str(data[ 0+edge, -1-edge])
    meta['LAT_REF3'] = str(data[-1-edge,  0+edge])
    meta['LAT_REF4'] = str(data[-1-edge, -1-edge])

    lon_file = glob.glob(os.path.join(geom_dir, f'*_{rlooks}rlks_{alooks}alks.lon'))[0]
    data = np.memmap(lon_file, dtype='float64', mode='r', shape=(length, width))
    meta['LON_REF1'] = str(data[ 0+edge,  0+edge])
    meta['LON_REF2'] = str(data[ 0+edge, -1-edge])
    meta['LON_REF3'] = str(data[-1-edge,  0+edge])
    meta['LON_REF4'] = str(data[-1-edge, -1-edge])

    # CENTER_INCIDENCE_ANGLE is optional
    los_files = glob.glob(os.path.join(geom_dir, f'*_{rlooks}rlks_{alooks}alks.los'))
    if len(los_files) > 0:
        data = np.memmap(los_files[0], dtype='float32', mode='r', shape=(length*2, width))[0:length*2:2, :]
        inc_angle = data[int(length/2), int(width/2)]
        meta['CENTER_INCIDENCE_ANGLE'] = str(inc_angle)

    pointingDirection = {'right': -1, 'left' :1}
    meta['ANTENNA_SIDE'] = str(pointingDirection[track.pointingDirection])

    return meta, track


def alos2_acquisition_modes():
    '''
    return ALOS-2 acquisition mode
    '''

    spotlightModes = ['SBS']
    # StripMap: Ultrafine [3 m], High sensitive [6 m], Fine [10 m]
    stripmapModes = ['UBS', 'UBD', 'HBS', 'HBD', 'HBQ', 'FBS', 'FBD', 'FBQ']
    scansarNominalModes = ['WBS', 'WBD', 'WWS', 'WWD']
    scansarWideModes = ['VBS', 'VBD']
    scansarModes = scansarNominalModes + scansarWideModes

    return (spotlightModes, stripmapModes, scansarNominalModes, scansarWideModes, scansarModes)


def extract_image_size_alosStack(geom_dir):
    import isce
    import isceobj

    # grab the number of looks in azimuth / range direction
    lats = glob.glob(os.path.join(geom_dir, '*_*rlks_*alks.lat'))
    rlooks = max(int(os.path.splitext(os.path.basename(x))[0].split('_')[1].strip('rlks')) for x in lats)
    alooks = max(int(os.path.splitext(os.path.basename(x))[0].split('_')[2].strip('alks')) for x in lats)

    # grab the number of rows / coluns
    lat = glob.glob(os.path.join(geom_dir, f'*_{rlooks}rlks_{alooks}alks.lat'))[0]
    img = isceobj.createImage()
    img.load(lat+'.xml')
    width = img.width
    length = img.length

    return (rlooks, alooks, width, length)


def load_track(trackDir, dateStr):
    '''Load the track using Product Manager.

    Parameters: trackDir - str, directory of the *.track.xml file
                dateStr  - str, date in YYMMDD format
    Returns:    track    - isceobj.Sensor.MultiMode.Track.Track object
    '''

    # read *.track.xml file
    track = load_product(os.path.join(trackDir, f'{dateStr}.track.xml'))

    # read *.frame.xml files
    track.frames = []
    fnames = sorted(glob.glob(os.path.join(trackDir, f'f*_*/{dateStr}.frame.xml')))
    for fname in fnames:
        track.frames.append(load_product(fname))

    return track



#####################################  geometry  #######################################
def extract_multilook_number(geom_dir, meta=dict(), fext_list=['.rdr','.geo','.rdr.full','.geo.full']):
    for fbase in ['hgt','lat','lon','los','shadowMask']:
        fbase = os.path.join(geom_dir, fbase)

        # get the file name of the geometry file of interest
        for fext in fext_list:
            fnames = glob.glob(fbase+fext)
            if len(fnames) > 0:
                fname = fnames[0]

                # get the file name of the full resolution metadata file
                full_meta_files = [f'{fname}.full.xml', f'{fname}.full.vrt']
                full_meta_files = [x for x in full_meta_files if os.path.isfile(x)]
                if len(full_meta_files) > 0:
                    full_meta_file = full_meta_files[0]

                    # calc A/RLOOKS
                    if full_meta_file.endswith('.xml'):
                        full_dict = readfile.read_isce_xml(full_meta_file)
                    else:
                        full_dict = readfile.read_gdal_vrt(full_meta_file)
                    mli_dict = readfile.read_attribute(fname)
                    meta['ALOOKS'] = int(int(full_dict['LENGTH']) / int(mli_dict['LENGTH']))
                    meta['RLOOKS'] = int(int(full_dict['WIDTH']) / int(mli_dict['WIDTH']))
                    break

    # default value
    for key in ['ALOOKS', 'RLOOKS']:
        if key not in meta:
            meta[key] = 1

    # NCORRLOOKS for coherence calibration
    rgfact = float(meta['rangeResolution']) / float(meta['rangePixelSize'])
    azfact = float(meta['azimuthResolution']) / float(meta['azimuthPixelSize'])
    meta['NCORRLOOKS'] = meta['RLOOKS'] * meta['ALOOKS'] / (rgfact * azfact)
    return meta


def extract_geometry_metadata(geom_dir, meta=dict(), box=None, fbase_list=['hgt','lat','lon','los'],
                              fext_list=['.rdr','.geo','.rdr.full','.geo.full']):
    """Extract / update metadata from geometry files

    extract LAT_REF1/2/3/4 from lat*
    extract LON_REF1/2/3/4 from lon*
    extract HEADING from los (azimuth angle)

    extract A/RLOOKS by comparing hgt.xml and hgt.full.xml file
    update azimuthPixelSize / rangePixelSize based on A/RLOOKS

    extract LENGTH/WIDTH from the first geom file
    update corresponding metadata if box is not None
    """

    def get_nonzero_row_number(data, buffer=2):
        """Find the first and last row number of rows without zero value
        for multiple swaths data
        """
        if np.all(data):
            r0, r1 = 0 + buffer, -1 - buffer
        else:
            row_flag = np.sum(data != 0., axis=1) == data.shape[1]
            row_idx = np.where(row_flag)[0]
            r0, r1 = row_idx[0] + buffer, row_idx[-1] - buffer
        return r0, r1

    # grab existing files
    geom_dir = os.path.abspath(geom_dir)
    for fext in fext_list:
        geom_files = [os.path.join(geom_dir, fbase+fext) for fbase in fbase_list]
        geom_files = [i for i in geom_files if os.path.isfile(i)]
        if len(geom_files) > 0:
            break

    # printout message
    if len(geom_files) == 0:
        msg = 'WARNING: No geometry files found with the following pattern!'
        msg += f'\n    file basenme: {fbase_list}'
        msg += f'\n    file extension: {fext_list}'
        print(msg)
        return meta

    print(f'extract metadata from geometry files: {[os.path.basename(i) for i in geom_files]}')

    # get A/RLOOKS
    meta = extract_multilook_number(geom_dir, meta, fext_list=fext_list)

    # update pixel_size for multilooked data
    meta['rangePixelSize'] *= meta['RLOOKS']
    meta['azimuthPixelSize'] *= meta['ALOOKS']

    # get LENGTH/WIDTH
    atr = readfile.read_attribute(geom_files[0])
    meta['LENGTH'] = atr['LENGTH']
    meta['WIDTH'] = atr['WIDTH']

    # update due to subset
    if box:
        meta = attr.update_attribute4subset(meta, box)

    # get LAT/LON_REF1/2/3/4 into metadata
    for geom_file in geom_files:
        if 'lat' in os.path.basename(geom_file):
            data = readfile.read(geom_file, box=box)[0]
            r0, r1 = get_nonzero_row_number(data)
            meta['LAT_REF1'] = str(data[r0, 0])
            meta['LAT_REF2'] = str(data[r0, -1])
            meta['LAT_REF3'] = str(data[r1, 0])
            meta['LAT_REF4'] = str(data[r1, -1])

        if 'lon' in os.path.basename(geom_file):
            data = readfile.read(geom_file, box=box)[0]
            r0, r1 = get_nonzero_row_number(data)
            meta['LON_REF1'] = str(data[r0, 0])
            meta['LON_REF2'] = str(data[r0, -1])
            meta['LON_REF3'] = str(data[r1, 0])
            meta['LON_REF4'] = str(data[r1, -1])

        if 'los' in os.path.basename(geom_file):
            # CENTER_INCIDENCE_ANGLE
            data = readfile.read(geom_file, datasetName='inc', box=box)[0]
            data[data == 0.] = np.nan
            inc_angle = data[int(data.shape[0]/2), int(data.shape[1]/2)]
            meta['CENTER_INCIDENCE_ANGLE'] = str(inc_angle)
    return meta



#####################################  baseline  #######################################
def read_tops_baseline(baseline_file):
    """Read baseline file generated by ISCE/topsStack processor.
    Example:
    baselines/20141213_20160418/20141213_20160418.txt:
        swath: IW1
        Bperp (average): 62.62863491739495
        Bpar (average): -29.435602419751426
        swath: IW2
        Bperp (average): 60.562020649374034
        Bpar (average): -34.56105358031081
    """
    bperps = []
    with open(baseline_file) as f:
        for line in f:
            c = line.split(":")
            if c[0] == "Bperp (average)":
                bperps.append(float(c[1]))
    bperp_top = np.mean(bperps)
    bperp_bottom = np.mean(bperps)
    return [bperp_top, bperp_bottom]


def read_stripmap_baseline(baseline_file):
    """Read baseline file generated by ISCE/stripmapStack processor.

    Example:
    baselines/20200111_20200125.txt
        PERP_BASELINE_BOTTOM 173.97914535263297
        PERP_BASELINE_TOP 174.05612879066618
    """
    fDict = readfile.read_template(baseline_file, delimiter=' ')
    bperp_top = float(fDict['PERP_BASELINE_TOP'])
    bperp_bottom = float(fDict['PERP_BASELINE_BOTTOM'])
    return [bperp_top, bperp_bottom]


def read_alosStack_baseline(baseline_file):
    '''read baseline file generated by alosStack
    '''
    bDict = {}
    with open(baseline_file) as f:
        lines = [line for line in f if line.strip() != '']
        for x in lines[2:]:
            blist = x.split()
            #to fit into the format of other processors, all alos satellites are after 2000
            blist[0] = '20' + blist[0]
            blist[1] = '20' + blist[1]
            bDict[blist[1]] = [float(blist[3]), float(blist[3])]
        bDict[blist[0]] = [0, 0]

    return bDict, blist[0]


def read_baseline_timeseries(baseline_dir, processor='tops', ref_date=None):
    """Read bperp time-series from files in baselines directory
    Parameters: baseline_dir : str, path to the baselines directory
                processor    : str, tops     for Sentinel-1/TOPS
                                    stripmap for StripMap data
                ref_date     : str, reference date in (YY)YYMMDD
    Returns:    bDict : dict, in the following format:
                    {'20141213': [0.0, 0.0],
                     '20141225': [104.6, 110.1],
                     ...
                    }
    """

    # grab all existed baseline files
    print(f'read perp baseline time-series from {baseline_dir}')
    if processor == 'tops':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*/*.txt')))
    elif processor == 'stripmap':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*.txt')))
    elif processor == 'alosStack':
        # all baselines are in baseline_center.txt
        bFiles = glob.glob(os.path.join(baseline_dir, 'baseline_center.txt'))
    else:
        raise ValueError(f'Un-recognized ISCE stack processor: {processor}')

    if len(bFiles) == 0:
        print(f'WARNING: no baseline text file found in dir {os.path.abspath(baseline_dir)}')
        return None

    if processor in ['tops', 'stripmap']:
        # ignore files with different date1
        # when re-run with different reference date
        date1s = [os.path.basename(i).split('_')[0] for i in bFiles]
        date1c = ut.most_common(date1s)
        bFiles = [i for i in bFiles if os.path.basename(i).split('_')[0] == date1c]

        # ignore empty files
        bFiles = [i for i in bFiles if os.path.getsize(i) > 0]

        # read files into dict
        bDict = {}
        for bFile in bFiles:
            dates = os.path.basename(bFile).split('.txt')[0].split('_')
            if processor == 'tops':
                bDict[dates[1]] = read_tops_baseline(bFile)
            else:
                bDict[dates[1]] = read_stripmap_baseline(bFile)
        bDict[dates[0]] = [0, 0]
        ref_date0 = dates[0]

    elif processor == 'alosStack':
        bDict, ref_date0 = read_alosStack_baseline(bFiles[0])

    else:
        raise ValueError(f'Un-recognized ISCE stack processor: {processor}')

    # change reference date
    if ref_date is not None and ref_date != ref_date0:
        ref_date = ptime.yyyymmdd(ref_date)
        print(f'change reference date to {ref_date}')
        ref_bperp = bDict[ref_date]

        for key in bDict.keys():
            bDict[key][0] -= ref_bperp[0]
            bDict[key][1] -= ref_bperp[1]

    return bDict



#####################################  multilook  #######################################

def multilook_number2resolution(meta_file, az_looks, rg_looks):
    # get full resolution info
    az_pixel_size, az_spacing, rg_pixel_size, rg_spacing = get_full_resolution(meta_file)

    # print out message
    print(f'Azimuth     pixel size : {az_pixel_size:.1f}')
    print(f'Azimuth ground spacing : {az_spacing:.1f}')
    print(f'Azimuth ground spacing : {az_spacing*az_looks:.1f} after multilooking by {az_looks}')

    print(f'Range       pixel size : {rg_pixel_size:.1f}')
    print(f'Range   ground spacing : {rg_spacing:.1f}')
    print(f'Range   ground spacing : {rg_spacing*rg_looks:.1f} after multilooking by {rg_looks}')
    return


def resolution2multilook_number(meta_file, resolution):
    """
    Calculate multilook number for InSAR processing given a desired output resolution on the ground

    Parameters: meta_file   : str, path of ISCE metadata file, i.e. IW1.xml, data.dat
                resolution  : float, target output resolution on the ground in meters
    Returns:    az/rg_looks : int, number of looks in azimuth / range direction
    """
    # get full resolution info
    az_pixel_size, az_spacing, rg_pixel_size, rg_spacing = get_full_resolution(meta_file)

    # calculate number of looks
    # 1. adjust the final resolution in one direction closest to the input resolution
    az_looks = resolution / az_spacing
    rg_looks = resolution / rg_spacing
    az_round_frac = abs(az_looks - np.rint(az_looks))
    rg_round_frac = abs(rg_looks - np.rint(rg_looks))
    if az_round_frac < rg_round_frac:
        resolution = np.rint(az_looks) * az_spacing
    else:
        resolution = np.rint(rg_looks) * rg_spacing

    # 2. calculate the multilook number based on the adjusted resolution
    az_looks = np.rint(resolution / az_spacing).astype(int)
    rg_looks = np.rint(resolution / rg_spacing).astype(int)

    # print out message
    print(f'Azimuth     pixel size : {az_pixel_size:.1f}')
    print(f'Azimuth ground spacing : {az_spacing:.1f}')
    print(f'Azimuth ground spacing : {az_spacing*az_looks:.1f} after multilooking by {az_looks}')

    print(f'Range       pixel size : {rg_pixel_size:.1f}')
    print(f'Range   ground spacing : {rg_spacing:.1f}')
    print(f'Range   ground spacing : {rg_spacing*rg_looks:.1f} after multilooking by {rg_looks}')

    return az_looks, rg_looks


def get_full_resolution(meta_file):
    """
    Grab the full resolution in terms of pixel_size and ground spacing
    """
    # check metadata file extension: only ISCE format is supported.
    fext = os.path.splitext(meta_file)[1]
    if fext not in ['.xml', '.dat']:
        raise ValueError(f'input ISCE metadata file extension "{fext}" not in [.xml, .dat]')

    # get middle sub-swath xml file for Sentinel-1 data
    if meta_file.endswith('.xml'):
        meta_files = glob.glob(meta_file)
        mid_idx = int(len(meta_files) / 2)
        meta_file = meta_files[mid_idx]

    # extract metadata
    meta, frame = extract_isce_metadata(meta_file, update_mode=False)
    meta['WIDTH'] = frame.numberOfSamples

    # calculate the full azimuth/range ground resolution
    az_pixel_size = float(meta['AZIMUTH_PIXEL_SIZE'])  #azimuth pixel size on the orbit
    rg_pixel_size = float(meta['RANGE_PIXEL_SIZE'])    #range   pixel size in LOS direction
    az_pixel_size /= int(meta.get('ALOOKS', 1))
    rg_pixel_size /= int(meta.get('RLOOKS', 1))

    height = float(meta['HEIGHT'])
    inc_angle = ut.incidence_angle(meta, dimension=0)

    az_spacing = az_pixel_size * EARTH_RADIUS / (EARTH_RADIUS + height)  #azimuth pixel size on the ground
    rg_spacing = rg_pixel_size / np.sin(inc_angle / 180. * np.pi)        #range   pixel size on the ground

    return az_pixel_size, az_spacing, rg_pixel_size, rg_spacing



#####################################  miscellaneous  #######################################

def get_IPF(proj_dir, ts_file):
    """Grab the IPF version number of each sub-swatch for Sentinel-1 time-series

    Parameters: proj_dir    - str, path of the project directory
                              E.g.: ~/data/AtacamaSenDT149
                ts_file     - str, path of HDF5 file for time-series
    Returns:    date_list   - list of str, dates in YYYYMMDD format
                IFP_IW1/2/3 - list of str, IFP version number
    """
    from mintpy.objects import timeseries

    s_dir = os.path.join(proj_dir, 'secondarys')
    m_dir = os.path.join(proj_dir, 'reference')

    # date list
    date_list = timeseries(ts_file).get_date_list()
    num_date = len(date_list)
    # reference date
    m_date = [i for i in date_list if not os.path.isdir(os.path.join(s_dir, i))][0]

    # grab IPF number
    IPF_IW1, IPF_IW2, IPF_IW3 = [], [], []
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        date_str = date_list[i]

        # get xml_dir
        if date_str == m_date:
            xml_dir = m_dir
        else:
            xml_dir = os.path.join(s_dir, date_str)

        # grab IPF version number
        for j, IPF_IW in enumerate([IPF_IW1, IPF_IW2, IPF_IW3]):
            xml_file = os.path.join(xml_dir, f'IW{j+1}.xml')
            IPFv = load_product(xml_file).processingSoftwareVersion
            IPF_IW.append(f'{float(IPFv):.02f}')

        prog_bar.update(i+1, suffix=f'{date_str} IW1/2/3')
    prog_bar.close()
    return date_list, IPF_IW1, IPF_IW2, IPF_IW3


def get_sensing_datetime_list(proj_dir, date_list=None):
    """Get the sensing datetime objects from ISCE stack results.
    It assumes the default directory structure from topsStack, as below:
    /proj_dir
        /reference/IW*.xml
        /secondarys
            /20150521/IW*.xml
            /20150614/IW*.xml
            ...
            /20210113/IW*.xml

    Parameters: proj_dir     - str, path to the root directory of stack processing
    Returns:    sensingMid   - list of datetime.datetime.obj
                sensingStart - list of datetime.datetime.obj
                sensingStop  - list of datetime.datetime.obj
    """
    # determine xml file basename
    ref_fname = glob.glob(os.path.join(proj_dir, 'reference', 'IW*.xml'))[0]
    fbase = os.path.basename(ref_fname)

    # get xml files for all acquisitions
    sec_fnames = sorted(glob.glob(os.path.join(proj_dir, 'secondarys', '*', fbase)))
    fnames = [ref_fname] + sec_fnames
    num_file = len(fnames)

    # loop to read file one by one
    sensingStart = []
    sensingStop = []
    for i, fname in enumerate(fnames):
        print(f'[{i+1}/{num_file}] read {fname}')
        obj = load_product(fname)
        sensingStart.append(obj.bursts[0].sensingStart)
        sensingStop.append(obj.bursts[-1].sensingStop)

    sensingStart = sorted(sensingStart)
    sensingStop  = sorted(sensingStop)

    # sensingStart/Stop --> sensingMid
    sensingMid = [i + (j - i)/2 for i, j in zip(sensingStart, sensingStop)]

    # round to the nearest second
    print('round sensingStart/Stop/Mid to the nearest second.')
    sensingStart = [ptime.round_seconds(i) for i in sensingStart]
    sensingStop  = [ptime.round_seconds(i) for i in sensingStop]
    sensingMid   = [ptime.round_seconds(i) for i in sensingMid]

    if date_list is not None:
        date_str_format = ptime.get_date_str_format(date_list[0])
        date_list_out = [i.strftime(date_str_format) for i in sensingMid]

        # check possible missing dates
        dates_missing = [i for i in date_list if i not in date_list_out]
        if dates_missing:
            raise ValueError(f'The following dates are missing:\n{dates_missing}')

        # prune dates not-needed
        flag = np.array([i in date_list for i in date_list_out], dtype=np.bool_)
        if np.sum(flag) > 0:
            sensingMid    = np.array(sensingMid)[flag].tolist()
            sensingStart  = np.array(sensingStart)[flag].tolist()
            sensingStop   = np.array(sensingStop)[flag].tolist()
            dates_removed = np.array(date_list_out)[~flag].tolist()
            print(f'The following dates are not needed and removed:\n{dates_removed}')

    return sensingMid, sensingStart, sensingStop



############################## Standard Processing ###########################################

def gaussian_kernel(sx, sy, sig_x, sig_y):
    '''Generate a Gaussian kernel (with all elements sum to 1).

    Parameters: sx/y    - int, dimensions of kernel
                sig_x/y - float, standard deviation of the Gaussian distribution
    '''
    # ensure sx/y are odd number
    sx += 1 if np.mod(sx, 2) == 0 else 0
    sy += 1 if np.mod(sy, 2) == 0 else 0

    x, y = np.meshgrid(np.arange(sx), np.arange(sy))
    x += 1
    y += 1

    xc = (sx + 1) / 2
    yc = (sy + 1) / 2
    fx = ((x-xc)**2.) / (2.*sig_x**2.)
    fy = ((y-yc)**2.) / (2.*sig_y**2.)

    k = np.exp(-1.0 * (fx+fy))
    a = 1./np.sum(k)
    k = a * k

    return k


def convolve(data, kernel):
    '''Convolve / filter the complex data based on the given kernel.

    Parameters: data   - 2D np.ndarray in complex
                kernel - 2D np.ndarray in float, convolution kernel
    '''
    real = ndimage.convolve(data.real, kernel, mode='constant', cval=0.0)
    imag = ndimage.convolve(data.imag, kernel, mode='constant', cval=0.0)
    return real + 1J * imag


def filter_goldstein(int_file, filt_file, filt_strength=0.2):
    """Filter wrapped interferogram with the power-spectral filter via isce2.

    Modified from ISCE-2/topsStack/FilterAndCoherence.py
    Reference: Goldstein, R. M., & Werner, C. L. (1998). Radar interferogram
        filtering for geophysical applications. Geophysical Research Letters,
        25(21), 4035-4038. doi:10.1029/1998GL900033

    Parameters: int_file      - str, path of wrapped interferogram
                filt_file     - str, path of filtered wrapped interferogram
                filt_strength - float, filtering strength between 0 and 1
    Returns:    filt_file     - str, path of filtered wrapped interferogram
    """

    import isce
    import isceobj
    from mroipac.filter.Filter import Filter
    print(f"Applying power-spectral filter (strength={filt_strength})...")

    # initialize the flattened interferogram
    int_img = isceobj.createIntImage()
    int_img.load(int_file + '.xml')
    int_img.setAccessMode('read')
    int_img.createImage()

    # create the filtered interferogram
    filt_img = isceobj.createIntImage()
    filt_img.setFilename(filt_file)
    filt_img.setWidth(int_img.getWidth())
    filt_img.setAccessMode('write')
    filt_img.createImage()

    # filter
    filt_obj = Filter()
    filt_obj.wireInputPort(name='interferogram', object=int_img)
    filt_obj.wireOutputPort(name='filtered interferogram', object=filt_img)
    filt_obj.goldsteinWerner(alpha=filt_strength)

    # close
    int_img.finalizeImage()
    filt_img.finalizeImage()

    return filt_file


def estimate_coherence(intfile, corfile):
    '''Estimate the spatial coherence (phase sigma) of the wrapped interferogram.

    Parameters: intfile - str, path to the *.int file
                corfile - str, path to the output correlation file
    '''
    import isce
    import isceobj
    from mroipac.icu.Icu import Icu

    # create filt interferogram file object
    filtImage = isceobj.createIntImage()
    filtImage.load(intfile + '.xml')
    filtImage.setAccessMode('read')
    filtImage.createImage()

    # create phase sigma correlation file object
    phsigImage = isceobj.createImage()
    phsigImage.dataType='FLOAT'
    phsigImage.bands = 1
    phsigImage.setWidth(filtImage.getWidth())
    phsigImage.setFilename(corfile)
    phsigImage.setAccessMode('write')
    phsigImage.createImage()

    # setup Icu() object
    icuObj = Icu(name='sentinel_filter_icu')
    icuObj.configure()
    icuObj.unwrappingFlag = False
    icuObj.useAmplitudeFlag = False
    #icuObj.correlationType = 'NOSLOPE'

    # run
    icuObj.icu(intImage=filtImage, phsigImage=phsigImage)
    phsigImage.renderHdr()

    # close
    filtImage.finalizeImage()
    phsigImage.finalizeImage()

    return


def unwrap_snaphu(int_file, cor_file, unw_file, max_defo=2.0, max_comp=32,
                  init_only=True, init_method='MCF', cost_mode='SMOOTH'):
    '''Unwrap interferograms using SNAPHU via isce2.

    Modified from ISCE-2/topsStack/unwrap.py
    Notes from Piyush:
        SNAPHU is an iterative solver, starting from the initial solution. It can get
            stuck in an infinite loop.
        The initial solution is created using MCF or MST method. The MST initial solution
            typically require lots of iterations and may not be a good starting point.
        DEFO cost mode requires geometry info for the program to interpret the coherence
            correctly and setup costs based on that. DEFO always sounds more theoretical
            to me, but I haven not fully explored it. TOPO cost mode requires spatial baseline.
            SMOOTH cost mode is purely data driven.
        Amplitude of zero is a mask in all cost modes. For TOPO mode, amplitude is used to find
            layover; for SMOOTH mode, only non-zero amplitude matters.

    Default configurations in ISCE-2/topsStack:
        init_only = True
        init_method = 'MCF'
        cost_mode = 'SMOOTH'
    Default configurations in FRInGE:
        init_only = False
        init_method = 'MST'
        cost_mode = 'DEFO'

    Parameters: int_file    - str, path to the wrapped interferogram file
                cor_file    - str, path to the correlation file: phase sigma or complex correlation
                unw_file    - str, path to the output unwrapped interferogram file
                max_defo    - float, maximum number of cycles for the deformation phase
                max_comp    - int, maximum number of connected components
                init_only   - bool, initlize-only mode
                init_method - str, algo used for initialization: MCF, MST
                cost_mode   - str, statistical-cost mode: TOPO, DEFO, SMOOTH, NOSTATCOSTS
    Returns:    unw_file    - str, path to the output unwrapped interferogram file
    '''
    import isce
    from contrib.Snaphu.Snaphu import Snaphu

    start_time = time.time()

    # configurations - atr
    atr = readfile.read_attribute(int_file)
    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    altitude = float(atr['HEIGHT'])
    earth_radius = float(atr['EARTH_RADIUS'])
    wavelength = float(atr['WAVELENGTH'])
    rg_looks = int(atr['RLOOKS'])
    az_looks = int(atr['ALOOKS'])
    corr_looks = float(atr.get('NCORRLOOKS', rg_looks * az_looks / 1.94))

    ## setup SNAPHU
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu.conf.full
    # https://github.com/isce-framework/isce2/blob/main/contrib/Snaphu/Snaphu.py
    print('phase unwrapping with SNAPHU ...')
    print(f'SNAPHU cost mode: {cost_mode}')
    print(f'SNAPHU init only: {init_only}')
    print(f'SNAPHU init method: {init_method}')
    print(f'SNAPHU max number of connected components: {max_comp}')

    snp = Snaphu()

    # file IO
    snp.setInput(int_file)
    snp.setOutput(unw_file)
    snp.setCorrfile(cor_file)
    snp.setWidth(width)

    atr_cor = readfile.read_attribute(cor_file)
    if int(atr_cor.get('BANDS', 1)) == 1:
        snp.setCorFileFormat('FLOAT_DATA')

    # runtime options
    snp.setCostMode(cost_mode)
    snp.setInitOnly(init_only)
    snp.setInitMethod(init_method)

    # geometry parameters
    # baseline info is not used in deformation mode, but is very important in topography mode
    snp.setAltitude(altitude)
    snp.setEarthRadius(earth_radius)
    snp.setWavelength(wavelength)
    snp.setRangeLooks(rg_looks)
    snp.setAzimuthLooks(az_looks)
    snp.setCorrLooks(corr_looks)

    # deformation mode parameters
    snp.setDefoMaxCycles(max_defo)

    # connected component control
    # grow connectedc components if init_only is True
    # https://github.com/isce-framework/isce2/blob/main/contrib/Snaphu/Snaphu.py#L413
    snp.setMaxComponents(max_comp)

    ## run SNAPHU
    snp.prepare()
    snp.unwrap()
    print('finished SNAPHU running')

    # mask out wired values from SNAPHU
    # based on https://github.com/isce-framework/isce2/pull/326
    flag = np.fromfile(int_file, dtype=np.complex64).reshape(length, width)
    data = np.memmap(unw_file, dtype='float32', mode='r+', shape=(length*2, width))
    data[0:length*2:2, :][np.nonzero(flag == 0)] = 0
    data[1:length*2:2, :][np.nonzero(flag == 0)] = 0

    ## render metadata
    print(f'write metadata file: {unw_file}.xml')
    atr['FILE_TYPE'] = '.unw'
    atr['DATA_TYPE'] = 'float32'
    atr['INTERLEAVE'] = 'BIL'
    atr['BANDS'] = '2'
    writefile.write_isce_xml(atr, unw_file)

    if snp.dumpConnectedComponents:
        print(f'write metadata file: {unw_file}.conncomp.xml')
        atr['FILE_TYPE'] = '.conncomp'
        atr['DATA_TYPE'] = 'uint8'
        atr['INTERLEAVE'] = 'BIP'
        atr['BANDS'] = '1'
        writefile.write_isce_xml(atr, f'{unw_file}.conncomp')

    # time usage
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return unw_file


def unwrap_icu(int_file, unw_file):
    """Unwrap interferograms using ICU via isce2.

    Modified from ISCE-2/topsStack/unwrap.py.
    Parameters: int_file - str, path of   wrapped interferogram
                unw_file - str, path of unwrapped interferogram
    Returns:    unw_file - str, path of unwrapped interferogram
    """
    import isce
    import isceobj
    from mroipac.icu.Icu import Icu

    start_time = time.time()

    # get width
    img = isceobj.Image.createImage()
    img.load(int_file + '.xml')
    width = img.getWidth()

    # create image object for .int file
    int_img = isceobj.Image.createIntImage()
    int_img.initImage(int_file, 'read', width)
    int_img.createImage()

    # create image object for .unw file
    unw_img = isceobj.Image.createImage()
    unw_img.setFilename(unw_file)
    unw_img.setWidth(width)
    unw_img.imageType = 'unw'
    unw_img.bands = 2
    unw_img.scheme = 'BIL'
    unw_img.dataType = 'FLOAT'
    unw_img.setAccessMode('write')
    unw_img.createImage()

    # run ICU
    icu_obj = Icu()
    icu_obj.filteringFlag = False
    icu_obj.useAmplitudeFlag = False
    icu_obj.singlePatch = True
    icu_obj.initCorrThresdhold = 0.1
    icu_obj.icu(intImage=int_img, unwImage=unw_img)

    int_img.finalizeImage()
    unw_img.finalizeImage()
    unw_img.renderHdr()
    unw_img.renderVRT()

    # time usage
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return unw_file
