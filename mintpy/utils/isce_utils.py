############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, Apr 2020           #
############################################################
# 2020-07: Talib Oliver-Cabrera, add UAVSAR support w/in stripmapStack
# 2020-10: Cunren Liang, add alosStack support
# Group contents:
#     metadata
#     geometry
#     baseline
#     multilook
#     miscellaneous
# Recommend import:
#   from mintpy.utils import isce_utils


import os
import re
import glob
import shelve
import datetime
import numpy as np
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile, utils1 as ut

# suppress matplotlib DEBUG message
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


SPEED_OF_LIGHT = 299792458  # m/s
EARTH_RADIUS = 6378122.65   # m



def get_processor(meta_file):
    """
    Get the name of ISCE processor (imaging mode)
    """
    meta_dir = os.path.dirname(meta_file)
    tops_meta_file = os.path.join(meta_dir, 'IW*.xml')
    stripmap_meta_files = [os.path.join(meta_dir, i) for i in ['data.dat', 'data']]
    alosStack_meta_frame_files = glob.glob(os.path.join(meta_dir, 'f1_*', '*.frame.xml'))

    processor = None
    if len(glob.glob(tops_meta_file)) > 0:
        # topsStack
        processor = 'tops'

    elif any(os.path.isfile(i) for i in stripmap_meta_files):
        # stripmapStack
        processor = 'stripmap'

    elif alosStack_meta_frame_files != []:
        # alosStack
        processor = 'alosStack'

    elif meta_file.endswith('.xml'):
        # stripmapApp
        processor = 'stripmap'

    else:
        raise ValueError('Un-recognized ISCE processor for metadata file: {}'.format(meta_file))
    return processor


#####################################  metadata  #######################################
def load_product(xmlname):
    """Load the product using Product Manager."""
    import isce
    from iscesys.Component.ProductManager import ProductManager as PM
    pm = PM()
    pm.configure()
    obj = pm.loadProduct(xmlname)
    return obj


def extract_isce_metadata(meta_file, geom_dir=None, rsc_file=None, update_mode=True):
    """Extract metadata from ISCE stack products
    Parameters: meta_file : str, path of metadata file, reference/IW1.xml or referenceShelve/data.dat
                geom_dir  : str, path of geometry directory.
                rsc_file  : str, output file name of ROIPAC format rsc file. None for not write to disk.
    Returns:    meta      : dict
                frame     : object, isceobj.Scene.Frame.Frame / isceobj.Scene.Burst.Burst
    """
    # check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=meta_file, check_readable=False) == 'skip':
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
        raise ValueError('un-recognized isce/stripmap metadata file: {}'.format(meta_file))

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
    lat_file = glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.lat'.format(rlooks, alooks)))[0]
    img = isceobj.createImage()
    img.load(lat_file+'.xml')
    width = img.width
    length = img.length
    data = np.memmap(lat_file, dtype='float64', mode='r', shape=(length, width))
    meta['LAT_REF1'] = str(data[0+edge, 0+edge])
    meta['LAT_REF2'] = str(data[0+edge, -1-edge])
    meta['LAT_REF3'] = str(data[-1-edge, 0+edge])
    meta['LAT_REF4'] = str(data[-1-edge, -1-edge])

    lon_file = glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.lon'.format(rlooks, alooks)))[0]
    data = np.memmap(lon_file, dtype='float64', mode='r', shape=(length, width))
    meta['LON_REF1'] = str(data[0+edge, 0+edge])
    meta['LON_REF2'] = str(data[0+edge, -1-edge])
    meta['LON_REF3'] = str(data[-1-edge, 0+edge])
    meta['LON_REF4'] = str(data[-1-edge, -1-edge])

    los_file = glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.los'.format(rlooks, alooks)))[0]
    data = np.memmap(los_file, dtype='float32', mode='r', shape=(length*2, width))[0:length*2:2, :]
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
    stripmapModes = ['UBS', 'UBD', 'HBS', 'HBD', 'HBQ', 'FBS', 'FBD', 'FBQ']
    scansarNominalModes = ['WBS', 'WBD', 'WWS', 'WWD']
    scansarWideModes = ['VBS', 'VBD']
    scansarModes = ['WBS', 'WBD', 'WWS', 'WWD', 'VBS', 'VBD']

    return (spotlightModes, stripmapModes, scansarNominalModes, scansarWideModes, scansarModes)


def extract_image_size_alosStack(geom_dir):
    import isce
    import isceobj

    # grab the number of looks in azimuth / range direction
    lats = glob.glob(os.path.join(geom_dir, '*_*rlks_*alks.lat'))
    rlooks = max([int(os.path.splitext(os.path.basename(x))[0].split('_')[1].strip('rlks')) for x in lats])
    alooks = max([int(os.path.splitext(os.path.basename(x))[0].split('_')[2].strip('alks')) for x in lats])

    # grab the number of rows / coluns
    lat = glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.lat'.format(rlooks, alooks)))[0]
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
    track = load_product(os.path.join(trackDir, '{}.track.xml'.format(dateStr)))

    # read *.frame.xml files
    track.frames = []
    fnames = sorted(glob.glob(os.path.join(trackDir, 'f*_*/{}.frame.xml'.format(dateStr))))
    for fname in fnames:
        track.frames.append(load_product(fname))

    return track



#####################################  geometry  #######################################
def extract_multilook_number(geom_dir, meta=dict(), fext_list=['.rdr','.geo','.rdr.full','.geo.full']):
    for fbase in ['hgt','lat','lon','los']:
        fbase = os.path.join(geom_dir, fbase)
        for fext in fext_list:
            fnames = glob.glob(fbase+fext)
            if len(fnames) > 0:
                break

        if len(fnames) > 0:
            fullXmlFile = '{}.full.xml'.format(fnames[0])
            if os.path.isfile(fullXmlFile):
                fullXmlDict = readfile.read_isce_xml(fullXmlFile)
                xmlDict = readfile.read_attribute(fnames[0])
                meta['ALOOKS'] = int(int(fullXmlDict['LENGTH']) / int(xmlDict['LENGTH']))
                meta['RLOOKS'] = int(int(fullXmlDict['WIDTH']) / int(xmlDict['WIDTH']))
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
        msg += '\n    file basenme: {}'.format(fbase_list)
        msg += '\n    file extension: {}'.format(fext_list)
        print(msg)
        return meta

    print('extract metadata from geometry files: {}'.format([os.path.basename(i) for i in geom_files]))

    # get A/RLOOKS
    meta = extract_multilook_number(geom_dir, meta, fext_list=fext_list)

    # update pixel_size for multilooked data
    meta['rangePixelSize'] *= meta['RLOOKS']
    meta['azimuthPixelSize'] *= meta['ALOOKS']

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
    with open(baseline_file, 'r') as f:
        for line in f:
            l = line.split(":")
            if l[0] == "Bperp (average)":
                bperps.append(float(l[1]))
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
    with open(baseline_file, 'r') as f:
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
    print('read perp baseline time-series from {}'.format(baseline_dir))
    if processor == 'tops':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*/*.txt')))
    elif processor == 'stripmap':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*.txt')))
    elif processor == 'alosStack':
        # all baselines are in baseline_center.txt
        bFiles = glob.glob(os.path.join(baseline_dir, 'baseline_center.txt'))
    else:
        raise ValueError('Un-recognized ISCE stack processor: {}'.format(processor))

    if len(bFiles) == 0:
        print('WARNING: no baseline text file found in dir {}'.format(os.path.abspath(baseline_dir)))
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
        raise ValueError('Un-recognized ISCE stack processor: {}'.format(processor))

    # change reference date
    if ref_date is not None and ref_date != ref_date0:
        ref_date = ptime.yyyymmdd(ref_date)
        print('change reference date to {}'.format(ref_date))
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
    print('Azimuth     pixel size : {:.1f}'.format(az_pixel_size))
    print('Azimuth ground spacing : {:.1f}'.format(az_spacing))
    print('Azimuth ground spacing : {:.1f} after multilooking by {}'.format(az_spacing*az_looks, az_looks))

    print('Range       pixel size : {:.1f}'.format(rg_pixel_size))
    print('Range   ground spacing : {:.1f}'.format(rg_spacing))
    print('Range   ground spacing : {:.1f} after multilooking by {}'.format(rg_spacing*rg_looks, rg_looks))
    return


def resolution2multilook_number(meta_file, resolution):
    """
    Calculate multilook number for InSAR processing given a disired output resolution on the ground

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
    print('Azimuth     pixel size : {:.1f}'.format(az_pixel_size))
    print('Azimuth ground spacing : {:.1f}'.format(az_spacing))
    print('Azimuth ground spacing : {:.1f} after multilooking by {}'.format(az_spacing*az_looks, az_looks))

    print('Range       pixel size : {:.1f}'.format(rg_pixel_size))
    print('Range   ground spacing : {:.1f}'.format(rg_spacing))
    print('Range   ground spacing : {:.1f} after multilooking by {}'.format(rg_spacing*rg_looks, rg_looks))

    return az_looks, rg_looks


def get_full_resolution(meta_file):
    """
    Grab the full resolution in terms of pixel_size and ground spacing
    """
    # check metadata file extension: only ISCE format is supported.
    fext = os.path.splitext(meta_file)[1]
    if fext not in ['.xml', '.dat']:
        raise ValueError('input ISCE metadata file extension "{}" not in [.xml, .dat]'.format(fext))

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

    # grab IPF numver
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
            xml_file = os.path.join(xml_dir, 'IW{}.xml'.format(j+1))
            IPFv = load_product(xml_file).processingSoftwareVersion
            IPF_IW.append('{:.02f}'.format(float(IPFv)))

        prog_bar.update(i+1, suffix='{} IW1/2/3'.format(date_str))
    prog_bar.close()
    return date_list, IPF_IW1, IPF_IW2, IPF_IW3


def safe_list_file2sensor_list(safe_list_file, date_list=None, print_msg=True):
    """Get list of Sentinel-1 sensor names from txt file with SAFE file names.

    Parameters: safe_list_file - str, path of the text file with Sentinel-1 SAFE file path
                                 E.g. SAFE_files.txt
                date_list      - list of str in YYYYMMDD format, reference list of dates
    Returns:    sensor_list    - list of str in S1A or S1B
                date_list      - list of str in YYYYMMDD format
    Example:
        date_list = timeseries('timeseries.h5').get_date_list()
        sensor_list = safe_list_file2sensor_list('../SAFE_files.txt',
                                                 date_list=date_list,
                                                 print_msg=False)[0]
        s1b_dates = [i for i, j in zip(date_list, sensor_list) if j == 'S1B']
        np.savetxt('S1B_date.txt', np.array(s1b_dates).reshape(-1,1), fmt='%s')
    """
    # read txt file
    fc = np.loadtxt(safe_list_file, dtype=str).astype(str).tolist()
    safe_fnames = [os.path.basename(i) for i in fc]

    # get date_list
    date_list_out = [re.findall('_\d{8}T', i)[0][1:-1] for i in safe_fnames]
    date_list_out = sorted(list(set(date_list_out)))

    # get sensor_list
    sensor_list = []
    for d in date_list_out:
        safe_fname = [i for i in safe_fnames if d in i][0]
        sensor = safe_fname.split('_')[0]
        sensor_list.append(sensor)

    # update against date_list_ref
    if date_list is not None:
        # check possible missing dates
        dates_missing = [i for i in date_list if i not in date_list_out]
        if dates_missing:
            raise ValueError('The following dates are missing:\n{}'.format(dates_missing))

        # prune dates not-needed
        flag = np.array([i in date_list for i in date_list_out], dtype=np.bool_)
        if np.sum(flag) > 0:
            sensor_list = np.array(sensor_list)[flag].tolist()
            dates_removed = np.array(date_list_out)[~flag].tolist()
            date_list_out = np.array(date_list_out)[flag].tolist()
            if print_msg:
                print('The following dates are not needed and removed:\n{}'.format(dates_removed))

    return sensor_list, date_list


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
        print('[{}/{}] read {}'.format(i+1, num_file, fname))
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
            raise ValueError('The following dates are missing:\n{}'.format(dates_missing))

        # prune dates not-needed
        flag = np.array([i in date_list for i in date_list_out], dtype=np.bool_)
        if np.sum(flag) > 0:
            sensingMid    = np.array(sensingMid)[flag].tolist()
            sensingStart  = np.array(sensingStart)[flag].tolist()
            sensingStop   = np.array(sensingStop)[flag].tolist()
            dates_removed = np.array(date_list_out)[~flag].tolist()
            print('The following dates are not needed and removed:\n{}'.format(dates_removed))

    return sensingMid, sensingStart, sensingStop

