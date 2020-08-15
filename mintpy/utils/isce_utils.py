############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, Apr 2020           #
############################################################
# Recommend import:
#   from mintpy.utils import isce_utils


import os
import glob
import shelve
import numpy as np
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile, utils1 as ut

# suppress matplotlib DEBUG message
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


SPEED_OF_LIGHT = 299792458  # m/s
EARTH_RADIUS = 6378122.65   # m



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


def get_processor(meta_file):
    """
    Get the name of ISCE processor (imaging mode)
    """
    meta_dir = os.path.dirname(meta_file)
    tops_meta_file = os.path.join(meta_dir, 'IW*.xml')
    stripmap_meta_files = [os.path.join(meta_dir, i) for i in ['data.dat', 'data']]

    processor = None
    if len(glob.glob(tops_meta_file)) > 0:
        # topsStack
        processor = 'tops'

    elif any(os.path.isfile(i) for i in stripmap_meta_files):
        # stripmapStack
        processor = 'stripmap'

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
    Returns:    metadata  : dict
                frame     : object, isceobj.Scene.Frame.Frame / isceobj.Scene.Burst.Burst
    """
    # check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=meta_file, check_readable=False) == 'skip':
        return readfile.read_roipac_rsc(rsc_file), None

    # 1. read/extract metadata from XML / shelve file
    processor = get_processor(meta_file)
    if processor == 'tops':
        print('extract metadata from ISCE/topsStack xml file:', meta_file)
        metadata, frame = extract_tops_metadata(meta_file)

    else:
        print('extract metadata from ISCE/stripmapStack shelve file:', meta_file)
        metadata, frame = extract_stripmap_metadata(meta_file)

    # 2. extract metadata from geometry file
    if geom_dir:
        metadata = extract_geometry_metadata(geom_dir, metadata)

    # 3. common metadata
    metadata['PROCESSOR'] = 'isce'
    metadata['ANTENNA_SIDE'] = '-1'

    # convert all value to string format
    for key, value in metadata.items():
        metadata[key] = str(value)

    # write to .rsc file
    metadata = readfile.standardize_metadata(metadata)
    if rsc_file:
        print('writing ', rsc_file)
        writefile.write_roipac_rsc(metadata, rsc_file)
    return metadata, frame


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

    metadata = {}
    metadata['prf'] = burst.prf
    metadata['startUTC'] = burst.burstStartUTC
    metadata['stopUTC'] = burstEnd.burstStopUTC
    metadata['radarWavelength'] = burst.radarWavelength
    metadata['startingRange'] = burst.startingRange
    metadata['passDirection'] = burst.passDirection
    metadata['polarization'] = burst.polarization
    metadata['trackNumber'] = burst.trackNumber
    metadata['orbitNumber'] = burst.orbitNumber
    metadata['PLATFORM'] = sensor.standardize_sensor_name(obj.spacecraftName)

    time_seconds = (burst.burstStartUTC.hour * 3600.0 +
                    burst.burstStartUTC.minute * 60.0 +
                    burst.burstStartUTC.second)
    metadata['CENTER_LINE_UTC'] = time_seconds

    orbit = burst.orbit
    peg = orbit.interpolateOrbit(burst.sensingMid, method='hermite')

    # Sentinel-1 TOPS pixel spacing
    Vs = np.linalg.norm(peg.getVelocity())   #satellite speed
    metadata['azimuthPixelSize'] = Vs*burst.azimuthTimeInterval
    metadata['rangePixelSize'] = burst.rangePixelSize

    # Sentinel-1 TOPS spatial resolution
    iw_str = 'IW2'
    if os.path.basename(xml_file).startswith('IW'):
        iw_str = os.path.splitext(os.path.basename(xml_file))[0]
    metadata['azimuthResolution'] = sensor.SENSOR_DICT['sen'][iw_str]['azimuth_resolution']
    metadata['rangeResolution'] = sensor.SENSOR_DICT['sen'][iw_str]['range_resolution']

    refElp = Planet(pname='Earth').ellipsoid
    llh = refElp.xyz_to_llh(peg.getPosition())
    refElp.setSCH(llh[0], llh[1], orbit.getENUHeading(burst.sensingMid))
    metadata['earthRadius'] = refElp.pegRadCur
    metadata['altitude'] = llh[2]

    # for Sentinel-1
    metadata['beam_mode'] = 'IW'
    metadata['swathNumber'] = burst.swathNumber
    # 1. multipel subswaths
    xml_files = glob.glob(os.path.join(os.path.dirname(xml_file), 'IW*.xml'))
    if len(xml_files) > 1:
        swath_num = [load_product(fname).bursts[0].swathNumber for fname in xml_files]
        metadata['swathNumber'] = ''.join(str(i) for i in sorted(swath_num))

    # 2. calculate ASF frame number for Sentinel-1
    metadata['firstFrameNumber'] = int(0.2 * (burst.burstStartUTC - obj.ascendingNodeTime).total_seconds())
    metadata['lastFrameNumber'] = int(0.2 * (burstEnd.burstStopUTC - obj.ascendingNodeTime).total_seconds())
    return metadata, burst


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

    metadata = {}
    metadata['prf'] = frame.PRF
    metadata['startUTC'] = frame.sensingStart
    metadata['stopUTC'] = frame.sensingStop
    metadata['radarWavelength'] = frame.radarWavelegth
    metadata['startingRange'] = frame.startingRange
    metadata['polarization'] = str(frame.polarization).replace('/', '')
    if metadata['polarization'].startswith("b'"):
        metadata['polarization'] = metadata['polarization'][2:4]
    metadata['trackNumber'] = frame.trackNumber
    metadata['orbitNumber'] = frame.orbitNumber
    metadata['PLATFORM'] = sensor.standardize_sensor_name(frame.platform.getSpacecraftName())

    time_seconds = (frame.sensingStart.hour * 3600.0 + 
                    frame.sensingStart.minute * 60.0 + 
                    frame.sensingStart.second)
    metadata['CENTER_LINE_UTC'] = time_seconds

    orbit = frame.orbit
    peg = orbit.interpolateOrbit(frame.sensingMid, method='hermite')

    Vs = np.linalg.norm(peg.getVelocity())  #satellite speed
    metadata['azimuthResolution'] = frame.platform.antennaLength / 2.0
    metadata['azimuthPixelSize'] = Vs / frame.PRF

    frame.getInstrument()
    rgBandwidth = frame.instrument.pulseLength * frame.instrument.chirpSlope
    metadata['rangeResolution'] = abs(SPEED_OF_LIGHT / (2.0 * rgBandwidth))
    metadata['rangePixelSize'] = frame.instrument.rangePixelSize

    refElp = Planet(pname='Earth').ellipsoid
    llh = refElp.xyz_to_llh(peg.getPosition())
    refElp.setSCH(llh[0], llh[1], orbit.getENUHeading(frame.sensingMid))
    metadata['earthRadius'] = refElp.pegRadCur
    metadata['altitude'] = llh[2]

    # for StripMap
    metadata['beam_mode'] = 'SM'
    return metadata, frame


#####################################  geometry  #######################################
def extract_multilook_number(geom_dir, metadata=dict(), fext_list=['.rdr','.geo','.rdr.full','.geo.full']):
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
                metadata['ALOOKS'] = int(int(fullXmlDict['LENGTH']) / int(xmlDict['LENGTH']))
                metadata['RLOOKS'] = int(int(fullXmlDict['WIDTH']) / int(xmlDict['WIDTH']))
                break

    # default value
    for key in ['ALOOKS', 'RLOOKS']:
        if key not in metadata:
            metadata[key] = 1

    # NCORRLOOKS for coherence calibration
    rgfact = float(metadata['rangeResolution']) / float(metadata['rangePixelSize'])
    azfact = float(metadata['azimuthResolution']) / float(metadata['azimuthPixelSize'])
    metadata['NCORRLOOKS'] = metadata['RLOOKS'] * metadata['ALOOKS'] / (rgfact * azfact)
    return metadata


def extract_geometry_metadata(geom_dir, metadata=dict(), box=None, fbase_list=['hgt','lat','lon','los'],
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
        return metadata

    print('extract metadata from geometry files: {}'.format([os.path.basename(i) for i in geom_files]))

    # get A/RLOOKS
    metadata = extract_multilook_number(geom_dir, metadata, fext_list=fext_list)

    # update pixel_size for multilooked data
    metadata['rangePixelSize'] *= metadata['RLOOKS']
    metadata['azimuthPixelSize'] *= metadata['ALOOKS']

    # get LAT/LON_REF1/2/3/4 and HEADING into metadata
    for geom_file in geom_files:
        if 'lat' in os.path.basename(geom_file):
            data = readfile.read(geom_file, box=box)[0]
            r0, r1 = get_nonzero_row_number(data)
            metadata['LAT_REF1'] = str(data[r0, 0])
            metadata['LAT_REF2'] = str(data[r0, -1])
            metadata['LAT_REF3'] = str(data[r1, 0])
            metadata['LAT_REF4'] = str(data[r1, -1])

        if 'lon' in os.path.basename(geom_file):
            data = readfile.read(geom_file, box=box)[0]
            r0, r1 = get_nonzero_row_number(data)
            metadata['LON_REF1'] = str(data[r0, 0])
            metadata['LON_REF2'] = str(data[r0, -1])
            metadata['LON_REF3'] = str(data[r1, 0])
            metadata['LON_REF4'] = str(data[r1, -1])

        if 'los' in os.path.basename(geom_file):
            data = readfile.read(geom_file, datasetName='az', box=box)[0]
            # HEADING
            data[data == 0.] = np.nan
            az_angle = np.nanmean(data)
            # convert isce azimuth angle to roipac orbit heading angle
            head_angle = -1 * (270 + az_angle)
            head_angle -= np.round(head_angle / 360.) * 360.
            metadata['HEADING'] = str(head_angle)

            # CENTER_INCIDENCE_ANGLE
            data = readfile.read(geom_file, datasetName='inc', box=box)[0]
            data[data == 0.] = np.nan
            inc_angle = data[int(data.shape[0]/2), int(data.shape[1]/2)]
            metadata['CENTER_INCIDENCE_ANGLE'] = str(inc_angle)
    return metadata


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


def read_baseline_timeseries(baseline_dir, processor='tops', ref_date=None):
    """Read bperp time-series from files in baselines directory
    Parameters: baseline_dir : str, path to the baselines directory
                processor    : str, tops     for Sentinel-1/TOPS
                                    stripmap for StripMap data
    Returns:    bDict : dict, in the following format:
                    {'20141213': [0.0, 0.0],
                     '20141225': [104.6, 110.1],
                     ...
                    }
    """

    print('read perp baseline time-series from {}'.format(baseline_dir))
    # grab all existed baseline files
    if processor == 'tops':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*/*.txt')))
    elif processor == 'stripmap':
        bFiles = sorted(glob.glob(os.path.join(baseline_dir, '*.txt')))
    else:
        raise ValueError('Un-recognized ISCE stack processor: {}'.format(processor))
    if len(bFiles) == 0:
        print('WARNING: no baseline text file found in dir {}'.format(os.path.abspath(baseline_dir)))
        return None

    # ignore files with different date1
    # when re-run with different reference date
    date1s = [os.path.basename(i).split('_')[0] for i in bFiles]
    date1 = ut.most_common(date1s)
    bFiles = [i for i in bFiles if os.path.basename(i).split('_')[0] == date1]

    # read files into dict
    bDict = {}
    for bFile in bFiles:
        dates = os.path.basename(bFile).split('.txt')[0].split('_')
        if processor == 'tops':
            bDict[dates[1]] = read_tops_baseline(bFile)
        else:
            bDict[dates[1]] = read_stripmap_baseline(bFile)
    bDict[dates[0]] = [0, 0]

    # change reference date
    if ref_date is not None and ref_date != dates[0]:
        print('change reference date to {}'.format(ref_date))
        ref_bperp = bDict[ref_date]

        for key in bDict.keys():
            bDict[key][0] -= ref_bperp[0]
            bDict[key][1] -= ref_bperp[1]

    return bDict

