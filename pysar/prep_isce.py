#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
# Modified from prep4timeseries.py in ISCE-2.2.0/contrib/stack/topsStack
#


import argparse
import os
import glob
#import gdal
#from gdalconst import GA_ReadOnly
import numpy as np
#import isce
#from isceobj.Planet.Planet import Planet
from pysar.utils import readfile, writefile, utils as ut


EXAMPLE = """example:
  prep_isce.py -i ./merged/interferograms -x ./master/IW1.xml -b ./baselines -g ./merged/geom_master
  prep_isce.py -i ./merged/interferograms -x ./master/IW1.xml -b ./baselines
  prep_isce.py -i ./merged/interferograms -x ./master/IW1.xml
  prep_isce.py -g ./merged/geom_master
"""

def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ISCE files for PySAR time series analysis.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('-i', '--ifg-dir', dest='ifgramDir', type=str, default=None,
                        help='The directory which contains all pairs\n'+
                             'e.g.: $PROJECT_DIR/merged/interferograms')
    parser.add_argument('-f', '--file-pattern', nargs = '+', dest='ifgramFiles', type=str,
                        default=['filt_*.unw','filt_*.cor'],
                        help='A list of files that will be used in pysar e.g.: filt_fine.unw filt_fine.cor')
    parser.add_argument('-x', '--xml-file', dest='xmlFile', type=str, default=None,
                        help='An xml file to extract common metada for the stack: e.g.: master/IW3.xml')
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default=None,
                        help=' directory with baselines ')
    parser.add_argument('-g', '--geometry-dir', dest='geometryDir', type=str, default=None,
                        help=' directory with geometry files ')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if not inps.ifgramDir and not inps.geometryDir:
        parser.print_usage()
        raise SystemExit('error: at least one of the following arguments are required: -i, -g')
    return inps


#########################################################################
def extract_isce_metadata(xmlFile, rscFile=None):
    import isce
    from isceobj.Planet.Planet import Planet

    print('extract metadata from xml file:',xmlFile)
    metadata = {}

    master = load_product(xmlFile)
    metadata['spacecraftName'] = master.spacecraftName

    burst = master.bursts[0]
    burstEnd = master.bursts[-1]
    metadata['radarWavelength'] = burst.radarWavelength
    metadata['rangePixelSize'] = burst.rangePixelSize
    metadata['prf'] = burst.prf
    metadata['startUTC'] = burst.burstStartUTC
    metadata['stopUTC'] = burstEnd.burstStopUTC
    metadata['startingRange'] = burst.startingRange
    metadata['passDirection'] = burst.passDirection
    metadata['polarization'] = burst.polarization
    metadata['trackNumber'] = burst.trackNumber
    metadata['orbitNumber'] = burst.orbitNumber

    metadata['swathNumber'] = burst.swathNumber
    # swathNumber for multipel subswaths
    xml_files = glob.glob(os.path.join(os.path.dirname(xmlFile), 'IW*.xml'))
    if len(xml_files) > 1:
        swath_num = []
        for xml_file in xml_files:
            swath_num.append(load_product(xml_file).bursts[0].swathNumber)
        metadata['swathNumber'] = ''.join(str(i) for i in sorted(swath_num))

    # calculate ASF frame number for Sentinel-1
    metadata['firstFrameNumber'] = int(np.floor(0.2 * (burst.burstStartUTC -
                                                       master.ascendingNodeTime).total_seconds()))
    metadata['lastFrameNumber']  = int(np.floor(0.2 * (burstEnd.burstStopUTC -
                                                       master.ascendingNodeTime).total_seconds()))

    time_seconds = (burst.burstStartUTC.hour*3600.0 +
                    burst.burstStartUTC.minute*60.0 +
                    burst.burstStartUTC.second)
    metadata['CENTER_LINE_UTC'] = time_seconds

    Vs = np.linalg.norm(burst.orbit.interpolateOrbit(burst.sensingMid, method='hermite').getVelocity())
    metadata['satelliteSpeed'] = Vs
    metadata['azimuthTimeInterval'] = burst.azimuthTimeInterval
    metadata['azimuthPixelSize'] = Vs*burst.azimuthTimeInterval

    tstart = burst.sensingStart
    tend   = burstEnd.sensingStop
    tmid = tstart + 0.5*(tend - tstart)

    orbit = burst.orbit
    peg = orbit.interpolateOrbit(tmid, method='hermite')

    refElp = Planet(pname='Earth').ellipsoid
    llh = refElp.xyz_to_llh(peg.getPosition())
    hdg = orbit.getENUHeading(tmid)
    refElp.setSCH(llh[0], llh[1], hdg)

    metadata['earthRadius'] = refElp.pegRadCur
    metadata['altitude'] = llh[2]

    # make all value in string format
    for key, value in metadata.items():
        metadata[key] = str(value)

    metadata = readfile.standardize_metadata_isce(metadata)
    if rscFile:
        print('writing ', rscFile)
        writefile.write_roipac_rsc(metadata, rscFile)
    return metadata


def load_product(xmlname):
    """Load the product using Product Manager."""
    from iscesys.Component.ProductManager import ProductManager as PM
    pm = PM()
    pm.configure()
    obj = pm.loadProduct(xmlname)
    return obj


def read_baseline(baselineFile):
    bperp, bpar = [], []
    with open(baselineFile, 'r') as f:
        for line in f:
            l = line.split(":")
            if l[0] == "Bperp (average)":
                bperp.append(float(l[1]))
            elif l[0] == "Bpar (average)":
                bpar.append(float(l[1]))
    return np.mean(bperp), np.mean(bpar)


def baseline_timeseries(baselineDir):
    bFiles = sorted(glob.glob(os.path.join(baselineDir, '*/*.txt')))

    bDict={}
    bDict['bperp'] = {}
    bDict['bpar'] = {}
    mDate = None
    for bFile in bFiles:
        mDate, sDate = os.path.basename(bFile).split('.txt')[0].split('_')
        bperp, bpar = read_baseline(bFile)
        bDict['bperp'][sDate] = bperp
        bDict['bpar'][sDate] = bpar

    bDict['bperp'][mDate] = 0
    bDict['bpar'][mDate] = 0
    return bDict


def extract_multilook_number(metadata=dict(), geom_dir=None):
    for fbase in ['hgt','lat','lon','los']:
        fbase = os.path.join(geom_dir, fbase)
        fname = glob.glob('{}*.rdr'.format(fbase)) + glob.glob('{}*.geo'.format(fbase))
        fname = fname[0]
        fullXmlFile = '{}.full.xml'.format(fname)
        if os.path.isfile(fullXmlFile):
            fullXmlDict = readfile.read_isce_xml(fullXmlFile)
            xmlDict = readfile.read_attribute(fname)
            metadata['ALOOKS'] = str(int(int(fullXmlDict['LENGTH']) / int(xmlDict['LENGTH'])))
            metadata['RLOOKS'] = str(int(int(fullXmlDict['WIDTH']) / int(xmlDict['WIDTH'])))
            break
    return metadata


def prepare_geometry(geometryDir, metadata=dict()):
    """Prepare and extract metadata from geometry files"""

    def get_nonzero_row_number(data, buffer=2):
        """Find the first and last row number of rows without zero value, for multiple swaths data"""
        if np.all(data):
            r0, r1 = 0 + buffer, -1 - buffer
        else:
            row_flag = np.sum(data != 0., axis=1) == data.shape[1]
            row_idx = np.where(row_flag)[0]
            r0, r1 = row_idx[0] + buffer, row_idx[-1] - buffer
        return r0, r1

    #print('prepare .rsc file for geometry files')
    isceFiles = [os.path.join(os.path.abspath(geometryDir), '{}.rdr'.format(i)) 
                 for i in ['hgt','lat','lon','los','shadowMask']]
    isceFiles = [i for i in isceFiles if os.path.isfile(i)]

    # get A/RLOOKS
    metadata = extract_multilook_number(metadata, geometryDir)

    # get LAT/LON_REF1/2/3/4 and HEADING_DEG into metadata
    for isceFile in isceFiles:
        if 'lat' in os.path.basename(isceFile):
            data = readfile.read(isceFile)[0]
            r0, r1 = get_nonzero_row_number(data)
            metadata['LAT_REF1'] = str(data[r0, 0])
            metadata['LAT_REF2'] = str(data[r0, -1])
            metadata['LAT_REF3'] = str(data[r1, 0])
            metadata['LAT_REF4'] = str(data[r1, -1])

        if 'lon' in os.path.basename(isceFile):
            data = readfile.read(isceFile)[0]
            r0, r1 = get_nonzero_row_number(data)
            metadata['LON_REF1'] = str(data[r0, 0])
            metadata['LON_REF2'] = str(data[r0, -1])
            metadata['LON_REF3'] = str(data[r1, 0])
            metadata['LON_REF4'] = str(data[r1, -1])

        if 'los' in os.path.basename(isceFile):
            data = readfile.read(isceFile, datasetName='inc')[0]
            data[data == 0.] = np.nan
            az_angle = np.nanmean(data)
            # convert isce azimuth angle to roipac orbit heading angle
            head_angle = -1 * (270 + az_angle)
            head_angle -= np.round(head_angle / 360.) * 360.
            metadata['HEADING_DEG'] = str(head_angle)

    # write rsc file for each geometry file
    for isceFile in isceFiles:
        # prepare metadata for current file
        geom_metadata = readfile.read_attribute(isceFile, meta_ext='.xml')
        geom_metadata.update(metadata)
        geom_metadata = readfile.standardize_metadata_isce(geom_metadata)

        # write .rsc file
        rscFile = isceFile+'.rsc'
        rscDict = dict()
        if os.path.isfile(rscFile):
            rscDict = readfile.read_roipac_rsc(rscFile)
        # update .rsc file only if there are new metadata key/value
        if not set(geom_metadata.items()).issubset(set(rscDict.items())):
            rscDictOut = {**rscDict, **geom_metadata}
            print('writing ', rscFile)
            writefile.write_roipac_rsc(rscDictOut, rscFile)
    return metadata


def prepare_stack(inputDir, filePattern, metadata, baselineDict):
    #print('prepare .rsc file for ', filePattern)
    isceFiles = glob.glob(os.path.join(os.path.abspath(inputDir), '*', filePattern))
    if len(isceFiles) == 0:
        raise FileNotFoundError('no file found in pattern: {}'.format(filePattern))

    # write .rsc file for each interferogram file
    for isceFile in isceFiles:
        # prepare metadata for current file
        ifg_metadata = readfile.read_attribute(isceFile, meta_ext='.xml')
        ifg_metadata.update(metadata)
        dates = os.path.basename(os.path.dirname(isceFile)).split('_')
        ifg_metadata = readfile.standardize_metadata_isce(ifg_metadata, dates, baselineDict)

        # write .rsc file
        rscFile = isceFile+'.rsc'
        rscDict = dict()
        if os.path.isfile(rscFile):
            rscDict = readfile.read_roipac_rsc(rscFile)
        # update .rsc file only if there are new metadata key/value
        if not set(ifg_metadata.items()).issubset(set(rscDict.items())):
            rscDictOut = {**rscDict, **ifg_metadata}
            print('writing ', rscFile)
            writefile.write_roipac_rsc(rscDictOut, rscFile)
    return


#########################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    baselineDict = {}
    if inps.baselineDir:
        baselineDict = baseline_timeseries(inps.baselineDir)

    metadata = {}
    if inps.xmlFile:
        rscFile = os.path.join(os.path.dirname(inps.xmlFile), 'data.rsc')
        if ut.run_or_skip(rscFile, in_file=inps.xmlFile, check_readable=False) == 'run':
            metadata = extract_isce_metadata(inps.xmlFile, rscFile)
        else:
            metadata = readfile.read_roipac_rsc(rscFile)

    if inps.geometryDir:
        metadata = prepare_geometry(inps.geometryDir, metadata)

    if inps.ifgramDir and inps.ifgramFiles:
        for namePattern in inps.ifgramFiles:
            prepare_stack(inps.ifgramDir, namePattern, metadata, baselineDict)
    print('Done.')
    return


#########################################################################
if __name__ == '__main__':
    """Main driver."""
    main() 
