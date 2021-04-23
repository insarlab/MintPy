#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Marin Govorcin, April 2021                       #
############################################################


import sys
import argparse
import datetime as dt
import numpy as np
from mintpy.utils import readfile

try:
    from kite.scene import Scene, SceneConfig
except ImportError:
    raise ImportError('Can not import kite / pyrocko!')


d2r = np.pi / 180.
r2d = 180. / np.pi


EXAMPLE = """example:
  #for Velocity [use step to export coseismic displacement]:
  save_kite.py geo/geo_velocity.h5 -g geo/geo_geometry.h5
  save_kite.py geo/geo_velocity.h5 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  
  #for Timeseries:
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  
  # IMPORT TO kite
  spool outfile_name              % /do quadtree,covariance/aps and then File>Save Scene and it is ready for GROND or BEAT
"""

KITE_URL = 'https://github.com/pyrocko/kite'


def create_parser():
    parser = argparse.ArgumentParser(description=f'Generate KITE ({KITE_URL}) npz and yaml from MintPy HDF5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset',
                        help='velocity step, date of timeseries, or date12 of unw. interferograms to be converted')
    parser.add_argument('-g', '--geom', dest='geom',
                        help='insar geometry; incidence, azimuth, height')
    parser.add_argument('-m', '--mask', dest='mask',
                        help='mask file, use mask.py to create on')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output filename')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps


#########################################################################################################
def mintpy2kite(ifg, attr, date1, date2, inc_angle, az_angle, out_file):
    """Create KITE container.
    Parameters: fig       - 2D np.ndarray for displacement in meters
                attr      - dict, dictionary of mintpy metadata
                date1/2   - str, reference/secondary date in YYYYMMDD format
                inc_angle - 2D np.ndarray, incidence angle of the LOS vector in degree
                az_angle  - 2D np.ndarray, azimutth  angle of the LOS vector in degree
                out_file  - str, path of the output file name of the KITE container
    Returns:    scene     - kite.scene.Scene object
    """

    print('\n---------------PREPARING KITE CONTAINER-----------')
    # fill the Kite container
    config = SceneConfig()
    config.frame.llLat = float(attr['Y_FIRST']) + float(attr['Y_STEP']) * float(attr['LENGTH'])
    config.frame.llLon = float(attr['X_FIRST'])
    config.frame.dE = float(attr['X_STEP'])
    config.frame.dN = float(attr['Y_STEP']) * -1
    config.frame.spacing = 'degree'

    config.meta.scene_title = attr['PROJECT_NAME']
    config.meta.scene_id = attr['trackNumber']

    if date1 is not None:
        utc_sec = dt.timedelta(seconds=float(atr['CENTER_LINE_UTC']))
        config.meta.time_master = dt.datetime.strptime(date1, '%Y%m%d') + utc_sec
        config.meta.time_slave =  dt.datetime.strptime(date2, '%Y%m%d') + utc_sec

    config.meta.orbital_node =  'Ascending' if attr['ORBIT_DIRECTION'].upper().startswith('ASC') else 'Descending'
    config.meta.wavelength = attr['WAVELENGTH']
    config.meta.satellite_name = attr['PLATFORM']

    scene = Scene(
        theta = np.flipud(np.pi/2 - inc_angle * d2r),
        phi = np.flipud(az_angle*d2r) + np.pi/2,
        displacement = np.flipud(ifg),
        config = config,
    )

    print('\n---------------SAVING KITE CONTAINER-----------')
    print('Save KITE data in file: {}'.format(out_file))
    scene.save(out_file)

    # print out msg
    urLat = config.frame.llLat + config.frame.dN * float(attr['LENGTH'])
    urLon = config.frame.llLon + config.frame.dE * float(attr['WIDTH'])
    print('Kite Scene info:')
    print('Scene title: {}'.format(config.meta.scene_title))
    print('Scene id: {}'.format(config.meta.scene_id))
    print('Scene orbit: {}'.format(config.meta.orbital_node))
    print('Scene platform: {}'.format(config.meta.satellite_name))
    print('Scene wavelength: {}'.format(config.meta.wavelength))
    print('Scene frame cols / rows: {0:g} / {1:g}'.format(float(attr['WIDTH']), float(attr['LENGTH'])))
    print('Scene frame first / last LAT: {0:.2f} / {1:.2f} °'.format(config.frame.llLat, urLat))
    print('Scene frame first / last LON: {0:.2f} / {1:.2f} °'.format(config.frame.llLon, urLon))
    print('Scene frame spacing in x / y: {} / {} {}'.format(config.frame.dE, config.frame.dN, config.frame.spacing))
    print('Scene min / mean / max displacement    : {0:.2f} / {1:.2f} / {2:.2f} m'.format(np.nanmin(scene.displacement),
                                                                                          np.nanmean(scene.displacement),
                                                                                          np.nanmax(scene.displacement)))
    print('Scene min / mean / max incidence angle : {0:.2f} / {1:.2f} / {2:.2f} °'.format(np.nanmin(scene.theta)*r2d,
                                                                                          np.nanmean(scene.theta)*r2d,
                                                                                          np.nanmax(scene.theta)*r2d))
    print('Scene min / mean / max azimuth angle   : {0:.2f} / {1:.2f} / {2:.2f} °'.format(np.nanmin(scene.phi)*r2d,
                                                                                          np.nanmean(scene.phi)*r2d,
                                                                                          np.nanmax(scene.phi)*r2d))

    return scene


#########################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    print('\n-------------------READ INPUTS -------------------')
    # check reference date
    print('Read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)
    if attr['FILE_TYPE'] == 'timeseries' and inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date, inps.dset = None, None

    # read data
    print('Read data from file: {}'.format(inps.file))
    dis, attr = readfile.read(inps.file, datasetName=inps.dset)
    if attr['FILE_TYPE'] == 'timeseries' and inps.ref_date:
        dis -= readfile.read(inps.file, datasetName=inps.ref_date)[0]

    #Mask data
    if inps.mask is not None:
        mask = readfile.read(inps.mask)[0]
        print('Masking data')
        dis[mask==0] = np.nan

    if inps.ref_date is not None:
        print('\nFirst  InSAR date: {}'.format(inps.ref_date))
        print('Second InSAR date: {}'.format(inps.dset))

    # output filename
    if not inps.outfile:
        inps.outfile = attr['PROJECT_NAME']

    # read geometry incidence / azimuth angle
    print('\nread incidence / azimuth angle from file: {}'.format(inps.geom))
    inc_angle = readfile.read(inps.geom, datasetName='incidenceAngle')[0]
    az_angle = readfile.read(inps.geom, datasetName='azimuthAngle')[0]
    print('Mean satellite incidence angle; {0:.2f}°'.format(np.nanmean(inc_angle)))
    print('Mean satellite heading angle: {0:.2f}°'.format(90 - np.nanmean(az_angle)))

    # create kite container
    scene = mintpy2kite(dis, attr, inps.ref_date, inps.dset, inc_angle, az_angle, out_file=inps.out_file)

    return


#########################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
