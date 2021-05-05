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
from mintpy.utils import ptime, readfile, arg_group, attribute
from mintpy import subset

d2r = np.pi / 180.
r2d = 180. / np.pi

EXAMPLE = """example:
  ## displacement [event-type inversion]
  # option 1: use velocity file with step estimation from timeseries2velocity.py for co-seismic displacement
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  # option 2: use time-series / ifgramStack file with date1_date2 for the transient displacement:
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  save_kite.py geo/geo_ifgramStack.h5     -d unwrapPhase-20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  ## velocity [interseismic or tensile dislocation inversion]
  # https://pyrocko.org/beat/docs/current/examples/Rectangular_tensile.html
  save_kite.py geo/geo_velocity.h5 -d velocity -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc

  ## import to kite
  spool outfile_name    % /do quadtree,covariance/aps and then File>Save Scene and it is ready for GROND or BEAT
"""

KITE_URL = 'https://github.com/pyrocko/kite'

def create_parser():
    parser = argparse.ArgumentParser(description=f'Generate KITE ({KITE_URL}) npz and yaml from MintPy HDF5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', type=str, help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset', type=str, required=True,
                        help='dataset of interest to be converted.\n'+
                             'e.g.: velocity / stepYYYYMMDD for velocity HDF5 file,\n'+
                             '      date12 in YYYYMMDD_YYYYMMDD for time-series HDF5 file,\n'+
                             '      date12 in unwrapPhase-YYYYMMDD_YYYYMMDD for ifgramStack HDF5 file.')
    parser.add_argument('-g', '--geom', dest='geom_file', type=str,
                        help='geometry file for incidence /azimuth angle and height.')
    parser.add_argument('-m', '--mask', dest='mask_file', type=str,
                        help='mask file, or run mask.py to mask the input file beforehand.')
    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        help='output filename')
    parser = arg_group.add_subset_argument(parser)
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
    try:
        from kite.scene import Scene, SceneConfig
    except ImportError:
        raise ImportError('Can not import kite / pyrocko!')

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
        utc_sec = dt.timedelta(seconds=float(attr['CENTER_LINE_UTC']))
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

    #Save kite container
    scene.save(out_file)

    # print out msg
    urLat = config.frame.llLat + config.frame.dN * float(attr['LENGTH'])
    urLon = config.frame.llLon + config.frame.dE * float(attr['WIDTH'])
    print('\nKite Scene info:')
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
    print('\n---------------SAVING KITE CONTAINER-----------')
    print('Save KITE data in file: {0}.npz {0}.yaml'.format(out_file))
    print('Import to KITE: spool {}'.format(out_file))

    return scene

#########################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    print('\n-------------------READ INPUTS -------------------')
    print('Read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)

    #Extract subset if defined
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), attr)

    # output filename
    if not inps.outfile:
        inps.outfile = attr['PROJECT_NAME']

    # date1/2
    if attr['FILE_TYPE'] in ['timeseries', 'HDFEOS']:
        date1, date2 = inps.dset.split('_')
        inps.dset = date2

    elif attr['FILE_TYPE'] == 'ifgramStack':
        date1, date2 = inps.dset.split('-')[1].split('_')

    else:
        # velocity and *.unw files
        date1, date2 = ptime.yyyymmdd(attr['DATE12'].replace('_','-').split('-'))
        if inps.dset.startswith('step'):
            date1 = inps.dset.split('step')[-1]
            date2 = date1
    print('First  InSAR date: {}'.format(date1))
    print('Second InSAR date: {}'.format(date2))

    # read data
    print('Read {} from file: {}'.format(inps.dset, inps.file))
    dis, attr = readfile.read(inps.file, datasetName=inps.dset, box=inps.pix_box)

    if attr['FILE_TYPE'] == 'timeseries':
        print('Read {} from file: {}'.format(date1, inps.file))
        dis -= readfile.read(inps.file, datasetName=date1, box=inps.pix_box)[0]

    # mask data
    if inps.mask_file is not None:
        mask = readfile.read(inps.mask_file, box=inps.pix_box)[0]
        print('Set data to NaN for pixels with zero value in file: {}'.format(inps.mask_file))
        dis[mask==0] = np.nan

    # read geometry incidence / azimuth angle
    print('\nread incidence / azimuth angle from file: {}'.format(inps.geom_file))
    inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle', box=inps.pix_box)[0]
    az_angle = readfile.read(inps.geom_file, datasetName='azimuthAngle', box=inps.pix_box)[0]
    print('Mean satellite incidence angle: {0:.2f}°'.format(np.nanmean(inc_angle)))
    print('Mean satellite heading   angle: {0:.2f}°\n'.format(90 - np.nanmean(az_angle)))

    # Update attributes
    if inps.subset_lat != None or inps.subset_x != None:
        attr = attribute.update_attribute4subset(attr, inps.pix_box)

    # create kite container
    scene = mintpy2kite(dis, attr, date1, date2, inc_angle, az_angle, out_file=inps.outfile)

    return scene

#########################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
