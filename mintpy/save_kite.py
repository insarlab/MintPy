############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Marin Govorcin, April 2021                       #
############################################################


import datetime as dt

import numpy as np

from mintpy import subset
from mintpy.utils import attribute as attr, ptime, readfile

d2r = np.pi / 180.
r2d = 180. / np.pi


#########################################################################################################
def create_kite_container(dis, atr, date1, date2, inc_angle, az_angle, out_file):
    """Create KITE container.

    Parameters: dis       - 2D np.ndarray for displacement in meters
                atr       - dict, dictionary of mintpy metadata
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
    config.frame.llLat = float(atr['Y_FIRST']) + float(atr['Y_STEP']) * float(atr['LENGTH'])
    config.frame.llLon = float(atr['X_FIRST'])
    config.frame.dE = float(atr['X_STEP'])
    config.frame.dN = float(atr['Y_STEP']) * -1
    config.frame.spacing = 'degree'

    config.meta.scene_title = atr['PROJECT_NAME']
    config.meta.scene_id = atr['trackNumber']

    if date1 is not None:
        utc_sec = dt.timedelta(seconds=float(atr['CENTER_LINE_UTC']))
        config.meta.time_master = dt.datetime.strptime(date1, '%Y%m%d') + utc_sec
        config.meta.time_slave =  dt.datetime.strptime(date2, '%Y%m%d') + utc_sec

    config.meta.orbital_node = atr['ORBIT_DIRECTION'].captilize()
    config.meta.wavelength = atr['WAVELENGTH']
    config.meta.satellite_name = atr['PLATFORM']

    scene = Scene(
        theta = np.flipud(np.pi/2 - inc_angle * d2r),
        phi = np.flipud(az_angle*d2r) + np.pi/2,
        displacement = np.flipud(dis),
        config = config,
    )

    # save kite container
    scene.save(out_file)

    # print out msg
    urLat = config.frame.llLat + config.frame.dN * float(atr['LENGTH'])
    urLon = config.frame.llLon + config.frame.dE * float(atr['WIDTH'])
    print('\nKite Scene info:')
    print('Scene title: {}'.format(config.meta.scene_title))
    print('Scene id: {}'.format(config.meta.scene_id))
    print('Scene orbit: {}'.format(config.meta.orbital_node))
    print('Scene platform: {}'.format(config.meta.satellite_name))
    print('Scene wavelength: {}'.format(config.meta.wavelength))
    print('Scene frame cols / rows: {0:g} / {1:g}'.format(float(atr['WIDTH']), float(atr['LENGTH'])))
    print('Scene frame first / last LAT: {0:.2f} / {1:.2f} °'.format(config.frame.llLat, urLat))
    print('Scene frame first / last LON: {0:.2f} / {1:.2f} °'.format(config.frame.llLon, urLon))
    print('Scene frame spacing in x / y: {} / {} degree'.format(config.frame.dE, config.frame.dN))

    data_list = [scene.displacement, scene.theta * r2d, scene.phi * r2d]
    ds_names = ['displacement', 'incidence angle', 'azimuth angle']
    ds_units = ['m', '°', '°']
    for data, ds_name, ds_unit in zip(data_list, ds_names, ds_units):
        dmin, dmean, dmax = np.nanmin(data), np.nanmean(data), np.nanmax(data)
        print(f'Scene min / mean / max {ds_name:18s}: {dmin:.2f} / {dmean:.2f} / {dmax:.2f} {ds_unit}')

    print('\n---------------SAVING KITE CONTAINER-----------')
    print('Save KITE data in file: {0}.npz {0}.yaml'.format(out_file))
    print('Import to KITE: spool {}'.format(out_file))

    return scene


def save_kite(inps):

    print('\n-------------------READ INPUTS -------------------')
    print('Read metadata from file: {}'.format(inps.file))
    atr = readfile.read_attribute(inps.file)
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), atr)
    inps.outfile = inps.outfile if inps.outfile else atr['PROJECT_NAME']

    # date1/2 and dset
    ftype = atr['FILE_TYPE']
    if ftype in ['timeseries', 'HDFEOS']:
        date1, date2 = inps.dset.split('_')
        inps.dset = date2
    elif ftype == 'ifgramStack':
        date1, date2 = inps.dset.split('-')[1].split('_')
    else:
        # velocity, unw
        date1, date2 = ptime.yyyymmdd(atr['DATE12'].replace('_','-').split('-'))
        if inps.dset.startswith('step'):
            date1 = date2 = inps.dset.split('step')[-1]
    print(f'InSAR start / end date: {date1} / {date2}')

    ## read data
    print('Read {} from file: {}'.format(inps.dset, inps.file))
    dis, atr = readfile.read(inps.file, datasetName=inps.dset, box=inps.pix_box)

    if ftype == 'timeseries':
        print('Read {} from file: {}'.format(date1, inps.file))
        dis -= readfile.read(inps.file, datasetName=date1, box=inps.pix_box)[0]

    # convert radians to meters
    if atr['UNIT'] == 'radian':
        phase2range = float(atr['WAVELENGTH']) / (-4 * np.pi)
        dis *= phase2range

    # mask
    if inps.mask_file is not None:
        mask = readfile.read(inps.mask_file, box=inps.pix_box)[0]
        print('Set data to NaN for pixels with zero value in file: {}'.format(inps.mask_file))
        dis[mask==0] = np.nan

    # read geometry incidence / azimuth angle
    print('\nread incidence / azimuth angle from file: {}'.format(inps.geom_file))
    inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle', box=inps.pix_box)[0]
    az_angle  = readfile.read(inps.geom_file, datasetName='azimuthAngle',   box=inps.pix_box)[0]
    print('Mean LOS incidence angle: {0:.2f}°'.format(np.nanmean(inc_angle)))
    print('Mean LOS azimuth   angle: {0:.2f}°'.format(np.nanmean(az_angle)))

    # update attributes
    if inps.subset_lat is not None or inps.subset_x is not None:
        atr = attr.update_attribute4subset(atr, inps.pix_box)

    ## create kite container
    create_kite_container(
        dis,
        atr,
        date1,
        date2,
        inc_angle,
        az_angle,
        out_file=inps.outfile,
    )

    return
