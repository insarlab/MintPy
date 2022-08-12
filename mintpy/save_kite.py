############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Marin Govorcin, April 2021                       #
############################################################


import datetime as dt
import numpy as np


d2r = np.pi / 180.
r2d = 180. / np.pi


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
