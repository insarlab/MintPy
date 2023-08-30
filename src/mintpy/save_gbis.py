############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################


import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, utils as ut

# suppress UserWarning from matplotlib
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


##############################################################################
def read_data(inps):
    """
    Returns: defo: 2D np.array with in-valid/masked-out pixel in NaN
    """
    # metadata
    inps.metadata = readfile.read_attribute(inps.file)
    k = inps.metadata['FILE_TYPE']
    inps.range2phase =  -4. * np.pi / float(inps.metadata['WAVELENGTH'])

    # mask
    if inps.mask_file:
        inps.mask = readfile.read(inps.mask_file)[0]
    else:
        inps.mask = np.ones((int(inps.metadata['LENGTH']),
                             int(inps.metadata['WIDTH'])), dtype=np.bool_)

    # data
    if k in ['.unw','velocity']:
        inps.phase = readfile.read(inps.file)[0]
        if k == 'velocity':
            if not inps.dset:
                inps.dset = 'velocity'
                print('No selected dataset, assuming "velocity" and continue.')
            inps.phase, atr = readfile.read(inps.file, datasetName=inps.dset)

            # velocity to displacement
            if inps.dset == 'velocity':
                print('convert velocity to displacement for {}'.format(atr['DATE12']))
                date1, date2 = inps.metadata['DATE12'].split('_')
                dt1, dt2 = ptime.date_list2vector([date1, date2])[0]
                inps.phase *= (dt2 - dt1).days / 365.25

            # meter to phase
            if atr.get('UNIT', 'm/year').startswith('m'):
                print('convert the unit from meter to radian')
                inps.phase *= inps.range2phase
                atr['UNIT'] = 'radian'

        # update mask to exclude pixel with NaN value
        inps.mask *= ~np.isnan(inps.phase)

    else:
        raise ValueError(f"input file not support yet: {k}")
    print(f'number of pixels: {np.sum(inps.mask)}')

    # change reference point
    if inps.ref_lalo:
        coord = ut.coordinate(inps.metadata)
        ref_lat, ref_lon = inps.ref_lalo
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]
        # update data
        inps.phase -= inps.phase[ref_y, ref_x]
        # update metadata
        inps.metadata['REF_LAT'] = ref_lat
        inps.metadata['REF_LON'] = ref_lon
        inps.metadata['REF_Y'] = ref_y
        inps.metadata['REF_X'] = ref_x

    # mask out pixels with zero phase value
    ref_y = int(inps.metadata['REF_Y'])
    ref_x = int(inps.metadata['REF_X'])
    inps.mask *= inps.phase != 0
    inps.mask[ref_y, ref_x] = 1
    print(f'number of pixels after excluding zero phase value: {np.sum(inps.mask)}')

    # read geometry
    inps.lat, inps.lon = ut.get_lat_lon(inps.metadata, geom_file=inps.geom_file)
    inps.inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle')[0]
    inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    inps.height = readfile.read(inps.geom_file, datasetName='height')[0]

    # convert the height of ellipsoid to geoid (mean sea level)
    # ref: https://github.com/vandry/geoidheight
    if inps.ellipsoid2geoid:
        # import geoid module
        try:
            import geoid
        except ImportError:
            raise ImportError('Can not import geoidheight (https://github.com/vandry/geoidheight.git)! ')

        # calculate offset and correct height
        egm_file = os.path.join(os.path.dirname(geoid.__file__), 'geoids/egm2008-1.pgm')
        gh_obj = geoid.GeoidHeight(egm_file)
        h_offset = gh_obj.get(lat=np.nanmean(inps.lat), lon=np.nanmean(inps.lon))
        inps.height -= h_offset

        # print message
        msg = 'convert height from ellipsoid to geoid'
        msg += f'\n\tby subtracting a constant offset of {h_offset:.2f} m'
        print(msg)

    # masking
    inps.phase[inps.mask==0] = np.nan
    inps.lat[inps.mask==0] = np.nan
    inps.lon[inps.mask==0] = np.nan
    inps.inc_angle[inps.mask==0] = np.nan
    inps.head_angle[inps.mask==0] = np.nan
    inps.height[inps.mask==0] = np.nan

    # output filename
    if not inps.outfile:
        proj_name = sensor.project_name2sensor_name(inps.file)[1]
        if not proj_name:
            raise ValueError('No custom/auto output filename found.')
        inps.outfile = '{}_{}.mat'.format(proj_name, inps.metadata['DATE12'])

        if not inps.outdir:
            inps.outdir = os.path.dirname(inps.file)
        inps.outfile = os.path.join(inps.outdir, inps.outfile)
    inps.outfile = os.path.abspath(inps.outfile)
    return


def plot_data(inps):
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=[14, 7])
    axs = axs.flatten()

    # plot deformation
    defo = inps.phase / inps.range2phase * 100. #convert to deformation in cm
    dmin, dmax = np.nanmin(defo), np.nanmax(defo)
    dlim = max(abs(dmin), abs(dmax))
    im = axs[0].imshow(defo, vmin=-dlim, vmax=dlim, cmap='jet', interpolation='nearest')
    # reference point
    axs[0].plot(int(inps.metadata['REF_X']),
                int(inps.metadata['REF_Y']), 'ks', ms=6)
    axs[0].set_title(f'Phase [{dmin:.1f}, {dmax:.1f}] um');
    # colorbar
    cbar = fig.colorbar(im, ax=axs[0]);
    cbar.set_label('cm')

    # plot geometry
    for ax, data, title in zip(axs[1:],
                               [inps.lat, inps.lon, inps.inc_angle, inps.head_angle, inps.height],
                               ['Latitude', 'Longitude', 'Incidence Angle', 'Head Angle', 'Height']):
        im = ax.imshow(data, cmap='jet', interpolation='nearest')
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax)
        if title == 'Height':
            cbar.set_label('m')
        else:
            cbar.set_label('degree')

    # save figure to file
    out_fig = f'{os.path.splitext(inps.outfile)[0]}.png'
    plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)
    print(f'saved figure to {out_fig}')
    return


def save2mat(inps):
    """write mat file"""
    mdict = {}
    # required by GBIS
    mdict['Heading'] = inps.head_angle[inps.mask].reshape(-1,1)
    mdict['Inc'] = inps.inc_angle[inps.mask].reshape(-1,1)
    mdict['Lat'] = inps.lat[inps.mask].reshape(-1,1)
    mdict['Lon'] = inps.lon[inps.mask].reshape(-1,1)
    mdict['Phase'] = inps.phase[inps.mask].reshape(-1,1)
    # optional
    mdict['Height'] = inps.height[inps.mask].reshape(-1,1)
    mdict['Mask'] = inps.mask
    mdict['Metadata'] = inps.metadata
    # save to mat file
    sio.savemat(inps.outfile, mdict, long_field_names=True)
    print(f'save to file: {os.path.abspath(inps.outfile)}')


##############################################################################
def save_gbis(inps):

    # matplotlib backend setting
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    inps.file = os.path.abspath(inps.file)

    read_data(inps)
    plot_data(inps)
    save2mat(inps)

    if inps.disp_fig:
        print('showing...')
        plt.show()

    return
