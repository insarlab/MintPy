#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2019, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################



import os
import argparse
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, utils as ut


EXAMPLE = """example:
  save_gbis.py velocity.h5 -g inputs/geometryGeo.h5 -o AlosDT73_20081012_20100302.mat
  save_gbis.py 20150223_20161031_msk.unw -g inputs/geometryGeo.h5 -o Alos2DT23_20150223_20161031.mat
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Convert MintPy product to GBIS .mat format.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='deformation file.')
    parser.add_argument('dset', nargs='?',
                        help='date/date12 of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-g','--geometry', dest='geom_file', required=True, help='geometry file')
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file.')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2,
                        help='custom reference pixel in lat/lon')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def read_data(inps):
    """
    Returns: defo: 2D np.array with in-valid/masked-out pixel in NaN
    """
    # metadata
    inps.metadata = readfile.read_attribute(inps.file)
    k = inps.metadata['FILE_TYPE']
    range2phase =  -4. * np.pi / float(inps.metadata['WAVELENGTH'])
    ext = os.path.splitext(inps.file)[1]

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
            # velocity to displacement
            date1, date2 = inps.metadata['DATE12'].split('_')
            dt1, dt2 = ptime.date_list2vector([date1, date2])[0]
            tdiff = (dt2 - dt1).days / 365.25
            inps.phase *= tdiff
            # displacement to phase
            inps.phase *= range2phase

        # update mask to exclude pixel with NaN value
        inps.mask *= ~np.isnan(inps.phase)
        # set all masked out pixel to NaN
        inps.phase[inps.mask==0] = np.nan
    else:
        raise ValueError("input file not support yet: {}".format(k))
    print('number of pixels: {}'.format(np.sum(inps.mask)))

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

    # read geometry
    inps.lat, inps.lon = ut.get_lat_lon(inps.metadata)
    inps.inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle')[0]
    inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    inps.lat[inps.mask==0] = np.nan
    inps.lon[inps.mask==0] = np.nan
    inps.inc_angle[inps.mask==0] = np.nan
    inps.head_angle[inps.mask==0] = np.nan

    # output filename
    if not inps.outfile:
        SAT = inps.metadata['PLATFORM'].lower().capitalize()
        if SAT.lower() not in sensor.sensorNames:
            raise ValueError('un-recognized sensor name: {}'.format(SAT))
        if inps.metadata['ORBIT_DIRECTION'].lower().startswith('asc'):
            ORBIT = 'A'
        else:
            ORBIT = 'D'
        TRACK = inps.metadata['trackNumber']
        inps.outfile = '{}{}T{}_{}.mat'.format(SAT, ORBIT, TRACK, inps.metadata['DATE12'])
    inps.outfile = os.path.abspath(inps.outfile)
    return


def plot_data(inps):
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=[14, 7])
    axs = axs.flatten()

    im = axs[0].imshow(ut.wrap(inps.phase), vmin=-np.pi, vmax=np.pi, cmap='jet');
    axs[0].set_title('Phase (wrapped for display)');
    cbar = fig.colorbar(im, ax=axs[0]);
    cbar.set_label('radian')

    for ax, data, title in zip(axs[1:5],
                               [inps.lat, inps.lon, inps.inc_angle, inps.head_angle],
                               ['Latitude', 'Longitude', 'Incidence Angle', 'Head Angle']):
        im = ax.imshow(data, cmap='jet')
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('degree')

    axs[5].axis('off')

    # save figure to file
    out_fig = '{}.png'.format(os.path.splitext(inps.outfile)[0])
    plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)
    print('saved figure to {}'.format(out_fig))
    return


def save2mat(inps):

    # 4. write mat file
    mdict = {}
    mdict['Lon'] = inps.lon[inps.mask].reshape(-1,1)
    mdict['Lat'] = inps.lat[inps.mask].reshape(-1,1)
    mdict['Phase'] = inps.phase[inps.mask].reshape(-1,1)
    mdict['Inc'] = inps.inc_angle[inps.mask].reshape(-1,1)
    mdict['Heading'] = inps.head_angle[inps.mask].reshape(-1,1)
    mdict['metadata'] = inps.metadata
    sio.savemat(inps.outfile, mdict, long_field_names=True)
    print('save to file: {}.mat'.format(os.path.abspath(inps.outfile)))
    return


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    read_data(inps)

    plot_data(inps)

    save2mat(inps)

    plt.show()
    return inps.outfile


##########################################################################
if __name__ == '__main__':
    main()
