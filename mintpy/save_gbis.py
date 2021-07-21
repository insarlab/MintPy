#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################



import os
import sys
import argparse
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
# suppress UserWarning from matplotlib
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, utils as ut


EXAMPLE = """example:
  save_gbis.py velocity.h5 -g inputs/geometryGeo.h5 -o AlosDT73_20081012_20100302.mat
  save_gbis.py 20150223_20161031_msk.unw -g inputs/geometryGeo.h5 -o Alos2DT23_20150223_20161031.mat
  save_gbis.py 20150223_20161031.unw -g inputs/geometryGeo.h5 --out-data ../Model/data --ellipsoid2geoid
"""

REFERENCE = """references:
  Bagnardi, M., and A. Hooper (2018), Inversion of Surface Deformation Data for Rapid Estimates of Source 
  Parameters and Uncertainties: A Bayesian Approach, Geochemistry, Geophysics, Geosystems, 19, 
  doi:10.1029/2018GC007585.

  Yunjun, Z., Amelung, F., & Aoki, Y. (2021), Imaging the hydrothermal system of Kirishima volcanic complex 
  with L-band InSAR time series, Geophysical Research Letters, 48(11), e2021GL092879. doi:10.1029/2021GL092879
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Convert MintPy product to GBIS .mat format.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('file', help='deformation file.')
    parser.add_argument('dset', nargs='?',
                        help='date/date12 of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-g','--geometry', dest='geom_file', required=True, help='geometry file')
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file.')

    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2,
                        help='custom reference pixel in lat/lon')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                        help='do not display the figure')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('--out-dir', dest='outdir',
                        help='custom output directory, ONLY IF --output is not specified.')
    parser.add_argument('--ellipsoid2geoid', action='store_true',
                        help='Convert the height of ellipsoid to geoid using "geoidheight" module\n'+
                             'Download & install geoidheight as below:\n'+
                             'https://github.com/geodesymiami/2021_Kirishima')
    return parser


def cmd_line_parse(iargs=None):
    print('{} {}'.format(os.path.basename(__file__), ' '.join(iargs)))
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = os.path.abspath(inps.file)

    # Backend setting
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


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
            # velocity to displacement
            date1, date2 = inps.metadata['DATE12'].split('_')
            dt1, dt2 = ptime.date_list2vector([date1, date2])[0]
            inps.phase *= (dt2 - dt1).days / 365.25
            # displacement to phase
            inps.phase *= inps.range2phase

        # update mask to exclude pixel with NaN value
        inps.mask *= ~np.isnan(inps.phase)
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

    # mask out pixels with zero phase value
    ref_y = int(inps.metadata['REF_Y'])
    ref_x = int(inps.metadata['REF_X'])
    inps.mask *= inps.phase != 0
    inps.mask[ref_y, ref_x] = 1
    print('number of pixels after excluding zero phase value: {}'.format(np.sum(inps.mask)))

    # read geometry
    inps.lat, inps.lon = ut.get_lat_lon(inps.metadata)
    inps.inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle')[0]
    inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    inps.height = readfile.read(inps.geom_file, datasetName='height')[0]

    # convert the height of ellipsoid to geoid (mean sea level)
    # ref: https://github.com/vandry/geoidheight
    if inps.ellipsoid2geoid:
        # import geoid module
        try:
            import geoid
        except:
            raise ImportError('Can not import geoidheight!')

        # calculate offset and correct height
        egm_file = os.path.join(os.path.dirname(geoid.__file__), 'geoids/egm2008-1.pgm')
        gh_obj = geoid.GeoidHeight(egm_file)
        h_offset = gh_obj.get(lat=np.nanmean(inps.lat), lon=np.nanmean(inps.lon))
        inps.height -= h_offset

        # print message
        msg = 'convert height from ellipsoid to geoid'
        msg += '\n\tby subtracting a constant offset of {:.2f} m'.format(h_offset)
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
    axs[0].set_title('Phase [{:.1f}, {:.1f}] um'.format(dmin, dmax));
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
    out_fig = '{}.png'.format(os.path.splitext(inps.outfile)[0])
    plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)
    print('saved figure to {}'.format(out_fig))
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
    print('save to file: {}'.format(os.path.abspath(inps.outfile)))
    return


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    read_data(inps)

    plot_data(inps)

    save2mat(inps)

    if inps.disp_fig:
        print('showing...')
        plt.show()
    return inps.outfile


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
