#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author:  Heresh Fattahi, 2013                            #
############################################################
# Based on the empirical model of Marinkovic and Larsen (2013),
# the LOD correction for each pixel is given in meters as:
#     C = 3.87e-7 * x * range_pixel_size_in_meter * time_difference_in_year
# where x is the pixel count in range direction


import os
import sys
import numpy as np

from mintpy.objects import timeseries
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, writefile, ptime
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################################
TEMPLATE = get_template_content('correct_LOD')

REFERENCE = """reference:
  Marinkovic, P., and Y. Larsen (2013), Consequences of long-term ASAR local oscillator 
  frequency decay - An empirical study of 10 years of data, in Living Planet Symposium,
  Edinburgh, U.K.
"""

EXAMPLE = """example:
  local_oscilator_drift.py  timeseries.h5                 inputs/geometryRadar.h5
  local_oscilator_drift.py  filt_101020_110220_4rlks.unw  inputs/geometryRadar.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Local Oscilator Drift (LOD) correction of Envisat'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument(dest='file', help='timeseries / interferograms file, i.e. timeseries.h5')
    parser.add_argument(dest='range_dist_file',
                        help='Slant range distance file, i.e. inputs/geometryRadar.h5, inputs/geometryGeo.h5\n' +
                        'or use range_distance.py to generate it.')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='Output file name for corrected file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#########################################################################################
def get_relative_range_distance(meta):
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    range_dist_1d = float(meta['RANGE_PIXEL_SIZE']) * np.linspace(0, width-1, width)
    range_dist = np.tile(range_dist_1d, (length, 1))
    range_dist -= range_dist[int(meta['REF_Y']), int(meta['REF_X'])]
    return range_dist


def correct_local_oscilator_drift(fname, rg_dist_file=None, out_file=None):
    print('-'*50)
    print('correct Local Oscilator Drift for Envisat using an empirical model (Marinkovic and Larsen, 2013)')
    print('-'*50)
    atr = readfile.read_attribute(fname)

    # Check Sensor Type
    platform = atr['PLATFORM']
    print('platform: '+platform)
    if not platform.lower() in ['env', 'envisat']:
        print('No need to correct LOD for '+platform)
        return

    # output file name
    if not out_file:
        out_file = '{}_LOD{}'.format(os.path.splitext(fname)[0], os.path.splitext(fname)[1])

    # Get LOD ramp rate from empirical model
    if not rg_dist_file:
        print('calculate range distance from file metadata')
        rg_dist = get_relative_range_distance(atr)
    else:
        print('read range distance from file: %s' % (rg_dist_file))
        rg_dist = readfile.read(rg_dist_file, datasetName='slantRangeDistance', print_msg=False)[0]
        rg_dist -= rg_dist[int(atr['REF_Y']), int(atr['REF_X'])]

    ramp_rate = np.array(rg_dist * 3.87e-7, np.float32)

    # Correct LOD Ramp for Input fname
    range2phase = -4*np.pi / float(atr['WAVELENGTH'])
    k = atr['FILE_TYPE']
    if k == 'timeseries':
        # read
        obj = timeseries(fname)
        obj.open()
        data = obj.read()

        # correct LOD
        diff_year = np.array(obj.yearList)
        diff_year -= diff_year[obj.refIndex]
        for i in range(data.shape[0]):
            data[i, :, :] -= ramp_rate * diff_year[i]

        # write
        obj_out = timeseries(out_file)
        obj_out.write2hdf5(data, refFile=fname)

    elif k in ['.unw']:
        data, atr = readfile.read(fname)

        dates = ptime.yyyymmdd2years(ptime.yyyymmdd(atr['DATE12'].split('-')))
        dt = dates[1] - dates[0]
        data -= ramp_rate * range2phase * dt

        writefile.write(data, out_file=out_file, metadata=atr)
    else:
        print('No need to correct for LOD for %s file' % (k))
    return out_file


#########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    inps.outfile = correct_local_oscilator_drift(inps.file, inps.range_dist_file, inps.outfile)
    print('Done.')
    return


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
