#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import argparse
import random
import sys
from datetime import datetime as dt

import h5py
import numpy as np

from mintpy.utils import ptime, readfile

##############################################################################################
EXAMPLE = """example:
  ifgram_simulation.py  unwrapIfgram.h5  velocity.h5
  ifgram_simulation.py  unwrapIfgram.h5  velocity.h5  -p 0.2  -m mask_aoi.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Simulating a set of interferograms based on ' +
                                                 'real interferograms and an existing displacement velocity field.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument(
        'ifgram_file', help='real interferograms file, e.g. unwrapIfgram.h5')
    parser.add_argument('velocity_file', help='velocity file')
    parser.add_argument('-o', '--output', dest='outfile', default='simulated_unwrapIfgram.h5',
                        help='output filename for simulated interferograms file. Default: simulated_unwrapIfgram.h5')
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2, metavar=('XMIN', 'XMAX'),
                        help='subset for simulation in x/cross-track/range direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2, metavar=('YMIN', 'YMAX'),
                        help='subset for simulation in y/along-track/azimuth direction')

    unwErr = parser.add_argument_group(
        'Add unwrapping error to the simulation.')
    unwErr.add_argument('-p', '--percentage', type=float, default=0.0,
                        help='percentage of unwrapping error, [0.0-1.0]. Default: 0.0')
    unwErr.add_argument('-m', '--mask', dest='mask_file', default='mask.h5',
                        help='mask for pixels with unwrapping error. Default: mask.h5')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if not 0.0 <= inps.percentage <= 1.0:
        raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % inps.percentage)
    return inps


##############################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    atr = readfile.read_attribute(inps.velocity_file)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    # Check subset input
    if inps.subset_y:
        inps.subset_y = sorted(inps.subset_y)
        print('subset in y/azimuth direction: '+str(inps.subset_y))
    else:
        inps.subset_y = [0, length]

    if inps.subset_x:
        inps.subset_x = sorted(inps.subset_x)
        print('subset in x/range direction: '+str(inps.subset_x))
    else:
        inps.subset_x = [0, width]
    y0, y1 = inps.subset_y
    x0, x1 = inps.subset_x

    # Read velocity/rate
    velocity = readfile.read(inps.velocity_file)[0]
    print('read velocity file: '+inps.velocity_file)

    k = 'interferograms'
    h5 = h5py.File(inps.ifgram_file, 'r')
    ifgram_list = sorted(h5[k].keys())
    ifgram_num = len(ifgram_list)
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    print('number of interferograms: '+str(ifgram_num))

    # Select interferograms with unwrapping error
    if inps.percentage > 0.0:
        mask = readfile.read(inps.mask_file, datasetName='mask')[0]
        print('read mask for pixels with unwrapping error from file: '+inps.mask_file)

        unw_err_ifgram_num = int(np.rint(inps.percentage*ifgram_num))
        unw_err_ifgram_idx = random.sample(list(range(ifgram_num)), unw_err_ifgram_num)
        unw_err_ifgram_list = [ifgram_list[i] for i in unw_err_ifgram_idx]
        unw_err_date12_list = [date12_list[i] for i in unw_err_ifgram_idx]
        print('randomly choose the following %d interferograms with unwrapping error' % unw_err_ifgram_num)
        print(unw_err_date12_list)

        unit_unw_err = 2.0*np.pi*mask
    else:
        unw_err_ifgram_list = []

    # Generate simulated interferograms
    m_dates = ptime.yyyymmdd([i.split('_')[0] for i in date12_list])
    s_dates = ptime.yyyymmdd([i.split('_')[1] for i in date12_list])
    range2phase = -4.0*np.pi/float(atr['WAVELENGTH'])

    print('writing simulated interferograms file: '+inps.outfile)
    h5out = h5py.File(inps.outfile, 'w')
    group = h5out.create_group('interferograms')
    for i in range(ifgram_num):
        ifgram = ifgram_list[i]
        # Get temporal baseline in years
        t1 = dt.strptime(m_dates[i], "%Y%m%d")
        t2 = dt.strptime(s_dates[i], "%Y%m%d")
        dt = (t2 - t1)
        dt = float(dt.days) / 365.25

        # Simuated interferograms with unwrap error
        unw = velocity*dt*range2phase
        if ifgram in unw_err_ifgram_list:
            rand_int = random.sample(list(range(1, 10)), 1)[0]
            unw += rand_int * unit_unw_err
            print(ifgram+'  - add unwrapping error of %d*2*pi' % rand_int)
        else:
            print(ifgram)

        gg = group.create_group(ifgram)
        dset = gg.create_dataset(ifgram, data=unw[y0:y1, x0:x1])

        for key, value in h5[k][ifgram].attrs.items():
            gg.attrs[key] = value
        if ifgram in unw_err_ifgram_list:
            gg.attrs['unwrap_error'] = 'yes'
        else:
            gg.attrs['unwrap_error'] = 'no'
        gg.attrs['LENGTH'] = y1-y0
        gg.attrs['WIDTH'] = x1-x0
    h5.close()
    h5out.close()
    print('Done.')
    return inps.outfile


##############################################################################################
if __name__ == '__main__':
    main()
