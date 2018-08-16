#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
from pysar.utils import readfile, writefile, utils as ut, deramp


######################################
TEMPLATE = """template:
## remove phase ramp for each epoch, useful to check localized deformation, i.e. volcanic, land subsidence, etc.
## [linear, quadratic]
pysar.deramp          = auto  #[no / linear / quadratic], auto for no - no ramp will be removed
pysar.deramp.maskFile = auto  #[filename / no], auto for maskTempCoh.h5, mask file for ramp estimation
"""

EXAMPLE = """example:
  remove_ramp.py  timeseries.h5      -m Mask.h5
  remove_ramp.py  timeseries.h5      -m Mask.h5         -s quadratic
  remove_ramp.py  090214_101120.unw  -m Mask_tempCoh.h5 -s quadratic  -y 0,2400,2000,6843
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Remove phase ramp',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE + '\n' + EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) for ramp removal')
    parser.add_argument('-m', '--mask', dest='mask_file', default='maskTempCoh.h5',
                        help='mask for pixels used in ramp estimation\n' +
                             'default - maskTempCoh.h5\n' +
                             'no - use the whole area')
    parser.add_argument('-s', dest='surface_type', default='linear',
                        choices={'linear', 'quadratic',
                                 'linear_range', 'linear_azimuth',
                                 'quadratic_range', 'quadratic_azimuth'},
                        help='type of surface/ramp to remove, linear by default')

    parser.add_argument('-y', dest='ysub', type=int, nargs='*',
                        help='subset in azimuth/row direction\h' +
                        ' for multiple surface removal within one track, i.e.:\n' +
                        '0,2400,2000,6843')

    parser.add_argument('-o', '--outfile', help='Output file name. Disabled when more than 1 input files')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file)

    if inps.ysub and not len(inps.ysub) % 2 == 0:
        raise Exception('ERROR: -y input has to have even length!')
    return inps


######################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('input file: ({})\n{}'.format(len(inps.file), inps.file))

    for File in inps.file:
        print('------------------------------------------')
        deramp.remove_surface(File, inps.surface_type, inps.mask_file, ysub=inps.ysub)

    print('Done.')
    return


###########################################################################################
if __name__ == '__main__':
    main()
