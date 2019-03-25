#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import argparse
from pysar.utils import utils as ut


###########################################################################################
TEMPLATE = """template:
## remove phase ramp for each epoch, useful to check localized deformation, i.e. volcanic, land subsidence, etc.
## [linear, quadratic]
pysar.deramp          = auto  #[no / linear / quadratic], auto for no - no ramp will be removed
pysar.deramp.maskFile = auto  #[filename / no], auto for maskTempCoh.h5, mask file for ramp estimation
"""

EXAMPLE = """example:
  remove_ramp.py  timeseries.h5      -m maskTempCoh.h5
  remove_ramp.py  ifgramStack.h5     -m maskTempCoh.h5  -d unwrapPhase_bridging
  remove_ramp.py  090214_101120.unw  -m maskTempCoh.h5  -s quadratic
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
    parser.add_argument('-d','--dset', dest='dset', 
                        help='dataset name to be derampped in ifgramStack file\n' + 
                             'e.g.: unwrapPhase\n' +
                             '      unwrapPhase_bridging')
    parser.add_argument('-o', '--outfile', help='Output file name. Disabled when more than 1 input files')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file)
    return inps


###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('input file: ({})\n{}'.format(len(inps.file), inps.file))

    for fname in inps.file:
        print('------------------------------------------')
        ut.run_deramp(fname,
                      ramp_type=inps.surface_type,
                      mask_file=inps.mask_file,
                      out_file=inps.outfile,
                      datasetName=inps.dset)
    return


###########################################################################################
if __name__ == '__main__':
    main()
