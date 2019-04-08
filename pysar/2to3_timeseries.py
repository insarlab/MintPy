#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018                              #
############################################################


import h5py
import argparse
import numpy as np
from pysar.objects import timeseries
from pysar.utils import ptime, readfile


################################################################################
EXAMPLE = """example:
  2to3_timeseries.py  timeseries_ECMWF_demErr_refDate_plane.h5  timeseries_ECMWF_demErr_ramp.h5
"""


def create_parser():
    """ Command line parser """
    parser = argparse.ArgumentParser(description='Convert time-series file from py2 PYSAR to py3 PYSAR format.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted')
    parser.add_argument('-o', '--output', dest='outfile', required=True, help='output file name')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################################
def run_2to3_timeseries(py2_file, py3_file):
    """Convert timeseries file from py2 pysar format to py3 pysar format"""
    # read data from py2_file
    atr = readfile.read_attribute(py2_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    with h5py.File(py2_file, 'r') as f:
        date_list = list(f['timeseries'].keys())
        num_date = len(date_list)
        ts_data = np.zeros((num_date, length, width), np.float32)
        print('reading time-series ...')
        prog_bar = ptime.progressBar(maxValue=num_date)
        for i in range(num_date):
            ts_data[i, :, :] = f['timeseries/{}'.format(date_list[i])][:]
            prog_bar.update(i+1, suffix=date_list[i])
        prog_bar.close()

    # prepare metadata
    bperp = np.array([float(i) for i in atr['P_BASELINE_TIMESERIES'].split()], dtype=np.float32)
    dates = np.array(date_list, np.string_)
    atr['REF_DATE'] = date_list[0]
    for key in ['P_BASELINE_TIMESERIES', 
                'P_BASELINE_TOP_TIMESERIES',
                'P_BASELINE_BOTTOM_TIMESERIES']:
        try:
            atr.pop(key)
        except:
            pass

    # write to py3_file
    ts_obj = timeseries(py3_file)
    ts_obj.write2hdf5(data=ts_data, dates=dates, bperp=bperp, metadata=atr)
    return py3_file


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    run_2to3_timeseries(inps.file, inps.outfile)
    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main()
