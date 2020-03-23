#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import time
import shutil
import argparse
import h5py
import numpy as np
from mintpy.objects import timeseries
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, writefile, ptime, utils as ut


##################################################################
TEMPLATE = get_template_content('reference_date')

EXAMPLE = """example:
  reference_date.py timeseries.h5 timeseries_ERA5.h5 timeseries_ERA5_demErr.h5 --template smallbaselineApp.cfg
  reference_date.py timeseries_ERA5_demErr.h5 --ref-date 20050107
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Change reference date of timeseries.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', nargs='+', help='timeseries file(s)')
    parser.add_argument('-r', '--ref-date', dest='refDate', default='minRMS',
                        help='reference date or method, default: auto. e.g.\n' +
                             '20101120\n' +
                             'reference_date.txt - text file with date in YYYYMMDD format in it\n' +
                             'minRMS             - choose date with min residual standard deviation')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    atr = readfile.read_attribute(inps.timeseries_file[0])
    if 'timeseries' not in atr['FILE_TYPE'].lower():
        raise ValueError('input file type: {} is not timeseries.'.format(atr['FILE_TYPE']))
    return inps


def read_template2inps(templateFile, inps=None):
    """Update inps with options from templateFile"""
    if not inps:
        inps = cmd_line_parse()
    template = readfile.read_template(templateFile)
    template = ut.check_template_auto_value(template)

    key = 'mintpy.reference.date'
    if key in template.keys() and template[key]:
        inps.refDate = template[key]
    return inps


def read_ref_date(inps):
    # check input reference date
    if not inps.refDate:
        print('No reference date input, skip this step.')
        return None
    # string in digit
    elif inps.refDate.isdigit():
        pass
    else:
        # txt file
        if os.path.isfile(inps.refDate):
            print('read reference date from file: ' + inps.refDate)
            inps.refDate = ptime.read_date_txt(inps.refDate)[0]
        else:
            print('input file {} does not exist, skip this step'.format(inps.refDate))
            return None
    inps.refDate = ptime.yyyymmdd(inps.refDate)
    print('input reference date: {}'.format(inps.refDate))

    # check available dates
    date_list = timeseries(inps.timeseries_file[0]).get_date_list()
    if inps.refDate not in date_list:
        msg = 'input reference date: {} is not found.'.format(inps.refDate)
        msg += '\nAll available dates:\n{}'.format(date_list)
        raise Exception(msg)
    return inps.refDate


##################################################################
def ref_date_file(ts_file, ref_date, outfile=None):
    """Change input file reference date to a different one.
    Parameters: ts_file : str, timeseries file to be changed
                ref_date : str, date in YYYYMMDD format
                outfile  : if str, save to a different file
                           if None, modify the data value in the existing input file
    """
    print('-'*50)
    print('change reference date for file: {}'.format(ts_file))
    atr = readfile.read_attribute(ts_file)
    if ref_date == atr['REF_DATE']:
        print('same reference date chosen as existing reference date.')
        if not outfile:
            print('Nothing to be done.')
            return ts_file
        else:
            print('Copy {} to {}'.format(ts_file, outfile))
            shutil.copy2(ts_file, outfile)
            return outfile
    else:
        obj = timeseries(ts_file)
        obj.open(print_msg=False)
        ref_idx = obj.dateList.index(ref_date)
        print('reading data ...')
        ts_data = readfile.read(ts_file)[0]

        ts_data -= np.tile(ts_data[ref_idx, :, :].reshape(1, obj.length, obj.width), (obj.numDate, 1, 1))

        if not outfile:
            print('open {} with r+ mode'.format(ts_file))
            with h5py.File(ts_file, 'r+') as f:
                print("update /timeseries dataset and 'REF_DATE' attribute value")
                f['timeseries'][:] = ts_data
                f.attrs['REF_DATE'] = ref_date
            print('close {}'.format(ts_file))
        else:
            atr['REF_DATE'] = ref_date
            writefile.write(ts_data, outfile, metadata=atr, ref_file=ts_file)
    return outfile


##################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    inps.refDate = read_ref_date(inps)

    if inps.refDate:
        for ts_file in inps.timeseries_file:
            ref_date_file(ts_file, inps.refDate, inps.outfile)
            time.sleep(1)   #to distinguish the modification time of input files
    return


##################################################################
if __name__ == '__main__':
    main()
