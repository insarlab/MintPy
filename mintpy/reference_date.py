#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import time
import shutil
import argparse
import h5py
import numpy as np

from mintpy.objects import timeseries
from mintpy.objects.cluster import split_box2sub_boxes
from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_group, readfile, writefile, ptime, utils as ut


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
                             'time-series HDF5 file with REF_DATE in its attributes\n' +
                             'reference_date.txt - text file with date in YYYYMMDD format in it\n' +
                             'minRMS             - choose date with min residual standard deviation')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('--force', action='store_true',
                        help='Force updating the data matrix.')

    # computing
    parser = arg_group.add_memory_argument(parser)

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

    key = 'mintpy.compute.maxMemory'
    if key in template.keys() and template[key]:
        inps.maxMemory = float(template[key])

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
        if os.path.isfile(inps.refDate):
            print('read reference date from file: ' + inps.refDate)
            if inps.refDate.endswith('.h5'):
                # HDF5 file
                atr = readfile.read_attribute(inps.refDate)
                inps.refDate = atr['REF_DATE']
            else:
                # txt file
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
def change_timeseries_ref_date(ts_file, ref_date, outfile=None, max_memory=4.0, force=False):
    """Change input file reference date to a different one.
    Parameters: ts_file : str, timeseries file to be changed
                ref_date : str, date in YYYYMMDD format
                outfile  : if str, save to a different file
                           if None, modify the data value in the existing input file
    """
    ts_file = os.path.abspath(ts_file)
    if not outfile:
        outfile = ts_file
    outfile = os.path.abspath(outfile)

    print('-'*50)
    print('change reference date for file: {}'.format(ts_file))
    atr = readfile.read_attribute(ts_file)
    dsName = atr['FILE_TYPE']

    # if the input reference date is the same as the existing one.
    if ref_date == atr.get('REF_DATE', None) and not force:
        print('input refDate is the same as the existing REF_DATE.')
        if outfile == ts_file:
            print('Nothing to be done.')
            return ts_file
        else:
            print('Copy {} to {}'.format(ts_file, outfile))
            shutil.copy2(ts_file, outfile)
            return outfile

    # basic info
    obj = timeseries(ts_file)
    obj.open(print_msg=False)
    num_date = obj.numDate
    length = obj.length
    width = obj.width
    ref_idx = obj.dateList.index(ref_date)

    # get list of boxes for block-by-block IO
    num_box = int(np.ceil((num_date * length * width * 4 * 2) / (max_memory * 1024**3)))
    box_list = split_box2sub_boxes(box=(0, 0, width, length),
                                   num_split=num_box,
                                   dimension='y',
                                   print_msg=True)

    # updating existing file or write new file
    if outfile == ts_file:
        mode = 'r+'
    else:
        mode = 'a'
        # instantiate output file
        writefile.layout_hdf5(outfile, ref_file=ts_file)

    # loop for block-by-block IO
    for i, box in enumerate(box_list):
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_width))
            print('box length: {}'.format(box_length))

        # reading
        print('reading data ...')
        ts_data = readfile.read(ts_file, box=box)[0]

        print('referencing in time ...')
        dshape = ts_data.shape
        ts_data -= np.tile(ts_data[ref_idx, :, :].reshape(1, dshape[1], dshape[2]), (dshape[0], 1, 1))

        # writing
        block = (0, num_date, box[1], box[3], box[0], box[2])
        writefile.write_hdf5_block(outfile,
                                   data=ts_data,
                                   datasetName=dsName,
                                   block=block,
                                   mode=mode)

    # update metadata
    print('update "REF_DATE" attribute value to {}'.format(ref_date))
    with h5py.File(outfile, 'r+') as f:
        f.attrs['REF_DATE'] = ref_date
        f.attrs['FILE_PATH'] = outfile

    return outfile


##################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # read reference date
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    inps.refDate = read_ref_date(inps)

    # run referencing in time
    if inps.refDate:
        for ts_file in inps.timeseries_file:
            change_timeseries_ref_date(ts_file,
                                       ref_date=inps.refDate,
                                       outfile=inps.outfile,
                                       max_memory=inps.maxMemory,
                                       force=inps.force)

            #to distinguish the modification time of input files
            time.sleep(1)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    return


##################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
