############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import shutil
import time

import h5py
import numpy as np

from mintpy.objects import timeseries
from mintpy.objects.cluster import split_box2sub_boxes
from mintpy.utils import ptime, readfile, writefile


##################################################################
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
            print(f'input file {inps.refDate} does not exist, skip this step')
            return None
    inps.refDate = ptime.yyyymmdd(inps.refDate)
    print(f'input reference date: {inps.refDate}')

    # check available dates
    date_list = timeseries(inps.timeseries_file[0]).get_date_list()
    if inps.refDate not in date_list:
        msg = f'input reference date: {inps.refDate} is not found.'
        msg += f'\nAll available dates:\n{date_list}'
        raise Exception(msg)

    return inps.refDate


##################################################################
def change_timeseries_ref_date(ts_file, ref_date, outfile=None, max_memory=4.0, force=False):
    """Change input file reference date to a different one.
    Parameters: ts_file  - str, timeseries file to be changed
                ref_date - str, date in YYYYMMDD format
                outfile  - if str, save to a different file
                           if None, modify the data value in the existing input file
    """
    ts_file = os.path.abspath(ts_file)
    outfile = outfile if outfile else ts_file
    outfile = os.path.abspath(outfile)

    print('-'*50)
    print(f'change reference date for file: {ts_file}')
    atr = readfile.read_attribute(ts_file)
    dsName = atr['FILE_TYPE']

    # if the input reference date is the same as the existing one.
    if ref_date == atr.get('REF_DATE', None) and not force:
        print('input refDate is the same as the existing REF_DATE.')
        if outfile == ts_file:
            print('Nothing to be done.')
            return ts_file
        else:
            print(f'Copy {ts_file} to {outfile}')
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
    box_list, num_box = split_box2sub_boxes(
        box=(0, 0, width, length),
        num_split=num_box,
        dimension='y',
        print_msg=True,
    )

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
            print(f'\n------- processing patch {i+1} out of {num_box} --------------')
            print(f'box width:  {box_width}')
            print(f'box length: {box_length}')

        # reading
        print('reading data ...')
        ts_data = readfile.read(ts_file, box=box)[0]

        print('referencing in time ...')
        dshape = ts_data.shape
        ts_data -= np.tile(ts_data[ref_idx, :, :].reshape(1, dshape[1], dshape[2]),
                           (dshape[0], 1, 1))

        # writing
        block = (0, num_date, box[1], box[3], box[0], box[2])
        writefile.write_hdf5_block(
            outfile,
            data=ts_data,
            datasetName=dsName,
            block=block,
            mode=mode,
        )

    # update metadata
    print(f'update "REF_DATE" attribute value to {ref_date}')
    with h5py.File(outfile, 'r+') as f:
        f.attrs['REF_DATE'] = ref_date
        f.attrs['FILE_PATH'] = outfile

    return outfile


##################################################################
def run_reference_date(inps):
    start_time = time.time()

    # read reference date
    inps.refDate = read_ref_date(inps)

    # apply temporal referencing
    if inps.refDate:
        for ts_file in inps.timeseries_file:
            change_timeseries_ref_date(
                ts_file,
                ref_date=inps.refDate,
                outfile=inps.outfile,
                max_memory=inps.maxMemory,
                force=inps.force,
            )

            # pause to distinguish the modification time of input files
            time.sleep(1)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return
