############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os

import h5py
import numpy as np

from mintpy.objects import (
    HDFEOS,
    giantIfgramStack,
    giantTimeseries,
    ifgramStack,
    timeseries,
)
from mintpy.utils import ptime, readfile


############################################################
def attributes2string(atr, sorting=True, max_meta_num=200):
    ## Get Dictionary of Attributes
    digits = max([len(key) for key in list(atr.keys())] + [0])
    atr_str = ''
    i = 0
    for key, value in sorted(atr.items(), key=lambda x: x[0]):
        i += 1
        if i > max_meta_num:
            atr_str += '  ...\n'
            break
        else:
            # format metadata key/value
            try:
                value = value.decode('utf8')
            except:
                pass
            atr_str += f'  {key:<{digits}}    {value}\n'

    return atr_str


def print_attributes(atr, max_meta_num=200):
    atr_str = attributes2string(atr, max_meta_num=max_meta_num)
    print(atr_str)


def print_hdf5_structure(fname, max_meta_num=200):
    # generate string
    global h5_str, maxDigit
    h5_str = ''

    def hdf5_structure2string(name, obj):
        global h5_str, maxDigit
        if isinstance(obj, h5py.Group):
            h5_str += f'HDF5 group   "/{name}"\n'
        elif isinstance(obj, h5py.Dataset):
            h5_str += f'HDF5 dataset "/{name:<{maxDigit}}": shape={str(obj.shape):<20}, '
            h5_str += f'dtype={str(obj.dtype):<10}, compression={obj.compression}\n'

        atr = dict(obj.attrs)
        if len(atr) > 0:
            h5_str += attributes2string(atr, max_meta_num=max_meta_num)

    with h5py.File(fname, 'r') as f:

        # grab metadata in root level as it will be missed in hdf5_structure2string()
        atr = dict(f.attrs)
        if len(atr) > 0:
            h5_str += 'Attributes in / level:\n'
            h5_str += attributes2string(atr, max_meta_num=max_meta_num)+'\n'

        # get maxDigit value
        maxDigit = max(len(i) for i in f.keys())
        maxDigit = max(20, maxDigit+1)
        if atr.get('FILE_TYPE', 'timeseries') == 'HDFEOS':
            maxDigit += 35

        # get structure string
        f.visititems(hdf5_structure2string)

    # print string
    print(h5_str)

    return h5_str


############################################################
def print_timseries_date_stat(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print(f'Start Date: {dateList[0]}')
    print(f'End   Date: {dateList[-1]}')
    print(f'Number of dates  : {len(dateList)}')
    print(f'STD of datetimes : {np.std(datevector):.2f} years')
    #print('----------------------')
    #print('List of dates:\n{}'.format(dateList))
    #print('----------------------')
    #print('List of dates in years:\n{}'.format(datevector))
    return


def print_date_list(fname, disp_ifgram='all', disp_num=False, max_num=1e4, print_msg=False):
    """Print time/date info of file"""
    k = readfile.read_attribute(fname)['FILE_TYPE']
    dateList = None
    if k in ['timeseries']:
        dateList = timeseries(fname).get_date_list()

    elif k == 'HDFEOS':
        dateList = HDFEOS(fname).get_date_list()

    elif k == 'giantTimeseries':
        dateList = giantTimeseries(fname).get_date_list()


    elif k in ['giantIfgramStack']:
        dateList = giantIfgramStack(fname).get_date12_list()

    elif k in ['ifgramStack']:
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        dateListAll = obj.get_date12_list(dropIfgram=False)
        dateListKept = obj.get_date12_list(dropIfgram=True)

        # show dropped ifgram or not
        if disp_ifgram == 'all':
            dateList = list(dateListAll)
        elif disp_ifgram == 'kept':
            dateList = list(dateListKept)
        else:
            dateList = sorted(list(set(dateListAll) - set(dateListKept)))

    else:
        print(f'--date option can not be applied to {k} file, ignore it.')

    # print list info
    max_num = int(max_num)
    if print_msg and dateList is not None:
        for d in dateList[:max_num]:
            if disp_num:
                if k in ['ifgramStack']:
                    num = dateListAll.index(d)
                else:
                    num = dateList.index(d)
                msg = f'{d}\t{num}'
            else:
                msg = d
            print(msg)

        # add ... at the end if --compact
        if max_num < len(dateList):
            print('...\n')

    return dateList


def print_slice_list(fname, print_msg=False):
    """Print slice info of file"""
    slice_list = readfile.get_slice_list(fname)
    if print_msg:
        for slice_name in slice_list:
            print(slice_name)
    return slice_list


def print_aux_info(fname):
    try:
        atr = readfile.read_attribute(fname)
        k = atr['FILE_TYPE']
        print('{} {:*<40}'.format('*'*20, 'Basic File Info '))
        print('file name: '+atr['FILE_PATH'])
        print('file type: '+atr['FILE_TYPE'])
        if 'Y_FIRST' in atr.keys():
            print('coordinates : GEO')
        else:
            print('coordinates : RADAR')
        if k in ['timeseries']:
            dateList = print_date_list(fname)
            print('\n{} {:*<40}'.format('*'*20, 'Date Stat Info '))
            print_timseries_date_stat(dateList)
    except:
        pass
    return


def print_dataset(fname, dsName):
    # get available dataset list
    global dsNames

    def get_hdf5_dataset(name, obj):
        global dsNames
        if isinstance(obj, h5py.Dataset):
            dsNames.append(name)

    dsNames = []
    with h5py.File(fname, 'r') as f:
        f.visititems(get_hdf5_dataset)

    # check input dataset
    if dsName not in dsNames:
        msg = f'input dataset {dsName} not found!'
        msg += f'\navailable datasets: {dsNames}'
        raise ValueError(msg)

    # print dataset values
    with h5py.File(fname, 'r') as f:
        data = f[dsName][:]
        print(data)

    # data stats
    print(f'dataset size: {data.shape}')
    print(f'dataset min / max: {np.nanmin(data)} / {np.nanmax(data)}')
    print(f'number of pixels in NaN: {np.sum(np.isnan(data))}')

    return


def print_info(inps):
    """Extract and print the input file structure information."""

    # --date/--num option
    if inps.disp_date or inps.disp_num:
        print_date_list(
            inps.file,
            disp_ifgram=inps.disp_ifgram,
            disp_num=inps.disp_num,
            max_num=inps.max_meta_num,
            print_msg=True,
        )
        return

    # --slice option
    if inps.disp_slice:
        if inps.disp_ifgram != 'all':
            raise ValueError('--show-ifgram option is not applicable to --slice.')
        print_slice_list(inps.file, print_msg=True)
        return

    # --dset option
    if inps.dset:
        print_dataset(inps.file, dsName=inps.dset)
        return

    # Basic info
    print_aux_info(inps.file)

    # Generic Attribute/Structure of all files
    fext = os.path.splitext(inps.file)[1].lower()
    if fext in ['.h5', '.he5']:
        print('\n{} {:*<40}'.format('*'*20, 'HDF5 File Structure '))
        print_hdf5_structure(
            inps.file,
            max_meta_num=inps.max_meta_num,
        )

    else:
        print('\n{} {:*<40}'.format('*'*20, 'Binary File Attributes '))
        atr = readfile.read_attribute(inps.file)
        print_attributes(
            atr,
            max_meta_num=inps.max_meta_num,
        )

    return
