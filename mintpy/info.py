#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import h5py
import numpy as np

from mintpy.utils import readfile, ptime
from mintpy.utils.arg_utils import create_argument_parser
from mintpy.objects import (
    giantIfgramStack,
    giantTimeseries,
    ifgramStack,
    timeseries,
    HDFEOS,
)


############################################################
EXAMPLE = """example:
  info.py timeseries.h5
  info.py velocity.h5
  info.py ifgramStack.h5

  # Display dataset
  info.py timeseries.py --dset date
  info.py timeseries.py --dset bperp

  # Time / Date Info
  info.py ifgramStack.h5 --date                 #print date1_date2 info for all  interferograms
  info.py timeseries.h5  --num                  #print date list of timeseries with its number
  info.py ifgramStack.h5 --date --show kept     #print date1_date2 info for kept interferograms
  info.py ifgramStack.h5 --date --show dropped  #print date1_date2 info for dropped/excluded interferograms
  info.py LS-PARAMS.h5   --date > date_list.txt #print date list of timeseries and save it to txt file.
  info.py S1_IW12_128_0593_0597_20141213_20180619.h5 --date

  # save date1_date2 info of interferograms to a text file
  info.py ifgramStack.h5 --date --show kept > date12_list.txt  

  # Slice / Dataset Info
  info.py timeseries.h5                              --slice
  info.py inputs/ifgramStack.h5                      --slice
  info.py S1_IW12_128_0593_0597_20141213_20180619.h5 --slice
  info.py LS-PARAMS.h5                               --slice
"""


def create_parser(subparsers=None):
    """Create command line parser."""
    synopsis = 'Display Metadata / Structure information of ANY File'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='File to check')
    parser.add_argument('--compact', action='store_true',
                        help='show compact info by displaying only the top 20 metadata')

    parser.add_argument('--dset', type=str, help='Show dataset')

    par_list = parser.add_argument_group('List','list date/slice info')
    par_list.add_argument('--date', dest='disp_date', action='store_true',
                          help='Show date/date12 info of input file')
    par_list.add_argument('--num', dest='disp_num', action='store_true',
                          help='Show date/date12 info with numbers')
    par_list.add_argument('--slice', dest='disp_slice', action='store_true',
                          help='Show slice list of the file')
    par_list.add_argument('--show','--show-ifgram', dest='disp_ifgram',
                          choices={'all','kept','dropped'}, default='all',
                          help='Show all / kept / dropped interferograms only. Default: all.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if not os.path.isfile(inps.file):
        raise FileNotFoundError(inps.file)

    inps.max_meta_num = 200
    if inps.compact:
        inps.max_meta_num = 20
    return inps


############################################################
def attributes2string(atr, sorting=True, max_meta_num=200):
    ## Get Dictionary of Attributes
    digits = max([len(key) for key in list(atr.keys())] + [0])
    atr_string = ''
    i = 0
    for key, value in sorted(atr.items(), key=lambda x: x[0]):
        i += 1
        if i > max_meta_num:
            atr_string += '  ...\n'
            break
        else:
            # format metadata key/value
            try:
                value = value.decode('utf8')
            except:
                pass
            atr_string += '  {k:<{d}}    {v}\n'.format(k=key,
                                                       d=digits,
                                                       v=value)
    return atr_string


def print_attributes(atr, max_meta_num=200):
    atr_string = attributes2string(atr, max_meta_num=max_meta_num)
    print(atr_string)


def print_hdf5_structure(fname, max_meta_num=200):
    # generate string
    global h5_string, maxDigit
    h5_string = ''

    def hdf5_structure2string(name, obj):
        global h5_string, maxDigit
        if isinstance(obj, h5py.Group):
            h5_string += 'HDF5 group   "/{n}"\n'.format(n=name)
        elif isinstance(obj, h5py.Dataset):
            h5_string += ('HDF5 dataset "/{n:<{w}}": shape {s:<20}, '
                          'dtype <{t}>\n').format(n=name,
                                                  w=maxDigit,
                                                  s=str(obj.shape),
                                                  t=obj.dtype)
        atr = dict(obj.attrs)
        if len(atr) > 0:
            h5_string += attributes2string(atr, max_meta_num=max_meta_num)+"\n"

    f = h5py.File(fname, 'r')
    # grab metadata in root level as it will be missed in hdf5_structure2string()
    atr = dict(f.attrs)
    if len(atr) > 0:
        h5_string += 'Attributes in / level:\n'
        h5_string += attributes2string(atr, max_meta_num=max_meta_num)+'\n'

    # get maxDigit value 
    maxDigit = max([len(i) for i in f.keys()])
    maxDigit = max(20, maxDigit+1)
    if atr.get('FILE_TYPE', 'timeseries') == 'HDFEOS':
        maxDigit += 35

    # get structure string
    f.visititems(hdf5_structure2string)
    f.close()

    # print string
    print(h5_string)
    return h5_string


############################################################
def print_timseries_date_stat(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print('Start Date: {}'.format(dateList[0]))
    print('End   Date: {}'.format(dateList[-1]))
    print('Number of dates  : {}'.format(len(dateList)))
    print('STD of datetimes : {:.2f} years'.format(np.std(datevector)))
    #print('----------------------')
    #print('List of dates:\n{}'.format(dateList))
    #print('----------------------')
    #print('List of dates in years:\n{}'.format(datevector))
    return


def print_date_list(fname, disp_ifgram='all', disp_num=False, print_msg=False):
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
        print('--date option can not be applied to {} file, ignore it.'.format(k))

    # print list info
    if print_msg and dateList is not None:
        for d in dateList:
            if disp_num:
                if k in ['ifgramStack']:
                    num = dateListAll.index(d)
                else:
                    num = dateList.index(d)
                msg = '{}\t{}'.format(d, num)
            else:
                msg = d
            print(msg)
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
        msg = 'input dataset {} not found!'.format(dsName)
        msg += '\navailable datasets: {}'.format(dsNames)
        raise ValueError(msg)

    # print dataset values
    with h5py.File(fname, 'r') as f:
        data = f[dsName][:]
        print(data)

    # data stats
    print('dataset size: {}'.format(data.shape))
    print('dataset min / max: {} / {}'.format(np.nanmin(data), np.nanmax(data)))
    print('number of pixels in NaN: {}'.format(np.sum(np.isnan(data))))
    return


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    ext = os.path.splitext(inps.file)[1].lower()

    # --date/--num option
    if inps.disp_date or inps.disp_num:
        print_date_list(inps.file,
                        disp_ifgram=inps.disp_ifgram,
                        disp_num=inps.disp_num,
                        print_msg=True)
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
    if ext in ['.h5', '.he5']:
        print('\n{} {:*<40}'.format('*'*20, 'HDF5 File Structure '))
        print_hdf5_structure(inps.file, max_meta_num=inps.max_meta_num)
    else:
        print('\n{} {:*<40}'.format('*'*20, 'Binary File Attributes '))
        atr = readfile.read_attribute(inps.file)
        print_attributes(atr, max_meta_num=inps.max_meta_num)

    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
