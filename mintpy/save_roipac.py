#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import numpy as np
from mintpy.objects import timeseries, HDFEOS
from mintpy.utils import readfile, writefile, ptime, utils as ut
from mintpy.utils.arg_utils import create_argument_parser
from mintpy import view


##############################################################################
EXAMPLE = """example:
  #----- unwrapped phase
  #for velocity: output an interferogram with temporal baseline in DATE12 metadata
  save_roipac.py  velocity.h5
  save_roipac.py  velocity.h5 -m maskTempCoh.h5 maskAoiShinmoe.h5

  #for time-series: specify (date1_)date2
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  #use the last date
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  20050601
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  20040728_20050601

  #for HDF-EOS5: specify displacement-date1_date2
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5  displacement-20170904_20170916
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5  displacement-20170916

  #for ifgramStack: specify date1_date2
  save_roipac.py  inputs/ifgramStack.h5  unwrapPhase-20091225_20100723
  save_roipac.py  inputs/ifgramStack.h5  unwrapPhase-20091225_20100723  --ref-yx 640 810

  #----- coherence
  save_roipac.py  inputs/ifgramStack.h5  coherence-20091225_20100723
  save_roipac.py  temporalCoherence.h5
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5 temporalCoherence -o 20170904_20170916.cor

  #----- DEM
  save_roipac.py  geo_geometryRadar.h5  height -o srtm1.dem
  save_roipac.py  geo_geometryRadar.h5  height -o srtm1.hgt
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5 height -o srtm1.dem
"""


def create_parser(subparsers=None):
    synopsis = 'Convert MintPy HDF5 file to ROI_PAC format.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='HDF5 file to be converted.')
    parser.add_argument('dset', nargs='?', help='date/date12 of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('-m','--mask', dest='mask_file', nargs='+', help='mask file')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='custom reference pixel in y/x')
    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2, help='custom reference pixel in lat/lon')
    parser.add_argument('--keep-all-metadata', dest='keepAllMetadata', action='store_true', help='Do not clean the metadata as ROIPAC format')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default dset
    if not inps.dset:
        atr = readfile.read_attribute(inps.file)
        k = atr['FILE_TYPE']
        if k in ['ifgramStack', 'HDFEOS']:
            raise Exception("NO input dataset! It's required for {} file".format(k))

        #for time-series
        if k == 'timeseries':
            inps.dset = timeseries(inps.file).get_date_list()[-1]
            print('NO date specified >>> continue with the last date: {}'.format(inps.dset))
    return inps


def read_data(inps):
    # metadata
    atr = readfile.read_attribute(inps.file)

    if 'WAVELENGTH' in atr.keys():
        range2phase = -4 * np.pi / float(atr['WAVELENGTH'])

    # change reference pixel
    if inps.ref_lalo:
        if 'Y_FIRST' in atr.keys():
            coord = ut.coordinate(atr)
            ref_y, ref_x = coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1])[0:2]
            inps.ref_yx = [ref_y, ref_x]
        else:
            raise ValueError("input file is not geocoded --> reference point in lat/lon is NOT support")

    if inps.ref_yx:
        atr['REF_Y'] = inps.ref_yx[0]
        atr['REF_X'] = inps.ref_yx[1]
        if 'Y_FIRST' in atr.keys():
            coord = ut.coordinate(atr)
            ref_lat, ref_lon = coord.radar2geo(inps.ref_yx[0], inps.ref_yx[1])[0:2]
            atr['REF_LAT'] = ref_lat
            atr['REF_LON'] = ref_lon
        print('change reference point to y/x: {}'.format(inps.ref_yx))

    # various file types
    print('read {} from file {}'.format(inps.dset, inps.file))
    k = atr['FILE_TYPE']
    if k == 'velocity':
        # read/prepare data
        data = readfile.read(inps.file)[0]

        # velocity to displacement
        print('convert velocity to displacement for {}'.format(atr['DATE12']))
        date1, date2 = atr['DATE12'].split('_')
        dt1, dt2 = ptime.date_list2vector([date1, date2])[0]
        data *= (dt2 - dt1).days / 365.25

        # displacement to phase
        print('convert displacement to phase in radian')
        data *= range2phase

        if inps.ref_yx:
            data -= data[inps.ref_yx[0], inps.ref_yx[1]]

        # metadata
        atr['FILE_TYPE'] = '.unw'
        atr['UNIT'] = 'radian'

        # output filename
        if not inps.outfile:
            inps.outfile = os.path.join(os.path.dirname(inps.file), '{}.unw'.format(atr['DATE12']))

    elif k == 'timeseries':
        # date1 and date2
        if '_' in inps.dset:
            date1, date2 = ptime.yyyymmdd(inps.dset.split('_'))
        else:
            date1 = atr['REF_DATE']
            date2 = ptime.yyyymmdd(inps.dset)

        # read/prepare data
        data = readfile.read(inps.file, datasetName=date2)[0]
        data -= readfile.read(inps.file, datasetName=date1)[0]
        print('converting range to phase')
        data *= range2phase
        if inps.ref_yx:
            data -= data[inps.ref_yx[0], inps.ref_yx[1]]

        # metadata
        atr['DATE'] = date1[2:8]
        atr['DATE12'] = '{}-{}'.format(date1[2:8], date2[2:8])
        atr['FILE_TYPE'] = '.unw'
        atr['UNIT'] = 'radian'

        # output filename
        if not inps.outfile:
            inps.outfile = '{}_{}.unw'.format(date1, date2)
            if inps.file.startswith('geo_'):
                inps.outfile = 'geo_'+inps.outfile

    elif k == 'HDFEOS':
        dname = inps.dset.split('-')[0]

        # date1 and date2
        if dname == 'displacement':
            if '-' in inps.dset:
                suffix = inps.dset.split('-')[1]
                if '_' in suffix:
                    date1, date2 = ptime.yyyymmdd(suffix.split('_'))
                else:
                    date1 = atr['REF_DATE']
                    date2 = ptime.yyyymmdd(suffix)
            else:
                raise ValueError("No '-' in input dataset! It is required for {}".format(dname))
        else:
            date_list = HDFEOS(inps.file).get_date_list()
            date1 = date_list[0]
            date2 = date_list[-1]
        date12 = '{}_{}'.format(date1, date2)

        # read / prepare data
        slice_list = readfile.get_slice_list(inps.file)
        if 'displacement' in inps.dset:
            # read/prepare data
            slice_name1 = view.search_dataset_input(slice_list, '{}-{}'.format(dname, date1))[0][0]
            slice_name2 = view.search_dataset_input(slice_list, '{}-{}'.format(dname, date2))[0][0]
            data = readfile.read(inps.file, datasetName=slice_name1)[0]
            data -= readfile.read(inps.file, datasetName=slice_name2)[0]
            print('converting range to phase')
            data *= range2phase
            if inps.ref_yx:
                data -= data[inps.ref_yx[0], inps.ref_yx[1]]
        else:
            slice_name = view.search_dataset_input(slice_list, inps.dset)[0][0]
            data = readfile.read(inps.file, datasetName=slice_name)[0]

        # metadata
        atr['DATE'] = date1[2:8]
        atr['DATE12'] = '{}-{}'.format(date1[2:8], date2[2:8])
        if dname == 'displacement':
            atr['FILE_TYPE'] = '.unw'
            atr['UNIT'] = 'radian'
        elif 'coherence' in dname.lower():
            atr['FILE_TYPE'] = '.cor'
            atr['UNIT'] = '1'
        elif dname == 'height':
            atr['FILE_TYPE'] = '.dem'
            atr['DATA_TYPE'] = 'int16'
        else:
            raise ValueError('unrecognized input dataset type: {}'.format(inps.dset))

        # output filename
        if not inps.outfile:
            inps.outfile = '{}{}'.format(date12, atr['FILE_TYPE'])

    elif k == 'ifgramStack':
        dname, date12 = inps.dset.split('-')
        date1, date2 = date12.split('_')

        # read / prepare data
        data = readfile.read(inps.file, datasetName=inps.dset)[0]
        if dname.startswith('unwrapPhase'):
            if 'REF_X' in atr.keys():
                data -= data[int(atr['REF_Y']), int(atr['REF_X'])]
                print('consider reference pixel in y/x: ({}, {})'.format(atr['REF_Y'], atr['REF_X']))
            else:
                print('No REF_Y/X found.')

        # metadata
        atr['DATE'] = date1[2:8]
        atr['DATE12'] = '{}-{}'.format(date1[2:8], date2[2:8])
        if dname.startswith('unwrapPhase'):
            atr['FILE_TYPE'] = '.unw'
            atr['UNIT'] = 'radian'
        elif dname == 'coherence':
            atr['FILE_TYPE'] = '.cor'
            atr['UNIT'] = '1'
        elif dname == 'wrapPhase':
            atr['FILE_TYPE'] = '.int'
            atr['UNIT'] = 'radian'
        elif dname == 'connectComponent':
            atr['FILE_TYPE'] = '.conncomp'
            atr['UNIT'] = '1'
            atr['DATA_TYPE'] = 'byte'
        else:
            raise ValueError('unrecognized dataset type: {}'.format(inps.dset))

        # output filename
        if not inps.outfile:
            inps.outfile = '{}{}'.format(date12, atr['FILE_TYPE'])
            if inps.file.startswith('geo_'):
                inps.outfile = 'geo_'+inps.outfile

    else:
        # read data
        data = readfile.read(inps.file, datasetName=inps.dset)[0]

        if inps.outfile:
            # get file extension
            fbase, fext = os.path.splitext(inps.outfile)
            # ignore certain meaningless file extensions
            while fext in ['.geo', '.rdr', '.full', '.wgs84']:
                fbase, fext = os.path.splitext(fbase)
            if not fext:
                fext = os.path.basename(fbase)

            atr['FILE_TYPE'] = fext
        else:
            # metadata
            if 'coherence' in k.lower():
                atr['FILE_TYPE'] = '.cor'
            elif k in ['mask']:
                atr['FILE_TYPE'] = '.msk'    
            elif k in ['geometry'] and inps.dset == 'height':
                if 'Y_FIRST' in atr.keys():
                    atr['FILE_TYPE'] = '.dem'
                else:
                    atr['FILE_TYPE'] = '.hgt'
                atr['UNIT'] = 'm'
            else:
                atr['FILE_TYPE'] = '.unw'

            inps.outfile = '{}{}'.format(os.path.splitext(inps.file)[0], atr['FILE_TYPE'])

    # mask
    if inps.mask_file:
        for m_file in inps.mask_file:
            print('mask data based on input file: {}'.format(m_file))
            mask = readfile.read(m_file)[0]
            mask *= ~np.isnan(data)
            data[mask==0] = np.nan

    # get rid of starting . if output as hdf5 file
    if inps.outfile.endswith('.h5'):
        if atr['FILE_TYPE'].startswith('.'):
            atr['FILE_TYPE'] = atr['FILE_TYPE'][1:]

    atr['PROCESSOR'] = 'roipac'
    return data, atr, inps.outfile


def clean_metadata4roipac(atr_in):
    atr = {}
    for key, value in atr_in.items():
        atr[key] = str(value)

    # drop the following keys
    key_list = ['width', 'Width', 'samples', 'length', 'lines',
                'SUBSET_XMIN','SUBSET_XMAX','SUBSET_YMIN','SUBSET_YMAX',
               ]
    for key in key_list:
        if key in atr.keys():
            atr.pop(key)

    # drop all keys that are not all UPPER_CASE
    key_list = list(atr.keys())
    for key in key_list:
        if not key.isupper():
            atr.pop(key)

    atr['FILE_LENGTH'] = atr['LENGTH']
    return atr


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    data, atr, out_file = read_data(inps)

    if not inps.keepAllMetadata:
        atr = clean_metadata4roipac(atr)

    writefile.write(data, out_file=out_file, metadata=atr)

    return


##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
