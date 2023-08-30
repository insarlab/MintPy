############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os

import numpy as np

from mintpy import view
from mintpy.objects import HDFEOS
from mintpy.utils import ptime, readfile, utils as ut, writefile


##############################################################################
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
        print(f'change reference point to y/x: {inps.ref_yx}')

    # various file types
    print(f'read {inps.dset} from file {inps.file}')
    k = atr['FILE_TYPE']
    if k == 'velocity':
        # read/prepare data
        if not inps.dset:
            inps.dset = 'velocity'
            print('No selected dataset, assuming "velocity" and continue.')
        data, atr = readfile.read(inps.file, datasetName=inps.dset)

        # convert velocity to cumulative displacement
        if inps.dset == 'velocity':
            print('convert velocity to displacement for {}'.format(atr['DATE12']))
            dt1, dt2 = ptime.date_list2vector(atr['DATE12'].split('_'))[0]
            data *= (dt2 - dt1).days / 365.25

        # convert data from the unit of meter to radian
        if atr.get('UNIT', 'm/year').startswith('m'):
            print('convert the unit from meter to radian')
            data *= range2phase
            atr['UNIT'] = 'radian'

        # apply the custom spatial referencing
        if inps.ref_yx:
            data -= data[inps.ref_yx[0], inps.ref_yx[1]]

        # metadata
        atr['FILE_TYPE'] = '.unw'

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
        atr['DATE12'] = f'{date1[2:8]}-{date2[2:8]}'
        atr['FILE_TYPE'] = '.unw'
        atr['UNIT'] = 'radian'

        # output filename
        if not inps.outfile:
            inps.outfile = f'{date1}_{date2}.unw'
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
                raise ValueError(f"No '-' in input dataset! It is required for {dname}")
        else:
            date_list = HDFEOS(inps.file).get_date_list()
            date1 = date_list[0]
            date2 = date_list[-1]
        date12 = f'{date1}_{date2}'

        # read / prepare data
        slice_list = readfile.get_slice_list(inps.file)
        if 'displacement' in inps.dset:
            # read/prepare data
            slice_name1 = view.search_dataset_input(slice_list, f'{dname}-{date1}')[0][0]
            slice_name2 = view.search_dataset_input(slice_list, f'{dname}-{date2}')[0][0]
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
        atr['DATE12'] = f'{date1[2:8]}-{date2[2:8]}'
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
            raise ValueError(f'unrecognized input dataset type: {inps.dset}')

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
        atr['DATE12'] = f'{date1[2:8]}-{date2[2:8]}'
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
            raise ValueError(f'unrecognized dataset type: {inps.dset}')

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
            print(f'mask data based on input file: {m_file}')
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
def save_roipac(inps):

    # read data and metadata
    data, atr, out_file = read_data(inps)

    # remove non-roipac metadata
    if not inps.keepAllMetadata:
        atr = clean_metadata4roipac(atr)

    # write
    writefile.write(data, out_file=out_file, metadata=atr)

    return
