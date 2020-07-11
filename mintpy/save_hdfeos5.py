#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import re
import argparse
import datetime as dt
import h5py
import numpy as np
from mintpy.objects import timeseries, geometry, sensor
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile
from mintpy import info


BOOL_ZERO = np.bool_(0)
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)
compression = 'lzf'


################################################################
TEMPALTE = TEMPLATE = get_template_content('hdfeos5')

EXAMPLE = """example:
  save_hdfeos5.py geo_timeseries_ERA5_ramp_demErr.h5 -c geo_temporalCoherence.h5 -m geo_maskTempCoh.h5 -g geo_geometryRadar.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Convert MintPy timeseries product into HDF-EOS5 format\n' +
                                     'https://earthdata.nasa.gov/esdis/eso/standards-and-references/hdf-eos5',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=TEMPALTE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', default='timeseries.h5', help='Timeseries file')
    parser.add_argument('-t', '--template', dest='template_file', help='Template file')

    parser.add_argument('-c', '--coherence', dest='coherence_file', required=True, 
                        help='Coherence/correlation file, i.e. avgSpatialCoh.h5, temporalCoherence.h5')
    parser.add_argument('-m', '--mask', dest='mask_file', required=True, help='Mask file')
    parser.add_argument('-g', '--geometry', dest='geom_file', required=True, help='geometry file')

    parser.add_argument('--update', action='store_true',
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')
    parser.add_argument('--subset', action='store_true',
                        help='Enable subset mode, a.k.a. put suffix _N31700_N32100_E130500_E131100')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################
def metadata_mintpy2unavco(meta_dict_in, dateList):
    # Extract UNAVCO format metadata from MintPy attributes dictionary and dateList
    meta_dict = {}
    for key in meta_dict_in.keys():
        meta_dict[key] = meta_dict_in[key]
        for prefix in ['unavco.', 'hdfeos5.']:
            if prefix in key.lower():
                key2 = key.lower().split(prefix)[1]
                meta_dict[key2] = meta_dict_in[key]

    unavco_meta_dict = dict()
    #################################
    # Required metadata
    #################################
    # Given manually
    # mission
    # ERS,ENV,S1,RS1,RS2,CSK,TSX,JERS,ALOS,ALOS2
    try:
        unavco_meta_dict['mission'] = sensor.get_unavco_mission_name(meta_dict)
    except ValueError:
        print('Missing required attribute: mission')

    # beam_mode/swath
    unavco_meta_dict['beam_mode'] = meta_dict['beam_mode']
    unavco_meta_dict['beam_swath'] = int(meta_dict.get('beam_swath', '0'))

    # relative_orbit, or track number
    #atr_dict['relative_orbit'] = int(re.match(r'(\w+)T([0-9+])',atr['PROJECT_NAME']).groups()[1])
    unavco_meta_dict['relative_orbit'] = int(meta_dict['relative_orbit'])

    # processing info
    unavco_meta_dict['processing_type'] = meta_dict.get('processing_type', 'LOS_TIMESERIES')
    #unavco_meta_dict['processing_software'] = meta_dict['processing_software']

    # Grabbed by script
    # date info
    unavco_meta_dict['first_date'] = dt.datetime.strptime(dateList[0], '%Y%m%d').isoformat()[0:10]
    unavco_meta_dict['last_date'] = dt.datetime.strptime(dateList[-1], '%Y%m%d').isoformat()[0:10]

    # footprint
    lons = [meta_dict['LON_REF1'],
            meta_dict['LON_REF3'],
            meta_dict['LON_REF4'],
            meta_dict['LON_REF2'],
            meta_dict['LON_REF1']]

    lats = [meta_dict['LAT_REF1'],
            meta_dict['LAT_REF3'],
            meta_dict['LAT_REF4'],
            meta_dict['LAT_REF2'],
            meta_dict['LAT_REF1']]

    unavco_meta_dict['scene_footprint'] = "POLYGON((" + ",".join(
        [lon+' '+lat for lon, lat in zip(lons, lats)]) + "))"

    unavco_meta_dict['history'] = dt.datetime.utcnow().isoformat()[0:10]

    #################################
    # Recommended metadata
    #################################
    # Given manually
    if 'frame' in meta_dict.keys():
        unavco_meta_dict['frame'] = int(meta_dict['frame'])
    elif 'first_frame' in meta_dict.keys():
        unavco_meta_dict['frame'] = int(meta_dict['first_frame'])
    else:
        unavco_meta_dict['frame'] = 0

    unavco_meta_dict['atmos_correct_method']   = meta_dict.get('atmos_correct_method', 'None')
    unavco_meta_dict['post_processing_method'] = meta_dict.get('post_processing_method', 'MintPy')
    unavco_meta_dict['processing_dem'] = meta_dict.get('processing_dem', 'Unknown')
    unavco_meta_dict['unwrap_method']  = meta_dict.get('unwrap_method', 'Unknown')

    # Grabbed by script
    unavco_meta_dict['flight_direction'] = meta_dict.get('ORBIT_DIRECTION', 'Unknown')[0].upper()

    if meta_dict['ANTENNA_SIDE'] == '-1':
        unavco_meta_dict['look_direction'] = 'R'
    else:
        unavco_meta_dict['look_direction'] = 'L'

    unavco_meta_dict['polarization'] = meta_dict.get('POLARIZATION', 'Unknown')
    unavco_meta_dict['prf'] = float(meta_dict.get('PRF', '0'))
    unavco_meta_dict['wavelength'] = float(meta_dict['WAVELENGTH'])

    #################################
    # insarmaps metadata
    #################################
    # footprint for data coverage
    if 'X_FIRST' in meta_dict.keys():
        lon0 = float(meta_dict['X_FIRST'])
        lat0 = float(meta_dict['Y_FIRST'])
        lon1 = lon0 + float(meta_dict['X_STEP'])*int(meta_dict['WIDTH'])
        lat1 = lat0 + float(meta_dict['Y_STEP'])*int(meta_dict['LENGTH'])
        lons = [str(lon0), str(lon1), str(lon1), str(lon0), str(lon0)]
        lats = [str(lat0), str(lat0), str(lat1), str(lat1), str(lat0)]
        unavco_meta_dict['data_footprint'] = "POLYGON((" + ",".join(
            [lon+' '+lat for lon, lat in zip(lons, lats)]) + "))"
    else:
        print('Input file is not geocoded, no data_footprint without X/Y_FIRST/STEP info.')

    return unavco_meta_dict


def prep_metadata(ts_file, print_msg=True):
    """Prepare metadata for HDF-EOS5 file"""
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    unavco_meta_dict = metadata_mintpy2unavco(ts_obj.metadata, ts_obj.dateList)
    if print_msg:
        print('## UNAVCO Metadata:')
        print('-----------------------------------------')
        info.print_attributes(unavco_meta_dict)
        print('-----------------------------------------')

    meta_dict = dict(ts_obj.metadata)
    meta_dict.update(unavco_meta_dict)
    meta_dict['FILE_TYPE'] = 'HDFEOS'
    return meta_dict


def get_output_filename(metadata, update_mode=False, subset_mode=False):
    """Get output file name of HDF-EOS5 time-series file"""
    SAT = metadata['mission']
    SW = metadata['beam_mode']
    if metadata['beam_swath']:
        SW += str(metadata['beam_swath'])
    RELORB = "%03d" % (int(metadata['relative_orbit']))

    # Frist and/or Last Frame
    frame1 = int(metadata['frame'])
    key = 'first_frame'
    if key in metadata.keys():
        frame1 = int(metadata[key])
    FRAME = "%04d" % (frame1)
    key = 'last_frame'
    if key in metadata.keys():
        frame2 = int(metadata[key])
        if frame2 != frame1:
            FRAME += "_%04d" % (frame2)

    DATE1 = dt.datetime.strptime(metadata['first_date'], '%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.datetime.strptime(metadata['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    if update_mode:
        print('Update mode is ON, put endDate as XXXXXXXX.')
        DATE2 = 'XXXXXXXX'

    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'_'+DATE2+'.he5'

    if subset_mode:
        print('Subset mode is enabled, put subset range info in output filename.')
        lat1 = float(metadata['Y_FIRST'])
        lon0 = float(metadata['X_FIRST'])
        lat0 = lat1 + float(metadata['Y_STEP']) * int(metadata['LENGTH'])
        lon1 = lon0 + float(metadata['X_STEP']) * int(metadata['WIDTH'])

        lat0Str = 'N%05d' % (round(lat0*1e3))
        lat1Str = 'N%05d' % (round(lat1*1e3))
        lon0Str = 'E%06d' % (round(lon0*1e3))
        lon1Str = 'E%06d' % (round(lon1*1e3))
        if lat0 < 0.0: lat0Str = 'S%05d' % (round(abs(lat0)*1e3))
        if lat1 < 0.0: lat1Str = 'S%05d' % (round(abs(lat1)*1e3))
        if lon0 < 0.0: lon0Str = 'W%06d' % (round(abs(lon0)*1e3))
        if lon1 < 0.0: lon1Str = 'W%06d' % (round(abs(lon1)*1e3))

        SUB = '_%s_%s_%s_%s' % (lat0Str, lat1Str, lon0Str, lon1Str)
        outName = '{}{}{}'.format(os.path.splitext(outName)[0],
                                  SUB,
                                  os.path.splitext(outName)[1])
    return outName


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()

    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(template_file)

    # Coherence-based network modification
    prefix = 'mintpy.save.hdfEos5.'

    key = prefix+'update'
    if key in template.keys() and template[key] == 'yes':
        inps.update = True

    key = prefix+'subset'
    if key in template.keys() and template[key] == 'yes':
        inps.subset = True

    return inps


def write2hdf5(out_file, ts_file, coh_file, mask_file, geom_file, metadata):
    """Write HDF5 file in HDF-EOS5 format"""
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    dateList = ts_obj.dateList

    # Open HDF5 File
    f = h5py.File(out_file, 'w')
    print('create HDF5 file: {} with w mode'.format(out_file))
    maxDigit = 55

    # Write Observation - Displacement
    gName = 'HDFEOS/GRIDS/timeseries/observation'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    dsName = 'displacement'
    data = ts_obj.read(print_msg=False)
    print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
           ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                           w=maxDigit,
                                           t=str(data.dtype),
                                           s=data.shape,
                                           c=compression))
    dset = group.create_dataset(dsName,
                                data=data,
                                dtype=np.float32,
                                chunks=True,
                                compression=compression)
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'meters'

    dsName = 'date'
    data = np.array(dateList, dtype=np.string_)
    group.create_dataset(dsName, data=data)
    print('create dataset /{d:<{w}} of {t:<10} in size of {s}'.format(d='{}/{}'.format(gName, dsName),
                                                                      w=maxDigit,
                                                                      t=str(data.dtype),
                                                                      s=data.shape))

    dsName = 'bperp'
    data = np.array(ts_obj.pbase, dtype=np.float32)
    group.create_dataset(dsName, data=data)
    print('create dataset /{d:<{w}} of {t:<10} in size of {s}'.format(d='{}/{}'.format(gName, dsName),
                                                                      w=maxDigit,
                                                                      t=str(data.dtype),
                                                                      s=data.shape))

    # Write Quality
    gName = 'HDFEOS/GRIDS/timeseries/quality'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    ## 1 - temporalCoherence
    dsName = 'temporalCoherence'
    data = readfile.read(coh_file)[0]
    print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
           ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                           w=maxDigit,
                                           t=str(data.dtype),
                                           s=data.shape,
                                           c=compression))
    dset = group.create_dataset(dsName,
                                data=data,
                                chunks=True,
                                compression=compression)
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'

    ## 2 - mask
    dsName = 'mask'
    data = readfile.read(mask_file, datasetName='mask')[0]
    print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
           ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                           w=maxDigit,
                                           t=str(data.dtype),
                                           s=data.shape,
                                           c=compression))
    dset = group.create_dataset(dsName,
                                data=data,
                                chunks=True,
                                compression=compression)
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = BOOL_ZERO
    dset.attrs['_FillValue'] = BOOL_ZERO
    dset.attrs['Units'] = '1'

    # Write Geometry
    # Required: height, incidenceAngle
    # Optional: rangeCoord, azimuthCoord, azimuthAngle, slantRangeDistance, waterMask, shadowMask
    gName = 'HDFEOS/GRIDS/timeseries/geometry'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    geom_obj = geometry(geom_file)
    geom_obj.open(print_msg=False)
    for dsName in geom_obj.datasetNames:
        data = geom_obj.read(datasetName=dsName, print_msg=False)
        print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
               ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                               w=maxDigit,
                                               t=str(data.dtype),
                                               s=data.shape,
                                               c=compression))
        dset = group.create_dataset(dsName,
                                    data=data,
                                    chunks=True,
                                    compression=compression)

        dset.attrs['Title'] = dsName
        if dsName in ['height',
                      'slantRangeDistance',
                      'bperp']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = 'meters'

        elif dsName in ['incidenceAngle',
                        'azimuthAngle',
                        'latitude',
                        'longitude']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = 'degrees'

        elif dsName in ['rangeCoord', 'azimuthCoord']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = '1'

        elif dsName in ['waterMask', 'shadowMask']:
            dset.attrs['MissingValue'] = BOOL_ZERO
            dset.attrs['_FillValue'] = BOOL_ZERO
            dset.attrs['Units'] = '1'

    # Write Attributes to the HDF File
    print('write metadata to root level')
    for key, value in iter(metadata.items()):
        f.attrs[key] = value
    f.close()
    print('finished writing to {}'.format(out_file))
    return out_file


################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # Prepare Metadata
    meta_dict = prep_metadata(ts_file=inps.timeseries_file, print_msg=True)

    # Get output filename
    outName = get_output_filename(metadata=meta_dict,
                                  update_mode=inps.update,
                                  subset_mode=inps.subset)

    # Open HDF5 File
    write2hdf5(out_file=outName,
               ts_file=inps.timeseries_file,
               coh_file=inps.coherence_file,
               mask_file=inps.mask_file,
               geom_file=inps.geom_file,
               metadata=meta_dict)
    return outName


################################################################
if __name__ == '__main__':
    main()
