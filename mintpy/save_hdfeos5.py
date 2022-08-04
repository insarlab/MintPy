#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import sys
import datetime as dt
import h5py
import numpy as np

from mintpy.objects import timeseries, geometry, sensor
from mintpy.defaults.template import get_template_content
from mintpy.utils import ptime, readfile
from mintpy.utils.arg_utils import create_argument_parser
from mintpy import info


BOOL_ZERO = np.bool_(0)
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)
COMPRESSION = 'lzf'


################################################################
TEMPALTE = TEMPLATE = get_template_content('hdfeos5')

EXAMPLE = """example:
  save_hdfeos5.py geo/geo_timeseries_ERA5_ramp_demErr.h5
  save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5 --tc temporalCoherence.h5 --asc avgSpatialCoh.h5 -m maskTempCoh.h5 -g inputs/geometryGeo.h5
"""

NOTE = """
  https://earthdata.nasa.gov/esdis/eso/standards-and-references/hdf-eos5
  https://mintpy.readthedocs.io/en/latest/hdfeos5/
"""

def create_parser(subparsers=None):
    synopsis = 'Convert MintPy timeseries product into HDF-EOS5 format'
    epilog = TEMPALTE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file', default='timeseries.h5', help='Timeseries file')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='Template file for 1) arguments/options and 2) missing metadata')

    parser.add_argument('--tc','--temp-coh', dest='tcoh_file',
                        help='Coherence/correlation file, i.e. temporalCoherence.h5')
    parser.add_argument('--asc','--avg-spatial-coh', dest='scoh_file',
                        help='Average spatial coherence file, i.e. avgSpatialCoh.h5')
    parser.add_argument('-m', '--mask', dest='mask_file', help='Mask file')
    parser.add_argument('-g', '--geometry', dest='geom_file', help='geometry file')
    parser.add_argument('--suffix', dest='suffix', help='suffix to be appended to file name (e.g. PS).')

    parser.add_argument('--update', action='store_true',
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')
    parser.add_argument('--subset', action='store_true',
                        help='Enable subset mode, a.k.a. put suffix _N31700_N32100_E130500_E131100')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default filenames
    ts_dir = os.path.dirname(inps.ts_file)
    if os.path.basename(inps.ts_file).startswith('geo_'):
        tcoh_file = os.path.join(ts_dir, 'geo_temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'geo_avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'geo_maskTempCoh.h5')
        geom_file = os.path.join(ts_dir, 'geo_geometryRadar.h5')
    else:
        tcoh_file = os.path.join(ts_dir, 'temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'maskTempCoh.h5')
        geom_file = os.path.join(ts_dir, 'inputs/geometryGeo.h5')

    if not inps.tcoh_file:  inps.tcoh_file = tcoh_file
    if not inps.scoh_file:  inps.scoh_file = scoh_file
    if not inps.mask_file:  inps.mask_file = mask_file
    if not inps.geom_file:  inps.geom_file = geom_file

    # check file existence
    for fname in [inps.ts_file, inps.tcoh_file, inps.scoh_file, inps.mask_file, inps.geom_file]:
        if not os.path.isfile(fname):
            raise FileNotFoundError(fname)

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()

    if not template_file:
        return inps, None

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

    return inps, template


################################################################
def prep_metadata(ts_file, template=None, print_msg=True):
    """Prepare metadata for HDF-EOS5 file."""
    # read metadata from ts_file
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    meta = dict(ts_obj.metadata)

    # read metadata from template_file
    for key, value in template.items():
        if not key.startswith(('mintpy', 'isce')):
            meta[key] = value

    # grab unavco metadata
    unavco_meta = metadata_mintpy2unavco(meta, ts_obj.dateList)
    if print_msg:
        print('## UNAVCO Metadata:')
        print('-----------------------------------------')
        info.print_attributes(unavco_meta)
        print('-----------------------------------------')

    # update metadata from unavco metadata
    meta.update(unavco_meta)
    meta['FILE_TYPE'] = 'HDFEOS'

    return meta


def metadata_mintpy2unavco(meta_in, dateList):
    """Convert metadata from mintpy format into unavco format."""
    # Extract UNAVCO format metadata from MintPy attributes dictionary and dateList
    meta = {}
    for key in meta_in.keys():
        meta[key] = meta_in[key]
        for prefix in ['unavco.', 'hdfeos5.']:
            if prefix in key.lower():
                key2 = key.lower().split(prefix)[1]
                meta[key2] = meta_in[key]

    unavco_meta = dict()
    #################################
    # Required metadata
    #################################
    # Given manually
    # mission
    # ERS,ENV,S1,RS1,RS2,CSK,TSX,JERS,ALOS,ALOS2
    try:
        unavco_meta['mission'] = sensor.get_unavco_mission_name(meta)
    except ValueError:
        print('Missing required attribute: mission')

    # beam_mode/swath
    unavco_meta['beam_mode'] = meta['beam_mode']
    unavco_meta['beam_swath'] = int(meta.get('beam_swath', '0'))

    # relative_orbit, or track number
    unavco_meta['relative_orbit'] = int(meta['relative_orbit'])

    # processing info
    unavco_meta['processing_type'] = 'LOS_TIMESERIES'
    unavco_meta['processing_software'] = meta.get('PROCESSOR', 'isce')

    # Grabbed by script
    # date info
    unavco_meta['first_date'] = dt.datetime.strptime(dateList[0], '%Y%m%d').isoformat()[0:10]
    unavco_meta['last_date'] = dt.datetime.strptime(dateList[-1], '%Y%m%d').isoformat()[0:10]

    # footprint
    lons = [meta['LON_REF1'],
            meta['LON_REF3'],
            meta['LON_REF4'],
            meta['LON_REF2'],
            meta['LON_REF1']]

    lats = [meta['LAT_REF1'],
            meta['LAT_REF3'],
            meta['LAT_REF4'],
            meta['LAT_REF2'],
            meta['LAT_REF1']]

    unavco_meta['scene_footprint'] = "POLYGON((" + ",".join(
        [lon+' '+lat for lon, lat in zip(lons, lats)]) + "))"

    unavco_meta['history'] = dt.datetime.utcnow().isoformat()[0:10]

    #################################
    # Recommended metadata
    #################################
    unavco_meta['first_frame'] = int(meta.get('first_frame', 0))
    unavco_meta['last_frame'] = int(meta.get('last_frame', unavco_meta['first_frame']))

    unavco_meta['atmos_correct_method']   = meta.get('atmos_correct_method', 'None')
    unavco_meta['post_processing_method'] = 'MintPy'
    unavco_meta['processing_dem'] = meta.get('processing_dem', 'Unknown')
    unavco_meta['unwrap_method']  = meta.get('unwrap_method', 'Unknown')

    # Grabbed by script
    unavco_meta['flight_direction'] = meta.get('ORBIT_DIRECTION', 'Unknown')[0].upper()

    if meta['ANTENNA_SIDE'] == '-1':
        unavco_meta['look_direction'] = 'R'
    else:
        unavco_meta['look_direction'] = 'L'

    unavco_meta['polarization'] = meta.get('POLARIZATION', 'Unknown')
    unavco_meta['prf'] = float(meta.get('PRF', '0'))
    unavco_meta['wavelength'] = float(meta['WAVELENGTH'])

    #################################
    # insarmaps metadata
    #################################
    # footprint for data coverage
    if 'X_FIRST' in meta.keys():
        lon0 = float(meta['X_FIRST'])
        lat0 = float(meta['Y_FIRST'])
        lon1 = lon0 + float(meta['X_STEP'])*int(meta['WIDTH'])
        lat1 = lat0 + float(meta['Y_STEP'])*int(meta['LENGTH'])
        lons = [str(lon0), str(lon1), str(lon1), str(lon0), str(lon0)]
        lats = [str(lat0), str(lat0), str(lat1), str(lat1), str(lat0)]
        unavco_meta['data_footprint'] = "POLYGON((" + ",".join(
            [lon+' '+lat for lon, lat in zip(lons, lats)]) + "))"
    else:
        print('Input file is not geocoded, no data_footprint without X/Y_FIRST/STEP info.')

    return unavco_meta


def get_output_filename(metadata, suffix=None, update_mode=False, subset_mode=False):
    """Get output file name of HDF-EOS5 time-series file."""
    SAT = metadata['mission']
    SW = metadata['beam_mode']
    if metadata['beam_swath']:
        SW += str(metadata['beam_swath'])
    RELORB = "%03d" % (int(metadata['relative_orbit']))

    # Frist and/or Last Frame
    frame1 = metadata['first_frame']
    frame2 = metadata['last_frame']
    FRAME = "%04d" % (int(frame1))
    if frame2 != frame1:
        FRAME += "_%04d" % (frame2)

    DATE1 = dt.datetime.strptime(metadata['first_date'], '%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.datetime.strptime(metadata['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    if update_mode:
        print('Update mode is ON, put endDate as XXXXXXXX.')
        DATE2 = 'XXXXXXXX'

    if not suffix:
       outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'_'+DATE2+'.he5'
    else:
       outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'_'+DATE2+'_'+suffix+'.he5'

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


def create_hdf5_dataset(group, dsName, data, max_digit=55, compression=COMPRESSION):
    """Create HDF5 dataset and print out message."""
    msg = 'create dataset {d:<{w}}'.format(d='{}/{}'.format(group.name, dsName), w=max_digit)
    msg += ' of {t:<10} in size of {s}'.format(t=str(data.dtype), s=data.shape)
    msg += ' with compression={c}'.format(c=compression)
    print(msg)

    if data.ndim == 1:
        dset = group.create_dataset(dsName,
                                    data=data,
                                    compression=compression)
    elif data.ndim == 2:
        dset = group.create_dataset(dsName,
                                    data=data,
                                    chunks=True,
                                    compression=compression)
    return dset


def write_hdf5_file(metadata, out_file, ts_file, tcoh_file, scoh_file, mask_file, geom_file):
    """Write HDF5 file in HDF-EOS5 format."""
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    dateList = ts_obj.dateList
    numDate = len(dateList)

    # Open HDF5 File
    f = h5py.File(out_file, 'w')
    print('create HDF5 file: {} with w mode'.format(out_file))
    max_digit = 55

    ##### Group - Observation
    gName = 'HDFEOS/GRIDS/timeseries/observation'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    ## O1 - displacement
    dsName = 'displacement'
    dsShape = (numDate, ts_obj.length, ts_obj.width)
    dsDataType = np.float32
    print(('create dataset /{d:<{w}} of {t:<10} in size of {s}'
           ' with compression={c}').format(d='{}/{}'.format(gName, dsName),
                                           w=max_digit,
                                           t='float32',
                                           s=dsShape,
                                           c=COMPRESSION))
    dset = group.create_dataset(dsName,
                                shape=dsShape,
                                maxshape=(None, dsShape[1], dsShape[2]),
                                dtype=dsDataType,
                                chunks=True,
                                compression=COMPRESSION)

    print('write data acquition by acquition ...')
    prog_bar = ptime.progressBar(maxValue=numDate)
    for i in range(numDate):
        dset[i, :, :] = readfile.read(ts_file, datasetName=dateList[i])[0]
        prog_bar.update(i+1, suffix='{}/{} {}'.format(i+1, numDate, dateList[i]))
    prog_bar.close()

    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'meters'

    ## O2 - date
    dsName = 'date'
    data = np.array(dateList, dtype=np.string_)
    dset = create_hdf5_dataset(group, dsName, data)

    ## O3 - perp baseline
    dsName = 'bperp'
    data = np.array(ts_obj.pbase, dtype=np.float32)
    dset = create_hdf5_dataset(group, dsName, data)

    ##### Group - Quality
    gName = 'HDFEOS/GRIDS/timeseries/quality'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    ## Q1 - temporalCoherence
    dsName = 'temporalCoherence'
    # read
    data = readfile.read(tcoh_file)[0]
    # write
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'

    ## Q2 - avgSpatialCoherence
    dsName = 'avgSpatialCoherence'
    # read
    data = readfile.read(scoh_file)[0]
    # write
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['_FillValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'

    ## Q3 - mask
    dsName = 'mask'
    # read
    data = readfile.read(mask_file, datasetName='mask')[0]
    # write
    dset = create_hdf5_dataset(group, dsName, data)
    # attributes
    dset.attrs['Title'] = dsName
    dset.attrs['MissingValue'] = BOOL_ZERO
    dset.attrs['_FillValue'] = BOOL_ZERO
    dset.attrs['Units'] = '1'

    ##### Group - Write Geometry
    # Required: height, incidenceAngle
    # Optional: rangeCoord, azimuthCoord, azimuthAngle, slantRangeDistance, waterMask, shadowMask
    gName = 'HDFEOS/GRIDS/timeseries/geometry'
    print('create group   /{}'.format(gName))
    group = f.create_group(gName)

    geom_obj = geometry(geom_file)
    geom_obj.open(print_msg=False)
    for dsName in geom_obj.datasetNames:
        # read
        data = geom_obj.read(datasetName=dsName, print_msg=False)
        # write
        dset = create_hdf5_dataset(group, dsName, data)

        # attributes
        dset.attrs['Title'] = dsName
        if dsName in ['height', 'slantRangeDistance', 'bperp']:
            dset.attrs['MissingValue'] = FLOAT_ZERO
            dset.attrs['_FillValue'] = FLOAT_ZERO
            dset.attrs['Units'] = 'meters'

        elif dsName in ['incidenceAngle', 'azimuthAngle', 'latitude', 'longitude']:
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
    inps, template = read_template2inps(inps.template_file, inps)

    # Prepare Metadata
    meta = prep_metadata(
        ts_file=inps.ts_file,
        template=template,
        print_msg=True)

    # Get output filename
    out_file = get_output_filename(
        metadata=meta,
        suffix=inps.suffix,
        update_mode=inps.update,
        subset_mode=inps.subset)

    # Open HDF5 File
    write_hdf5_file(
        metadata=meta,
        out_file=out_file,
        ts_file=inps.ts_file,
        tcoh_file=inps.tcoh_file,
        scoh_file=inps.scoh_file,
        mask_file=inps.mask_file,
        geom_file=inps.geom_file)

    return


################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

