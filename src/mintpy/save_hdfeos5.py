############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import datetime as dt
import os

import h5py
import numpy as np

from mintpy import info
from mintpy.objects import geometry, sensor, timeseries
from mintpy.utils import ptime, readfile, utils as ut

BOOL_ZERO = np.bool_(0)
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)
COMPRESSION = 'lzf'


################################################################
def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""

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
def prep_metadata(ts_file, geom_file, template=None, print_msg=True):
    """Prepare metadata for HDF-EOS5 file."""
    # read metadata from ts_file
    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    meta = dict(ts_obj.metadata)

    # read metadata from template_file
    if template:
        for key, value in template.items():
            if not key.startswith(('mintpy', 'isce')):
                meta[key] = value

    # grab unavco metadata
    unavco_meta = metadata_mintpy2unavco(meta, ts_obj.dateList, geom_file)
    if print_msg:
        print('## UNAVCO Metadata:')
        print('-----------------------------------------')
        info.print_attributes(unavco_meta)
        print('-----------------------------------------')

    # update metadata from unavco metadata
    meta.update(unavco_meta)
    meta['FILE_TYPE'] = 'HDFEOS'

    return meta


def metadata_mintpy2unavco(meta_in, dateList, geom_file):
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
    # footprint for actual data coverage in lat/lon bounding box.
    if 'Y_FIRST' in meta.keys():
        # time-series in geo-coordinates
        N = float(meta['Y_FIRST'])
        W = float(meta['X_FIRST'])
        S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
        E = W + float(meta['X_STEP']) * int(meta['WIDTH'])
        unavco_meta['data_footprint'] = ut.snwe_to_wkt_polygon([S, N, W, E])

    else:
        # time-series in radar-coordinates
        geom_meta = readfile.read_attribute(geom_file)
        geom_dset_list = readfile.get_dataset_list(geom_file)
        # potential extra geometry file (for roipac/gamma)
        geo_geom_file = os.path.join(os.path.dirname(geom_file), 'geometryGeo.h5')

        if 'Y_FIRST' not in geom_meta.keys() and 'latitude' in geom_dset_list:
            # geometry in radar-coodinates (isce/doris)
            lat_data = readfile.read(geom_file, datasetName='latitude')[0]
            lon_data = readfile.read(geom_file, datasetName='longitude')[0]

            # set pixels with invalid value or zero to nan
            lat_data[np.abs(lat_data) == 90] = np.nan
            lat_data[lat_data == 0] = np.nan
            lon_data[lon_data == 0] = np.nan

            S, N = np.nanmin(lat_data), np.nanmax(lat_data)
            W, E = np.nanmin(lon_data), np.nanmax(lon_data)
            unavco_meta['data_footprint'] = ut.snwe_to_wkt_polygon([S, N, W, E])

        elif os.path.isfile(geo_geom_file):
            # geometry in geo-coordinates (roipac/gamma)
            geom_meta = readfile.read_attribute(geo_geom_file)

            N = float(geom_meta['Y_FIRST'])
            W = float(geom_meta['X_FIRST'])
            S = N + float(geom_meta['Y_STEP']) * int(geom_meta['LENGTH'])
            E = W + float(geom_meta['X_STEP']) * int(geom_meta['WIDTH'])
            unavco_meta['data_footprint'] = ut.snwe_to_wkt_polygon([S, N, W, E])

        else:
            msg = 'WARNING: "data_footprint" is NOT assigned, '
            msg += 'due to the lack of X/Y_FIRST attributes and latitude/longitde datasets.'
            print(msg)

    return unavco_meta


def get_output_filename(metadata, suffix=None, update_mode=False, subset_mode=False):
    """Get output file name of HDF-EOS5 time-series file."""
    SAT = metadata['mission']
    SW = metadata['beam_mode']
    if metadata['beam_swath']:
        SW += str(metadata['beam_swath'])
    RELORB = "{:03d}".format(int(metadata['relative_orbit']))

    # First and/or Last Frame
    frame1 = metadata['first_frame']
    frame2 = metadata['last_frame']
    FRAME = f"{int(frame1):04d}"
    if frame2 != frame1:
        FRAME += f"_{frame2:04d}"

    DATE1 = dt.datetime.strptime(metadata['first_date'], '%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.datetime.strptime(metadata['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    if update_mode:
        print('Update mode is ON, put endDate as XXXXXXXX.')
        DATE2 = 'XXXXXXXX'

    if suffix:
        outName = f'{SAT}_{SW}_{RELORB}_{FRAME}_{DATE1}_{DATE2}_{suffix}.he5'
    else:
        outName = f'{SAT}_{SW}_{RELORB}_{FRAME}_{DATE1}_{DATE2}.he5'

    if subset_mode:
        print('Subset mode is enabled, put subset range info in output filename.')
        lat1 = float(metadata['Y_FIRST'])
        lon0 = float(metadata['X_FIRST'])
        lat0 = lat1 + float(metadata['Y_STEP']) * int(metadata['LENGTH'])
        lon1 = lon0 + float(metadata['X_STEP']) * int(metadata['WIDTH'])

        lat0Str = f'N{round(lat0*1e3):05d}'
        lat1Str = f'N{round(lat1*1e3):05d}'
        lon0Str = f'E{round(lon0*1e3):06d}'
        lon1Str = f'E{round(lon1*1e3):06d}'
        if lat0 < 0.0: lat0Str = f'S{round(abs(lat0)*1e3):05d}'
        if lat1 < 0.0: lat1Str = f'S{round(abs(lat1)*1e3):05d}'
        if lon0 < 0.0: lon0Str = f'W{round(abs(lon0)*1e3):06d}'
        if lon1 < 0.0: lon1Str = f'W{round(abs(lon1)*1e3):06d}'

        SUB = f'_{lat0Str}_{lat1Str}_{lon0Str}_{lon1Str}'
        fbase, fext = os.path.splitext(outName)
        outName = f'{fbase}{SUB}{fext}'

    return outName


def create_hdf5_dataset(group, dsName, data, max_digit=55, compression=COMPRESSION):
    """Create HDF5 dataset and print out message."""

    msg = 'create dataset {d:<{w}}'.format(d=f'{group.name}/{dsName}', w=max_digit)
    msg += f' of {str(data.dtype):<10} in size of {data.shape} with compression={compression}'
    print(msg)

    if data.ndim == 1:
        dset = group.create_dataset(
            dsName,
            data=data,
            compression=compression,
        )

    elif data.ndim == 2:
        dset = group.create_dataset(
            dsName,
            data=data,
            chunks=True,
            compression=compression,
        )

    return dset


def write_hdf5_file(metadata, out_file, ts_file, tcoh_file, scoh_file, mask_file, geom_file):
    """Write HDF5 file in HDF-EOS5 format."""

    ts_obj = timeseries(ts_file)
    ts_obj.open(print_msg=False)
    dateList = ts_obj.dateList
    numDate = len(dateList)

    # Open HDF5 File
    print(f'create HDF5 file: {out_file} with w mode')
    max_digit = 55

    with h5py.File(out_file, 'w') as f:

        ##### Group - Observation
        gName = 'HDFEOS/GRIDS/timeseries/observation'
        print(f'create group   /{gName}')
        group = f.create_group(gName)

        ## O1 - displacement
        dsName = 'displacement'
        dsShape = (numDate, ts_obj.length, ts_obj.width)
        dsDataType = np.float32
        msg = 'create dataset /{d:<{w}}'.format(d=f'{gName}/{dsName}', w=max_digit)
        msg += f' of {"float32":<10} in size of {dsShape} with compression={COMPRESSION}'
        print(msg)

        dset = group.create_dataset(
            dsName,
            shape=dsShape,
            maxshape=(None, dsShape[1], dsShape[2]),
            dtype=dsDataType,
            chunks=True,
            compression=COMPRESSION,
        )

        print('write data acquition by acquition ...')
        prog_bar = ptime.progressBar(maxValue=numDate)
        for i in range(numDate):
            dset[i, :, :] = readfile.read(ts_file, datasetName=dateList[i])[0]
            prog_bar.update(i+1, suffix=f'{i+1}/{numDate} {dateList[i]}')
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
        print(f'create group   /{gName}')
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
        # Optional: rangeCoord, azimuthCoord, azimuthAngle, slantRangeDistance,
        #           waterMask, shadowMask
        gName = 'HDFEOS/GRIDS/timeseries/geometry'
        print(f'create group   /{gName}')
        group = f.create_group(gName)

        geom_obj = geometry(geom_file)
        geom_obj.open(print_msg=False)

        # add latitude/longitude if missing, e.g. ARIA/HyP3
        dsNames = geom_obj.datasetNames + ['latitude', 'longitude']
        dsNames = list(set(dsNames))

        for dsName in dsNames:
            # read
            if dsName in geom_obj.datasetNames:
                data = geom_obj.read(datasetName=dsName, print_msg=False)
            elif dsName == 'latitude':
                data = ut.get_lat_lon(metadata, dimension=2)[0]
            elif dsName == 'longitude':
                data = ut.get_lat_lon(metadata, dimension=2)[1]
            else:
                raise ValueError(f'Un-recognized dataset name: {dsName}!')

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

    print(f'finished writing to {out_file}')

    return out_file


################################################################
def save_hdfeos5(inps):

    inps, template = read_template2inps(inps.template_file, inps)

    # prepare metadata
    meta = prep_metadata(
        ts_file=inps.ts_file,
        geom_file=inps.geom_file,
        template=template,
        print_msg=True)

    # get output filename
    out_file = get_output_filename(
        metadata=meta,
        suffix=inps.suffix,
        update_mode=inps.update,
        subset_mode=inps.subset)

    # write HDF5 File
    write_hdf5_file(
        metadata=meta,
        out_file=out_file,
        ts_file=inps.ts_file,
        tcoh_file=inps.tcoh_file,
        scoh_file=inps.scoh_file,
        mask_file=inps.mask_file,
        geom_file=inps.geom_file)

    return
