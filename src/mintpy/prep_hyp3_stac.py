############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################
import datetime as dt
import time

import h5py
import numpy as np
import pystac
import stackstac
from osgeo import osr

from mintpy.constants import SPEED_OF_LIGHT
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile
from mintpy.utils import utils1 as ut


osr.UseExceptions()


#########################################################################
def get_metadata(dataset):
    keys = list(dataset.coords.keys())
    hyp3_meta = {}
    for key in keys:
        if key in ['time', 'x', 'y', 'band']:
            continue

        value = dataset.coords[key].values
        if value.shape == ():
            value = value.item()
        else:
            # value = list(dict.fromkeys(value))
            value = list(value)

        hyp3_meta[key] = value

    # Add geospatial metadata
    meta = {}
    n_dates, n_bands, meta['LENGTH'], meta['WIDTH'] = dataset.shape
    example_image = dataset.isel(time=0)
    meta['X_STEP'], _, meta['X_FIRST'], _, meta['Y_STEP'], meta['Y_FIRST'], *_ = dataset.attrs['transform']
    meta['DATA_TYPE'] = example_image['data_type'].values.item()
    meta['EPSG'] = example_image['epsg'].values.item()
    meta['X_UNIT'] = 'meters'
    meta['Y_UNIT'] = 'meters'
    meta['NoDataValue'] = example_image['nodata'].values.item()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(meta['EPSG'])
    meta['UTM_ZONE'] = srs.GetName().split(' ')[-1]

    # add universal hyp3 metadata
    meta['PROCESSOR'] = 'hyp3'
    meta['ALOOKS'] = hyp3_meta['azimuth_looks']
    meta['RLOOKS'] = hyp3_meta['range_looks']
    meta['EARTH_RADIUS'] = hyp3_meta['earth_radius_at_nadir']
    meta['HEIGHT'] = hyp3_meta['spacecraft_height']
    meta['STARTING_RANGE'] = np.mean(hyp3_meta['slant_range_near'])
    meta['CENTER_LINE_UTC'] = np.mean(hyp3_meta['utc_time'])
    meta['HEADING'] = np.mean(hyp3_meta['heading']) % 360.0 - 360.0  # ensure negative value for the heading angle

    # add LAT/LON_REF1/2/3/4 based on whether satellite ascending or descending
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

    # convert UTM to lat/lon
    N, W = ut.utm2latlon(meta, W, N)
    S, E = ut.utm2latlon(meta, E, S)

    meta['ORBIT_DIRECTION'] = hyp3_meta['reference_orbit_direction'].upper()
    if meta['ORBIT_DIRECTION'] == 'ASCENDING':
        meta['LAT_REF1'] = str(S)
        meta['LAT_REF2'] = str(S)
        meta['LAT_REF3'] = str(N)
        meta['LAT_REF4'] = str(N)
        meta['LON_REF1'] = str(W)
        meta['LON_REF2'] = str(E)
        meta['LON_REF3'] = str(W)
        meta['LON_REF4'] = str(E)
    else:
        meta['LAT_REF1'] = str(N)
        meta['LAT_REF2'] = str(N)
        meta['LAT_REF3'] = str(S)
        meta['LAT_REF4'] = str(S)
        meta['LON_REF1'] = str(E)
        meta['LON_REF2'] = str(W)
        meta['LON_REF3'] = str(E)
        meta['LON_REF4'] = str(W)

    if hyp3_meta['reference_granule'][0].startswith('S1'):
        meta['PLATFORM'] = 'Sen'
        meta['ANTENNA_SIDE'] = -1
        meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
        meta['RANGE_PIXEL_SIZE'] = sensor.SEN['range_pixel_size'] * int(meta['RLOOKS'])
        meta['AZIMUTH_PIXEL_SIZE'] = sensor.SEN['azimuth_pixel_size'] * int(meta['ALOOKS'])
    else:
        raise NotImplementedError('Only Sentinel-1 data is currently supported')

    meta = readfile.standardize_metadata(meta)

    date1s = [dt.datetime.fromisoformat(x).strftime('%Y%m%d') for x in hyp3_meta['start_datetime']]
    date2s = [dt.datetime.fromisoformat(x).strftime('%Y%m%d') for x in hyp3_meta['end_datetime']]
    date12s = [f'{d1}_{d2}' for d1, d2 in zip(date1s, date2s)]

    perp_baseline = np.abs(hyp3_meta['baseline'])
    return meta, date12s, perp_baseline


def write_ifgram_stack(outfile, dataset, date12s, perp_baselines):
    with h5py.File(outfile, 'a') as f:
        f['date'][:, 0] = [d1.split('_')[0].encode('utf-8') for d1 in date12s]
        f['date'][:, 1] = [d2.split('_')[1].encode('utf-8') for d2 in date12s]
        f['dropIfgram'][:] = True
        f['bperp'][:] = perp_baselines
        f['unwrapPhase'][:, :, :] = dataset.sel(band='unw_phase').to_numpy().astype(np.float32)
        f['coherence'][:, :, :] = dataset.sel(band='corr').to_numpy().astype(np.float32)
        f['connectComponent'][:, :, :] = dataset.sel(band='conncomp').to_numpy().astype(int)

        # add MODIFICATION_TIME metadata to each 3D dataset
        for dsName in ['unwrapPhase', 'coherence', 'connectComponent']:
            f[dsName].attrs['MODIFICATION_TIME'] = str(time.time())


def write_geometry(outfile, dastaset):
    first_product = dastaset.isel(time=0)

    # Convert from hyp3/gamma to mintpy/isce2 convention
    incidence_angle = first_product.sel(band='lv_theta').to_numpy().astype(np.float32)
    incidence_angle[incidence_angle == 0] = np.nan
    incidence_angle = 90 - (incidence_angle * 180 / np.pi)

    # Calculate Slant Range distance
    atr = {
        'HEIGHT': first_product.coords['spacecraft_height'].item(),
        'EARTH_RADIUS': first_product.coords['earth_radius_at_nadir'].item(),
    }
    slant_range_distance = ut.incidence_angle2slant_range_distance(atr, incidence_angle)

    # Convert from hyp3/gamma to mintpy/isce2 convention
    azimuth_angle = first_product.sel(band='lv_phi').to_numpy().astype(np.float32)
    azimuth_angle[azimuth_angle == 0] = np.nan
    azimuth_angle = azimuth_angle * 180 / np.pi - 90  # hyp3/gamma to mintpy/isce2 convention
    azimuth_angle = ut.wrap(azimuth_angle, wrap_range=[-180, 180])  # rewrap within -180 to 180

    with h5py.File(outfile, 'a') as f:
        f['height'][:, :] = first_product.sel(band='dem').to_numpy().astype(np.float32)
        f['incidenceAngle'][:, :] = incidence_angle
        f['slantRangeDistance'][:, :] = slant_range_distance
        f['azimuthAngle'][:, :] = azimuth_angle
        f['waterMask'][:, :] = first_product.sel(band='water_mask').to_numpy().astype(np.bool_)


def load_hyp3_stac(
    stac_file,
    subset_yx=None,
    subset_geo=None,
    compression=None,
    ifg_outfile='./inputs/ifgramStack.h5',
    geom_outfile='./inputs/geometryGeo.h5',
):
    collection = pystac.Collection.from_file(stac_file)
    items = list(collection.get_all_items())
    dataset = stackstac.stack(items[:3])

    if subset_geo and subset_yx:
        print('Both geographic and index subsets were provided. Using geographic subset method.')

    if subset_geo:
        dataset = dataset.sel(y=slice(subset_geo[0], subset_geo[1]), x=slice(subset_geo[2], subset_geo[3]))
    elif subset_yx:
        dataset = dataset.isel(y=slice(subset_yx[0], subset_yx[1]), x=slice(subset_yx[2], subset_yx[3]))

    meta, date12s, perp_baselines = get_metadata(dataset)

    num_pair, _, length, width = dataset.shape

    # create ifgramStack.h5 file
    ds_name_dict = {
        'date': (np.dtype('S8'), (num_pair, 2)),
        'dropIfgram': (np.bool_, (num_pair,)),
        'bperp': (np.float32, (num_pair,)),
        'unwrapPhase': (np.float32, (num_pair, length, width)),
        'coherence': (np.float32, (num_pair, length, width)),
        'connectComponent': (int, (num_pair, length, width)),
    }
    meta['FILE_TYPE'] = 'ifgramStack'
    writefile.layout_hdf5(ifg_outfile, ds_name_dict, metadata=meta, compression=compression)
    write_ifgram_stack(ifg_outfile, dataset, date12s, perp_baselines)

    # create geometryGeo.h5 file
    ds_name_dict = {
        'height': (np.float32, (length, width)),
        'incidenceAngle': (np.float32, (length, width)),
        'azimuthAngle': (np.float32, (length, width)),
        'slantRangeDistance': (np.float32, (length, width)),
        'waterMask': (np.bool_, (length, width)),
    }
    meta['FILE_TYPE'] = 'geometry'
    writefile.layout_hdf5(geom_outfile, ds_name_dict, metadata=meta, compression=compression)
    write_geometry(geom_outfile, dataset)


def prep_hyp3_stac(inps):
    load_hyp3_stac(inps.stacFile, inps.yx, inps.lalo, inps.compression, inps.outfile[0], inps.outfile[1])
