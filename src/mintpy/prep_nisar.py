############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee, Jul 2023                           #
############################################################

import datetime
import glob
import os
from pathlib import Path

import h5py
import numpy as np
from osgeo import gdal
from pyproj.transformer import Transformer
from scipy.interpolate import RegularGridInterpolator

from mintpy.constants import EARTH_RADIUS, SPEED_OF_LIGHT
from mintpy.utils import ptime, writefile

DATASET_ROOT_UNW = '/science/LSAR/GUNW/grids/frequencyA/interferogram/unwrapped'
IDENTIFICATION = '/science/LSAR/identification'
RADARGRID_ROOT = 'science/LSAR/GUNW/metadata/radarGrid'
DATASETS = {
    'xcoord'           : f"{DATASET_ROOT_UNW}/xCoordinates",
    'ycoord'           : f"{DATASET_ROOT_UNW}/yCoordinates",
    'unw'              : f"{DATASET_ROOT_UNW}/POL/unwrappedPhase",
    'cor'              : f"{DATASET_ROOT_UNW}/POL/coherenceMagnitude",
    'connComp'         : f"{DATASET_ROOT_UNW}/POL/connectedComponents",
    'layoverShadowMask': f"{DATASET_ROOT_UNW}/layoverShadowMask",
    'waterMask'        : f"{DATASET_ROOT_UNW}/waterMask",
    'epsg'             : f"{DATASET_ROOT_UNW}/projection",
    'xSpacing'         : f"{DATASET_ROOT_UNW}/xCoordinateSpacing",
    'ySpacing'         : f"{DATASET_ROOT_UNW}/yCoordinateSpacing",
    'polarization'     : f"{DATASET_ROOT_UNW}/listOfPolarizations",
    'range_look'       : f"{DATASET_ROOT_UNW}/numberOfRangeLooks",
    'azimuth_look'     : f"{DATASET_ROOT_UNW}/numberOfAzimuthLooks",
}
PROCESSINFO = {
    'centerFrequency': "/science/LSAR/GUNW/grids/frequencyA/centerFrequency",
    'orbit_direction': f"{IDENTIFICATION}/orbitPassDirection",
    'platform'       : f"{IDENTIFICATION}/missionId",
    'start_time'     : f"{IDENTIFICATION}/referenceZeroDopplerStartTime",
    'end_time'       : f"{IDENTIFICATION}/referenceZeroDopplerEndTime",
    'rdr_xcoord'     : f"{RADARGRID_ROOT}/xCoordinates",
    'rdr_ycoord'     : f"{RADARGRID_ROOT}/yCoordinates",
    'rdr_slant_range': f"{RADARGRID_ROOT}/slantRange",
    'rdr_height'     : f"{RADARGRID_ROOT}/heightAboveEllipsoid",
    'rdr_incidence'  : f"{RADARGRID_ROOT}/incidenceAngle",
    'bperp'          : f"{RADARGRID_ROOT}/perpendicularBaseline",
}

####################################################################################

def load_nisar(inps):
    """Prepare and load NISAR data and metadata into HDF5/MintPy format."""

    print(f'update mode: {inps.update_mode}')

    input_files = sorted(glob.glob(inps.input_glob))
    print(f"Found {len(input_files)} unwrapped files")

    if inps.subset_lat:
        bbox = (inps.subset_lat[0], inps.subset_lon[0], inps.subset_lat[1], inps.subset_lon[1])
    else:
        bbox=None

    # extract metadata
    metadata, bounds = extract_metadata(input_files, bbox)

    # output filename
    stack_file = os.path.join(inps.out_dir, "inputs/ifgramStack.h5")
    geometry_file = os.path.join(inps.out_dir, "inputs/geometryGeo.h5")

    # get date info: date12_list
    date12_list = _get_date_pairs(input_files)

    metadata = prepare_geometry(
        outfile=geometry_file,
        metaFile=input_files[0],
        bbox=bounds,
        metadata=metadata,
        demFile=inps.dem_file
        )

    prepare_stack(
        outfile=stack_file,
        inp_files=input_files,
        metadata=metadata,
        bbox=bounds,
        date12_list=date12_list,
    )

    print("Done.")

    return


def extract_metadata(input_files, bbox=None, polarization='HH'):
    """Extract NISAR metadata for MintPy."""
    meta_file = input_files[0]
    meta = {}

    # update dataset
    for key, value in DATASETS.items():
        if value:
            value = value.replace("POL", polarization)
            DATASETS[key] = value

    with h5py.File(meta_file, 'r') as ds:
        pixel_height = ds[DATASETS['ySpacing']][()]
        pixel_width = ds[DATASETS['xSpacing']][()]
        x_origin = min(ds[DATASETS['xcoord']][()])
        y_origin = max(ds[DATASETS['ycoord']][()])
        xcoord = ds[DATASETS['xcoord']][()]
        ycoord = ds[DATASETS['ycoord']][()]
        meta["EPSG"] = ds[DATASETS['epsg']][()]
        meta['WAVELENGTH'] = SPEED_OF_LIGHT / ds[PROCESSINFO['centerFrequency']][()]
        meta["ORBIT_DIRECTION"] = ds[PROCESSINFO['orbit_direction']][()].decode('utf-8')
        meta['POLARIZATION'] = polarization
        meta["ALOOKS"] = ds[DATASETS['azimuth_look']][()]
        meta["RLOOKS"] = ds[DATASETS['range_look']][()]
        meta['PLATFORM'] = ds[PROCESSINFO['platform']][()].decode('utf-8')
        meta['STARTING_RANGE'] = min(ds[PROCESSINFO['rdr_slant_range']][()].flatten())
        start_time = datetime.datetime.strptime(ds[PROCESSINFO['start_time']][()].decode('utf-8')[0:26],
                                                "%Y-%m-%dT%H:%M:%S.%f")
        end_time = datetime.datetime.strptime(ds[PROCESSINFO['end_time']][()].decode('utf-8')[0:26],
                                              "%Y-%m-%dT%H:%M:%S.%f")

    t_mid = start_time + (end_time - start_time) / 2.0
    meta["CENTER_LINE_UTC"] = (
            t_mid - datetime.datetime(t_mid.year, t_mid.month, t_mid.day)
    ).total_seconds()

    meta["X_FIRST"] = x_origin - pixel_width//2
    meta["Y_FIRST"] = y_origin - pixel_height//2
    meta["X_STEP"] = pixel_width
    meta["Y_STEP"] = pixel_height
    meta["X_UNIT"] = meta["Y_UNIT"] = "meters"
    meta["EARTH_RADIUS"] = EARTH_RADIUS

    # NISAR Altitude
    meta['HEIGHT'] = 747000

    # Range and Azimuth pixel size need revision, values are just to fill in
    meta["RANGE_PIXEL_SIZE"] = abs(pixel_width)
    meta["AZIMUTH_PIXEL_SIZE"] = abs(pixel_height)

    # get the common raster bound among input files
    if bbox:
        # assuming bbox is in lat/lon coordinates
        epsg_src = 4326
        utm_bbox = bbox_to_utm(bbox, meta["EPSG"], epsg_src)
    else:
        utm_bbox = None
    bounds = common_raster_bound(input_files, utm_bbox)
    meta['bbox'] = ",".join([str(b) for b in bounds])

    col1, row1, col2, row2 = get_rows_cols(xcoord, ycoord, bounds)
    length = row2 - row1
    width = col2 - col1
    meta['LENGTH'] = length
    meta['WIDTH'] = width

    return meta, bounds


def get_rows_cols(xcoord, ycoord, bounds):
    """Get row and cols of the bounding box to subset"""
    xindex = np.where(np.logical_and(xcoord >= bounds[0], xcoord <= bounds[2]))[0]
    yindex = np.where(np.logical_and(ycoord >= bounds[1], ycoord <= bounds[3]))[0]
    row1, row2 = min(yindex), max(yindex)
    col1, col2 = min(xindex), max(xindex)
    return col1, row1, col2, row2


def get_raster_corners(input_file):
    """Get the (west, south, east, north) bounds of the image."""
    with h5py.File(input_file, 'r') as ds:
        xcoord = ds[DATASETS['xcoord']][:]
        ycoord = ds[DATASETS['ycoord']][:]
        west  = max(min(ds[PROCESSINFO['rdr_xcoord']][:]), min(xcoord))
        east  = min(max(ds[PROCESSINFO['rdr_xcoord']][:]), max(xcoord))
        north = min(max(ds[PROCESSINFO['rdr_ycoord']][:]), max(ycoord))
        south = max(min(ds[PROCESSINFO['rdr_ycoord']][:]), min(ycoord))
    return west, south, east, north


def common_raster_bound(input_files, utm_bbox=None):
    """Get common bounds among all data"""
    x_bounds = []
    y_bounds = []
    for file in input_files:
        west, south, east, north = get_raster_corners(file)
        x_bounds.append([west, east])
        y_bounds.append([south, north])
    if not utm_bbox is None:
        x_bounds.append([utm_bbox[0], utm_bbox[2]])
        y_bounds.append([utm_bbox[1], utm_bbox[3]])

    bounds = max(x_bounds)[0], max(y_bounds)[0], min(x_bounds)[1], min(y_bounds)[1]
    return bounds


def bbox_to_utm(bbox, epsg_dst, epsg_src=4326):
    """Convert a list of points to a specified UTM coordinate system.
        If epsg_src is 4326 (lat/lon), assumes points_xy are in degrees.
    """
    xmin, ymin, xmax, ymax = bbox
    t = Transformer.from_crs(epsg_src, epsg_dst, always_xy=True)
    xs = [xmin, xmax]
    ys = [ymin, ymax]
    xt, yt = t.transform(xs, ys)
    xys = list(zip(xt, yt))
    return *xys[0], *xys[1]


def read_subset(inp_file, bbox, geometry=False):
    """Read a subset of data using bounding box in rows and cols"""
    dataset = {}
    with h5py.File(inp_file, 'r') as ds:
        xcoord = ds[DATASETS['xcoord']][:]
        ycoord = ds[DATASETS['ycoord']][:]
        col1, row1, col2, row2 = get_rows_cols(xcoord, ycoord, bbox)

        if geometry:
            # Set all values to 1 temporarily because water mask is zero
            dataset['water_mask'] = ds[DATASETS['waterMask']][row1:row2, col1:col2] * 0 + 1
            dataset['layover_shadow_mask'] = ds[DATASETS['layoverShadowMask']][row1:row2, col1:col2]
            dataset['xybbox'] = (col1, row1, col2, row2)
        else:
            dataset['unw_data'] = ds[DATASETS['unw']][row1:row2, col1:col2]
            dataset['cor_data'] = ds[DATASETS['cor']][row1:row2, col1:col2]
            dataset['conn_comp'] = (ds[DATASETS['connComp']][row1:row2, col1:col2]).astype(np.float32)
            dataset['conn_comp'][dataset['conn_comp'] > 254] = np.nan
            dataset['pbase'] = np.nanmean(ds[PROCESSINFO['bperp']][()])
    return dataset


def read_and_interpolate_geometry(gunw_file, dem_file, xybbox):
    """Read the DEM, change projection and interpolate to data grid, interpolate slant range and incidence angle to data grid"""
    dataset = gdal.Open(dem_file, gdal.GA_ReadOnly)
    geotransform = dataset.GetGeoTransform()
    proj = gdal.osr.SpatialReference(wkt=dataset.GetProjection())
    src_epsg = int(proj.GetAttrValue('AUTHORITY', 1))
    raster_array = dataset.ReadAsArray()
    rdr_coords = {}

    with h5py.File(gunw_file, 'r') as ds:
        dst_epsg = ds[DATASETS['epsg']][()]
        xcoord = ds[DATASETS['xcoord']][xybbox[0]:xybbox[2]]
        ycoord = ds[DATASETS['ycoord']][xybbox[1]:xybbox[3]]

        rdr_coords['xcoord_radar_grid'] = ds[PROCESSINFO['rdr_xcoord']][()]
        rdr_coords['ycoord_radar_grid'] = ds[PROCESSINFO['rdr_ycoord']][()]
        rdr_coords['height_radar_grid'] = ds[PROCESSINFO['rdr_height']][()]
        rdr_coords['slant_range'] = ds[PROCESSINFO['rdr_slant_range']][()]
        rdr_coords['perp_baseline'] = ds[PROCESSINFO['bperp']][()]
        rdr_coords['incidence_angle'] = ds[PROCESSINFO['rdr_incidence']][()]

    subset_rows = len(ycoord)
    subset_cols = len(xcoord)
    # dem_subset_array = np.zeros((subset_rows, subset_cols), dtype=raster_array.dtype)
    Y_2d, X_2d = np.meshgrid(ycoord, xcoord, indexing='ij')

    if not src_epsg == dst_epsg:
        coord_transform = Transformer.from_crs(dst_epsg, src_epsg, always_xy=True)
        x_dem, y_dem = coord_transform.transform(X_2d.flatten(), Y_2d.flatten())
    else:
        x_dem, y_dem = X_2d.flatten(), Y_2d.flatten()

    cols = ((y_dem - geotransform[3]) / geotransform[5]).astype(int)
    rows = ((x_dem - geotransform[0]) / geotransform[1]).astype(int)

    dem_subset_array = raster_array[cols.reshape(subset_rows, subset_cols), rows.reshape(subset_rows, subset_cols)]

    slant_range, incidence_angle = interpolate_geometry(X_2d, Y_2d, dem_subset_array, rdr_coords)

    return dem_subset_array, slant_range, incidence_angle


def interpolate_geometry(X_2d, Y_2d, dem, rdr_coords):
    """Interpolate slant range and incidence angle"""
    pnts = np.stack((dem.flatten(), Y_2d.flatten(), X_2d.flatten()), axis=-1)
    length, width = Y_2d.shape

    # build the interpolator
    interpolator = RegularGridInterpolator((rdr_coords['height_radar_grid'],
                                            rdr_coords['ycoord_radar_grid'],
                                            rdr_coords['xcoord_radar_grid']),
                                            rdr_coords['slant_range'],
                                            method='linear')

    interpolated_slant_range = interpolator(pnts)

    interpolator = RegularGridInterpolator((rdr_coords['height_radar_grid'],
                                            rdr_coords['ycoord_radar_grid'],
                                            rdr_coords['xcoord_radar_grid']),
                                            rdr_coords['incidence_angle'],
                                            method='linear')
    interpolated_incidence_angle = interpolator(pnts)

    return interpolated_slant_range.reshape(length, width), interpolated_incidence_angle.reshape(length, width)


######################################

def _get_date_pairs(filenames):
    str_list = [Path(f).stem for f in filenames]
    return [str(f.split('_')[13]) + '_' + str(f.split('_')[11]) for f in str_list]


def prepare_geometry(
        outfile,
        metaFile,
        metadata,
        bbox,
        demFile
):
    """Prepare the geometry file."""
    print("-" * 50)
    print(f"preparing geometry file: {outfile}")

    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # Read waterMask, LayoverShadowMask, xybbox:
    geo_ds = read_subset(metaFile, bbox, geometry=True)
    dem_subset_array, slant_range, incidence_angle = read_and_interpolate_geometry(metaFile, demFile, geo_ds['xybbox'])

    length, width = dem_subset_array.shape

    ds_name_dict = {
        "height": [np.float32, (length, width), dem_subset_array],
        "incidenceAngle": [np.float32, (length, width), incidence_angle],
        "slantRangeDistance": [np.float32, (length, width), slant_range],
        "shadowMask": [np.bool_, (length, width), geo_ds['layover_shadow_mask']],
        "waterMask": [np.bool_, (length, width), geo_ds['water_mask']],
    }

    # initiate HDF5 file
    meta["FILE_TYPE"] = "geometry"
    meta['STARTING_RANGE'] = np.nanmin(slant_range)
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    return meta


def prepare_stack(
    outfile,
    inp_files,
    metadata,
    bbox,
    date12_list
):
    """Prepare the input unw stack."""
    print("-" * 50)
    print(f"preparing ifgramStack file: {outfile}")
    # copy metadata to meta
    meta = {key: value for key, value in metadata.items()}

    # get list of *.unw file
    num_pair = len(inp_files)

    print(f"number of inputs/unwrapped interferograms: {num_pair}")

    # read baseline data
    pbase = np.zeros(num_pair, dtype=np.float32)

    # size info
    cols, rows = meta['WIDTH'], meta['LENGTH']

    # define (and fill out some) dataset structure
    date12_arr = np.array([x.split("_") for x in date12_list], dtype=np.string_)
    drop_ifgram = np.ones(num_pair, dtype=np.bool_)
    ds_name_dict = {
        "date": [date12_arr.dtype, (num_pair, 2), date12_arr],
        "bperp": [np.float32, (num_pair,), pbase],
        "dropIfgram": [np.bool_, (num_pair,), drop_ifgram],
        "unwrapPhase": [np.float32, (num_pair, rows, cols), None],
        "coherence": [np.float32, (num_pair, rows, cols), None],
        "connectComponent": [np.float32, (num_pair, rows, cols), None],
    }

    # initiate HDF5 file
    meta["FILE_TYPE"] = "ifgramStack"
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print(f"writing data to HDF5 file {outfile} with a mode ...")

    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_pair)
        for i, file in enumerate(inp_files):
            dataset = read_subset(file, bbox)

            # read/write *.unw file
            f["unwrapPhase"][i] = dataset['unw_data']

            # read/write *.cor file
            f["coherence"][i] = dataset['cor_data']

            # read/write *.unw.conncomp file
            f["connectComponent"][i] = dataset['conn_comp']

            # read/write perpendicular baseline file
            f['bperp'][i] = dataset['pbase']

            prog_bar.update(i + 1, suffix=date12_list[i])
        prog_bar.close()

    print(f"finished writing to HDF5 file: {outfile}")
    return outfile
