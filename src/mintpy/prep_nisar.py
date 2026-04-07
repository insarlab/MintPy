#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee, Jul 2023                           #
#         Emre Havazli, Apr 2026                           #
############################################################

import datetime
import glob
import os
from pathlib import Path

import h5py
import numpy as np
from osgeo import gdal, osr
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

from mintpy.constants import EARTH_RADIUS, SPEED_OF_LIGHT
from mintpy.utils import attribute as attr, ptime, writefile

# ---------------------------------------------------------------------
# Constants / HDF5 paths (GUNW frequencyA, unwrappedInterferogram)
# ---------------------------------------------------------------------
DATASET_ROOT_UNW = "/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"
PARAMETERS = (
    "/science/LSAR/GUNW/metadata/processingInformation/parameters/"
    "unwrappedInterferogram/frequencyA"
)
IDENTIFICATION = "/science/LSAR/identification"
RADARGRID_ROOT = "science/LSAR/GUNW/metadata/radarGrid"

DATASETS = {
    "xcoord": f"{DATASET_ROOT_UNW}/POL/xCoordinates",
    "ycoord": f"{DATASET_ROOT_UNW}/POL/yCoordinates",
    "unw": f"{DATASET_ROOT_UNW}/POL/unwrappedPhase",
    "cor": f"{DATASET_ROOT_UNW}/POL/coherenceMagnitude",
    "connComp": f"{DATASET_ROOT_UNW}/POL/connectedComponents",
    "ion": f"{DATASET_ROOT_UNW}/POL/ionospherePhaseScreen",
    "epsg": f"{DATASET_ROOT_UNW}/POL/projection",
    "xSpacing": f"{DATASET_ROOT_UNW}/POL/xCoordinateSpacing",
    "ySpacing": f"{DATASET_ROOT_UNW}/POL/yCoordinateSpacing",
    "polarization": "/science/LSAR/GUNW/grids/frequencyA/listOfPolarizations",
    "range_look": f"{PARAMETERS}/numberOfRangeLooks",
    "azimuth_look": f"{PARAMETERS}/numberOfAzimuthLooks",
}

PROCESSINFO = {
    "centerFrequency": "/science/LSAR/GUNW/grids/frequencyA/centerFrequency",
    "orbit_direction": f"{IDENTIFICATION}/orbitPassDirection",
    "platform": f"{IDENTIFICATION}/missionId",
    "start_time": f"{IDENTIFICATION}/referenceZeroDopplerStartTime",
    "end_time": f"{IDENTIFICATION}/referenceZeroDopplerEndTime",
    "rdr_xcoord": f"{RADARGRID_ROOT}/xCoordinates",
    "rdr_ycoord": f"{RADARGRID_ROOT}/yCoordinates",
    "rdr_slant_range": f"{RADARGRID_ROOT}/referenceSlantRange",
    "rdr_height": f"{RADARGRID_ROOT}/heightAboveEllipsoid",
    "rdr_incidence": f"{RADARGRID_ROOT}/incidenceAngle",
    "rdr_los_x": f"{RADARGRID_ROOT}/losUnitVectorX",
    "rdr_los_y": f"{RADARGRID_ROOT}/losUnitVectorY",
    "rdr_wet_tropo": f"{RADARGRID_ROOT}/wetTroposphericPhaseScreen",
    "rdr_hs_tropo": f"{RADARGRID_ROOT}/hydrostaticTroposphericPhaseScreen",
    "rdr_SET": f"{RADARGRID_ROOT}/slantRangeSolidEarthTidesPhase",
    "bperp": f"{RADARGRID_ROOT}/perpendicularBaseline",
}

STACK_TYPES = {"ifgram", "ion", "tropo", "set"}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def _datasets_for_pol(polarization: str) -> dict:
    """Return a per-call datasets dict without mutating the global DATASETS."""
    out = {}
    for k, v in DATASETS.items():
        out[k] = (
            v.replace("POL", polarization) if isinstance(v, str) and "POL" in v else v
        )
    return out


def _grid_bounds_from_xy(xcoord: np.ndarray, ycoord: np.ndarray):
    """
    Compute pixel-edge bounds aligned to the xcoord/ycoord grid (pixel centers).
    Returns bounds in (minx, miny, maxx, maxy) and dx, dy (signed spacings).
    """
    if xcoord.size < 2 or ycoord.size < 2:
        raise ValueError(
            "xcoord/ycoord must have at least 2 elements to infer spacing."
        )

    dx = float(xcoord[1] - xcoord[0])
    dy = float(ycoord[1] - ycoord[0])

    left = float(xcoord[0] - dx / 2.0)
    right = float(xcoord[-1] + dx / 2.0)

    top = float(ycoord[0] - dy / 2.0)
    bottom = float(ycoord[-1] + dy / 2.0)

    miny, maxy = (bottom, top) if bottom < top else (top, bottom)
    minx, maxx = (left, right) if left < right else (right, left)
    return (minx, miny, maxx, maxy), dx, dy


def _warp_to_grid_mem(
    *,
    src_path: str,
    src_epsg: int,
    dst_epsg: int,
    xcoord: np.ndarray,
    ycoord: np.ndarray,
    resample_alg: str,
):
    """
    Warp a raster to the exact xcoord/ycoord grid using MEM output.
    Uses bounds derived from pixel-edge and xRes/yRes with targetAlignedPixels.
    """
    bounds, dx, dy = _grid_bounds_from_xy(xcoord, ycoord)

    warp_opts = gdal.WarpOptions(
        format="MEM",
        outputBounds=bounds,
        srcSRS=f"EPSG:{src_epsg}",
        dstSRS=f"EPSG:{dst_epsg}",
        xRes=abs(dx),
        yRes=abs(dy),
        targetAlignedPixels=True,
        resampleAlg=resample_alg,
    )
    dst = gdal.Warp("", src_path, options=warp_opts)
    if dst is None:
        raise RuntimeError(f"GDAL Warp failed for {src_path}")
    arr = dst.ReadAsArray()
    if arr is None:
        raise RuntimeError(f"Failed reading warped array for {src_path}")
    return arr


def _read_raster_epsg(path: str) -> int:
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    if ds is None:
        raise OSError(f"Cannot open raster: {path}")
    projection = ds.GetProjection()
    if not projection:
        raise ValueError(f"Raster has no projection metadata: {path}")

    srs = osr.SpatialReference()
    if srs.ImportFromWkt(projection) != 0:
        raise ValueError(
            f"Could not parse raster projection WKT for {path}: {projection!r}"
        )

    srs.AutoIdentifyEPSG()
    epsg = srs.GetAuthorityCode(None)
    if epsg is None:
        for authority_node in ["PROJCS", "GEOGCS"]:
            epsg = srs.GetAuthorityCode(authority_node)
            if epsg is not None:
                break

    if epsg is None:
        raise ValueError(
            f"Could not determine EPSG from raster projection for {path}: {projection!r}"
        )

    return int(epsg)


def _make_rgi(grid_axes, values, method="linear"):
    """
    RegularGridInterpolator wrapper:
      - flips decreasing axes (SciPy requires increasing)
      - bounds_error=False + fill_value=np.nan to avoid crashing
    """
    axes = [np.asarray(a) for a in grid_axes]
    vals = values
    for dim, ax in enumerate(axes):
        if ax[0] > ax[-1]:
            axes[dim] = ax[::-1]
            vals = np.flip(vals, axis=dim)

    return RegularGridInterpolator(
        tuple(axes),
        vals,
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )


def _coerce_subset_metadata_types(meta):
    """Keep subset-updated metadata numeric for downstream array sizing."""
    for key in [
        "LENGTH",
        "WIDTH",
        "XMAX",
        "YMAX",
        "SUBSET_XMIN",
        "SUBSET_XMAX",
        "SUBSET_YMIN",
        "SUBSET_YMAX",
    ]:
        if key in meta:
            meta[key] = int(meta[key])
    for key in ["X_FIRST", "Y_FIRST", "X_STEP", "Y_STEP"]:
        if key in meta:
            meta[key] = float(meta[key])
    return meta


def _read_valid_unw_mask(gunw_file: str, xybbox, pol: str):
    """
    Validity mask is ALWAYS based on finite unwrappedPhase (+ _FillValue check),
    using:
      /science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/{pol}/unwrappedPhase
    """
    path = f"{DATASET_ROOT_UNW}/{pol}/unwrappedPhase"
    with h5py.File(gunw_file, "r") as ds:
        dset = ds[path]
        unw = dset[xybbox[1] : xybbox[3], xybbox[0] : xybbox[2]]
        fill = dset.attrs.get("_FillValue", None)

    valid = np.isfinite(unw)
    if fill is not None:
        valid &= unw != fill
    return valid


def _read_perpendicular_baseline(gunw_file: str) -> np.float32:
    """Read the NISAR perpendicular baseline as one finite mean value."""
    with h5py.File(gunw_file, "r") as ds:
        dset = ds[PROCESSINFO["bperp"]]
        bperp = np.asarray(dset[()], dtype=np.float64).reshape(-1)
        fill = dset.attrs.get("_FillValue", None)

    if fill is not None:
        bperp = np.where(bperp == fill, np.nan, bperp)

    bperp = bperp[np.isfinite(bperp)]
    if bperp.size == 0:
        raise ValueError(
            f"No finite perpendicular baseline values found in {gunw_file}"
        )

    pbase = np.mean(bperp)
    return np.float32(pbase)


def _read_target_grid(gunw_file: str, xybbox, polarization: str):
    """Read the destination EPSG and subset grid axes from a GUNW file."""
    datasets = _datasets_for_pol(polarization)
    with h5py.File(gunw_file, "r") as ds:
        return (
            int(ds[datasets["epsg"]][()]),
            ds[datasets["xcoord"]][xybbox[0] : xybbox[2]],
            ds[datasets["ycoord"]][xybbox[1] : xybbox[3]],
        )


def _read_radar_grid_fields(gunw_file: str, field_map: dict):
    """Read radar-grid interpolation axes plus the requested data fields."""
    rdr_coords = {}
    with h5py.File(gunw_file, "r") as ds:
        rdr_coords["xcoord_radar_grid"] = ds[PROCESSINFO["rdr_xcoord"]][()]
        rdr_coords["ycoord_radar_grid"] = ds[PROCESSINFO["rdr_ycoord"]][()]
        rdr_coords["height_radar_grid"] = ds[PROCESSINFO["rdr_height"]][()]
        for out_key, process_key in field_map.items():
            rdr_coords[out_key] = ds[PROCESSINFO[process_key]][()]
    return rdr_coords


def _prepare_radar_grid_interpolation(
    gunw_file, dem_file, xybbox, polarization, field_map
):
    """Build the common DEM/grid/valid-mask context for radar-grid interpolation."""
    dem_src_epsg = _read_raster_epsg(dem_file)
    dst_epsg, xcoord, ycoord = _read_target_grid(gunw_file, xybbox, polarization)
    rdr_coords = _read_radar_grid_fields(gunw_file, field_map)

    dem_subset_array = _warp_to_grid_mem(
        src_path=dem_file,
        src_epsg=dem_src_epsg,
        dst_epsg=dst_epsg,
        xcoord=xcoord,
        ycoord=ycoord,
        resample_alg="bilinear",
    )

    y_2d, x_2d = np.meshgrid(ycoord, xcoord, indexing="ij")
    valid_mask = _read_valid_unw_mask(gunw_file, xybbox, polarization)

    return {
        "dst_epsg": dst_epsg,
        "xcoord": xcoord,
        "ycoord": ycoord,
        "x_2d": x_2d,
        "y_2d": y_2d,
        "dem": dem_subset_array,
        "valid_mask": valid_mask,
        "rdr_coords": rdr_coords,
    }


def _prepare_valid_interp_points(x_2d, y_2d, dem, valid_mask):
    """Return output shape plus 3D interpolation points for valid pixels only."""
    shape = y_2d.shape
    ii, jj = np.where(valid_mask)
    if ii.size == 0:
        return shape, ii, jj, None

    pts = np.column_stack(
        [
            dem[ii, jj].astype(np.float64),
            y_2d[ii, jj].astype(np.float64),
            x_2d[ii, jj].astype(np.float64),
        ]
    )
    return shape, ii, jj, pts


def _interpolate_radar_grid_field(rdr_coords, field_name, pts):
    """Interpolate one radar-grid field at valid pixel locations."""
    grid = (
        rdr_coords["height_radar_grid"],
        rdr_coords["ycoord_radar_grid"],
        rdr_coords["xcoord_radar_grid"],
    )
    interpolator = _make_rgi(grid, rdr_coords[field_name], method="linear")
    return interpolator(pts)


def _empty_interp_array(shape):
    """Allocate a float32 array initialized with NaNs for interpolation output."""
    return np.full(shape, np.nan, dtype=np.float32)


def _resolve_stack_type(stack_type, outfile):
    """Prefer explicit stack types while keeping legacy filename inference."""
    if stack_type is not None:
        if stack_type not in STACK_TYPES:
            raise ValueError(
                f"Unsupported stack_type {stack_type!r}; expected one of {sorted(STACK_TYPES)}"
            )
        return stack_type

    legacy_names = {
        "inputs/ifgramStack.h5": "ifgram",
        "inputs/ionStack.h5": "ion",
        "inputs/tropoStack.h5": "tropo",
        "inputs/setStack.h5": "set",
    }
    for legacy_name, inferred_type in legacy_names.items():
        if legacy_name in outfile:
            return inferred_type

    raise ValueError(
        f"Unable to infer stack_type from outfile {outfile!r}. "
        "Please pass stack_type explicitly."
    )


def _read_stack_observation(file, stack_type, bbox, dem_file, polarization):
    """Read one observation for the requested stack type."""
    pbase = _read_perpendicular_baseline(file)

    if stack_type in {"ifgram", "ion"}:
        dataset = read_subset(file, bbox, polarization=polarization)
        unwrap_key = "unw_data" if stack_type == "ifgram" else "ion_data"
        return {
            "unwrap_phase": dataset[unwrap_key],
            "coherence": dataset["cor_data"],
            "connect_component": dataset["conn_comp"],
            "pbase": pbase,
        }

    geo_ds = read_subset(file, bbox, polarization=polarization, geometry=True)
    if stack_type == "tropo":
        unwrap_phase = read_and_interpolate_troposphere(
            file, dem_file, geo_ds["xybbox"], polarization=polarization
        )
    else:
        unwrap_phase = read_and_interpolate_SET(
            file, dem_file, geo_ds["xybbox"], polarization=polarization
        )

    return {"unwrap_phase": unwrap_phase, "pbase": pbase}


# ---------------------------------------------------------------------
# Primary workflow
# ---------------------------------------------------------------------
def load_nisar(inps):
    """Prepare and load NISAR data and metadata into HDF5/MintPy format."""
    print(f"update mode: {inps.update_mode}")

    input_files = sorted(glob.glob(inps.input_glob))
    if not input_files:
        raise FileNotFoundError(
            f"No NISAR GUNW files found for input pattern: {inps.input_glob}"
        )

    if str(inps.dem_file).lower() in ["auto", "none", "no", ""]:
        raise ValueError(
            f"A real DEM path is required for prep_nisar.py; got {inps.dem_file!r}"
        )
    if not os.path.isfile(inps.dem_file):
        raise FileNotFoundError(f"DEM file not found: {inps.dem_file}")

    print(f"Found {len(input_files)} unwrapped files")

    if inps.subset_lat:
        bbox = (
            inps.subset_lon[0],
            inps.subset_lat[0],
            inps.subset_lon[1],
            inps.subset_lat[1],
        )
    else:
        bbox = None

    # extract metadata
    pol = getattr(inps, "polarization", "HH")
    metadata, bounds = extract_metadata(input_files, bbox=bbox, polarization=pol)

    # output filename
    stack_file = os.path.join(inps.out_dir, "inputs/ifgramStack.h5")
    geometry_file = os.path.join(inps.out_dir, "inputs/geometryGeo.h5")
    ion_stack_file = os.path.join(inps.out_dir, "inputs/ionStack.h5")
    tropo_stack_file = os.path.join(inps.out_dir, "inputs/tropoStack.h5")
    set_stack_file = os.path.join(inps.out_dir, "inputs/setStack.h5")

    # date pairs
    date12_list = _get_date_pairs(input_files)

    # geometry
    metadata = prepare_geometry(
        outfile=geometry_file,
        metaFile=input_files[0],
        bbox=bounds,
        metadata=metadata,
        demFile=inps.dem_file,
        maskFile=inps.mask_file,
        polarization=pol,
    )

    # standalone water mask (MintPy format)
    if getattr(inps, "mask_file", None) not in [None, "None", "auto"]:
        water_mask_file = os.path.join(inps.out_dir, "waterMask.h5")
        prepare_water_mask(
            outfile=water_mask_file,
            metaFile=input_files[0],
            metadata=metadata,
            bbox=bounds,
            maskFile=inps.mask_file,
            polarization=pol,
        )

    # ifgram stack
    prepare_stack(
        outfile=stack_file,
        inp_files=input_files,
        metadata=metadata,
        demFile=inps.dem_file,
        bbox=bounds,
        date12_list=date12_list,
        polarization=pol,
        stack_type="ifgram",
    )

    # ionosphere stack
    prepare_stack(
        outfile=ion_stack_file,
        inp_files=input_files,
        metadata=metadata,
        demFile=inps.dem_file,
        bbox=bounds,
        date12_list=date12_list,
        polarization=pol,
        stack_type="ion",
    )

    # troposphere stack
    prepare_stack(
        outfile=tropo_stack_file,
        inp_files=input_files,
        metadata=metadata,
        demFile=inps.dem_file,
        bbox=bounds,
        date12_list=date12_list,
        polarization=pol,
        stack_type="tropo",
    )
    print("Done.")

    # SET stack
    prepare_stack(
        outfile=set_stack_file,
        inp_files=input_files,
        metadata=metadata,
        demFile=inps.dem_file,
        bbox=bounds,
        date12_list=date12_list,
        polarization=pol,
        stack_type="set",
    )
    print("Done.")
    return


# ---------------------------------------------------------------------
# Metadata / subset utilities
# ---------------------------------------------------------------------
def extract_metadata(input_files, bbox=None, polarization="HH"):
    """Extract NISAR metadata for MintPy."""
    meta_file = input_files[0]
    meta = {}

    datasets = _datasets_for_pol(polarization)

    with h5py.File(meta_file, "r") as ds:
        pixel_height = ds[datasets["ySpacing"]][()]
        pixel_width = ds[datasets["xSpacing"]][()]
        x_origin = float(np.min(ds[datasets["xcoord"]][()]))
        y_origin = float(np.max(ds[datasets["ycoord"]][()]))
        xcoord = ds[datasets["xcoord"]][()]
        ycoord = ds[datasets["ycoord"]][()]
        meta["EPSG"] = int(ds[datasets["epsg"]][()])
        meta["WAVELENGTH"] = SPEED_OF_LIGHT / ds[PROCESSINFO["centerFrequency"]][()]
        meta["ORBIT_DIRECTION"] = ds[PROCESSINFO["orbit_direction"]][()].decode("utf-8")
        meta["POLARIZATION"] = polarization
        meta["ALOOKS"] = ds[datasets["azimuth_look"]][()]
        meta["RLOOKS"] = ds[datasets["range_look"]][()]
        meta["PLATFORM"] = ds[PROCESSINFO["platform"]][()].decode("utf-8")
        meta["STARTING_RANGE"] = float(
            np.min(ds[PROCESSINFO["rdr_slant_range"]][()].flatten())
        )

        start_time = datetime.datetime.strptime(
            ds[PROCESSINFO["start_time"]][()].decode("utf-8")[0:26],
            "%Y-%m-%dT%H:%M:%S.%f",
        )
        end_time = datetime.datetime.strptime(
            ds[PROCESSINFO["end_time"]][()].decode("utf-8")[0:26],
            "%Y-%m-%dT%H:%M:%S.%f",
        )

    t_mid = start_time + (end_time - start_time) / 2.0
    meta["CENTER_LINE_UTC"] = (
        t_mid - datetime.datetime(t_mid.year, t_mid.month, t_mid.day)
    ).total_seconds()

    # These were previously using //2 (integer) which is wrong for float spacing.
    meta["X_FIRST"] = x_origin - float(pixel_width) / 2.0
    meta["Y_FIRST"] = y_origin - float(pixel_height) / 2.0
    meta["X_STEP"] = float(pixel_width)
    meta["Y_STEP"] = float(pixel_height)

    if meta["EPSG"] == 4326:
        meta["X_UNIT"] = meta["Y_UNIT"] = "degree"
    else:
        meta["X_UNIT"] = meta["Y_UNIT"] = "meters"
        if str(meta["EPSG"]).startswith("326"):
            meta["UTM_ZONE"] = str(meta["EPSG"])[3:] + "N"
        else:
            meta["UTM_ZONE"] = str(meta["EPSG"])[3:] + "S"
    meta["EARTH_RADIUS"] = EARTH_RADIUS

    # NISAR altitude
    meta["HEIGHT"] = 747000

    # placeholder pixel sizes
    meta["RANGE_PIXEL_SIZE"] = abs(float(pixel_width))
    meta["AZIMUTH_PIXEL_SIZE"] = abs(float(pixel_height))

    # keep full-scene dimensions/origin first, then apply subset offsets below
    meta["LENGTH"] = int(ycoord.size)
    meta["WIDTH"] = int(xcoord.size)

    # bbox handling
    if bbox:
        epsg_src = 4326
        utm_bbox = bbox_to_utm(bbox, meta["EPSG"], epsg_src)
    else:
        utm_bbox = None

    bounds = common_raster_bound(input_files, utm_bbox, polarization=polarization)
    meta["bbox"] = ",".join([str(b) for b in bounds])

    col1, row1, col2, row2 = get_rows_cols(xcoord, ycoord, bounds)
    meta = attr.update_attribute4subset(meta, (col1, row1, col2, row2), print_msg=False)
    meta = _coerce_subset_metadata_types(meta)

    return meta, bounds


def get_rows_cols(xcoord, ycoord, bounds):
    """
    Get (col1, row1, col2, row2) for subsetting given bounds=(xmin,ymin,xmax,ymax).

    Robust to:
      - bounds slightly outside coordinate extent
      - collapsed/invalid overlap bounds
      - empty index selections
    Returns indices suitable for Python slicing: arr[row1:row2, col1:col2]
    """
    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)

    xmin, ymin, xmax, ymax = map(float, bounds)

    # Ensure bounds ordering
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    # Clip bounds to coordinate extent
    x_minc, x_maxc = float(np.nanmin(xcoord)), float(np.nanmax(xcoord))
    y_minc, y_maxc = float(np.nanmin(ycoord)), float(np.nanmax(ycoord))
    xmin = max(xmin, x_minc)
    xmax = min(xmax, x_maxc)
    ymin = max(ymin, y_minc)
    ymax = min(ymax, y_maxc)

    # If bounds collapse after clipping, fall back to full extent
    if xmin >= xmax or ymin >= ymax:
        xmin, xmax = x_minc, x_maxc
        ymin, ymax = y_minc, y_maxc

    # Primary selection
    xindex = np.where((xcoord >= xmin) & (xcoord <= xmax))[0]
    yindex = np.where((ycoord >= ymin) & (ycoord <= ymax))[0]

    # If empty (can happen due to floating roundoff), fall back to nearest
    if xindex.size == 0:
        i1 = int(np.nanargmin(np.abs(xcoord - xmin)))
        i2 = int(np.nanargmin(np.abs(xcoord - xmax)))
        col1, col2 = (i1, i2) if i1 <= i2 else (i2, i1)
    else:
        col1, col2 = int(np.min(xindex)), int(np.max(xindex))

    if yindex.size == 0:
        j1 = int(np.nanargmin(np.abs(ycoord - ymin)))
        j2 = int(np.nanargmin(np.abs(ycoord - ymax)))
        row1, row2 = (j1, j2) if j1 <= j2 else (j2, j1)
    else:
        row1, row2 = int(np.min(yindex)), int(np.max(yindex))

    # Make slice end exclusive (+1)
    col2 = min(col2 + 1, xcoord.size)
    row2 = min(row2 + 1, ycoord.size)

    return col1, row1, col2, row2


def get_raster_corners(input_file, polarization="HH"):
    """Get the (west, south, east, north) bounds of the image."""
    datasets = _datasets_for_pol(polarization)
    with h5py.File(input_file, "r") as ds:
        xcoord = ds[datasets["xcoord"]][:]
        ycoord = ds[datasets["ycoord"]][:]
        west = max(np.min(ds[PROCESSINFO["rdr_xcoord"]][:]), np.min(xcoord))
        east = min(np.max(ds[PROCESSINFO["rdr_xcoord"]][:]), np.max(xcoord))
        north = min(np.max(ds[PROCESSINFO["rdr_ycoord"]][:]), np.max(ycoord))
        south = max(np.min(ds[PROCESSINFO["rdr_ycoord"]][:]), np.min(ycoord))
    return float(west), float(south), float(east), float(north)


def common_raster_bound(input_files, utm_bbox=None, polarization="HH"):
    """Get common bounds among all data in (xmin, ymin, xmax, ymax)."""
    wests = []
    souths = []
    easts = []
    norths = []
    for file in input_files:
        west, south, east, north = get_raster_corners(file, polarization=polarization)
        wests.append(west)
        souths.append(south)
        easts.append(east)
        norths.append(north)

    common = [
        max(wests),
        max(souths),
        min(easts),
        min(norths),
    ]

    if utm_bbox:
        common = [
            max(common[0], utm_bbox[0]),
            max(common[1], utm_bbox[1]),
            min(common[2], utm_bbox[2]),
            min(common[3], utm_bbox[3]),
        ]

    if common[0] >= common[2] or common[1] >= common[3]:
        raise ValueError(
            f"No common overlap found among input files within bounds: {common}"
        )

    return common


def bbox_to_utm(bbox, dst_epsg, src_epsg=4326):
    """Convert a bounding box into the destination CRS.

    Use transform_bounds instead of projecting only two diagonal corners.
    For projected grids such as UTM, a lat/lon-aligned box is not guaranteed
    to remain axis-aligned after reprojection, so the diagonal-corner approach
    can clip valid data near the other two corners.
    """
    xmin, xmax = sorted((float(bbox[0]), float(bbox[2])))
    ymin, ymax = sorted((float(bbox[1]), float(bbox[3])))

    if int(dst_epsg) == int(src_epsg):
        return (xmin, ymin, xmax, ymax)

    transformer = Transformer.from_crs(
        f"EPSG:{src_epsg}", f"EPSG:{dst_epsg}", always_xy=True
    )
    return transformer.transform_bounds(xmin, ymin, xmax, ymax, densify_pts=21)


def read_subset(gunw_file, bbox, polarization="HH", geometry=False):
    """
    Read subset for unwrapped interferogram products.
    If geometry=True, returns bbox indices only (xybbox) without reading data arrays.
    """
    datasets = _datasets_for_pol(polarization)
    with h5py.File(gunw_file, "r") as ds:
        xcoord = ds[datasets["xcoord"]][()]
        ycoord = ds[datasets["ycoord"]][()]
        if bbox:
            col1, row1, col2, row2 = get_rows_cols(xcoord, ycoord, bbox)
        else:
            row1, row2 = 0, len(ycoord)
            col1, col2 = 0, len(xcoord)

        xybbox = (col1, row1, col2, row2)

        if geometry:
            return {"xybbox": xybbox}

        unw_dset = ds[datasets["unw"]]
        cor_dset = ds[datasets["cor"]]
        cc_dset = ds[datasets["connComp"]]
        ion_dset = ds[datasets["ion"]]

        unw_data = unw_dset[row1:row2, col1:col2].astype(np.float32)
        cor_data = cor_dset[row1:row2, col1:col2].astype(np.float32)
        conn_comp = cc_dset[row1:row2, col1:col2].astype(np.float32)
        ion_data = ion_dset[row1:row2, col1:col2].astype(np.float32)

        # fill handling
        fill_unw = unw_dset.attrs.get("_FillValue", None)
        fill_cor = cor_dset.attrs.get("_FillValue", None)
        fill_cc = cc_dset.attrs.get("_FillValue", None)
        fill_ion = ion_dset.attrs.get("_FillValue", None)

    if fill_unw is not None:
        unw_data[unw_data == fill_unw] = np.nan
    if fill_cor is not None:
        cor_data[cor_data == fill_cor] = np.nan
    if fill_cc is not None:
        conn_comp[conn_comp == fill_cc] = np.nan
    if fill_ion is not None:
        ion_data[ion_data == fill_ion] = np.nan

    return {
        "unw_data": unw_data,
        "cor_data": cor_data,
        "conn_comp": conn_comp,
        "ion_data": ion_data,
        "xybbox": xybbox,
    }


# ---------------------------------------------------------------------
# Geometry (DEM warp + 3D interpolation at valid pixels)
# ---------------------------------------------------------------------
def read_and_interpolate_geometry(
    gunw_file, dem_file, xybbox, polarization="HH", mask_file=None
):
    """
    Warp DEM to the interferogram grid (aligned), then interpolate slant range & incidence.
    Interpolation is evaluated at valid pixels only (validity from unwrappedPhase finite + _FillValue).
    """
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        {
            "slant_range": "rdr_slant_range",
            "incidence_angle": "rdr_incidence",
            "los_x": "rdr_los_x",
            "los_y": "rdr_los_y",
        },
    )
    slant_range, incidence_angle, azimuth_angle = interpolate_geometry(
        interp_ctx["x_2d"],
        interp_ctx["y_2d"],
        interp_ctx["dem"],
        interp_ctx["rdr_coords"],
        interp_ctx["valid_mask"],
    )

    # Mask handling (optional external mask warped to grid; otherwise ones)
    if mask_file in ["auto", "None", None, "no", ""]:
        mask_subset_array = np.ones(interp_ctx["dem"].shape, dtype="byte")
    else:
        mask_src_epsg = _read_raster_epsg(mask_file)
        mask_subset_array = _warp_to_grid_mem(
            src_path=mask_file,
            src_epsg=mask_src_epsg,
            dst_epsg=interp_ctx["dst_epsg"],
            xcoord=interp_ctx["xcoord"],
            ycoord=interp_ctx["ycoord"],
            resample_alg="near",
        ).astype("byte")

    return (
        interp_ctx["dem"],
        slant_range,
        incidence_angle,
        azimuth_angle,
        mask_subset_array,
    )


def interpolate_geometry(X_2d, Y_2d, dem, rdr_coords, valid_mask):
    """Interpolate slant range, incidence angle, and azimuth angle at valid pixels only."""
    shape, ii, jj, pts = _prepare_valid_interp_points(X_2d, Y_2d, dem, valid_mask)
    out_slant = _empty_interp_array(shape)
    out_incid = _empty_interp_array(shape)
    out_az = _empty_interp_array(shape)

    if pts is None:
        return out_slant, out_incid, out_az

    sl = _interpolate_radar_grid_field(rdr_coords, "slant_range", pts)
    inc = _interpolate_radar_grid_field(rdr_coords, "incidence_angle", pts)
    losx = _interpolate_radar_grid_field(rdr_coords, "los_x", pts)
    losy = _interpolate_radar_grid_field(rdr_coords, "los_y", pts)

    # Azimuth angle from horizontal LOS unit vector components.
    az = np.degrees(np.arctan2(-losy, -losx))

    out_slant[ii, jj] = sl.astype(np.float32)
    out_incid[ii, jj] = inc.astype(np.float32)
    out_az[ii, jj] = az.astype(np.float32)
    return out_slant, out_incid, out_az


def read_and_interpolate_troposphere(gunw_file, dem_file, xybbox, polarization="HH"):
    """Warp DEM to aligned grid and interpolate combined tropo at valid pixels only."""
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        {
            "wet_tropo": "rdr_wet_tropo",
            "hydrostatic_tropo": "rdr_hs_tropo",
        },
    )
    total_tropo = interpolate_troposphere(
        interp_ctx["x_2d"],
        interp_ctx["y_2d"],
        interp_ctx["dem"],
        interp_ctx["rdr_coords"],
        interp_ctx["valid_mask"],
    )
    return total_tropo


def interpolate_troposphere(X_2d, Y_2d, dem, rdr_coords, valid_mask):
    """Interpolate total tropo (hydrostatic + wet) at valid pixels only."""
    shape, ii, jj, pts = _prepare_valid_interp_points(X_2d, Y_2d, dem, valid_mask)
    out = _empty_interp_array(shape)

    if pts is None:
        return out

    rdr_coords = dict(rdr_coords)
    rdr_coords["total_tropo"] = (
        rdr_coords["hydrostatic_tropo"] + rdr_coords["wet_tropo"]
    )
    val = _interpolate_radar_grid_field(rdr_coords, "total_tropo", pts)
    out[ii, jj] = val.astype(np.float32)
    return out


def read_and_interpolate_SET(gunw_file, dem_file, xybbox, polarization="HH"):
    """Warp DEM to aligned grid and interpolate SET phase at valid pixels only."""
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        {"rdr_SET": "rdr_SET"},
    )
    set_phase = interpolate_set(
        interp_ctx["x_2d"],
        interp_ctx["y_2d"],
        interp_ctx["dem"],
        interp_ctx["rdr_coords"],
        interp_ctx["valid_mask"],
    )
    return set_phase


def interpolate_set(X_2d, Y_2d, dem, rdr_coords, valid_mask):
    """Interpolate SET phase at valid pixels only."""
    shape, ii, jj, pts = _prepare_valid_interp_points(X_2d, Y_2d, dem, valid_mask)
    out = _empty_interp_array(shape)

    if pts is None:
        return out

    val = _interpolate_radar_grid_field(rdr_coords, "rdr_SET", pts)
    out[ii, jj] = val.astype(np.float32)
    return out


# ---------------------------------------------------------------------
# MintPy file builders
# ---------------------------------------------------------------------
def _get_date_pairs(filenames):
    str_list = [Path(f).stem for f in filenames]
    return [
        str(f.split("_")[11].split("T")[0]) + "_" + str(f.split("_")[13].split("T")[0])
        for f in str_list
    ]


def prepare_geometry(
    outfile, metaFile, metadata, bbox, demFile, maskFile, polarization="HH"
):
    """Prepare the geometry file."""
    print("-" * 50)
    print(f"preparing geometry file: {outfile}")

    meta = {key: value for key, value in metadata.items()}

    geo_ds = read_subset(metaFile, bbox, polarization=polarization, geometry=True)
    dem_subset_array, slant_range, incidence_angle, azimuth_angle, mask = (
        read_and_interpolate_geometry(
            metaFile,
            demFile,
            geo_ds["xybbox"],
            polarization=polarization,
            mask_file=maskFile,
        )
    )

    length, width = dem_subset_array.shape
    ds_name_dict = {
        "height": [np.float32, (length, width), dem_subset_array],
        "incidenceAngle": [np.float32, (length, width), incidence_angle],
        "slantRangeDistance": [np.float32, (length, width), slant_range],
        "azimuthAngle": [np.float32, (length, width), azimuth_angle],
    }
    if maskFile:
        valid = _read_valid_unw_mask(metaFile, geo_ds["xybbox"], polarization)
        ds_name_dict["waterMask"] = [
            np.bool_,
            (length, width),
            mask.astype(bool) & valid,
        ]

    meta["FILE_TYPE"] = "geometry"
    meta["STARTING_RANGE"] = float(np.nanmin(slant_range))
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)
    return meta


def prepare_water_mask(outfile, metaFile, metadata, bbox, maskFile, polarization="HH"):
    """Prepare a standalone MintPy waterMask.h5 aligned to the NISAR grid."""
    print("-" * 50)
    print(f"preparing water mask file: {outfile}")

    if not maskFile or maskFile in ["auto", "None", None]:
        raise ValueError("maskFile must be a raster path (e.g., waterMask.msk)")

    meta = {key: value for key, value in metadata.items()}

    # get subset indices
    geo_ds = read_subset(metaFile, bbox, polarization=polarization, geometry=True)
    xybbox = geo_ds["xybbox"]

    # get target grid axes + EPSG from the NISAR file
    dst_epsg, xcoord, ycoord = _read_target_grid(metaFile, xybbox, polarization)

    # warp mask raster onto NISAR grid
    mask_src_epsg = _read_raster_epsg(maskFile)
    mask_arr = _warp_to_grid_mem(
        src_path=maskFile,
        src_epsg=mask_src_epsg,
        dst_epsg=dst_epsg,
        xcoord=xcoord,
        ycoord=ycoord,
        resample_alg="near",
    ).astype("byte")

    # Convention in this script: nonzero => True (valid), 0 => False (masked)
    water_mask_bool = mask_arr.astype(bool)

    # constrain to valid NISAR pixels (finite/unfilled unwrappedPhase)
    valid = _read_valid_unw_mask(metaFile, xybbox, polarization)
    water_mask_bool &= valid

    length, width = water_mask_bool.shape
    ds_name_dict = {"waterMask": [np.bool_, (length, width), water_mask_bool]}

    meta["FILE_TYPE"] = "waterMask"
    meta["LENGTH"] = int(length)
    meta["WIDTH"] = int(width)

    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)
    return outfile


def prepare_stack(
    outfile,
    inp_files,
    metadata,
    demFile,
    bbox,
    date12_list,
    polarization="HH",
    stack_type=None,
):
    """Prepare the input stacks."""
    effective_stack_type = _resolve_stack_type(stack_type, outfile)
    print("-" * 50)
    print(f"preparing {effective_stack_type} stack file: {outfile}")

    meta = {key: value for key, value in metadata.items()}
    num_pair = len(inp_files)
    print(f"number of inputs/unwrapped interferograms: {num_pair}")

    pbase = np.zeros(num_pair, dtype=np.float32)
    cols = int(meta["WIDTH"])
    rows = int(meta["LENGTH"])

    date12_arr = np.array([x.split("_") for x in date12_list], dtype=np.bytes_)
    drop_ifgram = np.ones(num_pair, dtype=np.bool_)

    ds_name_dict = {
        "date": [date12_arr.dtype, (num_pair, 2), date12_arr],
        "bperp": [np.float32, (num_pair,), pbase],
        "dropIfgram": [np.bool_, (num_pair,), drop_ifgram],
        "unwrapPhase": [np.float32, (num_pair, rows, cols), None],
        "coherence": [np.float32, (num_pair, rows, cols), None],
        "connectComponent": [np.float32, (num_pair, rows, cols), None],
    }

    meta["FILE_TYPE"] = "ifgramStack"

    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    print(f"writing data to HDF5 file {outfile} with a mode ...")
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_pair)
        for i, file in enumerate(inp_files):
            obs = _read_stack_observation(
                file, effective_stack_type, bbox, demFile, polarization
            )
            f["unwrapPhase"][i] = obs["unwrap_phase"]

            if "coherence" in obs:
                f["coherence"][i] = obs["coherence"]
                f["connectComponent"][i] = obs["connect_component"]

            if "pbase" in obs:
                f["bperp"][i] = obs["pbase"]

            prog_bar.update(i + 1, suffix=date12_list[i])
        prog_bar.close()

    print(f"finished writing to HDF5 file: {outfile}")
    return outfile
