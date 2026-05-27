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
# Constants / HDF5 paths
# ---------------------------------------------------------------------
FREQUENCY_MAP = {
    "A": "frequencyA",
    "B": "frequencyB",
    "frequencyA": "frequencyA",
    "frequencyB": "frequencyB",
}
IDENTIFICATION = "/science/LSAR/identification"
RADARGRID_ROOT = "science/LSAR/GUNW/metadata/radarGrid"

PROCESSINFO = {
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
def _normalize_frequency(frequency) -> str:
    """Return the GUNW frequency group name for CLI values auto/A/B."""
    if frequency is None or str(frequency).lower() == "auto":
        return "frequencyA"

    normalized = FREQUENCY_MAP.get(str(frequency))
    if normalized is None:
        raise ValueError("frequency must be one of: auto, A, B")
    return normalized


def _dataset_root_unw(frequency: str) -> str:
    return f"/science/LSAR/GUNW/grids/{frequency}/unwrappedInterferogram"


def _parameters_root(frequency: str) -> str:
    return (
        "/science/LSAR/GUNW/metadata/processingInformation/parameters/"
        f"unwrappedInterferogram/{frequency}"
    )


def _center_frequency_path(frequency: str) -> str:
    return f"/science/LSAR/GUNW/grids/{frequency}/centerFrequency"


def _datasets_for_pol(polarization: str, frequency: str) -> dict:
    """Return per-call dataset paths for the selected frequency/polarization."""
    root = _dataset_root_unw(frequency)
    parameters = _parameters_root(frequency)
    return {
        "xcoord": f"{root}/{polarization}/xCoordinates",
        "ycoord": f"{root}/{polarization}/yCoordinates",
        "unw": f"{root}/{polarization}/unwrappedPhase",
        "mask": f"{root}/mask",
        "cor": f"{root}/{polarization}/coherenceMagnitude",
        "connComp": f"{root}/{polarization}/connectedComponents",
        "ion": f"{root}/{polarization}/ionospherePhaseScreen",
        "epsg": f"{root}/{polarization}/projection",
        "xSpacing": f"{root}/{polarization}/xCoordinateSpacing",
        "ySpacing": f"{root}/{polarization}/yCoordinateSpacing",
        "polarization": f"/science/LSAR/GUNW/grids/{frequency}/listOfPolarizations",
        "range_look": f"{parameters}/numberOfRangeLooks",
        "azimuth_look": f"{parameters}/numberOfAzimuthLooks",
    }


def _resolve_frequency(gunw_file: str, frequency, polarization: str) -> str:
    """Resolve and validate the requested NISAR frequency."""
    resolved = _normalize_frequency(frequency)
    datasets = _datasets_for_pol(polarization, resolved)

    with h5py.File(gunw_file, "r") as ds:
        required_paths = [
            _dataset_root_unw(resolved),
            datasets["unw"],
            _center_frequency_path(resolved),
        ]
        missing = [path for path in required_paths if path not in ds]

    if missing:
        requested_frequency = "auto" if frequency is None else str(frequency)
        requested = (
            "auto (frequencyA)"
            if requested_frequency.lower() == "auto"
            else requested_frequency
        )
        if requested_frequency in ["auto", "A", "frequencyA"]:
            hint = "Use --frequency B for frequencyB products."
        else:
            hint = (
                "Check that the input file contains frequencyB for this polarization."
            )
        raise ValueError(
            f"NISAR {requested} data for polarization {polarization!r} was not found "
            f"in {gunw_file}. Missing path: {missing[0]}. {hint}"
        )

    return resolved


def _grid_bounds_from_xy(xcoord: np.ndarray, ycoord: np.ndarray):
    """Return pixel-edge bounds and signed spacing for x/y pixel centers."""
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
    """Warp a raster to the exact xcoord/ycoord grid using MEM output."""
    bounds, _, _ = _grid_bounds_from_xy(xcoord, ycoord)

    warp_opts = gdal.WarpOptions(
        format="MEM",
        outputBounds=bounds,
        srcSRS=f"EPSG:{src_epsg}",
        dstSRS=f"EPSG:{dst_epsg}",
        width=int(xcoord.size),
        height=int(ycoord.size),
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
    """Wrap scipy interpolation with axis-order and out-of-bounds handling."""
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


def _read_unwrapped_phase_valid_mask(gunw_file: str, xybbox, pol: str, frequency: str):
    """Fallback validity mask based on finite unwrappedPhase (+ _FillValue check)."""
    datasets = _datasets_for_pol(pol, frequency)
    path = datasets["unw"]
    with h5py.File(gunw_file, "r") as ds:
        dset = ds[path]
        unw = dset[xybbox[1] : xybbox[3], xybbox[0] : xybbox[2]]
        fill = dset.attrs.get("_FillValue", None)

    valid = np.isfinite(unw)
    if fill is not None:
        valid &= unw != fill
    return valid


def _read_is_land_and_valid_mask(gunw_file: str, xybbox, pol: str, frequency: str):
    """Decode the native GUNW mask into MintPy's keep-mask convention."""
    datasets = _datasets_for_pol(pol, frequency)
    path = datasets["mask"]

    with h5py.File(gunw_file, "r") as ds:
        if path not in ds:
            return _read_unwrapped_phase_valid_mask(gunw_file, xybbox, pol, frequency)

        dset = ds[path]
        mask_bits = dset[xybbox[1] : xybbox[3], xybbox[0] : xybbox[2]]
        fill = dset.attrs.get("_FillValue", None)

    valid_samples = np.isfinite(mask_bits)
    if fill is not None:
        valid_samples &= mask_bits != fill

    # Bits 0-7: subswath and water encoding; higher bits carry other flags.
    mask = np.where(valid_samples, mask_bits & 0xFF, 0).astype(np.int64, copy=False)
    water_mask = (mask // 100) == 1
    ref_subswath = (mask // 10) % 10
    sec_subswath = mask % 10
    is_valid = valid_samples & (ref_subswath > 0) & (sec_subswath > 0)
    return is_valid & ~water_mask


def _read_common_is_land_and_valid_mask(input_files, bbox, pol: str, frequency: str):
    """Return the common keep-mask across all input GUNW products."""
    common_mask = None
    for gunw_file in input_files:
        geo_ds = read_subset(
            gunw_file, bbox, polarization=pol, frequency=frequency, geometry=True
        )
        mask = _read_is_land_and_valid_mask(
            gunw_file, geo_ds["xybbox"], pol, frequency
        ).astype(bool, copy=False)

        if common_mask is None:
            common_mask = mask.copy()
            continue

        if mask.shape != common_mask.shape:
            raise ValueError(
                "Cannot combine NISAR masks with different subset shapes: "
                f"{common_mask.shape} versus {mask.shape} in {gunw_file}"
            )
        common_mask &= mask

    if common_mask is None:
        raise ValueError("No NISAR input files found for common mask generation.")

    return common_mask


def _external_mask_is_set(external_mask_file):
    return external_mask_file not in ["auto", "None", None, "no", ""]


def _apply_external_mask(
    keep_mask,
    gunw_file: str,
    xybbox,
    external_mask_file,
    polarization: str,
    frequency: str,
):
    """Refine a keep-mask with an optional external raster mask."""
    if not _external_mask_is_set(external_mask_file):
        return keep_mask

    dst_epsg, xcoord, ycoord = _read_target_grid(
        gunw_file, xybbox, polarization, frequency
    )
    mask_src_epsg = _read_raster_epsg(external_mask_file)
    external_mask = _warp_to_grid_mem(
        src_path=external_mask_file,
        src_epsg=mask_src_epsg,
        dst_epsg=dst_epsg,
        xcoord=xcoord,
        ycoord=ycoord,
        resample_alg="near",
    ).astype(bool)
    return keep_mask & external_mask


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


def _read_target_grid(gunw_file: str, xybbox, polarization: str, frequency: str):
    """Read the destination EPSG and subset grid axes from a GUNW file."""
    datasets = _datasets_for_pol(polarization, frequency)
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
    gunw_file, dem_file, xybbox, polarization, frequency, field_map, valid_mask=None
):
    """Build the common DEM/grid/valid-mask context for radar-grid interpolation."""
    dem_src_epsg = _read_raster_epsg(dem_file)
    dst_epsg, xcoord, ycoord = _read_target_grid(
        gunw_file, xybbox, polarization, frequency
    )
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
    if valid_mask is None:
        valid_mask = _read_is_land_and_valid_mask(
            gunw_file, xybbox, polarization, frequency
        )
    else:
        valid_mask = np.asarray(valid_mask, dtype=np.bool_)
        if valid_mask.shape != x_2d.shape:
            raise ValueError(
                "Provided NISAR valid mask shape does not match the target grid: "
                f"{valid_mask.shape} versus {x_2d.shape}"
            )

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


def _required_paths_for_stack_type(stack_type, polarization, frequency):
    """Return HDF5 source datasets needed to build the requested stack."""
    datasets = _datasets_for_pol(polarization, frequency)
    if stack_type == "ifgram":
        return [datasets["unw"], datasets["cor"], datasets["connComp"]]
    if stack_type == "ion":
        return [datasets["ion"], datasets["cor"], datasets["connComp"]]
    if stack_type == "tropo":
        return [PROCESSINFO["rdr_wet_tropo"], PROCESSINFO["rdr_hs_tropo"]]
    if stack_type == "set":
        return [PROCESSINFO["rdr_SET"]]

    raise ValueError(
        f"Unsupported stack_type {stack_type!r}; expected one of {sorted(STACK_TYPES)}"
    )


def _missing_required_paths(inp_files, stack_type, polarization, frequency):
    """Return missing required HDF5 source paths as (file, path) pairs."""
    required_paths = _required_paths_for_stack_type(stack_type, polarization, frequency)
    missing = []

    for file in inp_files:
        with h5py.File(file, "r") as ds:
            missing.extend((file, path) for path in required_paths if path not in ds)

    return missing


def _read_stack_observation(file, stack_type, bbox, dem_file, polarization, frequency):
    """Read one observation for the requested stack type."""
    pbase = _read_perpendicular_baseline(file)

    if stack_type in {"ifgram", "ion"}:
        dataset = read_subset(
            file, bbox, polarization=polarization, frequency=frequency
        )
        unwrap_key = "unw_data" if stack_type == "ifgram" else "ion_data"
        return {
            "unwrap_phase": dataset[unwrap_key],
            "coherence": dataset["cor_data"],
            "connect_component": dataset["conn_comp"],
            "pbase": pbase,
        }

    geo_ds = read_subset(
        file, bbox, polarization=polarization, frequency=frequency, geometry=True
    )
    if stack_type == "tropo":
        unwrap_phase = read_and_interpolate_troposphere(
            file,
            dem_file,
            geo_ds["xybbox"],
            polarization=polarization,
            frequency=frequency,
        )
    else:
        unwrap_phase = read_and_interpolate_SET(
            file,
            dem_file,
            geo_ds["xybbox"],
            polarization=polarization,
            frequency=frequency,
        )

    return {"unwrap_phase": unwrap_phase, "pbase": pbase}


def _apply_common_mask_to_observation(obs, common_mask):
    """Mask all 2D observation layers with the stack-wide common keep-mask."""
    if common_mask is None:
        return obs

    common_mask = np.asarray(common_mask, dtype=np.bool_)
    invalid = ~common_mask
    masked = {}
    for key, value in obs.items():
        if not isinstance(value, np.ndarray):
            masked[key] = value
            continue

        if value.shape != common_mask.shape:
            raise ValueError(
                f"Cannot apply NISAR common mask with shape {common_mask.shape} "
                f"to {key} layer with shape {value.shape}"
            )

        data = value.copy()
        if np.issubdtype(data.dtype, np.floating):
            data[invalid] = np.nan
        else:
            data[invalid] = 0
        masked[key] = data

    return masked


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
    frequency = _resolve_frequency(
        input_files[0], getattr(inps, "frequency", "auto"), pol
    )
    print(f"Using NISAR {frequency}")
    metadata, bounds = extract_metadata(
        input_files, bbox=bbox, polarization=pol, frequency=frequency
    )
    common_mask = _read_common_is_land_and_valid_mask(
        input_files, bounds, pol=pol, frequency=frequency
    )
    first_geo_ds = read_subset(
        input_files[0], bounds, polarization=pol, frequency=frequency, geometry=True
    )
    common_mask = _apply_external_mask(
        common_mask,
        input_files[0],
        first_geo_ds["xybbox"],
        inps.mask_file,
        pol,
        frequency,
    )
    print(
        "Common valid land pixels from all NISAR masks: "
        f"{np.sum(common_mask)} / {common_mask.size}"
    )

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
        externalMaskFile=None,
        commonMask=common_mask,
        polarization=pol,
        frequency=frequency,
    )

    # standalone water mask (MintPy format)
    water_mask_file = os.path.join(inps.out_dir, "waterMask.h5")
    prepare_water_mask(
        outfile=water_mask_file,
        metaFile=input_files[0],
        metadata=metadata,
        bbox=bounds,
        externalMaskFile=None,
        commonMask=common_mask,
        polarization=pol,
        frequency=frequency,
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
        frequency=frequency,
        stack_type="ifgram",
        commonMask=common_mask,
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
        frequency=frequency,
        stack_type="ion",
        commonMask=common_mask,
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
        frequency=frequency,
        stack_type="tropo",
        commonMask=common_mask,
    )

    # SET stack
    prepare_stack(
        outfile=set_stack_file,
        inp_files=input_files,
        metadata=metadata,
        demFile=inps.dem_file,
        bbox=bounds,
        date12_list=date12_list,
        polarization=pol,
        frequency=frequency,
        stack_type="set",
        commonMask=common_mask,
    )
    print("Done.")
    return


# ---------------------------------------------------------------------
# Metadata / subset utilities
# ---------------------------------------------------------------------
def extract_metadata(input_files, bbox=None, polarization="HH", frequency="frequencyA"):
    """Extract NISAR metadata for MintPy."""
    meta_file = input_files[0]
    meta = {}

    datasets = _datasets_for_pol(polarization, frequency)

    with h5py.File(meta_file, "r") as ds:
        pixel_height = ds[datasets["ySpacing"]][()]
        pixel_width = ds[datasets["xSpacing"]][()]
        x_origin = float(np.min(ds[datasets["xcoord"]][()]))
        y_origin = float(np.max(ds[datasets["ycoord"]][()]))
        xcoord = ds[datasets["xcoord"]][()]
        ycoord = ds[datasets["ycoord"]][()]
        meta["EPSG"] = int(ds[datasets["epsg"]][()])
        meta["WAVELENGTH"] = SPEED_OF_LIGHT / ds[_center_frequency_path(frequency)][()]
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

    bounds = common_raster_bound(
        input_files, utm_bbox, polarization=polarization, frequency=frequency
    )
    meta["bbox"] = ",".join([str(b) for b in bounds])

    col1, row1, col2, row2 = get_rows_cols(xcoord, ycoord, bounds)
    meta = attr.update_attribute4subset(meta, (col1, row1, col2, row2), print_msg=False)
    meta = _coerce_subset_metadata_types(meta)

    return meta, bounds


def get_rows_cols(xcoord, ycoord, bounds):
    """Get subset indices for bounds=(xmin, ymin, xmax, ymax)."""
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


def get_raster_corners(input_file, polarization="HH", frequency="frequencyA"):
    """Get the (west, south, east, north) bounds of the image."""
    datasets = _datasets_for_pol(polarization, frequency)
    with h5py.File(input_file, "r") as ds:
        xcoord = ds[datasets["xcoord"]][:]
        ycoord = ds[datasets["ycoord"]][:]
        west = max(np.min(ds[PROCESSINFO["rdr_xcoord"]][:]), np.min(xcoord))
        east = min(np.max(ds[PROCESSINFO["rdr_xcoord"]][:]), np.max(xcoord))
        north = min(np.max(ds[PROCESSINFO["rdr_ycoord"]][:]), np.max(ycoord))
        south = max(np.min(ds[PROCESSINFO["rdr_ycoord"]][:]), np.min(ycoord))
    return float(west), float(south), float(east), float(north)


def common_raster_bound(
    input_files, utm_bbox=None, polarization="HH", frequency="frequencyA"
):
    """Get common bounds among all data in (xmin, ymin, xmax, ymax)."""
    wests = []
    souths = []
    easts = []
    norths = []
    for file in input_files:
        west, south, east, north = get_raster_corners(
            file, polarization=polarization, frequency=frequency
        )
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
    """Convert a bounding box into the destination CRS."""
    # Use transform_bounds instead of projecting only two diagonal corners.
    # Projected lat/lon-aligned boxes are not guaranteed to remain axis-aligned.
    xmin, xmax = sorted((float(bbox[0]), float(bbox[2])))
    ymin, ymax = sorted((float(bbox[1]), float(bbox[3])))

    if int(dst_epsg) == int(src_epsg):
        return (xmin, ymin, xmax, ymax)

    transformer = Transformer.from_crs(
        f"EPSG:{src_epsg}", f"EPSG:{dst_epsg}", always_xy=True
    )
    return transformer.transform_bounds(xmin, ymin, xmax, ymax, densify_pts=21)


def read_subset(
    gunw_file, bbox, polarization="HH", frequency="frequencyA", geometry=False
):
    """Read subset arrays or only geometry bounds for unwrapped products."""
    datasets = _datasets_for_pol(polarization, frequency)
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
    gunw_file,
    dem_file,
    xybbox,
    polarization="HH",
    frequency="frequencyA",
    external_mask_file=None,
    valid_mask=None,
):
    """Warp DEM to the interferogram grid and interpolate geometry layers."""
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        frequency,
        {
            "slant_range": "rdr_slant_range",
            "incidence_angle": "rdr_incidence",
            "los_x": "rdr_los_x",
            "los_y": "rdr_los_y",
        },
        valid_mask=valid_mask,
    )
    slant_range, incidence_angle, azimuth_angle = interpolate_geometry(
        interp_ctx["x_2d"],
        interp_ctx["y_2d"],
        interp_ctx["dem"],
        interp_ctx["rdr_coords"],
        interp_ctx["valid_mask"],
    )

    # Base mask comes from the native/common GUNW mask; external masks only refine it.
    mask_subset_array = interp_ctx["valid_mask"].astype(bool)
    if _external_mask_is_set(external_mask_file):
        mask_subset_array = _apply_external_mask(
            mask_subset_array,
            gunw_file,
            xybbox,
            external_mask_file,
            polarization,
            frequency,
        )

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


def read_and_interpolate_troposphere(
    gunw_file, dem_file, xybbox, polarization="HH", frequency="frequencyA"
):
    """Warp DEM to aligned grid and interpolate combined tropo at valid pixels only."""
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        frequency,
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


def read_and_interpolate_SET(
    gunw_file, dem_file, xybbox, polarization="HH", frequency="frequencyA"
):
    """Warp DEM to aligned grid and interpolate SET phase at valid pixels only."""
    interp_ctx = _prepare_radar_grid_interpolation(
        gunw_file,
        dem_file,
        xybbox,
        polarization,
        frequency,
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
    """Return reference_secondary date pairs in YYYYMMDD_YYYYMMDD format."""
    date12_list = []
    for filename in filenames:
        with h5py.File(filename, "r") as ds:
            if (
                f"{IDENTIFICATION}/referenceZeroDopplerStartTime" in ds
                and f"{IDENTIFICATION}/secondaryZeroDopplerStartTime" in ds
            ):
                ref_time = ds[f"{IDENTIFICATION}/referenceZeroDopplerStartTime"][()]
                sec_time = ds[f"{IDENTIFICATION}/secondaryZeroDopplerStartTime"][()]
                ref_date = ref_time.decode("utf-8").split("T")[0].replace("-", "")
                sec_date = sec_time.decode("utf-8").split("T")[0].replace("-", "")
                date12_list.append(f"{ref_date}_{sec_date}")
                continue

        parts = Path(filename).stem.split("_")
        if len(parts) > 13:
            date12_list.append(f"{parts[11].split('T')[0]}_{parts[13].split('T')[0]}")
            continue

        raise ValueError(
            f"Could not determine reference/secondary dates from {filename}. "
            "Expected NISAR identification zero-Doppler start times or an "
            "OPERA-style filename."
        )

    return date12_list


def prepare_geometry(
    outfile,
    metaFile,
    metadata,
    bbox,
    demFile,
    externalMaskFile,
    polarization="HH",
    frequency="frequencyA",
    commonMask=None,
):
    """Prepare the geometry file."""
    print("-" * 50)
    print(f"preparing geometry file: {outfile}")

    meta = {key: value for key, value in metadata.items()}

    geo_ds = read_subset(
        metaFile, bbox, polarization=polarization, frequency=frequency, geometry=True
    )
    dem_subset_array, slant_range, incidence_angle, azimuth_angle, mask = (
        read_and_interpolate_geometry(
            metaFile,
            demFile,
            geo_ds["xybbox"],
            polarization=polarization,
            frequency=frequency,
            external_mask_file=externalMaskFile,
            valid_mask=commonMask,
        )
    )

    length, width = dem_subset_array.shape
    ds_name_dict = {
        "height": [np.float32, (length, width), dem_subset_array],
        "incidenceAngle": [np.float32, (length, width), incidence_angle],
        "slantRangeDistance": [np.float32, (length, width), slant_range],
        "azimuthAngle": [np.float32, (length, width), azimuth_angle],
    }
    ds_name_dict["waterMask"] = [
        np.bool_,
        (length, width),
        mask.astype(bool),
    ]

    meta["FILE_TYPE"] = "geometry"
    meta["STARTING_RANGE"] = float(np.nanmin(slant_range))
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)
    return meta


def prepare_water_mask(
    outfile,
    metaFile,
    metadata,
    bbox,
    externalMaskFile,
    polarization="HH",
    frequency="frequencyA",
    commonMask=None,
):
    """Prepare a standalone MintPy waterMask.h5 from the GUNW mask."""
    print("-" * 50)
    print(f"preparing water mask file: {outfile}")

    meta = {key: value for key, value in metadata.items()}

    # get subset indices
    geo_ds = read_subset(
        metaFile, bbox, polarization=polarization, frequency=frequency, geometry=True
    )
    xybbox = geo_ds["xybbox"]

    if commonMask is None:
        water_mask_bool = _read_is_land_and_valid_mask(
            metaFile, xybbox, polarization, frequency
        )
    else:
        water_mask_bool = np.asarray(commonMask, dtype=np.bool_)
        length = xybbox[3] - xybbox[1]
        width = xybbox[2] - xybbox[0]
        if water_mask_bool.shape != (length, width):
            raise ValueError(
                "Provided NISAR common mask shape does not match the water mask "
                f"grid: {water_mask_bool.shape} versus {(length, width)}"
            )

    water_mask_bool = _apply_external_mask(
        water_mask_bool,
        metaFile,
        xybbox,
        externalMaskFile,
        polarization,
        frequency,
    )

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
    frequency="frequencyA",
    stack_type=None,
    commonMask=None,
):
    """Prepare the input stacks."""
    effective_stack_type = _resolve_stack_type(stack_type, outfile)
    print("-" * 50)
    print(f"preparing {effective_stack_type} stack file: {outfile}")

    meta = {key: value for key, value in metadata.items()}
    num_pair = len(inp_files)
    print(f"number of inputs/unwrapped interferograms: {num_pair}")

    missing = _missing_required_paths(
        inp_files, effective_stack_type, polarization, frequency
    )
    if missing:
        first_file, first_path = missing[0]
        message = (
            f"required NISAR {effective_stack_type} layer is missing: "
            f"{first_path} in {first_file}"
        )
        if effective_stack_type == "ifgram":
            raise FileNotFoundError(message)

        print(f"WARNING: skipping {effective_stack_type} stack because {message}")
        for missing_file, missing_path in missing[1:]:
            print(
                f"WARNING: skipping {effective_stack_type} stack because "
                f"required NISAR {effective_stack_type} layer is missing: "
                f"{missing_path} in {missing_file}"
            )
        if os.path.exists(outfile):
            print(f"WARNING: existing stack file was not updated: {outfile}")
        return None

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
                file, effective_stack_type, bbox, demFile, polarization, frequency
            )
            obs = _apply_common_mask_to_observation(obs, commonMask)
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
