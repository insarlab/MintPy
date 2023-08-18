"""Utilities to make HDF5 files CF-compliant and readable by GDAL."""
from __future__ import annotations

import datetime

import h5py
import numpy as np
import pyproj

from mintpy.utils import readfile

__all__ = ["write_coordinate_system"]


def write_coordinate_system(
    filename: str,
    dset_name: str,
    xy_dim_names: tuple[str, str] = ("x", "y"),
    grid_mapping_dset: str = "spatial_ref",
) -> None:
    """
    Write the coordinate system CF metadata to an existing HDF5 file.

    Parameters
    ----------
    filename : str
        File path.
    dset_name : str
        Dataset name within the HDF5 file.
    xy_dim_names : tuple[str, str], optional
        x and y dimension names, by default ("x", "y").
    grid_mapping_dset : str, optional
        Dataset name for grid mapping HDF5 attribute
        By default "spatial_ref" (matching `rioxarray`).
    """
    atr = readfile.read_attribute(filename)
    epsg = int(atr.get("EPSG", 4326))

    with h5py.File(filename, "a") as hf:
        crs = pyproj.CRS.from_epsg(epsg)
        dset = hf[dset_name]

        # Setup the dataset holding the SRS information
        srs_dset = hf.require_dataset(grid_mapping_dset, shape=(), dtype=int)
        srs_dset.attrs.update(crs.to_cf())
        dset.attrs["grid_mapping"] = grid_mapping_dset

        _setup_time_dimension(hf, dset)
        _setup_xy_dimensions(hf, dset, atr, crs, xy_dim_names)


def _get_xy_arrays(atr: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate x and y arrays from attribute dictionary.

    Parameters
    ----------
    atr : dict
        MintPy attribute dictionary containing ROIPAC raster attributes.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        x and y arrays.
    """
    x0 = float(atr["X_FIRST"])
    y0 = float(atr["Y_FIRST"])
    x_step = float(atr["X_STEP"])
    y_step = float(atr["Y_STEP"])
    rows = int(atr["LENGTH"])
    cols = int(atr["WIDTH"])

    x_arr = x0 + x_step * np.arange(cols)
    y_arr = y0 + y_step * np.arange(rows)

    # Shift by half pixel to get the centers
    x_arr += x_step / 2
    y_arr += y_step / 2
    return x_arr, y_arr


def _get_coordinate_metadata(crs: pyproj.CRS) -> tuple[dict, dict]:
    """Get HDF5 coordinate metadata based on CRS type.

    Parameters
    ----------
    crs : pyproj.CRS
        Coordinate Reference System.

    Returns
    -------
    tuple[dict, dict]
        Dictionary of metadata for x and y respectively.
    """
    x_coord_attrs = {"axis": "X"}
    y_coord_attrs = {"axis": "Y"}

    if crs.is_projected:
        units = "meters"
        # X metadata
        x_coord_attrs.update(
            {
                "long_name": "x coordinate of projection",
                "standard_name": "projection_x_coordinate",
                "units": units,
            }
        )
        # Y metadata
        y_coord_attrs.update(
            {
                "long_name": "y coordinate of projection",
                "standard_name": "projection_y_coordinate",
                "units": units,
            }
        )
    elif crs.is_geographic:
        # X metadata
        x_coord_attrs.update(
            {
                "long_name": "longitude",
                "standard_name": "longitude",
                "units": "degrees_east",
            }
        )
        # Y metadata
        y_coord_attrs.update(
            {
                "long_name": "latitude",
                "standard_name": "latitude",
                "units": "degrees_north",
            }
        )
    return x_coord_attrs, y_coord_attrs


def _setup_xy_dimensions(
    hf: h5py.File,
    dset: h5py.Dataset,
    atr,
    crs: pyproj.CRS,
    xy_dim_names: tuple[str, str] = ("x", "y"),
):
    """
    Setup time dimension in the HDF5 file.

    Parameters
    ----------
    hf : h5py.File
        HDF5 file object.
    dset : h5py.Dataset
        Dataset within the HDF5 file.
    atr : dict
        MintPy attribute dictionary containing ROIPAC raster attributes.
    xy_dim_names : tuple[str, str], optional
        x and y dimension names, by default ("x", "y").

    Returns
    -------
    Optional[list[datetime.datetime]]
        list of dates if they exist, otherwise None.
    """
    x_dim_name, y_dim_name = xy_dim_names
    # add metadata to x, y coordinates
    x_arr, y_arr = _get_xy_arrays(atr)
    x_dim_dset = hf.create_dataset(x_dim_name, data=x_arr)
    x_dim_dset.make_scale(x_dim_name)
    y_dim_dset = hf.create_dataset(y_dim_name, data=y_arr)
    y_dim_dset.make_scale(y_dim_name)

    x_coord_attrs, y_coord_attrs = _get_coordinate_metadata(crs)

    y_dim_dset.attrs.update(y_coord_attrs)
    x_dim_dset.attrs.update(x_coord_attrs)

    ndim = dset.ndim
    dset.dims[ndim - 1].attach_scale(x_dim_dset)
    dset.dims[ndim - 2].attach_scale(y_dim_dset)
    dset.dims[ndim - 1].label = x_dim_name
    dset.dims[ndim - 2].label = y_dim_name


def _setup_time_dimension(hf: h5py.File, dset: h5py.Dataset):
    """Setup time dimension in the HDF5 file.

    Parameters
    ----------
    hf : h5py.File
        HDF5 file object.
    dset : h5py.Dataset
        Dataset within the HDF5 file.
    """
    if "date" not in hf:
        # If we want to do something other than time as a 3rd dimension
        # (e.g. for ifg date pairs) we'll need to figure out what other valid
        # dims there are. otherwise, you can just do `phony_dims="sort"` in xarray
        return None

    date_arr = [
        datetime.datetime.strptime(ds, "%Y%m%d") for ds in hf["date"][()].astype(str)
    ]
    days_since = [(d - date_arr[0]).days for d in date_arr]
    dt_dim = hf.create_dataset("time", data=days_since)
    dt_dim.make_scale()
    cf_attrs = dict(
        units=f"days since {str(date_arr[0])}", calendar="proleptic_gregorian"
    )
    dt_dim.attrs.update(cf_attrs)
    dset.dims[0].attach_scale(dt_dim)
    dset.dims[0].label = "time"
