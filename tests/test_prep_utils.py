from datetime import date, timedelta

import h5py
import numpy as np
import pytest
from osgeo import gdal

from mintpy.utils import prep_utils, writefile

gdal.UseExceptions()


@pytest.fixture
def metadata():
    return {
        "EPSG": "32611",
        "X_FIRST": "480540.0",
        "X_STEP": "30.0",
        "X_UNIT": "meters",
        "Y_FIRST": "3902670.0",
        "Y_STEP": "-30.0",
        "Y_UNIT": "meters",
        "LENGTH": "100",
        "WIDTH": "100",
        "REF_X": "50",
        "REF_Y": "50",
        # Sample S1 metadata
        "ALOOKS": "1",
        "AZIMUTH_PIXEL_SIZE": "14.1",
        "CENTER_LINE_UTC": "49436.654675",
        "EARTH_RADIUS": "6371000.0",
        "HEIGHT": "750000.0",
        "ORBIT_DIRECTION": "Descending",
        "PLATFORM": "S1A",
        "RANGE_PIXEL_SIZE": "2.329562114715323",
        "RLOOKS": "1",
        "STARTING_RANGE": "845960.7488998839",
        "UNIT": "m",
        "WAVELENGTH": "0.05546576",
        "DATA_TYPE": "float32",
        "PROCESSOR": "mintpy",
        "NO_DATA_VALUE": "none",
    }


@pytest.fixture
def data2d(metadata):
    rows = int(metadata["LENGTH"])
    cols = int(metadata["WIDTH"])
    return np.random.randn(rows, cols).astype("float32")


@pytest.fixture
def dates():
    num_dates = 10
    date_list = [date(2020, 1, 1) + i * timedelta(days=12) for i in range(num_dates)]
    date_strs = [d.strftime("%Y%m%d") for d in date_list]
    return np.array(date_strs, dtype=np.string_)


@pytest.fixture
def data3d(dates, metadata):
    rows = int(metadata["LENGTH"])
    cols = int(metadata["WIDTH"])
    return np.random.randn(len(dates), rows, cols).astype("float32")


def _check_gdal_metadata(gdal_str, meta):
    ds = gdal.Open(gdal_str)
    assert ds is not None
    assert int(ds.GetSpatialRef().GetAuthorityCode(None)) == int(meta["EPSG"])

    gt = ds.GetGeoTransform()
    # GEOTRANSFORM = [204500.0, 5.0, 0.0, 2151300.0, 0.0, -10.0]
    assert gt[0] == float(meta["X_FIRST"])
    assert gt[3] == float(meta["Y_FIRST"])
    assert gt[1] == float(meta["X_STEP"])
    assert gt[5] == float(meta["Y_STEP"])
    ds = None


def test_2d(metadata, data2d, tmp_path):
    rows, cols = data2d.shape
    ds_name_dict = {
        "velocity": [np.float32, (rows, cols), data2d],
    }

    # initiate HDF5 file
    metadata["FILE_TYPE"] = "velocity"
    # meta["REF_DATE"] = ref_date # might not be the first date!
    outfile = tmp_path / "velocity.h5"
    writefile.layout_hdf5(str(outfile), ds_name_dict, metadata=metadata)

    # Check the metadata was written
    with h5py.File(outfile) as hf:
        assert "x" in hf
        assert "y" in hf
        assert prep_utils.TIME_DSET_NAME not in hf
        assert "spatial_ref" in hf

    _check_gdal_metadata(f"NETCDF:{outfile}:velocity", metadata)


def test_3d(metadata, dates, data3d, tmp_path):
    num_dates, rows, cols = data3d.shape
    ds_name_dict = {
        "date": [dates.dtype, (num_dates,), dates],
        "timeseries": [np.float32, (num_dates, rows, cols), data3d],
    }

    # initiate HDF5 file
    metadata["FILE_TYPE"] = "timeseries"
    # meta["REF_DATE"] = ref_date # might not be the first date!
    outfile = tmp_path / "timeseries.h5"
    writefile.layout_hdf5(str(outfile), ds_name_dict, metadata=metadata)

    # Check the metadata was written
    with h5py.File(outfile) as hf:
        assert "x" in hf
        assert "y" in hf
        assert prep_utils.TIME_DSET_NAME in hf
        assert "spatial_ref" in hf

    _check_gdal_metadata(f"NETCDF:{outfile}:timeseries", metadata)
