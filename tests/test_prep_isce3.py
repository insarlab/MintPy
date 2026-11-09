"""Unit tests for the ISCE3/Dolphin interface."""


import numpy as np
from osgeo import gdal, osr

from mintpy.utils import readfile


def _write_synthetic_geotiff(fname, shape=(8, 10), gt=(500000.0, 10.0, 0.0, 4000000.0, 0.0, -10.0),
                             epsg=32605, nodata=0.0):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(str(fname), shape[1], shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    data = np.ones(shape, dtype=np.float32) * 1.5
    data[0, :] = nodata
    band.WriteArray(data)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    ds = None
    return fname


def test_read_attribute_isce3_unw_tif(tmp_path):
    # Dolphin-style name with the date pair in the filename
    fname = str(tmp_path / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['DATE12'] == '220105-220117'
    assert atr['FILE_TYPE'] == '.unw'
    assert atr['LENGTH'] == 8
    assert atr['WIDTH'] == 10
    assert atr['EPSG'] == '32605'
    assert atr['NO_DATA_VALUE'] == '0.0'


def test_read_attribute_isce3_fullres_parent_dir(tmp_path):
    # Dolphin "fullres.unw.tif" layout: date pair only in the parent dir name
    ifg_dir = tmp_path / '20220105_20220117'
    ifg_dir.mkdir()
    fname = str(ifg_dir / 'fullres.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['DATE12'] == '220105-220117'
    assert atr['FILE_TYPE'] == '.unw'


def test_read_attribute_isce3_rsc_priority(tmp_path):
    # a .rsc sidecar written by prep_isce3 must take priority over the
    # attributes derived from the GeoTIFF itself
    fname = str(tmp_path / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    with open(fname + '.rsc', 'w') as f:
        f.write('PROCESSOR isce3\n')
        f.write('DATE12 220105-220117\n')
        f.write('P_BASELINE_TOP_HDR -2.5048807779724456\n')
        f.write('P_BASELINE_BOTTOM_HDR -2.5048807779724456\n')

    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['P_BASELINE_TOP_HDR'] == '-2.5048807779724456'
    assert 'FILE_PATH' in atr


def test_read_isce3_geotiff_no_nodata(tmp_path):
    # a band without nodata must not produce the string 'None'
    fname = str(tmp_path / '20220105_20220117.int.tif')
    _write_synthetic_geotiff(fname, nodata=None)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'isce3'
    assert 'NO_DATA_VALUE' not in atr or atr['NO_DATA_VALUE'] != 'None'


#########################################################################
# water mask auto-align (geometryDict._warp_water_mask)
#########################################################################
