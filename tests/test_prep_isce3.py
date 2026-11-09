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
    atr = readfile.read_attribute(fname, processor='isce3')
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['DATE12'] == '220105-220117'
    assert atr['FILE_TYPE'] == '.unw'
    assert atr['LENGTH'] == 8
    assert atr['WIDTH'] == 10
    assert atr['EPSG'] == '32605'
    assert atr['NO_DATA_VALUE'] == '0.0'


def test_read_attribute_isce3_per_pair_dir(tmp_path):
    # per-pair-directory layout (the date pair is in the directory name):
    # ifgrams/YYYYMMDD_YYYYMMDD/fullres.unw.tif, date pair in the dir name.
    # (Dolphin itself writes flat YYYYMMDD_YYYYMMDD.int.tif, covered above.)
    ifg_dir = tmp_path / '20220105_20220117'
    ifg_dir.mkdir()
    fname = str(ifg_dir / 'fullres.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname, processor='isce3')
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['DATE12'] == '220105-220117'
    assert atr['FILE_TYPE'] == '.unw'


def test_read_attribute_isce3_rsc_priority(tmp_path):
    # a .rsc sidecar written by prep_isce3 must take priority over the
    # attributes derived from the GeoTIFF itself: write sidecar values that
    # CONFLICT with what the GeoTIFF/filename would produce, so that a
    # reversed precedence order fails this test
    fname = str(tmp_path / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    with open(fname + '.rsc', 'w') as f:
        f.write('PROCESSOR isce3\n')
        f.write('DATE12 990101-990102\n')   # the filename says 220105-220117
        f.write('LENGTH 8\n')
        f.write('WIDTH 10\n')
        f.write('P_BASELINE_TOP_HDR -2.5048807779724456\n')
        f.write('P_BASELINE_BOTTOM_HDR -2.5048807779724456\n')

    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['DATE12'] == '990101-990102'   # sidecar wins over the filename
    assert atr['P_BASELINE_TOP_HDR'] == '-2.5048807779724456'
    assert 'FILE_PATH' in atr


def test_read_isce3_geotiff_no_nodata(tmp_path):
    # a band without nodata must not invent a NO_DATA_VALUE: check the
    # GeoTIFF reader itself (so that an invented '0.0'/'nan'/'none' fails),
    # plus the public no-data-value lookup used by the readers.
    fname = str(tmp_path / '20220105_20220117.int.tif')
    _write_synthetic_geotiff(fname, nodata=None)
    assert 'NO_DATA_VALUE' not in readfile.read_isce3_geotiff(fname)
    assert readfile.read_attribute(fname, processor='isce3')['PROCESSOR'] == 'isce3'
    assert readfile.get_no_data_value(fname) is None


def test_read_attribute_isce3_date12_priority(tmp_path):
    # the interferogram date pair must win over a date pair in an ancestor dir
    root = tmp_path / '20190101_20191231_project'
    ifg_dir = root / 'ifgrams'
    ifg_dir.mkdir(parents=True)

    fname = str(ifg_dir / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    assert readfile.read_attribute(fname, processor='isce3')['DATE12'] == '220105-220117'

    # per-pair-directory layout: date pair in the immediate parent dir
    pair_dir = ifg_dir / '20220205_20220217'
    pair_dir.mkdir()
    fname = str(pair_dir / 'fullres.unw.tif')
    _write_synthetic_geotiff(fname)
    assert readfile.read_attribute(fname, processor='isce3')['DATE12'] == '220205-220217'

    # Dolphin flat .int.tif layout, again under an ancestor with a different
    # date pair
    int_dir = tmp_path / 'S1_20190101_20191231' / 'dolphin' / 'interferograms'
    int_dir.mkdir(parents=True)
    fname = str(int_dir / '20220305_20220317.int.tif')
    _write_synthetic_geotiff(fname)
    assert readfile.read_attribute(fname, processor='isce3')['DATE12'] == '220305-220317'

    # without the hint the file keeps the generic GDAL attributes; with the
    # hint it is read as isce3, but no date pair is invented for it
    fname = str(ifg_dir / 'fullres.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'gdal'
    assert 'DATE12' not in atr

    atr = readfile.read_attribute(fname, processor='isce3')
    assert atr['PROCESSOR'] == 'isce3'
    assert 'DATE12' not in atr


def test_read_attribute_tif_default_is_gdal(tmp_path):
    # GeoTIFFs are not promoted to isce3 based on their name: without an
    # explicit processor hint they keep the generic GDAL attributes
    fname = str(tmp_path / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'gdal'
    assert atr['FILE_TYPE'] != '.unw'
    assert 'DATE12' not in atr

    # the same for a LiCSAR-style name
    fname = str(tmp_path / '20220105_20220117.geo.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'gdal'
    assert atr['FILE_TYPE'] != '.unw'

    # an explicitly configured non-isce3 processor also disables isce3 handling
    fname2 = str(tmp_path / '20220205_20220217.unw.tif')
    _write_synthetic_geotiff(fname2)
    atr = readfile.read_attribute(fname2, processor='roipac')
    assert atr['PROCESSOR'] == 'gdal'
    assert 'DATE12' not in atr


def test_read_attribute_isce3_processor_hint(tmp_path):
    # an explicit processor='isce3' from the configuration is authoritative
    fname = str(tmp_path / '20220105_20220117.geo.unw.tif')
    _write_synthetic_geotiff(fname)
    atr = readfile.read_attribute(fname, processor='isce3')
    assert atr['PROCESSOR'] == 'isce3'
    assert atr['FILE_TYPE'] == '.unw'
    assert atr['DATE12'] == '220105-220117'

    # read() forwards the hint to the attribute reader used for the data
    _, atr = readfile.read(fname, processor='isce3')
    assert atr['PROCESSOR'] == 'isce3'


def test_read_attribute_rsc_sidecar_short_circuit(tmp_path):
    # with a .rsc sidecar, read_attribute() takes the .rsc branch and does not
    # call read_isce3_geotiff(); this is the pre-existing upstream
    # ".rsc > file itself" priority, so a partial sidecar is NOT completed
    # from the GeoTIFF (same behavior as main for any .tif + .rsc)
    fname = str(tmp_path / '20220105_20220117.unw.tif')
    _write_synthetic_geotiff(fname)
    with open(fname + '.rsc', 'w') as f:
        f.write('LENGTH 8\n')
        f.write('WIDTH 10\n')
    atr = readfile.read_attribute(fname)
    assert atr['PROCESSOR'] == 'roipac'
    assert 'DATE12' not in atr
