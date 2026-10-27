"""Unit tests for the ISCE3/Dolphin interface (prep_isce3, isce3_utils, readfile)."""

import os

import numpy as np
import pytest
from osgeo import gdal, osr

from mintpy.objects.stackDict import geometryDict
from mintpy.prep_isce3 import add_ifgram_metadata
from mintpy.utils import isce3_utils, readfile


#########################################################################
# baseline time series
#########################################################################
def _write_baseline_file(dir_path, date_pair, bperp, fmt='Bperp (m)'):
    fname = os.path.join(str(dir_path), f'{date_pair}.txt')
    with open(fname, 'w') as f:
        f.write(f'{fmt}: {bperp}\n')
        f.write('Bpar average (m): 12.345\n')
    return fname


def test_read_baseline_timeseries_isce3(tmp_path):
    # plain directory layout
    _write_baseline_file(tmp_path, '20210104_20210110', -85.65046462265644)
    _write_baseline_file(tmp_path, '20210104_20210116', -19.081015753493432,
                         fmt='Bperp average (m)')
    # a pair NOT starting with the reference date must be ignored
    _write_baseline_file(tmp_path, '20210110_20210116', 99.9)

    bdict = isce3_utils.read_baseline_timeseries_isce3(str(tmp_path))
    assert bdict['20210104'] == [0.0, 0.0]
    assert bdict['20210110'] == pytest.approx([-85.65046462265644] * 2)
    assert bdict['20210116'] == pytest.approx([-19.081015753493432] * 2)
    assert len(bdict) == 3


def test_read_baseline_timeseries_isce3_nested_glob(tmp_path):
    # Dolphin layout: baselines/t124_xxx/YYYYMMDD_YYYYMMDD.txt
    burst_dir = tmp_path / 'baselines' / 't124_264305_iw2'
    burst_dir.mkdir(parents=True)
    _write_baseline_file(burst_dir, '20210104_20210110', -85.65)
    _write_baseline_file(burst_dir, '20210104_20210116', -19.08)

    # glob in the baseline dir itself must be expanded
    bdict = isce3_utils.read_baseline_timeseries_isce3(
        str(tmp_path / 'baselines' / 't124_*/'))
    assert bdict['20210104'] == [0.0, 0.0]
    assert bdict['20210110'] == pytest.approx([-85.65] * 2)

    # plain dir without glob must fall back to a recursive search
    bdict2 = isce3_utils.read_baseline_timeseries_isce3(str(tmp_path / 'baselines'))
    assert bdict2['20210110'] == pytest.approx([-85.65] * 2)

    # a non-existent dir must NOT crash, just return an empty dict
    assert isce3_utils.read_baseline_timeseries_isce3(
        str(tmp_path / 'nope')) == {}


def test_add_ifgram_metadata():
    meta = add_ifgram_metadata(
        {'PROCESSOR': 'isce3'},
        dates=['20210104', '20210110'],
        baseline_dict={'20210104': [0.0, 0.0], '20210110': [-85.65, -85.65]},
    )
    assert meta['DATE12'] == '210104-210110'
    assert meta['P_BASELINE_TOP_HDR'] == '-85.65'
    assert meta['P_BASELINE_BOTTOM_HDR'] == '-85.65'


#########################################################################
# burst XML metadata extraction (CENTER_INCIDENCE_ANGLE regression)
#########################################################################
def _write_burst_xml(path, attrs):
    isce3_utils.save_burst_attributes_to_xml(attrs, str(path))
    return path


def test_extract_isce3_metadata_center_incidence_angle(tmp_path):
    # realistic S1 IW2 values (Hawaii T124, from an actual dolphin run)
    attrs = {
        'prf': 1451.627112193990,
        'burstStartUTC': '2021-01-04 04:30:18.346686',
        'burstStopUTC': '2021-01-04 04:30:21.430020',
        'radarWavelength': 0.05546576,
        'startingRange': 845516.028030508198,
        'passDirection': 'Ascending',
        'polarization': 'VV',
        'trackNumber': 124,
        'orbitNumber': 25000,
        'sensingMid': '2021-01-04 04:30:19.888353',
        'azimuthTimeInterval': 0.0020555563,
        'rangePixelSize': 2.329562114715,
        'swathNumber': 2,
        'satelliteSpeed': 7597.883725532828,
        'HEADING': 347.728905505086,
        'earthRadius': 6346500.407216606,
        'altitude': 697692.716961,
    }
    xml_file = _write_burst_xml(tmp_path / 'IW2.burst.xml', attrs)
    meta = isce3_utils.extract_isce3_metadata(str(xml_file))

    # the old formula returned exactly ~90 deg; the fixed one must be the
    # realistic near-range incidence angle (~36 deg for these values)
    inc = float(meta['CENTER_INCIDENCE_ANGLE'])
    assert 30.0 < inc < 45.0
    assert abs(inc - 90.0) > 45.0

    # sanity of a few other keys
    assert meta['PROCESSOR'] == 'isce3'
    assert meta['PLATFORM'] == 'sen'
    assert float(meta['CENTER_LINE_UTC']) == pytest.approx(4 * 3600 + 30 * 60 + 19.888353)
    assert float(meta['HEIGHT']) == pytest.approx(697692.716961)


def test_read_baseline_missing_dir_warns(tmp_path, capsys):
    isce3_utils.read_baseline_timeseries_isce3(str(tmp_path / 'empty'))
    assert 'no baseline text files found' in capsys.readouterr().out


#########################################################################
# read_attribute on ISCE3/Dolphin GeoTIFF products
#########################################################################
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
def test_warp_water_mask(tmp_path):
    # reference grid: 10 x 8 pixels at 10 m
    ref_fname = str(tmp_path / 'ref.tif')
    _write_synthetic_geotiff(ref_fname, shape=(8, 10))

    # water mask on a coarser grid: 5 x 4 pixels at 20 m, same CRS/extent
    src_fname = str(tmp_path / 'water_mask.tif')
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(src_fname, 5, 4, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((500000.0, 20.0, 0.0, 4000000.0, 0.0, -20.0))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32605)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    mask = np.ones((4, 5), dtype=np.float32)   # 1 = land
    mask[0, :] = 0.0                            # 0 = water
    band.WriteArray(mask)
    ds = None

    geom_obj = geometryDict(
        processor='isce3',
        datasetDict={'height': ref_fname, 'waterMask': src_fname},
        extraMetadata={'dummy': 'x'},
    )
    result = geom_obj._warp_water_mask('waterMask', 8, 10)
    assert result.shape == (8, 10)
    assert set(np.unique(result)).issubset({0.0, 1.0})
    # the top row of the source (water) stays water after the warp
    assert np.all(result[0, :] == 0.0)
    assert np.all(result[-1, :] == 1.0)
