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


#########################################################################
# geometry merging / angle computation (synthetic data, no real data needed)
#########################################################################
def _write_geotiff(fname, data, gt=(500000.0, 10.0, 0.0, 4000000.0, 0.0, -10.0), epsg=32605):
    driver = gdal.GetDriverByName('GTiff')
    data = np.asarray(data, dtype=np.float32)
    ds = driver.Create(str(fname), data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(data)
    ds = None
    return fname


def test_warp_water_mask_box_and_multilook(tmp_path):
    # reference grid 8 x 10 at 10 m; water mask on the same grid: left half water, right half land
    ref_fname = str(tmp_path / 'ref.tif')
    _write_synthetic_geotiff(ref_fname, shape=(8, 10))
    mask = np.zeros((8, 10), dtype=np.float32)
    mask[:, 5:] = 1.0
    mask_fname = str(tmp_path / 'water_mask.tif')
    _write_geotiff(mask_fname, mask)

    geom_obj = geometryDict(
        processor='isce3',
        datasetDict={'height': ref_fname, 'waterMask': mask_fname},
        extraMetadata={'dummy': 'x'},
    )

    # box=(x0, y0, x1, y1) = the lower-right land half, at half resolution
    result = geom_obj._warp_water_mask('waterMask', 4, 5, box=(5, 4, 10, 8), xstep=1, ystep=1)
    assert result.shape == (4, 5)
    assert np.all(result == 1.0)

    # full grid is returned unchanged when no box/step is given
    result_full = geom_obj._warp_water_mask('waterMask', 8, 10)
    assert np.all(result_full[:, :5] == 0.0)
    assert np.all(result_full[:, 5:] == 1.0)


def test_build_vrt_from_h5_requires_coordinate_arrays(tmp_path):
    import h5py

    # missing x/y coordinates
    h5_file = tmp_path / 'static_layers_missing.h5'
    with h5py.File(h5_file, 'w') as f:
        grp = f.create_group('data')
        grp.create_dataset('los_east', data=np.zeros((2, 2), dtype=np.float32))
    with pytest.raises(KeyError):
        isce3_utils.build_vrt_from_h5(h5_file, tmp_path / 'vrt', ['los_east.tif'])

    # degenerate 1-pixel coordinate arrays
    h5_file = tmp_path / 'static_layers_1px.h5'
    with h5py.File(h5_file, 'w') as f:
        grp = f.create_group('data')
        grp.create_dataset('x_coordinates', data=np.array([0.0], dtype=np.float64))
        grp.create_dataset('y_coordinates', data=np.array([0.0], dtype=np.float64))
        grp.create_dataset('los_east', data=np.zeros((1, 1), dtype=np.float32))
    with pytest.raises(ValueError):
        isce3_utils.build_vrt_from_h5(h5_file, tmp_path / 'vrt', ['los_east.tif'])


def test_compute_angles_with_finite_nodata(tmp_path):
    nodata = -9999.0
    east = np.full((4, 5), 0.1, dtype=np.float32)
    north = np.full((4, 5), 0.2, dtype=np.float32)
    east[0, :] = nodata
    north[0, :] = nodata
    east_file = _write_geotiff(tmp_path / 'los_east.tif', east)
    north_file = _write_geotiff(tmp_path / 'los_north.tif', north)

    inc_file = tmp_path / 'incidenceAngle.tif'
    isce3_utils.compute_incidence_angle(east_file, north_file, inc_file, nodata=nodata)
    inc = gdal.Open(str(inc_file)).ReadAsArray()
    assert np.all(inc[0, :] == nodata)
    assert np.all(np.isfinite(inc[1:, :]))

    az_file = tmp_path / 'azimuthAngle.tif'
    isce3_utils.compute_azimuth_angle(east_file, north_file, az_file, nodata=nodata)
    az = gdal.Open(str(az_file)).ReadAsArray()
    assert np.all(az[0, :] == nodata)
    assert np.all((az[1:, :] >= 0) & (az[1:, :] < 360))


#########################################################################
# prep_isce3.py: RSC generation and geometry preparation
#########################################################################
def _write_unw_tif(fname):
    return _write_synthetic_geotiff(fname, shape=(8, 10))


def test_prepare_stack_isce3_writes_rsc(tmp_path):
    from mintpy.prep_isce3 import prepare_stack_isce3

    _write_unw_tif(tmp_path / '20210104_20210110.unw.tif')
    _write_unw_tif(tmp_path / '20210110_20210116.unw.tif')
    _write_unw_tif(tmp_path / 'no_dates.unw.tif')

    baseline_dict = {'20210104': [0.0, 0.0], '20210110': [-85.65, -85.65]}
    prepare_stack_isce3(
        str(tmp_path / '*.unw.tif'),
        metadata={'PROCESSOR': 'isce3'},
        baseline_dict=baseline_dict,
        update_mode=False,
    )

    rsc = readfile.read_roipac_rsc(str(tmp_path / '20210104_20210110.unw.tif.rsc'))
    assert rsc['PROCESSOR'] == 'isce3'
    assert rsc['DATE12'] == '210104-210110'
    assert float(rsc['P_BASELINE_TOP_HDR']) == pytest.approx(-85.65)
    assert not os.path.isfile(str(tmp_path / 'no_dates.unw.tif.rsc'))


def test_prepare_stack_isce3_date_pair_from_parent_dir(tmp_path):
    from mintpy.prep_isce3 import prepare_stack_isce3

    ifg_dir = tmp_path / '20220105_20220117'
    ifg_dir.mkdir()
    _write_unw_tif(ifg_dir / 'fullres.unw.tif')

    prepare_stack_isce3(str(ifg_dir / '*.unw.tif'), metadata={'PROCESSOR': 'isce3'}, update_mode=False)
    rsc = readfile.read_roipac_rsc(str(ifg_dir / 'fullres.unw.tif.rsc'))
    assert rsc['DATE12'] == '220105-220117'


def test_prepare_stack_isce3_empty_pattern():
    from mintpy.prep_isce3 import prepare_stack_isce3

    with pytest.raises(FileNotFoundError):
        prepare_stack_isce3('/no/such/dir/*.unw.tif')


def _write_static_layers_h5(fname, length=8, width=10, dx=10.0):
    import h5py

    with h5py.File(fname, 'w') as f:
        grp = f.create_group('data')
        grp.create_dataset('x_coordinates', data=np.arange(width) * dx)
        grp.create_dataset('y_coordinates', data=np.arange(length) * dx)
        grp.create_dataset('projection', data=np.int32(32605))
        grp.create_dataset('z', data=np.zeros((length, width), dtype=np.float32))
        grp.create_dataset('los_east', data=np.zeros((length, width), dtype=np.float32))
        grp.create_dataset('los_north', data=np.zeros((length, width), dtype=np.float32))
    return fname


def test_prepare_geometry_isce3_writes_rsc(tmp_path, monkeypatch):
    from mintpy.prep_isce3 import prepare_geometry_isce3

    geom_dir = tmp_path / 'geom' / 't124_264305_iw2'
    geom_dir.mkdir(parents=True)
    _write_static_layers_h5(geom_dir / 'static_layers.h5')

    ref_int = _write_unw_tif(tmp_path / 'ref.unw.tif')
    (tmp_path / 'merged').mkdir()
    height = _write_synthetic_geotiff(tmp_path / 'merged' / 'height.tif', shape=(8, 10))

    def fake_extract_merge_geometry(**kwargs):
        return {'height.tif': height}

    monkeypatch.setattr(isce3_utils, 'extract_merge_geometry', fake_extract_merge_geometry)

    metadata = {'PROCESSOR': 'isce3'}
    metadata = prepare_geometry_isce3(
        geom_dir=str(tmp_path / 'geom'),
        out_dir=str(tmp_path / 'merged'),
        metadata=metadata,
        update_mode=False,
        ref_int_file=ref_int,
    )

    assert os.path.isfile(str(height) + '.rsc')
    assert metadata['LENGTH'] == '8'
    assert metadata['WIDTH'] == '10'
    assert metadata['ALOOKS'] == '1'
    assert metadata['RLOOKS'] == '1'


def test_add_ifgram_metadata_invalid_dates():
    from mintpy.prep_isce3 import add_ifgram_metadata

    with pytest.raises(ValueError):
        add_ifgram_metadata({'PROCESSOR': 'isce3'}, dates=['20210104'])


#########################################################################
# cli/prep_isce3.py
#########################################################################
def test_cli_prep_isce3_defaults(tmp_path):
    from mintpy.cli.prep_isce3 import cmd_line_parse

    geom_dir = tmp_path / 'geom'
    geom_dir.mkdir()
    inps = cmd_line_parse(['-g', str(geom_dir)])

    assert inps.out_dir == str(tmp_path / 'merged_geom')
    assert inps.geom_files == isce3_utils.GEOMETRY_FILENAMES
    assert inps.update_mode is True
    assert inps.obs_files is None
    assert inps.geom_dirs == [str(geom_dir)]

    inps = cmd_line_parse(['-g', str(geom_dir), '--force', '--out-dir', str(tmp_path / 'out')])
    assert inps.update_mode is False
    assert inps.out_dir == str(tmp_path / 'out')


def test_cli_prep_isce3_geom_glob(tmp_path):
    from mintpy.cli.prep_isce3 import cmd_line_parse

    for name in ['t124_a', 't124_b']:
        burst_dir = tmp_path / name / '20210104'
        burst_dir.mkdir(parents=True)
        _write_static_layers_h5(burst_dir / 'static_layers.h5')

    inps = cmd_line_parse(['-g', str(tmp_path / 't124*')])
    assert inps.geom_dir == str(tmp_path / 't124_a')
    assert len(inps.geom_dirs) == 2

    with pytest.raises(FileNotFoundError):
        cmd_line_parse(['-g', str(tmp_path / 'nope*')])
