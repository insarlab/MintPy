import numpy as np
import pytest
from pyproj import Transformer

from mintpy.tropo_gacos import calculate_delay_timeseries, get_delay_geo
from mintpy.utils import readfile, writefile


def _write_ztd(ztd_path, data, x_first, y_first, x_step, y_step):
    """Write a synthetic GACOS .ztd binary file and its .ztd.rsc metadata."""
    data = np.asarray(data, dtype=np.float32)
    length, width = data.shape
    data.tofile(ztd_path)
    rsc_content = (
        f'WIDTH         {width}\n'
        f'FILE_LENGTH   {length}\n'
        f'X_FIRST       {x_first}\n'
        f'Y_FIRST       {y_first}\n'
        f'X_STEP        {x_step}\n'
        f'Y_STEP        {y_step}\n'
        f'X_UNIT        degrees\n'
        f'Y_UNIT        degrees\n'
        f'Z_OFFSET      0\n'
        f'Z_SCALE       1\n'
        f'PROJECTION    LATLON\n'
        f'DATA_TYPE     float32\n'
    )
    with open(f'{ztd_path}.rsc', 'w') as f:
        f.write(rsc_content)


def _make_linear_ztd_grid(
    ztd_path,
    lat_north=40.0,
    lon_west=118.0,
    dlat=-0.01,
    dlon=0.01,
    nrows=200,
    ncols=250,
    offset=0.5,
    lat_coeff=0.02,
    lon_coeff=-0.015,
):
    """Create a GACOS .ztd grid with exact linear function of lat/lon."""
    lats = lat_north + dlat * (np.arange(nrows) + 0.5)
    lons = lon_west + dlon * (np.arange(ncols) + 0.5)
    lon2d, lat2d = np.meshgrid(lons, lats)
    ztd = offset + lat_coeff * (lat2d - 39.0) + lon_coeff * (lon2d - 119.0)
    _write_ztd(ztd_path, ztd, lon_west, lat_north, dlon, dlat)
    return offset, lat_coeff, lon_coeff


def test_get_delay_geo_latlon(tmp_path):
    """Regression test for standard WGS84 lat/lon geocoded products."""
    ztd_file = tmp_path / '20220101.ztd'
    offset, c_lat, c_lon = _make_linear_ztd_grid(ztd_file)

    length, width = 20, 25
    atr = {
        'LENGTH': length,
        'WIDTH': width,
        'X_FIRST': 118.5,
        'Y_FIRST': 39.5,
        'X_STEP': 0.005,
        'Y_STEP': -0.005,
        'X_UNIT': 'degrees',
        'Y_UNIT': 'degrees',
    }
    cos_inc = np.full((length, width), np.cos(np.deg2rad(35.0)), dtype=np.float32)

    delay = get_delay_geo(str(ztd_file), atr, cos_inc)

    lats = atr['Y_FIRST'] + atr['Y_STEP'] * (np.arange(length) + 0.5)
    lons = atr['X_FIRST'] + atr['X_STEP'] * (np.arange(width) + 0.5)
    lon2d, lat2d = np.meshgrid(lons, lats)
    expected = -(offset + c_lat * (lat2d - 39.0) + c_lon * (lon2d - 119.0)) / cos_inc
    assert delay.shape == (length, width)
    np.testing.assert_allclose(delay[1:-1, 1:-1], expected[1:-1, 1:-1], atol=1e-4)


def test_get_delay_geo_utm_epsg(tmp_path):
    """HyP3-style UTM product with EPSG and UTM_ZONE (near zone edge rotation)."""
    ztd_file = tmp_path / '20220101.ztd'
    offset, c_lat, c_lon = _make_linear_ztd_grid(ztd_file)

    # UTM Zone 50N near eastern edge (~119.4E, 39.1N) where grid convergence > 1.5 deg
    length, width = 30, 40
    atr = {
        'LENGTH': length,
        'WIDTH': width,
        'X_FIRST': 680000.0,
        'Y_FIRST': 4350000.0,
        'X_STEP': 1000.0,
        'Y_STEP': -1000.0,
        'X_UNIT': 'meters',
        'Y_UNIT': 'meters',
        'EPSG': '32650',
        'UTM_ZONE': '50N',
    }
    cos_inc = np.full((length, width), np.cos(np.deg2rad(38.0)), dtype=np.float32)

    delay = get_delay_geo(str(ztd_file), atr, cos_inc)

    xs = atr['X_FIRST'] + atr['X_STEP'] * (np.arange(width) + 0.5)
    ys = atr['Y_FIRST'] + atr['Y_STEP'] * (np.arange(length) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    lons, lats = Transformer.from_crs('EPSG:32650', 'EPSG:4326', always_xy=True).transform(xx, yy)
    expected = -(offset + c_lat * (lats - 39.0) + c_lon * (lons - 119.0)) / cos_inc

    assert delay.shape == (length, width)
    np.testing.assert_allclose(delay, expected, atol=1e-4)


def test_get_delay_geo_utm_zone_without_epsg(tmp_path):
    """UTM product having UTM_ZONE without EPSG (e.g. GAMMA or stripped metadata)."""
    ztd_file = tmp_path / '20220101.ztd'
    offset, c_lat, c_lon = _make_linear_ztd_grid(
        ztd_file,
        lat_north=-20.0,
        lon_west=32.0,
        dlat=-0.01,
        dlon=0.01,
        nrows=150,
        ncols=150,
    )

    # Southern hemisphere UTM Zone 36S
    length, width = 20, 25
    atr = {
        'LENGTH': length,
        'WIDTH': width,
        'X_FIRST': 450000.0,
        'Y_FIRST': 7750000.0,
        'X_STEP': 1000.0,
        'Y_STEP': -1000.0,
        'X_UNIT': 'meter',
        'Y_UNIT': 'meter',
        'UTM_ZONE': '36S',
    }
    cos_inc = np.ones((length, width), dtype=np.float32)

    delay = get_delay_geo(str(ztd_file), atr, cos_inc)

    xs = atr['X_FIRST'] + atr['X_STEP'] * (np.arange(width) + 0.5)
    ys = atr['Y_FIRST'] + atr['Y_STEP'] * (np.arange(length) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    lons, lats = Transformer.from_crs('EPSG:32736', 'EPSG:4326', always_xy=True).transform(xx, yy)
    expected = -(offset + c_lat * (lats - 39.0) + c_lon * (lons - 119.0))

    assert delay.shape == (length, width)
    np.testing.assert_allclose(delay, expected, atol=1e-4)


def test_get_delay_geo_projected_non_utm_epsg(tmp_path):
    """Projected non-UTM EPSG (Polar Stereographic EPSG:3413) without UTM_ZONE."""
    ztd_file = tmp_path / '20220101.ztd'
    offset, c_lat, c_lon = _make_linear_ztd_grid(
        ztd_file,
        lat_north=72.0,
        lon_west=-52.0,
        dlat=-0.02,
        dlon=0.04,
        nrows=150,
        ncols=200,
    )

    length, width = 16, 20
    atr = {
        'LENGTH': length,
        'WIDTH': width,
        'X_FIRST': -200000.0,
        'Y_FIRST': -2100000.0,
        'X_STEP': 1000.0,
        'Y_STEP': -1000.0,
        'X_UNIT': 'meters',
        'Y_UNIT': 'meters',
        'EPSG': '3413',
    }
    cos_inc = np.full((length, width), 0.8, dtype=np.float32)

    delay = get_delay_geo(str(ztd_file), atr, cos_inc)

    xs = atr['X_FIRST'] + atr['X_STEP'] * (np.arange(width) + 0.5)
    ys = atr['Y_FIRST'] + atr['Y_STEP'] * (np.arange(length) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    lons, lats = Transformer.from_crs('EPSG:3413', 'EPSG:4326', always_xy=True).transform(xx, yy)
    expected = -(offset + c_lat * (lats - 39.0) + c_lon * (lons - 119.0)) / cos_inc

    assert delay.shape == (length, width)
    np.testing.assert_allclose(delay, expected, atol=1e-4)


def test_get_delay_geo_utm_non_overlapping_raises(tmp_path):
    """Non-overlapping ZTD coverage should raise a descriptive ValueError."""
    ztd_file = tmp_path / '20220101.ztd'
    _make_linear_ztd_grid(ztd_file, lat_north=10.0, lon_west=10.0, nrows=50, ncols=50)

    atr = {
        'LENGTH': 10,
        'WIDTH': 10,
        'X_FIRST': 680000.0,
        'Y_FIRST': 4350000.0,
        'X_STEP': 1000.0,
        'Y_STEP': -1000.0,
        'X_UNIT': 'meters',
        'Y_UNIT': 'meters',
        'EPSG': '32650',
        'UTM_ZONE': '50N',
    }
    cos_inc = np.ones((10, 10), dtype=np.float32)

    with pytest.raises(ValueError, match='does NOT overlap'):
        get_delay_geo(str(ztd_file), atr, cos_inc)


def test_calculate_delay_timeseries_utm(tmp_path):
    """End-to-end calculate_delay_timeseries on HyP3 UTM timeseries & geometry."""
    gacos_dir = tmp_path / 'GACOS'
    gacos_dir.mkdir()
    dates = ['20220101', '20220113']
    coeffs = []
    for idx, date_str in enumerate(dates):
        ztd_file = gacos_dir / f'{date_str}.ztd'
        coeffs.append(_make_linear_ztd_grid(ztd_file, offset=0.4 + 0.2 * idx))

    length, width = 16, 20
    base_atr = {
        'LENGTH': str(length),
        'WIDTH': str(width),
        'X_FIRST': '680000.0',
        'Y_FIRST': '4350000.0',
        'X_STEP': '1000.0',
        'Y_STEP': '-1000.0',
        'X_UNIT': 'meters',
        'Y_UNIT': 'meters',
        'EPSG': '32650',
        'UTM_ZONE': '50N',
        'PROCESSOR': 'hyp3',
        'UNIT': 'm',
    }

    dis_file = tmp_path / 'timeseries.h5'
    ts_atr = {**base_atr, 'FILE_TYPE': 'timeseries'}
    dates_arr = np.array(dates, dtype=np.bytes_)
    writefile.layout_hdf5(
        str(dis_file),
        {
            'date': [dates_arr.dtype, (len(dates),), dates_arr],
            'timeseries': [np.float32, (len(dates), length, width), np.zeros((len(dates), length, width), np.float32)],
        },
        metadata=ts_atr,
    )

    geom_file = tmp_path / 'geometryGeo.h5'
    geom_atr = {**base_atr, 'FILE_TYPE': 'geometry'}
    inc_angle = np.full((length, width), 36.0, dtype=np.float32)
    writefile.layout_hdf5(
        str(geom_file),
        {'incidenceAngle': [np.float32, (length, width), inc_angle]},
        metadata=geom_atr,
    )

    tropo_file = tmp_path / 'GACOS.h5'
    out = calculate_delay_timeseries(str(tropo_file), str(dis_file), str(geom_file), str(gacos_dir))
    assert out == str(tropo_file)

    data, out_atr = readfile.read(str(tropo_file))
    assert data.shape == (len(dates), length, width)
    assert out_atr['UTM_ZONE'] == '50N'

    xs = float(base_atr['X_FIRST']) + float(base_atr['X_STEP']) * (np.arange(width) + 0.5)
    ys = float(base_atr['Y_FIRST']) + float(base_atr['Y_STEP']) * (np.arange(length) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    lons, lats = Transformer.from_crs('EPSG:32650', 'EPSG:4326', always_xy=True).transform(xx, yy)
    cos_inc = np.cos(np.deg2rad(inc_angle))
    for idx, (offset, c_lat, c_lon) in enumerate(coeffs):
        expected = -(offset + c_lat * (lats - 39.0) + c_lon * (lons - 119.0)) / cos_inc
        np.testing.assert_allclose(data[idx], expected, atol=1e-4)
