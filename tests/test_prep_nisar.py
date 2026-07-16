"""Unit coverage for NISAR GUNW metadata parsing in ``mintpy.prep_nisar``.

The tests build minimal synthetic GUNW HDF5 files with ``h5py`` (a few KB,
no GDAL warp and no network) and exercise the pure-Python / h5py parsing
layer of ``prep_nisar``: frequency resolution, required-path discovery and
metadata extraction.

They lock in regression coverage for the metadata reader, most importantly
issue #1485: the reader must read ``.../radarGrid/referenceSlantRange`` and
not a bare ``slantRange`` dataset (which does not exist in real GUNW files
and used to raise ``KeyError``).

The GDAL-backed data path (``read_subset`` -> ``_warp_to_grid_mem`` and the
radarGrid interpolation) is intentionally out of scope here; it is validated
against a real sample GUNW product separately.
"""

import datetime

import h5py
import numpy as np
import pytest

from mintpy import prep_nisar
from mintpy.constants import EARTH_RADIUS, SPEED_OF_LIGHT

GUNW_ROOT = "/science/LSAR/GUNW"
IDENT = "/science/LSAR/identification"
RADARGRID = f"{GUNW_ROOT}/metadata/radarGrid"


def _unw_root(freq):
    return f"{GUNW_ROOT}/grids/{freq}/unwrappedInterferogram"


def _params_root(freq):
    return (
        f"{GUNW_ROOT}/metadata/processingInformation/parameters/"
        f"unwrappedInterferogram/{freq}"
    )


def _build_gunw(
    path,
    *,
    frequencies=("frequencyA",),
    polarization="HH",
    epsg=32612,
    x_start=400000.0,
    y_start=4300000.0,
    x_step=80.0,
    y_step=-80.0,
    nx=5,
    ny=4,
    center_frequency=1.2575e9,
    slant_range_name="referenceSlantRange",
    orbit_direction="ASCENDING",
    mission_id="NISAR",
    start_time="2008-10-12T06:09:11.000000",
    end_time="2008-10-12T06:09:25.000000",
    range_looks=20,
    azimuth_looks=4,
    include_ifgram=True,
    include_ion=True,
    include_tropo=True,
    include_set=True,
):
    """Write a minimal synthetic NISAR GUNW HDF5 product to ``path``.

    Datasets mirror the paths that the GDAL-free readers touch. String fields
    are stored as byte strings to exercise the ``.decode('utf-8')`` path, and
    spacing / looks / EPSG are stored as scalar datasets on purpose.
    """
    xcoord = x_start + np.arange(nx) * x_step
    ycoord = y_start + np.arange(ny) * y_step
    shape = (ny, nx)

    with h5py.File(path, "w") as f:
        # identification group -- byte strings exercise .decode('utf-8')
        f[f"{IDENT}/orbitPassDirection"] = np.bytes_(orbit_direction)
        f[f"{IDENT}/missionId"] = np.bytes_(mission_id)
        f[f"{IDENT}/referenceZeroDopplerStartTime"] = np.bytes_(start_time)
        f[f"{IDENT}/referenceZeroDopplerEndTime"] = np.bytes_(end_time)

        # radarGrid cube: coords + reference slant range (issue #1485 dataset)
        f[f"{RADARGRID}/xCoordinates"] = xcoord
        f[f"{RADARGRID}/yCoordinates"] = ycoord
        slant_range = np.full(shape, 8.0e5, dtype=np.float64) + np.arange(nx)
        f[f"{RADARGRID}/{slant_range_name}"] = slant_range
        if include_tropo:
            f[f"{RADARGRID}/wetTroposphericPhaseScreen"] = np.zeros(shape, np.float32)
            f[f"{RADARGRID}/hydrostaticTroposphericPhaseScreen"] = np.zeros(shape, np.float32)
        if include_set:
            f[f"{RADARGRID}/slantRangeSolidEarthTidesPhase"] = np.zeros(shape, np.float32)

        for freq in frequencies:
            grids = f"{GUNW_ROOT}/grids/{freq}"
            # scalar center frequency
            f[f"{grids}/centerFrequency"] = np.float64(center_frequency)

            unw = f"{_unw_root(freq)}/{polarization}"
            f[f"{unw}/xCoordinates"] = xcoord
            f[f"{unw}/yCoordinates"] = ycoord
            # scalar spacing datasets
            f[f"{unw}/xCoordinateSpacing"] = np.float64(x_step)
            f[f"{unw}/yCoordinateSpacing"] = np.float64(y_step)
            # scalar integer EPSG code
            f[f"{unw}/projection"] = np.int32(epsg)
            if include_ifgram:
                f[f"{unw}/unwrappedPhase"] = np.zeros(shape, np.float32)
                f[f"{unw}/coherenceMagnitude"] = np.ones(shape, np.float32)
                f[f"{unw}/connectedComponents"] = np.zeros(shape, np.int32)
            if include_ion:
                f[f"{unw}/ionospherePhaseScreen"] = np.zeros(shape, np.float32)

            params = _params_root(freq)
            # scalar look counts
            f[f"{params}/numberOfRangeLooks"] = np.int32(range_looks)
            f[f"{params}/numberOfAzimuthLooks"] = np.int32(azimuth_looks)

    return str(path)


# ---------------------------------------------------------------------------
# Tier A -- pure helpers (no HDF5)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "value, expected",
    [
        (None, "frequencyA"),
        ("auto", "frequencyA"),
        ("A", "frequencyA"),
        ("B", "frequencyB"),
        ("frequencyA", "frequencyA"),
        ("frequencyB", "frequencyB"),
    ],
)
def test_normalize_frequency_valid(value, expected):
    assert prep_nisar._normalize_frequency(value) == expected


def test_normalize_frequency_invalid():
    with pytest.raises(ValueError, match="auto, A, B"):
        prep_nisar._normalize_frequency("C")


def test_dataset_path_builders():
    assert prep_nisar._dataset_root_unw("frequencyA").endswith(
        "grids/frequencyA/unwrappedInterferogram"
    )
    assert prep_nisar._center_frequency_path("frequencyB").endswith(
        "grids/frequencyB/centerFrequency"
    )
    datasets = prep_nisar._datasets_for_pol("HH", "frequencyA")
    assert datasets["unw"].endswith("frequencyA/unwrappedInterferogram/HH/unwrappedPhase")
    assert datasets["cor"].endswith("/HH/coherenceMagnitude")
    assert datasets["connComp"].endswith("/HH/connectedComponents")
    assert datasets["epsg"].endswith("/HH/projection")


@pytest.mark.parametrize(
    "stack_type, tail",
    [
        ("ifgram", "unwrappedPhase"),
        ("ion", "ionospherePhaseScreen"),
        ("tropo", "wetTroposphericPhaseScreen"),
        ("set", "slantRangeSolidEarthTidesPhase"),
    ],
)
def test_required_paths_for_stack_type(stack_type, tail):
    paths = prep_nisar._required_paths_for_stack_type(stack_type, "HH", "frequencyA")
    assert any(p.endswith(tail) for p in paths)


def test_required_paths_for_stack_type_unknown():
    with pytest.raises(ValueError, match="Unsupported stack_type"):
        prep_nisar._required_paths_for_stack_type("bogus", "HH", "frequencyA")


def test_processinfo_reads_reference_slant_range():
    """Regression sentinel for #1485: the reader targets referenceSlantRange."""
    slant_path = prep_nisar.PROCESSINFO["rdr_slant_range"]
    assert slant_path.endswith("/referenceSlantRange")
    # guard against a regression back to the bare (non-existent) dataset name
    assert not slant_path.endswith("/slantRange")


# ---------------------------------------------------------------------------
# Tier B -- HDF5 traversal: frequency resolution and required-path discovery
# ---------------------------------------------------------------------------
def test_resolve_frequency_auto_and_explicit(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_A.h5", frequencies=("frequencyA",))
    assert prep_nisar._resolve_frequency(gunw, None, "HH") == "frequencyA"
    assert prep_nisar._resolve_frequency(gunw, "auto", "HH") == "frequencyA"
    assert prep_nisar._resolve_frequency(gunw, "A", "HH") == "frequencyA"


def test_resolve_frequency_missing_b_raises(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_A.h5", frequencies=("frequencyA",))
    with pytest.raises(ValueError, match="frequencyB"):
        prep_nisar._resolve_frequency(gunw, "B", "HH")


def test_resolve_frequency_missing_polarization_raises(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_hh.h5", polarization="HH")
    with pytest.raises(ValueError, match="VV"):
        prep_nisar._resolve_frequency(gunw, "A", "VV")


def test_missing_required_paths_complete(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_full.h5")
    for stack_type in ["ifgram", "ion", "tropo", "set"]:
        missing = prep_nisar._missing_required_paths(
            [gunw], stack_type, "HH", "frequencyA"
        )
        assert missing == [], f"{stack_type} should have no missing paths"


def test_missing_required_paths_detects_missing(tmp_path):
    # build a file whose ifgram layer lacks connectedComponents
    gunw = tmp_path / "gunw_partial.h5"
    _build_gunw(gunw, include_ifgram=True)
    with h5py.File(gunw, "a") as f:
        del f[f"{_unw_root('frequencyA')}/HH/connectedComponents"]

    missing = prep_nisar._missing_required_paths(
        [str(gunw)], "ifgram", "HH", "frequencyA"
    )
    assert len(missing) == 1
    assert missing[0][1].endswith("/HH/connectedComponents")


# ---------------------------------------------------------------------------
# Tier C -- extract_metadata (the #1485 regression target)
# ---------------------------------------------------------------------------
def test_extract_metadata_utm(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_utm.h5", epsg=32612, nx=5, ny=4)
    meta, bounds = prep_nisar.extract_metadata([gunw])

    assert meta["EPSG"] == 32612
    assert meta["X_UNIT"] == "meters"
    assert meta["Y_UNIT"] == "meters"
    assert meta["UTM_ZONE"] == "12N"
    assert meta["ORBIT_DIRECTION"] == "ASCENDING"
    assert meta["POLARIZATION"] == "HH"
    assert meta["PLATFORM"] == "NISAR"
    assert int(meta["RLOOKS"]) == 20
    assert int(meta["ALOOKS"]) == 4
    assert int(meta["LENGTH"]) == 4
    assert int(meta["WIDTH"]) == 5
    assert meta["X_STEP"] == pytest.approx(80.0)
    assert meta["Y_STEP"] == pytest.approx(-80.0)
    assert meta["HEIGHT"] == 747000
    assert float(meta["EARTH_RADIUS"]) == pytest.approx(EARTH_RADIUS)
    assert float(meta["WAVELENGTH"]) == pytest.approx(SPEED_OF_LIGHT / 1.2575e9)
    # STARTING_RANGE comes from radarGrid/referenceSlantRange (min = 8.0e5)
    assert float(meta["STARTING_RANGE"]) == pytest.approx(8.0e5)
    # mid of 06:09:11 -> 06:09:25 is 06:09:18 = 22158 s past midnight
    assert float(meta["CENTER_LINE_UTC"]) == pytest.approx(22158.0)
    assert len(bounds) == 4


def test_extract_metadata_southern_utm_zone(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_utm_s.h5", epsg=32712)
    meta, _ = prep_nisar.extract_metadata([gunw])
    assert meta["UTM_ZONE"] == "12S"


def test_extract_metadata_geographic_4326(tmp_path):
    gunw = _build_gunw(
        tmp_path / "gunw_geo.h5",
        epsg=4326,
        x_start=-78.0,
        y_start=0.98,
        x_step=0.001,
        y_step=-0.001,
    )
    meta, _ = prep_nisar.extract_metadata([gunw])
    assert meta["EPSG"] == 4326
    assert meta["X_UNIT"] == "degree"
    assert meta["Y_UNIT"] == "degree"
    assert "UTM_ZONE" not in meta


def test_extract_metadata_ascending_y(tmp_path):
    # ascending y coordinates (positive y_step) must still parse cleanly
    gunw = _build_gunw(tmp_path / "gunw_asc_y.h5", y_start=4200000.0, y_step=80.0)
    meta, _ = prep_nisar.extract_metadata([gunw])
    assert meta["Y_STEP"] == pytest.approx(80.0)
    assert int(meta["LENGTH"]) == 4
    assert int(meta["WIDTH"]) == 5


def test_extract_metadata_requires_reference_slant_range(tmp_path):
    """Guard #1485: a file with a bare ``slantRange`` (and no
    ``referenceSlantRange``) must fail -- the old buggy behavior."""
    gunw = _build_gunw(
        tmp_path / "gunw_bad_slant.h5", slant_range_name="slantRange"
    )
    # match on the dataset name so an unrelated KeyError cannot pass this guard
    with pytest.raises(KeyError, match="referenceSlantRange"):
        prep_nisar.extract_metadata([gunw])


def test_extract_metadata_center_line_utc_matches_manual(tmp_path):
    start = "2020-01-01T00:00:00.000000"
    end = "2020-01-01T00:01:00.000000"
    gunw = _build_gunw(
        tmp_path / "gunw_utc.h5", start_time=start, end_time=end
    )
    meta, _ = prep_nisar.extract_metadata([gunw])
    t0 = datetime.datetime.fromisoformat(start)
    t1 = datetime.datetime.fromisoformat(end)
    t_mid = t0 + (t1 - t0) / 2.0
    expected = (t_mid - datetime.datetime(2020, 1, 1)).total_seconds()
    assert float(meta["CENTER_LINE_UTC"]) == pytest.approx(expected)


def test_extract_metadata_origin_and_bounds_values(tmp_path):
    """Pin the half-pixel origin correction and common bounds numerically."""
    gunw = _build_gunw(tmp_path / "gunw_origin.h5")  # default 5x4 grid
    meta, bounds = prep_nisar.extract_metadata([gunw])
    # X_FIRST = min(x) - x_step/2 = 400000 - 40; Y_FIRST = max(y) - y_step/2
    assert float(meta["X_FIRST"]) == pytest.approx(399960.0)
    assert float(meta["Y_FIRST"]) == pytest.approx(4300040.0)
    # bounds = (west, south, east, north) at pixel centers
    assert tuple(pytest.approx(v) for v in bounds) == (
        pytest.approx(400000.0),
        pytest.approx(4299760.0),
        pytest.approx(400320.0),
        pytest.approx(4300000.0),
    )


# ---------------------------------------------------------------------------
# frequencyB positive path (not just normalization / missing-B)
# ---------------------------------------------------------------------------
def test_resolve_frequency_b_positive(tmp_path):
    gunw = _build_gunw(
        tmp_path / "gunw_ab.h5", frequencies=("frequencyA", "frequencyB")
    )
    assert prep_nisar._resolve_frequency(gunw, "B", "HH") == "frequencyB"
    # auto still prefers frequencyA even when B exists
    assert prep_nisar._resolve_frequency(gunw, "auto", "HH") == "frequencyA"


def test_extract_metadata_frequency_b(tmp_path):
    gunw = _build_gunw(
        tmp_path / "gunw_ab.h5", frequencies=("frequencyA", "frequencyB")
    )
    meta, _ = prep_nisar.extract_metadata([gunw], frequency="frequencyB")
    assert meta["POLARIZATION"] == "HH"
    assert int(meta["WIDTH"]) == 5
    assert float(meta["STARTING_RANGE"]) == pytest.approx(8.0e5)


# ---------------------------------------------------------------------------
# Multi-file common bounds
# ---------------------------------------------------------------------------
def test_get_raster_corners(tmp_path):
    gunw = _build_gunw(tmp_path / "gunw_corners.h5")
    west, south, east, north = prep_nisar.get_raster_corners(gunw)
    assert (west, south, east, north) == (
        pytest.approx(400000.0),
        pytest.approx(4299760.0),
        pytest.approx(400320.0),
        pytest.approx(4300000.0),
    )


def test_common_raster_bound_intersection(tmp_path):
    # second scene shifted east by two pixels (160 m) -> overlapping range
    f1 = _build_gunw(tmp_path / "gunw_a.h5", x_start=400000.0)
    f2 = _build_gunw(tmp_path / "gunw_b.h5", x_start=400160.0)
    west, south, east, north = prep_nisar.common_raster_bound([f1, f2])
    assert west == pytest.approx(400160.0)  # max of the two wests
    assert east == pytest.approx(400320.0)  # min of the two easts
    assert south == pytest.approx(4299760.0)
    assert north == pytest.approx(4300000.0)


def test_common_raster_bound_no_overlap_raises(tmp_path):
    f1 = _build_gunw(tmp_path / "gunw_a.h5", x_start=400000.0)
    f2 = _build_gunw(tmp_path / "gunw_far.h5", x_start=500000.0)
    with pytest.raises(ValueError, match="No common overlap"):
        prep_nisar.common_raster_bound([f1, f2])
