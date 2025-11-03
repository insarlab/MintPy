import pytest

from mintpy.prep_hyp3 import add_hyp3_metadata, _get_product_name_and_type


def test_get_product_name_and_type():
    assert _get_product_name_and_type(
        'S1_136231_IW2_20200604_20200616_VV_INT80_10C1_foo.tif'
    ) == ('S1_136231_IW2_20200604_20200616_VV_INT80_10C1', 'INSAR_ISCE_BURST')

    assert _get_product_name_and_type(
        'S1_064_000000s1n00-136231s2n02-000000s3n00_IW_20200604_20200616_VV_INT80_77F1_foo.tif'
    ) == ('S1_064_000000s1n00-136231s2n02-000000s3n00_IW_20200604_20200616_VV_INT80_77F1', 'INSAR_ISCE_MULTI_BURST')

    assert _get_product_name_and_type(
        'S1AA_20150504T120217_20150621T120220_VVP048_INT80_G_ueF_5CED_foo.tif'
    ) == ('S1AA_20150504T120217_20150621T120220_VVP048_INT80_G_ueF_5CED', 'INSAR_GAMMA')

    # Old INSAR_ISCE_MULTI_BURST naming convention
    with pytest.raises(
        ValueError,
        match=r'^Failed to parse product name from filename: '
              r'S1A_064_E053_1_N27_3_E054_1_N27_8_20200604_20200616_VV_INT80_3FBF_foo\.tif$',
    ):
        _get_product_name_and_type(
            'S1A_064_E053_1_N27_3_E054_1_N27_8_20200604_20200616_VV_INT80_3FBF_foo.tif'
        )

    with pytest.raises(ValueError, match=r'^Failed to parse product name from filename: foo$'):
        _get_product_name_and_type('foo')


def test_add_hyp3_metadata_insar_isce_burst(test_data_dir):
    assert add_hyp3_metadata(
        fname=str(test_data_dir / 'S1_056072_IW2_20220814_20220907_VV_INT80_E09B_corr.tif'),
        meta={
            'WIDTH': 1335,
            'LENGTH': 485,
            'X_STEP': 80.0,
            'Y_STEP': -80.0,
            'X_FIRST': 445520.0,
            'Y_FIRST': 4289840.0
        },
        is_ifg=True,
    ) == {
       'WIDTH': 1335,
       'LENGTH': 485,
       'X_STEP': 80.0,
       'Y_STEP': -80.0,
       'X_FIRST': 445520.0,
       'Y_FIRST': 4289840.0,
       'PROCESSOR': 'hyp3',
       'CENTER_LINE_UTC': '46709.112304',
       'ALOOKS': '4',
       'RLOOKS': '20',
       'EARTH_RADIUS': '6337286.638938101',
       'HEIGHT': '693000.0',
       'STARTING_RANGE': '846099.1914484155',
       'HEADING': -13.310112268820319,
       'ORBIT_DIRECTION': 'ASCENDING',
       'LAT_REF1': '4251040.0',
       'LAT_REF2': '4251040.0',
       'LAT_REF3': '4289840.0',
       'LAT_REF4': '4289840.0',
       'LON_REF1': '445520.0',
       'LON_REF2': '552320.0',
       'LON_REF3': '445520.0',
       'LON_REF4': '552320.0',
       'PLATFORM': 'Sen',
       'ANTENNA_SIDE': -1,
       'WAVELENGTH': 0.055465764662349676,
       'RANGE_PIXEL_SIZE': 46.0,
       'AZIMUTH_PIXEL_SIZE': 56.4,
       'DATE12': '220814-220907',
       'P_BASELINE_TOP_HDR': '158.57439820410497',
       'P_BASELINE_BOTTOM_HDR': '158.57439820410497',
       'beam_mode': 'IW',
       'beam_swath': '2',
       'unwrap_method': 'snaphu_mcf'
   }


def test_add_hyp3_metadata_insar_isce_multi_burst(test_data_dir):
    assert add_hyp3_metadata(
        fname=str(test_data_dir / 'S1_044_000000s1n00-093117s2n01-093118s3n01_IW_20250718_20250730_VV_INT80_B4FA_unw_phase.tif'),
        meta={
            'WIDTH': 2314,
            'LENGTH': 718,
            'X_STEP': 80.0,
            'Y_STEP': -80.0,
            'X_FIRST': 660960.0,
            'Y_FIRST': 5950880.0,
        },
        is_ifg=True,
    ) == {
       'WIDTH': 2314,
       'LENGTH': 718,
       'X_STEP': 80.0,
       'Y_STEP': -80.0,
       'X_FIRST': 660960.0,
       'Y_FIRST': 5950880.0,
       'PROCESSOR': 'hyp3',
       'CENTER_LINE_UTC': '62482.201944',
       'ALOOKS': '4',
       'RLOOKS': '20',
       'EARTH_RADIUS': '6337286.638938101',
       'HEIGHT': '693000.0',
       'STARTING_RANGE': '849199.1730339259',
       'HEADING': -164.275514873119,
       'ORBIT_DIRECTION': 'DESCENDING',
       'LAT_REF1': '5950880.0',
       'LAT_REF2': '5950880.0',
       'LAT_REF3': '5893440.0',
       'LAT_REF4': '5893440.0',
       'LON_REF1': '846080.0',
       'LON_REF2': '660960.0',
       'LON_REF3': '846080.0',
       'LON_REF4': '660960.0',
       'PLATFORM': 'Sen',
       'ANTENNA_SIDE': -1,
       'WAVELENGTH': 0.055465764662349676,
       'RANGE_PIXEL_SIZE': 46.0,
       'AZIMUTH_PIXEL_SIZE': 56.4,
       'DATE12': '250718-250730',
       'P_BASELINE_TOP_HDR': '33.88969537726557',
       'P_BASELINE_BOTTOM_HDR': '33.88969537726557',
       'beam_mode': 'IW',
       'beam_swath': '23',
       'unwrap_method': 'snaphu_mcf',
   }


def test_add_hyp3_metadata_insar_gamma(test_data_dir):
    assert add_hyp3_metadata(
        fname=str(test_data_dir / 'S1AC_20251001T204513_20251007T204359_HHR006_INT40_G_ueF_1DBE_dem.tif'),
        meta={
            'WIDTH': 6829,
            'LENGTH': 3735,
            'X_STEP': 40.0,
            'Y_STEP': -40.0,
            'X_FIRST': 282180.0,
            'Y_FIRST': 6802380.0,
        },
        is_ifg=False,
    ) == {
       'WIDTH': 6829,
       'LENGTH': 3735,
       'X_STEP': 40.0,
       'Y_STEP': -40.0,
       'X_FIRST': 282180.0,
       'Y_FIRST': 6802380.0,
       'PROCESSOR': 'hyp3',
       'CENTER_LINE_UTC': '74714.869576',
       'ALOOKS': '2',
       'RLOOKS': '10',
       'EARTH_RADIUS': '6362008.3766',
       'HEIGHT': '705239.4029000001',
       'STARTING_RANGE': '803591.1723',
       'HEADING': -18.150511100000017,
       'ORBIT_DIRECTION': 'ASCENDING',
       'LAT_REF1': '6652980.0',
       'LAT_REF2': '6652980.0',
       'LAT_REF3': '6802380.0',
       'LAT_REF4': '6802380.0',
       'LON_REF1': '282180.0',
       'LON_REF2': '555340.0',
       'LON_REF3': '282180.0',
       'LON_REF4': '555340.0',
       'PLATFORM': 'Sen',
       'ANTENNA_SIDE': -1,
       'WAVELENGTH': 0.055465764662349676,
       'RANGE_PIXEL_SIZE': 23.0,
       'AZIMUTH_PIXEL_SIZE': 28.2,
       'beam_mode': 'IW',
       'beam_swath': '123',
       'relative_orbit': 90,
       'startUTC': '2025-10-01 20:45:13.000000',
       'stopUTC': '2025-10-01 20:45:40.000000',
       'unwrap_method': 'mcf'
   }
