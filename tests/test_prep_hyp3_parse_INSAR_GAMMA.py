from mintpy.prep_hyp3 import _get_product_name_and_type


def test_get_product_name_and_type_gamma_3digit_days():
    fname = "S1CC_20250705T053629_20250717T053629_VVP012_INT80_G_weF_F288_unw_phase_clip.tif"
    product_name, job_type = _get_product_name_and_type(fname)

    assert product_name == "S1CC_20250705T053629_20250717T053629_VVP012_INT80_G_weF_F288"
    assert job_type == "INSAR_GAMMA"


def test_get_product_name_and_type_gamma_4digit_days():
    fname = "S1BC_20210807T053645_20250729T053630_VVP1452_INT80_G_weF_9C74_unw_phase_clip.tif"
    product_name, job_type = _get_product_name_and_type(fname)

    assert product_name == "S1BC_20210807T053645_20250729T053630_VVP1452_INT80_G_weF_9C74"
    assert job_type == "INSAR_GAMMA"
