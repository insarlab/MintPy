import numpy as np
import pytest

from mintpy.utils.utils0 import median_abs_deviation


@pytest.mark.parametrize('data', [
    [[0., 1., 2.], [10., 11., 12.]],
    [[0., 2.], [10., 12.]],
    [[0., np.nan, 2.], [10., 12., np.nan]],
], ids=['rectangular', 'square', 'nan'])
@pytest.mark.parametrize('center, expected', [
    ([1., 11.], [1., 1.]),
    ([-1., 14.], [2., 3.]),
], ids=['median', 'nonmedian'])
@pytest.mark.parametrize('scale', [1., 0.67449])
def test_median_abs_deviation_row_centers(data, center, expected, scale):
    actual = median_abs_deviation(data, center=np.array(center), scale=scale)
    np.testing.assert_allclose(actual, np.array(expected) / scale)


@pytest.mark.parametrize('center, expected', [
    (None, [1., 1.]),
    (0., [1., 11.]),
    (np.array(0.), [1., 11.]),
])
def test_median_abs_deviation_default_and_scalar_centers(center, expected):
    data = np.array([[0., 1., 2.], [10., 11., 12.]])
    actual = median_abs_deviation(data, center=center, scale=1.)
    np.testing.assert_allclose(actual, expected)


@pytest.mark.parametrize('center, expected', [
    (None, 1.),
    (1., 1.),
    (np.array(1.), 1.),
    (-1., 2.),
    (np.array(-1.), 2.),
])
def test_median_abs_deviation_one_dimensional(center, expected):
    actual = median_abs_deviation([0., 1., 2.], center=center, scale=1.)
    np.testing.assert_allclose(actual, expected)
