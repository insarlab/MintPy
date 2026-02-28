import numpy as np
import pytest

from mintpy.utils.utils0 import (
    flatten_for_resample,
    move_spatial_dimension,
    restore_from_resample,
)


def test_flatten_for_resample_2d():
    data = np.random.rand(20, 30)

    out = flatten_for_resample(data)

    assert out.shape == (20, 30)
    np.testing.assert_array_equal(out, data)

def test_flatten_for_resample_3d():
    data = np.random.rand(20, 30, 5)

    out = flatten_for_resample(data)

    assert out.shape == (20, 30, 5)
    np.testing.assert_array_equal(out, data)

def test_flatten_for_resample_4d():
    data = np.random.rand(20, 30, 5, 4)

    out = flatten_for_resample(data)

    assert out.shape == (20, 30, 20)  # 5 * 4
    np.testing.assert_array_equal(out.reshape(20, 30, 5, 4), data)

def test_restore_from_resample_2d():
    rows, cols = 20, 30
    data = np.random.rand(rows, cols)

    out = restore_from_resample(
        data,
        non_spatial_shape=()
    )

    assert out.shape == (rows, cols)
    np.testing.assert_array_equal(out, data)

def test_restore_from_resample_3d():
    rows, cols = 20, 30
    non_spatial_shape = (5,)
    data = np.random.rand(rows, cols, 5)

    out = restore_from_resample(
        data,
        non_spatial_shape=non_spatial_shape
    )

    assert out.shape == (20, 30, 5)
    np.testing.assert_array_equal(out, data)

def test_restore_from_resample_4d():
    rows, cols = 20, 30
    non_spatial_shape = (5, 4)
    data = np.random.rand(rows, cols, 20)  # 5 * 4

    out = restore_from_resample(
        data,
        non_spatial_shape=non_spatial_shape
    )

    assert out.shape == (20, 30, 5, 4)
    np.testing.assert_array_equal(out.reshape(20, 30, 20), data)

@pytest.mark.parametrize(
    "shape",
    [
        (20, 30),           # 2D
        (20, 30, 5),        # 3D
        (20, 30, 5, 4),     # 4D
        (20, 30, 2, 3, 4),  # 5D
    ],
)
def test_flatten_restore_roundtrip(shape):
    data = np.random.rand(*shape)
    non_spatial_shape = shape[2:]

    flat = flatten_for_resample(data)
    restored = restore_from_resample(
        flat,
        non_spatial_shape=non_spatial_shape
    )

    np.testing.assert_array_equal(restored, data)

def test_restore_from_resample_invalid_shape():
    rows, cols = 20, 30
    data = np.random.rand(rows, cols, 10)

    with pytest.raises(ValueError):
        restore_from_resample(
            data,
            non_spatial_shape=(3, 4)  # 12 != 10
        )

def test_move_spatial_dimension_2d_front_and_back():
    # shape: (time, row, col)
    arr = np.zeros((20, 30))

    front = move_spatial_dimension(arr, to_front=True)
    assert front.shape == (20, 30)

    back = move_spatial_dimension(front, to_front=False)
    assert back.shape == arr.shape

    # ensure data integrity
    np.testing.assert_array_equal(back, arr)

def test_move_spatial_dimension_3d_front_and_back():
    # shape: (time, row, col)
    arr = np.zeros((5, 20, 30))

    front = move_spatial_dimension(arr, to_front=True)
    assert front.shape == (20, 30, 5)

    back = move_spatial_dimension(front, to_front=False)
    assert back.shape == arr.shape

    # ensure data integrity
    np.testing.assert_array_equal(back, arr)

def test_move_spatial_dimension_4d_front_and_back():
    # shape: (time, band, row, col)
    arr = np.random.rand(3, 4, 20, 30)

    front = move_spatial_dimension(arr, to_front=True)
    assert front.shape == (20, 30, 3, 4)

    back = move_spatial_dimension(front, to_front=False)
    assert back.shape == arr.shape

    np.testing.assert_array_equal(back, arr)

@pytest.mark.parametrize(
    "shape",
    [
        (20, 30),           # 2D
        (20, 30, 5),        # 3D
        (20, 30, 5, 4),     # 4D
        (20, 30, 2, 3, 4),  # 5D
    ],
)
@pytest.mark.parametrize("to_front", [True, False])
def test_move_spatial_dimension_roundtrip(shape, to_front):
    data = np.random.rand(*shape)

    # Move spatial dimensions
    moved = move_spatial_dimension(data, to_front=to_front)

    # Invert the move
    restored = move_spatial_dimension(moved, to_front=not to_front)

    # Shape should match original
    assert restored.shape == data.shape

    # Values should match original exactly
    np.testing.assert_array_equal(restored, data)
