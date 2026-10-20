"""Numerical-equivalence tests for src/mintpy/ifgram_inversion_gpu.py.

Compare the GPU-batched solver against the per-pixel CPU reference
(scipy.linalg.lstsq via mintpy.ifgram_inversion.estimate_timeseries) on
synthetic SBAS-like networks. Tests are skipped automatically when
PyTorch / CUDA are unavailable.
"""

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from mintpy.ifgram_inversion import estimate_timeseries
from mintpy.ifgram_inversion_gpu import estimate_timeseries_batch

requires_cuda = pytest.mark.skipif(
    not torch.cuda.is_available(),
    reason="CUDA-capable GPU required for ifgram_inversion_gpu tests",
)


def make_redundant_network(num_date, num_pair, *, max_span=4, seed=0):
    """Synthetic SBAS-like network with multiple ifgrams per adjacent date
    interval.

    All adjacent pairs (i, i+1) are forced first to guarantee minimal
    connectivity, then short-baseline pairs (i, j) with j - i <= max_span
    are sampled until ``num_pair`` is reached. This mirrors a well-connected
    Sentinel-1 SBAS network and stays full-rank under partial NaN masking
    of individual interferograms.
    """
    rng = np.random.default_rng(seed)
    candidate = [(i, j) for i in range(num_date)
                        for j in range(i + 1, min(i + max_span + 1, num_date))]
    if num_pair > len(candidate):
        raise ValueError(
            f"requested num_pair={num_pair} exceeds {len(candidate)} candidates "
            f"for num_date={num_date}, max_span={max_span}"
        )
    forced = [(i, i + 1) for i in range(num_date - 1)]
    forced_set = set(forced)
    optional = [p for p in candidate if p not in forced_set]
    rng.shuffle(optional)
    selected = sorted(forced + optional[: num_pair - len(forced)])

    A = np.zeros((num_pair, num_date - 1), dtype=np.float32)
    for k, (i, j) in enumerate(selected):
        A[k, i:j] = 1.0
    tbase = np.cumsum(rng.uniform(0.05, 0.2, size=num_date - 1).astype(np.float32))
    tbase = np.concatenate([[0.0], tbase])
    tbase_diff = np.diff(tbase).reshape(-1, 1).astype(np.float32)
    B = A * tbase_diff.T
    return A, B, tbase_diff


def synthesize_observations(A, B, num_pixel, *, nan_frac=0.0, weighted=True, seed=0):
    rng = np.random.default_rng(seed)
    num_pair, num_unknown = A.shape
    X_true = rng.normal(0, 1, size=(num_unknown, num_pixel)).astype(np.float32)
    y = (B @ X_true) + rng.normal(0, 0.05, size=(num_pair, num_pixel)).astype(np.float32)
    weight_sqrt = (rng.uniform(0.5, 2.0, size=(num_pair, num_pixel)).astype(np.float32)
                   if weighted else None)
    if nan_frac > 0:
        nan_mask = rng.random(y.shape) < nan_frac
        y[nan_mask] = np.nan
    return y, weight_sqrt


def cpu_reference(A, B, y, weight_sqrt, tbase_diff, *, min_norm_velocity=True):
    """Per-pixel CPU reference via mintpy.ifgram_inversion.estimate_timeseries."""
    num_pair, num_unknown = A.shape
    num_pixel = y.shape[1]
    num_date = num_unknown + 1
    ts = np.zeros((num_date, num_pixel), dtype=np.float32)
    tcoh = np.zeros(num_pixel, dtype=np.float32)
    nobs = np.zeros(num_pixel, dtype=np.int16)
    for k in range(num_pixel):
        wk = None if weight_sqrt is None else weight_sqrt[:, k]
        tsi, qi, ni = estimate_timeseries(
            A=A, B=B, y=y[:, k], tbase_diff=tbase_diff,
            weight_sqrt=wk,
            min_norm_velocity=min_norm_velocity,
            rcond=1e-5, min_redundancy=1.0,
            inv_quality_name='temporalCoherence',
            print_msg=False,
        )
        ts[:, k] = tsi.flatten()
        tcoh[k] = qi
        nobs[k] = ni
    return ts, tcoh, nobs


def assert_equivalent(cpu, gpu, *, ts_rel_tol, tcoh_abs_tol):
    ts_cpu, tcoh_cpu, nobs_cpu = cpu
    ts_gpu, tcoh_gpu, nobs_gpu = gpu

    assert np.all(np.isfinite(ts_gpu)), 'GPU returned non-finite timeseries'
    assert np.all(np.isfinite(tcoh_gpu)), 'GPU returned non-finite temporalCoh'

    ts_scale = max(float(np.abs(ts_cpu).max()), 1e-6)
    ts_rms = float(np.sqrt(np.mean((ts_cpu - ts_gpu) ** 2)))
    assert ts_rms < ts_rel_tol * ts_scale, (
        f'timeseries rms={ts_rms:.3e} > {ts_rel_tol} * scale={ts_scale:.3f}'
    )

    tcoh_rms = float(np.sqrt(np.mean((tcoh_cpu - tcoh_gpu) ** 2)))
    assert tcoh_rms < tcoh_abs_tol, (
        f'temporalCoherence rms={tcoh_rms:.3e} > {tcoh_abs_tol:.3e}'
    )

    np.testing.assert_array_equal(nobs_cpu, nobs_gpu)


@pytest.fixture(scope='module')
def network():
    """Shared SBAS-like network sized to match FernandinaSenDT128
    (num_date=98, num_pair=288). max_span=6 reproduces the typical
    short-baseline coverage of a Sentinel-1 SBAS network and provides
    enough redundancy that small fractions of NaN observations do not
    push individual pixels into rank-deficiency.
    """
    return make_redundant_network(num_date=98, num_pair=288, max_span=6, seed=0)


def _all_pixels_full_rank(A, y):
    """Return True if every pixel's design matrix (after dropping NaN rows)
    is still full-rank. Used to keep tests off the rank-deficient edge case
    (which is handled separately at runtime by a CPU fallback path).
    """
    n_unknown = A.shape[1]
    for k in range(y.shape[1]):
        valid = ~np.isnan(y[:, k])
        if np.linalg.matrix_rank(A[valid]) < n_unknown:
            return False
    return True


@requires_cuda
def test_wls_no_nan(network):
    """WLS, no NaN observations — expect ~ float32 round-off match."""
    A, B, tbase_diff = network
    y, w = synthesize_observations(A, B, num_pixel=64, nan_frac=0.0, seed=1)
    cpu = cpu_reference(A, B, y, w, tbase_diff)
    gpu = estimate_timeseries_batch(
        A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=w,
        min_norm_velocity=True,
        chunk_size=64, solver='torch', print_msg=False,
    )
    assert_equivalent(cpu, gpu, ts_rel_tol=1e-5, tcoh_abs_tol=1e-5)


@requires_cuda
def test_wls_with_nan_redundant(network):
    """WLS with low NaN rate on a redundant network — float32 round-off match.

    NaN fraction is kept at 3% so each pixel still has enough observations to
    keep its (NaN-masked) design matrix full-rank. Higher NaN rates exercise
    the rank-deficient edge case, which is handled by a separate CPU fallback
    path (not yet implemented at the time of writing).
    """
    A, B, tbase_diff = network
    y, w = synthesize_observations(A, B, num_pixel=64, nan_frac=0.03, seed=2)
    assert _all_pixels_full_rank(A, y), \
        'fixture broke: a pixel is rank-deficient after NaN masking'
    cpu = cpu_reference(A, B, y, w, tbase_diff)
    gpu = estimate_timeseries_batch(
        A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=w,
        min_norm_velocity=True,
        chunk_size=32, solver='torch', print_msg=False,
    )
    assert_equivalent(cpu, gpu, ts_rel_tol=1e-4, tcoh_abs_tol=1e-4)


@requires_cuda
def test_ols_no_nan(network):
    """OLS path (weight_sqrt=None)."""
    A, B, tbase_diff = network
    y, _ = synthesize_observations(A, B, num_pixel=64, nan_frac=0.0,
                                   weighted=False, seed=3)
    cpu = cpu_reference(A, B, y, None, tbase_diff)
    gpu = estimate_timeseries_batch(
        A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=None,
        min_norm_velocity=True,
        chunk_size=32, solver='torch', print_msg=False,
    )
    assert_equivalent(cpu, gpu, ts_rel_tol=1e-5, tcoh_abs_tol=1e-5)


@requires_cuda
def test_min_norm_phase(network):
    """min_norm_velocity=False solves on A directly; ts[1:] = X."""
    A, B, tbase_diff = network
    y, w = synthesize_observations(A, B, num_pixel=64, nan_frac=0.0, seed=4)
    cpu = cpu_reference(A, B, y, w, tbase_diff, min_norm_velocity=False)
    gpu = estimate_timeseries_batch(
        A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=w,
        min_norm_velocity=False,
        chunk_size=32, solver='torch', print_msg=False,
    )
    assert_equivalent(cpu, gpu, ts_rel_tol=1e-5, tcoh_abs_tol=1e-5)


@requires_cuda
def test_chunk_size_invariance(network):
    """Output must be effectively identical across chunk sizes."""
    A, B, tbase_diff = network
    y, w = synthesize_observations(A, B, num_pixel=64, nan_frac=0.03, seed=5)
    assert _all_pixels_full_rank(A, y)
    common = dict(A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=w,
                  min_norm_velocity=True, solver='torch', print_msg=False)
    ts_a, tcoh_a, nobs_a = estimate_timeseries_batch(chunk_size=16, **common)
    ts_b, tcoh_b, nobs_b = estimate_timeseries_batch(chunk_size=64, **common)
    np.testing.assert_allclose(ts_a, ts_b, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(tcoh_a, tcoh_b, rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(nobs_a, nobs_b)


@requires_cuda
def test_unsupported_solver_raises(network):
    A, B, tbase_diff = network
    y, w = synthesize_observations(A, B, num_pixel=8, seed=6)
    with pytest.raises(ValueError, match='unsupported solver'):
        estimate_timeseries_batch(
            A=A, B=B, y=y, tbase_diff=tbase_diff, weight_sqrt=w,
            solver='cupy', print_msg=False,
        )
