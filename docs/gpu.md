# Configure GPU acceleration for the network inversion #

The `invert_network` step (in `ifgram_inversion.py`) ships an opt-in GPU solver that batches the per-pixel weighted least-squares inversion as normal-equations + Cholesky on a CUDA device via PyTorch. This is a partial GPU implementation: only `invert_network` is offloaded to the GPU; every other step in `smallbaselineApp.py` continues to run on the CPU. The solver is opt-in — the default `mintpy.networkInversion.solver = auto` resolves to `cpu`, so existing setups are unaffected.

The `torch` solver is orthogonal to Dask parallel processing (see [dask.md](./dask.md)): the former replaces the per-pixel CPU loop with a single batched Cholesky on one CUDA device, the latter distributes that same per-pixel loop across multiple worker processes. The two paths are not currently combined; pick one.

## 1. Setup ##

See [installation.md](./installation.md) section 2.4 for installing the `[gpu]` extras with the matching CUDA wheel index. Selecting `solver = torch` on a host without a visible CUDA device is a hard error (no silent CPU fallback).

## 2. Enable ##

#### 2.1 via command line ####

Run the following in the terminal:

```bash
ifgram_inversion.py inputs/ifgramStack.h5 --solver torch
ifgram_inversion.py inputs/ifgramStack.h5 --solver torch --gpu-chunk-size 20000
```

`--gpu-chunk-size 0` (the default) auto-sizes the per-chunk pixel count from free VRAM; pass a positive integer to override.

#### 2.2 via template file ####

Adjust options in the template file:

```cfg
mintpy.networkInversion.solver       = torch  #[cpu / torch], auto for cpu
mintpy.networkInversion.gpuChunkSize = auto   #[int >= 0], auto for 0 (auto-size from free VRAM)
```

and feed the template file to the script:

```bash
ifgram_inversion.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg
smallbaselineApp.py smallbaselineApp.cfg
```

#### 2.3 Testing using example data ####

Download and run the FernandinaSenDT128 example data; then run with and without the GPU solver:

```bash
cd FernandinaSenDT128/mintpy
ifgram_inversion.py inputs/ifgramStack.h5 -w no --solver cpu
ifgram_inversion.py inputs/ifgramStack.h5 -w no --solver torch
```

The two outputs should agree to float32 round-off (RMS on the order of 1e-5).

## 3. Behavior notes ##

+ **VRAM auto-sizing.** `gpuChunkSize = 0` (auto) probes free VRAM at runtime and chooses a per-chunk pixel count with a fixed headroom factor. Set an explicit integer to override (e.g. for reproducible chunking across hosts with different VRAM).

+ **Rank-deficient pixels.** Detected via `torch.linalg.cholesky_ex` info codes; their solution is set to zero so NaN/Inf never propagate downstream. A warning line reports the count per chunk.

+ **Per-pixel NaN observations.** Handled by zeroing the corresponding row weight, which is mathematically equivalent to dropping that row from the WLS system.

+ **No silent CPU fallback.** Selecting `solver = torch` on a host without a visible CUDA device raises immediately rather than silently falling back to CPU; this keeps performance regressions visible.

## 4. Performance ##

Indicative numbers below were measured on an NVIDIA RTX 5080 (Blackwell sm_120, CUDA 12.8, PyTorch 2.11) at the time this feature was submitted. Speedup depends on scene size, GPU class, and chunk-size tuning, so reproduce on your own data and hardware before drawing conclusions.

+ **Tutorial-scale** (FernandinaSenDT128: 270k pixels, 288 ifgs) — `invert_network` runs roughly **16×** faster internally and **4.5×** faster end-to-end versus the CPU path.
+ **Large-scene** (GalapagosSenDT128: 3.4M pixels, 475 ifgs; ~12.6× pixels and 1.65× ifgs over Fernandina) — roughly **44×** internal and **36×** step-wall speedup on `invert_network` (CPU 6189 s → torch 170 s on the same machine), confirming the speedup grows at scale.
+ **Numerical equivalence** between the `cpu` and `torch` solvers holds to float32 round-off: RMS on the order of `1e-5` on the tutorial case, with absolute RMS at most ~16 µm on the large-scene case.
