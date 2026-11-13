# Numerical contract

For each oriented displacement, residual is endpoint minus origin minus measured
displacement. The Gaussian objective is the sum of residual quadratic forms in
supplied precision, including soft controls. Observations remain distinct,
including parallels, reversals, self-loops, and fixed–fixed factors. Exact
relationships do not remove the residual contributions of contracted factors.

Covariance is finite, symmetric, and positive definite. Horizontal solving
marginalizes with the XY covariance, before inversion. First-order polar and
underwater Jacobian helpers use radians and positive-down depths. The underwater
helper rejects impossible or caller-defined unreliable near-vertical geometry.
These helpers do not track correlations themselves: normalized factors carry raw
reading identities. Shared readings are rejected unless an explicit independent
factor approximation is requested; that approximation disables calibrated statistics.

The reduced design has full column rank under the supported full-dimensional
interstation factors and axis-subset unary controls after reported gauges.
Stochastic redundancy is scalar row count minus free scalar unknown count. The
variance factor is absent at zero redundancy or with relative weights. Calibrated
statistics eligibility additionally assumes the caller's covariance and Gaussian
model are valid. No residual is a probability of blunder. Series chains can have
indistinguishable standardized residuals; controls create redundancy even on trees.

Every solve validates the reconstructed normal-system residual and computes
original observation residuals, objective, hard-constraint discrepancy, actual
kernel size, iterations, and scaled infinity-norm backward error. The default
backward error target is numerical, not an accuracy guarantee. Poor conditioning,
large coordinate offsets, and floating-point cancellation can limit accuracy.
Unsafe elimination pivots remain in the kernel. Nonpositive curvature, overflow,
iteration exhaustion, stagnation, and cancellation never receive a converged
status. Strict convenience calls reject incomplete reports.

Reduction uses full small blocks and Cholesky solves. Off-diagonal Schur blocks
are not assumed symmetric. Degree-limited elimination has linear-size storage
for bounded dimension; sorting adds overhead. PCG costs iterations times sparse
system size. Direct factorization can incur substantial fill. No performance
budget or cross-platform bit identity is promised. Deterministic ID ordering
makes input permutation immaterial within a build.

Huber, general correlated and latent factors, nonlinear raw-reading adjustment,
sequential fundamental-cycle corrections, and influence diagnostics are outside
the initial released capability. There is no geological prior. Magnetic models,
CRS transformations, epochs, and vertical datums belong to adapters. Adapters
record transformation provenance and must apply declination and convergence once.

The Python generator `tests/generate_graphs.py` creates seeded independent
multigraphs at configurable sizes, component counts, and extra-edge counts.
It emits chains, rings, and random spanning trees with additional chords, full
noncommuting covariance blocks, parallel/reversed observations, self-loops, and
multiple fixed stations. Each reading identity is distinct. Small cases use an
independent pivoted dense QR reference; larger cases can use independent sparse
LSMR of the whitened rectangular system. Reference nonconvergence is an error.
Generated measurements and reference values never come from the production solver.

```sh
python tests/generate_graphs.py --stations 48 2500 --components 5 \
  --seeds 3 17 101 --topologies random ring --reference auto \
  --solver-results --output-dir tests/fixtures
```

The persistent corpus includes independent references and results from every
implemented linear backend. See [the fixture contract](fixtures.md) for generation,
coverage, numerical tolerances, and tests that check fresh and saved results.

The original-system stopping test uses the true residual, not a recursively
updated PCG residual. Its absolute-plus-relative scale is
`||H||inf * ||q||inf + ||g||inf`. This accounts for cancellation when large
coordinates satisfy a small RHS. The separate backward-error target remains
explicit. The kernel is solved more tightly before reconstruction; up to three
original-system refinement corrections are available by default. Refinement uses
the same reduction and PCG backend, reports its step count, and respects the total
iteration budget. A small backward residual still does not establish coordinate
accuracy in a poorly conditioned survey.

The initial local macOS arm64 release-build benchmark used five repetitions per
case. These are observations for the generated chain/ring model, not release
gates or a comparison against a different solver's workload:

| Stations | Topology | Median preparation (ms) | Median solve (ms) |
| -------: | -------- | ----------------------: | ----------------: |
|    1,000 | Chain    |                   0.806 |             2.575 |
|    1,000 | Ring     |                   0.544 |             2.080 |
|   10,000 | Chain    |                   4.364 |            25.567 |
|   10,000 | Ring     |                   4.221 |            14.861 |
|   24,000 | Chain    |                  10.282 |            51.918 |
|   24,000 | Ring     |                  10.040 |            35.091 |

Run `cargo bench --locked --bench networks` to reproduce the workload. The
benchmark checks the analytic objective on every repetition. Use
`tests/benchmark_graphs.py` for generated chorded random cases, including separate
JSON preparation/solve timing and independent-reference objective and coordinate
differences. Debug-wheel timings must not be compared with release benchmarks.

`coordinate_buffer(result)` exports a read-only owned copy with three native
doubles per sorted station. Its memoryview retains a bytes owner independently
of result lifetime. This interface deliberately makes no zero-copy promise.
