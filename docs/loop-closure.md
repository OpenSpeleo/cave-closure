# Measurement-based loop closure

`cave_closure::loop_closure` owns measurement-to-vector conversion, uncertainty
propagation, station substitutions, component controls, floating translations,
sparse graph assembly, convergence and fallback solves. Consumer libraries only
map their data and units to this API and map the result to their output formats.
No consumer or reader dependency enters the Rust implementation.

Rust callers pass an owned `LoopClosureInput` to `prepare_loop_closure`, then
solve the resulting immutable `LoopClosurePreparation`. Python callers use
`prepare_loop_closure(request)` for the same reusable preparation, or
`solve_loop_closure_json(request)` for a single computation boundary. The latter
copies the entire request once, releases the interpreter, prepares and solves in
Rust, then returns the complete serialized `SolveReport`. It includes every
coordinate, residual, status and diagnostic; no intermediate matrix, per-axis
solve or gauge correction returns to Python. Cancellation and Python signals use
the same worker mechanism as the Gaussian API.

## Owned request

The Python request is a mapping with these seven fields:

- `station_ids`: unique nonempty opaque station identifiers in input order.
- `measurements`: rows
  `(id, start_index, stop_index, length_m, bearing_deg, vertical_delta_m, secondary_bearing_deg_or_None, quality, hard_offset)`.
- `initial_positions`: local ENU triples in metres, one for each station.
- `controls`: rows `(id, station_index, axes, values)` for exact controls.
- `relationships`: rows `(id, kind, start_index, stop_index, value)`, using
  `equivalent` or `hard_offset`. An equivalence carries no displacement.
- `source_order`: a permutation of station indexes for deterministic component
  anchor selection.
- `options`: every field of the Rust `LoopClosureOptions` structure, including
  mode, numerical algorithm, weighting, graph method, standard deviations, graph
  policies and convergence settings.

Python row records are tuples; contained axis/value sequences and arrays may be
ordinary lists. Preparation validates shapes, indexes, identities, finite
values, units-dependent uncertainty choices and incompatible configurations.
Inputs are copied and preparations own all their evidence. Concurrent solves use
independent workspaces. `Prepared.normalized_json()` exposes the exact
normalized numerical model, with one entry per independent instrument axis when
applicable.

## Numerical models and sources

`least_squares` supports OpenSpeleo's independent Cartesian Gaussian placement
model and capped inverse-length relative precision. It also supports
`length_azimuth_diagonal`: Stalagmight's independent horizontal axis variances,
with section-quality-scaled length/angle sigmas and optional averaged forward
bearings. Its source is Stalagmight's
`src-tauri/crates/models/src/project/lsq.rs`. Zero-variance axes omit their
observations; floating translations retain the measured centroid. The same
native preparation and linear backends used by the Gaussian API solve these
models. Relative weights do not claim calibrated instrument statistics.

`reference_cg` reproduces the horizontal Java/Rust reference in
`docs/loop_closure/ariane_solver.java` and `graph_solver.rs`, also present in
`compass_lib/loop_closure`. This is the solver inside Ariane's Line: capped
inverse-length precision, unreduced unpreconditioned CG, measured-position warm
starts, an absolute residual threshold and a near-zero denominator exit. It
returns the current iterate at its historical exits, with a diagnostic warning
when the residual tolerance was not reached.

`traverse_quality_cg` reproduces `compass_lib/compass_lib/solver/ariane.py`:
inverse-squared-length precision, anchor-pair BFS traverses blocking
intermediate anchors, best misclosure-per-shot quality per endpoint pair,
centered sparse normal equations, independent XYZ CG, and sparse direct
fallback. Fewer than two merged anchors retains the initial positions. The
source attributes its traverse-quality strategy to Larry Fish. This
implementation does not require Python, NumPy or SciPy. With `sparse-direct`,
fallback uses the library's native sparse factorization support.

The graph methods retain scalar-power, compensated-sum and fused-product
rounding conventions of the original implementation. On Apple platforms,
reference-compatible dot reduction uses the native system BLAS. Existing
consumer GeoJSON goldens remain the byte-level migration contract; the consumer
regression suite reads those files without regenerating them.

## Verification

`tests/loop_closure.rs` contains independent analytical solutions for each
weighting and graph policy, control and identity regressions, malformed-input
checks, cancellation, concurrency, spatial depth behavior and fallback tests.
Native unit tests check anchor-blocked BFS selection, shared traverses, parallel
edges, zero quality, length floors, compensated norms and fused-product
rounding. `tests/core.rs` exercises the underlying Gaussian model and numerical
backends. Both default and all-features Cargo test configurations must pass.

## Inverse-length weighted L1 loop closure

`cave_closure::l1` supplies least-absolute-deviations loop closure. Unlike the
Gaussian models, it minimizes `sum(weight * abs(residual))` independently on X,
Y and Z. The relative weight of a leg is
`1 / max(measured_displacement_length_m, length_floor_m)`. Errors can
concentrate on a few weaker legs instead of being spread as squared residuals.

This is the formulation previously used by Compass's `SparseSolver` in
`compass_lib/solver/sparse.py`, before its Rust migration. That implementation
used `scipy.optimize.linprog`; the library now assembles a sparse column matrix
and calls HiGHS directly from Rust. Its backend dependency is pinned in
`Cargo.toml` and `Cargo.lock`. There is no Python optimizer in the solve path.
The native model retains the source's centered coordinates, free-station order,
positive/negative slack-column order, presolve and dual-simplex settings. All
anchor coordinates stay exact. Failed or interrupted solves raise an error; they
never substitute guessed geographic coordinates.

Rust callers use `prepare_l1_loop_closure(L1Input)` and reuse the immutable
`L1Prepared`. Python exposes `prepare_l1_loop_closure(request)` and
`solve_l1_loop_closure_json(request)`. The latter copies, validates, assembles,
solves every axis and serializes the complete result in one library call. It
releases the interpreter while computing and supports `CancellationToken` and
Python signals through HiGHS interrupt callbacks. Each solve owns its own model;
shared preparations can be solved concurrently.

The Python request has five fields:

- `station_ids`: unique nonempty opaque identifiers.
- `observations`: tuples `(id, start_index, stop_index, displacement_or_None)`.
  Displacements are ENU triples in metres. `None` explicitly requests the
  initial endpoint-coordinate difference.
- `initial_positions`: owned ENU triples in metres.
- `anchor_indices`: stations fixed at their initial coordinates.
- `options`: `mode`, `length_floor_m`, `minimum_anchor_count`,
  `deduplicate_edges`, and `max_iterations` (an integer or `None`).

Every option is explicit. The Compass compatibility configuration uses spatial
mode, a 0.1 metre length floor, a minimum of two anchors, first-observation
undirected-edge deduplication, and no iteration limit. With fewer anchors,
initial positions are retained. Horizontal mode adjusts XY and preserves each
initial Z. Deduplication can be disabled to retain parallel observations. The
source dense-assignment behavior of self-observations is retained for
compatibility. This graph API takes displacement observations; it does not infer
raw instrument calibration or protected-shot relationships. Consumers must
choose the graph model matching their format semantics.

The objective is an uncalibrated L1 objective. Gaussian variance, redundancy and
backward-error statistics are not reported for it; zero placeholders in those
fields are accompanied by an explicit diagnostic warning. Native builds enable
this model with the `l1` Cargo feature; Python builds include that feature by
default. Building the backend requires CMake, a C++ compiler and libclang for
the C bindings. Wheels contain the compiled backend.

`tests/loop_closure_l1.rs` checks analytical optima distinct from least squares,
XYZ controls, depth preservation, reversed and repeated legs, zero lengths,
source-compatible self-observations, malformed inputs, cancellation and
concurrent reuse. `tests/test_loop_closure_l1.py` checks the owned Python API,
single-call/prepared agreement, exception types and immutable source evidence.
