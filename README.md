# cave-closure

Independent Gaussian weighted least-squares correction of oriented survey
multigraphs in local ENU metres. The Rust core has no consumer, geographic,
Tauri, or Python dependency. This source tree is standalone; upstream and
registry identities must be supplied before publication or submodule
registration.

Inverse-length weighted L1 (least-absolute-deviations) loop closure is also
available through the owned measurement API; see
[measurement-based loop closure](docs/loop-closure.md). All numerical work runs
in the native library, including sparse LP assembly and the HiGHS solve.

The caller selects spatial or marginal-horizontal loop closure. Horizontal mode
requires an explicit preserved-Z value for every station. Covariance models
include supplied SPD matrices, equal isotropic variance, length-proportional
variance, and uncalibrated relative weights. No variance floor is inserted.

Exact equivalences and offsets are contracted before solving. Hard controls may
select individual axes. Soft controls identify translation directions without
being replaced by artificial gauges. Uncontrolled directions receive reported,
deterministic gauges unless `RequireControl` behavior is requested.

The default backend eliminates unknown station blocks with at most two
neighbors, solves the remaining system with block-Jacobi PCG, and reconstructs
coordinates. The optional `sparse-direct` feature adds maintained sparse QR of
the whitened design and sparse Cholesky of the normal system. Backends are
explicit choices. There is no automatic fallback or regularization.

```rust
use cave_closure::{json, prepare, PrepareOptions, SolveConfig, Workspace, SolveControl};
let problem = json::from_json(include_str!("tests/fixtures/controlled_tree.json"))?;
let prepared = prepare(&problem, &PrepareOptions::default())?;
let result = prepared.solve(&SolveConfig::default(), &mut Workspace::default(), &SolveControl::default())?;
result.require_converged()?;
# Ok::<(), cave_closure::Error>(())
```

```python
from cave_closure import prepare_json, coordinates

prepared = prepare_json(input_json)
result = prepared.solve()
result.require_converged()
positions = coordinates(result)
```

`prepare_arrays` is the ordinary owned bulk API. `prepare_buffers` accepts
C-contiguous native unsigned 64-bit endpoint indexes and float64 displacement
and covariance buffers. Both take `mode=LoopClosureMode.SPATIAL` or
`mode=LoopClosureMode.HORIZONTAL`, with `LoopClosureMode` imported from
`cave_closure`. Horizontal mode also requires `preserved_z`. `LoopClosureMode`
is a `StrEnum` with values `"spatial"` and `"horizontal"`; existing string
arguments remain accepted. NumPy is optional. Preparation copies input;
exporters must not be modified concurrently while the copy runs. Prepared
objects and results are immutable. Concurrent solves have separate workspaces
and may share a preparation. `CancellationToken` supports worker cancellation;
main-thread Python signals retain their exception type. Direct factorization
observes cancellation before and after its external call.

Both bulk APIs accept `algorithm=LoopClosureAlgorithm.LSE_REDUCED_PCG` (the
default), `LSE_UNREDUCED_PCG`, `LSE_SPARSE_QR`, or `LSE_SPARSE_CHOLESKY`. Import
`LoopClosureAlgorithm` from `cave_closure`; it is a Rust enum exposed through
PyO3. Standard Python builds include all four algorithms. Spatial and horizontal
modes select loop-closure dimensions independently of the algorithm.

The repository uses one Cargo package with Rust sources in `src_rust/` and
Python sources in `src_python/`. The optional `python` feature enables the
bindings; Maturin enables it and all numerical backends for Python builds. Rust
consumers use the same library without enabling Python.

CI and release use the same native platform matrix in
`.github/workflows/workflow.yml` and `.github/workflows/release.yml`: Linux,
Windows, and macOS on both x86-64 and ARM64. CI tests every configured ordinary
and free-threaded Python interpreter. Release builds the corresponding wheel ABI
variants. Builds run directly on those runners, without containers or
cross-compilation.

Release wheels use three ABI paths: `abi3` for ordinary CPython, `abi3t` for
free-threaded interpreters supporting the stable ABI, and an
interpreter-specific legacy free-threaded wheel. Build floors are selected by
`MINIMUM_PYTHON_VERSION` and `MINIMUM_FT_PYTHON_VERSION` in the release
workflow; the corresponding ABI floors are declared in `Cargo.toml`. Every wheel
enables both numerical backends, `sparse-direct` and `l1`. Each platform builds
one ordinary stable-ABI wheel at the minimum interpreter, covering compatible
ordinary Python versions without rebuilding for each interpreter. See
[artifact qualification](docs/artifacts.md) for the runner matrix and Linux
compatibility limits.

Local checks, CI, and release builds use the exact Rust toolchain declared in
`rust-toolchain.toml`, including its Clippy and rustfmt components. Run commands
from this repository directory so rustup selects that pin. Update the toolchain
file deliberately and validate both Clippy configurations when changing it;
`rustup update` alone must not change the version used by this project. CI logs
the selected compiler and Clippy versions before validation.

Run checks directly from this tree:

```sh
cargo test --locked --all-features
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
uv run maturin build --locked --release --features abi3,sparse-direct,l1,pyo3/extension-module --out dist
uv run pytest
uv run --frozen python scripts/generate_stubs.py --check
```

Native Python stubs in `src_python/cave_closure/_native/__init__.pyi` are
generated from Rust using `pyo3-stub-gen`. Run
`uv run --frozen python scripts/generate_stubs.py` after changing a binding to
regenerate and format them. Add `--check` to compare without modifying files.
The `generate-stubs` pre-commit hook runs this check, including in CI.
`src_python/cave_closure/__init__.pyi` re-exports the generated API and
describes the Python helpers.

The Rust-only `horizontal` API also accepts ordered, precontracted groups and
interpreted XY observations for scalar warm-start compatibility. Its numerical
contract and caller responsibilities are documented with the measurement API.

See [measurement-based loop closure](docs/loop-closure.md),
[the numerical contract](docs/numerics.md),
[interchange contract](docs/interchange.md), and
[artifact qualification](docs/artifacts.md). The
[persisted fixture corpus](docs/fixtures.md) checks every solver against
independent references. AGPL terms are in [LICENSE](LICENSE).
