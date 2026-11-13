# cave_closure

Independent Gaussian weighted least-squares adjustment of oriented survey
multigraphs in local ENU metres. The Rust core has no consumer, geographic,
Tauri, or Python dependency. This source tree is standalone; upstream and
registry identities must be supplied before publication or submodule
registration.

The caller selects spatial or marginal-horizontal adjustment. Horizontal mode
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
and covariance buffers. NumPy is optional. Preparation copies input; exporters
must not be modified concurrently while the copy runs. Prepared objects and
results are immutable. Concurrent solves have separate workspaces and may share
a preparation. `CancellationToken` supports worker cancellation; main-thread
Python signals retain their exception type. Direct factorization observes
cancellation before and after its external call.

Run checks directly from this tree:

```sh
cargo test --locked --all-features
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
uv run maturin build --features abi3 --out dist
uv run pytest
```

See [the numerical contract](docs/numerics.md),
[interchange contract](docs/interchange.md), and
[artifact qualification](docs/artifacts.md). The
[persisted fixture corpus](docs/fixtures.md) checks every solver against
independent references. AGPL terms are in [LICENSE](LICENSE).
