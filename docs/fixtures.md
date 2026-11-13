# Persisted numerical corpus

`tests/fixtures/corpus.json` names every input, independent reference, and solver
result. The corpus contains the analytic controlled tree and twelve generated
multigraphs: three seeds for both random and ring topologies at 48 and 2,500
stations. Each generated graph has five components and two fixed stations per
component, full noncommuting covariance blocks, parallel and reversed readings,
self-loops, and additional chords. The larger cases contain over 7,000 factors.

Each input has a `.reference.json` containing coordinates, original observation
residuals, objective, stationarity, and scalar redundancy. Small references use
SciPy's pivoted dense QR; larger references use sparse LSMR on an independently
assembled whitened rectangular design. The offline `auto` reference policy uses
dense QR for at most 600 unknowns to bound reference-tool memory and records the
chosen method. It does not select a production backend.

These references describe the noisy problem's unique least-squares optimum.
The original latent coordinates used to sample noisy measurements generally do
not minimize that noisy objective and are not treated as the expected solution.
The production solver never generates its own reference coordinates.

Four `.result.json` files accompany every input: reduced PCG, unreduced PCG,
sparse QR, and sparse Cholesky. Each envelope names the solver and reduction
setting and contains the complete finite, converged `SolveReport`. All four solve
the same Gaussian objective. The unreduced PCG configuration supplies a further
check on reduction and reconstruction.

Run from the standalone source root to regenerate the inputs and results:

```sh
python tests/generate_graphs.py --stations 48 2500 --components 5 \
  --seeds 3 17 101 --topologies random ring --reference auto \
  --solver-results --output-dir tests/fixtures
```

The generator retains the analytic input and includes every existing spatial
input in the destination when rebuilding references and the manifest. The current
independent reference tool supports full XYZ hard controls and displacement
covariances or isotropic uncertainties. It rejects models outside that scope.
Generating solver results requires Cargo and builds the `fixture_results`
example with `sparse-direct`. Generation fails if any solve is incomplete or
disagrees with its independent reference.

```sh
python -m pytest tests/test_fixture_corpus.py -q
```

Tests verify saved results and freshly computed results against the reference,
then compare fresh and saved coordinates. They check fixed controls, all original
residuals and corrections, the independently recomputed objective, rank-based
redundancy, and normal-gradient stationarity. Fresh native results go to the test
runner's temporary directory. Both PCG configurations are also exercised through
the installed Python interface. Tests never regenerate persistent expectations.

Coordinate comparisons use `1e-9 m + 1e-10 * S`, where `S` is the maximum absolute
reference coordinate. Objective comparisons use `1e-12 + 1e-9 * abs(Q)`. These
scale-aware tolerances permit floating-point differences across backends and
platforms; iteration counts and report bytes need not be identical. Regression
tests deliberately corrupt coordinates, residuals, and objective values to prove
that incorrect saved results fail verification.

The sdist includes every manifest member. Artifact inspection rejects a missing
input, reference, or solver result. The Rust crate retains its compact analytic
unit fixture; the complete multi-backend integration corpus belongs to the
standalone source distribution and uses that crate's packaged native example.
