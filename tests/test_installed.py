"""Common installed-artifact corpus; all reference systems are assembled independently."""

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest
from cave_closure import (
    CancellationToken,
    LoopClosureAlgorithm,
    LoopClosureMode,
    coordinate_buffer,
    coordinates,
    prepare_arrays,
    prepare_buffers,
    prepare_json,
)

ALGORITHMS = (
    LoopClosureAlgorithm.LSE_REDUCED_PCG,
    LoopClosureAlgorithm.LSE_UNREDUCED_PCG,
    LoopClosureAlgorithm.LSE_SPARSE_QR,
    LoopClosureAlgorithm.LSE_SPARSE_CHOLESKY,
)

FIXTURE = Path(__file__).parent / "fixtures" / "controlled_tree.json"


def test_analytic_fixture():
    result = prepare_json(FIXTURE.read_text()).solve()
    result.require_converged()
    assert coordinates(result)["B"] == pytest.approx((10.2, 0, 0), abs=1e-12)
    assert json.loads(result.to_json())["diagnostics"]["objective"] == pytest.approx(
        0.2
    )
    with pytest.raises(TypeError):
        coordinates(result)["B"] = (0, 0, 0)
    with pytest.raises(AttributeError):
        result.status = "Converged"


def arrays():
    return dict(
        station_ids=["A", "B", "C"],
        observation_ids=["ab", "bc"],
        endpoints=np.array([[0, 1], [1, 2]], dtype=np.uint64),
        displacements=np.array([[10, 0, 0], [10, 0, 0]], dtype=np.float64),
        covariances=np.array([np.eye(3), 4 * np.eye(3)]),
        mode=LoopClosureMode.SPATIAL,
        controls=[
            ("a", 0, [0, 1, 2], [0.0, 0.0, 0.0], None),
            ("c", 2, [0, 1, 2], [21.0, 0.0, 0.0], None),
        ],
    )


@pytest.mark.parametrize("prepare", [prepare_arrays, prepare_buffers])
@pytest.mark.parametrize("mode", list(LoopClosureMode))
@pytest.mark.parametrize("algorithm", ALGORITHMS)
@pytest.mark.parametrize("as_string", [False, True])
def test_adjustment_modes(prepare, mode, algorithm, as_string):
    args = arrays()
    args["mode"] = mode.value if as_string else mode
    args["algorithm"] = algorithm
    args["controls"] = []
    if mode is LoopClosureMode.HORIZONTAL:
        args["preserved_z"] = [5.0, 7.0, 9.0]
    result = prepare(**args).solve()
    result.require_converged()
    positions = coordinates(result)
    assert positions["B"][:2] == pytest.approx((10.0, 0.0))
    expected_z = (
        (5.0, 7.0, 9.0) if mode is LoopClosureMode.HORIZONTAL else (0.0, 0.0, 0.0)
    )
    assert tuple(positions[station][2] for station in ["A", "B", "C"]) == expected_z


def test_buffers_are_owned():
    args = arrays()
    prepared = prepare_buffers(**args)
    args["displacements"][:] = 999
    args["covariances"][:] = 0
    del args
    result = prepared.solve()
    del prepared
    assert coordinates(result)["B"][0] == pytest.approx(10.2)


@pytest.mark.parametrize("prepare", [prepare_arrays, prepare_buffers])
@pytest.mark.parametrize("algorithm", ALGORITHMS)
def test_algorithm_selects_solver(prepare, algorithm):
    result = prepare(**arrays(), algorithm=algorithm).solve()
    result.require_converged()
    diagnostics = json.loads(result.to_json())["diagnostics"]
    reduced = algorithm == LoopClosureAlgorithm.LSE_REDUCED_PCG
    iterative = algorithm == LoopClosureAlgorithm.LSE_UNREDUCED_PCG
    assert diagnostics["kernel_unknowns"] == (0 if reduced else 3)
    assert diagnostics["iterations"] == (1 if iterative else 0)
    assert coordinates(result)["B"][0] == pytest.approx(10.2)


@pytest.mark.parametrize("prepare", [prepare_arrays, prepare_buffers])
def test_algorithm_requires_native_enum(prepare):
    with pytest.raises(TypeError):
        prepare(**arrays(), algorithm="lse_sparse_qr")


@pytest.mark.parametrize(
    "field,mutation",
    [
        ("endpoints", lambda a: a.astype(np.int64)),
        ("displacements", lambda a: a.astype(np.float32)),
        ("displacements", lambda a: a[:, ::-1]),
        ("displacements", lambda a: a.byteswap().view(a.dtype.newbyteorder())),
        ("covariances", lambda a: a.reshape(2, 9)),
        ("endpoints", lambda a: np.array([[0, 3], [1, 2]], dtype=np.uint64)),
    ],
)
def test_bad_buffers(field, mutation):
    args = arrays()
    args[field] = mutation(args[field])
    with pytest.raises((ValueError, TypeError, OverflowError)):
        prepare_buffers(**args)


def test_concurrent_and_cancelled_solves():
    prepared = prepare_json(FIXTURE.read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(
            pool.map(lambda _: coordinates(prepared.solve())["B"], range(16))
        )
    assert all(v == results[0] for v in results)
    token = CancellationToken()
    token.cancel()
    with pytest.raises(RuntimeError):
        prepared.solve(cancellation=token)


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("algorithm", ALGORITHMS)
def test_random_full_covariance_against_independent_dense_qr(seed, algorithm):
    rng = np.random.default_rng(seed)
    n = 12
    edges = [(i, i + 1) for i in range(n - 1)] + [
        tuple(rng.choice(n, 2, replace=False)) for _ in range(12)
    ]
    observed = rng.normal(size=(len(edges), 3))
    matrices = rng.normal(size=(len(edges), 3, 3))
    covariance = matrices @ matrices.transpose(0, 2, 1) + 0.5 * np.eye(3)
    rows, rhs = [], []
    for (u, v), delta, cov in zip(edges, observed, covariance):
        a = np.zeros((3, 3 * (n - 1)))
        if u:
            a[:, 3 * (u - 1) : 3 * u] -= np.eye(3)
        if v:
            a[:, 3 * (v - 1) : 3 * v] += np.eye(3)
        root = np.linalg.cholesky(cov)
        rows.append(np.linalg.solve(root, a))
        rhs.append(np.linalg.solve(root, delta))
    a, y = np.vstack(rows), np.concatenate(rhs)
    optimum = np.linalg.lstsq(a, y, rcond=None)[0]
    result = prepare_arrays(
        [str(i) for i in range(n)],
        [str(i) for i in range(len(edges))],
        edges,
        observed.tolist(),
        covariance.tolist(),
        mode=LoopClosureMode.SPATIAL,
        algorithm=algorithm,
        controls=[("anchor", 0, [0, 1, 2], [0, 0, 0], None)],
    ).solve()
    result.require_converged()
    actual = coordinates(result)
    values = np.array([actual[str(i)] for i in range(1, n)]).ravel()
    np.testing.assert_allclose(values, optimum, atol=1e-9, rtol=1e-10)
    expected_objective = np.linalg.norm(a @ optimum - y) ** 2
    assert json.loads(result.to_json())["diagnostics"]["objective"] == pytest.approx(
        expected_objective, abs=1e-12, rel=1e-9
    )
    assert np.linalg.norm(a.T @ (a @ values - y), ord=np.inf) < 1e-8


def test_iteration_exhaustion_rejected_by_strict_api():
    result = prepare_json(
        FIXTURE.read_text(), algorithm=LoopClosureAlgorithm.LSE_UNREDUCED_PCG
    ).solve(max_iterations=0)
    assert result.status == "IterationLimit"
    with pytest.raises(RuntimeError):
        result.require_converged()


def test_readonly_result_buffer_outlives_report():
    result = prepare_json(FIXTURE.read_text()).solve()
    view = coordinate_buffer(result)
    del result
    assert view.readonly
    assert list(view)[3:6] == pytest.approx([10.2, 0.0, 0.0])
    with pytest.raises(TypeError):
        view[0] = 42.0


def test_standard_library_buffers_need_no_numpy_import():
    script = """
import builtins
original = builtins.__import__
def without_numpy(name, *args, **kwargs):
    if name.split('.')[0] in {'numpy', 'scipy'}:
        raise ImportError('optional scientific package blocked by qualification test')
    return original(name, *args, **kwargs)
builtins.__import__ = without_numpy
from array import array
from cave_closure import LoopClosureAlgorithm, LoopClosureMode, prepare_buffers, coordinates
endpoints = memoryview(array('Q', [0,1])).cast('B').cast('Q', shape=(1,2))
displacements = memoryview(array('d', [10.,0.,0.])).cast('B').cast('d', shape=(1,3))
covariances = memoryview(array('d', [1.,0.,0.,0.,1.,0.,0.,0.,1.])).cast('B').cast('d', shape=(1,3,3))
prepared = prepare_buffers(station_ids=['A','B'], observation_ids=['ab'], endpoints=endpoints,
                           displacements=displacements, covariances=covariances, mode=LoopClosureMode.SPATIAL)
result = prepared.solve()
result.require_converged()
assert coordinates(result)['B'] == (10.,0.,0.)
"""
    subprocess.run([sys.executable, "-c", script], check=True)
