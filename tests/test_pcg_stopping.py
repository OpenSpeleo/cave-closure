"""Cancellation in true residuals must not make a solved network diverge."""

import json
import math

import numpy as np
import pytest
from cave_closure import (
    LoopClosureAlgorithm,
    LoopClosureMode,
    coordinates,
    prepare_arrays,
)


@pytest.mark.parametrize("count", [100, 250])
@pytest.mark.parametrize(
    "algorithm",
    [LoopClosureAlgorithm.LSE_REDUCED_PCG, LoopClosureAlgorithm.LSE_UNREDUCED_PCG],
)
def test_long_weighted_network_stops_at_scaled_true_residual(count, algorithm):
    positions = np.array([[i * 13.7, math.sin(i / 7) * 40, 0] for i in range(count)])
    edges = [(i, i + 1) for i in range(count - 1)] + [
        (i, i + 19) for i in range(0, count - 19, 7)
    ]
    deltas = np.array(
        [
            positions[v] - positions[u] + [math.sin(i) * 0.1, math.cos(i) * 0.2, 0]
            for i, (u, v) in enumerate(edges)
        ]
    )
    lengths = np.linalg.norm(deltas, axis=1)
    # Independently assemble the whitened rectangular design with station 0 fixed.
    design = np.zeros((len(edges), count - 1))
    for i, (u, v) in enumerate(edges):
        if u:
            design[i, u - 1] = -1 / math.sqrt(lengths[i])
        design[i, v - 1] = 1 / math.sqrt(lengths[i])
    rhs = deltas[:, :2] / np.sqrt(lengths[:, None])
    expected = np.linalg.lstsq(design, rhs, rcond=None)[0]
    report = prepare_arrays(
        [str(i) for i in range(count)],
        [str(i) for i in range(len(edges))],
        edges,
        deltas.tolist(),
        [(np.eye(3) * length).tolist() for length in lengths],
        mode=LoopClosureMode.HORIZONTAL,
        preserved_z=[0] * count,
        algorithm=algorithm,
        controls=[("start", 0, [0, 1], [0.0, 0.0], None)],
    ).solve()
    report.require_converged()
    actual = coordinates(report)
    np.testing.assert_allclose(
        [actual[str(i)][:2] for i in range(1, count)], expected, atol=1e-7, rtol=1e-10
    )
    diagnostics = json.loads(report.to_json())["diagnostics"]
    assert diagnostics["objective"] == pytest.approx(
        np.linalg.norm(design @ expected - rhs) ** 2, rel=1e-8, abs=1e-12
    )
    assert diagnostics["backward_error"] < 1e-10
