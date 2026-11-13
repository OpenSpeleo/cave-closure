"""Generated surveys, not production-generated golden coordinates."""

import json

import numpy as np
import pytest
from cave_closure import LoopClosureAlgorithm, coordinates, prepare_json
from generate_graphs import design, generate_graph, reference


@pytest.mark.parametrize("seed", [3, 17, 101, 617])
@pytest.mark.parametrize("topology", ["random", "ring", "chain"])
def test_generated_graph_against_independent_pivoted_qr(seed, topology):
    problem = generate_graph(
        seed=seed, stations=45, extra_edges=70, components=3, topology=topology
    )
    expected = reference(problem)
    prepared = prepare_json(json.dumps(problem))
    reduced = prepared.solve()
    reduced.require_converged()
    ordinary = prepare_json(
        json.dumps(problem), algorithm=LoopClosureAlgorithm.LSE_UNREDUCED_PCG
    ).solve()
    ordinary.require_converged()
    for actual in [coordinates(reduced), coordinates(ordinary)]:
        for station, value in expected["coordinates"].items():
            np.testing.assert_allclose(actual[station], value, atol=1e-9, rtol=1e-10)
    diagnostics = json.loads(reduced.to_json())["diagnostics"]
    assert diagnostics["objective"] == pytest.approx(
        expected["objective"], abs=1e-12, rel=1e-9
    )
    assert diagnostics["stochastic_redundancy"] == expected["scalar_redundancy"]


@pytest.mark.parametrize("seed", [5, 41, 509])
def test_large_random_network_against_independent_sparse_lsmr(seed):
    problem = generate_graph(seed=seed, stations=1200, extra_edges=2400, components=4)
    expected = reference(problem, method="sparse_lsmr")
    result = prepare_json(json.dumps(problem)).solve()
    result.require_converged()
    actual = coordinates(result)
    for station, value in expected["coordinates"].items():
        np.testing.assert_allclose(actual[station], value, atol=1e-8, rtol=1e-10)
    a, y, free, _ = design(problem)
    q = np.array([actual[key] for key in free]).ravel()
    assert np.linalg.norm(a.T @ (a @ q - y), ord=np.inf) < 1e-7
    assert json.loads(result.to_json())["diagnostics"]["objective"] == pytest.approx(
        expected["objective"], abs=1e-10, rel=1e-9
    )


def test_random_input_permutation_and_id_relabeling():
    problem = generate_graph(seed=71, stations=60, extra_edges=90)
    first = coordinates(prepare_json(json.dumps(problem)).solve())
    rng = np.random.default_rng(19)
    rng.shuffle(problem["stations"])
    rng.shuffle(problem["observations"])
    for obs in problem["observations"][::2]:
        obs["from"], obs["to"] = obs["to"], obs["from"]
        obs["value"] = [-v for v in obs["value"]]
    second = coordinates(prepare_json(json.dumps(problem)).solve())
    for station in first:
        np.testing.assert_allclose(
            first[station], second[station], atol=1e-9, rtol=1e-10
        )
