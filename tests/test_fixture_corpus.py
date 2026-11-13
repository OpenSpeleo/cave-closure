"""Every persisted solver result and fresh solve must match an independent optimum."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
from cave_closure import LoopClosureAlgorithm, prepare_json
from fixture_corpus import (
    ALGORITHMS,
    coordinate_tolerance,
    generate_solver_results,
    verify_result,
)
from generate_graphs import design, reference

FIXTURES = Path(__file__).parent / "fixtures"


def cases():
    corpus = json.loads((FIXTURES / "corpus.json").read_text())
    assert corpus["schema_version"] == 1
    assert corpus["cases"], "the persisted fixture corpus must not be empty"
    return corpus["cases"]


CASES = cases()


@pytest.fixture(scope="module", params=CASES, ids=lambda case: case["name"])
def fixture_case(request):
    case = request.param
    problem = json.loads((FIXTURES / case["input"]).read_text())
    expected = json.loads((FIXTURES / case["reference"]).read_text())
    return case, problem, expected, design(problem)


@pytest.fixture(scope="module")
def fresh_native_results(tmp_path_factory):
    # Build once per suite, then run all four backends once per input. This also
    # exercises optional direct solvers when the installed Python wheel omits them.
    output = tmp_path_factory.mktemp("cave-closure-native-corpus")
    generate_solver_results(FIXTURES, output)
    return output


def test_corpus_is_complete():
    assert len({case["name"] for case in CASES}) == len(CASES)
    inputs, artifacts = set(), set()
    for case in CASES:
        assert set(case["results"]) == set(ALGORITHMS)
        inputs.add(case["input"])
        for name in (case["input"], case["reference"], *case["results"].values()):
            assert Path(name).name == name
            assert name not in artifacts, f"corpus cases share an artifact: {name}"
            artifacts.add(name)
            assert (FIXTURES / name).is_file(), f"missing corpus artifact: {name}"
    actual_inputs = set()
    for path in FIXTURES.glob("*.json"):
        if path.name == "corpus.json" or path.name.endswith(
            (".reference.json", ".result.json")
        ):
            continue
        data = json.loads(path.read_text())
        if "stations" in data and "observations" in data:
            actual_inputs.add(path.name)
    assert inputs == actual_inputs


def test_reference_satisfies_independent_optimality_checks(fixture_case):
    _, problem, expected, system = fixture_case
    assert expected["schema_version"] == 1
    a, y, free, fixed = system
    values = np.array([expected["coordinates"][station] for station in free]).ravel()
    residual = a @ values - y
    assert np.max(np.abs(a.T @ residual), initial=0.0) <= 1e-7
    assert float(residual @ residual) == pytest.approx(
        expected["objective"], abs=1e-12, rel=1e-9
    )
    assert expected["scalar_redundancy"] == a.shape[0] - a.shape[1]
    for station, value in fixed.items():
        np.testing.assert_array_equal(expected["coordinates"][station], value)
    if a.shape[1] <= 600:
        # Cross-check small persisted references through a fresh rectangular QR.
        qr = reference(problem, method="dense_qr")
        for station in expected["coordinates"]:
            np.testing.assert_allclose(
                expected["coordinates"][station],
                qr["coordinates"][station],
                atol=coordinate_tolerance(qr),
                rtol=0,
            )


@pytest.mark.parametrize("solver", ALGORITHMS)
def test_saved_solver_result_matches_ground_truth(fixture_case, solver):
    case, problem, expected, system = fixture_case
    result = json.loads((FIXTURES / case["results"][solver]).read_text())
    verify_result(problem, expected, result, solver, system=system)


@pytest.mark.parametrize("solver", ALGORITHMS)
def test_fresh_native_solver_matches_ground_truth(
    fixture_case, solver, fresh_native_results
):
    case, problem, expected, system = fixture_case
    result = json.loads((fresh_native_results / case["results"][solver]).read_text())
    verify_result(problem, expected, result, solver, system=system)
    saved = json.loads((FIXTURES / case["results"][solver]).read_text())
    actual = result["report"]["coordinates"]
    previous = saved["report"]["coordinates"]
    assert [value["station"] for value in actual] == [
        value["station"] for value in previous
    ]
    np.testing.assert_allclose(
        [value["value"] for value in actual],
        [value["value"] for value in previous],
        atol=coordinate_tolerance(expected),
        rtol=0,
    )


@pytest.mark.parametrize(
    "reduce", [False, True], ids=["LSE_UNREDUCED_PCG", "LSE_REDUCED_PCG"]
)
def test_installed_python_solver_matches_ground_truth(fixture_case, reduce):
    _, problem, expected, system = fixture_case
    report = prepare_json(
        json.dumps(problem),
        algorithm=(
            LoopClosureAlgorithm.LSE_REDUCED_PCG
            if reduce
            else LoopClosureAlgorithm.LSE_UNREDUCED_PCG
        ),
    ).solve()
    report.require_converged()
    solver = "LSE_REDUCED_PCG" if reduce else "LSE_UNREDUCED_PCG"
    result = {
        "schema_version": 1,
        "algorithm": solver,
        "report": json.loads(report.to_json()),
    }
    verify_result(problem, expected, result, solver, system=system)


@pytest.mark.parametrize("field", ["coordinates", "residual", "objective"])
def test_comparison_rejects_incorrect_results(field):
    problem = json.loads((FIXTURES / "controlled_tree.json").read_text())
    expected = json.loads((FIXTURES / "controlled_tree.reference.json").read_text())
    result = copy.deepcopy(
        json.loads(
            (FIXTURES / "controlled_tree.LSE_REDUCED_PCG.result.json").read_text()
        )
    )
    if field == "coordinates":
        result["report"]["coordinates"][1]["value"][0] += 0.01
    elif field == "residual":
        result["report"]["observations"][0]["residual"][0] += 0.01
    else:
        result["report"]["diagnostics"]["objective"] += 1.0
    with pytest.raises(AssertionError):
        verify_result(problem, expected, result, "LSE_REDUCED_PCG")
