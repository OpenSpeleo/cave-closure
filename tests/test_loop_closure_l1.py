"""Owned L1 preparation and single-call computation through the Rust bindings."""

import copy
import json
from concurrent.futures import ThreadPoolExecutor

import pytest
from cave_closure import (
    CancellationToken,
    IncompleteResultError,
    ValidationError,
    prepare_l1_loop_closure,
    solve_l1_loop_closure_json,
)


def request():
    return {
        "station_ids": ["A", "B", "C"],
        "observations": [
            ("AB", 0, 1, (5.0, 0.0, 0.0)),
            ("BC", 1, 2, (15.0, 0.0, 0.0)),
        ],
        "initial_positions": [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (24.0, 0.0, 0.0)],
        "anchor_indices": [0, 2],
        "options": {
            "mode": "spatial",
            "length_floor_m": 0.1,
            "minimum_anchor_count": 2,
            "deduplicate_edges": True,
            "max_iterations": None,
        },
    }


def test_single_call_matches_prepared_solution():
    source = request()
    before = copy.deepcopy(source)
    direct = json.loads(solve_l1_loop_closure_json(source))
    prepared = json.loads(prepare_l1_loop_closure(source).solve().to_json())
    assert direct == prepared
    points = {p["station"]: p["value"] for p in direct["coordinates"]}
    assert points == {"A": [0.0, 0.0, 0.0], "B": [5.0, 0.0, 0.0], "C": [24.0, 0.0, 0.0]}
    assert direct["diagnostics"]["objective"] == pytest.approx(4 / 15)
    assert not direct["diagnostics"]["calibrated_statistics_eligible"]
    assert source == before


def test_preparation_owns_mutable_inputs():
    source = request()
    prepared = prepare_l1_loop_closure(source)
    source["initial_positions"][2] = (1e9, 0.0, 0.0)
    source["observations"].clear()
    report = prepared.solve()
    report.require_converged()
    assert dict(report.coordinates)["C"] == [24.0, 0.0, 0.0]
    evidence = json.loads(prepared.normalized_json()[0])
    assert len(evidence["observations"]) == 2


def test_shared_preparation_solves_concurrently():
    prepared = prepare_l1_loop_closure(request())
    with ThreadPoolExecutor(max_workers=4) as pool:
        outputs = list(pool.map(lambda _: prepared.solve().to_json(), range(8)))
    assert len(set(outputs)) == 1


@pytest.mark.parametrize("direct", [False, True])
def test_cancellation_keeps_the_standard_exception_contract(direct):
    token = CancellationToken()
    token.cancel()
    with pytest.raises(IncompleteResultError, match="Cancelled"):
        if direct:
            solve_l1_loop_closure_json(request(), cancellation=token)
        else:
            prepare_l1_loop_closure(request()).solve(cancellation=token)


@pytest.mark.parametrize("field", ["station_ids", "observations", "anchor_indices"])
def test_missing_request_fields_are_rejected(field):
    source = request()
    source.pop(field)
    with pytest.raises(ValidationError):
        prepare_l1_loop_closure(source)


def test_unknown_options_are_rejected():
    source = request()
    source["options"]["typo"] = True
    with pytest.raises(ValidationError, match="unexpected"):
        solve_l1_loop_closure_json(source)


def test_nonfinite_measurement_is_rejected():
    source = request()
    source["observations"][0] = ("AB", 0, 1, (float("nan"), 0.0, 0.0))
    with pytest.raises(ValidationError, match="nonfinite"):
        solve_l1_loop_closure_json(source)


def test_horizontal_request_preserves_depths():
    source = request()
    source["options"]["mode"] = "horizontal"
    source["initial_positions"][1] = (5.0, 0.0, -3.0)
    source["initial_positions"][2] = (24.0, 0.0, -8.0)
    report = json.loads(solve_l1_loop_closure_json(source))
    points = {p["station"]: p["value"] for p in report["coordinates"]}
    assert points["B"] == [5.0, 0.0, -3.0]
    assert points["C"] == [24.0, 0.0, -8.0]
