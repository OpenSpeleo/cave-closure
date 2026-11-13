"""Shared corpus orchestration and independent numerical acceptance checks.

Only the explicit generator writes persistent fixtures. Tests compute fresh
native results in temporary directories and never rewrite the saved expectations.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np

SOLVERS = ("reduced_pcg", "unreduced_pcg", "sparse_qr", "sparse_cholesky")


def generate_solver_results(fixtures: Path, output: Path):
    root = Path(__file__).resolve().parent.parent
    completed = subprocess.run(
        [
            "cargo",
            "build",
            "--release",
            "--locked",
            "--features",
            "sparse-direct",
            "--example",
            "fixture_results",
            "--message-format=json",
        ],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    executables = [
        record["executable"]
        for line in completed.stdout.splitlines()
        if (record := json.loads(line)).get("reason") == "compiler-artifact"
        and record.get("target", {}).get("name") == "fixture_results"
        and record.get("executable")
    ]
    if len(executables) != 1:
        raise RuntimeError("could not identify the native fixture result generator")
    subprocess.run(
        [executables[0], str(fixtures.resolve()), str(output.resolve())],
        cwd=root,
        check=True,
    )


def coordinate_tolerance(reference):
    scale = max(
        abs(value) for point in reference["coordinates"].values() for value in point
    )
    return 1e-9 + 1e-10 * scale


def verify_result(problem, reference, result, solver, *, system=None):
    """Check coordinates, original residuals, objective, rank and stationarity.

    The reference is the noisy problem's optimum, not its latent sampled geometry.
    Tolerances are numerical acceptance bounds, never golden-file bit equality.
    """
    # generate_graphs imports this module while defining the reference generator.
    from generate_graphs import design  # noqa: PLC0415 - circular import

    assert result["schema_version"] == 1
    assert result["solver"] == solver
    assert result["reduce"] == (solver == "reduced_pcg")
    report = result["report"]
    assert report["schema_version"] == 1
    assert report["status"] == "converged"
    actual = {entry["station"]: entry["value"] for entry in report["coordinates"]}
    assert len(actual) == len(report["coordinates"]) == len(problem["stations"])
    assert actual.keys() == reference["coordinates"].keys()
    tolerance = coordinate_tolerance(reference)
    station_ids = list(reference["coordinates"])
    actual_values = np.asarray([actual[station] for station in station_ids])
    expected_values = np.asarray(
        [reference["coordinates"][station] for station in station_ids]
    )
    assert np.all(np.isfinite(actual_values)) and np.all(np.isfinite(expected_values))
    np.testing.assert_allclose(
        actual_values,
        expected_values,
        atol=tolerance,
        rtol=0,
        err_msg=f"{solver}: station coordinates",
    )
    np.testing.assert_allclose(
        [actual[control["station"]] for control in problem["controls"]],
        [control["value"] for control in problem["controls"]],
        atol=tolerance,
        rtol=0,
    )
    residuals = {entry["id"]: entry for entry in report["observations"]}
    assert len(residuals) == len(report["observations"]) == len(problem["observations"])
    assert residuals.keys() == reference["observation_residuals"].keys()
    observation_ids = [observation["id"] for observation in problem["observations"]]
    measured_residuals = np.asarray(
        [
            np.subtract(actual[observation["to"]], actual[observation["from"]])
            - observation["value"]
            for observation in problem["observations"]
        ]
    )
    for field in ("residual", "correction"):
        values = np.asarray([residuals[key][field] for key in observation_ids])
        assert np.all(np.isfinite(values))
        np.testing.assert_allclose(
            values,
            measured_residuals,
            atol=2 * tolerance,
            rtol=0,
            err_msg=f"{solver}: original-factor {field}",
        )
    np.testing.assert_allclose(
        [residuals[key]["residual"] for key in observation_ids],
        [reference["observation_residuals"][key] for key in observation_ids],
        atol=2 * tolerance,
        rtol=0,
    )
    a, y, free, _ = system if system is not None else design(problem)
    values = np.array([actual[key] for key in free]).ravel()
    residual = a @ values - y
    objective = float(residual @ residual)
    diagnostics = report["diagnostics"]
    objective_tolerance = 1e-12 + 1e-9 * abs(reference["objective"])
    assert abs(objective - reference["objective"]) <= objective_tolerance
    assert abs(diagnostics["objective"] - objective) <= objective_tolerance
    assert (
        diagnostics["stochastic_redundancy"]
        == reference["scalar_redundancy"]
        == a.shape[0] - a.shape[1]
    )
    gradient = float(np.max(np.abs(a.T @ residual), initial=0.0))
    assert np.isfinite(gradient) and gradient <= 1e-7, (
        f"{solver}: nonstationary result {gradient}"
    )
    assert diagnostics["artificial_gauges"] == []
    assert diagnostics["warnings"] == []
