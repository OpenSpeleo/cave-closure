"""Seeded independent random multigraph corpus and dense/sparse reference solvers.

Run directly with --output-dir to generate persistent or temporary JSON cases.
The production solver is never used to construct measurements or reference values.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from fixture_corpus import SOLVERS, generate_solver_results, verify_result
from scipy import linalg, sparse
from scipy.sparse.linalg import lsmr


def generate_graph(
    *,
    seed: int,
    stations: int,
    extra_edges: int,
    components: int = 1,
    topology: str = "random",
    controls_per_component: int = 2,
):
    if (
        stations < components
        or components < 1
        or extra_edges < 0
        or controls_per_component < 1
    ):
        raise ValueError("invalid graph generation sizes")
    rng = np.random.default_rng(seed)
    positions = np.cumsum(rng.normal(scale=5.0, size=(stations, 3)), axis=0)
    groups = [
        list(map(int, group))
        for group in np.array_split(np.arange(stations), components)
    ]
    edges, controls = [], []
    for group in groups:
        for i in range(1, len(group)):
            parent = (
                group[i - 1]
                if topology in {"chain", "ring"}
                else int(rng.choice(group[:i]))
            )
            edges.append((parent, group[i]))
        if topology == "ring" and len(group) > 2:
            edges.append((group[-1], group[0]))
        count = min(controls_per_component, len(group))
        chosen = [group[0]] + list(rng.choice(group[1:], count - 1, replace=False))
        controls.extend(chosen)
    for _ in range(extra_edges):
        group = groups[int(rng.integers(components))]
        u, v = map(int, rng.choice(group, 2, replace=True))
        edges.append((u, v))
    # Explicit parallel, reversed and self-loop records exercise multigraph
    # assembly; they carry distinct stochastic evidence and observation IDs.
    if edges:
        u, v = edges[int(rng.integers(len(edges)))]
        edges.extend([(u, v), (v, u), (u, u)])
    observations = []
    for i, (u, v) in enumerate(edges):
        matrix = rng.normal(size=(3, 3))
        covariance = matrix @ matrix.T + rng.uniform(0.25, 1.5) * np.eye(3)
        noise = rng.multivariate_normal(np.zeros(3), covariance)
        observations.append(
            {
                "id": f"observation:{i:08d}",
                "kind": "displacement",
                "from": f"station:{u:08d}",
                "to": f"station:{v:08d}",
                "value": (positions[v] - positions[u] + noise).tolist(),
                "uncertainty": {"kind": "covariance", "matrix": covariance.tolist()},
                "reading_ids": [f"independent-reading:{i}"],
            }
        )
    return {
        "schema_version": 1,
        "frame": {
            "kind": "local_enu",
            "length_unit": "m",
            "provenance": [
                f"Synthetic independent Gaussian multigraph; seed={seed}; topology={topology}"
            ],
        },
        "mode": {"kind": "spatial"},
        "stations": [{"id": f"station:{i:08d}"} for i in range(stations)],
        "observations": observations,
        "controls": [
            {
                "id": f"control:{i:08d}",
                "kind": "hard",
                "station": f"station:{i:08d}",
                "axes": ["x", "y", "z"],
                "value": positions[i].tolist(),
            }
            for i in sorted(controls)
        ],
        "relationships": [],
        "algorithm": {"objective": "gaussian", "backend": "reduced_pcg"},
    }


def design(problem):
    """Independent scalar rectangular design; hard XYZ controls are substituted."""
    if problem["relationships"] or problem["mode"] != {"kind": "spatial"}:
        raise ValueError(
            "reference corpus requires spatial factors without exact relationships"
        )
    if any(
        c["kind"] != "hard" or c["axes"] != ["x", "y", "z"] for c in problem["controls"]
    ):
        raise ValueError("reference corpus requires full XYZ hard controls")
    ids = [s["id"] for s in problem["stations"]]
    fixed = {c["station"]: np.array(c["value"]) for c in problem["controls"]}
    free = [identity for identity in ids if identity not in fixed]
    columns = {identity: 3 * i for i, identity in enumerate(free)}
    rows, cols, data, targets = [], [], [], []
    for i, observation in enumerate(problem["observations"]):
        u, v = observation["from"], observation["to"]
        uncertainty = observation["uncertainty"]
        if uncertainty["kind"] == "covariance":
            covariance = np.array(uncertainty["matrix"])
        elif uncertainty["kind"] == "isotropic":
            covariance = uncertainty["sigma_m"] ** 2 * np.eye(3)
        else:
            raise ValueError("unsupported reference uncertainty model")
        whitening = linalg.solve_triangular(
            linalg.cholesky(covariance, lower=True), np.eye(3), lower=True
        )
        target = np.array(observation["value"], dtype=float)
        if u in fixed:
            target += fixed[u]
        if v in fixed:
            target -= fixed[v]
        targets.extend(whitening @ target)
        for station, sign in ((u, -1.0), (v, 1.0)):
            if station in fixed:
                continue
            for a in range(3):
                for b in range(3):
                    rows.append(3 * i + a)
                    cols.append(columns[station] + b)
                    data.append(sign * whitening[a, b])
    a = sparse.coo_matrix(
        (data, (rows, cols)), shape=(len(targets), len(free) * 3)
    ).tocsr()
    return a, np.array(targets), free, fixed


def reference(problem, *, method="dense_qr"):
    a, y, free, fixed = design(problem)
    if not free:
        x = np.empty(0)
    elif method == "dense_qr":
        x, _, rank, _ = linalg.lstsq(a.toarray(), y, lapack_driver="gelsy")
        if rank != a.shape[1]:
            raise RuntimeError("reference design is rank deficient")
    elif method == "sparse_lsmr":
        result = lsmr(
            a, y, atol=1e-13, btol=1e-13, maxiter=max(1000, 4 * a.shape[1]), conlim=1e13
        )
        x = result[0]
        if result[1] not in {1, 2, 4, 5}:
            raise RuntimeError(f"independent LSMR did not converge: {result[1]}")
    else:
        raise ValueError("unknown reference method")
    residual = a @ x - y
    coordinates = {
        **{key: value.tolist() for key, value in fixed.items()},
        **{key: x[3 * i : 3 * i + 3].tolist() for i, key in enumerate(free)},
    }
    return {
        "schema_version": 1,
        "method": method,
        "coordinates": coordinates,
        "observation_residuals": {
            observation["id"]: (
                np.array(coordinates[observation["to"]])
                - np.array(coordinates[observation["from"]])
                - np.array(observation["value"])
            ).tolist()
            for observation in problem["observations"]
        },
        "objective": float(residual @ residual),
        "stationarity_infinity": float(np.max(np.abs(a.T @ residual), initial=0.0)),
        "scalar_redundancy": a.shape[0] - a.shape[1],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, nargs="+", default=[3, 17, 101])
    parser.add_argument("--stations", type=int, nargs="+", default=[1000])
    parser.add_argument(
        "--extra-edges", type=int, help="defaults to twice the station count"
    )
    parser.add_argument("--components", type=int, default=1)
    parser.add_argument(
        "--topologies",
        nargs="+",
        choices=["random", "chain", "ring"],
        default=["random", "ring"],
    )
    parser.add_argument("--reference", choices=["dense_qr", "sparse_lsmr", "auto"])
    parser.add_argument(
        "--solver-results",
        action="store_true",
        help="generate and verify every native solver result; requires --reference",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.solver_results and args.reference is None:
        parser.error("--solver-results requires --reference")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for stations in args.stations:
        for topology in args.topologies:
            for seed in args.seeds:
                problem = generate_graph(
                    seed=seed,
                    stations=stations,
                    extra_edges=args.extra_edges
                    if args.extra_edges is not None
                    else 2 * stations,
                    components=args.components,
                    topology=topology,
                )
                name = f"{topology}-{stations}-{seed}"
                (args.output_dir / f"{name}.json").write_text(
                    json.dumps(problem, allow_nan=False)
                )
                print(
                    f"{name}: {len(problem['stations'])} stations, {len(problem['observations'])} independent observations"
                )
    if args.reference:
        cases = []
        for path in sorted(args.output_dir.glob("*.json")):
            if path.name == "corpus.json" or path.name.endswith(
                (".reference.json", ".result.json")
            ):
                continue
            problem = json.loads(path.read_text())
            if (
                not isinstance(problem, dict)
                or "stations" not in problem
                or "observations" not in problem
            ):
                continue
            method = args.reference
            if method == "auto":
                # Offline reference memory policy, recorded in every reference;
                # it never selects or changes a production solver backend.
                unknowns = 3 * (len(problem["stations"]) - len(problem["controls"]))
                method = "dense_qr" if unknowns <= 600 else "sparse_lsmr"
            oracle = reference(problem, method=method)
            reference_name = f"{path.stem}.reference.json"
            (args.output_dir / reference_name).write_text(
                json.dumps(oracle, allow_nan=False)
            )
            cases.append(
                {
                    "name": path.stem,
                    "input": path.name,
                    "reference": reference_name,
                    "results": {
                        solver: f"{path.stem}.{solver}.result.json"
                        for solver in SOLVERS
                    },
                }
            )
        (args.output_dir / "corpus.json").write_text(
            json.dumps({"schema_version": 1, "cases": cases}, indent=2) + "\n"
        )
        if args.solver_results:
            generate_solver_results(args.output_dir, args.output_dir)
            for case in cases:
                problem = json.loads((args.output_dir / case["input"]).read_text())
                oracle = json.loads((args.output_dir / case["reference"]).read_text())
                for solver, filename in case["results"].items():
                    verify_result(
                        problem,
                        oracle,
                        json.loads((args.output_dir / filename).read_text()),
                        solver,
                    )
            print(
                f"Verified {len(cases)} fixtures against independent references with all {len(SOLVERS)} solvers"
            )


if __name__ == "__main__":
    main()
