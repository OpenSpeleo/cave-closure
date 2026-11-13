"""Measure generated cases at declared accuracy, without performance pass/fail gates."""

import argparse
import json
import statistics
import time
from pathlib import Path

from cave_closure import coordinates, prepare_json


def benchmark(path, repeats):
    text = path.read_text()
    prepare_times, solve_times = [], []
    for _ in range(repeats):
        started = time.perf_counter()
        prepared = prepare_json(text)
        prepared_at = time.perf_counter()
        result = prepared.solve()
        solved_at = time.perf_counter()
        result.require_converged()
        prepare_times.append(1000 * (prepared_at - started))
        solve_times.append(1000 * (solved_at - prepared_at))
    report = json.loads(result.to_json())
    reference_path = path.with_suffix(".reference.json")
    record = {
        "case": path.name,
        "repeats": repeats,
        "preparation_median_ms": statistics.median(prepare_times),
        "solve_median_ms": statistics.median(solve_times),
        **{
            key: report["diagnostics"][key]
            for key in (
                "objective",
                "iterations",
                "kernel_unknowns",
                "original_unknowns",
                "backward_error",
            )
        },
    }
    if reference_path.exists():
        reference = json.loads(reference_path.read_text())
        actual = coordinates(result)
        record["objective_gap"] = (
            report["diagnostics"]["objective"] - reference["objective"]
        )
        record["maximum_coordinate_difference_m"] = max(
            abs(actual[key][axis] - value[axis])
            for key, value in reference["coordinates"].items()
            for axis in range(3)
        )
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for path in sorted(args.directory.glob("*.json")):
        if path.name.endswith(".reference.json"):
            continue
        record = benchmark(path, args.repeats)
        records.append(record)
        print(json.dumps(record, allow_nan=False))
    args.output.write_text(json.dumps(records, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
