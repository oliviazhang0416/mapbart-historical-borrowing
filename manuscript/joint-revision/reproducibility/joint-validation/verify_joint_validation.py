#!/usr/bin/env python3
"""Recompute the published joint-validation summaries from canonical.csv.

This is a no-fit check.  It uses only the canonical per-method/per-replicate
columns required by the validation summary and does not require R.
"""
from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "09-method-exploration" / "validation" / "results"
CANONICAL = RESULTS / "canonical.csv"
SUMMARY = RESULTS / "summary.csv"
CHECKCOLS = {
    "method", "scenario", "rep", "ate_mean", "truth", "lower", "upper",
    "max_rhat", "min_bulk_ess",
}
EXPECTED_FLAGS = {
    ("JOINT-LRC-01", "compatible"): 0,
    ("JOINT-SOURCE-01", "compatible"): 1,
    ("SEPARATED-SOURCE-01", "compatible"): 1,
    ("JOINT-LRC-01", "one"): 3,
    ("JOINT-SOURCE-01", "one"): 2,
    ("SEPARATED-SOURCE-01", "one"): 0,
    ("JOINT-LRC-01", "two"): 3,
    ("JOINT-SOURCE-01", "two"): 2,
    ("SEPARATED-SOURCE-01", "two"): 7,
    ("JOINT-LRC-01", "heterogeneous"): 8,
    ("JOINT-SOURCE-01", "heterogeneous"): 4,
    ("SEPARATED-SOURCE-01", "heterogeneous"): 5,
    ("JOINT-LRC-01", "null"): 4,
    ("JOINT-SOURCE-01", "null"): 3,
    ("SEPARATED-SOURCE-01", "null"): 8,
}


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def close(a: float, b: float, tol: float = 1e-12) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def main() -> None:
    rows = read_csv(CANONICAL)
    if len(rows) != 600:
        raise SystemExit(f"expected 600 canonical rows, found {len(rows)}")
    missing = CHECKCOLS.difference(rows[0])
    if missing:
        raise SystemExit(f"canonical.csv missing check columns: {sorted(missing)}")
    groups = defaultdict(list)
    for row in rows:
        for col in CHECKCOLS - {"method", "scenario"}:
            value = float(row[col])
            if not math.isfinite(value):
                raise SystemExit(f"non-finite {col} in {row['method']}/{row['scenario']}/{row['rep']}")
        groups[(row["method"], row["scenario"])].append(row)
    if set(groups) != set(EXPECTED_FLAGS):
        raise SystemExit("method/scenario groups differ from the frozen 15-group design")

    summary = {(r["method"], r["scenario"]): r for r in read_csv(SUMMARY)}
    if set(summary) != set(groups):
        raise SystemExit("summary.csv and canonical.csv have different groups")

    print("method,scenario,n,rmse,rmse_mcse,coverage_count,flags")
    for key in sorted(groups):
        group = groups[key]
        if len(group) != 40:
            raise SystemExit(f"expected 40 rows for {key}, found {len(group)}")
        errors = [float(r["ate_mean"]) - float(r["truth"]) for r in group]
        squared = [e * e for e in errors]
        rmse = math.sqrt(sum(squared) / len(squared))
        rmse_mcse = math.sqrt(statistics.variance(squared) /
                              (4.0 * len(squared) * rmse * rmse))
        coverage = sum(float(r["lower"]) <= float(r["truth"]) <= float(r["upper"])
                       for r in group)
        flags = sum(float(r["max_rhat"]) > 1.05 or
                    float(r["min_bulk_ess"]) < 400.0 for r in group)
        published = summary[key]
        if not close(rmse, float(published["ate_rmse"]), 1e-11):
            raise SystemExit(f"RMSE mismatch for {key}: {rmse} != {published['ate_rmse']}")
        if coverage != int(published["coverage_count"]):
            raise SystemExit(f"coverage mismatch for {key}")
        if flags != EXPECTED_FLAGS[key]:
            raise SystemExit(f"flag mismatch for {key}: {flags} != {EXPECTED_FLAGS[key]}")
        print(f"{key[0]},{key[1]},{len(group)},{rmse:.12f},{rmse_mcse:.12f},{coverage},{flags}")

    print("PASS: canonical.csv recomputes RMSE, RMSE MCSE, coverage, and ATE flags")


if __name__ == "__main__":
    main()
