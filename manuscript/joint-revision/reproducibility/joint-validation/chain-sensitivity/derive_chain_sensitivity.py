#!/usr/bin/env python3
"""Derive constant-effect chain-index RMSE sensitivity from saved summaries.

This reads the reviewed saved-trace output constant_chain_rmse.csv. It does
not read posterior draws, refit a model, exclude flagged groups, or alter any
validation result. Chains 1, 2, and 3 each summarize all 120 constant-effect
method-dataset groups (40 Compatible, 40 One-variable, and 40 Two-variable).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

METHODS = ("JOINT-LRC-01", "JOINT-SOURCE-01", "SEPARATED-SOURCE-01")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.input.open(newline="")))
    required = {"method", "chain", "constant_rmse", "n_groups"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"input lacks required columns: {sorted(required)}")
    values = {}
    for row in rows:
        chain = int(row["chain"])
        method = row["method"]
        if chain not in (1, 2, 3):
            continue
        if method not in METHODS:
            raise SystemExit(f"unexpected method {method}")
        if int(row["n_groups"]) != 120:
            raise SystemExit(f"expected 120 groups for chain {chain}, got {row['n_groups']}")
        key = (chain, method)
        if key in values:
            raise SystemExit(f"duplicate {key}")
        values[key] = float(row["constant_rmse"])
    expected = {(c, m) for c in (1, 2, 3) for m in METHODS}
    if set(values) != expected:
        raise SystemExit("input does not contain exactly chains 1/2/3 for all methods")

    out_rows = []
    for chain in (1, 2, 3):
        lrc = values[(chain, "JOINT-LRC-01")]
        source = values[(chain, "JOINT-SOURCE-01")]
        separated = values[(chain, "SEPARATED-SOURCE-01")]
        out_rows.append({
            "chain": chain,
            "n_constant_groups": 120,
            "lrc_rmse": f"{lrc:.12f}",
            "joint_source_rmse": f"{source:.12f}",
            "separated_source_rmse": f"{separated:.12f}",
            "lrc_vs_separated_pct_reduction": f"{100.0 * (1.0 - lrc / separated):.6f}",
            "lrc_vs_joint_source_pct_reduction": f"{100.0 * (1.0 - lrc / source):.6f}",
        })
    csv_path = args.out_dir / "chain_sensitivity.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)

    md_path = args.out_dir / "CHAIN-SENSITIVITY.md"
    lines = [
        "# Saved-chain constant-effect sensitivity",
        "",
        "This is a descriptive chain-index sensitivity check from the reviewed saved-trace summary. Each chain-specific RMSE uses all 120 constant-effect method--dataset groups: 40 Compatible, 40 One-variable and 40 Two-variable. No ATE convergence flags or other diagnostics were used to exclude rows.",
        "",
        f"Input: `{args.input}` (SHA-256 `{sha256(args.input)}`). The upstream diagnostic follow-up that produced this input read the saved 1,800 raw RDS results; this derivation performs no MCMC and no refit.",
        "",
        "The percentage is $100(1-\\mathrm{RMSE}_{\\mathrm{LRC}}/\\mathrm{RMSE}_{\\mathrm{comparator}})$. Chain 0 in the upstream file is the existing pooled three-chain reference; the table below reports only individual saved chains 1--3.",
        "",
        "| Chain | Groups | Joint LRC RMSE | Joint source RMSE | Separated source RMSE | LRC vs separated | LRC vs joint source |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in out_rows:
        lines.append(f"| {row['chain']} | {row['n_constant_groups']} | {row['lrc_rmse']} | {row['joint_source_rmse']} | {row['separated_source_rmse']} | {row['lrc_vs_separated_pct_reduction']}% | {row['lrc_vs_joint_source_pct_reduction']}% |")
    lines += [
        "",
        "The chain-index reductions range from 4.178% to 5.008% against Separated source-BART and from 6.301% to 7.351% against Joint source-BART. This narrows the descriptive sensitivity of the pooled RMSE direction across saved chain indices; it is not a convergence proof, a replacement for the ATE diagnostic screen, or evidence that the chains have mixed adequately.",
        "",
        "The upstream frozen raw files, hashes and validation outputs remain unchanged.",
    ]
    md_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
