#!/usr/bin/env python3
"""Run and validate the small synthetic reviewer demo.

This demo exercises the central taxonomic-versus-functional beta-diversity
workflow used for Supplementary Fig. 3 / Table S5 without requiring manuscript
data. It is a software smoke test, not a biological result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd


REQUIRED_OUTPUTS = [
    "distance_summary.csv",
    "paired_wilcoxon.csv",
    "sample_mean_dissimilarity.csv",
    "taxonomic_abundance.csv",
    "guild_abundance.csv",
    "primary_abundance.csv",
    "metadata_analysis_samples.csv",
    "supplementary_fig3_beta_diversity.png",
    "supplementary_fig3_beta_diversity.svg",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory (default: outputs/demo)",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    demo = root / "demo"
    outdir = Path(args.outdir) if args.outdir else root / "outputs" / "demo"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(root / "scripts" / "02_taxonomic_functional_beta_diversity.py"),
        "--abundance", str(demo / "mag_abundance.tsv"),
        "--guild-membership", str(demo / "guild_membership.tsv"),
        "--primary-membership", str(demo / "primary_membership.tsv"),
        "--metadata", str(demo / "metadata.tsv"),
        "--outdir", str(outdir),
    ]

    start = time.perf_counter()
    subprocess.run(cmd, check=True)
    elapsed = time.perf_counter() - start

    missing = [name for name in REQUIRED_OUTPUTS if not (outdir / name).exists()]
    if missing:
        raise RuntimeError(f"Demo did not create expected files: {missing}")

    expected = json.loads((demo / "expected_metrics.json").read_text())
    summary = pd.read_csv(outdir / "distance_summary.csv")
    observed_means = {
        row["representation"]: float(row["mean"])
        for _, row in summary.loc[summary["subset"] == "all"].iterrows()
    }

    for key, target in expected["mean_bray_curtis"].items():
        if key not in observed_means or not np.isclose(observed_means[key], target, rtol=1e-10, atol=1e-12):
            raise AssertionError(
                f"Unexpected {key} mean Bray-Curtis: {observed_means.get(key)} != {target}"
            )

    tests = pd.read_csv(outdir / "paired_wilcoxon.csv")
    observed_p = dict(zip(tests["comparison"], tests["P"].astype(float)))
    for key, target in expected["paired_wilcoxon_p"].items():
        if key not in observed_p or not np.isclose(observed_p[key], target, rtol=0, atol=1e-12):
            raise AssertionError(f"Unexpected Wilcoxon P for {key}: {observed_p.get(key)} != {target}")

    print("PASS: synthetic reviewer demo reproduced expected outputs")
    print(f"  taxonomic mean Bray-Curtis: {observed_means['taxonomic']:.12f}")
    print(f"  guild mean Bray-Curtis:     {observed_means['guild']:.12f}")
    print(f"  primary mean Bray-Curtis:   {observed_means['primary']:.12f}")
    print(f"  paired Wilcoxon P values:   {observed_p}")
    print(f"  runtime: {elapsed:.2f} s")
    print(f"  outputs: {outdir}")


if __name__ == "__main__":
    main()
