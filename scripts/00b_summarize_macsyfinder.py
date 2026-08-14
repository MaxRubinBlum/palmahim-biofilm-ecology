#!/usr/bin/env python3
"""Summarize the supplied combined MacSyFinder `all_systems` output.

The archived origin file may contain repeated header rows created while
concatenating per-MAG outputs. This script removes them explicitly and writes:
1. one row per MAG/model/system with completeness and component-gene evidence;
2. a compact MAG x model count matrix suitable for provenance inspection.

This is a provenance helper. Final ecological trait calls remain the curated
assignments in Supplementary Table S4, applying the operational definitions in S1.
"""
import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="combined MacSyFinder ALL_systems TSV")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(a.input, sep="\t", dtype=str, low_memory=False)

    # Remove concatenated header rows and malformed rows conservatively.
    d = d[d["MAG"].notna() & d["sys_id"].notna()].copy()
    d = d[(d["MAG"] != "MAG") & (d["sys_id"] != "sys_id")]
    d["sys_wholeness_num"] = pd.to_numeric(d["sys_wholeness"], errors="coerce")

    def genes(x):
        return ";".join(sorted(set(x.dropna().astype(str))))

    group_cols = ["MAG", "Tool", "model_fqn", "sys_id"]
    system = d.groupby(group_cols, dropna=False).agg(
        sys_wholeness=("sys_wholeness_num", "max"),
        n_hits=("hit_id", "count"),
        n_unique_genes=("gene_name", "nunique"),
        component_genes=("gene_name", genes),
    ).reset_index()

    mandatory = (d[d["hit_status"].eq("mandatory")]
                 .groupby(group_cols)["gene_name"].apply(genes).rename("mandatory_genes"))
    accessory = (d[d["hit_status"].eq("accessory")]
                 .groupby(group_cols)["gene_name"].apply(genes).rename("accessory_genes"))
    system = system.join(mandatory, on=group_cols).join(accessory, on=group_cols)
    system.to_csv(out / "macsyfinder_system_summary.tsv", sep="\t", index=False)

    counts = system.groupby(["MAG", "model_fqn"]).size().unstack(fill_value=0)
    counts.to_csv(out / "macsyfinder_MAG_system_matrix.tsv", sep="\t")
    print(f"Retained {len(d)} hit rows, {len(system)} systems, {len(counts)} MAGs")

if __name__ == "__main__":
    main()
