#!/usr/bin/env python3
"""Build Supplementary Tables 6–8 summaries from an integrated MAG profile.

Input should contain MAG, Order, sample abundance columns, transporter/fermentation
columns, MEROPS/CAZyme hit and family-count columns. Optional category columns
(e.g. MEROPS_* or CAZyme_*) are averaged per order if present.

The final Supplementary Table S6 uses abundance in percentage points, so
`--abundance-units percent` is the default. Use `fraction` for 0–1 matrices.
"""
import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--sample-columns", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--abundance-units", choices=("percent", "fraction"), default="percent")
    a = ap.parse_args()

    sep = "\t" if a.profiles.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(a.profiles, sep=sep)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    A = df.set_index("MAG")[a.sample_columns].astype(float)

    threshold_1 = 1.0 if a.abundance_units == "percent" else 0.01
    threshold_5 = 5.0 if a.abundance_units == "percent" else 0.05
    base = pd.DataFrame({
        "MAG": A.index,
        "mean_relative_abundance": A.mean(axis=1),
        "max_relative_abundance": A.max(axis=1),
        "samples_above_1pct": (A >= threshold_1).sum(axis=1),
        "samples_above_5pct": (A >= threshold_5).sum(axis=1),
        "dominant_sample": A.idxmax(axis=1),
    }).set_index("MAG")

    s6 = df.set_index("MAG").join(base, how="left").reset_index()
    s6.to_csv(out / "supp_table6_heterotroph_profiles.tsv", sep="\t", index=False)

    numeric = [c for c in df.columns
               if c.startswith(("MEROPS_", "CAZyme_"))
               or c in ["MEROPS_hits", "MEROPS_family_count", "CAZyme_hits", "CAZyme_family_count"]]
    if numeric:
        agg = df.groupby("Order")[numeric].mean(numeric_only=True)
        n = df.groupby("Order").MAG.nunique().rename("n_MAGs")
        summ = n.to_frame().join(agg)
        mer = [c for c in summ.columns if c.startswith("MEROPS")]
        caz = [c for c in summ.columns if c.startswith("CAZyme")]
        summ[["n_MAGs"] + mer].to_csv(out / "supp_table7_merops_order_summary.tsv", sep="\t")
        summ[["n_MAGs"] + caz].to_csv(out / "supp_table8_cazyme_order_summary.tsv", sep="\t")

if __name__ == "__main__":
    main()
