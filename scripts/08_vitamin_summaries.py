#!/usr/bin/env python3
"""Abundance-weighted vitamin-provider summaries (Supplementary Figs. 7–8).

Input abundance is samples x MAGs; metadata contains sample, habitat; traits are
MAG rows with complete vitamin-pathway 0/1 columns, Order, optional hierarchical
functional category, and optional B12_steps (0–7).
"""
import argparse
from pathlib import Path
import pandas as pd


def b12_status(n):
    return "complete" if n == 7 else "near-complete" if n == 6 else "partial" if n >= 1 else "absent"


def read(p, index_col=0):
    return pd.read_csv(p, sep="\t" if p.endswith((".tsv", ".txt")) else ",", index_col=index_col)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--vitamin-columns", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--display-threshold-pct", type=float, default=5.0)
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    A = read(a.abundance).astype(float)
    T = read(a.traits)
    meta = read(a.metadata)
    samples = [s for s in A.index if s in meta.index]
    mags = [m for m in A.columns if m in T.index]
    A = A.loc[samples, mags]; T = T.loc[mags]; meta = meta.loc[samples]

    # Habitat-specific mean MAG abundance, matching the Supplementary Methods.
    H = A.groupby(meta.habitat).mean()

    rows = []
    for vit in a.vitamin_columns:
        carriers = T[vit].fillna(0).astype(bool)
        for habitat in H.index:
            vals = H.loc[habitat, carriers.index[carriers]]
            total = float(vals.sum())
            orders = T.loc[vals.index].assign(abundance=vals.values).groupby("Order").abundance.sum()
            for order, x in orders.items():
                rows.append({"vitamin": vit, "habitat": habitat, "Order": order,
                             "provider_abundance": x,
                             "provider_fraction_pct": 100 * x / total if total else 0})
    tax = pd.DataFrame(rows)
    keep = tax.groupby("Order").provider_fraction_pct.max()
    keep = set(keep[keep >= a.display_threshold_pct].index)
    tax["display_order"] = tax.Order.where(tax.Order.isin(keep), "Other")
    tax.to_csv(out / "vitamin_providers_by_order.tsv", sep="\t", index=False)

    if "B12_steps" in T.columns:
        b = T[["B12_steps"]].copy()
        b["B12_status"] = b.B12_steps.fillna(0).astype(int).map(b12_status)
        b.to_csv(out / "b12_pathway_completeness.tsv", sep="\t")

        # Supplementary Fig. 8c: percentage of total habitat abundance assigned
        # to complete, near-complete, partial, or absent M00122 profiles.
        comp = []
        for habitat in H.index:
            denom = float(H.loc[habitat].sum())
            for status in ("complete", "near-complete", "partial", "absent"):
                smags = b.index[b.B12_status.eq(status)]
                x = float(H.loc[habitat, smags].sum()) if len(smags) else 0.0
                comp.append({"habitat": habitat, "B12_status": status,
                             "abundance": x,
                             "total_habitat_abundance": denom,
                             "habitat_fraction_pct": 100 * x / denom if denom else 0})
        pd.DataFrame(comp).to_csv(out / "b12_completeness_by_habitat.tsv", sep="\t", index=False)

        if "functional_category" in T.columns:
            complete = b.index[b.B12_status.eq("complete")]
            fr = []
            for habitat in H.index:
                tmp = T.loc[complete, ["functional_category"]].assign(abundance=H.loc[habitat, complete].values)
                z = tmp.groupby("functional_category").abundance.sum()
                total = float(z.sum())
                for cat, x in z.items():
                    fr.append({"habitat": habitat, "functional_category": cat,
                               "abundance": x,
                               "complete_provider_fraction_pct": 100 * x / total if total else 0})
            pd.DataFrame(fr).to_csv(out / "b12_complete_producers_by_function.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
