#!/usr/bin/env python3
"""Taxonomic versus functional beta-diversity (Supplementary Fig. 3 / Table S5).

Inputs
------
--abundance        samples x MAG relative abundance (fractions or percentages)
--guild-membership MAG x strict ecological guilds (0/1)
--primary-membership MAG x primary-production traits (0/1)
--metadata         columns: sample, habitat; optional inclusion flag

If metadata contains `include_taxonomic_functional_beta_diversity`, rows with 0
are excluded. In the manuscript this removes AnemPM22 and retains 18 biofilms.
The taxonomic matrix always uses all MAGs in the abundance table. Missing MAGs
in functional membership tables are treated as zero rather than being dropped.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt


def read_table(path: str, index_col=0):
    sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep, index_col=index_col)


def hellinger(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    totals = x.sum(axis=1, keepdims=True)
    rel = np.divide(x, totals, out=np.zeros_like(x), where=totals != 0)
    return np.sqrt(rel)


def bray_curtis(x: np.ndarray) -> np.ndarray:
    return squareform(pdist(x, metric="braycurtis"))


def mean_distance_to_others(d: np.ndarray) -> np.ndarray:
    return d.sum(axis=1) / (d.shape[0] - 1)


def distance_long(d, samples, metadata):
    rows = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            s1, s2 = samples[i], samples[j]
            h1, h2 = metadata.loc[s1, "habitat"], metadata.loc[s2, "habitat"]
            rows.append({"sample1": s1, "sample2": s2, "distance": d[i, j],
                         "habitat1": h1, "habitat2": h2,
                         "comparison": "within" if h1 == h2 else "between"})
    return pd.DataFrame(rows)


def summarize_distances(long_df, representation):
    def one(x, label):
        q1, med, q3 = np.quantile(x, [0.25, 0.50, 0.75])
        return {"representation": representation, "subset": label,
                "n_pairs": len(x), "mean": np.mean(x), "median": med,
                "q1": q1, "q3": q3}
    out = [one(long_df.distance.to_numpy(), "all")]
    for label in ("within", "between"):
        x = long_df.loc[long_df.comparison == label, "distance"].to_numpy()
        if len(x): out.append(one(x, label))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--guild-membership", required=True)
    ap.add_argument("--primary-membership", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--include-column", default="include_taxonomic_functional_beta_diversity")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    A = read_table(args.abundance).astype(float)
    G = read_table(args.guild_membership).astype(float)
    P = read_table(args.primary_membership).astype(float)
    sep = "\t" if args.metadata.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(args.metadata, sep=sep).set_index("sample")

    samples = [s for s in A.index if s in meta.index]
    if args.include_column in meta.columns:
        samples = [s for s in samples if bool(meta.loc[s, args.include_column])]
    A = A.loc[samples]
    meta = meta.loc[samples]

    # Taxonomic distances use every MAG. Functional tables are aligned to the
    # same MAG set, with absent assignments explicitly filled as zero.
    mags = list(A.columns)
    G = G.reindex(mags).fillna(0.0)
    P = P.reindex(mags).fillna(0.0)

    taxonomic = A.to_numpy()
    guild = taxonomic @ G.to_numpy()
    primary = taxonomic @ P.to_numpy()
    matrices = {"taxonomic": taxonomic, "guild": guild, "primary": primary}

    # Export the exact matrices consumed by the companion vegan script.
    pd.DataFrame(taxonomic, index=samples, columns=mags).to_csv(outdir / "taxonomic_abundance.csv")
    pd.DataFrame(guild, index=samples, columns=G.columns).to_csv(outdir / "guild_abundance.csv")
    pd.DataFrame(primary, index=samples, columns=P.columns).to_csv(outdir / "primary_abundance.csv")
    meta[["habitat"]].to_csv(outdir / "metadata_analysis_samples.csv")

    distances = {k: bray_curtis(hellinger(v)) for k, v in matrices.items()}
    means = pd.DataFrame(index=samples)
    summaries = []
    for name, d in distances.items():
        pd.DataFrame(d, index=samples, columns=samples).to_csv(outdir / f"distance_{name}.csv")
        means[name] = mean_distance_to_others(d)
        long_df = distance_long(d, samples, meta)
        long_df.to_csv(outdir / f"pairwise_{name}.csv", index=False)
        summaries.append(summarize_distances(long_df, name))
    means.index.name = "sample"
    means.join(meta[["habitat"]]).to_csv(outdir / "sample_mean_dissimilarity.csv")
    pd.concat(summaries, ignore_index=True).to_csv(outdir / "distance_summary.csv", index=False)

    tests = []
    for other in ("guild", "primary"):
        stat = wilcoxon(means.taxonomic, means[other], alternative="greater")
        tests.append({"comparison": f"taxonomic > {other}", "W": stat.statistic,
                      "P": stat.pvalue, "n_samples": len(means)})
    pd.DataFrame(tests).to_csv(outdir / "paired_wilcoxon.csv", index=False)

    order = ["taxonomic", "guild", "primary"]
    vals = [means[x].to_numpy() for x in order]
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.boxplot(vals, tick_labels=order, showfliers=False, showmeans=True)
    rng = np.random.default_rng(1)
    for i, x in enumerate(vals, 1):
        ax.scatter(np.full(len(x), i) + rng.normal(0, 0.035, len(x)), x, s=22, alpha=0.75)
    ax.set_ylabel("Mean Bray–Curtis dissimilarity to other biofilms")
    ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(outdir / "supplementary_fig3_beta_diversity.png", dpi=300)
    fig.savefig(outdir / "supplementary_fig3_beta_diversity.svg")

if __name__ == "__main__":
    main()
