#!/usr/bin/env python3
"""Validate repository machine-readable inputs against manuscript S5 values."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon


def hellinger(x):
    x = np.asarray(x, dtype=float)
    totals = x.sum(axis=1, keepdims=True)
    rel = np.divide(x, totals, out=np.zeros_like(x), where=totals != 0)
    return np.sqrt(rel)


def bray_curtis(x):
    return squareform(pdist(hellinger(x), metric="braycurtis"))


def mean_pairwise(d):
    return d[np.triu_indices_from(d, k=1)].mean()


def sample_mean_distance(d):
    return d.sum(axis=1) / (d.shape[0] - 1)


def main(root):
    root = Path(root)
    abundance = pd.read_csv(root / "data/mag_abundance.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv(root / "data/sample_metadata.tsv", sep="\t").set_index("sample")
    guild = pd.read_csv(root / "data/strict_guild_membership.tsv", sep="\t", index_col=0)
    primary = pd.read_csv(root / "data/primary_production_membership.tsv", sep="\t", index_col=0)
    traits = pd.read_csv(root / "data/curated_mag_traits.tsv", sep="\t", index_col=0)

    keep = metadata.index[metadata["include_taxonomic_functional_beta_diversity"].astype(int) == 1]
    a = abundance.loc[keep]
    common = a.columns.intersection(guild.index).intersection(primary.index)
    if len(common) != a.shape[1]:
        missing = sorted(set(a.columns) - set(common))
        raise ValueError(f"MAGs missing from membership tables: {missing[:10]}")

    d_tax = bray_curtis(a.values)
    d_guild = bray_curtis(a[common].values @ guild.loc[common].values)
    d_primary = bray_curtis(a[common].values @ primary.loc[common].values)

    observed = {
        "taxonomic": mean_pairwise(d_tax),
        "guild": mean_pairwise(d_guild),
        "primary": mean_pairwise(d_primary),
    }
    expected = {
        "taxonomic": 0.6210749038914138,
        "guild": 0.13066208127467652,
        "primary": 0.069498470943168,
    }
    for key in expected:
        if not np.isclose(observed[key], expected[key], rtol=0, atol=1e-12):
            raise AssertionError(f"{key}: observed {observed[key]} != expected {expected[key]}")

    p_guild = wilcoxon(sample_mean_distance(d_tax), sample_mean_distance(d_guild), alternative="greater").pvalue
    p_primary = wilcoxon(sample_mean_distance(d_tax), sample_mean_distance(d_primary), alternative="greater").pvalue
    expected_p = 3.814697265625e-06
    if not np.isclose(p_guild, expected_p, rtol=0, atol=1e-15):
        raise AssertionError(f"guild Wilcoxon P mismatch: {p_guild}")
    if not np.isclose(p_primary, expected_p, rtol=0, atol=1e-15):
        raise AssertionError(f"primary Wilcoxon P mismatch: {p_primary}")

    redundancy_targets = {
        "Methanotroph": (56, 5.182640584838682),
        "Sulfur_oxidizing_autotroph": (68, 5.2501957192348785),
        "Autotrophic_carbon_fixation": (133, 6.978115142936405),
        "CBB_I": (95, 6.519453620189055),
        "CBB_II": (27, 2.67307830145243),
        "rTCA": (19, 3.575476265458539),
    }
    for trait, (n_expected, eff_expected) in redundancy_targets.items():
        carriers = traits.index[traits[trait].fillna(0).astype(int) == 1]
        if len(carriers) != n_expected:
            raise AssertionError(f"{trait}: carrier count {len(carriers)} != {n_expected}")
        eff = []
        for row in a[carriers].values:
            z = row[row > 0]
            if len(z):
                p = z / z.sum()
                eff.append(1.0 / np.sum(p * p))
        med = float(np.median(eff))
        if not np.isclose(med, eff_expected, rtol=0, atol=1e-12):
            raise AssertionError(f"{trait}: median effective MAG number {med} != {eff_expected}")

    print("PASS: repository inputs reproduce Supplementary Table S5 targets.")
    for key, value in observed.items():
        print(f"{key}\tmean_Bray_Curtis={value:.16g}")
    print(f"guild\tWilcoxon_P={p_guild:.16g}")
    print(f"primary\tWilcoxon_P={p_primary:.16g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    args = parser.parse_args()
    main(args.root)
