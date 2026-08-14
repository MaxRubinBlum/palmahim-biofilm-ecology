#!/usr/bin/env python3
"""Validate reconstructed inputs against published Supplementary Table S5.

Run after `00_prepare_reproducibility_inputs.py`. The script fails if the key
beta-diversity and functional-redundancy values no longer reproduce the final
manuscript table within tight numerical tolerance.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon

EXPECTED_MEAN_BC = {
    "taxonomic": 0.6210749038914138,
    "guild": 0.13066208127467652,
    "primary": 0.069498470943168,
}
EXPECTED_WILCOXON_P = 3.814697265625e-06
EXPECTED_REDUNDANCY = {
    "Methanotroph": (56, 0.10557377850262348, 0.34186264956462364, 0.6365795298467646, 18, 5.182640584838683),
    "Sulfur_oxidizing_autotroph": (68, 0.0380177505757026, 0.09912832057477684, 0.18342028806508898, 18, 5.25019571923488),
    "Autotrophic_carbon_fixation": (133, 0.05030632875825378, 0.24354850495100763, 0.35090368365433233, 18, 6.978115142936405),
    "CBB_I": (95, 0.04260993323682298, 0.16682002091061932, 0.23692374197071858, 18, 6.519453620189055),
    "CBB_II": (27, 0.009311855034296717, 0.05605683354418451, 0.1428109408151841, 17, 2.67307830145243),
    "rTCA": (19, 0.0009895365670411023, 0.005903263292698565, 0.07223099382228745, 6, 3.575476265458539),
}


def read(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def hellinger(x):
    x = np.asarray(x, float)
    totals = x.sum(axis=1, keepdims=True)
    return np.sqrt(np.divide(x, totals, out=np.zeros_like(x), where=totals != 0))


def dist(x):
    return squareform(pdist(hellinger(x), metric="braycurtis"))


def mean_upper(d):
    return float(d[np.triu_indices_from(d, k=1)].mean())


def sample_means(d):
    return d.sum(axis=1) / (d.shape[0] - 1)


def effective_number(row):
    x = np.asarray(row, float)
    x = x[x > 0]
    if not len(x):
        return np.nan
    p = x / x.sum()
    return 1.0 / np.sum(p ** 2)


def assert_close(name, observed, expected, atol=1e-12):
    if not np.isclose(observed, expected, rtol=0, atol=atol):
        raise AssertionError(f"{name}: observed {observed!r}, expected {expected!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="data/generated")
    args = ap.parse_args()
    ddir = Path(args.datadir)

    A = read(ddir / "mag_abundance.tsv").astype(float)
    G = read(ddir / "strict_guild_membership.tsv").astype(float)
    P = read(ddir / "primary_production_membership.tsv").astype(float)
    T = read(ddir / "curated_mag_traits.tsv")
    meta = read(ddir / "sample_metadata.tsv")

    samples = [s for s in A.index if bool(meta.loc[s, "include_taxonomic_functional_beta_diversity"])]
    if len(samples) != 18:
        raise AssertionError(f"Expected 18 S5 samples, found {len(samples)}")
    A = A.loc[samples]
    mags = list(A.columns)
    G = G.reindex(mags).fillna(0)
    P = P.reindex(mags).fillna(0)
    T = T.reindex(mags)

    Dtax = dist(A.to_numpy())
    Dguild = dist(A.to_numpy() @ G.to_numpy())
    Dprimary = dist(A.to_numpy() @ P.to_numpy())
    ds = {"taxonomic": Dtax, "guild": Dguild, "primary": Dprimary}

    for name, expected in EXPECTED_MEAN_BC.items():
        observed = mean_upper(ds[name])
        assert_close(f"mean Bray-Curtis {name}", observed, expected)
        print(f"PASS mean Bray-Curtis {name}: {observed:.16g}")

    for name, other in (("guild", Dguild), ("primary", Dprimary)):
        p = float(wilcoxon(sample_means(Dtax), sample_means(other), alternative="greater").pvalue)
        assert_close(f"Wilcoxon P taxonomic>{name}", p, EXPECTED_WILCOXON_P, atol=1e-15)
        print(f"PASS Wilcoxon taxonomic>{name}: P={p:.16g}")

    for trait, expected in EXPECTED_REDUNDANCY.items():
        carriers = pd.to_numeric(T[trait], errors="coerce").fillna(0).astype(bool)
        X = A.loc[:, carriers.index[carriers]].to_numpy()
        summed = X.sum(axis=1)
        eff = np.array([effective_number(row) for row in X])
        observed = (int(carriers.sum()), float(summed.min()), float(np.median(summed)),
                    float(summed.max()), int((summed >= 0.01).sum()), float(np.nanmedian(eff)))
        if observed[0] != expected[0] or observed[4] != expected[4]:
            raise AssertionError(f"{trait}: count/core mismatch {observed} vs {expected}")
        for label, obs, exp in zip(("min", "median", "max", "effective"),
                                   (observed[1], observed[2], observed[3], observed[5]),
                                   (expected[1], expected[2], expected[3], expected[5])):
            assert_close(f"{trait} {label}", obs, exp)
        print(f"PASS redundancy {trait}: carriers={observed[0]}, samples>=1%={observed[4]}, effective={observed[5]:.12g}")

    print("\nAll Supplementary Table S5 validation checks passed.")

if __name__ == "__main__":
    main()
