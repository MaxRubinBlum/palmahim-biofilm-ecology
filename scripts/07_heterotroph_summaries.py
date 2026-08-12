"""Summaries underlying Supplementary Tables 6–8."""

import pandas as pd


def abundance_summary(abundance):
    return pd.DataFrame({
        "mean_relative_abundance": abundance.mean(axis=1),
        "max_relative_abundance": abundance.max(axis=1),
        "samples_above_1pct": (abundance >= 1).sum(axis=1),
        "samples_above_5pct": (abundance >= 5).sum(axis=1),
        "dominant_sample": abundance.idxmax(axis=1),
    })


def merops_order_summary(mag_profiles):
    return mag_profiles.groupby("Order").agg(
        n_MAGs=("MAG", "nunique"),
        mean_MEROPS_genes=("MEROPS_hits", "mean"),
        mean_MEROPS_families=("MEROPS_family_count", "mean"),
    )


def cazyme_order_summary(mag_profiles):
    return mag_profiles.groupby("Order").agg(
        n_MAGs=("MAG", "nunique"),
        mean_CAZyme_genes=("CAZyme_hits", "mean"),
        mean_CAZyme_families=("CAZyme_family_count", "mean"),
    )
