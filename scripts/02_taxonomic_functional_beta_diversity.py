"""Taxonomic versus functional beta-diversity.

Data required
-------------
A                  : samples x MAG relative-abundance matrix
guild_membership   : MAGs x ecological guilds membership matrix
primary_membership : MAGs x primary-production trait matrix

The functional matrices are abundance weighted, not sample-level presence/absence.
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon


def hellinger(x):
    x = np.asarray(x, dtype=float)
    totals = x.sum(axis=1, keepdims=True)
    rel = np.divide(x, totals, out=np.zeros_like(x), where=totals != 0)
    return np.sqrt(rel)


def bray_curtis(x):
    return squareform(pdist(x, metric="braycurtis"))


def mean_distance_to_others(d):
    return d.sum(axis=1) / (d.shape[0] - 1)


def compare(A, guild_membership, primary_membership):
    guild_abundance = A @ guild_membership
    primary_abundance = A @ primary_membership

    d_tax = bray_curtis(hellinger(A))
    d_guild = bray_curtis(hellinger(guild_abundance))
    d_primary = bray_curtis(hellinger(primary_abundance))

    tax_mean = mean_distance_to_others(d_tax)
    guild_mean = mean_distance_to_others(d_guild)
    primary_mean = mean_distance_to_others(d_primary)

    tests = {
        "taxonomic_gt_guild": wilcoxon(tax_mean, guild_mean, alternative="greater"),
        "taxonomic_gt_primary": wilcoxon(tax_mean, primary_mean, alternative="greater"),
    }

    return {
        "distance_taxonomic": d_tax,
        "distance_guild": d_guild,
        "distance_primary": d_primary,
        "sample_mean_taxonomic": tax_mean,
        "sample_mean_guild": guild_mean,
        "sample_mean_primary": primary_mean,
        "tests": tests,
    }
