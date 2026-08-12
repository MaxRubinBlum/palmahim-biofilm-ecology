"""Abundance-weighted vitamin-provider summaries."""

import pandas as pd


def b12_status(n_steps):
    if n_steps == 7:
        return "complete"
    if n_steps == 6:
        return "near-complete"
    if n_steps >= 1:
        return "partial"
    return "absent"


def provider_abundance(habitat_mean_abundance, complete_provider_mask):
    """Sum abundance of MAGs carrying a complete pathway."""
    return habitat_mean_abundance.loc[:, complete_provider_mask].sum(axis=1)


def aggregate_providers_by_order(provider_table, abundance_col="abundance"):
    return provider_table.groupby("Order")[abundance_col].sum().sort_values(ascending=False)
