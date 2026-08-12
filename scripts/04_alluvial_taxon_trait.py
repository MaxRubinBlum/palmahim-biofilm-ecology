"""Prepare abundance-weighted taxon-trait links for Fig. 3.

Data required
-------------
One row per MAG containing:
    Order
    abundance
    binary trait columns
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def taxon_trait_statistics(mag_table, taxa, traits):
    rows = []

    for taxon in taxa:
        in_taxon = mag_table["Order"].eq(taxon)

        for trait in traits:
            has_trait = mag_table[trait].astype(bool)

            a = int((in_taxon & has_trait).sum())
            b = int((in_taxon & ~has_trait).sum())
            c = int((~in_taxon & has_trait).sum())
            d = int((~in_taxon & ~has_trait).sum())

            _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")

            denom = np.sqrt((a+b)*(c+d)*(a+c)*(b+d))
            phi = ((a*d) - (b*c)) / denom if denom else np.nan

            weighted_abundance = mag_table.loc[in_taxon & has_trait, "abundance"].sum()

            rows.append((taxon, trait, phi, p, weighted_abundance))

    out = pd.DataFrame(rows, columns=["taxon", "trait", "phi", "p", "abundance"])
    out["q"] = multipletests(out["p"], method="fdr_bh")[1]

    # Fig. 3 shows positive relationships; biological traits are curated upstream.
    out = out[out["phi"] > 0].copy()
    out["ribbon_width"] = np.sqrt(out["abundance"])
    out["ribbon_opacity"] = out["phi"]
    return out
