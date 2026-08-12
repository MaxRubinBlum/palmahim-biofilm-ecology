"""Gene- and MAG-level TPM for Figs. 6 and 7.

Data required: gene_id, MAG, length_bp, read_count.
"""

import pandas as pd


def calculate_tpm(genes):
    genes = genes.copy()
    genes["length_kb"] = genes["length_bp"] / 1000.0
    genes["RPK"] = genes["read_count"] / genes["length_kb"]
    scale = genes["RPK"].sum() / 1e6
    genes["TPM"] = genes["RPK"] / scale
    return genes


def mag_level_tpm(genes_with_tpm):
    mag_tpm = genes_with_tpm.groupby("MAG")["TPM"].sum()
    relative_mag_tpm = 100 * mag_tpm / mag_tpm.sum()
    return pd.DataFrame({"MAG_TPM": mag_tpm, "relative_MAG_TPM_pct": relative_mag_tpm})
