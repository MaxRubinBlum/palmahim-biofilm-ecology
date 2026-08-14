#!/usr/bin/env python3
"""Build machine-readable analysis inputs from Supplementary Tables S3 and S4.

This script makes the published supplementary tables the canonical public input
for downstream custom analyses. It avoids committing duplicate ~100-kB wide
matrices while preserving an exact, executable transformation.

Required
--------
--s3  Supplementary Table S3 (ATLAS summary; sheet `genome_quality`)
--s4  Supplementary Table S4 (ecological annotation)

Outputs
-------
mag_abundance.tsv                 samples x MAG relative abundance, fractions
mag_taxonomy.tsv                  MAG taxonomy
curated_mag_traits.tsv            final curated S4 trait table
strict_guild_membership.tsv       four mutually exclusive ecological guilds
primary_production_membership.tsv three abundance-weighted production traits
sample_metadata.tsv               habitat assignments and S5 inclusion flag

The taxonomic-functional analysis uses 18 biofilms: AnemPM22 is retained in the
published abundance table but marked for exclusion by the metadata flag.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


SAMPLE_HABITATS = {
    "AnemPM22": "anemone",
    "Bio1": "carbonate", "Bio2": "carbonate", "Bio3": "carbonate",
    "Bio4": "carbonate", "Bio6": "carbonate", "Bio7": "carbonate",
    "Bio8": "carbonate", "Bio9": "carbonate",
    "CrabPM21": "crab_carapace",
    "EGGPM35": "shark_egg_capsule", "EGGPM36": "shark_egg_capsule",
    "PBB": "ARMS", "PBT": "ARMS",
    "PM27": "plastic", "PM28": "plastic", "PM29": "plastic",
    "PM30": "plastic", "PM39": "plastic",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s3", required=True, help="Supplementary Table S3 xlsx")
    ap.add_argument("--s4", required=True, help="Supplementary Table S4 xlsx")
    ap.add_argument("--outdir", default="data")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # S3 contains one row per MAG. Relative-abundance columns are explicit (%)
    # and are converted to fractions for the repository analysis convention.
    s3 = pd.read_excel(args.s3, sheet_name="genome_quality")
    if "MAG" not in s3.columns:
        raise ValueError("S3 must contain a MAG column")
    s3 = s3.set_index("MAG", drop=False)

    abundance_cols = [c for c in s3.columns if str(c).endswith(" relative abundance (%)")]
    if len(abundance_cols) != 19:
        raise ValueError(f"Expected 19 relative-abundance columns in S3, found {len(abundance_cols)}")
    samples = [str(c).removesuffix(" relative abundance (%)") for c in abundance_cols]
    unknown = [s for s in samples if s not in SAMPLE_HABITATS]
    if unknown:
        raise ValueError(f"Unrecognized samples in S3: {unknown}")

    abundance_pct = s3[abundance_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    col_sums = abundance_pct.sum(axis=0).to_numpy()
    if not np.allclose(col_sums, 100.0, atol=1e-6):
        raise ValueError(f"S3 normalized abundance columns do not sum to 100%: {col_sums}")
    abundance = (abundance_pct / 100.0).T
    abundance.index = samples
    abundance.columns = s3.index
    abundance.index.name = "sample"
    abundance.to_csv(out / "mag_abundance.tsv", sep="\t")

    tax_map = {
        "Domain": "Domain", "Phylum": "phylum", "Class": "class",
        "Order": "Order", "Family": "family", "Genus": "genus", "Species": "species"
    }
    tax = s3[["MAG"] + list(tax_map)].rename(columns=tax_map).set_index("MAG")
    tax.to_csv(out / "mag_taxonomy.tsv", sep="\t")

    # S4 has a category row followed by the true column-name row.
    s4 = pd.read_excel(args.s4, sheet_name=0, header=1)
    s4 = s4.loc[:, ~s4.columns.astype(str).str.startswith("Unnamed:")]
    if "MAG" not in s4.columns or "Ecological_guild" not in s4.columns:
        raise ValueError("S4 header was not recognized; expected MAG and Ecological_guild")
    s4 = s4.dropna(subset=["MAG"]).set_index("MAG", drop=False)
    s4.to_csv(out / "curated_mag_traits.tsv", sep="\t", index=False)

    # All S3 MAGs should have an S4 row; fail loudly rather than silently
    # changing the taxonomic matrix when a trait table is incomplete.
    missing = sorted(set(s3.index) - set(s4.index))
    if missing:
        raise ValueError(f"S4 is missing {len(missing)} S3 MAGs; first examples: {missing[:10]}")

    guild = pd.DataFrame(index=s4.index)
    guild["methanotroph"] = pd.to_numeric(s4["Methanotroph"], errors="coerce").fillna(0).astype(int)
    guild["sulfur_oxidizing_autotroph"] = pd.to_numeric(s4["Sulfur_oxidizing_autotroph"], errors="coerce").fillna(0).astype(int)
    guild["other_autotroph"] = s4["Ecological_guild"].eq("Other autotroph").astype(int)
    guild["other_community_member"] = s4["Ecological_guild"].eq("Other").astype(int)
    guild.index.name = "MAG"
    if not (guild.sum(axis=1) == 1).all():
        bad = guild.index[guild.sum(axis=1) != 1].tolist()[:10]
        raise ValueError(f"Strict ecological guilds are not mutually exclusive for: {bad}")
    guild.to_csv(out / "strict_guild_membership.tsv", sep="\t")

    primary = pd.DataFrame(index=s4.index)
    primary["methanotrophy"] = pd.to_numeric(s4["Methanotroph"], errors="coerce").fillna(0).astype(int)
    primary["sulfur_oxidizing_autotrophy"] = pd.to_numeric(s4["Sulfur_oxidizing_autotroph"], errors="coerce").fillna(0).astype(int)
    primary["autotrophic_carbon_fixation"] = pd.to_numeric(s4["Autotrophic_carbon_fixation"], errors="coerce").fillna(0).astype(int)
    primary.index.name = "MAG"
    primary.to_csv(out / "primary_production_membership.tsv", sep="\t")

    meta = pd.DataFrame({
        "sample": samples,
        "habitat": [SAMPLE_HABITATS[s] for s in samples],
        "include_taxonomic_functional_beta_diversity": [0 if s == "AnemPM22" else 1 for s in samples],
    })
    meta.to_csv(out / "sample_metadata.tsv", sep="\t", index=False)

    print(f"Wrote reproducibility inputs for {len(s3)} MAGs and {len(samples)} samples to {out}")
    print("Taxonomic-functional comparison samples:", int(meta.include_taxonomic_functional_beta_diversity.sum()))

if __name__ == "__main__":
    main()
