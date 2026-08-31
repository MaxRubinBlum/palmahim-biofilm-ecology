# Synthetic reviewer demo

This directory contains a **small simulated dataset** for testing the public analysis code without access to the manuscript data. It exercises the central taxonomic-versus-functional beta-diversity workflow used for Supplementary Fig. 3 / Table S5.

The demo is only a software/reproducibility test. Its values have no biological interpretation and are not manuscript results.

## Inputs

- `mag_abundance.tsv` - six simulated samples x eight simulated MAGs; rows sum to 1.
- `guild_membership.tsv` - mutually exclusive strict ecological-guild assignments.
- `primary_membership.tsv` - primary-production trait assignments.
- `metadata.tsv` - three simulated habitat classes with two samples each.
- `expected_metrics.json` - numerical values used to validate the demo.

## Run

From the repository root, after installing the Python dependencies:

```bash
python scripts/run_demo.py
```

By default, outputs are written to `outputs/demo/`.

## Expected result

A successful run prints:

```text
PASS: synthetic reviewer demo reproduced expected outputs
```

and validates these summary values:

| Representation | Mean Bray-Curtis dissimilarity |
|---|---:|
| Taxonomic MAG composition | 0.11866717634372449 |
| Strict ecological guilds | 0.04663683547885813 |
| Primary-production traits | 0.03864436134709344 |

The one-sided paired Wilcoxon tests give `P = 0.015625` for both `taxonomic > guild` and `taxonomic > primary` in this simulated dataset.

The demo also writes the analysis matrices, distance tables, Wilcoxon table, and PNG/SVG beta-diversity figure produced by `02_taxonomic_functional_beta_diversity.py`.

## Optional R companion test

If R and `vegan` are installed, the demo-generated matrices can also be passed to the companion permutation workflow:

```bash
Rscript scripts/02b_beta_diversity_permutation_tests.R \
  outputs/demo/taxonomic_abundance.csv \
  outputs/demo/guild_abundance.csv \
  outputs/demo/primary_abundance.csv \
  outputs/demo/metadata_analysis_samples.csv \
  outputs/demo
```

GitHub Actions runs the Python demo automatically across supported Python versions and also runs the R companion workflow on Ubuntu.
