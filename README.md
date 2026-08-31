# Palmahim seep biofilm ecology analyses

[![reproducibility-smoke-test](https://github.com/MaxRubinBlum/palmahim-biofilm-ecology/actions/workflows/reproducibility.yml/badge.svg)](https://github.com/MaxRubinBlum/palmahim-biofilm-ecology/actions/workflows/reproducibility.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11-3.13](https://img.shields.io/badge/Python-3.11--3.13-blue.svg)](requirements.txt)

Custom analysis and figure-preparation code accompanying the manuscript **“Hydrocarbon seep biofilms share a common functional organization across diverse substrates.”**

This repository records the manuscript-specific statistical, aggregation, network and figure-preparation steps that connect curated genome annotations and abundance profiles to the reported figures and supplementary tables. It intentionally does **not** copy the source code of established third-party bioinformatic tools; their historical versions and parameters are reported in the manuscript Methods and Supplementary Table S2.

> **For editors and reviewers:** the repository includes source code, a self-contained simulated demo dataset, expected demo outputs, environment specifications, an MIT license, figure-by-figure provenance, numerical validation against Supplementary Table S5, and automated smoke tests. See [`docs/nature_software_checklist.md`](docs/nature_software_checklist.md) for a direct mapping to the Nature Research Code and Software Submission Checklist.

## Reviewer quick start

The quickest way to test the custom analysis code does **not** require the manuscript data.

```bash
git clone https://github.com/MaxRubinBlum/palmahim-biofilm-ecology.git
cd palmahim-biofilm-ecology

python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell:
# .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/run_demo.py
```

A successful run prints:

```text
PASS: synthetic reviewer demo reproduced expected outputs
```

and reports the measured runtime. On a current standard desktop/laptop, the synthetic demo is expected to complete in **well under one minute** after dependencies are installed. The automated GitHub Actions workflow runs the same demo on Ubuntu 24.04 under Python 3.11, 3.12 and 3.13.

The demo dataset and expected numerical values are documented in [`demo/README.md`](demo/README.md). It is simulated solely to test the code and has no biological interpretation.

## System requirements

### Reviewer/demo environment

- **Operating system:** automated tests use **Ubuntu 24.04**. The Python analysis scripts use standard cross-platform libraries and can also be run on current 64-bit Linux, macOS or Windows installations. Bash reference workflows are intended for Unix-like shells or WSL.
- **Python:** 3.11-3.13.
- **Python dependencies:** version-bounded requirements are listed in [`requirements.txt`](requirements.txt) and [`environment.yml`](environment.yml): NumPy 1.26-<3, pandas 2-<3, SciPy 1.11-<2, statsmodels 0.14-<1, Matplotlib 3.9-<4 and openpyxl 3.1-<4.
- **R companion analyses:** R 4.3-<5 with vegan 2.6-<3; ggplot2 3.4-<4 is included for figure preparation. The CI companion workflow tests R 4.4 on Ubuntu.
- **Hardware:** no non-standard hardware is required for the reviewer demo or downstream custom analyses. A conventional desktop/laptop is sufficient. Upstream metagenomic assembly, read mapping and genome annotation can require substantially more CPU, memory and storage, but those established third-party workflows are outside the small reviewer demo.

The reviewer/test environment above is provided for reproducibility of the public custom code. It should not be confused with **historical software provenance**: the exact versions used for the manuscript analyses are reported in Supplementary Table S2 and the Supplementary Methods.

## Installation

### Fast Python-only installation

```bash
python -m venv .venv
source .venv/bin/activate            # Linux/macOS
# .venv\Scripts\Activate.ps1         # Windows PowerShell
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Typical installation time is approximately **1-5 minutes** on a standard desktop with a normal broadband connection, depending mainly on package download speed.

### Full Conda environment, including R

```bash
conda env create -f environment.yml
conda activate palmahim-biofilm-ecology
```

A full Conda solve/install including the R environment typically takes approximately **5-15 minutes**, depending on platform, solver cache and network speed.

## Demo

The repository contains a six-sample, eight-MAG simulated dataset under [`demo/`](demo/) that exercises the central taxonomic-versus-functional beta-diversity workflow.

```bash
python scripts/run_demo.py
```

The wrapper executes `scripts/02_taxonomic_functional_beta_diversity.py`, verifies that all expected output files were produced, and checks the numerical summaries against [`demo/expected_metrics.json`](demo/expected_metrics.json).

Expected demo summary:

| Representation | Mean Bray-Curtis dissimilarity |
|---|---:|
| Taxonomic MAG composition | 0.11866717634372449 |
| Strict ecological guilds | 0.04663683547885813 |
| Primary-production traits | 0.03864436134709344 |

Both one-sided paired Wilcoxon comparisons in the simulated dataset return `P = 0.015625`.

With R and `vegan` installed, reviewers can also run the companion permutation workflow on the demo-generated matrices:

```bash
Rscript scripts/02b_beta_diversity_permutation_tests.R \
  outputs/demo/taxonomic_abundance.csv \
  outputs/demo/guild_abundance.csv \
  outputs/demo/primary_abundance.csv \
  outputs/demo/metadata_analysis_samples.csv \
  outputs/demo
```

## Reproducing the central manuscript analysis

The central custom analyses can be reconstructed from manuscript Supplementary Tables S3 and S4. During peer review these tables accompany the submission; after publication they are the canonical public data source.

Generate the machine-readable inputs with:

```bash
python scripts/00_prepare_reproducibility_inputs.py \
  --s3 Supplementary_Table_S3_ATLAS_Summary.xlsx \
  --s4 Supplementary_Table_S4_Ecological_Annotation.xlsx \
  --outdir data/generated
```

This produces:

```text
data/generated/mag_abundance.tsv
data/generated/mag_taxonomy.tsv
data/generated/curated_mag_traits.tsv
data/generated/strict_guild_membership.tsv
data/generated/primary_production_membership.tsv
data/generated/sample_metadata.tsv
```

`mag_abundance.tsv` uses fractions (0-1). The input builder verifies that the 19 abundance columns in S3 each sum to 100% before conversion, checks that all S3 MAGs occur in S4, and reconstructs the curated ecological-guild and primary-production matrices used downstream.

Reproduce Supplementary Fig. 3 / Table S5:

```bash
python scripts/02_taxonomic_functional_beta_diversity.py \
  --abundance data/generated/mag_abundance.tsv \
  --guild-membership data/generated/strict_guild_membership.tsv \
  --primary-membership data/generated/primary_production_membership.tsv \
  --metadata data/generated/sample_metadata.tsv \
  --outdir outputs/taxonomic_functional

Rscript scripts/02b_beta_diversity_permutation_tests.R \
  outputs/taxonomic_functional/taxonomic_abundance.csv \
  outputs/taxonomic_functional/guild_abundance.csv \
  outputs/taxonomic_functional/primary_abundance.csv \
  outputs/taxonomic_functional/metadata_analysis_samples.csv \
  outputs/taxonomic_functional

python scripts/03_functional_redundancy.py \
  --abundance data/generated/mag_abundance.tsv \
  --traits data/generated/curated_mag_traits.tsv \
  --taxonomy data/generated/mag_taxonomy.tsv \
  --metadata data/generated/sample_metadata.tsv \
  --out outputs/taxonomic_functional/redundancy.csv
```

For the exact S5 redundancy rows, select the curated columns `Methanotroph`, `Sulfur_oxidizing_autotroph`, `Autotrophic_carbon_fixation`, `CBB_I`, `CBB_II` and `rTCA` from `curated_mag_traits.tsv`.

See [`docs/reproducibility_validation.md`](docs/reproducibility_validation.md) for the numerical cross-check against the manuscript output.

## Analysis map

| Script / provenance item | Analysis | Manuscript output |
|---|---|---|
| `00_prepare_reproducibility_inputs.py` | Reconstruct machine-readable abundance, taxonomy and curated membership matrices from Supplementary Tables S3-S4 | Shared downstream inputs |
| `00b_summarize_macsyfinder.py` | Summarize archived MacSyFinder system evidence | Annotation provenance |
| `01a_fig2_order_abundance.py` | MAG -> order aggregation and top-20 selection | Fig. 2a |
| `01_fig2_community_analysis.R` | Shannon diversity, Hellinger/Bray-Curtis PCoA, PERMANOVA, PERMDISP | Fig. 2b,c |
| `02_taxonomic_functional_beta_diversity.py` | Taxonomic/guild/primary-production Bray-Curtis, paired Wilcoxon, distance summaries, boxplot | Supp. Fig. 3; Table S5 |
| `02b_beta_diversity_permutation_tests.R` | Replicated-substrate PERMANOVA, Mantel, Procrustes/PROTEST | Table S5 |
| `03_functional_redundancy.py` | Carrier counts, inverse-Simpson effective MAG number, core >=1% classification | Table S5 |
| `04_alluvial_taxon_trait.py` | Top taxa, Fisher tests, phi, BH correction, abundance-weighted links | Fig. 3 |
| `05_c1_network.py` | C1 compatibility, shared habitat, interaction score, CLR support, Cytoscape node/edge tables | Fig. 4 |
| `09_metatranscriptome_mapping.sh` + `06_metatranscriptome_tpm.py` | RNA mapping/count reference workflow; gene/MAG TPM and integrated trait matrix | Figs. 6-7 |
| `07_heterotroph_summaries.py` | MAG abundance and order-level MEROPS/CAZyme summaries | Tables S6-S8 |
| `08_vitamin_summaries.py` | Habitat-weighted vitamin providers and B12 completeness/function | Supp. Figs. 7-8 |
| `10_phylogenomics_commands.sh` | Exact lineage-specific GToTree commands; GTDB-Tk provenance documented separately | Fig. 5; Supp. Figs. 2, 4-6 |
| `11_ctd_figure_prep.py` | Processed CTD/ROV timeline plotting | Supp. Fig. 1 |

A figure-by-figure provenance map is provided in [`docs/analysis_manifest.md`](docs/analysis_manifest.md).

## Data and provenance

Detailed input schemas are in [`data/README.md`](data/README.md). The final curated ecological table is the authoritative downstream input; its assignments were assembled from the operational definitions in Supplementary Table S1 using METABOLIC, additional functional HMMs, QSAP, MacSyFinder, MEROPS and dbCAN/CAZyme evidence, with expert review where appropriate.

The repository does **not** represent this curation step as a fully automated classifier. Instead it separates:

```text
raw / third-party annotation outputs
        |
        v
operational definitions + documented curation
        |
        v
Supplementary Tables S3-S4
        |
        v
versioned input-construction script
        |
        v
custom statistical / aggregation scripts
        |
        v
figure- and table-ready outputs
```

`provenance/source_files_manifest.tsv` records source filenames, sizes and SHA-256 checksums. Small inspectable provenance derivatives are retained in GitHub. Larger annotation and transcriptomic intermediates are intended for the archival publication deposit rather than duplication in the code repository.

Historical featureCounts commands recovered from original output headers are stored in `provenance/metatranscriptomics/featurecounts_commands.txt`. Because historical libraries used differing command-line flags, `scripts/09_metatranscriptome_mapping.sh` is explicitly a **reference workflow**, not a claim that every library used identical featureCounts settings.

## Important interpretation rules

- Functional matrices are abundance weighted, not sample-level presence/absence.
- Supplementary Fig. 3/Table S5 uses 18 biofilms; `AnemPM22` is excluded by an explicit metadata flag.
- The singleton crab-carapace biofilm remains in the 18-sample distance/Wilcoxon analysis but is excluded from replicated-substrate PERMANOVA.
- Sulfur-oxidizing autotrophy uses the strict curated ecological classification, not generic sulfur-metabolism genes.
- The C1 network represents **potential metabolic handoffs**, not measured metabolite flux.
- CLR correlations annotate support for C1 edges but do not define biological edge inclusion.
- Supplementary Fig. 3 uses one mean dissimilarity value per biofilm for paired Wilcoxon tests, avoiding pseudoreplication of all pairwise distances.

## Automated validation

The GitHub Actions workflow [`.github/workflows/reproducibility.yml`](.github/workflows/reproducibility.yml) provides a reviewer-facing smoke test. It:

1. installs the documented Python dependencies;
2. compiles all Python scripts;
3. runs and numerically validates the simulated demo under Python 3.11, 3.12 and 3.13 on Ubuntu 24.04;
4. syntax-checks the shell and R scripts;
5. runs the R 4.4 / `vegan` companion permutation workflow on the demo matrices.

This CI validates software execution. It does not replace the manuscript-specific numerical validation documented in [`docs/reproducibility_validation.md`](docs/reproducibility_validation.md).

## License

Code is released under the [MIT License](LICENSE). Data files retain the terms specified by their original repositories and the manuscript data-availability statement.

## Citation

GitHub citation metadata are provided in [`CITATION.cff`](CITATION.cff). Until the manuscript is published, please cite the manuscript title and this repository; the accepted release will be archived and assigned a persistent DOI.

## Archival release

At acceptance/publication, create a tagged release matching the accepted manuscript and archive that release together with the larger provenance files (for example through Zenodo). Add the archival DOI to the final Code availability statement.
