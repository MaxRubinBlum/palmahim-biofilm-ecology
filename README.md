# Palmahim seep biofilm ecology analyses

Custom analysis and figure-preparation code accompanying the manuscript **“Hydrocarbon seep biofilms share a common functional organization across diverse substrates.”**

The repository records the custom statistical, aggregation and figure-preparation steps that connect curated genome annotations and abundance profiles to the manuscript figures and supplementary tables. It intentionally does **not** reproduce the source code of established third-party bioinformatic tools; their versions and parameters are reported in the Methods and Supplementary Table S2.

## Analysis map

| Script | Analysis | Manuscript output |
|---|---|---|
| `00_prepare_reproducibility_inputs.py` | Reconstruct machine-readable abundance, taxonomy and curated membership matrices from published Supplementary Tables S3–S4 | Shared downstream inputs |
| `00b_summarize_macsyfinder.py` | Summarize archived MacSyFinder system evidence while removing concatenated header rows | Annotation provenance |
| `01a_fig2_order_abundance.py` | MAG → order aggregation and top-20 selection | Fig. 2a |
| `01_fig2_community_analysis.R` | Shannon diversity, Hellinger/Bray–Curtis PCoA, PERMANOVA, PERMDISP | Fig. 2b,c |
| `02_taxonomic_functional_beta_diversity.py` | Taxonomic/guild/primary-production Bray–Curtis, paired Wilcoxon, distance summaries, boxplot | Supp. Fig. 3; Table S5 |
| `02b_beta_diversity_permutation_tests.R` | Replicated-substrate PERMANOVA, Mantel, Procrustes/PROTEST | Table S5 |
| `03_functional_redundancy.py` | Carrier counts, inverse-Simpson effective MAG number, core ≥1% classification | Table S5 |
| `04_alluvial_taxon_trait.py` | Top taxa, Fisher tests, phi, BH correction, abundance-weighted links | Fig. 3 |
| `05_c1_network.py` | C1 compatibility, shared habitat, interaction score, CLR support, Cytoscape node/edge tables | Fig. 4 |
| `09_metatranscriptome_mapping.sh` + `06_metatranscriptome_tpm.py` | RNA mapping/count reference workflow; gene/MAG TPM and integrated trait matrix | Figs. 6–7 |
| `07_heterotroph_summaries.py` | MAG abundance and order-level MEROPS/CAZyme summaries | Tables S6–S8 |
| `08_vitamin_summaries.py` | Habitat-weighted vitamin providers and B12 completeness/function | Supp. Figs. 7–8 |
| `10_phylogenomics_commands.sh` | Exact lineage-specific GToTree commands; GTDB-Tk tree provenance documented separately | Fig. 5; Supp. Figs. 2, 4–6 |
| `11_ctd_figure_prep.py` | Plotting a processed CTD/ROV timeline | Supp. Fig. 1 |

See [`docs/analysis_manifest.md`](docs/analysis_manifest.md) for a figure-by-figure provenance map and [`docs/reproducibility_validation.md`](docs/reproducibility_validation.md) for numerical validation against Supplementary Table S5.

## Canonical public inputs

The central custom analyses can be reconstructed from the manuscript's published Supplementary Tables S3 and S4 rather than from project-specific intermediate filenames. S3 contains the final MAG catalogue, GTDB taxonomy and sample-level relative abundance; S4 contains the final curated ecological assignments used throughout the downstream analyses.

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

`mag_abundance.tsv` uses fractions (0–1). The input builder verifies that the 19 abundance columns in S3 each sum to 100% before conversion. `sample_metadata.tsv` retains all samples but marks `AnemPM22` for exclusion from the 18-biofilm taxonomic-functional comparison.

These generated matrices are deliberately not duplicated in version control: the published S3/S4 workbooks are the canonical data source and the transformation into analysis-ready TSV files is executable and versioned here.

## Environment

```bash
conda env create -f environment.yml
conda activate palmahim-biofilm-ecology
```

R analyses require `vegan`; figure preparation may use `ggplot2`.

```r
install.packages(c("vegan", "ggplot2"))
```

## Reproduce Supplementary Fig. 3 / Table S5

```bash
python scripts/02_taxonomic_functional_beta_diversity.py \
  --abundance data/generated/mag_abundance.tsv \
  --guild-membership data/generated/strict_guild_membership.tsv \
  --primary-membership data/generated/primary_production_membership.tsv \
  --metadata data/generated/sample_metadata.tsv \
  --outdir outputs/taxonomic_functional
```

The Python analysis exports the exact taxonomic, guild and primary-production matrices consumed by the companion `vegan` workflow:

```bash
Rscript scripts/02b_beta_diversity_permutation_tests.R \
  outputs/taxonomic_functional/taxonomic_abundance.csv \
  outputs/taxonomic_functional/guild_abundance.csv \
  outputs/taxonomic_functional/primary_abundance.csv \
  outputs/taxonomic_functional/metadata_analysis_samples.csv \
  outputs/taxonomic_functional
```

Functional redundancy/core metrics can be regenerated from the same inputs:

```bash
python scripts/03_functional_redundancy.py \
  --abundance data/generated/mag_abundance.tsv \
  --traits data/generated/curated_mag_traits.tsv \
  --taxonomy data/generated/mag_taxonomy.tsv \
  --metadata data/generated/sample_metadata.tsv \
  --out outputs/taxonomic_functional/redundancy.csv
```

For the exact S5 redundancy rows, select the curated columns `Methanotroph`, `Sulfur_oxidizing_autotroph`, `Autotrophic_carbon_fixation`, `CBB_I`, `CBB_II` and `rTCA` from `curated_mag_traits.tsv`.

## Provenance of upstream annotations

The curated ecological table is the authoritative downstream input. Its assignments were assembled using the operational definitions in Supplementary Table S1 from METABOLIC, additional functional HMMs, QSAP, MacSyFinder, MEROPS and dbCAN/CAZyme evidence, with expert review where appropriate. The repository does not falsely represent this curation step as a completely automated classifier.

`provenance/source_files_manifest.tsv` records the exact supplied origin filenames, file sizes and SHA-256 checksums. Small, inspectable provenance outputs are retained in GitHub; large annotation and transcriptomic origin files are intended for the archival publication deposit (for example Zenodo) rather than duplication in the code repository. A compact MAG-level QSAP class matrix is included, and `00b_summarize_macsyfinder.py` reconstructs compact MacSyFinder summaries from the archived combined output.

Historical featureCounts commands recovered from the original count-file headers are stored in `provenance/metatranscriptomics/featurecounts_commands.txt`. Because the historical commands differed between libraries, `scripts/09_metatranscriptome_mapping.sh` is explicitly a reference workflow rather than a claim that all libraries used identical featureCounts flags.

## Reproducibility principle

```text
raw/third-party annotation outputs
        ↓
operational definitions + documented curation
        ↓
published Supplementary Tables S3–S4
        ↓
versioned input-construction script
        ↓
custom statistical / aggregation scripts
        ↓
figure- and table-ready outputs
```

## Phylogenomic provenance

Phylogenies have two distinct provenance routes. The broad bacterial tree in Supplementary Fig. 2 is the GTDB-Tk bac120 phylogeny, and the Methylococcales tree in Fig. 5 is a subtree pruned from that GTDB-Tk tree. The exact pruned Fig. 5 Newick file is included in `provenance/phylogenomics/gtdbtk/`. The lineage-specific *Methyloprofundus*, QPIN01/MMG2 and CAJXQU01 trees in Supplementary Figs. 4–6 were reconstructed with GToTree v1.7.10; their original run logs and exact commands are retained in `provenance/phylogenomics/`.

## C1 network interpretation

The C1 network represents **potential metabolic handoffs** supported by genome-encoded metabolic complementarity and ecological co-occurrence. It does not directly measure metabolite exchange or carbon flux. CLR correlations are supporting annotations rather than edge-inclusion criteria.

## Archival release

At acceptance/publication, create a tagged release matching the accepted manuscript and archive that release together with the larger provenance files. Add the resulting archival DOI to the manuscript Code availability statement.

## License

Code is released under the MIT License. Data files retain the terms specified by their original repositories and the manuscript data-availability statement.
