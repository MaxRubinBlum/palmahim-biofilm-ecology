# Palmahim seep biofilm ecology analyses

Custom analysis and figure-preparation code accompanying the manuscript **“Hydrocarbon seep biofilms share a common functional organization across diverse substrates.”**

The repository records the custom statistical, aggregation and figure-preparation steps that connect curated genome annotations and abundance profiles to the manuscript figures and supplementary tables. It intentionally does **not** reproduce the source code of established third-party bioinformatic tools; their versions and parameters are reported in the Methods and Supplementary Table 2.

## Analysis map

| Script | Analysis | Manuscript output |
|---|---|---|
| `01a_fig2_order_abundance.py` | MAG → order aggregation and top-20 selection | Fig. 2a |
| `01_fig2_community_analysis.R` | Shannon diversity, Hellinger/Bray–Curtis PCoA, PERMANOVA, PERMDISP | Fig. 2b,c |
| `02_taxonomic_functional_beta_diversity.py` | Taxonomic/guild/primary-production Bray–Curtis, paired Wilcoxon, distance summaries, boxplot | Supp. Fig. 3; Table 5 |
| `02b_beta_diversity_permutation_tests.R` | Replicated-substrate PERMANOVA, Mantel, Procrustes/PROTEST | Table 5 |
| `03_functional_redundancy.py` | Carrier counts, inverse-Simpson effective MAG number, core ≥1% classification | Table 5 |
| `04_alluvial_taxon_trait.py` | Top taxa, Fisher tests, phi, BH correction, abundance-weighted links | Fig. 3 |
| `05_c1_network.py` | C1 compatibility, shared habitat, exact interaction score, CLR support, Cytoscape node/edge tables | Fig. 4 |
| `09_metatranscriptome_mapping.sh` + `06_metatranscriptome_tpm.py` | RNA mapping/count reference commands; gene/MAG TPM and integrated trait matrix | Figs. 6–7 |
| `07_heterotroph_summaries.py` | MAG abundance and order-level MEROPS/CAZyme summaries | Tables 6–8 |
| `08_vitamin_summaries.py` | Habitat-weighted vitamin providers and B12 completeness/function | Supp. Figs. 7–8 |
| `10_phylogenomics_commands.sh` | Exact lineage-specific GToTree commands; GTDB-Tk tree provenance documented separately | Fig. 5; Supp. Figs. 2, 4–6 |
| `11_ctd_figure_prep.py` | Plotting a processed CTD/ROV timeline | Supp. Fig. 1 |

See [`docs/analysis_manifest.md`](docs/analysis_manifest.md) for a figure-by-figure provenance map.

## Data required

Scripts use biological data types rather than project-specific filenames:

- sample × MAG relative-abundance matrix;
- sample × MAG raw-abundance/coverage matrix where required;
- sample metadata with habitat/substrate assignments;
- MAG taxonomy;
- curated MAG ecological/metabolic trait assignments;
- gene-level RNA counts, gene lengths and gene-to-MAG mapping;
- MEROPS and CAZyme annotations;
- complete vitamin-pathway assignments;
- processed CTD/ROV timeline for Supplementary Fig. 1.

Biological trait assignments must follow the operational definitions in Supplementary Table 1. The scripts treat those classifications as input; they do not independently reinterpret genome annotations.

Detailed schemas are in [`data/README.md`](data/README.md).

## Environment

```bash
conda env create -f environment.yml
conda activate palmahim-biofilm-ecology
```

R analyses require `vegan`; figure preparation may use `ggplot2`.

```r
install.packages(c("vegan", "ggplot2"))
```

## Example workflow

```bash
python scripts/02_taxonomic_functional_beta_diversity.py \
  --abundance data/mag_abundance.csv \
  --guild-membership data/strict_guild_membership.csv \
  --primary-membership data/primary_production_membership.csv \
  --metadata data/sample_metadata.csv \
  --outdir outputs/taxonomic_functional

Rscript scripts/02b_beta_diversity_permutation_tests.R \
  data/mag_abundance.csv \
  data/guild_abundance.csv \
  data/primary_abundance.csv \
  data/sample_metadata.csv \
  outputs/taxonomic_functional
```

## Reproducibility principle

```text
genome annotations
        ↓
curated ecological classifications
        ↓
sample-specific MAG abundance
        ↓
statistical and aggregation scripts
        ↓
figure- and table-ready data
        ↓
visualization
```

## Phylogenomic provenance

Phylogenies have two distinct provenance routes. The broad bacterial tree in Supplementary Fig. 2 is the GTDB-Tk bac120 phylogeny, and the Methylococcales tree in Fig. 5 is a subtree pruned from that GTDB-Tk tree. The exact pruned Fig. 5 Newick file is included in `provenance/phylogenomics/gtdbtk/`. The lineage-specific *Methyloprofundus*, QPIN01/MMG2 and CAJXQU01 trees in Supplementary Figs. 4–6 were reconstructed with GToTree v1.7.10; their original run logs and exact commands are retained in `provenance/phylogenomics/`.

## C1 network interpretation

The C1 network represents **potential metabolic handoffs** supported by genome-encoded metabolic complementarity and ecological co-occurrence. It does not directly measure metabolite exchange or carbon flux. CLR correlations are supporting annotations rather than edge-inclusion criteria.

## Archival release

At acceptance/publication, create a tagged release matching the accepted manuscript and archive that release with Zenodo. Add the DOI to the manuscript Code availability statement.

## License

Code is released under the MIT License. Data files retain the terms specified by their original repositories and the manuscript data-availability statement.
