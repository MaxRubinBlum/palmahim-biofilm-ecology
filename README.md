# Palmahim seep biofilm ecology analyses

Custom analysis and figure-preparation code accompanying the manuscript **“Hydrocarbon seep biofilms share a common functional organization across diverse substrates.”**

This repository documents the custom analyses that are not fully reproduced by standard bioinformatics packages. Upstream metagenomic assembly, MAG reconstruction, annotation, read mapping and phylogenomics are described in the manuscript Methods and Supplementary Methods together with software versions and parameters.

## Repository scope

The repository covers:

| Script | Analysis | Manuscript output |
|---|---|---|
| `01_fig2_community_analysis.R` | Shannon diversity, Hellinger transformation, Bray–Curtis PCoA, PERMANOVA and PERMDISP | Fig. 2b,c |
| `02_taxonomic_functional_beta_diversity.py` | Taxonomic versus functional Bray–Curtis turnover and paired Wilcoxon tests | Supplementary Fig. 3; Supplementary Table 5 |
| `03_functional_redundancy.py` | Inverse-Simpson effective MAG number and core-function classification | Supplementary Table 5 |
| `04_alluvial_taxon_trait.py` | Taxon–trait Fisher tests, phi coefficients, abundance weighting and figure-ready links | Fig. 3 |
| `05_c1_network.py` | Metabolic compatibility, habitat overlap and C1 interaction scores; Cytoscape node/edge tables | Fig. 4 |
| `06_metatranscriptome_tpm.py` | Gene TPM and MAG-level transcript abundance | Figs. 6–7 |
| `07_heterotroph_summaries.py` | Abundance, MEROPS and CAZyme summaries | Supplementary Tables 6–8 |
| `08_vitamin_summaries.py` | Abundance-weighted vitamin providers and cobalamin completeness | Supplementary Figs. 7–8 |

## Data required

The scripts do not assume project-specific filenames. They require combinations of the following tabular data:

- sample × MAG relative-abundance matrix;
- sample × MAG raw-abundance matrix where needed;
- sample metadata with habitat/substrate assignments;
- MAG taxonomy;
- curated MAG ecological and metabolic trait assignments;
- gene-level RNA read counts, gene lengths and gene-to-MAG mapping for transcriptomic analyses;
- MAG-level MEROPS and CAZyme annotations;
- complete vitamin-pathway assignments.

Biological trait assignments should be generated using the operational definitions reported in Supplementary Table 1. The analysis scripts treat those curated classifications as input and do not reinterpret genome annotations.

## Minimal input schemas

See [`data/README.md`](data/README.md) for the required columns and table orientation for each analysis.

## Environment

Python dependencies are listed in [`environment.yml`](environment.yml). The community analysis uses R with `vegan`.

Create the Python environment with:

```bash
conda env create -f environment.yml
conda activate palmahim-biofilm-ecology
```

For R:

```r
install.packages("vegan")
```

## Reproducibility principle

```text
genome annotations
        ↓
curated ecological classifications
        ↓
sample-specific MAG abundance
        ↓
statistical / aggregation scripts
        ↓
figure- and table-ready data
        ↓
visualization
```

The same curated classifications should feed all manuscript figures, supplementary figures and supplementary tables.

## Network interpretation

The C1 network represents **potential metabolic handoffs** supported by genome-encoded metabolic complementarity and ecological co-occurrence. It does not directly measure metabolite exchange or carbon flux.

## Citation and archival release

Before publication, create a tagged release corresponding to the accepted manuscript and archive that release with Zenodo. The resulting DOI can then be added to the manuscript Code availability statement.

## License

Code is released under the MIT License. Data files retain the terms specified by their original repositories and the manuscript data-availability statement.
