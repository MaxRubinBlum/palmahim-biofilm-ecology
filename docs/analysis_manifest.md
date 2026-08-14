# Analysis-to-manuscript manifest

This file maps each custom or manuscript-specific analysis to the code and data provenance supplied in this repository. Standard third-party bioinformatic programs are cited and parameterized in the manuscript/Supplementary Methods; the repository records the manuscript-specific aggregation, statistical and figure-preparation steps.

## Canonical shared inputs

Published Supplementary Tables S3 and S4 are the canonical public source for the central abundance-weighted custom analyses. `00_prepare_reproducibility_inputs.py` reconstructs the sample × MAG abundance matrix, MAG taxonomy, curated trait table, strict ecological-guild membership, primary-production membership and sample metadata. It verifies that all 19 S3 relative-abundance columns sum to 100%, converts them to fractions, checks that all S3 MAGs occur in S4 and checks that strict ecological guilds are mutually exclusive.

| Manuscript item | Code | Required data | Main outputs |
|---|---|---|---|
| Shared input reconstruction | `00_prepare_reproducibility_inputs.py` | published Supplementary Tables S3–S4 | analysis-ready abundance, taxonomy, curated traits, guild/production membership, metadata |
| Fig. 2a | `01a_fig2_order_abundance.py` | MAG abundance + taxonomy | top-20 order matrix; bubble-plot table |
| Fig. 2b,c | `01_fig2_community_analysis.R` | sample × MAG abundance + habitat metadata | Shannon tests, PCoA, PERMANOVA, PERMDISP |
| Supp. Fig. 1 | `11_ctd_figure_prep.py` | processed CTD/ROV timeline | CTD plot |
| Supp. Fig. 3 / Table S5 | `02_taxonomic_functional_beta_diversity.py`, `02b_beta_diversity_permutation_tests.R`, `03_functional_redundancy.py` | S3/S4-derived abundance + strict guild/primary matrices + taxonomy | exact analysis matrices, distance matrices, paired Wilcoxon, replicated-substrate PERMANOVA, Mantel, Procrustes, redundancy/core |
| Fig. 3 | `04_alluvial_taxon_trait.py` | MAG order, abundance, curated S4 traits | taxon-trait link table with Fisher/phi/BH and visual scaling |
| Fig. 4 | `05_c1_network.py` | relative/raw abundance + curated C1 traits + habitat + taxonomy | complete edge table, display edge table, Cytoscape nodes |
| Figs. 6–7 | `09_metatranscriptome_mapping.sh`, `06_metatranscriptome_tpm.py` | RNA FASTQ/CDS counts + DNA abundance + traits | gene TPM, MAG TPM, integrated plotting matrix |
| Fig. 5 | GTDB-Tk bac120 provenance + supplied pruned Newick | broad GTDB-Tk bac120 tree | pruned Methylococcales subtree |
| Supp. Fig. 2 | GTDB-Tk bac120 workflow described in Methods | bacterial MAG FASTA + GTDB-Tk database | broad bacterial MAG tree |
| Supp. Fig. 4 | `10_phylogenomics_commands.sh` + original run log | 88 *Methyloprofundus* focal/reference genomes | GToTree tree; 87 genomes retained |
| Supp. Fig. 5 | `10_phylogenomics_commands.sh` + original run log | 55 QPIN01/MMG2 focal/reference genomes | GToTree tree; 53 genomes retained |
| Supp. Fig. 6 | `10_phylogenomics_commands.sh` + original run log | 20 CAJXQU01 focal/reference genomes | GToTree tree; 20 genomes retained |
| Tables S6–S8 | `07_heterotroph_summaries.py` | integrated heterotroph profile; percent abundance by default | MAG and order-level summaries |
| Supp. Figs. 7–8 | `08_vitamin_summaries.py` | habitat abundance + vitamin calls + taxonomy/function | provider tables, B12 complete-provider function, abundance-weighted B12 completeness classes |
| Annotation provenance | `00b_summarize_macsyfinder.py`, compact QSAP matrix, source manifest | archived raw annotation outputs | inspectable provenance derivatives and checksums |

## Numerically validated manuscript outputs

The S3/S4-derived inputs reproduce the Supplementary Table S5 mean Bray–Curtis dissimilarities exactly: 0.6210749038914138 (taxonomic), 0.13066208127467652 (strict ecological guilds) and 0.069498470943168 (strict primary-production traits). Both one-sided paired Wilcoxon comparisons give P = 3.814697265625e-06. Carrier counts, core status and median inverse-Simpson effective MAG numbers for Methanotroph, Sulfur_oxidizing_autotroph, Autotrophic_carbon_fixation, CBB_I, CBB_II and rTCA also reproduce S5 exactly. See `docs/reproducibility_validation.md`.

## Upstream provenance boundary

The final curated ecological assignments in Supplementary Table S4 are authoritative downstream inputs. Their derivation uses operational definitions from Supplementary Table S1 applied to METABOLIC, functional HMM, QSAP, MacSyFinder, MEROPS and dbCAN/CAZyme evidence, with expert review where appropriate. The repository therefore distinguishes annotation provenance from downstream computational reproducibility rather than implying that all S4 assignments were produced by one automated classifier.

`provenance/source_files_manifest.tsv` records exact supplied source filenames, sizes and SHA-256 checksums. Large raw annotation/transcriptomic files should accompany the archival publication release rather than being duplicated in GitHub.

## Items intentionally not reimplemented

ATLAS, metaSPAdes, VAMB, dRep, CheckM2, GTDB-Tk, Minimap2, Prodigal, RAST, eggNOG-mapper, METABOLIC, dbCAN2, MEROPS, MacSyFinder, QSAP, Bowtie2, SAMtools, featureCounts, GToTree and iTOL are third-party tools. Their versions and parameters belong in the Methods/Supplementary Table S2; this repository does not copy their source code.

## Important interpretation rules

- Functional matrices are abundance weighted, not sample-level presence/absence.
- Supplementary Fig. 3/Table S5 uses 18 biofilms; AnemPM22 is excluded by an explicit metadata flag.
- The singleton crab-carapace biofilm remains in the 18-sample distance/Wilcoxon analysis but is excluded from replicated-substrate PERMANOVA.
- Sulfur-oxidizing autotrophy uses the strict curated ecological classification, not generic sulfur-metabolism genes.
- The C1 network represents potential metabolic handoffs, not measured metabolite flux.
- CLR correlations annotate support for C1 edges but do not define biological edge inclusion.
- Supplementary Fig. 3 uses one mean dissimilarity value per biofilm for paired Wilcoxon tests, avoiding pseudoreplication of all pairwise distances.
