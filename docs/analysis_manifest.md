# Analysis-to-manuscript manifest

This file maps each custom or manuscript-specific analysis to the code supplied in this repository. Standard third-party bioinformatic programs are cited and parameterized in the manuscript/Supplementary Methods; the repository records the custom aggregation/statistical steps and reference command templates needed to connect those outputs to manuscript figures and tables.

| Manuscript item | Code | Required data | Main outputs |
|---|---|---|---|
| Fig. 2a | `01a_fig2_order_abundance.py` | MAG abundance + taxonomy | top-20 order matrix; bubble-plot table |
| Fig. 2b,c | `01_fig2_community_analysis.R` | sample × MAG abundance + habitat metadata | Shannon tests, PCoA, PERMANOVA, PERMDISP |
| Supp. Fig. 1 | `11_ctd_figure_prep.py` | processed CTD/ROV timeline | CTD plot |
| Supp. Fig. 3 / Table 5 | `02_...py`, `02b_...R`, `03_...py` | MAG abundance + strict guild/primary matrices + taxonomy | distance matrices, Wilcoxon, PERMANOVA, Mantel, Procrustes, redundancy/core |
| Fig. 3 | `04_alluvial_taxon_trait.py` | MAG order, abundance, curated traits | taxon-trait link table with Fisher/phi/BH and visual scaling |
| Fig. 4 | `05_c1_network.py` | relative/raw abundance + C1 traits + habitat + taxonomy | complete edge table, display edge table, Cytoscape nodes |
| Figs. 6–7 | `09_metatranscriptome_mapping.sh`, `06_metatranscriptome_tpm.py` | RNA FASTQ/CDS counts + DNA abundance + traits | gene TPM, MAG TPM, integrated plotting matrix |
| Fig. 5 | GTDB-Tk bac120 provenance + supplied pruned Newick | broad GTDB-Tk bac120 tree | pruned Methylococcales subtree |
| Supp. Fig. 2 | GTDB-Tk bac120 workflow described in Methods | bacterial MAG FASTA + GTDB-Tk database | broad bacterial MAG tree |
| Supp. Fig. 4 | `10_phylogenomics_commands.sh` + original run log | 88 *Methyloprofundus* focal/reference genomes | GToTree tree; 87 genomes retained |
| Supp. Fig. 5 | `10_phylogenomics_commands.sh` + original run log | 55 QPIN01/MMG2 focal/reference genomes | GToTree tree; 53 genomes retained |
| Supp. Fig. 6 | `10_phylogenomics_commands.sh` + original run log | 20 CAJXQU01 focal/reference genomes | GToTree tree; 20 genomes retained |
| Tables 6–8 | `07_heterotroph_summaries.py` | integrated heterotroph profile | MAG and order-level summaries |
| Supp. Figs. 7–8 | `08_vitamin_summaries.py` | habitat abundance + vitamin calls + taxonomy/function | provider tables, B12 completeness/function tables |

## Items intentionally not reimplemented

ATLAS, metaSPAdes, VAMB, dRep, CheckM2, GTDB-Tk, Minimap2, Prodigal, RAST, eggNOG-mapper, METABOLIC, dbCAN2, MEROPS, MacSyFinder, QSAP, Bowtie2, SAMtools, featureCounts, GToTree and iTOL are third-party tools. Their versions and parameters belong in the Methods/Supplementary Table 2; this repository does not copy their source code.

## Important interpretation rules

- Functional matrices are abundance weighted, not sample-level presence/absence.
- Sulfur-oxidizing autotrophy uses the strict curated ecological classification, not generic sulfur-metabolism genes.
- The C1 network represents potential metabolic handoffs, not measured metabolite flux.
- CLR correlations annotate support for C1 edges but do not define biological edge inclusion.
- Supplementary Fig. 3 uses one mean dissimilarity value per biofilm for paired Wilcoxon tests, avoiding pseudoreplication of all pairwise distances.
