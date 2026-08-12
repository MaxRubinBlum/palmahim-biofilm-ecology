#!/usr/bin/env bash
set -euo pipefail

# Phylogenomic provenance for the Palmahim seep biofilm manuscript.
#
# Fig. 5: Methylococcales subtree pruned from the broader GTDB-Tk bac120 tree
# used for Supplementary Fig. 2. The final pruned Newick tree is supplied under
# provenance/phylogenomics/gtdbtk/.
#
# Supplementary Figs. 4-6: lineage-specific GToTree v1.7.10 analyses using
# the Gammaproteobacteria HMM set (172 targets). The commands below are copied
# from the original run logs supplied with the manuscript analysis.

# Supplementary Fig. 5 — QPIN01/MMG2
# 55 genomes supplied; 53 retained after GToTree filtering.
GToTree -f filenames.txt -H Gammaproteobacteria -j 32 -o GtoTree

# Supplementary Fig. 6 — CAJXQU01
# 20 genomes supplied; all 20 retained.
GToTree -f filenames.txt -H Gammaproteobacteria -j 24 -o GtoTree2

# Supplementary Fig. 4 — Methyloprofundus
# 88 genomes supplied; 87 retained after GToTree filtering.
GToTree -f genomes.list -H Gammaproteobacteria -j 24 -o GtoTree

# GToTree default processing recorded in the run logs:
# - marker set: Gammaproteobacteria, 172 targets
# - gene-length filter: within 20% of the median for each gene set
# - genomes with <50% of targeted SCGs removed
# - concatenated amino-acid alignment
# - FastTreeMP v2.1.11, JTT + CAT (20 rate categories)
# - SH-like local support based on 1000 resamples
#
# Original run logs are retained under provenance/phylogenomics/gtotree/.
