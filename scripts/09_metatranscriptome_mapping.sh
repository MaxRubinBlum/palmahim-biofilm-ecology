#!/usr/bin/env bash
set -euo pipefail
# Reference workflow for upstream RNA mapping/count processing.
#
# IMPORTANT: this file documents the analysis logic, not a claim that every
# historical library used identical featureCounts flags. The exact featureCounts
# commands recovered from original output headers are preserved under
# provenance/metatranscriptomics/featurecounts_commands.txt.
#
# Data required: paired quality-filtered RNA FASTQ files, complete MAG catalog,
# and Prodigal CDS GFF. Replace placeholders with your paths/sample names.

bowtie2-build MAG_catalog.fasta MAG_catalog
bowtie2 --very-sensitive -x MAG_catalog \
  -1 SAMPLE_R1.fastq.gz -2 SAMPLE_R2.fastq.gz -S SAMPLE.sam

samtools view -b -q 20 -f 2 SAMPLE.sam | \
  samtools sort -o SAMPLE.q20.proper.bam
samtools index SAMPLE.q20.proper.bam

# Generic paired-end CDS counting example. Historical command lines differed
# between libraries; consult the provenance file before reporting exact flags.
featureCounts -p -t CDS -g ID \
  -a MAG_catalog.gff \
  -o SAMPLE.featureCounts.txt SAMPLE.q20.proper.bam

# TPM is calculated downstream with 06_metatranscriptome_tpm.py.
