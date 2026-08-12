#!/usr/bin/env bash
set -euo pipefail
# Reference command template for the upstream RNA mapping/count workflow.
# Data required: paired quality-filtered RNA FASTQ files, complete MAG catalog,
# Prodigal CDS GFF. Replace placeholders with your paths/sample names.

bowtie2-build MAG_catalog.fasta MAG_catalog
bowtie2 --very-sensitive -x MAG_catalog -1 SAMPLE_R1.fastq.gz -2 SAMPLE_R2.fastq.gz -S SAMPLE.sam
samtools view -b -q 20 -f 2 SAMPLE.sam | samtools sort -o SAMPLE.q20.proper.bam
samtools index SAMPLE.q20.proper.bam
featureCounts -p -B -C -t CDS -g ID -a MAG_catalog.gff -o SAMPLE.featureCounts.txt SAMPLE.q20.proper.bam

# TPM is calculated downstream with 06_metatranscriptome_tpm.py.
