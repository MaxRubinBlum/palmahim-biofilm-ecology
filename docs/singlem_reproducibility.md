# SingleM revision workflow

This workflow accompanies the assembly-independent prokaryotic community analysis of the 19 seep biofilms. It separates third-party marker profiling from custom recovery and diversity analyses. The scripts are a portable reconstruction of the workstation workflow, not verbatim copies of its files. Before tagging a release, verify the reconstructed commands against the original execution logs and compare the regenerated tables with the original results.

## Inputs and software

- Paired raw metagenomic reads: `SAMPLE_1.fastq.gz`, `SAMPLE_2.fastq.gz`. Use the same libraries as the ATLAS project; the raw reads are not the historical ATLAS QC-filtered reads. Confirm library provenance and record checksums in the archival manifest.
- The 19 original ATLAS assembly FASTAs named `SAMPLE.fasta` and the final 565 MAG FASTAs named `MAG*.fasta`. Use the explicit sample list in the shell script (no Bio5 or Bio10).
- SingleM 0.18.3, metapackage `S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb`. Obtain the appropriate metapackage separately; record its checksum and installed dependency versions. The workstation used an existing SingleM environment. See the [SingleM documentation](https://wwood.github.io/singlem/).
- Python with pandas and NumPy (including Generator.multivariate_hypergeometric); the repository requirements cover these dependencies. For the statistical step the historical environment used R 4.6.2, vegan 2.5-7 and permute 0.9-7. Record the exact Python environment used for the original rarefaction; fixed seeds do not guarantee identical random draws across software versions.
- Original manuscript abundance and taxonomy tables, generated from Supplementary Table S3 (and S4 using the existing input preparation script). Do not substitute the synthetic demo dataset.

From the repository root:

```bash
python scripts/00_prepare_reproducibility_inputs.py \
    --s3 /path/to/Supplementary_Table_S3_ATLAS_Summary.xlsx \
    --s4 /path/to/Supplementary_Table_S4_Ecological_Annotation.xlsx \
    --outdir data/generated

export SINGLEM_METAPACKAGE_PATH=/path/to/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb
export THREADS=8
bash scripts/singlem/run_singlem.sh /path/to/reads /path/to/assemblies /path/to/MAGs /path/to/singlem all

python scripts/singlem/analyze_singlem.py --out /path/to/singlem --mag-inputs data/generated
Rscript scripts/singlem/singlem_statistics.R /path/to/singlem
```

The raw-read wrapper is reconstructed from the documented input/output layout; check its flags against the original raw-read logs before using it for a release. The assembly, MAG and appraisal commands use the options exercised on the workstation. The shell script supports stages `reads`, `references`, `appraisal` and `all`. Completed jobs use sentinel files; remove the relevant sentinel if intentionally changing inputs or parameters. Use a new output directory for a clean reproducibility run. If the existing workstation outputs are available, skip upstream processing and run the Python and R commands against those tables (prefer a copy for parity checks). Python stages include recovery, alpha, rarefaction, beta, compare and diagnose. The last two require `--mag-inputs`; diagnose also requires targeted appraisal with found_in.

## Analytical definitions

**Appraisal.** Profile all raw metagenomes and obtain marker windows independently from each assembly and all final MAGs. Appraise against the MAG catalogue alone, pooled assemblies plus MAGs, and each sample's own assembly separately. Exact and imperfect modes are retained; the latter uses a requested identity of 0.9666666666666667 (maximum divergence two bases for a 60-base window). This is marker-window matching, not measured genome ANI.

**Off-target filtering.** Use SingleM summarise --exclude-off-target-hits before any custom calculation. Limit custom analyses to Bacteria and Archaea. Join on sample, marker gene and nucleotide window; require the three recovery classes to be mutually exclusive and exhaustive, and require the binned subset to be included in assembled markers for these data. A nonzero discrepancy is an error, not evidence of unassembled populations. Do not use the earlier unfiltered analysis for reporting.

**Recovery estimator.** For each sample, domain and marker, divide OTU coverage by the total marker coverage. Average category fractions across markers with equal marker weights. MAG-linked, assembled-only and unaccounted are mutually exclusive. Report this estimator separately from SingleM native appraise estimates, which have a different aggregation procedure. These are marker abundance fractions, not fractions of total metagenomic DNA reads. Pooled assembly matches do not imply recovery in the same sample. Domain estimates from condensed profiles sum assigned coverage within each domain and normalize across prokaryotic domains.

**Taxonomic representation.** Aggregate normalized marker fractions by rank and recovery category, retaining unresolved lineages; divide by the number of markers and compute within-lineage MAG representation. For the three-profile comparison, restrict original S3 MAG abundances to Bacteria and renormalize all profiles within their respective fractions. Retain the original GTDB names, including unresolved ancestors. Label matching across releases is not a validated taxonomic crosswalk; fine-rank disagreement can reflect marker resolution or taxonomy changes. The taxon labels are not a metabolic classifier.

**Alpha diversity.** Calculate Shannon entropy (natural logarithm), richness and effective Shannon diversity separately for marker sequences. Unrarefied coverage-weighted summaries and minimum-hit sensitivity are retained. Rarefy integer num_hits without replacement with 100 multivariate hypergeometric draws per sample, marker and depth (seed 42). Primary depth 1,000 hits; require at least 30 common bacterial markers (fallback 500), with all 35 retained in the actual dataset. Use the same retained marker set across samples and depths. Average each marker's Shannon values over draws, then take the median across markers. Primary analysis retains all observed variants; thresholds 2 and 5 after subsampling are sensitivity analyses, not sequencing-error correction. The 2,000-hit analysis is a depth sensitivity test. These are marker-window diversity estimates, not species richness. Inspect rarefaction curves before claiming saturation.

**Beta diversity.** For each bacterial marker, normalize coverage within sample, calculate Bray-Curtis distances, then average the 35 matrices. Minimum-hit thresholds 1 (primary), 2 and 5. Use square-root distances for PCoA, PERMANOVA and PERMDISP. Include all 19 samples in ordination; restrict inference to the four replicated categories (17 samples) and the carbonate-plastic contrast (13 samples). Use 9,999 free permutations, median-based bias-adjusted dispersion, seed 42. Alpha tests use Kruskal-Wallis and exact two-sided Wilcoxon rank-sum tests. The carbonate-plastic P values are unadjusted planned contrasts; sensitivity analyses are not independent confirmations. Free permutations assume exchangeable samples; check cruise/site/deployment dependence before final inference. Nonsignificant dispersion tests do not establish homogeneity.

**Targeted diagnostic.** Use exact/imperfect appraise with found_in for EGGPM35, EGGPM36 and PM29. Group reported MAG identifiers per marker window, and retain both exclusively reported and possible matches to MAG124 or UWMA-0217. Normalize either against all bacterial marker abundance or within the MAG-linked fraction. These are conditional on reported hits; they are not rigorous population abundance bounds or a proof of genome purity.

## Outputs and validation

- `analysis_qc/on_target_otu.tsv`: filtered marker observations.
- `appraisal/{exact,imperfect}/`: native summaries and category OTU tables.
- `analysis_filtered/`: recovery, taxonomic and unrarefied alpha tables; compressed per-OTU recovery.
- `rarefaction/`: marker depth, selection, per-marker curves, sample summaries.
- `beta_diversity/`: dissimilarities, ordination, statistical tests and R session information.
- `mag_comparison/`: three taxonomic profiles per rank, pairwise dissimilarities, unmatched labels and MAG124/UWMA-0217 diagnostic.

The recorded workstation run resolved all off-target discrepancies and retained 35 bacterial markers in every sample and 37 archaeal marker families overall. A clean reproduction should have no overlapping or missing appraisal categories. Compare sample-level recovery, rarefaction and profile dissimilarity tables against the original outputs and record differences. The original workstation analyses ran successfully; the consolidated repository versions have not yet been rerun end-to-end on those inputs. The existing repository smoke test is not a biological validation of SingleM.

Large raw reads, metapackages, assemblies, OTU tables and manuscript-specific abundance data are intentionally not duplicated in GitHub. Archive the exact inputs, execution logs, environment, hashes and final source-data tables with the revision release.
