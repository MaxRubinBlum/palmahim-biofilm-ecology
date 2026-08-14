# Input data schemas

## Canonical reconstruction from the published supplements

For the central abundance-weighted analyses, use `scripts/00_prepare_reproducibility_inputs.py` with the published Supplementary Tables S3 and S4. The script reconstructs the exact machine-readable inputs used by the custom analysis code and verifies basic invariants before writing them.

```bash
python scripts/00_prepare_reproducibility_inputs.py \
  --s3 Supplementary_Table_S3_ATLAS_Summary.xlsx \
  --s4 Supplementary_Table_S4_Ecological_Annotation.xlsx \
  --outdir data/generated
```

The generated files are deliberately not duplicated in Git because S3 and S4 are the canonical published data source.

## MAG abundance

`mag_abundance.tsv`: rows = samples; columns = MAG identifiers; values = relative-abundance **fractions (0–1)**. The input builder reads the explicit relative-abundance (%) columns from S3, verifies that each sample sums to 100%, and divides by 100.

## Sample metadata

`sample_metadata.tsv` contains:

```text
sample
habitat
include_taxonomic_functional_beta_diversity
```

The inclusion flag is 0 for `AnemPM22` and 1 for the 18 biofilms used for Supplementary Fig. 3/Table S5. The singleton crab-carapace biofilm remains part of the 18-sample distance/Wilcoxon analysis; the companion PERMANOVA script independently restricts substrate tests to replicated classes.

## MAG taxonomy

`mag_taxonomy.tsv`: index/first column = MAG. Required downstream field:

```text
Order
```

Other GTDB ranks are retained by the input builder.

## Curated ecological traits

`curated_mag_traits.tsv` is reconstructed from Supplementary Table S4 and retains the final curated MAG-level assignments. These are the authoritative trait calls for downstream custom analyses. The raw-annotation-to-curated-trait step includes the operational rules in Supplementary Table S1 and expert review where applicable; downstream plotting/statistical scripts should not reinterpret the raw annotation files independently.

## Strict ecological-guild membership

Rows = MAGs; mutually exclusive binary columns:

```text
methanotroph
sulfur_oxidizing_autotroph
other_autotroph
other_community_member
```

The input builder checks that each MAG belongs to exactly one strict ecological guild.

## Primary-production membership

Rows = MAGs; binary columns:

```text
methanotrophy
sulfur_oxidizing_autotrophy
autotrophic_carbon_fixation
```

These traits are not mutually exclusive. Functional beta-diversity remains abundance weighted: sample MAG abundance is multiplied by these membership matrices.

## Raw abundance / coverage

Rows = samples; columns = MAG identifiers; non-negative abundance/coverage values used in the C1 interaction score. These values are distinct from normalized relative abundance.

## C1 traits

Rows = MAGs; binary columns:

```text
methane_oxidizer
methanol_consumer
formate_consumer
formaldehyde_processing
```

Use the curated C1 classifications described in the Supplementary Methods rather than deriving these labels ad hoc from individual annotation hits.

## Fig. 3 MAG trait table

One row per MAG with:

```text
MAG
Order
abundance
<binary curated trait columns>
```

## Metatranscriptome gene table

```text
gene_id
MAG
length_bp
read_count
```

Optional functional annotation columns can be retained. Historical featureCounts invocations recovered from original count-file headers are in `provenance/metatranscriptomics/featurecounts_commands.txt`.

## Heterotroph profile / Supplementary Table S6

One row per MAG containing `MAG`, `Order`, sample abundance columns, transporter/fermentation fields, and MEROPS/CAZyme hit/family counts. The final integrated S6 abundance fields are expressed in **percentage points**, not 0–1 fractions. Accordingly, `07_heterotroph_summaries.py` defaults to `--abundance-units percent`; use `--abundance-units fraction` only when supplying a 0–1 matrix.

Curated category columns may begin with `MEROPS_` or `CAZyme_` and are averaged per order by the summary script.

## Vitamin traits

Rows = MAGs; required `Order`, complete-pathway 0/1 columns supplied to `--vitamin-columns`; optional:

```text
B12_steps
functional_category
```

`functional_category` should already be mutually exclusive/hierarchical to prevent double counting. `08_vitamin_summaries.py` additionally produces abundance-weighted habitat distributions of complete (7/7), near-complete (6/7), partial (1–5/7), and absent (0/7) B12 pathway status.

## CTD timeline

```text
timestamp
dissolved_oxygen
```

Optional `event` column marks synchronized ROV operations.

## Upstream provenance

`provenance/source_files_manifest.tsv` records the supplied origin filenames, sizes, SHA-256 checksums and repository/archive policy. Large third-party annotation and transcriptomic outputs should accompany the archival publication release rather than be duplicated in GitHub.
