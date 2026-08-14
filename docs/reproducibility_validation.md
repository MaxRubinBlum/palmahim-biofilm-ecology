# Reproducibility validation against manuscript supplementary outputs

The central custom analysis inputs were independently reconstructed from the final processed manuscript data and checked against Supplementary Table S5. For the public repository, the same inputs are generated reproducibly from the published Supplementary Tables S3 and S4 using `scripts/00_prepare_reproducibility_inputs.py`; duplicate wide TSV matrices are therefore not stored in version control.

## Source-data validation

The final normalized abundance block in the supplied `summary_biofilm(20260814-104958).xlsx` contained 19 sample columns. Each column summed to 100% before conversion to fractions. The corresponding published Supplementary Table S3 exposes these columns explicitly as `<sample> relative abundance (%)`, and the input builder verifies the same 100% invariant before conversion.

Supplementary Table S4 contains 565 curated MAG rows. The input builder checks that every S3 MAG is represented in S4 and reconstructs the four mutually exclusive strict ecological guilds and three primary-production traits used in the manuscript analyses. `AnemPM22` is retained in the source tables but marked for exclusion from the 18-biofilm taxonomic-functional comparison.

## Taxonomic-functional beta-diversity

Using the S3 abundance values and S4 curated assignments, the reconstructed analysis inputs reproduce Supplementary Table S5 exactly:

| Profile | Mean Bray-Curtis |
|---|---:|
| Taxonomic MAG composition | 0.6210749038914138 |
| Strict ecological guilds | 0.13066208127467652 |
| Strict primary-production traits | 0.069498470943168 |

The one-sided paired Wilcoxon signed-rank tests comparing each biofilm's mean taxonomic dissimilarity with its corresponding functional mean give `P = 3.814697265625e-06` for both comparisons, matching Supplementary Table S5 exactly.

The workflow uses all 18 biofilms for the distance summaries and paired tests. The singleton crab-carapace sample remains in these analyses; the companion `vegan` script separately restricts substrate PERMANOVA to replicated substrate classes.

## Functional redundancy/core validation

The following values calculated from the reconstructed abundance and curated trait matrices match Supplementary Table S5 exactly:

| Function | Carrier MAGs | Min abundance | Median abundance | Max abundance | Samples >=1% | Median effective MAG number |
|---|---:|---:|---:|---:|---:|---:|
| Methanotroph | 56 | 0.10557377850262348 | 0.34186264956462364 | 0.6365795298467646 | 18 | 5.182640584838683 |
| Sulfur-oxidizing autotroph | 68 | 0.0380177505757026 | 0.09912832057477684 | 0.18342028806508898 | 18 | 5.25019571923488 |
| Autotrophic carbon fixer | 133 | 0.05030632875825378 | 0.24354850495100763 | 0.35090368365433233 | 18 | 6.978115142936405 |
| CBB Form I | 95 | 0.04260993323682298 | 0.16682002091061932 | 0.23692374197071858 | 18 | 6.519453620189055 |
| CBB Form II | 27 | 0.009311855034296717 | 0.05605683354418451 | 0.1428109408151841 | 17 | 2.67307830145243 |
| rTCA | 19 | 0.0009895365670411023 | 0.005903263292698565 | 0.07223099382228745 | 6 | 3.575476265458539 |

For the published S5 redundancy panel, run `03_functional_redundancy.py` with the explicit trait columns `Methanotroph Sulfur_oxidizing_autotroph Autotrophic_carbon_fixation CBB_I CBB_II rTCA`.

## Provenance boundaries

The repository can reproduce the central custom abundance-weighted analyses from the published curated inputs. The conversion from third-party annotation outputs to the final expert-curated ecological assignments in S4 remains partly rule-based/manual. This boundary is stated explicitly rather than being represented as a fully automated classification step.

`provenance/source_files_manifest.tsv` records the exact supplied origin filenames, sizes and SHA-256 checksums. Small inspectable provenance derivatives are committed; larger raw annotation and transcriptomic files should accompany the archival publication release rather than be duplicated in GitHub.

Historical featureCounts commands recovered from the original count-file headers are retained in `provenance/metatranscriptomics/featurecounts_commands.txt`. These commands are the historical source of truth for featureCounts invocation; the generic shell workflow is illustrative because the original libraries used different command-line flags.
