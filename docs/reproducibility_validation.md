# Reproducibility validation against manuscript supplementary outputs

The committed machine-readable inputs were reconstructed directly from the final processed analysis files supplied with the manuscript.

## Taxonomic-functional beta-diversity

`data/mag_abundance.tsv` contains sample × MAG relative abundance as fractions (0–1), derived from the normalized abundance block of `summary_biofilm(20260814-104958).xlsx`. Each source sample summed to 100% before conversion to fractions. `data/sample_metadata.tsv` marks `AnemPM22` for exclusion from the 18-biofilm taxonomic-functional comparison.

Using `data/mag_abundance.tsv`, `data/strict_guild_membership.tsv` and `data/primary_production_membership.tsv`, the repository workflow reproduces Supplementary Table S5 exactly:

| Profile | Mean Bray-Curtis |
|---|---:|
| Taxonomic MAG composition | 0.6210749038914138 |
| Strict ecological guilds | 0.13066208127467652 |
| Strict primary-production traits | 0.069498470943168 |

The one-sided paired Wilcoxon signed-rank tests comparing taxonomic mean dissimilarity with each functional representation both give `P = 3.814697265625e-06`, matching Supplementary Table S5.

## Functional redundancy/core validation

The following values calculated from `data/mag_abundance.tsv` and `data/curated_mag_traits.tsv` match Supplementary Table S5 exactly:

| Function | Carrier MAGs | Min abundance | Median abundance | Max abundance | Samples >=1% | Median effective MAG number |
|---|---:|---:|---:|---:|---:|---:|
| Methanotroph | 56 | 0.10557377850262348 | 0.34186264956462364 | 0.6365795298467646 | 18 | 5.182640584838683 |
| Sulfur-oxidizing autotroph | 68 | 0.0380177505757026 | 0.09912832057477684 | 0.18342028806508898 | 18 | 5.25019571923488 |
| Autotrophic carbon fixer | 133 | 0.05030632875825378 | 0.24354850495100763 | 0.35090368365433233 | 18 | 6.978115142936405 |
| CBB Form I | 95 | 0.04260993323682298 | 0.16682002091061932 | 0.23692374197071858 | 18 | 6.519453620189055 |
| CBB Form II | 27 | 0.009311855034296717 | 0.05605683354418451 | 0.1428109408151841 | 17 | 2.67307830145243 |
| rTCA | 19 | 0.0009895365670411023 | 0.005903263292698565 | 0.07223099382228745 | 6 | 3.575476265458539 |

## Provenance boundaries

The repository can reproduce the central custom abundance-weighted analyses from curated inputs. The conversion from raw annotation-tool output to final expert-curated ecological assignments remains partly rule-based/manual. `provenance/source_files_manifest.tsv` records the exact supplied origin files and SHA-256 checksums. Small inspectable provenance derivatives are committed; larger raw annotation and transcriptomic files should be archived with the publication release rather than duplicated in GitHub.

Historical featureCounts commands recovered from the original count-file headers are retained in `provenance/metatranscriptomics/featurecounts_commands.txt`. These commands should be treated as the historical source of truth for featureCounts invocation; generic shell examples in the repository are illustrative unless explicitly stated otherwise.
