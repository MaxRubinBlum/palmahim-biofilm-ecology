# Functional-annotation provenance

The manuscript uses a curated MAG-level ecological annotation table as the authoritative input for downstream abundance-weighted analyses (`data/curated_mag_traits.tsv`).

The curated table was assembled from multiple annotation sources using the operational definitions in Supplementary Table S1. Relevant supplied origin files include METABOLIC results, additional functional HMM calls, QSAP signaling annotations, MacSyFinder system calls, MEROPS peptidase annotations and dbCAN/CAZyme annotations.

Small, inspectable provenance derivatives are retained here. Exact source filenames, sizes and SHA-256 checksums are recorded in `../source_files_manifest.tsv`. Large raw annotation outputs should be deposited with the archival publication release (for example Zenodo) and verified using those checksums rather than duplicated in GitHub.

The raw-annotation-to-curated-trait transition includes conservative operational rules and, where needed, expert review. The repository therefore distinguishes **annotation provenance** from **downstream computational reproducibility** rather than implying that every curated ecological assignment arose from a fully automated classifier.
