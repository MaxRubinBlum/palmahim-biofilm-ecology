# Nature Research Code and Software Submission Checklist mapping

This repository is a collection of manuscript-specific analysis and figure-generation scripts rather than a standalone software package. The table below maps the Nature Research software checklist to concrete repository resources so that editors and reviewers can quickly verify each item.

| Nature checklist item | Repository resource |
|---|---|
| Compiled standalone software and/or source code | Source code in [`scripts/`](../scripts/) |
| Small simulated or real dataset to demo the code | Synthetic reviewer dataset in [`demo/`](../demo/) |
| All software dependencies and operating systems, including version numbers | [`README.md`](../README.md), `System requirements`; [`requirements.txt`](../requirements.txt); [`environment.yml`](../environment.yml) |
| Versions the software has been tested on | GitHub Actions workflow [`.github/workflows/reproducibility.yml`](../.github/workflows/reproducibility.yml): Ubuntu 24.04, Python 3.11/3.12/3.13, R 4.4 companion test |
| Required non-standard hardware | None for the reviewer demo or downstream custom analyses; see `System requirements` in the README |
| Installation instructions | `Installation` in [`README.md`](../README.md) |
| Typical install time | `Installation` in [`README.md`](../README.md) |
| Demo instructions | [`demo/README.md`](../demo/README.md) and `Reviewer quick start` in the main README |
| Expected demo output | [`demo/expected_metrics.json`](../demo/expected_metrics.json) and [`demo/README.md`](../demo/README.md) |
| Expected demo run time | `Reviewer quick start` / `Demo` in [`README.md`](../README.md); runtime is also printed by `scripts/run_demo.py` |
| Instructions for use on data | [`data/README.md`](../data/README.md), [`docs/analysis_manifest.md`](analysis_manifest.md), and command examples in the main README |
| Reproduction instructions | [`docs/analysis_manifest.md`](analysis_manifest.md) and [`docs/reproducibility_validation.md`](reproducibility_validation.md) |
| Software license | [`LICENSE`](../LICENSE) - MIT License |
| Open-source repository link | https://github.com/MaxRubinBlum/palmahim-biofilm-ecology |
| Detailed description of code functionality | Manuscript Supplementary Methods; repository-side mapping in [`docs/analysis_manifest.md`](analysis_manifest.md) |

## Scope and provenance boundary

The public repository captures the **custom manuscript-specific computational steps** connecting curated genome annotations and abundance profiles to statistical results, tables and figures. Established third-party bioinformatic tools are not vendored. Their historical versions and parameters are reported in the manuscript Methods and Supplementary Table S2.

The central abundance-weighted analyses can be reconstructed from Supplementary Tables S3 and S4 using `scripts/00_prepare_reproducibility_inputs.py`. The small simulated dataset under `demo/` is provided specifically so reviewers can test the code without requiring access to those manuscript data files.

Large upstream annotation and transcriptomic intermediates are represented by provenance records and checksums rather than duplicated in this repository. This distinction is intentional and is described in the README and provenance documentation.
