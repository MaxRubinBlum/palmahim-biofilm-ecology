# Metatranscriptomic provenance

This directory records historical provenance recovered from the supplied transcriptomic analysis outputs.

`featurecounts_commands.txt` contains the actual featureCounts command lines embedded in the headers of the supplied count files. These commands differ between libraries and are the historical source of truth for featureCounts invocation.

Large count/expression tables are tracked by filename and SHA-256 checksum in `../source_files_manifest.tsv` and should be archived with the publication release rather than duplicated in GitHub.

The repository mapping shell script is a reference workflow, not a claim that every historical library was processed with identical featureCounts flags.
