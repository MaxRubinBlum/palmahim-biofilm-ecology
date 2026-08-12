# Input data schemas

The repository avoids hard-coded project filenames. Use CSV/TSV files with the following logical structures.

## MAG abundance
Rows = samples; columns = MAG identifiers; values = relative abundance. Use fractions (0–1) for threshold-based analyses unless a script explicitly expects percentages.

## Raw abundance / coverage
Rows = samples; columns = MAG identifiers; non-negative abundance/coverage values used in the C1 interaction score.

## Sample metadata
```text
sample
habitat
```

## MAG taxonomy
Index/first column = MAG. Required field:
```text
Order
```
Other GTDB ranks can be retained.

## Strict ecological-guild membership
Rows = MAGs; binary columns such as:
```text
methanotroph
sulfur_oxidizing_autotroph
other_autotroph
other_community_member
```

## Primary-production membership
Rows = MAGs; binary columns:
```text
methanotrophy
sulfur_oxidizing_autotrophy
autotrophic_carbon_fixation
```

## C1 traits
Rows = MAGs; binary columns:
```text
methane_oxidizer
methanol_consumer
formate_consumer
formaldehyde_processing
```

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
Optional functional annotation columns can be retained.

## Heterotroph profile
One row per MAG containing `MAG`, `Order`, sample abundance columns, transporter/fermentation fields, and MEROPS/CAZyme hit/family counts. Curated category columns may begin with `MEROPS_` or `CAZyme_` and will be averaged per order by the summary script.

## Vitamin traits
Rows = MAGs; required `Order`, complete-pathway 0/1 columns supplied to `--vitamin-columns`; optional:
```text
B12_steps
functional_category
```
`functional_category` should already be mutually exclusive/hierarchical to prevent double counting.

## CTD timeline
```text
timestamp
dissolved_oxygen
```
Optional `event` column marks synchronized ROV operations.
