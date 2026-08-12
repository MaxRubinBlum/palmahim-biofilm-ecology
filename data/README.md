# Input data schemas

The repository intentionally avoids hard-coded project filenames. Supply input tables in CSV/TSV/XLSX form and adapt the short loading block at the top of each script.

## Abundance matrix

Rows: samples  
Columns: MAG identifiers  
Values: relative abundance (fractions or percentages, but use one convention consistently).

## Raw-abundance matrix

Rows: samples  
Columns: MAG identifiers  
Values: non-negative abundance/coverage values used for C1 edge weighting.

## Sample metadata

Required columns:

```text
sample
habitat
```

Optional columns can include replicate identifiers or other environmental metadata.

## MAG taxonomy

Required columns:

```text
MAG
Order
```

Additional GTDB ranks may be retained.

## Ecological trait table

Rows: MAGs. Binary/categorical columns may include:

```text
MAG
Ecological_guild
Methanotroph
Sulfur_oxidizing_autotroph
Autotrophic_carbon_fixation
CBB_I
CBB_II
rTCA
Methane_oxidation
Methanol_consumption
Formate_consumption
Formaldehyde_processing
Fermentation
Protein_degradation
Polysaccharide_degradation
```

The exact trait names may be changed, provided the mapping in the relevant script is updated.

## Metatranscriptome gene table

Required columns:

```text
gene_id
MAG
length_bp
read_count
```

## Heterotroph profile table

One row per MAG with taxonomy, sample-level abundance, MEROPS hits/families, CAZyme hits/families and curated transporter/fermentation features.

## Vitamin table

One row per MAG with complete-pathway calls for each vitamin and, for cobalamin, the number of supported M00122 steps (0–7).
