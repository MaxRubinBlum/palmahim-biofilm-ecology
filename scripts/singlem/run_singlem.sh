#!/usr/bin/env bash
# Reconstructed reference workflow; raw pipe invocation should be checked against archived logs.
set -euo pipefail
: "${SINGLEM_METAPACKAGE_PATH:?Set path to S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb}"
: "${1:?Usage: run_singlem.sh READS_DIR ASSEMBLY_DIR MAG_DIR OUTPUT_DIR [stage]}"
READS=$1
ASM=$2
MAG=$3
OUT=$4
STAGE=${5:-all}
THREADS=${THREADS:-8}
samples=(AnemPM22 Bio1 Bio2 Bio3 Bio4 Bio6 Bio7 Bio8 Bio9 CrabPM21 EGGPM35 EGGPM36 PBB PBT PM27 PM28 PM29 PM30 PM39)
mkdir -p "$OUT"/{otu,archive,profiles,assemblies,genomes,logs,appraisal,analysis_qc}
singlem --version | tee "$OUT/logs/singlem_version.txt"
run() {
    local id=$1
    shift
    if [[ ! -f "$OUT/logs/$id.done" ]]; then
        printf '%q ' "$@" > "$OUT/logs/$id.command"
        printf '\n' >> "$OUT/logs/$id.command"
        "$@" > "$OUT/logs/$id.log" 2>&1
        touch "$OUT/logs/$id.done"
    fi
}
if [[ "$STAGE" == all || "$STAGE" == reads ]]; then
    for s in "${samples[@]}"; do
        run "reads_$s" singlem pipe --threads "$THREADS" \
            -1 "$READS/${s}_1.fastq.gz" -2 "$READS/${s}_2.fastq.gz" \
            --metapackage "$SINGLEM_METAPACKAGE_PATH" \
            --otu-table "$OUT/otu/${s}_otu.tsv" \
            --archive-otu-table "$OUT/archive/${s}_archive.json" \
            --taxonomic-profile "$OUT/profiles/${s}_profile.tsv"
    done
fi
if [[ "$STAGE" == all || "$STAGE" == references ]]; then
    genomes=("$MAG"/MAG*.fasta)
    [[ "${#genomes[@]}" == 565 ]]
    run MAG_catalogue singlem pipe --genome-fasta-files "${genomes[@]}" \
        --metapackage "$SINGLEM_METAPACKAGE_PATH" --threads "$THREADS" \
        --otu-table "$OUT/genomes/MAG_catalogue_otu.tsv"
    for s in "${samples[@]}"; do
        run "assembly_$s" singlem pipe --genome-fasta-files "$ASM/$s.fasta" \
            --metapackage "$SINGLEM_METAPACKAGE_PATH" --threads "$THREADS" \
            --otu-table "$OUT/assemblies/${s}_otu.tsv"
    done
fi
if [[ "$STAGE" == all || "$STAGE" == appraisal ]]; then
    reads=(); assemblies=()
    for s in "${samples[@]}"; do
        reads+=("$OUT/otu/${s}_otu.tsv")
        assemblies+=("$OUT/assemblies/${s}_otu.tsv")
    done
    genome="$OUT/genomes/MAG_catalogue_otu.tsv"
    for mode in exact imperfect; do
        args=()
        if [[ "$mode" == imperfect ]]; then
            args=(--imperfect --sequence-identity 0.9666666666666667)
        fi
        for kind in MAGs pooled; do
            d="$OUT/appraisal/$mode/$kind"
            mkdir -p "$d"
            if [[ -f "$d/.done" ]]; then continue; fi
            extra=()
            if [[ "$kind" == pooled ]]; then
                extra=(--assembly-otu-tables "${assemblies[@]}" \
                    --output-assembled-otu-table "$d/assembled.tsv" \
                    --output-unbinned-otu-table "$d/unbinned.tsv")
            fi
            singlem appraise --metapackage "$SINGLEM_METAPACKAGE_PATH" \
                --metagenome-otu-tables "${reads[@]}" --genome-otu-tables "$genome" \
                "${args[@]}" "${extra[@]}" \
                --output-binned-otu-table "$d/binned.tsv" \
                --output-unaccounted-for-otu-table "$d/unaccounted.tsv" \
                > "$d/summary.tsv" 2> "$d/appraise.log"
            touch "$d/.done"
        done
        for s in "${samples[@]}"; do
            d="$OUT/appraisal/$mode/individual/$s"
            mkdir -p "$d"
            if [[ -f "$d/.done" ]]; then continue; fi
            singlem appraise --metapackage "$SINGLEM_METAPACKAGE_PATH" \
                --metagenome-otu-tables "$OUT/otu/${s}_otu.tsv" \
                --assembly-otu-tables "$OUT/assemblies/${s}_otu.tsv" \
                "${args[@]}" --output-assembled-otu-table "$d/assembled.tsv" \
                --output-unaccounted-for-otu-table "$d/unaccounted.tsv" \
                > "$d/summary.tsv" 2> "$d/appraise.log"
            touch "$d/.done"
        done
        d="$OUT/mag_comparison/targeted_appraisal"
        mkdir -p "$d"
        if [[ ! -f "$d/$mode.done" ]]; then
            singlem appraise --metapackage "$SINGLEM_METAPACKAGE_PATH" \
                --metagenome-otu-tables "$OUT/otu/EGGPM35_otu.tsv" \
                    "$OUT/otu/EGGPM36_otu.tsv" "$OUT/otu/PM29_otu.tsv" \
                --genome-otu-tables "$genome" "${args[@]}" --output-found-in \
                --output-binned-otu-table "$d/${mode}_binned.tsv" \
                > "$d/${mode}_summary.tsv" 2> "$d/$mode.log"
            touch "$d/$mode.done"
        fi
    done
    singlem summarise --input-otu-tables "${reads[@]}" \
        --exclude-off-target-hits --metapackage "$SINGLEM_METAPACKAGE_PATH" \
        --output-otu-table "$OUT/analysis_qc/on_target_otu.tsv"
fi
