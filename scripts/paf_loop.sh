#!/bin/bash

set -euo pipefail

: "${BASE:?Need BASE}"
: "${GENOME_FASTA_DIR:?Need GENOME_FASTA_DIR}"

PAF_BASE="${MINIMAP_OUT:-${BASE}/minimap2_output}"
OUTBASE="${BASE}/extracted_genes"
SCRIPTS="$(dirname "$(realpath "$0")")"

mkdir -p "$OUTBASE"

for dir in "$PAF_BASE"/*_minimap_hits; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    lineage="${core##*_}"   # everything after the last _  → e.g. 34
    gene="${core%_*}"       # everything before the last _ → e.g. cyaA

    outdir="${OUTBASE}/${gene}_${lineage}_minimap_hits"
    mkdir -p "$outdir"

    echo "=== Processing $name (gene=$gene lineage=$lineage) ==="

    paf_count=0
    for paf in "$dir"/*.paf; do
        [ -f "$paf" ] || continue
        paf_count=$((paf_count + 1))

        asm=$(basename "$paf" .paf)

        # Resolve assembly path — try extensions in order
        genome=""
        for ext in fna fasta fna.gz; do
            candidate="${GENOME_FASTA_DIR}/${asm}.${ext}"
            if [ -f "$candidate" ]; then
                genome="$candidate"
                break
            fi
        done

        if [ -z "$genome" ]; then
            echo "WARNING: Assembly not found for $asm, skipping"
            continue
        fi

        out_fasta="${outdir}/${asm}.bakta.${gene}.fasta"

        python "$SCRIPTS/extract_gene_from_paf.py" \
            "$paf" \
            "$genome" \
            "$out_fasta"

        if [ ! -s "$out_fasta" ]; then
            echo "  Skipping empty extraction for $asm"
            rm -f "$out_fasta"
        fi
    done

    if [ "$paf_count" -eq 0 ]; then
        echo "WARNING: No .paf files found in $dir"
    fi

done
