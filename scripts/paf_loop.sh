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

    # ✅ SAFE PARSING
    if [[ "$core" == *_* ]]; then
        gene="${core%_*}"
        lineage="${core##*_}"
    else
        gene="$core"
        lineage="unknown"
    fi

    outdir="${OUTBASE}/${gene}_${lineage}_minimap_hits"
    mkdir -p "$outdir"

    echo "=== Processing $name (gene=$gene lineage=$lineage) ==="

    paf_count=0
    for paf in "$dir"/*.paf; do
        [ -f "$paf" ] || continue
        paf_count=$((paf_count + 1))

        asm=$(basename "$paf" .paf)

        # asm may be the full genome filename including extension
        # (e.g. 34A.fasta from 34A.fasta.paf) or just the stem (e.g. 34A).
        # Try the direct filename first, then fall back to appending extensions.
        genome=""
        if [ -f "${GENOME_FASTA_DIR}/${asm}" ]; then
            genome="${GENOME_FASTA_DIR}/${asm}"
        else
            for ext in fna fasta fna.gz; do
                candidate="${GENOME_FASTA_DIR}/${asm}.${ext}"
                if [ -f "$candidate" ]; then
                    genome="$candidate"
                    break
                fi
            done
        fi

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
