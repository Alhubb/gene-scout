#!/bin/bash
set -euo pipefail

: "${BASE:?Need BASE}"
: "${GENOME_FASTA_DIR:?Need GENOME_FASTA_DIR}"

PAF_BASE="${BASE}/minimap2_output"
OUTBASE="${BASE}/extracted_genes"
SCRIPTS="$(dirname "$(realpath "$0")")"

mkdir -p "$OUTBASE"

for dir in "$PAF_BASE"/*_minimap_hits; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    gene="${core%%_*}"
    lineage="${core#*_}"

    outdir="${OUTBASE}/${gene}_${lineage}_minimap_hits"
    mkdir -p "$outdir"

    for paf in "$dir"/*.paf; do
        asm=$(basename "$paf" .paf)

        for ext in fna fasta fna.gz; do
            genome="${GENOME_FASTA_DIR}/${asm}.${ext}"
            [ -f "$genome" ] && break
        done

        [ -f "$genome" ] || continue

        out_fasta="${outdir}/${asm}.bakta.${gene}.fasta"

        zcat -f "$genome" | \
        python "$SCRIPTS/extract_gene_from_paf.py" \
            "$paf" /dev/stdin "$out_fasta"

        [ -s "$out_fasta" ] || rm -f "$out_fasta"
    done
done
