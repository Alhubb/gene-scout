#!/bin/bash
set -euo pipefail

: "${BASE:?Need BASE}"
: "${GENOME_FASTA_DIR:?Need GENOME_FASTA_DIR}"
: "${WT_GENES:?Need WT_GENES}"

OUTBASE="${BASE}/minimap2_output"
mkdir -p "$OUTBASE"

for gene in "$WT_GENES"/*.fasta "$WT_GENES"/*.fna; do
    [ -f "$gene" ] || continue
    gene_base=$(basename "$gene")
    gene_base="${gene_base%.*}"

    outdir="${OUTBASE}/${gene_base}_minimap_hits"
    mkdir -p "$outdir"

    for g in "$GENOME_FASTA_DIR"/*.fna "$GENOME_FASTA_DIR"/*.fasta "$GENOME_FASTA_DIR"/*.fna.gz; do
        [ -f "$g" ] || continue
        asm=$(basename "$g")
        asm="${asm%.fna.gz}"
        asm="${asm%.fna}"
        asm="${asm%.fasta}"

        minimap2 -x asm5 -c "$g" "$gene" > "${outdir}/${asm}.paf"
    done
done
