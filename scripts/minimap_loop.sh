#!/bin/bash
# minimap_loop.sh
# Aligns each WT gene CDS against all genome assemblies, producing PAF files.
#
# Strategy: pre-build a minimap2 index (.mmi) for each genome once, then map
# all genes against each pre-built index in parallel.
#
# Why genome = reference:
#   PAF target columns (tname, tstart, tend, strand) always refer to the
#   reference sequence. extract_gene_from_paf.py uses these coordinates to
#   slice the CDS out of the genome assembly, so the genome must be the
#   reference and the gene must be the query.
#
# Required environment variables:
#   BASE             — project root directory
#   GENOME_FASTA_DIR — directory containing genome FASTA files (.fna / .fasta / .fna.gz)
#   WT_GENES         — directory containing WT gene CDS FASTA files
#
# Optional environment variables:
#   THREADS          — number of parallel mapping jobs (default: 4)
#   PRESET           — minimap2 preset (default: asm20)
#                      asm20 handles up to ~5% divergence; use a looser preset
#                      for more divergent datasets.
#
# Usage:
#   export BASE=... GENOME_FASTA_DIR=... WT_GENES=...
#   bash minimap_loop.sh

set -euo pipefail

: "${BASE:?Need BASE}"
: "${GENOME_FASTA_DIR:?Need GENOME_FASTA_DIR}"
: "${WT_GENES:?Need WT_GENES}"

OUTBASE="${BASE}/minimap2_output"
INDEX_DIR="${BASE}/minimap2_indices"
THREADS="${THREADS:-4}"
PRESET="${PRESET:-asm20}"

# ============================================================
# PHASE 1 — Build genome indices (once per genome, serial)
# ============================================================
# Indexing is I/O bound so is kept serial to avoid saturating the disk.
# Already-indexed genomes are skipped, making re-runs safe.

mkdir -p "$INDEX_DIR" "$OUTBASE"

echo "==============================="
echo "Phase 1: Indexing genomes"
echo "==============================="

for g in "$GENOME_FASTA_DIR"/*.fna \
          "$GENOME_FASTA_DIR"/*.fasta \
          "$GENOME_FASTA_DIR"/*.fna.gz; do
    [ -f "$g" ] || continue

    asm=$(basename "$g")
    asm="${asm%.fna.gz}"
    asm="${asm%.fna}"
    asm="${asm%.fasta}"

    idx="${INDEX_DIR}/${asm}.mmi"

    if [ -f "$idx" ]; then
        echo "  Skipping (already indexed): $asm"
        continue
    fi

    echo "  Indexing: $asm"
    minimap2 -x "$PRESET" -d "$idx" "$g"
done

echo "✅ Indexing complete"
echo

# ============================================================
# PHASE 2 — Map each gene against every pre-built index
#           Parallelised across THREADS cores with xargs -P
# ============================================================
# Each gene is processed in turn (outer loop).
# Within each gene, all genome mappings run in parallel — each
# minimap2 call is independent so there are no race conditions.
# Already-completed PAF files are skipped, making re-runs safe.

echo "==============================="
echo "Phase 2: Mapping genes ($THREADS parallel jobs)"
echo "==============================="

for gene in "$WT_GENES"/*.fasta "$WT_GENES"/*.fna; do
    [ -f "$gene" ] || continue

gene_file=$(basename "$gene")
gene_name="${gene_file%.*}"

# Optional: fallback if WT filename still contains lineage
if [[ "$gene_name" == *_* ]]; then
    gene="${gene_name%_*}"
    lineage="${gene_name##*_}"
else
    gene="$gene_name"
    lineage="unknown"
fi

outdir="${OUTBASE}/${gene}_${lineage}_minimap_hits"

    outdir="${OUTBASE}/${gene_base}_minimap_hits"
    mkdir -p "$outdir"

    echo "Processing gene: $gene_base"

    # Write index paths that still need mapping to a temp file.
    # Using a temp file avoids subshell variable-scoping issues
    # under set -euo pipefail when using process substitution.
    tmplist=$(mktemp)
    for idx in "$INDEX_DIR"/*.mmi; do
        [ -f "$idx" ] || continue
        asm=$(basename "$idx" .mmi)
        paf="${outdir}/${asm}.paf"
        [ -f "$paf" ] || echo "$idx"
    done > "$tmplist"

    total=$(wc -l < "$tmplist")
    echo "  $total mappings to run"

    # Run up to THREADS minimap2 processes simultaneously.
    xargs -a "$tmplist" -P "$THREADS" -I{} bash -c '
        idx="{}"
        asm=$(basename "$idx" .mmi)
        paf="'"$outdir"'/${asm}.paf"
        minimap2 -x "'"$PRESET"'" -c "$idx" "'"$gene"'" > "$paf"
    '

    rm -f "$tmplist"
    echo "✅ Finished gene: $gene_base"
    echo
done

echo "✅ Mapping complete"
