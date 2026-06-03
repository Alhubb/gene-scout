#!/bin/bash
# FIX: Added set -euo pipefail throughout for consistent error handling.
set -euo pipefail

BASE="/home/alasdair/3D_UHU_evo_analysis"
EXTRACTED="$BASE/extracted_genes"
WT="$BASE/WT_gene_fasta"
OUT="$BASE/per_gene_fastas"

mkdir -p "$OUT"

for dir in "$EXTRACTED"/*; do
    [ -d "$dir" ] || continue

    # e.g. cyaA_34
    name=$(basename "$dir")

    # FIX: Original used cut -d'_' -f1 / -f2, which silently breaks for any
    # gene name containing an underscore. We now split on the LAST underscore
    # so that multi-part gene names (e.g. rpo_B) survive intact, and the
    # lineage suffix (always a plain number) is cleanly separated.
    lineage="${name##*_}"      # everything after the last underscore  → 34
    gene="${name%_*}"          # everything before the last underscore → cyaA

    wt_fasta="$WT/${gene}_${lineage}.fasta"
    out_fasta="$OUT/${gene}_${lineage}.fasta"

    if [ ! -f "$wt_fasta" ]; then
        echo "⚠ Missing WT file for $gene $lineage, skipping"
        continue
    fi

    # Collect matching extracted CDS files
    shopt -s nullglob
    files=( "$dir"/*.bakta."$gene".fasta )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "⚠ No mutant CDS files found in $dir, skipping"
        continue
    fi

    echo "✅ Building $out_fasta (${#files[@]} mutant(s))"

    {
        # WT ancestor first
        cat "$wt_fasta"
        echo

        # All mutant CDS for this gene
        for f in "${files[@]}"; do
            cat "$f"
            echo
        done
    } > "$out_fasta"
done

echo "✅ Per-gene FASTAs written to $OUT"
