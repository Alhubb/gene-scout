#!/bin/bash
set -euo pipefail

BASE="/home/alasdair/3D_UHU_evo_analysis"
EXTRACTED="$BASE/extracted_genes"
WT="$BASE/WT_gene_fasta"
OUT="$BASE/per_gene_fastas"

mkdir -p "$OUT"

for dir in "$EXTRACTED"/*; do
    [ -d "$dir" ] || continue

    # Strip suffix correctly
    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    # ✅ robust parsing
    if [[ "$core" == *_* ]]; then
        gene="${core%_*}"
        lineage="${core##*_}"
    else
        gene="$core"
        lineage="unknown"
    fi

    wt_fasta="$WT/${gene}_${lineage}.fasta"

    # ✅ fallback logic
    if [ ! -f "$wt_fasta" ]; then
        fallback="$WT/${gene}.fasta"
        if [ -f "$fallback" ]; then
            wt_fasta="$fallback"
        else
            echo "⚠ Missing WT file for $gene ($lineage), skipping"
            continue
        fi
    fi

    out_fasta="$OUT/${gene}_${lineage}.fasta"

    shopt -s nullglob
    files=( "$dir"/*.bakta."$gene".fasta )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "⚠ No mutant CDS files found in $dir, skipping"
        continue
    fi

    echo "✅ Building $out_fasta (${#files[@]} mutant(s))"

    {
        cat "$wt_fasta"
        echo

        for f in "${files[@]}"; do
            cat "$f"
            echo
        done
    } > "$out_fasta"

done

echo "✅ Per-gene FASTAs written to $OUT"
