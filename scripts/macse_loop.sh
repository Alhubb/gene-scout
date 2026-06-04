#!/bin/bash
set -euo pipefail

BASE="/home/alasdair/3D_UHU_evo_analysis"
GENE_DIR="$BASE/extracted_genes"
ANC_DIR="$BASE/WT_gene_fasta"
OUT_DIR="$BASE/macse_output"

mkdir -p "$OUT_DIR"

for dir in "$GENE_DIR"/*; do
    [ -d "$dir" ] || continue

    # ✅ Safe parsing
    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    if [[ "$core" == *_* ]]; then
        gene="${core%_*}"
        anc="${core##*_}"
    else
        gene="$core"
        anc="unknown"
    fi

    anc_fasta="$ANC_DIR/${gene}_${anc}.fasta"

    # ✅ WT fallback
    if [ ! -s "$anc_fasta" ]; then
        fallback="$ANC_DIR/${gene}.fasta"
        if [ -s "$fallback" ]; then
            anc_fasta="$fallback"
        else
            echo "⚠ Missing or empty ancestor FASTA for ${gene} (${anc}), skipping"
            continue
        fi
    fi

    all_fasta="$OUT_DIR/${gene}_${anc}_all.fasta"
    aln_nt="$OUT_DIR/${gene}_${anc}_aligned_nt.fasta"
    aln_aa="$OUT_DIR/${gene}_${anc}_aligned_aa.fasta"

    # ✅ Skip if already done
    if [ -f "$aln_aa" ]; then
        echo "   ⏭ Already aligned: ${gene}_${anc}"
        continue
    fi

    shopt -s nullglob
    extracted=( "$dir"/*.bakta."$gene".fasta )
    shopt -u nullglob

    if [ ${#extracted[@]} -eq 0 ]; then
        echo "⚠ No extracted CDS files for ${gene}_${anc}, skipping"
        continue
    fi

    {
        cat "$anc_fasta"
        echo
        for f in "${extracted[@]}"; do
            cat "$f"
            echo
        done
    } > "$all_fasta"

    macse \
        -prog alignSequences \
        -seq "$all_fasta" \
        -out_NT "$aln_nt" \
        -out_AA "$aln_aa"

    echo "✅ MACSE finished: ${gene}_${anc}"
done
``
