#!/bin/bash
# FIX: Added set -euo pipefail so failures inside the loop are not silently
# swallowed. Without this a failed MACSE run produces no output but the loop
# continues, giving the impression all genes were processed.
set -euo pipefail

BASE="/home/alasdair/3D_UHU_evo_analysis"
GENE_DIR="$BASE/extracted_genes"
ANC_DIR="$BASE/WT_gene_fasta"
OUT_DIR="$BASE/macse_output"

mkdir -p "$OUT_DIR"

for dir in "$GENE_DIR"/*; do
    [ -d "$dir" ] || continue

    gene_anc=$(basename "$dir")

    # FIX: Same underscore-safe split as build_gene_fastas.sh —
    # split on the LAST underscore so gene names with underscores work.
    anc="${gene_anc##*_}"
    gene="${gene_anc%_*}"

    anc_fasta="$ANC_DIR/${gene}_${anc}.fasta"
    all_fasta="$OUT_DIR/${gene}_${anc}_all.fasta"
    aln_nt="$OUT_DIR/${gene}_${anc}_aligned_nt.fasta"
    aln_aa="$OUT_DIR/${gene}_${anc}_aligned_aa.fasta"

    # Skip if ancestor FASTA is missing or empty
    if [ ! -s "$anc_fasta" ]; then
        echo "⚠ Missing or empty ancestor FASTA for ${gene}_${anc}, skipping"
        continue
    fi

    # Collect extracted CDS files for this gene
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
