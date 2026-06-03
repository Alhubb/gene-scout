#!/bin/bash
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================

BASE="/home/alasdair/3D_UHU_evo_analysis"
SCRIPTS="$BASE/scripts"
MINIMAP_OUT="/media/alasdair/External/3D_UHU/minimap2_output"

EXTRACTED="$BASE/extracted_genes"
WT="$BASE/WT_gene_fasta"
PER_GENE_FASTAS="$BASE/per_gene_fastas"
PER_GENE_FASTAS_DEDUP="$BASE/per_gene_fastas_dedup"
MACSE_OUT="$BASE/macse_output"
STOP_CALLS="$BASE/stop_calls"

MACSE_SUMMARY="$BASE/mutation_cds_macse_strict_summary.tsv"
FRAMESHIFTS="$BASE/frameshifts.tsv"
INDELS="$BASE/indels.tsv"
STOP_GENES="$BASE/stop_genes.tsv"

# MUTATION_LIST is optional — if not set, screening step is skipped.
# Set via: export MUTATION_LIST=/path/to/mutation_list.tsv
# or pass --mutations to gene-scout.
MUTATION_LIST="${MUTATION_LIST:-}"

# ============================================================
# HELPERS
# ============================================================

step() { echo; echo "── Step $1: $2 ──────────────────────────"; }
done_() { echo "   ✅ $1"; }
warn() { echo "   ⚠  $1"; }

# ============================================================
# 1. BUILD PER-GENE FASTAs
# ============================================================

step "1" "Building per-gene FASTAs"
mkdir -p "$PER_GENE_FASTAS"

built=0; skipped=0
for dir in "$EXTRACTED"/*; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"
    lineage="${core##*_}"
    gene="${core%_*}"

    wt_fasta="$WT/${gene}_${lineage}.fasta"
    out_fasta="$PER_GENE_FASTAS/${gene}_${lineage}.fasta"

    if [ ! -f "$wt_fasta" ]; then
        warn "Missing WT for ${gene}_${lineage} — skipping"
        skipped=$((skipped + 1)); continue
    fi

    shopt -s nullglob
    files=( "$dir"/*.bakta."$gene".fasta )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        warn "No mutant CDS files in $dir — skipping"
        skipped=$((skipped + 1)); continue
    fi

    { cat "$wt_fasta"; echo
      for f in "${files[@]}"; do cat "$f"; echo; done
    } > "$out_fasta"
    built=$((built + 1))
done

done_ "$built gene/lineage FASTAs built${skipped:+, $skipped skipped}"

# ============================================================
# 1.5 DEDUPLICATE PER-GENE FASTAs
# ============================================================

step "1.5" "Deduplicating per-gene FASTAs"
mkdir -p "$PER_GENE_FASTAS_DEDUP"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f")
    out="$PER_GENE_FASTAS_DEDUP/$base"
    conda run -n python seqkit rmdup -s "$f" > "$out"
    n=$(grep -c "^>" "$out")
    echo "   ${base%.fasta}: $n unique CDS variants"
done

done_ "Deduplication complete"

# ============================================================
# 2. MACSE CODON-AWARE ALIGNMENT
# ============================================================

step "2" "Running MACSE alignment (4 cores)"
mkdir -p "$MACSE_OUT"

find "$PER_GENE_FASTAS_DEDUP" -maxdepth 1 -type f -name "*.fasta" -print0 | \
xargs -0 -r -n 1 -P 4 bash -c '
    f="$1"
    base=$(basename "$f" .fasta)
    if [ -f "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" ]; then
        echo "   ${base}: already aligned — skipping"
        exit 0
    fi
    echo "   Aligning: $base"
    conda run -n macse macse \
      -prog alignSequences \
      -seq "$f" \
      -out_AA "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" \
      -out_NT "'"$MACSE_OUT"'/${base}_aligned_nt.fasta"
' _

done_ "MACSE alignments complete"

# ============================================================
# 3. MACSE MUTATION SCAN
# ============================================================

step "3" "Scanning MACSE alignments for mutations"

conda run -n python python \
    "$SCRIPTS/scan_mutations_cds_macse.py" \
    "$MACSE_OUT" \
    > "$MACSE_SUMMARY"

done_ "MACSE summary → $(basename "$MACSE_SUMMARY")"

# ============================================================
# 4. FRAMESHIFT DETECTION (position-aware, MACSE NT alignments)
# ============================================================

step "4" "Detecting position-specific frameshifts"

if [ -z "${MUTATION_LIST}" ]; then
    warn "No MUTATION_LIST — skipping frameshift detection"
    > "$FRAMESHIFTS"
else
    conda run -n python python \
        "$SCRIPTS/detect_frameshifts_macse.py" \
        "$MACSE_OUT" \
        "$MUTATION_LIST" \
        "$FRAMESHIFTS"

    n=$(wc -l < "$FRAMESHIFTS")
    done_ "$n position-specific frameshift hits → $(basename "$FRAMESHIFTS")"
fi

# ============================================================
# 4.5 IN-FRAME INDEL DETECTION (three-tier, MACSE AA alignments)
# ============================================================

step "4.5" "Detecting in-frame indels (delins / inframe_deletion)"

if [ -z "${MUTATION_LIST}" ]; then
    warn "No MUTATION_LIST — skipping indel detection"
    > "$INDELS"
else
    conda run -n python python \
        "$SCRIPTS/detect_indels_macse.py" \
        "$MACSE_OUT" \
        "$MUTATION_LIST" \
        "$INDELS"

    n=$(wc -l < "$INDELS")
    done_ "$n indel hits → $(basename "$INDELS")"
fi

# ============================================================
# 5. STOP-GAIN / STOP-LOSS DETECTION
# ============================================================

step "5" "Detecting stop-gain / stop-loss mutations"
mkdir -p "$STOP_CALLS"
> "$STOP_CALLS/stop_calls.tsv"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    conda run -n python python \
        "$SCRIPTS/detect_stop_mutations.py" "$f" \
        >> "$STOP_CALLS/stop_calls.tsv"
done

awk 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' "$STOP_CALLS/stop_calls.tsv" \
    | sort -u > "$STOP_GENES"

n=$(wc -l < "$STOP_GENES")
done_ "$n stop events → $(basename "$STOP_GENES")"

# ============================================================
# 6. MUTATION SCREENING
# ============================================================

step "6" "Screening predefined mutations"

conda run -n python python \
    "$SCRIPTS/screen_mutations.py" \
    --mutations   "$MUTATION_LIST" \
    --macse       "$MACSE_SUMMARY" \
    --frameshifts "$FRAMESHIFTS" \
    --indels      "$INDELS" \
    --stops       "$STOP_GENES" \
    --output      "$BASE/mutation_screen_results.tsv"

done_ "Results → mutation_screen_results.tsv"

# ============================================================

echo
echo "══════════════════════════════════════════"
echo "  ✅  PIPELINE COMPLETE"
echo "  ➡   $BASE/mutation_screen_results.tsv"
echo "══════════════════════════════════════════"
echo
