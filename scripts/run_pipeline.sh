#!/bin/bash
set -euo pipefail

# ============================================================
# BASE DIRECTORY (auto-detected)
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${BASE:-$(cd "$SCRIPT_DIR/.." && pwd)}"

SCRIPTS="$SCRIPT_DIR"

EXTRACTED="${EXTRACTED:-$BASE/extracted_genes}"
WT="${WT_GENES:-${WT:-$BASE/WT_gene_fasta}}"
PER_GENE_FASTAS="$BASE/per_gene_fastas"
PER_GENE_FASTAS_DEDUP="$BASE/per_gene_fastas_dedup"
MACSE_OUT="$BASE/macse_output"
STOP_CALLS="$BASE/stop_calls"

MACSE_SUMMARY="$BASE/mutation_cds_macse_strict_summary.tsv"
FRAMESHIFTS="$BASE/frameshifts.tsv"
INDELS="$BASE/indels.tsv"
STOP_GENES="$BASE/stop_genes.tsv"

# ============================================================
# AUTO-DETECT MUTATION LIST
# ============================================================

if [ -z "${MUTATION_LIST:-}" ]; then
    for f in \
        "$BASE/mutation_list.tsv" \
        "$SCRIPT_DIR/../mutation_list.tsv" \
        "./mutation_list.tsv"
    do
        if [ -f "$f" ]; then
            MUTATION_LIST="$f"
            echo "   ✅ Found mutation list: $f"
            break
        fi
    done
fi

# Final validation
if [ -n "${MUTATION_LIST:-}" ] && [ ! -f "${MUTATION_LIST}" ]; then
    echo "   ⚠  MUTATION_LIST not found at ${MUTATION_LIST} — disabling mutation steps"
    MUTATION_LIST=""
fi

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
    files=( "$dir"/*.fasta )
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

done_ "$built FASTAs built${skipped:+, $skipped skipped}"

# ============================================================
# 1.5 DEDUPLICATION
# ============================================================

step "1.5" "Deduplicating FASTAs"
mkdir -p "$PER_GENE_FASTAS_DEDUP"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f")
    out="$PER_GENE_FASTAS_DEDUP/$base"

    seqkit rmdup -s "$f" > "$out"

    n=$(grep -c "^>" "$out")
    echo "   ${base%.fasta}: $n unique variants"
done

done_ "Deduplication complete"

# ============================================================
# 2. MACSE ALIGNMENT
# ============================================================

step "2" "Running MACSE alignment"
mkdir -p "$MACSE_OUT"

THREADS="${THREADS:-4}"

find "$PER_GENE_FASTAS_DEDUP" -maxdepth 1 -name "*.fasta" -print0 | \
xargs -0 -r -n 1 -P "$THREADS" bash -c '
    set -e
    f="$1"
    base=$(basename "$f" .fasta)

    if [ -f "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" ]; then
        echo "   ${base}: already aligned"
        exit 0
    fi

    echo "   Aligning: $base"

    macse \
        -prog alignSequences \
        -seq "$f" \
        -out_AA "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" \
        -out_NT "'"$MACSE_OUT"'/${base}_aligned_nt.fasta"
' _

done_ "MACSE complete"

# ============================================================
# 3. MUTATION SCAN
# ============================================================

step "3" "Scanning mutations"

python \
    "$SCRIPTS/scan_mutations_cds_macse.py" \
    "$MACSE_OUT" \
    > "$MACSE_SUMMARY"

done_ "Summary written"

# ============================================================
# 4. STOP MUTATIONS (always runs — needed by screening step)
# ============================================================

step "4" "Detecting stop mutations"
mkdir -p "$STOP_CALLS"
> "$STOP_CALLS/stop_calls.tsv"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    python \
        "$SCRIPTS/detect_stop_mutations.py" "$f" \
        >> "$STOP_CALLS/stop_calls.tsv"
done

awk 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' "$STOP_CALLS/stop_calls.tsv" \
    | sort -u > "$STOP_GENES"

done_ "Stop mutations complete"

# ============================================================
# 5–7 OPTIONAL STEPS (require mutation list)
# ============================================================

if [ -z "${MUTATION_LIST}" ]; then
    warn "No mutation list — skipping Steps 5–7"
    > "$FRAMESHIFTS"
    > "$INDELS"
    > "$BASE/mutation_screen_results.tsv"
else

    # ---------- Step 5 ----------
    step "5" "Frameshift detection"

    python \
        "$SCRIPTS/detect_frameshifts_macse.py" \
        "$MACSE_OUT" \
        "$MUTATION_LIST" \
        "$FRAMESHIFTS"

    done_ "$(wc -l < "$FRAMESHIFTS") frameshift hits → $(basename "$FRAMESHIFTS")"

    # ---------- Step 5.5 ----------
    step "5.5" "Indel detection"

    python \
        "$SCRIPTS/detect_indels_macse.py" \
        "$MACSE_OUT" \
        "$MUTATION_LIST" \
        "$INDELS"

    done_ "$(wc -l < "$INDELS") indel hits → $(basename "$INDELS")"

    # ---------- Step 6 ----------
    step "6" "Mutation screening"

    python \
        "$SCRIPTS/screen_mutations.py" \
        --mutations   "$MUTATION_LIST" \
        --macse       "$MACSE_SUMMARY" \
        --frameshifts "$FRAMESHIFTS" \
        --indels      "$INDELS" \
        --stops       "$STOP_GENES" \
        --output      "$BASE/mutation_screen_results.tsv"

    done_ "Results → mutation_screen_results.tsv"
fi

# ============================================================

echo
echo "══════════════════════════════════════════"
echo "  ✅ PIPELINE COMPLETE"
echo "  ➡  Output: $BASE/mutation_screen_results.tsv"
echo "══════════════════════════════════════════"
