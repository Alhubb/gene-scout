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
PER_GENE_FASTAS_CLEAN="$BASE/per_gene_fastas_clean"
MACSE_OUT="$BASE/macse_output"
MACSE_OUT_CLEAN="$BASE/macse_output_clean"
STOP_CALLS="$BASE/stop_calls"
CONFIRMED_SEQS="$BASE/confirmed_sequences"

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

if [ -n "${MUTATION_LIST:-}" ]; then
    echo "   ✅ Mutation list: $MUTATION_LIST"
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

    if [[ "$core" == *_* ]]; then
        gene="${core%_*}"
        lineage="${core##*_}"
    else
        gene="$core"
        lineage="unknown"
    fi

    out_fasta="$PER_GENE_FASTAS/${gene}_${lineage}.fasta"

    wt_fasta="$WT/${gene}_${lineage}.fasta"
    if [ ! -f "$wt_fasta" ]; then
        candidate="$WT/${gene}.fasta"
        if [ -f "$candidate" ]; then
            wt_fasta="$candidate"
        else
            warn "Missing WT for ${gene} (${lineage}) — skipping"
            skipped=$((skipped + 1))
            continue
        fi
    fi

    shopt -s nullglob
    files=( "$dir"/*bakta.${gene}.fasta )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        warn "No mutant CDS files in $dir — skipping"
        skipped=$((skipped + 1))
        continue
    fi

    {
        cat "$wt_fasta"; echo
        for f in "${files[@]}"; do
            cat "$f"; echo
        done
    } > "$out_fasta"

    built=$((built + 1))
done

done_ "$built FASTAs built${skipped:+, $skipped skipped}"

# ============================================================
# 1.5 DEDUPLICATION
# ============================================================

step "1.5" "Deduplicating FASTAs"
mkdir -p "$PER_GENE_FASTAS_DEDUP"

shopt -s nullglob
for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f")
    out="$PER_GENE_FASTAS_DEDUP/$base"

    seqkit rmdup -s "$f" > "$out"

    n=$(grep -c "^>" "$out")
    echo "   ${base%.fasta}: $n unique variants"
done
shopt -u nullglob

done_ "Deduplication complete"

# ============================================================
# 1.6 CLEAN PER-GENE FASTAs
# Removes sequences that are too short, mostly gaps, or contain
# non-IUPAC characters. Frameshifted and stop-codon sequences
# are kept — they are handled by the detection scripts downstream.
# ============================================================

step "1.6" "Cleaning per-gene FASTAs"
mkdir -p "$PER_GENE_FASTAS_CLEAN"

python "$SCRIPTS/clean_gene_fastas.py" \
    "$PER_GENE_FASTAS_DEDUP" \
    "$PER_GENE_FASTAS_CLEAN" \
    2>&1

done_ "Cleaning complete → $(basename "$PER_GENE_FASTAS_CLEAN")"

step "2" "Running MACSE alignment"
mkdir -p "$MACSE_OUT"

THREADS="${THREADS:-4}"

find "$PER_GENE_FASTAS_CLEAN" -maxdepth 1 -name "*.fasta" -print0 | \
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
# 2.5 CLEAN MACSE ALIGNMENTS
# Removes sequences that MACSE rendered as near-entirely gaps
# after codon-aware alignment (< 50% of ancestor ungapped length).
# For stop_lost genes, also removes sequences too short to reach
# the WT stop position. Both AA and NT files cleaned in lockstep.
# All downstream steps use MACSE_OUT_CLEAN.
# ============================================================

step "2.5" "Cleaning MACSE alignments"
mkdir -p "$MACSE_OUT_CLEAN"

# Derive stop_lost gene names from mutation list if available
STOP_LOST_GENES=""
if [ -n "${MUTATION_LIST:-}" ]; then
    STOP_LOST_GENES=$(awk -F'\t' '$4=="stop_lost" {print $1}' "$MUTATION_LIST" \
        | sort -u | paste -sd,)
fi

python "$SCRIPTS/clean_macse_alignments.py" \
    "$MACSE_OUT" \
    "$MACSE_OUT_CLEAN" \
    ${STOP_LOST_GENES:+--stop-lost-genes "$STOP_LOST_GENES"} \
    2>&1

done_ "MACSE cleaning complete → $(basename "$MACSE_OUT_CLEAN")"

# ============================================================
# 3. MUTATION SCAN
# ============================================================

step "3" "Scanning mutations"

if [ -z "${MUTATION_LIST}" ]; then
    warn "No mutation list — skipping MACSE scan"
    > "$MACSE_SUMMARY"
else
    python "$SCRIPTS/scan_mutations_cds_macse.py" \
        "$MACSE_OUT_CLEAN" \
        "$MUTATION_LIST" \
        > "$MACSE_SUMMARY"
fi

done_ "Summary written"

# ============================================================
# 4. STOP MUTATIONS
# ============================================================

step "4" "Detecting stop mutations"
mkdir -p "$STOP_CALLS"
> "$STOP_CALLS/stop_calls.tsv"

shopt -s nullglob
for f in "$PER_GENE_FASTAS_CLEAN"/*.fasta; do
    if [ -n "${MUTATION_LIST:-}" ]; then
        python "$SCRIPTS/detect_stop_mutations.py" "$f" "$MUTATION_LIST" \
            >> "$STOP_CALLS/stop_calls.tsv"
    else
        python "$SCRIPTS/detect_stop_mutations.py" "$f" \
            >> "$STOP_CALLS/stop_calls.tsv"
    fi
done
shopt -u nullglob

awk 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' \
    "$STOP_CALLS/stop_calls.tsv" | sort -u > "$STOP_GENES"

done_ "Stop mutations complete"

# ============================================================
# OPTIONAL STEPS
# ============================================================

if [ -z "${MUTATION_LIST}" ]; then
    warn "No mutation list — skipping Steps 5–7"
    > "$FRAMESHIFTS"
    > "$INDELS"
    > "$BASE/mutation_screen_results.tsv"
else

    step "5" "Frameshift detection"
    python "$SCRIPTS/detect_frameshifts_macse.py" \
        "$MACSE_OUT_CLEAN" "$MUTATION_LIST" "$FRAMESHIFTS"
    done_ "$(wc -l < "$FRAMESHIFTS") frameshift hits"

    step "5.5" "Indel detection"
    python "$SCRIPTS/detect_indels_macse.py" \
        "$MACSE_OUT_CLEAN" "$MUTATION_LIST" "$INDELS"
    done_ "$(wc -l < "$INDELS") indel hits"

    step "6" "Mutation screening"
    python "$SCRIPTS/screen_mutations.py" \
        --mutations   "$MUTATION_LIST" \
        --macse       "$MACSE_SUMMARY" \
        --frameshifts "$FRAMESHIFTS" \
        --indels      "$INDELS" \
        --stops       "$STOP_GENES" \
        --output      "$BASE/mutation_screen_results.tsv"
    done_ "Mutation screening complete"

    step "7" "Extracting confirmed sequences"
    mkdir -p "$CONFIRMED_SEQS"
    python "$SCRIPTS/extract_confirmed_sequences.py" \
        "$MACSE_OUT_CLEAN" \
        "$MACSE_SUMMARY" \
        "$MUTATION_LIST" \
        "$CONFIRMED_SEQS" \
        --frameshifts     "$FRAMESHIFTS" \
        --stops           "$STOP_GENES" \
        --per-gene-fastas "$PER_GENE_FASTAS_CLEAN" \
        --screen-results  "$BASE/mutation_screen_results.tsv" \
        --any
    done_ "Confirmed sequences → $(basename "$CONFIRMED_SEQS")"

    step "8" "Counting confirmed evidence"
    python "$SCRIPTS/count_confirmed_evidence.py" \
        "$CONFIRMED_SEQS" \
        "$CONFIRMED_SEQS/confirmed_evidence_counts.tsv"
    done_ "Evidence counts → confirmed_evidence_counts.tsv"

fi

# ============================================================

echo
echo "══════════════════════════════════════════"
echo "  ✅ PIPELINE COMPLETE"
echo "  ➡ Output: $BASE/mutation_screen_results.tsv"
echo "══════════════════════════════════════════"
