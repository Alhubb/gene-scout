#!/bin/bash
set -euo pipefail

# ============================================================
# BASE DIRECTORY (auto-detected)
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${BASE:-$(cd "$SCRIPT_DIR/.." && pwd)}"

SCRIPTS="$SCRIPT_DIR"

EXTRACTED="${EXTRACTED:-$BASE/extracted_genes}"
WT="${WT:-$BASE/WT_gene_fasta}"
PER_GENE_FASTAS="$BASE/per_gene_fastas"
PER_GENE_FASTAS_DEDUP="$BASE/per_gene_fastas_dedup"
MACSE_OUT="$BASE/macse_output"
STOP_CALLS="$BASE/stop_calls"

MACSE_SUMMARY="$BASE/mutation_cds_macse_strict_summary.tsv"
FRAMESHIFTS="$BASE/frameshifts.tsv"
INDELS="$BASE/indels.tsv"
STOP_GENES="$BASE/stop_genes.tsv"

# ============================================================
# CHECK DEPENDENCIES
# ============================================================

echo "Checking dependencies..."

for cmd in seqkit macse python; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "❌ ERROR: '$cmd' not found in PATH"
        echo "👉 Activate the gene-scout environment first"
        exit 1
    fi
done

echo "   ✅ All dependencies found"

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
# 1. BUILD FASTAS
# ============================================================

step "1" "Building per-gene FASTAs"
mkdir -p "$PER_GENE_FASTAS"

for dir in "$EXTRACTED"/*; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"
    lineage="${core##*_}"
    gene="${core%_*}"

    wt_fasta="$WT/${gene}_${lineage}.fasta"
    out_fasta="$PER_GENE_FASTAS/${gene}_${lineage}.fasta"

    if [ ! -f "$wt_fasta" ]; then
        warn "Missing WT: ${gene}_${lineage}"
        continue
    fi

    files=( "$dir"/*.bakta."$gene".fasta )
    [ ${#files[@]} -eq 0 ] && continue

    { cat "$wt_fasta"; echo
      for f in "${files[@]}"; do cat "$f"; echo; done
    } > "$out_fasta"
done

done_ "FASTAs built"

# ============================================================
# 1.5 DEDUPLICATION
# ============================================================

step "1.5" "Deduplicating"
mkdir -p "$PER_GENE_FASTAS_DEDUP"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f")
    out="$PER_GENE_FASTAS_DEDUP/$base"

    seqkit rmdup -s "$f" > "$out"

    n=$(grep -c "^>" "$out")
    echo "   ${base%.fasta}: $n variants"
done

done_ "Dedup complete"

# ============================================================
# 2. MACSE ALIGNMENT
# ============================================================

step "2" "Running MACSE"
mkdir -p "$MACSE_OUT"

for f in "$PER_GENE_FASTAS_DEDUP"/*.fasta; do
    base=$(basename "$f" .fasta)

    if [ -f "$MACSE_OUT/${base}_aligned_aa.fasta" ]; then
        echo "   ${base}: already aligned"
        continue
    fi

    echo "   Aligning: $base"

    macse \
        -prog alignSequences \
        -seq "$f" \
        -out_AA "$MACSE_OUT/${base}_aligned_aa.fasta" \
        -out_NT "$MACSE_OUT/${base}_aligned_nt.fasta"
done

done_ "MACSE complete"

# ============================================================
# 3. MUTATION SCAN
# ============================================================

step "3" "Scanning mutations"

python "$SCRIPTS/scan_mutations_cds_macse.py" \
    "$MACSE_OUT" \
    > "$MACSE_SUMMARY"

done_ "Summary written"

# ============================================================
# OPTIONAL STEPS
# ============================================================

if [ -z "${MUTATION_LIST}" ]; then
    warn "No mutation list — skipping Steps 4–6"
    > "$FRAMESHIFTS"
    > "$INDELS"
    > "$BASE/mutation_screen_results.tsv"
else

    step "4" "Frameshifts"
    python "$SCRIPTS/detect_frameshifts_macse.py" \
        "$MACSE_OUT" "$MUTATION_LIST" "$FRAMESHIFTS"

    step "4.5" "Indels"
    python "$SCRIPTS/detect_indels_macse.py" \
        "$MACSE_OUT" "$MUTATION_LIST" "$INDELS"

    step "6" "Screening"
    python "$SCRIPTS/screen_mutations.py" \
        --mutations "$MUTATION_LIST" \
        --macse "$MACSE_SUMMARY" \
        --frameshifts "$FRAMESHIFTS" \
        --indels "$INDELS" \
        --stops "$STOP_GENES" \
        --output "$BASE/mutation_screen_results.tsv"

    done_ "Mutation screening complete"
fi

# ============================================================
# 5. STOP MUTATIONS
# ============================================================

step "5" "Stop mutations"
mkdir -p "$STOP_CALLS"
> "$STOP_CALLS/stop_calls.tsv"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    python "$SCRIPTS/detect_stop_mutations.py" "$f" \
        >> "$STOP_CALLS/stop_calls.tsv"
done

awk 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' \
    "$STOP_CALLS/stop_calls.tsv" | sort -u > "$STOP_GENES"

done_ "Stop mutations complete"

# ============================================================

echo
echo "══════════════════════════════════════════"
echo "  ✅ PIPELINE COMPLETE"
echo "  ➡ Output: $BASE/mutation_screen_results.tsv"
echo "══════════════════════════════════════════"
