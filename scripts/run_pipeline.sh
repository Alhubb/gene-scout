#!/bin/bash

set -euo pipefail

###############################################################################
# REQUIRED ENVIRONMENT VARIABLES
###############################################################################

: "${BASE:?Need BASE}"
: "${WT_GENES:?Need WT_GENES}"

###############################################################################
# PATHS
###############################################################################

EXTRACTED="${BASE}/extracted_genes"
MINIMAP_OUT="${BASE}/minimap2_output"

PER_GENE_FASTAS="${BASE}/per_gene_fastas"
PER_GENE_FASTAS_DEDUP="${BASE}/per_gene_fastas_dedup"
MACSE_OUT="${BASE}/macse_output"

STOP_CALLS="${BASE}/stop_calls"
STOP_GENES="${BASE}/stop_genes.tsv"
FRAMESHIFTS="${BASE}/frameshifts.tsv"
MACSE_SUMMARY="${BASE}/mutation_cds_macse_strict_summary.tsv"
FINAL_OUT="${BASE}/mutation_final_evidence.tsv"

SCRIPTS="$(dirname "$(realpath "$0")")"

###############################################################################
# STEP 1 — BUILD PER-GENE FASTAs
###############################################################################

echo
echo "==============================="
echo "1️⃣  Building per-gene FASTAs"
echo "==============================="

mkdir -p "$PER_GENE_FASTAS"

for dir in "$EXTRACTED"/*_minimap_hits; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    # FIX: split on LAST underscore, not first.
    # Original: gene="${core%%_*}" / lineage="${core#*_}"
    # This silently truncates gene names that contain an underscore.
    # Splitting on the last underscore keeps multi-part names intact.
    lineage="${core##*_}"   # e.g. 34
    gene="${core%_*}"       # e.g. cyaA

    wt_fasta="${WT_GENES}/${gene}_${lineage}.fasta"
    out_fasta="${PER_GENE_FASTAS}/${gene}_${lineage}.fasta"

    if [ ! -f "$wt_fasta" ]; then
        echo "⚠ Missing WT ${gene}_${lineage}, skipping"
        continue
    fi

    shopt -s nullglob
    files=( "$dir"/*.bakta."$gene".fasta )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "⚠ No mutant CDS for ${gene}_${lineage}, skipping"
        continue
    fi

    echo "✅ Building $out_fasta"

    {
        cat "$wt_fasta"
        echo
        for f in "${files[@]}"; do
            cat "$f"
            echo
        done
    } > "$out_fasta"
done

###############################################################################
# STEP 2 — DEDUPLICATE PER-GENE FASTAs
###############################################################################

echo
echo "==============================="
echo "2️⃣  Deduplicating per-gene FASTAs"
echo "==============================="

mkdir -p "$PER_GENE_FASTAS_DEDUP"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f")
    out="${PER_GENE_FASTAS_DEDUP}/${base}"
    seqkit rmdup -s "$f" > "$out"
    n=$(grep -c "^>" "$out")
    echo "✅ $base → $n unique CDS variants"
done

###############################################################################
# STEP 3 — MACSE (CODON-AWARE ALIGNMENT)
###############################################################################

echo
echo "==============================="
echo "3️⃣  Running MACSE (codon-aware)"
echo "==============================="

mkdir -p "$MACSE_OUT"

find "$PER_GENE_FASTAS_DEDUP" -name "*.fasta" -print0 | \
xargs -0 -r -n 1 -P 4 bash -c '
    set -o noglob
    f="$1"
    base=$(basename "$f" .fasta)

    if [ -f "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" ]; then
        echo "Skipping $base (already aligned)"
        exit 0
    fi

    echo "Running MACSE on $base"

    macse -prog alignSequences \
      -seq "$f" \
      -out_AA "'"$MACSE_OUT"'/${base}_aligned_aa.fasta" \
      -out_NT "'"$MACSE_OUT"'/${base}_aligned_nt.fasta"
' _

###############################################################################
# STEP 4 — SCAN MACSE ALIGNMENTS FOR MISSENSE / IN-FRAME INDELS
###############################################################################

echo
echo "==============================="
echo "4️⃣  Scanning MACSE alignments"
echo "==============================="

python "$SCRIPTS/scan_mutations_cds_macse.py" "$MACSE_OUT" \
    > "$MACSE_SUMMARY"

###############################################################################
# STEP 5 — FRAMESHIFT DETECTION FROM CIGAR STRINGS
###############################################################################

echo
echo "==============================="
echo "5️⃣  Detecting frameshifts"
echo "==============================="

> "$FRAMESHIFTS"

for dir in "$MINIMAP_OUT"/*_minimap_hits; do
    [ -d "$dir" ] || continue

    dirname=$(basename "$dir")
    core="${dirname%_minimap_hits}"
    gene="${core%_*}"   # underscore-safe: split on last _

    for paf in "$dir"/*.paf; do
        [ -f "$paf" ] || continue

        awk -v gene="$gene" -v src="$paf" '
        {
            cigar = ""
            for (i = 1; i <= NF; i++)
                if ($i ~ /^cg:Z:/) cigar = substr($i, 6)
            if (cigar == "") next

            fs = 0
            tmp = cigar
            while (match(tmp, /[0-9]+[ID]/)) {
                op  = substr(tmp, RSTART, RLENGTH)
                n   = substr(op, 1, length(op)-1) + 0
                if (n % 3 != 0) fs = 1
                tmp = substr(tmp, RSTART + RLENGTH)
            }
            if (fs) print src "\t" gene
        }' "$paf" >> "$FRAMESHIFTS"
    done
done

# Deduplicate (same gene may appear from multiple PAF alignments)
sort -u "$FRAMESHIFTS" -o "$FRAMESHIFTS"

echo "✅ Frameshifts written to $FRAMESHIFTS"

###############################################################################
# STEP 6 — STOP-GAIN / STOP-LOSS
###############################################################################

echo
echo "==============================="
echo "6️⃣  Detecting stop-gain / stop-loss"
echo "==============================="

mkdir -p "$STOP_CALLS"
> "${STOP_CALLS}/stop_calls.tsv"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    python "$SCRIPTS/detect_stop_mutations.py" "$f" \
        >> "${STOP_CALLS}/stop_calls.tsv"
done

# stop_calls.tsv format: isolate_id <tab> event_type
# Deduplicate and write to STOP_GENES
awk 'NF >= 2 { print $1 "\t" $2 }' "${STOP_CALLS}/stop_calls.tsv" \
    | sort -u > "$STOP_GENES"

echo "✅ Stop mutations written to $STOP_GENES"

###############################################################################
# STEP 7 — FINAL MERGE
###############################################################################

echo
echo "==============================="
echo "7️⃣  Final merge"
echo "==============================="

python "$SCRIPTS/merge_mutation_evidence.py" \
    "$MACSE_SUMMARY" \
    "$FRAMESHIFTS" \
    "$STOP_GENES" \
    > "$FINAL_OUT"

###############################################################################
# STEP 8 — OPTIONAL MUTATION SCREENING
###############################################################################

if [ -n "${MUTATION_LIST:-}" ]; then
    echo
    echo "==============================="
    echo "8️⃣  Screening predefined mutations"
    echo "==============================="

    python "$SCRIPTS/screen_mutations.py" \
        --mutations "$MUTATION_LIST" \
        --macse     "$MACSE_SUMMARY" \
        --frameshifts "$FRAMESHIFTS" \
        --stops     "$STOP_GENES" \
        --output    "${BASE}/mutation_screen_results.tsv"

    echo "✅ Screening results: ${BASE}/mutation_screen_results.tsv"
else
    echo
    echo "ℹ No MUTATION_LIST supplied — discovery mode only"
fi

echo
echo "✅ PIPELINE COMPLETE"
echo "➡ Final evidence table: $FINAL_OUT"
