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

PROTEIN_FASTAS="${BASE}/protein_fastas"
STOP_CALLS="${BASE}/stop_calls"
STOP_GENES="${BASE}/stop_genes.tsv"

SCRIPTS="$(dirname "$(realpath "$0")")"

###############################################################################
# STEP 1 — BUILD PER-GENE FASTAs
###############################################################################

echo
echo "==============================="
echo "1️⃣ Building per-gene FASTAs"
echo "==============================="

mkdir -p "$PER_GENE_FASTAS"

for dir in "$EXTRACTED"/*_minimap_hits; do
    [ -d "$dir" ] || continue

    name=$(basename "$dir")
    core="${name%_minimap_hits}"

    gene="${core%%_*}"
    lineage="${core#*_}"

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
echo "2️⃣ Deduplicating per-gene FASTAs"
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
# STEP 3 — MACSE (FRAME-PRESERVING MUTATIONS)
###############################################################################

echo
echo "==============================="
echo "3️⃣ Running MACSE (codon-aware)"
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
# STEP 4 — SCAN MACSE ALIGNMENTS
###############################################################################

echo
echo "==============================="
echo "4️⃣ Scanning MACSE alignments"
echo "==============================="

python "$SCRIPTS/scan_mutations_cds_macse.py" "$MACSE_OUT" \
  > "${BASE}/mutation_cds_macse_strict_summary.tsv"

###############################################################################
# STEP 5 — FRAMESHIFT DETECTION (minimap2 evidence)
###############################################################################

echo
echo "==============================="
echo "5️⃣ Detecting frameshifts"
echo "==============================="

find "$MINIMAP_OUT" -name "*.paf" -exec awk '
{
    cigar=""; fs=0
    for (i=1;i<=NF;i++)
        if ($i ~ /^cg:Z:/) cigar=substr($i,6)

    while (match(cigar,/[0-9]+[ID]/)) {
        n=substr(cigar,RSTART,RLENGTH-1)
        if (n%3!=0) fs=1
        cigar=substr(cigar,RSTART+RLENGTH)
    }

    if (fs==1) {
        split(FILENAME,a,"/")
        split(a[length(a)-1],b,"_")
        print b[1]
    }
}' {} + | sort -u > "${BASE}/frameshifts.tsv"

###############################################################################
# STEP 6 — STOP-GAIN / STOP-LOSS
###############################################################################

echo
echo "==============================="
echo "6️⃣ Detecting stop-gain / stop-loss"
echo "==============================="

mkdir -p "$PROTEIN_FASTAS" "$STOP_CALLS"

> "${STOP_CALLS}/stop_calls.tsv"

for f in "$PER_GENE_FASTAS"/*.fasta; do
    base=$(basename "$f" .fasta)

    transeq -sequence "$f" \
      -outseq "${PROTEIN_FASTAS}/${base}.faa"

    python "$SCRIPTS/detect_stop_mutations.py" \
      "${PROTEIN_FASTAS}/${base}.faa" \
      >> "${STOP_CALLS}/stop_calls.tsv"
done

awk '{print $1}' "${STOP_CALLS}/stop_calls.tsv" | sort -u > "$STOP_GENES"

###############################################################################
# STEP 7 — FINAL MERGE
###############################################################################

echo
echo "==============================="
echo "7️⃣ Final merge"
echo "==============================="

python "$SCRIPTS/merge_mutation_evidence.py" \
  "${BASE}/mutation_cds_macse_strict_summary.tsv" \
  "${BASE}/frameshifts.tsv" \
  "$STOP_GENES" \
  > "${BASE}/mutation_final_evidence.tsv"

echo
echo "✅ MUTATION PIPELINE COMPLETE"
echo "➡ Output: ${BASE}/mutation_final_evidence.tsv"
``
