#!/bin/bash
# FIX: Added set -euo pipefail. Without it, a failed mv (e.g. cross-filesystem
# interruption) continues silently and the missing file is never flagged.
set -euo pipefail

SOURCE="/media/alasdair/External/3D_UHU/bakrep/ecoli"
DEST="/media/alasdair/External/3D_UHU/bakrep/fasta"

mkdir -p "$DEST"

# NOTE: This script deduplicates by *filename* only. Genome files with
# different names but identical content will both be moved and processed
# downstream. If you want to avoid redundant computation, run seqkit rmdup
# or a checksum-based dedup on the assembled FASTA sequences after this step.

find "$SOURCE" -type f -print0 | while IFS= read -r -d '' file; do
    filename=$(basename "$file")

    if [ -e "$DEST/$filename" ]; then
        base="${filename%.*}"
        ext="${filename##*.}"
        counter=1
        while [ -e "$DEST/${base}_${counter}.${ext}" ]; do
            ((counter++))
        done
        newname="${base}_${counter}.${ext}"
        echo "  Name collision: renaming to $newname"
        # FIX: Use cp + rm (copy-then-delete) rather than bare mv so that a
        # cross-filesystem transfer that is interrupted mid-way does not leave
        # a half-written file in DEST with the source already deleted.
        cp "$file" "$DEST/$newname"
        rm "$file"
    else
        cp "$file" "$DEST/$filename"
        rm "$file"
    fi
done

echo "✅ Files moved to $DEST"
