#!/bin/bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: bash move_ecoli.sh <source_dir> <dest_dir>" >&2
    exit 1
fi

SOURCE="$1"
DEST="$2"

if [ ! -d "$SOURCE" ]; then
    echo "ERROR: source directory not found: $SOURCE" >&2
    exit 1
fi

mkdir -p "$DEST"

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
        cp "$file" "$DEST/$newname"
        rm "$file"
    else
        cp "$file" "$DEST/$filename"
        rm "$file"
    fi
done

echo "✅ Files moved to $DEST"
