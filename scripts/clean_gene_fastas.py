#!/usr/bin/env python3
"""
clean_gene_fastas.py  <input_dir>  <output_dir>

Cleans per-gene deduplicated FASTA files before downstream screening.
Removes sequences that would cause false positives or confound detection:

  1. All-gap / near-empty sequences — contigs that mapped but contain
     almost no sequence (likely very short or failed assemblies).
     Threshold: sequence must have residues ≥ MIN_FRACTION of WT length.

  2. Non-IUPAC sequences — sequences containing characters other than
     ACGTN (corrupted or mis-assembled contigs).

  3. Sequences shorter than MIN_ABS_LENGTH nucleotides — absolute
     minimum to be a credible CDS extraction.

What is NOT removed at this stage:
  - Frameshifted sequences (length not divisible by 3) — needed for
    detect_frameshifts_macse.py
  - Sequences with early stop codons — needed for detect_stop_mutations.py
  - Sequences with inframe deletions — needed for detect_indels_macse.py

These are handled by type-specific detection scripts downstream.

Usage:
    python clean_gene_fastas.py <input_dir> <output_dir>

Output:
    One cleaned FASTA per gene in <output_dir>, plus a summary to stderr.
"""

import os
import sys
import glob

# ============================================================
# Configuration
# ============================================================

# Minimum fraction of WT sequence length an isolate must have
MIN_FRACTION = 0.5

# Absolute minimum nucleotide length
MIN_ABS_LENGTH = 100

# Valid nucleotide characters
VALID = set("ACGTNacgtn")


# ============================================================
# Helpers
# ============================================================

def parse_fasta(path):
    """Parse a FASTA file, return list of (header, seq) tuples."""
    records = []
    with open(path) as fh:
        current_header, current_seq = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            records.append((current_header, "".join(current_seq)))
    return records


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def clean_seq(seq):
    """Return the sequence stripped of gaps and MACSE markers."""
    return seq.replace("-", "").replace("!", "").upper()


def is_valid(seq, wt_len):
    """
    Return (keep, reason) — True if the sequence passes all quality filters.
    """
    raw = clean_seq(seq)

    # Absolute minimum length
    if len(raw) < MIN_ABS_LENGTH:
        return False, f"too_short ({len(raw)} nt)"

    # Minimum fraction of WT length
    if len(raw) < wt_len * MIN_FRACTION:
        return False, f"short_fraction ({len(raw)}/{wt_len} = {len(raw)/wt_len:.1%})"

    # Non-IUPAC characters
    invalid = set(raw) - VALID
    if invalid:
        return False, f"non_IUPAC ({invalid})"

    return True, "ok"


# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: clean_gene_fastas.py <input_dir> <output_dir>")

    in_dir  = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    fasta_files = sorted(glob.glob(os.path.join(in_dir, "*.fasta")))
    if not fasta_files:
        sys.exit(f"ERROR: no .fasta files found in {in_dir}")

    total_in = total_out = total_removed = 0

    for fasta_path in fasta_files:
        base     = os.path.basename(fasta_path)
        out_path = os.path.join(out_dir, base)

        records = parse_fasta(fasta_path)
        if len(records) < 2:
            # Only WT, nothing to filter — copy as-is
            write_fasta(records, out_path)
            continue

        # First record is always the WT ancestor
        wt_header, wt_seq = records[0]
        wt_len = len(clean_seq(wt_seq))

        kept    = []
        removed = []

        for header, seq in records[1:]:
            keep, reason = is_valid(seq, wt_len)
            if keep:
                kept.append((header, seq))
            else:
                removed.append((header, reason))

        n_in      = len(records) - 1
        n_kept    = len(kept)
        n_removed = len(removed)

        total_in      += n_in
        total_out     += n_kept
        total_removed += n_removed

        write_fasta([(wt_header, wt_seq)] + kept, out_path)

        gene_ref = base.replace(".fasta", "")
        if n_removed > 0:
            print(f"  {gene_ref}: {n_in} → {n_kept} "
                  f"(removed {n_removed})", file=sys.stderr)
            for h, r in removed[:3]:  # show first 3 removed
                print(f"    ✗ {h[:60]}: {r}", file=sys.stderr)
            if n_removed > 3:
                print(f"    ... and {n_removed - 3} more", file=sys.stderr)
        else:
            print(f"  {gene_ref}: {n_in} → {n_kept} ✅", file=sys.stderr)

    print(f"\n✅ Cleaning complete: "
          f"{total_in} sequences in → "
          f"{total_out} kept, "
          f"{total_removed} removed",
          file=sys.stderr)


if __name__ == "__main__":
    main()
