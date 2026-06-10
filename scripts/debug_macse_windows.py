#!/usr/bin/env python3
"""
debug_macse_windows.py  <macse_output_dir>  <mutation_list.tsv>

Prints the ancestor and isolate residues at each missense/delins window
defined in the mutation list so you can verify what amino acids are
present in the alignment before troubleshooting scan_mutations_cds_macse.py.

Usage:
    python debug_macse_windows.py <macse_output_dir> <mutation_list.tsv>
"""

import os
import sys
import csv
import glob
from collections import Counter
from Bio import SeqIO


def ungapped_to_alignment_index(seq, pos):
    count = 0
    for i, aa in enumerate(seq):
        if aa not in ("-", "!"):
            count += 1
            if count == pos:
                return i
    return None


def get_alignment_columns(anc_seq, start, end):
    col_start = ungapped_to_alignment_index(anc_seq, start)
    if col_start is None:
        return None, None
    span = end - start + 1
    count = 0
    col_end = col_start
    while col_end < len(anc_seq) and count < span:
        if anc_seq[col_end] not in ("-", "!"):
            count += 1
        col_end += 1
    return col_start, col_end


def load_inspect(mutation_file):
    """Load missense and delins entries from the mutation list."""
    inspect = []
    with open(mutation_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mtype = row.get("mutation_type", "").strip()
            if mtype not in ("missense", "delins"):
                continue
            gene = row.get("gene", "").strip()
            try:
                s = int(row["aa_start"])
                e = int(row["aa_end"])
            except (KeyError, ValueError):
                continue
            expected = row.get("expected", "").strip() or None
            inspect.append((gene, s, e, mtype, expected))
    return inspect


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: debug_macse_windows.py <macse_output_dir> <mutation_list.tsv>")

    macse_dir     = sys.argv[1]
    mutation_file = sys.argv[2]

    inspect = load_inspect(mutation_file)
    if not inspect:
        sys.exit("No missense or delins mutations found in mutation list.")

    for g, s, e, m, expected in inspect:
        print(f"\n{'='*60}")
        print(f"Gene={g}  positions={s}-{e}  type={m}  expected='{expected}'")
        print(f"{'='*60}")

        pattern = os.path.join(macse_dir, f"{g}*_aligned_aa.fasta")
        aln_files = glob.glob(pattern)

        if not aln_files:
            print(f"  ⚠  No alignment files found matching {pattern}")
            continue

        iso_windows = Counter()
        anc_window_shown = False
        n_isolates = 0

        for aln in aln_files:
            records = list(SeqIO.parse(aln, "fasta"))
            if len(records) < 2:
                continue

            anc = str(records[0].seq)
            c0, c1 = get_alignment_columns(anc, s, e)

            if c0 is None:
                print(f"  ⚠  Could not map position {s} in {os.path.basename(aln)}")
                continue

            if not anc_window_shown:
                anc_w = anc[c0:c1].replace("-", "").replace("!", "")
                ctx_start = max(0, c0 - 3)
                ctx_end   = min(len(anc), c1 + 3)
                anc_ctx = anc[ctx_start:ctx_end].replace("-", "").replace("!", "")
                print(f"  Ancestor residues at pos {s}-{e}: '{anc_w}'")
                print(f"  Ancestor context (±3 aa):       '{anc_ctx}'")
                anc_window_shown = True

            for rec in records[1:]:
                iso = str(rec.seq)
                iso_w = iso[c0:c1].replace("-", "").replace("!", "")
                iso_windows[iso_w] += 1
                n_isolates += 1

        print(f"\n  Isolate residues at this window ({n_isolates} total):")
        print(f"  {'Residues':<30} {'Count':>6}  {'%':>8}")
        print(f"  {'-'*50}")
        for residues, count in iso_windows.most_common(20):
            pct = 100 * count / n_isolates if n_isolates else 0
            marker = "  ← expected" if residues == expected else ""
            print(f"  {residues:<30} {count:>6}  {pct:>7.1f}%{marker}")

        if len(iso_windows) > 20:
            print(f"  ... and {len(iso_windows) - 20} more unique combinations")

        if m == "delins":
            gap_count = sum(c for w, c in iso_windows.items() if "-" in w or "!" in w)
            print(f"\n  Delins check: {gap_count}/{n_isolates} isolates have gaps at positions {s}-{e}")


if __name__ == "__main__":
    main()
