#!/usr/bin/env python3
"""
debug_macse_windows.py  <macse_output_dir>

Prints the ancestor and isolate residues at each missense/delins window
so you can see exactly what amino acids are present before deciding on
the expected values in scan_mutations_cds_macse.py.

Usage:
    python debug_macse_windows.py /home/alasdair/3D_UHU_evo_analysis/macse_output
"""

import os
import sys
import glob
from collections import Counter
from Bio import SeqIO

# Only the mutations where we need to inspect the actual residues.
# (gene, aa_start, aa_end, mutation_type, currently_expected)
INSPECT = [
    ("glpT", 164, 164, "missense", "V"),
    ("cyaA", 394, 394, "missense", "D"),
    ("glpT", 182, 187, "delins",   "C"),
]


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


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: debug_macse_windows.py <macse_output_dir>")

    macse_dir = sys.argv[1]

    for g, s, e, m, expected in INSPECT:
        print(f"\n{'='*60}")
        print(f"Gene={g}  positions={s}-{e}  type={m}  expected='{expected}'")
        print(f"{'='*60}")

        pattern = os.path.join(macse_dir, f"{g}*_aligned_aa.fasta")
        aln_files = glob.glob(pattern)

        if not aln_files:
            print(f"  ⚠  No alignment files found matching {pattern}")
            continue

        # Collect all isolate windows across all alignment files
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
                print(f"  Ancestor residues at pos {s}-{e}: '{anc_w}'")
                # Also show a slightly wider window for context
                ctx_start = max(0, c0 - 3)
                ctx_end   = min(len(anc), c1 + 3)
                anc_ctx = anc[ctx_start:ctx_end].replace("-","").replace("!","")
                print(f"  Ancestor context (±3 aa):       '{anc_ctx}'")
                anc_window_shown = True

            for rec in records[1:]:
                iso = str(rec.seq)
                iso_w = iso[c0:c1].replace("-", "").replace("!", "")
                iso_windows[iso_w] += 1
                n_isolates += 1

        print(f"\n  Isolate residues at this window ({n_isolates} total isolates):")
        print(f"  {'Residues':<30} {'Count':>6}  {'% of isolates':>14}")
        print(f"  {'-'*54}")
        for residues, count in iso_windows.most_common(20):
            pct = 100 * count / n_isolates if n_isolates else 0
            marker = "  ← expected" if residues == expected else ""
            print(f"  {residues:<30} {count:>6}  {pct:>13.1f}%{marker}")

        if n_isolates > 20:
            remaining = len(iso_windows) - 20
            if remaining > 0:
                print(f"  ... and {remaining} more unique residue combinations")

        # For delins: show how many isolates have gaps in the window
        if m == "delins":
            print(f"\n  Delins check — isolates with gaps in window:")
            gap_count = sum(
                c for w, c in iso_windows.items() if "-" in w or "!" in w
            )
            print(f"  {gap_count} / {n_isolates} isolates have gaps at positions {s}-{e}")


if __name__ == "__main__":
    main()
