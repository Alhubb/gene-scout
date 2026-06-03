#!/usr/bin/env python3
"""
scan_mutations_cds_macse.py  <macse_output_dir>

Scans MACSE amino-acid alignments for the predefined mutations and writes
a summary TSV to stdout. run_pipeline.sh redirects stdout to MACSE_SUMMARY.

FIX summary vs previous version:
  1. Missense confirmation uses a two-tier approach:
       - present_exact=1 if the specific expected residue is found.
       - present_any=1 if ANY amino acid difference is found at that position.
     Both written to output TSV.
  2. confirm_any_missense() now returns False if the isolate window is empty
     after stripping gaps/frameshift markers. Previously an empty window
     (e.g. due to a deletion or assembly gap) would compare as unequal to
     the ancestor and incorrectly trigger any_missense=1.
  3. confirm_indel() searches a wider flanking window (±5 aa) for delins.
  4. Debug mode is always on: per-isolate window content printed to stderr.
     Redirect stderr to suppress: python scan_mutations_cds_macse.py macse_output 2>/dev/null
"""

import os
import sys
import glob
from Bio import SeqIO

DEBUG = True

# ============================================================
# Mutation definitions
# ============================================================
# Columns: (gene, aa_start, aa_end, mutation_type, expected_aa_or_insert)
# aa_start / aa_end are 1-based positions in the WT ungapped protein.
# Frameshift and stop entries are listed here for the merge table but
# are NOT scanned by MACSE — they are handled by the frameshift and
# stop-codon modules respectively.

MUTATIONS = [
    ("glpT", 256, 268, "inframe_deletion", None),
    ("uhpC", 440, 440, "stop_lost",         None),   # handled by stop module
    ("glpT", 182, 187, "delins",           "C"),
    ("glpT", 164, 164, "missense",         "V"),
    ("uhpT", 7,   7,   "stop_gained",      None),   # handled by stop module
    ("ptsI", 500, 500, "frameshift",       None),   # handled by frameshift module
    ("cyaA", 394, 394, "missense",         "D"),
    ("uhpT", 423, 423, "frameshift",       None),   # handled by frameshift module
    ("glpT", 294, 294, "frameshift",       None),   # handled by frameshift module
    ("uhpA", 41,  43,  "inframe_deletion", None),
]

# Mutation types handled by MACSE scanning
MACSE_TYPES = {"missense", "inframe_deletion", "delins"}

# How many extra alignment columns to search either side of the defined
# window when looking for an inserted residue in a delins. Increase this
# if MACSE consistently places the insertion just outside the window.
DELINS_FLANK = 5


# ============================================================
# Alignment helpers
# ============================================================

def ungapped_to_alignment_index(seq, pos):
    """Return the alignment column index corresponding to ungapped position pos (1-based)."""
    count = 0
    for i, aa in enumerate(seq):
        if aa not in ("-", "!"):
            count += 1
            if count == pos:
                return i
    return None


def get_alignment_columns(anc_seq, start, end):
    """
    Return (col_start, col_end) alignment column indices spanning ungapped
    positions start..end (1-based, inclusive) in anc_seq.
    """
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


# ============================================================
# Confirmation logic
# ============================================================

def confirm_any_missense(anc, iso, start, end):
    """
    Check that ANY amino-acid difference exists in the window.
    Returns False if the isolate window is empty after stripping gaps/frameshifts
    — an empty window means the position is absent in this isolate (e.g. due to
    a deletion or assembly gap) and should not be called as a missense.
    """
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return False
    anc_w = anc[c0:c1].replace("-", "").replace("!", "")
    iso_w = iso[c0:c1].replace("-", "").replace("!", "")
    # Only call missense if the isolate actually has residues at this position
    if not iso_w:
        return False
    return anc_w != iso_w


def get_isolate_residue(anc, iso, start, end):
    """Return the isolate residue string at the window (for logging)."""
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return "?"
    return iso[c0:c1].replace("-", "").replace("!", "")


def confirm_indel(anc, iso, start, end, inserted=None):
    """
    Confirm an in-frame deletion (or delins) in the isolate.

    Logic:
      - Ancestor window must be mostly ungapped (≥80% of expected residues).
      - Isolate must have at least as many gaps as the expected deletion size.
      - For delins, the inserted residue is searched in the window PLUS a
        flanking region of DELINS_FLANK columns on each side, to account for
        MACSE placing the insertion just outside the strict deletion boundary.
    """
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return False

    anc_w = anc[c0:c1]
    iso_w = iso[c0:c1]

    anc_ungapped = len([a for a in anc_w if a not in ("-", "!")])
    iso_gaps     = sum(1 for a in iso_w if a in ("-", "!"))
    expected_del = end - start + 1

    # Skip if ancestor window is too gappy to be reliable
    if anc_ungapped < expected_del * 0.8:
        return False

    # Isolate must show at least the expected number of deleted residues
    if iso_gaps < expected_del:
        return False

    if inserted:
        # FIX: search a wider flanking window for the inserted residue
        flank_start = max(0, c0 - DELINS_FLANK)
        flank_end   = min(len(iso), c1 + DELINS_FLANK)
        iso_wide = iso[flank_start:flank_end].replace("-", "").replace("!", "")
        return inserted in iso_wide

    return True


# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: scan_mutations_cds_macse.py <macse_output_dir>")

    macse_dir = sys.argv[1]

    present_exact = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}
    present_any   = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}

    for g, s, e, m, x in MUTATIONS:
        if m not in MACSE_TYPES:
            continue

        pattern = os.path.join(macse_dir, f"{g}*_aligned_aa.fasta")
        for aln in glob.glob(pattern):
            records = list(SeqIO.parse(aln, "fasta"))
            if len(records) < 2:
                continue

            anc = str(records[0].seq)

            for r in records[1:]:
                iso = str(r.seq)

                if m == "missense":
                    res = get_isolate_residue(anc, iso, s, e)
                    anc_res = get_isolate_residue(anc, anc, s, e)

                    # Tier 1: exact match — isolate has the specific expected AA
                    if res == x:
                        present_exact[(g, s, e, m)] = 1

                    # Tier 2: any missense — isolate differs from ancestor at all
                    if confirm_any_missense(anc, iso, s, e):
                        present_any[(g, s, e, m)] = 1

                    if DEBUG:
                        print(
                            f"DEBUG missense {g} {s}: anc='{anc_res}' iso='{res}' expected='{x}'  "
                            f"{'EXACT MATCH' if res == x else 'any missense' if res != anc_res else 'no change'}",
                            file=sys.stderr,
                        )

                elif m in ("inframe_deletion", "delins"):
                    if confirm_indel(anc, iso, s, e, x):
                        present_exact[(g, s, e, m)] = 1
                        if DEBUG:
                            print(
                                f"DEBUG {m} {g} {s}-{e}: confirmed in {r.id}",
                                file=sys.stderr,
                            )

    # Write to stdout
    print("gene\taa_start\taa_end\tmutation_type\tpresent\tany_missense")
    for g, s, e, m, _ in MUTATIONS:
        print(
            f"{g}\t{s}\t{e}\t{m}\t"
            f"{present_exact.get((g,s,e,m), 0)}\t"
            f"{present_any.get((g,s,e,m), 0)}"
        )


if __name__ == "__main__":
    main()
