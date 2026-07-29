#!/usr/bin/env python3

import os
import sys
import csv
import glob
from Bio import SeqIO

DEBUG = False  # set True if you want verbose debug output

# ============================================================
# Load mutation list
# ============================================================

def load_mutations(tsv_file):
    """
    Load mutation list TSV. Expected columns:
        gene  aa_start  aa_end  mutation_type  [expected]

    The 'expected' column is optional — it should contain the single-letter
    amino acid expected in isolates for missense mutations, or the inserted
    residue for delins. If absent or blank, exact matching is skipped and
    any missense at that position counts as a hit.

    Uses csv.DictReader so column order doesn't matter and the 'expected'
    column is found by name regardless of position.
    """
    mutations = []

    with open(tsv_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene  = row.get("gene", "").strip()
            mtype = row.get("mutation_type", "").strip()
            if not gene or not mtype:
                continue
            try:
                start = int(row["aa_start"])
                end   = int(row["aa_end"])
            except (KeyError, ValueError):
                continue

            expected_raw = row.get("expected", "").strip()
            expected = None if expected_raw in ("", ".", "NA", "None") else expected_raw

            mutations.append((gene, start, end, mtype, expected))

    return mutations


# ============================================================
# Constants
# ============================================================

MACSE_TYPES = {"missense", "inframe_deletion", "delins"}
DELINS_FLANK = 5

# ============================================================
# Alignment helpers
# ============================================================

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


# ============================================================
# Mutation logic
# ============================================================

# Positional tolerance for missense detection — number of residues either
# side of the called position to search. Handles off-by-one shifts from
# MACSE inserting a spurious ! at alignment position 1 for some isolates,
# which shifts all downstream positions by 1 in the ungapped coordinate system.
POSITION_TOLERANCE = 2


def confirm_any_missense(anc, iso, start, end):
    """
    Return True if ANY amino-acid difference exists within ±POSITION_TOLERANCE
    of the called position. Handles MACSE leading ! artefacts that shift
    positions by 1. Returns False if no difference found in any position checked.
    """
    for offset in range(-POSITION_TOLERANCE, POSITION_TOLERANCE + 1):
        s = max(1, start + offset)
        e = max(1, end + offset)
        c0, c1 = get_alignment_columns(anc, s, e)
        if c0 is None:
            continue
        anc_w = anc[c0:c1].replace("-", "").replace("!", "")
        iso_w = iso[c0:c1].replace("-", "").replace("!", "")
        if not iso_w:
            continue
        if anc_w != iso_w:
            return True
    return False


def get_isolate_residue(anc, iso, start, end):
    """
    Return isolate residue at the called position.
    Falls back within ±POSITION_TOLERANCE if exact position maps to a gap.
    """
    for offset in range(0, POSITION_TOLERANCE + 1):
        for s, e in [(start + offset, end + offset),
                     (start - offset, end - offset)]:
            s = max(1, s)
            e = max(1, e)
            c0, c1 = get_alignment_columns(anc, s, e)
            if c0 is None:
                continue
            res = iso[c0:c1].replace("-", "").replace("!", "")
            if res:
                return res
    return "?"


def confirm_indel(anc, iso, start, end, inserted=None):
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return False

    anc_w = anc[c0:c1]
    iso_w = iso[c0:c1]

    anc_ungapped = len([a for a in anc_w if a not in ("-", "!")])
    iso_gaps = sum(1 for a in iso_w if a in ("-", "!"))
    expected_del = end - start + 1

    if anc_ungapped < expected_del * 0.8:
        return False

    if iso_gaps < expected_del:
        return False

    if inserted:
        flank_start = max(0, c0 - DELINS_FLANK)
        flank_end   = min(len(iso), c1 + DELINS_FLANK)

        iso_wide = iso[flank_start:flank_end].replace("-", "").replace("!", "")
        return inserted in iso_wide

    return True


# ============================================================
# Main
# ============================================================

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]

    if len(args) != 2:
        sys.exit("Usage: scan_mutations_cds_macse.py <macse_dir> <mutation_list.tsv> [--debug]")

    macse_dir     = args[0]
    mutation_file = args[1]
    debug         = "--debug" in flags

    MUTATIONS = load_mutations(mutation_file)

    if not MUTATIONS:
        print("[ERROR] No mutations loaded", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(MUTATIONS)} mutations", file=sys.stderr)

    # ✅ collect all alignments
    all_alignments = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))
    print(f"[INFO] Found {len(all_alignments)} alignment files", file=sys.stderr)

    present_exact = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}
    present_any   = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}

    # ============================================================
    # Scan mutations
    # ============================================================

    for g, s, e, m, x in MUTATIONS:

        if m not in MACSE_TYPES:
            continue

        # Match alignment files by gene name prefix with underscore boundary
        # e.g. 'cyaA' matches 'cyaA_34_aligned_aa.fasta' but NOT 'cyaABC_34...'
        # Falls back to case-insensitive substring only if no prefix hits found.
        matching_files = [
            f for f in all_alignments
            if os.path.basename(f).startswith(g + "_")
        ]
        if not matching_files:
            matching_files = [
                f for f in all_alignments
                if g.lower() in os.path.basename(f).lower()
            ]

        print(f"[INFO] {g}: matched {len(matching_files)} file(s)", file=sys.stderr)

        for aln in matching_files:
            records = list(SeqIO.parse(aln, "fasta"))

            if len(records) < 2:
                continue

            anc = str(records[0].seq)

            if debug:
                print(f"DEBUG {g}: using '{records[0].id}' as ancestor in {os.path.basename(aln)}",
                      file=sys.stderr)

            for r in records[1:]:
                iso = str(r.seq)

                if m == "missense":
                    # Single unified tolerance search — find the offset where
                    # the isolate differs from the ancestor, then use that same
                    # offset for both any_missense and exact_missense checks.
                    # This ensures both columns are consistent when a leading !
                    # artefact shifts all positions by 1.
                    best_res = None
                    found_any = False
                    for offset in range(-POSITION_TOLERANCE, POSITION_TOLERANCE + 1):
                        os_ = max(1, s + offset)
                        oe  = max(1, e + offset)
                        c0, c1 = get_alignment_columns(anc, os_, oe)
                        if c0 is None:
                            continue
                        anc_w = anc[c0:c1].replace("-", "").replace("!", "")
                        iso_w = iso[c0:c1].replace("-", "").replace("!", "")
                        if not iso_w:
                            continue
                        if anc_w != iso_w:
                            found_any = True
                            best_res = iso_w
                            break

                    if found_any:
                        present_any[(g, s, e, m)] = 1
                        if x is None or best_res == x:
                            present_exact[(g, s, e, m)] = 1

                    if debug:
                        anc_res = get_isolate_residue(anc, anc, s, e)
                        print(
                            f"DEBUG {g}:{s} WT='{anc_res}' iso='{best_res}' expected='{x}' "
                            f"any={found_any} exact={present_exact.get((g,s,e,m),0)}",
                            file=sys.stderr
                        )

                elif m in ("inframe_deletion", "delins"):
                    if confirm_indel(anc, iso, s, e, x):
                        # Indel confirmation stored in present_any for output
                        # as macse_indel_confirmed — separate from missense
                        present_any[(g, s, e, m)] = 1

    # ============================================================
    # Output
    # exact_missense: only populated for missense mutations
    # any_missense:   only populated for missense mutations
    # macse_indel_confirmed: only populated for inframe_deletion/delins
    # ============================================================

    print("gene\taa_start\taa_end\tmutation_type\texact_missense\tany_missense\tmacse_indel_confirmed")

    for g, s, e, m, _ in MUTATIONS:
        exact_missense        = present_exact.get((g,s,e,m), 0) if m == "missense" else 0
        any_missense          = present_any.get((g,s,e,m), 0)   if m == "missense" else 0
        macse_indel_confirmed = present_any.get((g,s,e,m), 0)   if m in ("inframe_deletion", "delins") else 0
        print(f"{g}\t{s}\t{e}\t{m}\t{exact_missense}\t{any_missense}\t{macse_indel_confirmed}")


# ============================================================

if __name__ == "__main__":
    main()
