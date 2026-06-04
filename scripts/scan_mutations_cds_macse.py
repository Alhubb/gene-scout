#!/usr/bin/env python3

import os
import sys
import glob
from Bio import SeqIO

DEBUG = True

# ============================================================
# Load mutation list
# ============================================================

def load_mutations(tsv_file):
    mutations = []

    with open(tsv_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")

            if len(fields) < 4:
                continue

            gene = fields[0]
            start = int(fields[1])
            end   = int(fields[2])
            mtype = fields[3]

            expected = None
            if len(fields) >= 5 and fields[4] not in (".", "", "NA"):
                expected = fields[4]

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

def confirm_any_missense(anc, iso, start, end):
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return False

    anc_w = anc[c0:c1].replace("-", "").replace("!", "")
    iso_w = iso[c0:c1].replace("-", "").replace("!", "")

    if not iso_w:
        return False

    return anc_w != iso_w


def get_isolate_residue(anc, iso, start, end):
    c0, c1 = get_alignment_columns(anc, start, end)
    if c0 is None:
        return "?"
    return iso[c0:c1].replace("-", "").replace("!", "")


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
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: scan_mutations_cds_macse.py "
            "<macse_output_dir> <mutation_list.tsv>"
        )

    macse_dir = sys.argv[1]
    mutation_file = sys.argv[2]

    MUTATIONS = load_mutations(mutation_file)

    if not MUTATIONS:
        print("[ERROR] No valid mutations loaded", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(MUTATIONS)} mutations", file=sys.stderr)

    # Collect ALL alignment files once
    all_alignments = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))

    print(f"[INFO] Found {len(all_alignments)} total alignment files", file=sys.stderr)

    present_exact = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}
    present_any   = {(g, s, e, m): 0 for g, s, e, m, _ in MUTATIONS}

    for g, s, e, m, x in MUTATIONS:
        if m not in MACSE_TYPES:
            continue

        # ✅ Robust gene matching
        matching_files = [
            f for f in all_alignments
            if g.lower() in os.path.basename(f).lower()
        ]

        print(f"[INFO] {g}: matched {len(matching_files)} alignment(s)", file=sys.stderr)

        for aln in matching_files:
            records = list(SeqIO.parse(aln, "fasta"))

            if len(records) < 2:
                continue

            anc = str(records[0].seq)

            for r in records[1:]:
                iso = str(r.seq)

                if m == "missense":
                    res = get_isolate_residue(anc, iso, s, e)
                    anc_res = get_isolate_residue(anc, anc, s, e)

                    # ✅ exact match logic (with fallback)
                    if x is None:
                        if confirm_any_missense(anc, iso, s, e):
                            present_exact[(g, s, e, m)] = 1
                    else:
                        if res == x:
                            present_exact[(g, s, e, m)] = 1

                    if confirm_any_missense(anc, iso, s, e):
                        present_any[(g, s, e, m)] = 1

                    if DEBUG:
                        print(
                            f"DEBUG {g} {s}: WT='{anc_res}' iso='{res}' expected='{x}'",
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

    # ============================================================
    # Output
    # ============================================================

    print("gene\taa_start\taa_end\tmutation_type\tpresent\tany_missense")

    for g, s, e, m, _ in MUTATIONS:
        print(
            f"{g}\t{s}\t{e}\t{m}\t"
            f"{present_exact.get((g,s,e,m), 0)}\t"
            f"{present_any.get((g,s,e,m), 0)}"
        )


if __name__ == "__main__":
    main()
