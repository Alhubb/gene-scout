#!/usr/bin/env python3
"""
convert_snippy_to_mutation_list.py  <snippy.tsv>  <output.tsv>

Converts a Snippy variant table into the mutation list format expected by
screen_mutations.py:
    gene  aa_start  aa_end  mutation_type

Mutation type mapping from Snippy EFFECT strings:
    *NNNfs / ins / del producing frameshift → frameshift
    NNN_NNNdel / NNNdelinsNNN              → inframe_deletion / delins
    NNNXaaYaa (single AA change)           → missense
    NNN* (stop codon gained)               → stop_gained
    Ter/Stop extended (stop codon lost)    → stop_lost

AA positions are parsed from the EFFECT field rather than AA_POS because
AA_POS only gives the start position, while range mutations (e.g.
Trp256_Leu268del) encode both start and end in the EFFECT string.
"""

import sys
import re
import csv


# ---------------------------------------------------------------------------
# Mutation type inference from Snippy EFFECT strings
# ---------------------------------------------------------------------------

def parse_effect(effect, aa_pos_field):
    """
    Return (aa_start, aa_end, mutation_type) parsed from the Snippy EFFECT
    string, falling back to aa_pos_field for the start position if needed.

    Snippy EFFECT examples:
        Trp256_Leu268del        → inframe_deletion,  start=256, end=268
        Trp182_His187delinsCys  → delins,             start=182, end=187
        Ala164Val               → missense,           start=164, end=164
        Gln7*                   → stop_gained,        start=7,   end=7
        Ter440Cysext*?          → stop_lost,          start=440, end=440
        Gly500fs                → frameshift,         start=500, end=500
        Ala423fs                → frameshift,         start=423, end=423
        Ala41_Leu43del          → inframe_deletion,   start=41,  end=43
    """

    # --- Frameshift: ends in 'fs' ---
    m = re.match(r'[A-Za-z*]+(\d+)fs', effect)
    if m:
        pos = int(m.group(1))
        return pos, pos, "frameshift"

    # --- Stop lost: starts with Ter/Stop or contains 'ext' ---
    if effect.startswith("Ter") or effect.startswith("Stop") or "ext" in effect:
        # Extract the number from the effect string
        m = re.search(r'(\d+)', effect)
        pos = int(m.group(1)) if m else _aa_pos_start(aa_pos_field)
        return pos, pos, "stop_lost"

    # --- Stop gained: ends with * (but not ext*?) ---
    m = re.match(r'[A-Za-z]+(\d+)\*$', effect)
    if m:
        pos = int(m.group(1))
        return pos, pos, "stop_gained"

    # --- Range delins: NNN_NNNdelinsNNN ---
    m = re.match(r'[A-Za-z*]+(\d+)_[A-Za-z*]+(\d+)delins', effect)
    if m:
        return int(m.group(1)), int(m.group(2)), "delins"

    # --- Range deletion: NNN_NNNdel ---
    m = re.match(r'[A-Za-z*]+(\d+)_[A-Za-z*]+(\d+)del$', effect)
    if m:
        return int(m.group(1)), int(m.group(2)), "inframe_deletion"

    # --- Single-position deletion: NNNdel ---
    m = re.match(r'[A-Za-z*]+(\d+)del$', effect)
    if m:
        pos = int(m.group(1))
        return pos, pos, "inframe_deletion"

    # --- Missense: AlaXXXVal style (two amino acids around a number) ---
    m = re.match(r'[A-Za-z]+(\d+)[A-Za-z]+$', effect)
    if m:
        pos = int(m.group(1))
        return pos, pos, "missense"

    # --- Fallback: use AA_POS start, unknown type ---
    pos = _aa_pos_start(aa_pos_field)
    print(
        f"WARNING: could not parse EFFECT '{effect}', "
        f"using AA_POS={aa_pos_field}, type=unknown",
        file=sys.stderr,
    )
    return pos, pos, "unknown"


def _aa_pos_start(aa_pos_field):
    """Extract the start position from Snippy AA_POS field (e.g. '256/452' → 256)."""
    try:
        return int(aa_pos_field.split("/")[0])
    except (ValueError, IndexError):
        return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: convert_snippy_to_mutation_list.py <snippy.tsv> <output.tsv>")

    in_path  = sys.argv[1]
    out_path = sys.argv[2]

    written = 0
    skipped = 0

    with open(in_path) as fh_in, open(out_path, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")

        # Validate required columns are present
        required = {"GENE", "AA_POS", "EFFECT", "FTYPE"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: input TSV is missing columns: {missing}")

        writer = csv.writer(fh_out, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "aa_start", "aa_end", "mutation_type"])

        for lineno, row in enumerate(reader, start=2):
            # Only process CDS rows
            if row.get("FTYPE", "").strip() != "CDS":
                skipped += 1
                continue

            gene   = row["GENE"].strip()
            effect = row["EFFECT"].strip()
            aa_pos = row["AA_POS"].strip()

            if not gene or not effect:
                print(
                    f"WARNING: line {lineno} missing GENE or EFFECT, skipping",
                    file=sys.stderr,
                )
                skipped += 1
                continue

            aa_start, aa_end, mut_type = parse_effect(effect, aa_pos)

            if mut_type == "unknown":
                skipped += 1
                continue

            writer.writerow([gene, aa_start, aa_end, mut_type])
            written += 1

    print(f"✅ Done — {written} mutations written to {out_path}")
    if skipped:
        print(f"⚠  {skipped} rows skipped (non-CDS or unparseable)")


if __name__ == "__main__":
    main()
