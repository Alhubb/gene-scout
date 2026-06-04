#!/usr/bin/env python3
"""
convert_snippy_to_mutation_list.py  <snippy.tsv>  <output.tsv>

Converts a Snippy variant table into the mutation list format expected by
scan_mutations_cds_macse.py and screen_mutations.py:
    gene  aa_start  aa_end  mutation_type  expected

The 'expected' column contains:
  - The single-letter mutant amino acid for missense (e.g. Ala164Val → V)
  - The single-letter inserted residue for delins (e.g. Trp182delinsCys → C)
  - Empty for all other mutation types (frameshift, stop, inframe_deletion)
"""

import sys
import re
import csv

# ---------------------------------------------------------------------------
# Three-letter to single-letter amino acid conversion
# ---------------------------------------------------------------------------

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*",
}


def three_to_one(aa3):
    """Convert three-letter amino acid code to single-letter. Returns None if unknown."""
    return THREE_TO_ONE.get(aa3.capitalize(), None)


# ---------------------------------------------------------------------------
# Mutation type and expected AA inference from Snippy EFFECT strings
# ---------------------------------------------------------------------------

def parse_effect(effect, aa_pos_field):
    """
    Return (aa_start, aa_end, mutation_type, expected_aa) parsed from the
    Snippy EFFECT string. expected_aa is the single-letter amino acid expected
    in isolates (mutant residue for missense, inserted residue for delins,
    None for all other types).

    Snippy EFFECT examples:
        Trp256_Leu268del        → inframe_deletion,  start=256, end=268, expected=None
        Trp182_His187delinsCys  → delins,             start=182, end=187, expected='C'
        Ala164Val               → missense,           start=164, end=164, expected='V'
        Gln7*                   → stop_gained,        start=7,   end=7,   expected=None
        Ter440Cysext*?          → stop_lost,          start=440, end=440, expected=None
        Gly500fs                → frameshift,         start=500, end=500, expected=None
        Ala41_Leu43del          → inframe_deletion,   start=41,  end=43,  expected=None
    """

    # --- Frameshift: ends in 'fs' ---
    m = re.match(r'[A-Za-z*]+(\d+)fs', effect)
    if m:
        return int(m.group(1)), int(m.group(1)), "frameshift", None

    # --- Stop lost: starts with Ter/Stop or contains 'ext' ---
    if effect.startswith("Ter") or effect.startswith("Stop") or "ext" in effect:
        m = re.search(r'(\d+)', effect)
        pos = int(m.group(1)) if m else _aa_pos_start(aa_pos_field)
        return pos, pos, "stop_lost", None

    # --- Stop gained: ends with * ---
    m = re.match(r'[A-Za-z]+(\d+)\*$', effect)
    if m:
        return int(m.group(1)), int(m.group(1)), "stop_gained", None

    # --- Range delins: NNN_NNNdelinsNNN — extract inserted residue ---
    m = re.match(r'[A-Za-z*]+(\d+)_[A-Za-z*]+(\d+)delins([A-Za-z]+)', effect)
    if m:
        inserted = three_to_one(m.group(3)[:3]) if len(m.group(3)) >= 3 else m.group(3)[0].upper()
        return int(m.group(1)), int(m.group(2)), "delins", inserted

    # --- Range deletion: NNN_NNNdel ---
    m = re.match(r'[A-Za-z*]+(\d+)_[A-Za-z*]+(\d+)del$', effect)
    if m:
        return int(m.group(1)), int(m.group(2)), "inframe_deletion", None

    # --- Single-position deletion: NNNdel ---
    m = re.match(r'[A-Za-z*]+(\d+)del$', effect)
    if m:
        return int(m.group(1)), int(m.group(1)), "inframe_deletion", None

    # --- Missense: AlaXXXVal — extract mutant (destination) residue ---
    m = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3})$', effect)
    if m:
        dest = three_to_one(m.group(3))
        return int(m.group(2)), int(m.group(2)), "missense", dest

    # --- Fallback ---
    pos = _aa_pos_start(aa_pos_field)
    print(
        f"WARNING: could not parse EFFECT '{effect}', "
        f"using AA_POS={aa_pos_field}, type=unknown",
        file=sys.stderr,
    )
    return pos, pos, "unknown", None


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

        required = {"GENE", "AA_POS", "EFFECT", "FTYPE"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: input TSV is missing columns: {missing}")

        writer = csv.writer(fh_out, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "aa_start", "aa_end", "mutation_type", "expected"])

        for lineno, row in enumerate(reader, start=2):
            if row.get("FTYPE", "").strip() != "CDS":
                skipped += 1
                continue

            gene   = row["GENE"].strip()
            effect = row["EFFECT"].strip()
            aa_pos = row["AA_POS"].strip()

            if not gene or not effect:
                print(f"WARNING: line {lineno} missing GENE or EFFECT, skipping",
                      file=sys.stderr)
                skipped += 1
                continue

            aa_start, aa_end, mut_type, expected = parse_effect(effect, aa_pos)

            if mut_type == "unknown":
                skipped += 1
                continue

            writer.writerow([gene, aa_start, aa_end, mut_type,
                             expected if expected is not None else ""])
            written += 1

    print(f"✅ Done — {written} mutations written to {out_path}")
    if skipped:
        print(f"⚠  {skipped} rows skipped (non-CDS or unparseable)")


if __name__ == "__main__":
    main()
