#!/usr/bin/env python3
"""
detect_indels_macse.py  <macse_output_dir>  <mutation_list.tsv>  <output.tsv>

Three-tier in-frame indel detection from MACSE amino-acid alignments,
for delins and inframe_deletion mutations.

This script provides regional and gene-level evidence for indel mutations
that were not confirmed exactly by scan_mutations_cds_macse.py. It is
intentionally parallel in structure to detect_frameshifts_macse.py.

Three tiers of evidence:

  Tier 1 — exact_deletion:
    Gaps in the isolate at the exact residue window from the mutation list.
    For inframe_deletion: confirms the exact deletion is present.
    For delins: confirms the deletion component is present but does NOT
    confirm the insertion — use macse_strict for full delins confirmation
    (deletion + specific inserted residue). This tier fires when clinical
    isolates show the same deletion as the evolved strain, even if they
    carry a different insertion or no insertion at all.

  Tier 2 — regional_indel:
    Any gap-containing window within ±REGIONAL_WINDOW residues of the
    expected position. Captures indels in the same functional region.

  Tier 3 — gene_indel:
    Any in-frame indel (gap) anywhere in the gene. Broadest evidence
    of selection pressure on that gene.

Output TSV columns (no header):
  gene  aa_position  tier  isolate_id

Where tier is one of: exact_deletion / regional / gene
"""

import os
import sys
import csv
import glob
from Bio import SeqIO

# Window (in residues) either side of expected position for Tier 2
REGIONAL_WINDOW = 20


# ============================================================
# Alignment helpers
# ============================================================

def ungapped_to_alignment_index(seq, pos):
    """Return alignment column for ungapped position pos (1-based)."""
    count = 0
    for i, aa in enumerate(seq):
        if aa not in ("-", "!"):
            count += 1
            if count == pos:
                return i
    return None


def get_alignment_columns(anc_seq, start, end):
    """Return (col_start, col_end) for ungapped positions start..end (1-based)."""
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


def get_gene_length(anc_seq):
    """Return number of ungapped residues in ancestor AA sequence."""
    return sum(1 for aa in anc_seq if aa not in ("-", "!"))


def indel_positions_in_isolate(anc_seq, iso_seq):
    """
    Return a sorted list of (aa_pos, gap_count) tuples where the isolate
    has gaps relative to the ancestor. aa_pos is the 1-based ancestor
    residue position at the start of each gap run.
    """
    positions = []
    anc_count = 0
    in_gap = False
    gap_start_pos = 0
    gap_count = 0

    for anc_aa, iso_aa in zip(anc_seq, iso_seq):
        if anc_aa not in ("-", "!"):
            anc_count += 1

        iso_is_gap = iso_aa in ("-", "!")

        if iso_is_gap and not in_gap:
            in_gap = True
            gap_start_pos = anc_count
            gap_count = 1
        elif iso_is_gap and in_gap:
            gap_count += 1
        elif not iso_is_gap and in_gap:
            positions.append((gap_start_pos, gap_count))
            in_gap = False
            gap_count = 0

    if in_gap:
        positions.append((gap_start_pos, gap_count))

    return positions


# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: detect_indels_macse.py "
            "<macse_output_dir> <mutation_list.tsv> <output.tsv>"
        )

    macse_dir     = sys.argv[1]
    mutation_file = sys.argv[2]
    output_file   = sys.argv[3]

    # Read delins and inframe_deletion targets from mutation list
    indel_targets = []
    try:
        with open(mutation_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("mutation_type", "").strip() in ("delins", "inframe_deletion"):
                    gene = row["gene"].strip()
                    start = int(row["aa_start"].strip())
                    end   = int(row["aa_end"].strip())
                    indel_targets.append((gene, start, end))
    except FileNotFoundError:
        sys.exit(f"ERROR: mutation list not found: {mutation_file}")

    if not indel_targets:
        print("No delins/inframe_deletion mutations in mutation list — nothing to do.",
              file=sys.stderr)
        sys.exit(0)

    unique_targets = sorted(set(indel_targets))
    print(
        f"Searching for {len(unique_targets)} unique indel target(s) "
        f"(regional window: ±{REGIONAL_WINDOW} residues):",
        file=sys.stderr,
    )
    for gene, start, end in unique_targets:
        print(f"  {gene} residues {start}–{end}", file=sys.stderr)

    # results: list of (gene, aa_position, tier, isolate_id)
    results = []

    for gene, expected_start, expected_end in unique_targets:
        pattern = os.path.join(macse_dir, f"{gene}*_aligned_aa.fasta")
        aln_files = glob.glob(pattern)

        if not aln_files:
            print(f"  WARNING: no AA alignment files found for gene '{gene}'",
                  file=sys.stderr)
            continue

        tier1 = tier2 = tier3 = 0

        for aln_path in sorted(aln_files):
            records = list(SeqIO.parse(aln_path, "fasta"))
            if len(records) < 2:
                continue

            anc_seq  = str(records[0].seq)
            gene_len = get_gene_length(anc_seq)

            region_lo = max(1, expected_start - REGIONAL_WINDOW)
            region_hi = min(gene_len, expected_end + REGIONAL_WINDOW)
            expected_centre = (expected_start + expected_end) // 2

            # Get alignment columns for exact window
            c0, c1 = get_alignment_columns(anc_seq, expected_start, expected_end)

            for rec in records[1:]:
                iso_seq = str(rec.seq)
                indel_pos = indel_positions_in_isolate(anc_seq, iso_seq)

                if not indel_pos:
                    continue

                # Tier 1 — exact: gap overlaps the expected window
                if c0 is not None:
                    iso_window = iso_seq[c0:c1]
                    expected_del = expected_end - expected_start + 1
                    iso_gaps = sum(1 for a in iso_window if a in ("-", "!"))
                    anc_ungapped = sum(1 for a in anc_seq[c0:c1]
                                       if a not in ("-", "!"))
                    if (anc_ungapped >= expected_del * 0.8
                            and iso_gaps >= expected_del):
                        results.append((gene, expected_centre, "exact_deletion", rec.id))
                        tier1 += 1
                        continue

                # Tier 2 — regional: any gap within ±REGIONAL_WINDOW
                regional = [(p, g) for p, g in indel_pos
                            if region_lo <= p <= region_hi]
                if regional:
                    closest_pos = min(regional,
                                      key=lambda pg: abs(pg[0] - expected_centre))
                    results.append((gene, closest_pos[0], "regional", rec.id))
                    tier2 += 1
                    continue

                # Tier 3 — gene: any gap anywhere
                results.append((gene, indel_pos[0][0], "gene", rec.id))
                tier3 += 1

        print(
            f"  {gene} {expected_start}–{expected_end}: "
            f"exact={tier1}  regional={tier2}  gene={tier3}",
            file=sys.stderr,
        )

    # Write output
    with open(output_file, "w") as out:
        for gene, pos, tier, isolate_id in results:
            out.write(f"{gene}\t{pos}\t{tier}\t{isolate_id}\n")

    total    = len(results)
    exact    = sum(1 for r in results if r[2] == "exact_deletion")
    regional = sum(1 for r in results if r[2] == "regional")
    gene_lvl = sum(1 for r in results if r[2] == "gene")
    print(
        f"\n✅ {total} indel hits written to {output_file}\n"
        f"   exact_deletion={exact}  regional={regional}  gene={gene_lvl}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
