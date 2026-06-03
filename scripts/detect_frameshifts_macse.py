#!/usr/bin/env python3
"""
detect_frameshifts_macse.py  <macse_output_dir>  <mutation_list.tsv>  <output.tsv>

Three-tier position-aware frameshift detection from MACSE NT alignments.

MACSE inserts '!' characters into NT alignments to mark frameshifting
insertions. For each frameshift in the mutation list, three tiers of
evidence are reported:

  Tier 1 — exact_frameshift:
    '!' at the exact codon position from the mutation list.

  Tier 2 — regional_frameshift:
    '!' within ±REGIONAL_WINDOW codons of the expected position.
    Captures frameshifts in the same functional region that may have
    the same phenotypic effect as the evolved strain mutation.

  Tier 3 — gene_frameshift:
    Any '!' anywhere in the gene. Broadest evidence of selection
    pressure on that gene in clinical isolates.

Output TSV columns:
  gene  aa_position  tier  isolate_id

Where:
  aa_position  — exact codon where '!' was found (all tiers)
  tier         — exact / regional / gene
  isolate_id   — FASTA record ID of the isolate

The output file is used by screen_mutations.py which reads the gene
and tier columns to populate exact_frameshift, regional_frameshift,
and gene_frameshift columns in mutation_screen_results.tsv.
"""

import os
import sys
import csv
import glob
from Bio import SeqIO

# Window size (in codons) either side of expected position for Tier 2
REGIONAL_WINDOW = 20


# ============================================================
# Alignment helpers
# ============================================================

def nt_to_alignment_column(anc_nt, nt_pos):
    """
    Return the alignment column index corresponding to nucleotide position
    nt_pos (1-based) in the ancestor NT sequence, skipping gaps and '!'.
    """
    count = 0
    for i, base in enumerate(anc_nt):
        if base not in ("-", "!"):
            count += 1
            if count == nt_pos:
                return i
    return None


def codon_alignment_columns(anc_nt, aa_pos):
    """
    Return (col_start, col_end) alignment column indices for the codon
    at amino acid position aa_pos (1-based) in the ancestor NT sequence.
    col_end is exclusive (suitable for slicing).
    Returns (None, None) if aa_pos is out of range.
    """
    nt_pos = (aa_pos - 1) * 3 + 1
    col_start = nt_to_alignment_column(anc_nt, nt_pos)
    if col_start is None:
        return None, None

    count = 0
    col_end = col_start
    while col_end < len(anc_nt) and count < 3:
        if anc_nt[col_end] not in ("-", "!"):
            count += 1
        col_end += 1

    return col_start, col_end


def get_codon_length(anc_nt):
    """Return the number of codons in the ancestor sequence."""
    ungapped = sum(1 for b in anc_nt if b not in ("-", "!"))
    return ungapped // 3


def frameshift_positions_in_isolate(anc_nt, iso_nt):
    """
    Return a sorted list of codon positions (1-based) where the isolate
    has a '!' character. Walks the ancestor to map alignment columns
    back to codon positions.
    """
    positions = []
    anc_nt_count = 0   # ungapped nt count in ancestor

    for col, (anc_base, iso_base) in enumerate(zip(anc_nt, iso_nt)):
        if anc_base not in ("-", "!"):
            anc_nt_count += 1
        if iso_base == "!":
            # Codon position based on ancestor NT count at this column
            codon_pos = (anc_nt_count - 1) // 3 + 1
            if not positions or positions[-1] != codon_pos:
                positions.append(codon_pos)

    return positions


# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: detect_frameshifts_macse.py "
            "<macse_output_dir> <mutation_list.tsv> <output.tsv>"
        )

    macse_dir     = sys.argv[1]
    mutation_file = sys.argv[2]
    output_file   = sys.argv[3]

    # Read frameshift targets from mutation list
    frameshift_targets = []
    try:
        with open(mutation_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("mutation_type", "").strip() == "frameshift":
                    gene = row["gene"].strip()
                    pos  = int(row["aa_start"].strip())
                    frameshift_targets.append((gene, pos))
    except FileNotFoundError:
        sys.exit(f"ERROR: mutation list not found: {mutation_file}")

    if not frameshift_targets:
        print("No frameshift mutations in mutation list — nothing to do.",
              file=sys.stderr)
        sys.exit(0)

    unique_targets = sorted(set(frameshift_targets))
    print(
        f"Searching for {len(unique_targets)} unique frameshift position(s) "
        f"(regional window: ±{REGIONAL_WINDOW} codons):",
        file=sys.stderr,
    )
    for gene, pos in unique_targets:
        print(f"  {gene} codon {pos}", file=sys.stderr)

    # results: list of (gene, found_pos, tier, isolate_id)
    results = []

    for gene, expected_pos in unique_targets:
        pattern = os.path.join(macse_dir, f"{gene}*_aligned_nt.fasta")
        aln_files = glob.glob(pattern)

        if not aln_files:
            print(f"  WARNING: no NT alignment files found for gene '{gene}'",
                  file=sys.stderr)
            continue

        tier1 = tier2 = tier3 = 0

        for aln_path in sorted(aln_files):
            records = list(SeqIO.parse(aln_path, "fasta"))
            if len(records) < 2:
                continue

            anc_nt    = str(records[0].seq)
            gene_len  = get_codon_length(anc_nt)
            region_lo = max(1, expected_pos - REGIONAL_WINDOW)
            region_hi = min(gene_len, expected_pos + REGIONAL_WINDOW)

            for rec in records[1:]:
                iso_nt = str(rec.seq)
                fs_positions = frameshift_positions_in_isolate(anc_nt, iso_nt)

                if not fs_positions:
                    continue

                # Tier 1 — exact position
                if expected_pos in fs_positions:
                    results.append((gene, expected_pos, "exact", rec.id))
                    tier1 += 1
                    continue   # already confirmed at highest tier

                # Tier 2 — regional (within ±REGIONAL_WINDOW codons)
                regional = [p for p in fs_positions if region_lo <= p <= region_hi]
                if regional:
                    # Record the closest frameshift position to expected
                    closest = min(regional, key=lambda p: abs(p - expected_pos))
                    results.append((gene, closest, "regional", rec.id))
                    tier2 += 1
                    continue

                # Tier 3 — anywhere in gene
                results.append((gene, fs_positions[0], "gene", rec.id))
                tier3 += 1

        print(
            f"  {gene} codon {expected_pos}: "
            f"exact={tier1}  regional={tier2}  gene={tier3}",
            file=sys.stderr,
        )

    # Write output
    with open(output_file, "w") as out:
        for gene, pos, tier, isolate_id in results:
            out.write(f"{gene}\t{pos}\t{tier}\t{isolate_id}\n")

    total = len(results)
    exact    = sum(1 for r in results if r[2] == "exact")
    regional = sum(1 for r in results if r[2] == "regional")
    gene_lvl = sum(1 for r in results if r[2] == "gene")
    print(
        f"\n✅ {total} frameshift hits written to {output_file}\n"
        f"   exact={exact}  regional={regional}  gene={gene_lvl}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
