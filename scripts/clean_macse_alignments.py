#!/usr/bin/env python3
"""
clean_macse_alignments.py  <macse_output_dir>  <cleaned_output_dir>

Cleans MACSE-aligned AA and NT FASTAs by removing isolate sequences that
are uninformative after alignment. Two failure modes are caught:

  1. Near-entirely-gap sequences — sequences that MACSE rendered as mostly
     gaps despite passing the pre-alignment nucleotide quality filter in
     clean_gene_fastas.py. A sequence can be long enough in nucleotides to
     pass that filter but MACSE may still collapse it into gaps after
     codon-aware alignment, contributing nothing to downstream detection.
     Threshold: ungapped AA residues must be ≥ MIN_FRACTION of ancestor length.

  2. Stop-lost sequences that are too short to reach the WT stop position —
     these sequences simply end before position of the WT stop codon due to
     truncated assemblies. They cannot confirm stop loss and would be
     misleading if included. The WT stop position is inferred from the
     ancestor AA sequence for genes whose name appears in STOP_LOST_GENES.

What is NOT removed:
  - The ancestor (first) sequence — always kept.
  - Sequences with early stop codons or frameshifts — informative and
    handled by detect_stop_mutations.py and detect_frameshifts_macse.py.
  - Sequences with in-frame deletions (gaps at specific positions) — these
    are exactly what detect_indels_macse.py is looking for.

Both the AA and NT alignment files are cleaned in lockstep so that record
sets remain consistent between the two files.

Usage:
    python clean_macse_alignments.py <macse_output_dir>  <cleaned_output_dir> \\
        [--stop-lost-genes gene1,gene2,...]

Output:
    One cleaned AA FASTA and one cleaned NT FASTA per gene in
    <cleaned_output_dir>, plus a summary to stderr.
"""

import os
import sys
import glob
import argparse

# ============================================================
# Configuration
# ============================================================

# Minimum fraction of ancestor ungapped residues an isolate must retain
# in the AA alignment to be kept.
MIN_FRACTION = 0.5


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


def ungapped_length(seq):
    """Count non-gap, non-! characters in an aligned sequence."""
    return sum(1 for aa in seq if aa not in ("-", "!"))


def wt_stop_position(anc_seq):
    """
    Return the 1-based position of the first stop codon in the ancestor
    AA sequence, or None if no stop codon is present.
    """
    ungapped = [aa for aa in anc_seq if aa not in ("-", "!")]
    return next((i + 1 for i, aa in enumerate(ungapped) if aa == "*"), None)


def records_by_id(records):
    """Return dict of header_id → (header, seq) for fast lookup."""
    return {h.split()[0]: (h, s) for h, s in records}


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Clean MACSE-aligned FASTAs (post-alignment gap filter)"
    )
    parser.add_argument("macse_output_dir",   help="Directory of MACSE aligned FASTAs")
    parser.add_argument("cleaned_output_dir", help="Output directory for cleaned FASTAs")
    parser.add_argument(
        "--stop-lost-genes",
        default="",
        help="Comma-separated gene names with stop_lost mutations. "
             "Isolates shorter than the WT stop position are removed. "
             "Example: --stop-lost-genes uhpC,rpoB"
    )
    args = parser.parse_args()

    in_dir  = args.macse_output_dir
    out_dir = args.cleaned_output_dir
    os.makedirs(out_dir, exist_ok=True)

    stop_lost_genes = {g.strip() for g in args.stop_lost_genes.split(",") if g.strip()}
    if stop_lost_genes:
        print(f"Stop-lost length filter active for: {', '.join(sorted(stop_lost_genes))}",
              file=sys.stderr)

    aa_files = sorted(glob.glob(os.path.join(in_dir, "*_aligned_aa.fasta")))
    if not aa_files:
        sys.exit(f"ERROR: no *_aligned_aa.fasta files found in {in_dir}")

    total_in = total_kept = total_removed = 0

    for aa_path in aa_files:
        base     = os.path.basename(aa_path)
        gene_ref = base.replace("_aligned_aa.fasta", "")
        # Extract gene name from gene_ref (e.g. uhpC_34 → uhpC)
        gene_name = gene_ref.rsplit("_", 1)[0]

        nt_path = os.path.join(in_dir, f"{gene_ref}_aligned_nt.fasta")
        if not os.path.exists(nt_path):
            print(f"  WARNING: no NT file for {gene_ref}, skipping", file=sys.stderr)
            continue

        aa_records = parse_fasta(aa_path)
        nt_records = parse_fasta(nt_path)

        if len(aa_records) < 2:
            write_fasta(aa_records, os.path.join(out_dir, base))
            write_fasta(nt_records, os.path.join(out_dir, f"{gene_ref}_aligned_nt.fasta"))
            continue

        anc_aa_header, anc_aa_seq = aa_records[0]
        anc_len  = ungapped_length(anc_aa_seq)
        wt_stop  = wt_stop_position(anc_aa_seq) if gene_name in stop_lost_genes else None

        if anc_len == 0:
            print(f"  WARNING: ancestor in {gene_ref} has zero ungapped residues, skipping",
                  file=sys.stderr)
            continue

        nt_by_id = records_by_id(nt_records)

        kept_aa  = []
        kept_nt  = []
        removed  = []

        for aa_header, aa_seq in aa_records[1:]:
            iso_len = ungapped_length(aa_seq)
            frac    = iso_len / anc_len

            # Filter 1: near-entirely-gap
            if frac < MIN_FRACTION:
                removed.append((aa_header, f"gap_fraction {frac:.1%} of ancestor"))
                continue

            # Filter 2: stop_lost gene — sequence too short to reach WT stop
            if wt_stop is not None and iso_len < wt_stop:
                removed.append((aa_header,
                                 f"too_short_for_stop_lost "
                                 f"({iso_len} aa < wt_stop pos {wt_stop})"))
                continue

            kept_aa.append((aa_header, aa_seq))
            iso_id = aa_header.split()[0]
            if iso_id in nt_by_id:
                kept_nt.append(nt_by_id[iso_id])
            else:
                # Fallback: match by ID field only (first whitespace-delimited
                # token), never by substring of the full header string —
                # short IDs like "10" would match unrelated records otherwise.
                matched = [(h, s) for h, s in nt_records[1:]
                           if h.split()[0] == iso_id]
                if matched:
                    kept_nt.append(matched[0])
                else:
                    print(f"  WARNING: no NT record found for {iso_id} in {gene_ref}",
                          file=sys.stderr)

        n_in      = len(aa_records) - 1
        n_kept    = len(kept_aa)
        n_removed = len(removed)

        total_in      += n_in
        total_kept    += n_kept
        total_removed += n_removed

        aa_out = os.path.join(out_dir, base)

        # Strip single leading ! MACSE artefacts.
        # MACSE occasionally inserts a spurious ! at alignment position 1 for
        # isolates that are otherwise intact — caused by global optimisation
        # rather than a genuine frameshift. This shifts all downstream AA
        # positions by 1, causing missense detection to look at the wrong
        # residue. Only a single leading ! is stripped; multiple consecutive
        # ! at the start indicate a genuine frameshift and are left intact.
        def strip_leading_bang(seq):
            if seq.startswith("!") and not seq.startswith("!!"):
                return seq[1:]
            return seq

        cleaned_aa = [(h, strip_leading_bang(s)) for h, s in kept_aa]
        cleaned_nt = [(h, strip_leading_bang(s)) for h, s in kept_nt]

        write_fasta([(anc_aa_header, anc_aa_seq)] + cleaned_aa, aa_out)

        anc_nt_header, anc_nt_seq = nt_records[0]
        nt_out = os.path.join(out_dir, f"{gene_ref}_aligned_nt.fasta")
        write_fasta([(anc_nt_header, anc_nt_seq)] + cleaned_nt, nt_out)

        if n_removed > 0:
            print(f"  {gene_ref}: {n_in} → {n_kept} (removed {n_removed})",
                  file=sys.stderr)
            for h, reason in removed[:3]:
                print(f"    ✗ {h[:70]}: {reason}", file=sys.stderr)
            if n_removed > 3:
                print(f"    ... and {n_removed - 3} more", file=sys.stderr)
        else:
            print(f"  {gene_ref}: {n_in} → {n_kept} ✅", file=sys.stderr)

    print(
        f"\n✅ MACSE alignment cleaning complete: "
        f"{total_in} sequences in → "
        f"{total_kept} kept, "
        f"{total_removed} removed",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
