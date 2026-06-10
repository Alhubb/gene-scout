#!/usr/bin/env python3
"""
count_mutation_types_per_gene.py

Produces a comprehensive per-gene mutation spectrum by combining two sources:

  1. Per-gene CDS FASTAs (per_gene_fastas_dedup/)
     Used for: frameshift, missense, synonymous, stop_gained, stop_lost
     These are reliably detected from raw CDS sequence length and translation.

  2. MACSE alignments (macse_output/)
     Used for: inframe_deletion, inframe_insertion, delins
     These require codon-aware alignment to detect correctly — raw CDS
     sequences from minimap2 extraction are all approximately the same
     length as the WT (alignment coordinates absorb indels), so length-based
     detection of inframe indels is not possible from the CDS FASTAs alone.
     Also used for: macse_frameshift (MACSE '!' markers, cross-check)

Strategy:
  - Classify each unique CDS variant from the dedup FASTAs into one of:
    frameshift, missense, synonymous, stop_gained, stop_lost, or unknown
  - For each gene, scan the MACSE AA alignment to count isolates with:
      gaps only        → inframe_deletion
      insertions only  → inframe_insertion
      both             → delins
  - Merge the two count tables into a single output TSV

Note on double-counting:
  MACSE classifies all sequences including those that are frameshifted or
  have stop codons. The MACSE indel counts are therefore independent of the
  CDS-based counts and represent the full set of unique variants with
  alignment-confirmed indels, regardless of other mutation types.
"""

import os
import sys
import glob
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

CATEGORIES = [
    "frameshift",
    "missense",
    "synonymous",
    "stop_gained",
    "stop_lost",
    "inframe_deletion",
    "inframe_insertion",
    "delins",
]


# ============================================================
# PASS 1 — CDS-based classification (frameshift, missense, etc.)
# ============================================================

def clean_nt(seq_str):
    """Strip gaps and MACSE frameshift markers. Does NOT trim."""
    return seq_str.upper().replace("-", "").replace("!", "")


def translate(nt_str):
    """Translate using bacterial genetic code (table 11)."""
    seq = nt_str if len(nt_str) % 3 == 0 else nt_str[:-(len(nt_str) % 3)]
    return str(Seq(seq).translate(table=11))


def classify_cds(wt_nt, wt_aa, wt_stop, iso_nt_raw):
    """
    Classify a single isolate CDS relative to WT.
    Returns one of: frameshift, stop_gained, stop_lost, missense,
                    synonymous, or unknown (fallback).
    Inframe indels are NOT reliably detected here — see MACSE pass.
    """
    iso_nt = clean_nt(iso_nt_raw)

    # Frameshift: check length BEFORE trimming
    if len(iso_nt) % 3 != 0:
        return "frameshift"

    iso_aa = translate(iso_nt)
    iso_stop = iso_aa.find("*")
    wt_stop_pos = wt_stop if wt_stop != -1 else len(wt_aa)

    if iso_stop != -1 and iso_stop < wt_stop_pos:
        return "stop_gained"

    if iso_stop == -1 or iso_stop > wt_stop_pos:
        return "stop_lost"

    if iso_aa == wt_aa:
        return "synonymous" if iso_nt != wt_nt else "synonymous"

    if iso_aa != wt_aa:
        return "missense"

    return "synonymous"


def pass1_cds(fasta_dir):
    """
    Returns dict: gene_ref → {mutation_type: count}
    for frameshift, missense, synonymous, stop_gained, stop_lost.
    """
    counts = defaultdict(lambda: defaultdict(int))

    fasta_files = sorted(f for f in os.listdir(fasta_dir) if f.endswith(".fasta"))
    if not fasta_files:
        print(f"⚠  No .fasta files found in {fasta_dir}", file=sys.stderr)
        return counts

    for fasta in fasta_files:
        gene_ref = fasta.replace(".fasta", "")
        path = os.path.join(fasta_dir, fasta)
        records = list(SeqIO.parse(path, "fasta"))
        if len(records) < 2:
            print(f"⚠  Skipping {fasta} — fewer than 2 records", file=sys.stderr)
            continue

        wt_nt  = clean_nt(str(records[0].seq))
        wt_aa  = translate(wt_nt)
        wt_stop = wt_aa.find("*")

        for rec in records[1:]:
            cat = classify_cds(wt_nt, wt_aa, wt_stop, str(rec.seq))
            counts[gene_ref][cat] += 1

    return counts


# ============================================================
# PASS 2 — MACSE AA alignment-based indel classification
# ============================================================

def get_ungapped_length(seq):
    """Count non-gap, non-! characters in an alignment row."""
    return sum(1 for aa in seq if aa not in ("-", "!"))


def classify_indel_aa(anc_seq, iso_seq):
    """
    Compare isolate to ancestor in the AA alignment.
    Returns: 'inframe_deletion', 'inframe_insertion', 'delins', or None.

    Logic:
      - Count ungapped residues in ancestor and isolate windows.
      - If isolate has gaps (fewer ungapped aa) → deletion component
      - If isolate has extra residues (longer ungapped) → insertion component
      - Both present → delins
      - Neither → no indel detected
    """
    anc_cols_with_residue = [i for i, aa in enumerate(anc_seq)
                              if aa not in ("-", "!")]
    iso_gaps_at_anc = sum(1 for i in anc_cols_with_residue
                          if i < len(iso_seq) and iso_seq[i] in ("-", "!"))

    # Extra residues: columns where isolate has residue but ancestor has gap
    anc_gap_cols = [i for i, aa in enumerate(anc_seq) if aa in ("-", "!")]
    iso_extra = sum(1 for i in anc_gap_cols
                    if i < len(iso_seq) and iso_seq[i] not in ("-", "!"))

    has_deletion  = iso_gaps_at_anc > 0
    has_insertion = iso_extra > 0

    if has_deletion and has_insertion:
        return "delins"
    elif has_deletion:
        return "inframe_deletion"
    elif has_insertion:
        return "inframe_insertion"
    return None


def pass2_macse(macse_dir):
    """
    Returns dict: gene_ref → {indel_type: count}
    for inframe_deletion, inframe_insertion, delins.
    gene_ref is derived from the alignment filename (e.g. glpT_34).
    """
    counts = defaultdict(lambda: defaultdict(int))

    aa_files = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))
    if not aa_files:
        print(f"⚠  No AA alignment files found in {macse_dir}", file=sys.stderr)
        return counts

    for aln_path in sorted(aa_files):
        base = os.path.basename(aln_path)
        # e.g. glpT_34_aligned_aa.fasta → glpT_34
        gene_ref = base.replace("_aligned_aa.fasta", "")

        records = list(SeqIO.parse(aln_path, "fasta"))
        if len(records) < 2:
            continue

        anc_seq = str(records[0].seq)

        for rec in records[1:]:
            iso_seq = str(rec.seq)
            indel_type = classify_indel_aa(anc_seq, iso_seq)
            if indel_type:
                counts[gene_ref][indel_type] += 1

    return counts


# ============================================================
# MAIN — merge and write
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Count mutation types per gene from CDS FASTAs and MACSE alignments"
    )
    parser.add_argument("fasta_dir", help="per_gene_fastas_dedup/ directory")
    parser.add_argument("macse_dir", help="macse_output/ directory")
    parser.add_argument("output",    help="Output TSV file path")
    args = parser.parse_args()

    print("Pass 1: classifying CDS sequences...", file=sys.stderr)
    cds_counts = pass1_cds(args.fasta_dir)

    print("Pass 2: scanning MACSE AA alignments for indels...", file=sys.stderr)
    macse_counts = pass2_macse(args.macse_dir)

    all_gene_refs = sorted(set(cds_counts.keys()) | set(macse_counts.keys()))

    with open(args.output, "w") as out:
        out.write("\t".join(["gene_ref"] + CATEGORIES) + "\n")
        for gene_ref in all_gene_refs:
            row = [gene_ref]
            for cat in CATEGORIES:
                if cat in ("inframe_deletion", "inframe_insertion", "delins"):
                    row.append(str(macse_counts[gene_ref].get(cat, 0)))
                else:
                    row.append(str(cds_counts[gene_ref].get(cat, 0)))
            out.write("\t".join(row) + "\n")

    print(f"✅ {args.output} written ({len(all_gene_refs)} gene/lineage combinations)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
