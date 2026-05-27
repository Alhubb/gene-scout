#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

FRAMESHIFT_CHAR = "!"


def scan_alignment(aln_file):

    records = list(SeqIO.parse(aln_file, "fasta"))
    if len(records) < 2:
        return []

    # FIX: derive gene name by stripping '_aligned_nt.fasta' suffix precisely
    gene = os.path.basename(aln_file).replace("_aligned_nt.fasta", "")

    wt = str(records[0].seq)
    mutations = []

    for rec in records[1:]:
        iso = str(rec.seq)
        if len(iso) != len(wt):
            # Alignment length mismatch — skip this record
            continue

        codon_pos = 0  # 1-based WT codon counter

        for i in range(0, len(wt) - 2, 3):
            wt_codon  = wt[i:i+3]
            iso_codon = iso[i:i+3]

            # FIX: skip codons containing MACSE frameshift characters — these
            # are frameshifts, not missense. Original treated '!' as a normal
            # nucleotide and reported them as missense substitutions.
            if FRAMESHIFT_CHAR in wt_codon or FRAMESHIFT_CHAR in iso_codon:
                # Still advance position if WT codon has no gaps
                if "-" not in wt_codon:
                    codon_pos += 1
                continue

            # Advance WT codon position only for non-gap WT codons
            if "-" not in wt_codon:
                codon_pos += 1

            if wt_codon == iso_codon:
                continue

            wt_has_gap  = "-" in wt_codon
            iso_has_gap = "-" in iso_codon

            if wt_has_gap or iso_has_gap:
                # In-frame indel — both positions reported as the same codon
                mutations.append((gene, codon_pos, codon_pos, "inframe_indel"))
            else:
                mutations.append((gene, codon_pos, codon_pos, "missense"))

    return mutations


def main(macse_dir):
    print("gene\tstart\tend\ttype")

    for fname in sorted(os.listdir(macse_dir)):
        
        if not fname.endswith("_aligned_nt.fasta"):
            continue

        aln_file = os.path.join(macse_dir, fname)
        for m in scan_alignment(aln_file):
            print("\t".join(map(str, m)))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: scan_mutations_cds_macse.py <macse_output_dir>")
    main(sys.argv[1])
