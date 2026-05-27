#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO

FRAMESHIFT_CHAR = "!"

def scan_alignment(aln_file):
    records = list(SeqIO.parse(aln_file, "fasta"))
    if len(records) < 2:
        return []

    gene = os.path.basename(aln_file).replace("_aligned_nt.fasta", "")
    wt = str(records[0].seq)
    mutations = []

    for rec in records[1:]:
        iso = str(rec.seq)
        if len(iso) != len(wt):
            continue

        codon_pos = 0  # 1-based WT codon counter

        for i in range(0, len(wt) - 2, 3):
            wt_codon  = wt[i:i+3]
            iso_codon = iso[i:i+3]

            # Skip codons containing MACSE frameshift characters —
            # these are frameshifts, not missense.
            if FRAMESHIFT_CHAR in wt_codon or FRAMESHIFT_CHAR in iso_codon:
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
                mutations.append((gene, codon_pos, codon_pos, "inframe_indel"))
            else:
                mutations.append((gene, codon_pos, codon_pos, "missense"))

    return mutations


def main(macse_dir):
    # FIX: column names match what screen_mutations.py expects:
    # aa_start, aa_end, mutation_type — not start, end, type.
    print("gene\taa_start\taa_end\tmutation_type")

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
