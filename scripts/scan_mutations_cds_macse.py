#!/usr/bin/env python3

import sys
import os
from Bio import AlignIO

def scan_alignment(aln_file):
    aln = AlignIO.read(aln_file, "fasta")
    wt = aln[0]  # WT is first by construction

    gene = os.path.basename(aln_file).split("_aligned")[0]

    mutations = []

    for rec in aln[1:]:
        pos = 0
        for i in range(0, len(wt.seq), 3):
            wt_codon = str(wt.seq[i:i+3])
            mut_codon = str(rec.seq[i:i+3])

            pos += 1

            if "-" in wt_codon or "-" in mut_codon:
                # in-frame gap
                if wt_codon != mut_codon:
                    mutations.append(
                        (gene, pos, pos, "inframe_deletion")
                    )
                continue

            if wt_codon != mut_codon:
                mutations.append(
                    (gene, pos, pos, "missense")
                )

    return mutations


def main(macse_dir):
    print("gene\tstart\tend\ttype", file=sys.stdout)

    for f in os.listdir(macse_dir):
        if not f.endswith("_aligned_nt.fasta"):
            continue

        aln_file = os.path.join(macse_dir, f)
        muts = scan_alignment(aln_file)

        for m in muts:
            print("\t".join(map(str, m)))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: scan_mutations_cds_macse.py <macse_output_dir>")
    main(sys.argv[1])
