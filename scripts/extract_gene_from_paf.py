#!/usr/bin/env python3

import sys
from Bio import SeqIO

def revcomp(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def parse_best_hit(paf):
    best = None
    best_len = 0
    with open(paf) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            strand = fields[4]
            tname = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            alen = abs(tend - tstart)
            if alen > best_len:
                best_len = alen
                best = (tname, strand, tstart, tend)
    return best

def main():
    paf, _, out = sys.argv[1:]
    hit = parse_best_hit(paf)
    if hit is None:
        sys.exit(0)

    tname, strand, start, end = hit

    seq = None
    for rec in SeqIO.parse(sys.stdin, "fasta"):
        if rec.id == tname or rec.id.split()[0] == tname:
            s = str(rec.seq[start:end])
            seq = revcomp(s) if strand == "-" else s
            break

    if seq:
        with open(out, "w") as fh:
            fh.write(f">{tname}:{start}-{end}({strand})\n")
            fh.write(seq + "\n")

if __name__ == "__main__":
    main()
