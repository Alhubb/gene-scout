#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main(faa_file):
    records = list(SeqIO.parse(faa_file, "fasta"))
    if len(records) < 2:
        return

    wt = records[0]
    wt_seq = str(wt.seq)

    wt_stop_pos = wt_seq.find("*")
    if wt_stop_pos == -1:
        wt_stop_pos = len(wt_seq)

    for rec in records[1:]:
        seq = str(rec.seq)
        stop_pos = seq.find("*")

        if stop_pos == -1:
            if len(seq) > wt_stop_pos:
                print(f"{wt.id}\tstop_lost")
        elif stop_pos < wt_stop_pos:
            print(f"{wt.id}\tstop_gained")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: detect_stop_mutations.py <protein_fasta>")
    main(sys.argv[1])
