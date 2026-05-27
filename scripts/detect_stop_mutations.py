#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq


def translate_nt(nt_str):

    seq = nt_str.upper().replace("-", "").replace("!", "")
    remainder = len(seq) % 3
    if remainder:
        seq = seq[:-remainder]
    return str(Seq(seq).translate(table=11))


def main(nt_fasta):
    records = list(SeqIO.parse(nt_fasta, "fasta"))
    if len(records) < 2:
        return

    wt_aa = translate_nt(str(records[0].seq))
    wt_stop_pos = wt_aa.find("*")
    if wt_stop_pos == -1:
        wt_stop_pos = len(wt_aa)

    for rec in records[1:]:
        aa = translate_nt(str(rec.seq))
        stop_pos = aa.find("*")

        if stop_pos != -1 and stop_pos < wt_stop_pos:
            # FIX: print rec.id (the isolate) not wt.id
            print(f"{rec.id}\tstop_gained")
            continue

        if stop_pos == -1 or stop_pos > wt_stop_pos:
            print(f"{rec.id}\tstop_lost")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: detect_stop_mutations.py <nucleotide_cds_fasta>")
    main(sys.argv[1])
