#!/usr/bin/env python3

import sys
import gzip
import io
from Bio import SeqIO

VALID = set("ACGTNacgtn")


def revcomp(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def parse_best_hit(paf_file):

    best = None
    best_matches = -1

    with open(paf_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            try:
                mapq    = int(fields[11])
                matches = int(fields[9])
                strand  = fields[4]
                tname   = fields[5]
                tstart  = int(fields[7])
                tend    = int(fields[8])
            except (ValueError, IndexError):
                continue

            # Skip ambiguous / unmapped hits
            if mapq == 0:
                continue

            if matches > best_matches:
                best_matches = matches
                best = (tname, strand, tstart, tend)

    return best


def read_assembly(assembly_file):

    if assembly_file == "/dev/stdin":
        raw = sys.stdin.buffer.read()
        handle = (
            io.TextIOWrapper(gzip.open(io.BytesIO(raw)))
            if raw[:2] == b"\x1f\x8b"
            else io.StringIO(raw.decode())
        )
    elif assembly_file.endswith(".gz"):
        handle = io.TextIOWrapper(gzip.open(assembly_file, "rb"))
    else:
        handle = open(assembly_file)

    return {rec.id: str(rec.seq) for rec in SeqIO.parse(handle, "fasta")}


def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: extract_gene_from_paf.py "
            "<paf> <assembly.fna[.gz] | /dev/stdin> <output.fasta>"
        )

    paf_file      = sys.argv[1]
    assembly_file = sys.argv[2]
    output_file   = sys.argv[3]

    records = read_assembly(assembly_file)

    hit = parse_best_hit(paf_file)
    if hit is None:
        print(f"[WARN] No valid hits in {paf_file}", file=sys.stderr)
        sys.exit(0)

    tname, strand, start, end = hit
    if end < start:
        start, end = end, start

    seq = records.get(tname) or records.get(tname.split()[0])
    if seq is None:
        print(f"[WARN] Contig {tname} not found in assembly", file=sys.stderr)
        sys.exit(0)

    extracted = str(seq[start:end])
    if strand == "-":
        extracted = revcomp(extracted)

    if not extracted:
        print(f"[WARN] Empty sequence extracted from {tname}", file=sys.stderr)
        sys.exit(0)

    if any(b not in VALID for b in extracted):
        print(
            f"[WARN] Non-IUPAC characters in sequence from {tname}",
            file=sys.stderr,
        )
        sys.exit(0)

    with open(output_file, "w") as fh:
        fh.write(f">{tname}:{start}-{end}({strand})\n{extracted}\n")


if __name__ == "__main__":
    main()
