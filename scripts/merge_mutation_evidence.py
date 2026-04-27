#!/usr/bin/env python3

import sys
from collections import defaultdict

def read_table(path):
    data = set()
    try:
        with open(path) as f:
            for line in f:
                if line.startswith("gene") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                data.add(parts[0])
    except FileNotFoundError:
        pass
    return data


def main(macse, frameshift, stop):
    macse_genes = read_table(macse)
    fs_genes = read_table(frameshift)
    stop_genes = read_table(stop)

    all_genes = macse_genes | fs_genes | stop_genes

    print("gene\tmacse\tframeshift\tstop\tfinal")

    for g in sorted(all_genes):
        m = int(g in macse_genes)
        f = int(g in fs_genes)
        s = int(g in stop_genes)
        final = int(m or f or s)

        print(f"{g}\t{m}\t{f}\t{s}\t{final}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: merge_mutation_evidence.py <macse> <frameshift> <stop>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
