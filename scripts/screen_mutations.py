#!/usr/bin/env python3

import argparse
import csv

def read_macse(macse_file):
    """Read MACSE strict mutation evidence."""
    macse_hits = set()
    with open(macse_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = (
                row["gene"],
                int(row["aa_start"]),
                int(row["aa_end"]),
                row["mutation_type"]
            )
            macse_hits.add(key)
    return macse_hits


def read_simple_list(path):
    """Read single-column gene list (frameshifts or stops)."""
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                genes.add(line)
    return genes


def main(args):
    macse_hits = read_macse(args.macse)
    frameshift_genes = read_simple_list(args.frameshifts)
    stop_genes = read_simple_list(args.stops)

    print(
        "gene\taa_start\taa_end\tmutation_type\t"
        "macse_strict\tminimap2_frameshift\tstop_confirmed\tfinal_confirmed"
    )

    with open(args.mutations) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene"]
            aa_start = int(row["aa_start"])
            aa_end = int(row["aa_end"])
            mut_type = row["mutation_type"]

            key = (gene, aa_start, aa_end, mut_type)

            macse_strict = 1 if key in macse_hits else 0
            minimap2_frameshift = 1 if (
                mut_type == "frameshift" and gene in frameshift_genes
            ) else 0
            stop_confirmed = 1 if (
                mut_type in ("stop_gain", "stop_loss") and gene in stop_genes
            ) else 0

            final_confirmed = 1 if (
                macse_strict or minimap2_frameshift or stop_confirmed
            ) else 0

            print(
                f"{gene}\t{aa_start}\t{aa_end}\t{mut_type}\t"
                f"{macse_strict}\t{minimap2_frameshift}\t"
                f"{stop_confirmed}\t{final_confirmed}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Screen predefined mutations against detected mutation evidence"
    )
    parser.add_argument("--mutations", required=True, help="Mutation list TSV")
    parser.add_argument("--macse", required=True, help="MACSE mutation summary TSV")
    parser.add_argument("--frameshifts", required=True, help="Frameshift gene list")
    parser.add_argument("--stops", required=True, help="Stop mutation gene list")

    main(parser.parse_args())
