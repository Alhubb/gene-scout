#!/usr/bin/env python3

import argparse
import csv
import sys


def read_macse_hits(macse_file):

    hits = set()
    try:
        with open(macse_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    key = (
                        row["gene"],
                        int(row["start"]),
                        int(row["end"]),
                        row["type"],
                    )
                    hits.add(key)
                except (KeyError, ValueError):
                    continue
    except FileNotFoundError:
        print(f"WARNING: {macse_file} not found", file=sys.stderr)
    return hits


def read_frameshift_genes(path):
    """Return set of gene names with frameshift evidence."""
    genes = set()
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    # Accept single-column or two-column (filename TAB gene)
                    parts = line.split("\t")
                    genes.add(parts[-1].strip())
    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)
    return genes


def read_stop_pairs(path):

    pairs = set()
    _normalise = {
        "stop_gain":   "stop_gained",
        "stop_gained": "stop_gained",
        "stop_loss":   "stop_lost",
        "stop_lost":   "stop_lost",
    }
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t", maxsplit=1)
                if len(parts) < 2:
                    continue
                event = _normalise.get(parts[1].strip(), parts[1].strip())
                pairs.add((parts[0].strip(), event))
    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)
    return pairs


def normalise_stop_type(mut_type):
    return {
        "stop_gain":   "stop_gained",
        "stop_gained": "stop_gained",
        "stop_loss":   "stop_lost",
        "stop_lost":   "stop_lost",
    }.get(mut_type, mut_type)


def main(args):
    macse_hits       = read_macse_hits(args.macse)
    frameshift_genes = read_frameshift_genes(args.frameshifts)
    stop_pairs       = read_stop_pairs(args.stops)

    out = open(args.output, "w") if args.output else sys.stdout

    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow([
        "gene", "aa_start", "aa_end", "mutation_type",
        "macse_strict", "minimap2_frameshift", "stop_confirmed", "final_confirmed",
    ])

    try:
        with open(args.mutations) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for lineno, row in enumerate(reader, start=2):
                try:
                    gene      = row["gene"]
                    aa_start  = int(row["aa_start"])
                    aa_end    = int(row["aa_end"])
                    mut_type  = row["mutation_type"]
                except (KeyError, ValueError) as e:
                    print(
                        f"WARNING: skipping malformed row {lineno} in "
                        f"{args.mutations}: {e}",
                        file=sys.stderr,
                    )
                    continue

                norm_type = normalise_stop_type(mut_type)

                # MACSE match: compare against (gene, start, end, type)
                macse_key    = (gene, aa_start, aa_end, mut_type)
                macse_strict = int(macse_key in macse_hits)

                minimap2_fs = int(
                    mut_type == "frameshift" and gene in frameshift_genes
                )

                stop_confirmed = int(
                    norm_type in ("stop_gained", "stop_lost")
                    and any(iso == gene and evt == norm_type for iso, evt in stop_pairs)
                )

                final = int(bool(macse_strict or minimap2_fs or stop_confirmed))

                writer.writerow([
                    gene, aa_start, aa_end, mut_type,
                    macse_strict, minimap2_fs, stop_confirmed, final,
                ])
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Screen predefined mutations against pipeline evidence"
    )
    parser.add_argument("--mutations",   required=True, help="Mutation list TSV (gene, aa_start, aa_end, mutation_type)")
    parser.add_argument("--macse",       required=True, help="MACSE mutation summary TSV")
    parser.add_argument("--frameshifts", required=True, help="Frameshift gene list TSV")
    parser.add_argument("--stops",       required=True, help="Stop mutation TSV (isolate_id, event_type)")
    parser.add_argument("--output",      default=None,  help="Output file (default: stdout)")
    main(parser.parse_args())
