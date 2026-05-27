#!/usr/bin/env python3

import sys


def read_macse_genes(path):

    genes = set()
    try:
        with open(path) as fh:
            header = fh.readline()  # skip header
            for lineno, line in enumerate(fh, start=2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    print(
                        f"WARNING: {path} line {lineno} has fewer than 4 fields, skipping",
                        file=sys.stderr,
                    )
                    continue
                genes.add(parts[0])
    except FileNotFoundError:
        print(f"WARNING: {path} not found — MACSE evidence will be absent", file=sys.stderr)
    return genes


def read_frameshift_genes(path):

    genes = set()
    try:
        with open(path) as fh:
            for lineno, line in enumerate(fh, start=1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                # Accept both single-column and two-column (filename TAB gene) formats
                gene = parts[-1].strip()
                if gene:
                    genes.add(gene)
    except FileNotFoundError:
        print(f"WARNING: {path} not found — frameshift evidence will be absent", file=sys.stderr)
    return genes


def read_stop_genes(path):

    stops = set()
    try:
        with open(path) as fh:
            for lineno, line in enumerate(fh, start=1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t", maxsplit=1)
                if len(parts) < 2:
                    print(
                        f"WARNING: {path} line {lineno} missing event type, skipping",
                        file=sys.stderr,
                    )
                    continue
                stops.add((parts[0].strip(), parts[1].strip()))
    except FileNotFoundError:
        print(f"WARNING: {path} not found — stop evidence will be absent", file=sys.stderr)
    return stops


def main(macse_path, frameshift_path, stop_path):
    macse_genes     = read_macse_genes(macse_path)
    frameshift_genes = read_frameshift_genes(frameshift_path)
    stop_pairs      = read_stop_genes(stop_path)

    # Collect all gene names that appear in any evidence source.
    # For stops we use the isolate ID as the key (matches what
    # detect_stop_mutations.py writes), so we report at isolate level.
    all_keys = set()
    all_keys.update(macse_genes)
    all_keys.update(frameshift_genes)
    all_keys.update(iso for iso, _ in stop_pairs)

    print("gene\tmacse\tframeshift\tstop\tfinal")

    for key in sorted(all_keys):
        m = int(key in macse_genes)
        f = int(key in frameshift_genes)
        # Stop: flag if this isolate/gene appears in any stop event
        s = int(any(iso == key for iso, _ in stop_pairs))
        final = int(bool(m or f or s))
        print(f"{key}\t{m}\t{f}\t{s}\t{final}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: merge_mutation_evidence.py "
            "<macse_summary.tsv> <frameshifts.tsv> <stop_genes.tsv>"
        )
    main(sys.argv[1], sys.argv[2], sys.argv[3])
