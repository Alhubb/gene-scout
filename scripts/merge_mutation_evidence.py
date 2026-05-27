#!/usr/bin/env python3
import sys


def read_macse_genes(path):
    genes = set()
    try:
        with open(path) as fh:
            fh.readline()  # skip header
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
                # Two-column format: paf_path TAB gene_name
                gene = parts[-1].strip()
                if gene:
                    genes.add(gene)
    except FileNotFoundError:
        print(f"WARNING: {path} not found — frameshift evidence will be absent", file=sys.stderr)
    return genes


def read_stop_genes(path):
    """
    Read stop_genes.tsv — three columns: gene_name  isolate_id  event_type.
    Returns a set of gene names that have any stop event.
    Column 0 is the gene name, written by detect_stop_mutations.py.
    """
    genes = set()
    try:
        with open(path) as fh:
            for lineno, line in enumerate(fh, start=1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    print(
                        f"WARNING: {path} line {lineno} has fewer than 3 fields, skipping",
                        file=sys.stderr,
                    )
                    continue
                genes.add(parts[0].strip())
    except FileNotFoundError:
        print(f"WARNING: {path} not found — stop evidence will be absent", file=sys.stderr)
    return genes


def main(macse_path, frameshift_path, stop_path):
    macse_genes      = read_macse_genes(macse_path)
    frameshift_genes = read_frameshift_genes(frameshift_path)
    stop_genes       = read_stop_genes(stop_path)

    all_genes = macse_genes | frameshift_genes | stop_genes

    print("gene\tmacse\tframeshift\tstop\tfinal")

    for gene in sorted(all_genes):
        m     = int(gene in macse_genes)
        f     = int(gene in frameshift_genes)
        s     = int(gene in stop_genes)
        final = int(bool(m or f or s))
        print(f"{gene}\t{m}\t{f}\t{s}\t{final}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: merge_mutation_evidence.py "
            "<macse_summary.tsv> <frameshifts.tsv> <stop_genes.tsv>"
        )
    main(sys.argv[1], sys.argv[2], sys.argv[3])
