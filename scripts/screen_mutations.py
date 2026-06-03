#!/usr/bin/env python3
"""
screen_mutations.py

Screens a predefined mutation list (from Snippy or similar) against the
evidence files produced by the pipeline, and reports per-mutation
confirmation status.

Fixes vs previous version:
  1. Three-tier frameshift detection: exact_frameshift (exact codon match),
     regional_frameshift (within ±REGIONAL_WINDOW codons), gene_frameshift
     (anywhere in gene). Populated from detect_frameshifts_macse.py output.
  2. macse_any_missense retained for missense mutations.
  3. final_confirmed includes all three frameshift tiers.
"""

import argparse
import csv
import sys


def read_macse_hits(macse_file):
    """
    Return two sets of (gene, aa_start, aa_end, mutation_type) keys:
      - exact_hits:  present=1 (exact residue match for missense, or indel confirmed)
      - any_hits:    any_missense=1 (any amino acid change at that position)
    """
    exact_hits = set()
    any_hits   = set()
    try:
        with open(macse_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    key = (
                        row["gene"],
                        int(row["aa_start"]),
                        int(row["aa_end"]),
                        row["mutation_type"],
                    )
                    if int(row.get("present", 0)):
                        exact_hits.add(key)
                    if int(row.get("any_missense", 0)):
                        any_hits.add(key)
                except (KeyError, ValueError):
                    continue
    except FileNotFoundError:
        print(f"WARNING: {macse_file} not found", file=sys.stderr)
    return exact_hits, any_hits


def read_frameshift_hits(path):
    """
    Read frameshifts.tsv — four columns: gene  aa_position  tier  isolate_id.
    Returns three sets of (gene, aa_position) tuples, one per tier:
      exact_hits    — exact codon match
      regional_hits — within ±REGIONAL_WINDOW codons
      gene_hits     — anywhere in gene (gene stored as (gene, 0))
    Written by detect_frameshifts_macse.py.
    """
    exact_hits    = set()
    regional_hits = set()
    gene_hits     = set()

    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                try:
                    pos = int(parts[1].strip())
                except ValueError:
                    continue
                tier = parts[2].strip() if len(parts) > 2 else "gene"

                if tier == "exact":
                    exact_hits.add((gene, pos))
                elif tier == "regional":
                    regional_hits.add((gene, pos))
                else:
                    gene_hits.add(gene)

    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)

    return exact_hits, regional_hits, gene_hits


def read_indel_hits(path):
    """
    Read indels.tsv — four columns: gene  aa_position  tier  isolate_id.
    Returns three sets for tier matching:
      exact_deletion_hits — set of (gene, aa_position) tuples where aa_position
                            is the centre of the confirmed deletion window.
                            Matched against (gene, (aa_start+aa_end)//2) from
                            the mutation list to avoid cross-contamination between
                            multiple indel targets in the same gene.
      regional_hits       — set of gene names with regional indel confirmation
      gene_hits           — set of gene names with gene-level indel confirmation
    Written by detect_indels_macse.py.
    """
    exact_deletion_hits = set()
    regional_hits       = set()
    gene_hits           = set()

    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                try:
                    pos = int(parts[1].strip())
                except ValueError:
                    continue
                tier = parts[2].strip()

                if tier == "exact_deletion":
                    exact_deletion_hits.add((gene, pos))
                elif tier == "regional":
                    regional_hits.add(gene)
                else:
                    gene_hits.add(gene)
    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)

    return exact_deletion_hits, regional_hits, gene_hits


def read_stop_pairs(path):
    """
    Read stop_genes.tsv — three columns: gene_name  isolate_id  event_type.
    Returns a set of (gene_name, event_type) tuples for matching against the
    mutation list. event_type is normalised to 'stop_gained' or 'stop_lost'.
    """
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
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                gene  = parts[0].strip()
                event = _normalise.get(parts[2].strip(), parts[2].strip())
                pairs.add((gene, event))
    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)
    return pairs


def normalise_stop_type(mut_type):
    """Normalise stop mutation type strings for consistent matching."""
    return {
        "stop_gain":   "stop_gained",
        "stop_gained": "stop_gained",
        "stop_loss":   "stop_lost",
        "stop_lost":   "stop_lost",
    }.get(mut_type, mut_type)


def main(args):
    macse_exact, macse_any          = read_macse_hits(args.macse)
    fs_exact, fs_regional, fs_gene  = read_frameshift_hits(args.frameshifts)
    ind_exact_del, ind_regional, ind_gene = read_indel_hits(args.indels)
    stop_pairs                      = read_stop_pairs(args.stops)

    out = open(args.output, "w") if args.output else sys.stdout

    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow([
        "gene", "aa_start", "aa_end", "mutation_type",
        "macse_strict", "macse_any_missense",
        "exact_frameshift", "regional_frameshift", "gene_frameshift",
        "exact_deletion", "regional_indel", "gene_indel",
        "stop_confirmed", "final_confirmed",
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
                macse_key = (gene, aa_start, aa_end, mut_type)

                # MACSE: exact residue match or indel confirmed
                macse_strict = int(macse_key in macse_exact)

                # MACSE: any missense at that position (missense only)
                macse_any_missense = int(
                    mut_type == "missense" and macse_key in macse_any
                )

                # Frameshift tiers (frameshift mutations only)
                is_fs = mut_type == "frameshift"
                exact_fs    = int(is_fs and (gene, aa_start) in fs_exact)
                regional_fs = int(is_fs and not exact_fs
                                  and any(g == gene for g, _ in fs_regional))
                gene_fs     = int(is_fs and not exact_fs and not regional_fs
                                  and gene in fs_gene)

                # Indel tiers (delins and inframe_deletion only)
                is_indel = mut_type in ("delins", "inframe_deletion")
                centre   = (aa_start + aa_end) // 2
                exact_del    = int(is_indel and (gene, centre) in ind_exact_del)
                regional_ind = int(is_indel and not exact_del
                                   and gene in ind_regional)
                gene_ind     = int(is_indel and not exact_del
                                   and not regional_ind
                                   and gene in ind_gene)

                # Stop: match on (gene_name, normalised_event_type)
                stop_confirmed = int(
                    norm_type in ("stop_gained", "stop_lost")
                    and (gene, norm_type) in stop_pairs
                )

                final = int(bool(
                    macse_strict or macse_any_missense
                    or exact_fs or regional_fs or gene_fs
                    or exact_del or regional_ind or gene_ind
                    or stop_confirmed
                ))

                writer.writerow([
                    gene, aa_start, aa_end, mut_type,
                    macse_strict, macse_any_missense,
                    exact_fs, regional_fs, gene_fs,
                    exact_del, regional_ind, gene_ind,
                    stop_confirmed, final,
                ])
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Screen predefined mutations against pipeline evidence"
    )
    parser.add_argument("--mutations",   required=True, help="Mutation list TSV")
    parser.add_argument("--macse",       required=True, help="MACSE mutation summary TSV")
    parser.add_argument("--frameshifts", required=True, help="Frameshift hits TSV (detect_frameshifts_macse.py)")
    parser.add_argument("--indels",      required=True, help="Indel hits TSV (detect_indels_macse.py)")
    parser.add_argument("--stops",       required=True, help="Stop mutation TSV")
    parser.add_argument("--output",      default=None,  help="Output file (default: stdout)")
    main(parser.parse_args())
