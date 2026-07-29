#!/usr/bin/env python3
"""
extract_confirmed_sequences.py

Extracts aligned sequences for all confirmed mutations, allowing visual
verification of results.

Handles three evidence types:
  - MACSE AA alignments: missense, inframe_deletion, delins
  - MACSE NT alignments: frameshifts
  - Per-gene CDS FASTAs: stop_gained, stop_lost

Usage:
    python extract_confirmed_sequences.py \\
        <macse_output_dir> \\
        <mutation_cds_macse_strict_summary.tsv> \\
        <mutation_list.tsv> \\
        <output_dir> \\
        [--frameshifts frameshifts.tsv] \\
        [--stops stop_genes.tsv] \\
        [--per-gene-fastas per_gene_fastas/] \\
        [--context N] \\
        [--any]
"""

import os
import sys
import csv
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


# ============================================================
# Loaders
# ============================================================

def load_confirmed(summary_file, include_any=False):
    """
    Returns (confirmed, any_hits, strict_only) where:
      confirmed   — keys to extract (exact_missense=1 or macse_indel_confirmed=1,
                    plus any_missense=1 if include_any)
      any_hits    — keys where any_missense=1
      strict_only — keys where exact_missense=1 or macse_indel_confirmed=1
                    but any_missense=0
    """
    strict_hits = set()
    any_hits    = set()
    try:
        with open(summary_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    key = (row["gene"].strip(), int(row["aa_start"]),
                           int(row["aa_end"]), row["mutation_type"].strip())
                    mut_type = row["mutation_type"].strip()
                    # exact_missense for missense, macse_indel_confirmed for indels
                    if mut_type == "missense":
                        if int(row.get("exact_missense", 0)):
                            strict_hits.add(key)
                    elif mut_type in ("inframe_deletion", "delins"):
                        if int(row.get("macse_indel_confirmed", 0)):
                            strict_hits.add(key)
                    if int(row.get("any_missense", 0)):
                        any_hits.add(key)
                except (KeyError, ValueError):
                    continue
    except FileNotFoundError:
        sys.exit(f"ERROR: summary file not found: {summary_file}")

    confirmed = strict_hits | any_hits if include_any else set(strict_hits)
    strict_only = strict_hits - any_hits
    return confirmed, any_hits, strict_only


def load_exact_deletion_hits(screen_results_file):
    """
    Read mutation_screen_results.tsv and return a set of
    (gene, aa_start, aa_end, mutation_type) tuples where exact_deletion=1.
    These are delins/inframe_deletion mutations where the deletion component
    is confirmed but the insertion may not be present.
    """
    exact_deletion = set()
    try:
        with open(screen_results_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    if int(row.get("exact_deletion", 0)):
                        key = (row["gene"].strip(), int(row["aa_start"]),
                               int(row["aa_end"]), row["mutation_type"].strip())
                        exact_deletion.add(key)
                except (KeyError, ValueError):
                    continue
    except FileNotFoundError:
        print(f"WARNING: {screen_results_file} not found — skipping exact_deletion extraction",
              file=sys.stderr)
    return exact_deletion


def load_regional_deletion_hits(screen_results_file):
    """
    Read mutation_screen_results.tsv and return a set of
    (gene, aa_start, aa_end, mutation_type) tuples where regional_indel=1
    but exact_deletion=0. These are delins mutations where a deletion was
    confirmed in the region but the specific called boundary was not matched —
    sequences are extracted and labelled as 'regional deletion confirmed'
    to distinguish them from exact confirmations.
    """
    regional_hits = set()
    try:
        with open(screen_results_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    if (int(row.get("regional_indel", 0))
                            and not int(row.get("exact_deletion", 0))):
                        key = (row["gene"].strip(), int(row["aa_start"]),
                               int(row["aa_end"]), row["mutation_type"].strip())
                        regional_hits.add(key)
                except (KeyError, ValueError):
                    continue
    except FileNotFoundError:
        print(f"WARNING: {screen_results_file} not found — skipping regional_deletion extraction",
              file=sys.stderr)
    return regional_hits


def load_mutation_list(tsv_file):
    mutations = []
    try:
        with open(tsv_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene  = row.get("gene", "").strip()
                mtype = row.get("mutation_type", "").strip()
                if not gene or not mtype:
                    continue
                try:
                    start = int(row["aa_start"])
                    end   = int(row["aa_end"])
                except (KeyError, ValueError):
                    continue
                expected_raw = row.get("expected", "").strip()
                expected = None if expected_raw in ("", ".", "NA", "None") else expected_raw
                mutations.append((gene, start, end, mtype, expected))
    except FileNotFoundError:
        sys.exit(f"ERROR: mutation list not found: {tsv_file}")
    return mutations


def load_frameshift_hits(path):
    hits = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                gene, pos_str, tier, isolate = parts[0], parts[1], parts[2], parts[3]
                try:
                    pos = int(pos_str)
                except ValueError:
                    continue
                hits.setdefault((gene.strip(), pos), []).append(
                    (tier.strip(), isolate.strip()))
    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)
    return hits


def load_stop_hits(path):
    """
    Read stop_genes.tsv. Handles two formats:

    New format (4 columns): gene  aa_position  tier  isolate_id
      Written by the updated detect_stop_mutations.py. Supports full
      exact/regional/gene tiering.

    Old format (3 columns): gene  isolate_id  event_type
      Written by the original detect_stop_mutations.py. All hits are
      assigned gene tier with position 0 since no position info is available.
      extract_stops will match these by gene name only.

    Returns dict of gene → list of (aa_position, tier, isolate_id).
    """
    hits = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")

                if len(parts) >= 4:
                    # New format: gene  aa_position  tier  isolate_id
                    gene = parts[0].strip()
                    try:
                        pos = int(parts[1].strip())
                    except ValueError:
                        continue
                    tier    = parts[2].strip()
                    isolate = parts[3].strip()
                    hits.setdefault(gene, []).append((pos, tier, isolate))

                elif len(parts) == 3:
                    # Old format: gene  isolate_id  event_type
                    # No position — use pos=0, tier=gene so extract_stops
                    # includes it under gene-level evidence only.
                    gene    = parts[0].strip()
                    isolate = parts[1].strip()
                    hits.setdefault(gene, []).append((0, "gene", isolate))

    except FileNotFoundError:
        print(f"WARNING: {path} not found", file=sys.stderr)
    return hits


# ============================================================
# Alignment helpers
# ============================================================

def ungapped_to_alignment_index(seq, pos):
    count = 0
    for i, aa in enumerate(seq):
        if aa not in ("-", "!"):
            count += 1
            if count == pos:
                return i
    return None


def get_window_columns(anc_seq, start, end, context):
    ctx_start  = max(1, start - context)
    ctx_end_aa = end + context
    ungapped   = sum(1 for aa in anc_seq if aa not in ("-", "!"))
    ctx_end_aa = min(ctx_end_aa, ungapped)

    col_lo = ungapped_to_alignment_index(anc_seq, ctx_start) or 0
    col_hi = ungapped_to_alignment_index(anc_seq, ctx_end_aa)
    col_hi = (col_hi + 1) if col_hi is not None else len(anc_seq)
    return col_lo, col_hi


def is_confirmed_isolate(anc, iso, start, end, mutation_type, expected):
    col_start = ungapped_to_alignment_index(anc, start)
    if col_start is None:
        return False
    span, count, col_end = end - start + 1, 0, col_start
    while col_end < len(anc) and count < span:
        if anc[col_end] not in ("-", "!"):
            count += 1
        col_end += 1

    if mutation_type == "missense":
        iso_w = iso[col_start:col_end].replace("-", "").replace("!", "")
        anc_w = anc[col_start:col_end].replace("-", "").replace("!", "")
        if not iso_w:
            return False
        return (iso_w == expected) if expected else (iso_w != anc_w)

    elif mutation_type in ("inframe_deletion", "delins"):
        anc_w = anc[col_start:col_end]
        iso_w = iso[col_start:col_end]
        anc_ug  = sum(1 for a in anc_w if a not in ("-", "!"))
        iso_gaps = sum(1 for a in iso_w if a in ("-", "!"))
        exp_del  = end - start + 1
        if anc_ug < exp_del * 0.8 or iso_gaps < exp_del:
            return False
        if mutation_type == "delins" and expected:
            fs = max(0, col_start - 5)
            fe = min(len(iso), col_end + 5)
            return expected in iso[fs:fe].replace("-", "").replace("!", "")
        return True
    return False


# ============================================================
# Sequence quality filter
# ============================================================

def is_clean_sequence(seq, mutation_type, anc_seq=None, full_seq=None,
                      wt_stop_pos=None, del_window=None):
    """
    Return True if the sequence carries genuine evidence for its mutation type.

    Validation is now per-mutation-type rather than a single generic filter,
    catching the specific failure modes identified in the evidence audit:

    ALL types:
      - All-gap / near-empty sequences are rejected (checked against the
        full sequence if supplied, else the windowed excerpt).
        Threshold: >10 ungapped characters required.

    frameshift:
      - Must carry a '!' marker in the aligned NT sequence OR have a length
        (ungapped nt) not divisible by 3. All-gap sequences have neither
        and are excluded.

    inframe_deletion / delins / delins_del:
      - Must have at least one gap ('-') in the window. Sequences with '!'
        but no '-' are frameshifted without a deletion — reject. Entirely
        gap sequences carry no boundary information — reject.

    stop_gained:
      - Translated sequence must contain a stop codon upstream of the WT
        stop position (wt_stop_pos, 1-based aa position).
      - Sequences with >3 stop codons are flagged as likely mis-assemblies.

    stop_lost:
      - Sequence must extend to at least the WT stop position. Sequences
        shorter than this cannot confirm stop loss — they may simply be
        truncated assemblies and are excluded.
      - The residue at the WT stop position must NOT be '*'.

    missense:
      - Handled by is_confirmed_isolate; only the generic length filter
        applies here.
    """
    # Generic: use full_seq for length check if available
    length_check_seq = full_seq if full_seq is not None else seq
    stripped = length_check_seq.replace("-", "").replace("!", "")
    if len(stripped) <= 10:
        return False

    # Frameshift
    if mutation_type == "frameshift":
        check     = full_seq if full_seq is not None else seq
        has_bang  = "!" in check
        clean_len = len(check.replace("-", "").replace("!", ""))
        if clean_len == 0:
            return False
        return has_bang or (clean_len % 3 != 0)

    # Inframe deletion / delins
    if mutation_type in ("inframe_deletion", "delins", "delins_del"):
        has_gap  = "-" in seq
        has_bang = "!" in seq
        all_gap  = all(c in ("-", "!") for c in seq) if seq else True
        if all_gap:
            return False
        if has_bang and not has_gap:
            return False
        if not has_gap and not has_bang:
            return False
        return True

    # Stop gained
    if mutation_type == "stop_gained":
        if seq.count("*") > 3:
            return False
        if wt_stop_pos is not None:
            first_stop = next((i + 1 for i, aa in enumerate(seq) if aa == "*"), None)
            if first_stop is None or first_stop >= wt_stop_pos:
                return False
        return True

    # Stop lost
    if mutation_type == "stop_lost":
        if wt_stop_pos is not None:
            ungapped_aa = [aa for aa in seq if aa not in ("-", "!")]
            if len(ungapped_aa) < wt_stop_pos:
                return False
            if ungapped_aa[wt_stop_pos - 1] == "*":
                return False
        return True

    # Missense — generic filter only
    return True


# ============================================================
# Extraction functions
# ============================================================

def extract_macse_aa(macse_dir, mut_list, confirmed, any_hits, strict_only, output_dir, context):
    aa_files   = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))
    macse_types = {"missense", "inframe_deletion", "delins"}

    # Extract if macse_strict=1 OR any_missense=1 — both are worth verifying
    to_extract = [(g, s, e, m, x) for g, s, e, m, x in mut_list
                  if m in macse_types and
                  ((g, s, e, m) in confirmed or (g, s, e, m) in any_hits)]

    print(f"MACSE AA: {len(to_extract)} mutation(s) to extract", file=sys.stderr)

    for gene, start, end, mtype, expected in to_extract:
        matching = [f for f in aa_files
                    if os.path.basename(f).startswith(gene + "_")]
        if not matching:
            matching = [f for f in aa_files
                        if gene.lower() in os.path.basename(f).lower()]
        if not matching:
            print(f"  ⚠  No AA alignment for {gene}", file=sys.stderr)
            continue

        label    = f"{gene}_{start}" + (f"_{end}" if end != start else "") + f"_{mtype}"
        out_path = os.path.join(output_dir, f"{label}_confirmed.fasta")
        written  = 0

        with open(out_path, "w") as fh:
            for aln in sorted(matching):
                records = list(SeqIO.parse(aln, "fasta"))
                if len(records) < 2:
                    continue
                anc            = str(records[0].seq)
                col_lo, col_hi = get_window_columns(anc, start, end, context)
                if written == 0:
                    fh.write(f">Ancestor | {gene} pos {start}-{end}\n")
                    fh.write(anc[col_lo:col_hi] + "\n")
                    written += 1
                for rec in records[1:]:
                    iso = str(rec.seq)
                    # is_any_only: the mutation has any_missense evidence but is
                    # NOT in the strict-only set (i.e. expected AA not matched).
                    # A mutation in both strict_hits and any_hits is NOT any_only
                    # — it has an exact residue match and must be checked as such.
                    key = (gene, start, end, mtype)
                    is_any_only = (key in any_hits) and (key not in strict_only)
                    # For any_missense-only mutations, accept any AA difference
                    # rather than requiring the specific expected AA
                    check_expected = None if is_any_only else expected
                    if is_confirmed_isolate(anc, iso, start, end, mtype, check_expected) \
                            and is_clean_sequence(iso[col_lo:col_hi], mtype, anc, full_seq=iso):
                        # Label distinguishes strict (exact AA) from any_missense
                        if key in strict_only or (key in any_hits and key in confirmed - any_hits):
                            ev_label = f"{mtype} confirmed"
                        else:
                            ev_label = f"{mtype} any_missense"
                        fh.write(f">{rec.id} | {ev_label}\n")
                        fh.write(iso[col_lo:col_hi] + "\n")
                        written += 1

        n = written - 1
        if n > 0:
            print(f"  ✅ {gene} {start}-{end} ({mtype}): {n} isolate(s) → {os.path.basename(out_path)}",
                  file=sys.stderr)
        else:
            os.remove(out_path)
            print(f"  ⚠  {gene} {start}-{end} ({mtype}): no sequences matched",
                  file=sys.stderr)


def extract_deletion_component(macse_dir, mut_list, exact_deletion_hits,
                               output_dir, context):
    """
    For delins mutations where exact_deletion=1 but macse_strict=0 —
    the deletion component is confirmed but the insertion is absent.
    Extracts isolates showing gaps at the deletion window, labelled
    as 'deletion component confirmed'.
    """
    aa_files = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))
    print(f"Deletion components: scanning hits...", file=sys.stderr)

    seen = set()
    for gene, start, end, mtype, expected in mut_list:
        if mtype != "delins":
            continue
        key = (gene, start, end, mtype)
        if key in seen or key not in exact_deletion_hits:
            continue
        seen.add(key)

        matching = [f for f in aa_files
                    if os.path.basename(f).startswith(gene + "_")]
        if not matching:
            matching = [f for f in aa_files
                        if gene.lower() in os.path.basename(f).lower()]
        if not matching:
            print(f"  ⚠  No AA alignment for {gene}", file=sys.stderr)
            continue

        label    = f"{gene}_{start}_{end}_delins_deletion_component"
        out_path = os.path.join(output_dir, f"{label}_confirmed.fasta")
        written  = 0

        with open(out_path, "w") as fh:
            for aln in sorted(matching):
                records = list(SeqIO.parse(aln, "fasta"))
                if len(records) < 2:
                    continue
                anc            = str(records[0].seq)
                col_lo, col_hi = get_window_columns(anc, start, end, context)

                if written == 0:
                    fh.write(f">Ancestor | {gene} pos {start}-{end} (delins — deletion component)\n")
                    fh.write(anc[col_lo:col_hi] + "\n")
                    written += 1

                # Check for deletion at the window (ignore insertion requirement)
                col_start = ungapped_to_alignment_index(anc, start)
                if col_start is None:
                    continue
                span, count, col_end = end - start + 1, 0, col_start
                while col_end < len(anc) and count < span:
                    if anc[col_end] not in ("-", "!"):
                        count += 1
                    col_end += 1

                for rec in records[1:]:
                    iso      = str(rec.seq)
                    anc_w    = anc[col_start:col_end]
                    iso_w    = iso[col_start:col_end]
                    anc_ug   = sum(1 for a in anc_w if a not in ("-", "!"))
                    iso_gaps = sum(1 for a in iso_w if a in ("-", "!"))
                    exp_del  = end - start + 1

                    # Require gaps to fall within the called window, not just
                    # regionally — this prevents larger upstream deletions from
                    # being counted as confirmation of the specific boundary.
                    # The gap must start at or before the first residue of the
                    # deletion window (col_start) and cover at least exp_del
                    # positions within the window itself.
                    first_gap_in_window = next(
                        (i for i, a in enumerate(iso_w) if a in ("-", "!")), None
                    )
                    gap_starts_at_boundary = (
                        first_gap_in_window is not None and first_gap_in_window == 0
                    )

                    if (anc_ug >= exp_del * 0.8
                            and iso_gaps >= exp_del
                            and gap_starts_at_boundary
                            and is_clean_sequence(iso[col_lo:col_hi], "delins_del",
                                                  anc, full_seq=iso)):
                        fh.write(f">{rec.id} | deletion component confirmed (no insertion)\n")
                        fh.write(iso[col_lo:col_hi] + "\n")
                        written += 1

        n = written - 1
        if n > 0:
            print(f"  ✅ {gene} {start}-{end} (delins deletion component): "
                  f"{n} isolate(s) → {os.path.basename(out_path)}", file=sys.stderr)
        else:
            os.remove(out_path)
            print(f"  ⚠  {gene} {start}-{end} (delins deletion component): "
                  f"no sequences matched", file=sys.stderr)


def extract_regional_deletion_component(macse_dir, mut_list, regional_hits,
                                         output_dir, context):
    """
    For delins mutations where regional_indel=1 but exact_deletion=0 —
    a deletion was confirmed in the region around the called window but the
    specific boundary was not matched. Sequences are extracted and labelled
    'regional deletion confirmed' to be honest that the exact event was not
    replicated in clinical isolates, only a deletion in the same region.
    """
    aa_files = glob.glob(os.path.join(macse_dir, "*_aligned_aa.fasta"))
    print(f"Regional deletion components: scanning hits...", file=sys.stderr)

    REGIONAL_WINDOW = 20

    seen = set()
    for gene, start, end, mtype, expected in mut_list:
        if mtype != "delins":
            continue
        key = (gene, start, end, mtype)
        if key in seen or key not in regional_hits:
            continue
        seen.add(key)

        matching = [f for f in aa_files
                    if os.path.basename(f).startswith(gene + "_")]
        if not matching:
            matching = [f for f in aa_files
                        if gene.lower() in os.path.basename(f).lower()]
        if not matching:
            print(f"  ⚠  No AA alignment for {gene}", file=sys.stderr)
            continue

        label    = f"{gene}_{start}_{end}_delins_regional_deletion"
        out_path = os.path.join(output_dir, f"{label}_confirmed.fasta")
        written  = 0

        region_lo = max(1, start - REGIONAL_WINDOW)
        region_hi = end + REGIONAL_WINDOW

        with open(out_path, "w") as fh:
            for aln in sorted(matching):
                records = list(SeqIO.parse(aln, "fasta"))
                if len(records) < 2:
                    continue
                anc            = str(records[0].seq)
                col_lo, col_hi = get_window_columns(anc, start, end, context)

                if written == 0:
                    fh.write(f">Ancestor | {gene} pos {start}-{end} "
                             f"(delins — regional deletion)\n")
                    fh.write(anc[col_lo:col_hi] + "\n")
                    written += 1

                # Get alignment columns for the regional window
                col_region_lo = ungapped_to_alignment_index(anc, region_lo) or 0
                col_region_hi_idx = ungapped_to_alignment_index(anc, region_hi)
                col_region_hi = (col_region_hi_idx + 1) if col_region_hi_idx else len(anc)

                for rec in records[1:]:
                    iso = str(rec.seq)

                    # Require the isolate to have a gap at a position where
                    # the ancestor has a RESIDUE (not an alignment gap).
                    # This prevents sequences that are identical to the ancestor
                    # — which also contain '-' characters from the alignment —
                    # from being incorrectly included.
                    anc_residue_cols = [i for i, a in enumerate(anc[col_region_lo:col_region_hi],
                                                                 start=col_region_lo)
                                        if a not in ("-", "!")]
                    gaps_at_residues = sum(
                        1 for col in anc_residue_cols
                        if col < len(iso) and iso[col] in ("-", "!")
                    )

                    if gaps_at_residues > 0 and is_clean_sequence(
                            iso[col_lo:col_hi], "delins_del", anc, full_seq=iso):
                        fh.write(f">{rec.id} | regional deletion confirmed\n")
                        fh.write(iso[col_lo:col_hi] + "\n")
                        written += 1

        n = written - 1
        if n > 0:
            print(f"  ✅ {gene} {start}-{end} (delins regional deletion): "
                  f"{n} isolate(s) → {os.path.basename(out_path)}", file=sys.stderr)
        else:
            os.remove(out_path)
            print(f"  ⚠  {gene} {start}-{end} (delins regional deletion): "
                  f"no sequences matched", file=sys.stderr)


def extract_frameshifts(macse_dir, mut_list, fs_hits, output_dir, context):
    nt_files = glob.glob(os.path.join(macse_dir, "*_aligned_nt.fasta"))
    print(f"Frameshifts: scanning hits...", file=sys.stderr)

    # Build gene-level lookup: gene → set of confirmed isolate IDs
    # This covers exact, regional, and gene-tier hits
    gene_to_isolates = {}
    for (gene, pos), hit_list in fs_hits.items():
        for tier, iso in hit_list:
            gene_to_isolates.setdefault(gene, set()).add(iso)

    # Deduplicate targets — one output file per unique (gene, start)
    seen = set()
    for gene, start, end, mtype, _ in mut_list:
        if mtype != "frameshift":
            continue
        if (gene, start) in seen:
            continue
        seen.add((gene, start))

        if gene not in gene_to_isolates:
            print(f"  ⚠  No frameshift hits for {gene}", file=sys.stderr)
            continue

        confirmed_ids = gene_to_isolates[gene]
        matching = [f for f in nt_files
                    if os.path.basename(f).startswith(gene + "_")]
        if not matching:
            matching = [f for f in nt_files
                        if gene.lower() in os.path.basename(f).lower()]
        if not matching:
            print(f"  ⚠  No NT alignment for {gene}", file=sys.stderr)
            continue

        out_path = os.path.join(output_dir, f"{gene}_{start}_frameshift_confirmed.fasta")
        written  = 0
        rejected = 0

        with open(out_path, "w") as fh:
            for aln in sorted(matching):
                records = list(SeqIO.parse(aln, "fasta"))
                if len(records) < 2:
                    continue
                anc = str(records[0].seq)

                # Find NT window: (start-context)*3 to (start+context)*3
                nt_lo = max(1, (start - context) * 3)
                nt_hi = (start + context) * 3
                col_lo, col_hi = 0, len(anc)
                count = 0
                for i, b in enumerate(anc):
                    if b not in ("-", "!"):
                        count += 1
                    if count == nt_lo:
                        col_lo = i
                    if count == nt_hi:
                        col_hi = i + 1
                        break

                if written == 0:
                    fh.write(f">Ancestor | {gene} codon {start} (frameshift)\n")
                    fh.write(anc[col_lo:col_hi] + "\n")
                    written += 1

                for rec in records[1:]:
                    if rec.id not in confirmed_ids:
                        continue
                    full_seq = str(rec.seq)
                    window   = full_seq[col_lo:col_hi]

                    # Check the window itself is not all-gap before the
                    # full-sequence frameshift check — a sequence can have
                    # residues outside the window but contribute nothing
                    # to the windowed excerpt if that region is all-gap.
                    window_ungapped = len(window.replace("-","").replace("!",""))
                    if window_ungapped == 0:
                        rejected += 1
                        print(f"    ✗ {rec.id}: all-gap in excerpt window",
                              file=sys.stderr)
                        continue

                    if not is_clean_sequence(window, "frameshift", full_seq=full_seq):
                        rejected += 1
                        print(f"    ✗ {rec.id}: failed frameshift validation "
                              f"(no '!' and length divisible by 3)",
                              file=sys.stderr)
                        continue
                    tier = next((t for (g, p), hits in fs_hits.items()
                                 if g == gene
                                 for t, i in hits
                                 if i == rec.id), "gene")
                    fh.write(f">{rec.id} | frameshift {tier}\n")
                    fh.write(window + "\n")
                    written += 1

        n = written - 1
        if n > 0:
            msg = f"{n} isolate(s)"
            if rejected:
                msg += f", {rejected} rejected (no frameshift evidence in sequence)"
            print(f"  ✅ {gene} {start} (frameshift): {msg} → {os.path.basename(out_path)}",
                  file=sys.stderr)
        else:
            os.remove(out_path)
            print(f"  ⚠  {gene} {start} (frameshift): no sequences matched",
                  file=sys.stderr)


def extract_stops(per_gene_fastas, mut_list, stop_hits, output_dir, context):
    """
    Extract confirmed stop_gained / stop_lost sequences.
    Only includes isolates where the stop event is at the exact called codon
    or within ±REGIONAL_WINDOW — consistent with frameshift and indel extraction.
    Gene-level stop events (stop anywhere in gene, unrelated to called position)
    are not written to the confirmed FASTA.
    """
    REGIONAL_WINDOW = 20

    fasta_files = glob.glob(os.path.join(per_gene_fastas, "*.fasta"))
    print(f"Stop mutations: scanning hits...", file=sys.stderr)

    def translate(s):
        s = s.upper().replace("-", "").replace("!", "")
        r = len(s) % 3
        if r:
            s = s[:-r]
        return str(Seq(s).translate(table=11))

    seen = set()
    for gene, start, end, mtype, _ in mut_list:
        norm = {"stop_gain":"stop_gained","stop_gained":"stop_gained",
                "stop_loss":"stop_lost","stop_lost":"stop_lost"}.get(mtype, mtype)
        if norm not in ("stop_gained", "stop_lost"):
            continue
        if (gene, norm) in seen:
            continue
        seen.add((gene, norm))

        gene_hits = stop_hits.get(gene, [])
        if not gene_hits:
            print(f"  ⚠  No stop hits for {gene}", file=sys.stderr)
            continue

        # Build confirmed_ids with tier info so we can label sequences correctly.
        # exact/regional → position-specific evidence
        # gene (pos>0)   → stop anywhere in gene, unrelated to called position
        # gene (pos=0)   → old format, no position info
        confirmed_ids        = {}   # isolate_id → label
        for pos, tier, isolate in gene_hits:
            if tier == "exact":
                dist = abs(pos - start)
                if dist == 0:
                    confirmed_ids[isolate] = f"{mtype} exact confirmed"
                elif dist <= REGIONAL_WINDOW:
                    confirmed_ids[isolate] = f"{mtype} regional confirmed"
            elif tier == "regional":
                dist = abs(pos - start)
                if dist <= REGIONAL_WINDOW:
                    confirmed_ids[isolate] = f"{mtype} regional confirmed"
            elif tier == "gene":
                if pos == 0:
                    # Old format — no position, include as gene-level
                    confirmed_ids[isolate] = f"{mtype} gene confirmed"
                else:
                    # New format gene-tier: stop at unrelated position — include
                    # but label clearly so it's distinguishable
                    confirmed_ids[isolate] = f"{mtype} gene confirmed (pos {pos})"

        if not confirmed_ids:
            print(f"  ⚠  No stop hits for {gene} {norm} at pos {start}",
                  file=sys.stderr)
            continue

        matching = [f for f in fasta_files
                    if os.path.basename(f).startswith(gene + "_")
                    or os.path.basename(f) == f"{gene}.fasta"]
        if not matching:
            print(f"  ⚠  No per-gene FASTA for {gene}", file=sys.stderr)
            continue

        out_path = os.path.join(output_dir, f"{gene}_{start}_{mtype}_confirmed.fasta")
        written  = 0
        rejected = 0

        with open(out_path, "w") as fh:
            for fasta_path in sorted(matching):
                records = list(SeqIO.parse(fasta_path, "fasta"))
                if len(records) < 2:
                    continue

                wt_aa       = translate(str(records[0].seq))
                wt_stop_pos = next((i+1 for i,aa in enumerate(wt_aa) if aa=="*"), None)

                if written == 0:
                    fh.write(f">Ancestor | {gene} (WT protein)\n")
                    fh.write(wt_aa + "\n")
                    written += 1

                for rec in records[1:]:
                    if rec.id not in confirmed_ids:
                        continue
                    aa = translate(str(rec.seq))
                    if not is_clean_sequence(aa, norm, wt_stop_pos=wt_stop_pos):
                        rejected += 1
                        print(f"    ✗ {rec.id}: failed sequence-level validation",
                              file=sys.stderr)
                        continue
                    ev_label = confirmed_ids[rec.id]
                    fh.write(f">{rec.id} | {ev_label}\n")
                    fh.write(aa + "\n")
                    written += 1

        n = written - 1
        if n > 0:
            msg = f"{n} isolate(s)"
            if rejected:
                msg += f", {rejected} rejected"
            print(f"  ✅ {gene} {start} ({mtype}): {msg} → {os.path.basename(out_path)}",
                  file=sys.stderr)
        else:
            os.remove(out_path)
            print(f"  ⚠  {gene} {start} ({mtype}): no sequences passed validation",
                  file=sys.stderr)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extract confirmed mutation sequences for visual verification"
    )
    parser.add_argument("macse_dir",         help="MACSE output directory")
    parser.add_argument("summary",           help="mutation_cds_macse_strict_summary.tsv")
    parser.add_argument("mutations",         help="mutation_list.tsv")
    parser.add_argument("output_dir",        help="Output directory")
    parser.add_argument("--frameshifts",     default=None,
                        help="frameshifts.tsv")
    parser.add_argument("--stops",           default=None,
                        help="stop_genes.tsv")
    parser.add_argument("--per-gene-fastas", default=None,
                        help="per_gene_fastas/ directory (required for stop extraction)")
    parser.add_argument("--screen-results",  default=None,
                        help="mutation_screen_results.tsv (for deletion component extraction)")
    parser.add_argument("--context",         type=int, default=10,
                        help="Residues either side of mutation window (default: 10)")
    parser.add_argument("--any",             action="store_true",
                        help="Also extract any_missense hits")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    confirmed, any_hits, strict_only = load_confirmed(args.summary, include_any=args.any)
    mut_list  = load_mutation_list(args.mutations)

    # MACSE AA — missense, inframe_deletion, delins
    extract_macse_aa(args.macse_dir, mut_list, confirmed, any_hits, strict_only,
                     args.output_dir, args.context)

    # Deletion component — delins where deletion confirmed but insertion absent
    if args.screen_results:
        exact_deletion_hits = load_exact_deletion_hits(args.screen_results)
        regional_deletion_hits = load_regional_deletion_hits(args.screen_results)

        # Exact deletion component — gap matches called window precisely
        delins_not_strict = exact_deletion_hits - confirmed
        if delins_not_strict:
            extract_deletion_component(args.macse_dir, mut_list, delins_not_strict,
                                       args.output_dir, args.context)

        # Regional deletion — gap confirmed in region but boundary not exact
        # Only for delins not already covered by exact or strict confirmation
        delins_regional_only = regional_deletion_hits - confirmed - exact_deletion_hits
        if delins_regional_only:
            extract_regional_deletion_component(args.macse_dir, mut_list,
                                                delins_regional_only,
                                                args.output_dir, args.context)
    else:
        print("⚠  --screen-results not supplied — skipping deletion component extraction",
              file=sys.stderr)

    # Frameshifts
    if args.frameshifts:
        fs_hits = load_frameshift_hits(args.frameshifts)
        extract_frameshifts(args.macse_dir, mut_list, fs_hits,
                            args.output_dir, args.context)
    else:
        print("\n⚠  --frameshifts not supplied — skipping frameshift extraction",
              file=sys.stderr)

    # Stop mutations
    if args.stops and args.per_gene_fastas:
        stop_hits = load_stop_hits(args.stops)
        extract_stops(args.per_gene_fastas, mut_list, stop_hits,
                      args.output_dir, args.context)
    else:
        print("⚠  --stops/--per-gene-fastas not supplied — skipping stop extraction",
              file=sys.stderr)

    print(f"\n✅ Done — output in {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
