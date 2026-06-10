#!/usr/bin/env python3
"""
count_confirmed_evidence.py  <confirmed_sequences_dir>

Scans all confirmed mutation FASTA files in the given directory and counts
how many sequences carry genuine sequence-level evidence for their called
mutation type.

The ancestor (first record in each FASTA) is excluded from all counts.

Evidence checked per mutation type:

  frameshift:
    Confirmed FASTAs for frameshifts are written from MACSE NT alignments
    by extract_confirmed_sequences.py, so sequences are nucleotide.
    Evidence requires either:
      - '!' marker anywhere in the aligned NT sequence (MACSE frameshift
        marker), OR
      - Ungapped NT length not divisible by 3.
    All-gap sequences (no ungapped characters) are excluded.
    NOTE: if a frameshift FASTA unexpectedly contains AA sequences (e.g.
    from a manual re-run), the length%3 check will be unreliable.

  inframe_deletion / delins / regional_deletion:
    - At least one gap ('-' or '!') at a position where the ancestor
      has a residue (i.e. a real deletion, not an alignment gap).
    Sequences identical to the ancestor (no real gap) are excluded.

  stop_gained:
    - A stop codon ('*') present upstream of the ancestor's stop position.

  stop_lost:
    - No stop codon at the ancestor's stop position (residue at that
      position is not '*').
    Sequences too short to reach the ancestor stop position are excluded.

  missense:
    - At least one amino acid difference from the ancestor in the
      ungapped sequence. All-gap sequences are excluded.

Output:
    Tab-separated table printed to stdout with columns:
      file  mutation  total_samples  unique_samples  with_evidence
      unique_with_evidence  evidence_rate

Usage:
    python count_confirmed_evidence.py <confirmed_sequences_dir>

    # Or for a specific directory:
    python count_confirmed_evidence.py /path/to/confirmed_sequences/
"""

import os
import sys
import re
import glob
from collections import Counter


# ============================================================
# FASTA parsing
# ============================================================

def parse_fasta(path):
    """Return list of (header, seq) tuples."""
    records = []
    with open(path) as fh:
        header, seq = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq)))
    return records


def ungapped_len(seq):
    return sum(1 for c in seq if c not in ("-", "!"))


# ============================================================
# Mutation type inference from filename
# ============================================================

# Filename pattern: gene_start[_end]_muttype[_qualifier]_confirmed.fasta
# e.g. glpT_294_frameshift_confirmed.fasta
#      uhpA_41_43_inframe_deletion_confirmed.fasta
#      glpT_182_187_delins_regional_deletion_confirmed.fasta
#      uhpC_440_stop_lost_confirmed.fasta
#      uhpT_7_stop_gained_confirmed.fasta

MTYPE_KEYWORDS = [
    "frameshift",
    "inframe_deletion",
    "delins_regional_deletion",
    "delins",
    "stop_gained",
    "stop_lost",
    "missense",
]

def infer_mutation_type(filename):
    """
    Infer mutation type from the filename.
    Returns one of: frameshift, inframe_deletion, delins,
                    regional_deletion, stop_gained, stop_lost,
                    missense, or unknown.
    """
    base = filename.lower().replace(".fasta", "")
    # Order matters — check more specific patterns first
    for kw in MTYPE_KEYWORDS:
        if kw in base:
            return kw
    return "unknown"


# ============================================================
# Evidence checks
# ============================================================

def real_gaps(anc_seq, iso_seq):
    """
    Count gaps in iso_seq at positions where anc_seq has a residue.
    These are genuine deletions relative to the ancestor.
    """
    return sum(
        1 for i, c in enumerate(anc_seq)
        if c not in ("-", "!") and i < len(iso_seq) and iso_seq[i] in ("-", "!")
    )


def is_nucleotide_seq(seq):
    """
    Return True if the sequence looks like nucleotide (contains only ACGTN-!
    characters after stripping gaps). AA sequences contain letters outside
    this set (e.g. D, E, F, H, I, K, L, M, P, Q, R, S, V, W, Y, *).
    """
    nt_chars = set("ACGTNacgtn-!")
    return all(c in nt_chars for c in seq if c not in ("-", "!"))


def has_evidence(anc_seq, iso_seq, mtype, wt_stop_pos):
    """
    Return (has_evidence: bool, reason: str) for a single isolate sequence.

    anc_seq      — ancestor sequence (AA or NT depending on mtype)
    iso_seq      — isolate sequence
    mtype        — inferred mutation type string
    wt_stop_pos  — 1-based position of WT stop codon in AA sequence,
                   or None if not applicable
    """
    ul = ungapped_len(iso_seq)

    # All-gap — no evidence regardless of type
    if ul == 0:
        return False, "all_gap"

    if mtype == "frameshift":
        has_bang = "!" in iso_seq
        if has_bang:
            return True, "bang_marker"
        # length%3 check is only meaningful for NT sequences.
        # Frameshift confirmed FASTAs are written from MACSE NT alignments,
        # so this should always be NT — but guard against AA sequences
        # (which would give spurious length%3 hits) by checking the alphabet.
        if not is_nucleotide_seq(iso_seq):
            return False, "aa_sequence_in_frameshift_fasta(cannot_check_length)"
        len_fs = ul % 3 != 0
        if len_fs:
            return True, "length_frameshift"
        return False, "no_frameshift_signal"

    if mtype in ("inframe_deletion", "delins",
                 "delins_regional_deletion", "regional_deletion"):
        n_gaps = real_gaps(anc_seq, iso_seq)
        if n_gaps > 0:
            return True, f"real_gaps({n_gaps})"
        return False, "no_real_gap"

    if mtype == "stop_gained":
        if wt_stop_pos is None:
            # No WT stop — check for any stop codon
            if "*" in iso_seq:
                return True, "stop_present"
            return False, "no_stop"
        first_stop = next((i + 1 for i, aa in enumerate(iso_seq) if aa == "*"), None)
        if first_stop is not None and first_stop < wt_stop_pos:
            return True, f"premature_stop_at_{first_stop}"
        return False, f"no_premature_stop(wt_pos={wt_stop_pos})"

    if mtype == "stop_lost":
        if wt_stop_pos is None:
            return False, "no_wt_stop_defined"
        if ul < wt_stop_pos:
            return False, f"too_short({ul}<{wt_stop_pos})"
        # Check the residue at the WT stop position
        ungapped_aa = [aa for aa in iso_seq if aa not in ("-", "!")]
        aa_at_wt = ungapped_aa[wt_stop_pos - 1] if len(ungapped_aa) >= wt_stop_pos else "?"
        if aa_at_wt != "*":
            return True, f"stop_gone(aa={aa_at_wt})"
        return False, "stop_retained"

    if mtype == "missense":
        # Compare ungapped AA sequences — any difference is evidence
        # (covers both macse_strict and any_missense labelled sequences)
        anc_ungapped = "".join(aa for aa in anc_seq if aa not in ("-", "!"))
        iso_ungapped = "".join(aa for aa in iso_seq if aa not in ("-", "!"))
        if not iso_ungapped:
            return False, "all_gap"
        min_len = min(len(anc_ungapped), len(iso_ungapped))
        if any(anc_ungapped[i] != iso_ungapped[i] for i in range(min_len)):
            return True, "aa_difference"
        return False, "identical_to_ancestor"

    return False, f"unknown_type({mtype})"


# ============================================================
# Process one FASTA file
# ============================================================

def process_fasta(path):
    """
    Returns dict with:
      filename, mutation_type, total_samples, unique_samples,
      with_evidence, unique_with_evidence,
      evidence_rate, issues (list of problem sequences)
    """
    filename = os.path.basename(path)
    mtype    = infer_mutation_type(filename)

    records = parse_fasta(path)
    if len(records) < 2:
        return {
            "filename": filename, "mutation_type": mtype,
            "total_samples": 0, "unique_samples": 0,
            "with_evidence": 0, "unique_with_evidence": 0,
            "evidence_rate": 0.0, "issues": ["no_samples_in_fasta"],
        }

    anc_header, anc_seq = records[0]
    samples = records[1:]

    # Infer WT stop position from ancestor (AA sequences only)
    wt_stop_pos = None
    if mtype in ("stop_gained", "stop_lost", "missense"):
        wt_stop_pos = next(
            (i + 1 for i, aa in enumerate(anc_seq) if aa == "*"), None
        )

    # Track unique sample IDs for dedup count
    seen_ids      = set()
    seen_ids_evid = set()
    total_evid    = 0
    issues        = []

    for header, seq in samples:
        sample_id = header.split("|")[0].strip()
        is_new    = sample_id not in seen_ids
        seen_ids.add(sample_id)

        ok, reason = has_evidence(anc_seq, seq, mtype, wt_stop_pos)
        if ok:
            total_evid += 1
            seen_ids_evid.add(sample_id)
        else:
            issues.append((sample_id, reason))

    total   = len(samples)
    unique  = len(seen_ids)
    u_evid  = len(seen_ids_evid)
    rate    = total_evid / total if total > 0 else 0.0

    return {
        "filename":              filename,
        "mutation_type":         mtype,
        "total_samples":         total,
        "unique_samples":        unique,
        "with_evidence":         total_evid,
        "unique_with_evidence":  u_evid,
        "evidence_rate":         rate,
        "issues":                issues,
    }


# ============================================================
# Main
# ============================================================

def main():
    if len(sys.argv) not in (2, 3):
        sys.exit("Usage: count_confirmed_evidence.py <confirmed_sequences_dir> [output.tsv]")

    confirmed_dir = sys.argv[1]
    output_path   = sys.argv[2] if len(sys.argv) == 3 \
                    else os.path.join(confirmed_dir, "confirmed_evidence_counts.tsv")

    # If the supplied output path is an existing directory, write the
    # default filename inside it rather than erroring.
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "confirmed_evidence_counts.tsv")

    if not os.path.isdir(confirmed_dir):
        sys.exit(f"ERROR: not a directory: {confirmed_dir}")

    fasta_files = sorted(glob.glob(os.path.join(confirmed_dir, "*_confirmed.fasta")))
    if not fasta_files:
        sys.exit(f"ERROR: no *_confirmed.fasta files found in {confirmed_dir}")

    print(f"Found {len(fasta_files)} confirmed FASTA file(s) in {confirmed_dir}")

    cols = [
        "file", "mutation_type",
        "total_samples", "unique_samples",
        "with_evidence", "unique_with_evidence",
        "evidence_rate",
    ]

    results = []
    for path in fasta_files:
        results.append(process_fasta(path))

    with open(output_path, "w") as out:
        out.write("\t".join(cols) + "\n")

        for r in results:
            out.write("\t".join([
                r["filename"],
                r["mutation_type"],
                str(r["total_samples"]),
                str(r["unique_samples"]),
                str(r["with_evidence"]),
                str(r["unique_with_evidence"]),
                f"{r['evidence_rate']:.1%}",
            ]) + "\n")

        total_files = len(results)
        files_any   = sum(1 for r in results if r["unique_with_evidence"] > 0)
        files_none  = total_files - files_any

        out.write(f"\n# Files checked:              {total_files}\n")
        out.write(f"# Files with ≥1 confirmed:    {files_any}\n")
        out.write(f"# Files with no evidence:     {files_none}\n")

        any_issues = [(r["filename"], r["issues"]) for r in results if r["issues"]]
        if any_issues:
            out.write("\n# Sequences lacking evidence:\n")
            for fname, issues in any_issues:
                out.write(f"#   {fname}:\n")
                for sid, reason in issues:
                    out.write(f"#     x {sid[:80]}: {reason}\n")
        else:
            out.write("\n# All sequences in all files carry evidence\n")

    print(f"Results written to {output_path}")
    print(f"\nFiles checked:           {total_files}")
    print(f"Files with ≥1 confirmed: {files_any}")
    print(f"Files with no evidence:  {files_none}")
    if any_issues:
        print("\nSequences lacking evidence:")
        for fname, issues in any_issues:
            print(f"  {fname}: {len(issues)} sequence(s) without evidence")


if __name__ == "__main__":
    main()

