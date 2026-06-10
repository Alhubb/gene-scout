#!/usr/bin/env python3
"""
detect_stop_mutations.py  <nucleotide_cds_fasta>  [mutation_list.tsv]

Translates CDS sequences and detects stop-gained / stop-lost mutations
relative to the first (wild-type) record.

Three-tier position-aware detection, consistent with detect_frameshifts_macse.py
and detect_indels_macse.py:

  Tier 1 — exact:
    stop_gained: premature stop at the exact called codon position.
    stop_lost:   stop codon absent at the exact WT stop position.

  Tier 2 — regional:
    stop_gained: premature stop within ±REGIONAL_WINDOW codons of called position.
    stop_lost:   stop absent and sequence extends past the WT stop position.

  Tier 3 — gene:
    Any premature stop anywhere in the gene (stop_gained), or any stop loss
    anywhere (stop_lost), regardless of position.

Output: four columns — gene_name  aa_position  tier  isolate_id

If no mutation_list is supplied, all hits are assigned to the gene tier.
"""

import sys
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq

REGIONAL_WINDOW = 20


def translate_nt(nt_str):
    seq = nt_str.upper().replace("-", "").replace("!", "")
    remainder = len(seq) % 3
    if remainder:
        seq = seq[:-remainder]
    return str(Seq(seq).translate(table=11))


def main(nt_fasta, called_positions=None):
    """
    called_positions: dict of {gene: {event_type: aa_position}} — the specific
    codon positions from the mutation list, used for exact/regional tiering.
    If not supplied, all hits are written as 'gene' tier.
    """
    records = list(SeqIO.parse(nt_fasta, "fasta"))
    if len(records) < 2:
        return

    gene = os.path.basename(nt_fasta)
    gene = gene.replace(".fasta", "").replace(".fna", "")
    gene = gene.rsplit("_", 1)[0]

    wt_aa = translate_nt(str(records[0].seq))
    # wt_stop_pos is 1-based
    wt_stop_pos = wt_aa.find("*") + 1 if "*" in wt_aa else len(wt_aa) + 1

    gene_calls    = called_positions.get(gene, {}) if called_positions else {}
    called_gain   = gene_calls.get("stop_gained")
    called_lost   = gene_calls.get("stop_lost")

    for rec in records[1:]:
        aa = translate_nt(str(rec.seq))
        iso_stop = aa.find("*") + 1 if "*" in aa else None

        # ── stop_gained ───────────────────────────────────────────────────
        if iso_stop is not None and iso_stop < wt_stop_pos:
            pos = iso_stop
            if called_gain is not None:
                dist = abs(pos - called_gain)
                tier = "exact" if dist == 0 else ("regional" if dist <= REGIONAL_WINDOW else "gene")
            else:
                tier = "gene"
            print(f"{gene}\t{pos}\t{tier}\t{rec.id}")
            continue

        # ── stop_lost ─────────────────────────────────────────────────────
        if iso_stop is None or iso_stop > wt_stop_pos:
            pos = wt_stop_pos
            if called_lost is not None:
                dist = abs(pos - called_lost)
                tier = "exact" if dist == 0 else ("regional" if dist <= REGIONAL_WINDOW else "gene")
            else:
                tier = "gene"
            print(f"{gene}\t{pos}\t{tier}\t{rec.id}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: detect_stop_mutations.py <nucleotide_cds_fasta> [mutation_list.tsv]")

    nt_fasta = sys.argv[1]
    called_positions = {}

    if len(sys.argv) >= 3:
        with open(sys.argv[2]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                mtype = row.get("mutation_type", "").strip()
                if mtype in ("stop_gained", "stop_lost"):
                    g = row.get("gene", "").strip()
                    try:
                        pos = int(row["aa_start"])
                    except (KeyError, ValueError):
                        continue
                    called_positions.setdefault(g, {})[mtype] = pos

    main(nt_fasta, called_positions)
