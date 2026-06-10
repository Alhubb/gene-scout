# GENE-SCOUT

GENE-SCOUT is a gene-centric bioinformatics pipeline for detecting coding-sequence (CDS) mutations across large collections of genome assemblies using nucleotide-resolved, tiered evidence.

⚠️ **Software status**
This repository is released as **v0.1** and has **not yet been formally validated**.
While the pipeline is functional and has been applied to real datasets, results should be interpreted with appropriate caution. Validation and benchmarking will be addressed in future releases.

---

## Overview

GENE-SCOUT is designed for population-scale analysis of bacterial genomes. Given a set of genome assemblies and a set of ancestral gene CDS sequences, it extracts and aligns the corresponding gene from every assembly, then screens for mutations using multiple, biologically appropriate evidence sources.

Each mutation is assigned evidence across three tiers:

- **Tier 1 — Exact:** the specific mutation called in the query (exact amino acid, exact frameshift position, exact stop codon position)
- **Tier 2 — Regional:** the same mutation type within ±20 amino acids of the called position
- **Tier 3 — Gene:** any mutation of the same type anywhere in the gene

For missense mutations only two tiers are used (exact residue and any amino acid change at the same position), as gene-level missense evidence is too non-specific.

---

## Mutation classes detected

- **Missense substitutions**
- **In-frame insertions and deletions**
- **Complex substitutions (delins)**
- **Frameshift mutations**
- **Stop-gain mutations** (premature termination)
- **Stop-loss mutations** (loss of the ancestral stop codon)

---

## Input requirements

### Genome assemblies (`--genomes`)

Directory containing genome assemblies as nucleotide FASTA files.

Supported formats: `.fasta`, `.fna`, `.fna.gz`

Assemblies may contain one or multiple contigs. GENE-SCOUT does not assume any specific genome source.

### Ancestor genes (`--ancestors`)

Directory containing ancestral gene coding sequences (CDS only) as nucleotide FASTA files.

Supported formats: `.fasta`, `.fna`

> ❗ Protein sequences are **not supported**. Each file must contain a single complete CDS (no UTRs, no introns).

### Mutation list (`--mutations`, optional)

A tab-separated file defining the mutations to screen for. Required columns:

```
gene    aa_start    aa_end    mutation_type    expected
```

Where `mutation_type` is one of: `missense`, `frameshift`, `inframe_deletion`, `delins`, `stop_gained`, `stop_lost`.

The `expected` column contains the mutant amino acid for missense, the inserted residue for delins, and is left blank for all other types.

If you have Snippy output, use the provided helper script to convert it:

```bash
python scripts/convert_snippy_to_mutation_list.py snippy_summary.tsv mutation_list.tsv
```

---

## Installation

```bash
git clone https://github.com/alhubb/gene-scout
cd gene-scout
conda env create -f environment.yml
conda activate gene-scout
chmod +x gene-scout
```

---

## Usage

### Basic — detect all mutation types across a genome collection

```bash
gene-scout \
  --genomes path/to/genomes/ \
  --ancestors path/to/WT_gene_fasta/ \
  --out results/
```

### With mutation screening — screen for specific predefined mutations

```bash
gene-scout \
  --genomes genomes/ \
  --ancestors WT_gene_fasta/ \
  --mutations mutation_list.tsv \
  --out results/
```

### All arguments

| Argument | Required | Description |
|---|---|---|
| `--genomes DIR` | ✅ | Directory of genome assemblies (`.fasta`, `.fna`, `.fna.gz`) |
| `--ancestors DIR` | ✅ | Directory of ancestral gene CDS FASTA files |
| `--out DIR` | ✅ | Output directory for all results and intermediate files |
| `--mutations TSV` | optional | Mutation list TSV to screen against |
| `--threads N` | optional | Parallel threads for minimap2 and MACSE (default: 4) |

---

## Pipeline steps

1. **Align** — each ancestral gene CDS is aligned against all genome assemblies using minimap2
2. **Extract** — gene sequences are extracted from assemblies using alignment coordinates
3. **Deduplicate** — exact duplicate CDS sequences are removed with seqkit
4. **Clean sequences** — short, near-empty, or non-IUPAC sequences are removed
5. **Align (MACSE)** — codon-aware multiple alignment is performed with MACSE
6. **Clean alignments** — sequences that collapse to near-entirely gaps after MACSE alignment are removed
7. **Scan for mutations** — MACSE amino-acid alignments are scanned for missense and in-frame indels
8. **Detect frameshifts** — MACSE `!` markers in NT alignments are mapped to codon positions with tiered evidence
9. **Detect indels** — in-frame gaps in AA alignments are detected with tiered evidence
10. **Detect stop mutations** — CDS sequences are translated and compared to the ancestor to detect stop-gain and stop-loss events with tiered evidence
11. **Screen mutations** — all evidence is merged and reported per mutation with tier assignments
12. **Extract confirmed sequences** — aligned sequences for each confirmed mutation are written to per-mutation FASTA files
13. **Count evidence** — unique isolates with genuine sequence-level evidence are counted per mutation

---

## Output files

All outputs are written to the directory specified by `--out`.

| File | Description |
|---|---|
| `mutation_screen_results.tsv` | Primary output — per-mutation evidence summary with tier flags and `final_confirmed` column |
| `confirmed_sequences/` | Per-mutation FASTA files containing the aligned sequences that support each call |
| `confirmed_sequences/confirmed_evidence_counts.tsv` | Counts of unique isolates with sequence-level evidence per mutation |
| `mutation_cds_macse_strict_summary.tsv` | MACSE alignment scan results (missense and indel calls) |
| `frameshifts.tsv` | Frameshift hits with tier and isolate ID |
| `indels.tsv` | In-frame indel hits with tier and isolate ID |
| `stop_genes.tsv` | Stop-gain and stop-loss hits with tier and isolate ID |
| `per_gene_fastas_dedup/` | Deduplicated CDS FASTAs |
| `per_gene_fastas_clean/` | Cleaned CDS FASTAs (input to MACSE) |
| `macse_output/` | Raw MACSE alignments |
| `macse_output_clean/` | Cleaned MACSE alignments (used by all detection steps) |

---

## Debugging missense and delins calls

If expected mutations are not being detected, use the debug helper to inspect what amino acids are present at each position in the MACSE alignment:

```bash
python scripts/debug_macse_windows.py \
  results/macse_output_clean/ \
  mutation_list.tsv
```

This prints the ancestor residues and the distribution of isolate residues at each missense and delins window, allowing you to verify the `expected` column in your mutation list.

---

## Dependencies

| Tool | Version | Citation |
|---|---|---|
| minimap2 | 2.30 | Li (2018) Bioinformatics |
| MACSE | 2.07 | Ranwez et al. (2011) PLoS ONE; Ranwez et al. (2018) MBE |
| seqkit | 2.13.0 | Shen et al. (2016) PLoS ONE |
| Biopython | 1.87 | Cock et al. (2009) Bioinformatics |
| Python | ≥3.9 | |

---

## License

GENE-SCOUT is released under the MIT License. See the LICENSE file for details.

---

## Development status

GENE-SCOUT is under active development. Future work will include formal validation, automated testing, additional input safeguards, and performance optimisation.

## Development status
GENE-SCOUT is under active development.
Future work will include formal validation, automated testing,
additional input safeguards, and performance optimisation.
