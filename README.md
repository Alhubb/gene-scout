# GENE-SCOUT

GENE-SCOUT is a gene-centric bioinformatics workflow for detecting
coding-sequence (CDS) mutations across genome assemblies using
nucleotide-resolved evidence.

⚠️ **Software status:**  
This repository is released as **v0.1** and has **not yet been formally validated**.
While the pipeline is functional, results should be interpreted with caution.
Validation and benchmarking will be addressed in future releases.

---

## Overview

GENE-SCOUT is designed for population-scale analysis of bacterial genomes.
It identifies mutations in specific genes by extracting coding sequences
directly from genome assemblies and analysing them using multiple,
biologically appropriate evidence sources.

The workflow explicitly distinguishes between:
- frame‑preserving mutations (e.g. missense, in‑frame indels), and
- frame‑disrupting mutations (frameshifts and stop codon changes).

---

## Mutation classes detected

GENE-SCOUT detects the following mutation types:

- **Missense substitutions**
- **In-frame insertions and deletions**
- **Frameshift mutations**
- **Stop-gain mutations** (premature termination)
- **Stop-loss mutations** (loss of the ancestral stop codon)

Each class is detected using a method appropriate to its underlying
molecular mechanism.

---

## Input requirements

### Genome assemblies (`--genomes`)

A directory containing genome assemblies as **nucleotide FASTA files**.
The following file types are supported:

- `.fasta`
- `.fna`
- `.fna.gz`

Assemblies may contain one or multiple contigs.
GENE-SCOUT does not assume any specific genome source
(e.g. BakRep, NCBI, local assemblies).

---

### Ancestor genes (`--ancestors`)

A directory containing **ancestral gene coding sequences (CDS only)**,
provided as **nucleotide FASTA files**:

- `.fasta`
- `.fna`

❗ **Important:**  
Protein sequences are **not supported**.

Each ancestor FASTA file should contain a single CDS corresponding to one
gene of interest and should represent the coding region only
(no UTRs or introns).

---

## Installation

GENE-SCOUT is implemented using Bash and Python and uses Conda to manage
dependencies.

Clone the repository:

```bash
git clone https://github.com/alhubb/gene-scout
cd gene-scout
``

## Creat and activate the Conda environment

conda env create -f environment.yml
conda activate gene-scout

## Make the main command executable 

chmod +x gene-s

## Useage

gene-scout \
  --genomes path/to/genomes/ \
  --ancestors path/to/WT_gene_fasta/ \
  --out results/
``
Arguments


--genomes
Directory containing genome FASTA files (.fasta, .fna, or .fna.gz)


--ancestors
Directory containing ancestral gene CDS FASTA files (.fasta or .fna)


--out
Output directory for all intermediate files and final results


All intermediate outputs (alignments, extracted CDS, deduplicated FASTAs)
are written inside the specified output directory.

## Output

The primary output file is:

mutation_final_evidence.tsv

This file provides a gene-level summary indicating whether each gene shows
evidence of mutation from any of the detection modules.
Intermediate files are retained to allow transparency, inspection,
and debugging.

## Methodology summary

At a high level, GENE-SCOUT performs the following steps:

Aligns ancestral gene CDS to genome assemblies using minimap2
Extracts CDS directly from genome sequences using alignment coordinates
Deduplicates identical CDS variants to improve scalability
Uses codon-aware alignment (MACSE) to detect frame-preserving mutations
Detects frameshifts directly from minimap2 CIGAR strings
Translates CDS to identify stop-gain and stop-loss mutations
Merges all evidence into a final mutation table


## License

GENE-SCOUT is released under the MIT License.
See the LICENSE file for details.

## Development status

GENE-SCOUT is under active development.
Future work will include formal validation, automated testing,
additional input safeguards, and performance optimisation.
