# GENE-SCOUT

GENE-SCOUT is a gene-centric bioinformatics workflow for detecting
coding-sequence (CDS) mutations across genome assemblies using
nucleotide-resolved evidence.

⚠️ **Software status**  
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
- frame-preserving mutations (e.g. missense, in-frame indels)
- frame-disrupting mutations (frameshifts and stop codon changes)

---

## Mutation classes detected

GENE-SCOUT detects the following mutation types:

- **Missense substitutions**
- **In-frame insertions and deletions**
- **Frameshift mutations**
- **Stop-gain mutations** (premature termination)
- **Stop-loss mutations** (loss of the ancestral stop codon)

Each mutation class is detected using evidence appropriate to its
underlying molecular mechanism.

---

## Input requirements

### Genome assemblies (`--genomes`)

Directory containing genome assemblies as **nucleotide FASTA files**.

Supported formats:
- `.fasta`
- `.fna`
- `.fna.gz`

Assemblies may contain one or multiple contigs.
GENE-SCOUT does not assume any specific genome source.

---

### Ancestor genes (`--ancestors`)

Directory containing **ancestral gene coding sequences (CDS only)**,
provided as **nucleotide FASTA files**.

Supported formats:
- `.fasta`
- `.fna`

❗ Protein sequences are **not supported**.

Each ancestor FASTA file must contain a single CDS corresponding to one
gene of interest and represent the coding region only
(no UTRs or introns).

---

## Installation

Clone the repository:

```bash
git clone https://github.com/alhubb/gene-scout
cd gene-scout
```
---

## Create and activate the Conda environment
```
conda env create -f environment.yml
conda activate gene-scout
```
---

## Make the main command executable
```
chmod +x gene-scout
```
## Useage
```
gene-scout \
  --genomes path/to/genomes/ \
  --ancestors path/to/WT_gene_fasta/ \
  --out results/
```

### Screening for predefined mutations

To screen a genome dataset for a user-defined list of mutations
(e.g. derived from Snippy), provide a mutation list:

```bash
gene-scout \
  --genomes genomes/ \
  --ancestors WT_gene_fasta/ \
  --mutations mutations.tsv \
  --out results/
```
---

## Arguments
--genomes
Directory containing genome FASTA files in nucleotide format.
Supported extensions:

.fasta
.fna
.fna.gz


--ancestors
Directory containing ancestral gene coding sequences (CDS only) in
nucleotide FASTA format.
Supported extensions:

.fasta
.fna

Protein sequences are not supported.

--out
Output directory where all intermediate files and final results are written.
All intermediate outputs (including minimap2 alignments, extracted CDS,
and deduplicated FASTAs) are retained to allow transparency,
inspection, and debugging.

---

## Output
The primary output file is:

mutation_final_evidence.tsv

This file provides a gene-level summary indicating whether each gene shows
evidence of mutation from any of the detection modules
(codon-aware alignment, frameshift detection, or stop codon analysis).
Intermediate files are retained by default to support inspection
and reproducibility.

---

## Methodology Summary
At a high level, GENE-SCOUT performs the following steps:

Aligns ancestral gene CDS to genome assemblies using minimap2
Extracts CDS directly from genome sequences using alignment coordinates
Deduplicates identical CDS variants to improve scalability
Uses codon-aware alignment (MACSE) to detect frame-preserving mutations
Detects frameshifts directly from minimap2 CIGAR strings
Translates CDS to identify stop-gain and stop-loss mutations
Merges all evidence into a final mutation table

---

## License
GENE-SCOUT is released under the MIT License.
See the LICENSE file for details.

---

## Development status
GENE-SCOUT is under active development.
Future work will include formal validation, automated testing,
additional input safeguards, and performance optimisation.
