# GENE-SCOUT

GENE-SCOUT is a gene-centric pipeline for detecting coding-sequence (CDS)
mutations across genome assemblies using nucleotide-resolved evidence.

⚠️ **Status**: This software is provided as **v0.1** and has **not yet been formally validated**.
Results should be interpreted with caution.

---

## What GENE-SCOUT does

GENE-SCOUT performs gene-level mutation analysis by:

- Aligning ancestral gene CDS to genome assemblies using **minimap2**
- Extracting CDS directly from genome sequences using PAF coordinates
- Identifying:
  - missense mutations
  - in-frame insertions/deletions
  - frameshifts
  - stop-gain and stop-loss mutations
- Integrating multiple evidence streams into a final mutation table

The pipeline is designed for population-scale analyses of bacterial genomes.

---

## Input requirements

### Genome assemblies (`--genomes`)

A directory containing genome FASTA files in **nucleotide format**:

- `.fasta`
- `.fna`
- `.fna.gz`

Assemblies may contain one or multiple contigs.

### Ancestor genes (`--ancestors`)

A directory containing **nucleotide coding sequences (CDS only)** for ancestral genes:

- `.fasta`
- `.fna`

❗ Protein sequences are **not supported**.

Each file should contain a single CDS corresponding to a gene of interest.

---

## Installation

GENE-SCOUT is distributed as a set of Bash and Python scripts and requires Conda.

```bash
git clone https://github.com/alhubb/gene-scout
cd gene-scout
conda env create -f environment.yml
conda activate gene-scout

