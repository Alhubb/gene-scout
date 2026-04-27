# GENE-SCOUT

GENE-SCOUT is a gene-centric pipeline for detecting coding-sequence
mutations across genome assemblies.

## Inputs

### Genomes
- `.fasta`
- `.fna`
- `.fna.gz`

### Ancestor genes
Nucleotide CDS ONLY:
- `.fasta`
- `.fna`

Protein sequences are not supported.

## Usage

```bash
gene-scout \
  --genomes genomes/ \
  --ancestors WT_gene_fasta/ \
  --out results/
