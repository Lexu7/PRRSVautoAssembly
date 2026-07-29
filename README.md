# PRRSVAssembly

Nanopore tiled-amplicon assembly workflow for PRRSV. Canu seed contigs are matched to the bundled full-genome PRRSV BLAST database, and the candidate with the highest cumulative BLAST bit score is used for two rounds of read-backed consensus generation.

## Reference selection

Only full-genome BLAST results are used. For each candidate reference, bit scores from all megablast hits against the Canu contigs are summed. The reference with the largest cumulative score is selected; score ties are resolved by reference ID in ascending lexical order. The complete ranked score table is written to `work/<sample>.best_ref_scores.tsv`.

When a diagnostic best-reference-filled FASTA is produced, zero-depth consensus positions are filled only from homologous reference coordinates established by a minimap2 alignment. Zero-depth inserted or unmapped positions are written as `N`; this diagnostic file is not the primary final consensus.

## Contents

- `PRRSVAssembly.sh`: main pipeline.
- `scripts/`: required Python helpers.
- `references/prrsv_refs.fasta`: reference sequences.
- `references/prrsv_blast.*`: prebuilt BLAST database for the bundled reference FASTA.

## Requirements

- Canu 2.2
- fastp
- NCBI BLAST+ (`blastn`)
- minimap2
- samtools
- bcftools
- seqkit
- Python 3

## Conda environment

Create and activate the workflow environment from the repository root:

```bash
conda env create -f environment.yml
conda activate prrsv-assembly
```

The environment provides Canu, fastp, NCBI BLAST+, minimap2, samtools, bcftools, seqkit and Python. The workflow discovers the active Conda environment automatically. If no environment is active, it attempts to locate the Conda base installation through `conda info --base`.

Verify that the required executables are available before processing samples:

```bash
command -v canu fastp blastn minimap2 samtools bcftools seqkit python3
```

## Run

```bash
bash PRRSVAssembly.sh -t 16 -o BC01_prrsv_assembly BC01.fastq.gz
```

Use `bash PRRSVAssembly.sh --help` for optional parameters. The packaged reference FASTA and BLAST database are used by default; override them with `--ref-fasta` and `--db-full` when necessary.

The pipeline refuses to reuse a non-empty output directory by default. Use `--resume` only when intentionally continuing an existing run with the same input and parameters.
