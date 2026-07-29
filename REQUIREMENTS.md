# Conda environment

Create the workflow environment from the repository root:

```bash
conda env create -f environment.yml
conda activate prrsv-assembly
```

The environment provides Canu, fastp, NCBI BLAST+, minimap2, samtools, bcftools, seqkit and Python. The workflow discovers the active Conda environment automatically. If no environment is active, it attempts to locate the Conda base installation through `conda info --base`.

Verify that the required executables are available before processing samples:

```bash
command -v canu fastp blastn minimap2 samtools bcftools seqkit python3
```
