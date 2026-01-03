<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logo_dark_bg.png">
  <source media="(prefers-color-scheme: light)" srcset="logo_white_bg.png">
  <img alt="GenomeToWindows Logo" src="logo_white_bg.png" width="400">
</picture>

<p>
  <strong>Automated Genome Downloading & Windowing Tool</strong><br>
  <em>Part of the MethylSense Package</em>
</p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2025.04.11.648151-b31b1b.svg)](https://www.biorxiv.org/content/10.1101/2025.04.11.648151v1.full)
[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](https://github.com/markusdrag/GenomeToWindows)

</div>

---

## Overview

A versatile Bash script that generates genomic windows of specified sizes from FASTA files. Can also automatically download genomes from Ensembl or UCSC databases. This tool is part of the **MethylSense** package preprocessing pipeline.

## Features

- **Direct FASTA input** - Process local genome files with `--genome` flag (v1.1)
- **Automatic genome downloading** from Ensembl or UCSC databases
- **Flexible window sizes** with sensible defaults
- **macOS compatible** - Works on both Linux and macOS (v1.1)
- **System tools fallback** - Use system bedtools/samtools with `--no-micromamba` (v1.1)
- **Progress tracking** with visual progress bars
- **Dry-run mode** to preview commands before execution

## Requirements

- **Bash** (3.2+, macOS compatible)
- **bedtools** and **samtools** (via system install OR micromamba)
- **curl** (for downloading genomes)

## Installation

### Quick Install (Recommended)

```bash
git clone https://github.com/markusdrag/GenomeToWindows.git
cd GenomeToWindows
chmod +x window_those_genomes.sh
```

If you have bedtools and samtools installed, you're ready to go with `--no-micromamba`.

### With Micromamba (Optional)

For automatic dependency management:
```bash
./install.sh
```

## Usage

### Process a Local FASTA File (Recommended)

```bash
# Process your local genome file
./window_those_genomes.sh \
  --genome /path/to/genome.fna \
  --windows 200,400,500,750,1000 \
  --output ./my_windows \
  --no-micromamba
```

### Download and Process

```bash
# Download human genome from Ensembl and generate default windows
./window_those_genomes.sh --species Homo_sapiens

# Download chicken genome with benchmarking window sizes
./window_those_genomes.sh --species Gallus_gallus --windows 200,400,500,750,1000

# Download mouse genome with custom window sizes
./window_those_genomes.sh --species Mus_musculus --windows 1000,5000,10000

# Use UCSC database instead of Ensembl
./window_those_genomes.sh --species Homo_sapiens --database ucsc

# Preview download commands without executing (dry-run)
./window_those_genomes.sh --species Gallus_gallus --dry-run
```

### Process Directory of Genomes

```bash
./window_those_genomes.sh --input ./my_genomes --output ./my_windows
```

### Command-Line Options

```
Options:
  -g, --genome FILE    Direct path to a FASTA file (.fna/.fa/.fasta)
                       This is the PREFERRED option for local files
  --fai FILE           Direct path to a FAI index file (skips samtools indexing)
  -s, --species SPECIES  Scientific name (e.g., "Homo_sapiens") - downloads genome
  -d, --database DB      Database: ensembl or ucsc (default: ensembl)
  -i, --input DIR        Input directory with .fna/.fa files (default: ./genomes)
  -o, --output DIR       Output directory for BED files (default: ./windowed_genomes)
  -w, --windows SIZES    Comma-separated window sizes in bp
                         Example: -w 200,400,500,750,1000
  --no-micromamba        Use system bedtools/samtools instead of micromamba
  --dry-run              Preview commands without executing
  -h, --help             Show help message
  -v, --version          Show version
```

## Examples

### Quick Start: Download and Window a Genome

The easiest way to get started is to let the script automatically download a genome from Ensembl. Just specify the species name and your desired window sizes:

```bash
# Download human genome from Ensembl and generate 1kb, 5kb, and 10kb windows
./window_those_genomes.sh --species Homo_sapiens --windows 1000,5000,10000

# Download chicken genome with MethylSense benchmarking sizes (200bp to 1000bp)
./window_those_genomes.sh --species Gallus_gallus --windows 200,400,500,750,1000

# Download mouse genome from UCSC instead of Ensembl
./window_those_genomes.sh --species Mus_musculus --database ucsc
```

The script will:
1. Download the genome FASTA from the selected database
2. Index it with `samtools faidx`
3. Generate BED files with non-overlapping windows for each specified size

### Using a Local Genome File

If you already have a genome FASTA file on your system, use `--genome` to process it directly. This is faster and avoids re-downloading:

```bash
./window_those_genomes.sh \
  --genome /path/to/Gallus_gallus.GRCg7b.dna.toplevel.fna \
  --windows 200,400,500,750,1000 \
  --output ./benchmark_windows \
  --no-micromamba
```

> **Tip:** Use `--no-micromamba` if you have bedtools and samtools already installed on your system.

### Skipping Indexing with Pre-existing FAI

If you already have a `.fai` index file, you can skip the indexing step entirely by providing it with `--fai`:

```bash
./window_those_genomes.sh \
  --genome /path/to/genome.fna \
  --fai /path/to/genome.fna.fai \
  --windows 1000,5000,10000
```

This is useful when working with very large genomes where indexing can take several minutes.

## Supported Species (for download)

### Ensembl
All species in Ensembl (100+ species). Common examples:
- `Homo_sapiens`, `Mus_musculus`, `Rattus_norvegicus`
- `Danio_rerio`, `Gallus_gallus`, `Bos_taurus`

### UCSC
Pre-configured mapping:
- `Homo_sapiens` → hg38
- `Mus_musculus` → mm39
- `Gallus_gallus` → galGal6
- And more...

## Output

BED files in `./windowed_genomes/` (or custom `--output`):
- `{genome_name}_{window_size}bp.bed`
- Example: `Gallus_gallus_500bp.bed`

## Changelog

### v1.1.0 (2026-01-02)
- **NEW:** `--genome` flag for direct FASTA input
- **NEW:** `--fai` flag for direct FAI index input
- **NEW:** `--no-micromamba` to use system tools
- **FIX:** macOS compatibility (removed GNU grep -P)
- **FIX:** Bash 3.x compatibility (removed associative arrays)

### v1.0.0
- Initial release

## License

MIT License - see LICENSE file for details

## Citation

```
Drag, M., et al. (2025). New high accuracy diagnostics for avian Aspergillus fumigatus
infection using Nanopore methylation sequencing of host cell-free DNA and machine
learning prediction. bioRxiv. https://doi.org/10.1101/2025.04.11.648151
```

## Contact

- Issues: [GitHub Issues](https://github.com/markusdrag/GenomeToWindows/issues)
- Email: markusdrag@gmail.com
