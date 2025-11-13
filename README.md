# Genome Window Generator

A versatile Bash script that automatically downloads genomes from Ensembl or UCSC databases and generates genomic windows of specified sizes. This tool is part of the **MethylSense** package preprocessing pipeline.

## Features

- **Automatic genome downloading** from Ensembl or UCSC databases
- **Flexible window sizes** with sensible defaults (1kb, 5kb, 10kb, 15kb, 20kb, 25kb)
- **Automatic fallback** between databases if download fails
- **Progress tracking** with visual progress bars
- **Dry-run mode** to preview commands before execution
- **Automated environment setup** using micromamba/conda

## Requirements

- **Bash** (4.0+)
- **micromamba** or **conda** (for managing dependencies)
- **curl** (for downloading genomes)
- **bedtools** (installed automatically via micromamba)
- **samtools** (installed automatically via micromamba)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/GenomeToWindows.git
cd GenomeToWindows
```

2. Ensure micromamba is installed:
```bash
# If you don't have micromamba, install it:
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

3. Make the script executable:
```bash
chmod +x window_those_genomes.sh
```

## Usage

### Quick Start

Download human genome and generate default windows:
```bash
./window_those_genomes.sh --species Homo_sapiens
```

### Common Examples

**Download mouse genome with custom window sizes:**
```bash
./window_those_genomes.sh --species Mus_musculus --windows 1000,5000,10000
```

**Use UCSC database instead of Ensembl:**
```bash
./window_those_genomes.sh --species Homo_sapiens --database ucsc
```

**Process existing genome files:**
```bash
./window_those_genomes.sh --input ./my_genomes --output ./my_windows
```

**Preview commands without executing (dry-run):**
```bash
./window_those_genomes.sh --species Homo_sapiens --dry-run
```

### Command-Line Options

```
Options:
  -s, --species SPECIES  Scientific name of species (e.g., "Homo_sapiens")
                         When specified, genome will be automatically downloaded
  -d, --database DB      Database to use: ensembl or ucsc (default: ensembl)
  -i, --input DIR        Input directory with .fna/.fa files (default: ./genomes)
  -o, --output DIR       Output directory for BED files (default: ./windowed_genomes)
  -w, --windows SIZES    Comma-separated window sizes in bp (default: 1000,5000,10000,15000,20000,25000)
                         Example: -w 1000,5000,10000
  --dry-run              Preview commands without executing them
  -h, --help             Show this help message
```

## Supported Species

### Ensembl
Supports all species available in Ensembl (100+ species). The script automatically finds the latest assembly.

Common examples:
- `Homo_sapiens` (Human)
- `Mus_musculus` (Mouse)
- `Rattus_norvegicus` (Rat)
- `Danio_rerio` (Zebrafish)
- `Drosophila_melanogaster` (Fruit fly)
- `Caenorhabditis_elegans` (C. elegans)
- And many more...

### UCSC
Pre-configured species mapping:
- `Homo_sapiens` → hg38
- `Mus_musculus` → mm39
- `Rattus_norvegicus` → rn7
- `Danio_rerio` → danRer11
- `Drosophila_melanogaster` → dm6
- `Caenorhabditis_elegans` → ce11
- `Saccharomyces_cerevisiae` → sacCer3
- `Gallus_gallus` → galGal6
- `Sus_scrofa` → susScr11
- `Bos_taurus` → bosTau9
- `Canis_familiaris` → canFam6

## Output

The script generates:
1. **Genome files** in `./genomes/` (if downloaded)
2. **BED window files** in `./windowed_genomes/` with naming format:
   - `{genome_name}_{window_size}bp.bed`
   - Example: `Homo_sapiens_1000bp.bed`

### BED File Format
Standard BED format with 3 columns:
```
chromosome    start    end
chr1          0        1000
chr1          1000     2000
chr1          2000     3000
...
```

## How It Works

1. **Download Phase** (if `--species` is specified):
   - Connects to Ensembl or UCSC FTP servers
   - Finds the latest genome assembly
   - Downloads and extracts the genome
   - Falls back to alternate database if download fails

2. **Windowing Phase**:
   - Creates micromamba environment with bedtools and samtools
   - Indexes genome files using samtools
   - Generates genomic windows using bedtools
   - Creates BED files for each specified window size

## Troubleshooting

### Download fails for both databases
- Check your internet connection
- Verify species name spelling (use underscores: `Homo_sapiens`)
- Try specifying a different database with `--database`

### micromamba not found
Install micromamba:
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### Permission denied
Make the script executable:
```bash
chmod +x window_those_genomes.sh
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see LICENSE file for details

## Citation

If you use this tool in your research, please cite:

```
Drag, M. (2025). Genome Window Generator: Automated genome downloading and windowing tool.
Part of the MethylSense package. GitHub: https://github.com/yourusername/GenomeToWindows
```

## Related Projects

- **MethylSense**: Main package for methylation analysis (link TBA)

## Contact

For questions or issues, please open an issue on GitHub or contact:
- Markus Drag - markusdrag@gmail.com

## Acknowledgments

- Ensembl and UCSC Genome Browser for providing genome data
- bedtools and samtools development teams
