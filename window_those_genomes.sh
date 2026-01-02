#!/bin/bash

################################################################################
# Genome Window Generator v1.1
# Part of the MethylSense package
#
# This script downloads genomes from Ensembl or UCSC and generates genomic
# windows of specified sizes using bedtools.
#
# v1.1 - Added --genome flag for direct FASTA file input
#      - Added --fai flag for direct FAI index input (skip samtools indexing)
#      - Fixed macOS compatibility (removed GNU grep -P)
#      - Added fallback to system bedtools/samtools if available
#
# Author: Markus Drag
# License: MIT
################################################################################

VERSION="1.1.0"

# === DEFAULTS ===
GENOME_DIR="./genomes"
OUTPUT_DIR="./windowed_genomes"
WINDOW_SIZES=(1000 5000 10000 15000 20000 25000)
DRY_RUN=false
ENV_NAME="bedtools_env"
SPECIES=""
DOWNLOAD_GENOME=false
DATABASE="ensembl"  # ensembl or ucsc
DIRECT_GENOME_FILE=""  # Direct path to a FASTA file
DIRECT_FAI_FILE=""     # Direct path to a FAI index file
USE_MICROMAMBA=true    # Whether to use micromamba environment

# Detect and set project root to script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

# === USAGE INFO ===
usage() {
    cat << EOF
Usage: $0 [options]

Genome Window Generator v${VERSION} - Part of the MethylSense package

This script can process existing genome files, use direct FASTA paths, or
automatically download genomes from Ensembl or UCSC databases.

Options:
  -g, --genome FILE    Direct path to a FASTA file (.fna/.fa/.fasta)
                       This is the PREFERRED option for local files
  --fai FILE           Direct path to a FAI index file (skips samtools indexing)
  -s, --species SPECIES  Scientific name of species (e.g., "Homo_sapiens")
                         When specified, genome will be automatically downloaded
  -d, --database DB      Database to use: ensembl or ucsc (default: ensembl)
  -i, --input DIR        Input directory with .fna/.fa files (default: ./genomes)
  -o, --output DIR       Output directory for BED files (default: ./windowed_genomes)
  -w, --windows SIZES    Comma-separated window sizes in bp (default: 1000,5000,10000,15000,20000,25000)
                         Example: -w 200,400,500,750,1000
  --no-micromamba        Use system bedtools/samtools instead of micromamba environment
  --dry-run              Preview commands without executing them
  -h, --help             Show this help message

Examples:
  # Process a specific FASTA file (RECOMMENDED for local files)
  $0 --genome /path/to/Gallus_gallus.GRCg7b.dna.toplevel.fna --windows 200,400,500,750,1000

  # Process existing genome files from a directory
  $0 --input ./my_genomes --output ./my_windows

  # Download human genome and create windows
  $0 --species Homo_sapiens

  # Download mouse genome with custom window sizes
  $0 --species Mus_musculus --windows 1000,5000,10000

  # Use UCSC database instead of Ensembl
  $0 --species Homo_sapiens --database ucsc

EOF
    exit 1
}

# === FLAG PARSING ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--genome)
            DIRECT_GENOME_FILE="$2"
            shift 2
            ;;
        --fai)
            DIRECT_FAI_FILE="$2"
            shift 2
            ;;
        -s|--species)
            SPECIES="$2"
            DOWNLOAD_GENOME=true
            shift 2
            ;;
        -d|--database)
            DATABASE="$2"
            shift 2
            ;;
        -i|--input)
            GENOME_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -w|--windows)
            IFS=',' read -ra WINDOW_SIZES <<< "$2"
            shift 2
            ;;
        --no-micromamba)
            USE_MICROMAMBA=false
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        -v|--version)
            echo "Genome Window Generator v${VERSION}"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate database option
if [[ "$DATABASE" != "ensembl" && "$DATABASE" != "ucsc" ]]; then
    echo "[ERROR] Database must be either 'ensembl' or 'ucsc'"
    exit 1
fi

# === GENOME DOWNLOAD FUNCTIONS ===
download_from_ensembl() {
    local species=$1
    local genome_dir=$2

    echo "[NET] Downloading genome for $species from Ensembl..."

    # Create genome directory if it doesn't exist
    mkdir -p "$genome_dir"

    # Convert species name format (e.g., Homo_sapiens to homo_sapiens for URL)
    local species_lower=$(echo "$species" | tr '[:upper:]' '[:lower:]')

    # Try to get the latest release number (macOS compatible - no grep -P)
    echo "[SEARCH] Finding latest Ensembl release..."
    local latest_release=$(curl -s "https://ftp.ensembl.org/pub/" | grep -o 'release-[0-9]*' | sed 's/release-//' | sort -n | tail -1)

    if [ -z "$latest_release" ]; then
        echo "[WARN]  Could not determine latest release, using current_fasta"
        local base_url="https://ftp.ensembl.org/pub/current_fasta"
    else
        echo "[PKG] Using Ensembl release $latest_release"
        local base_url="https://ftp.ensembl.org/pub/release-${latest_release}/fasta"
    fi

    local fasta_url="${base_url}/${species_lower}/dna/"

    echo "[SEARCH] Searching for primary assembly or toplevel DNA file..."

    # Try to find the primary assembly file first, then toplevel (macOS compatible)
    local genome_file=$(curl -s "$fasta_url" | grep -o "${species}[^\"]*\.dna\.\(primary_assembly\|toplevel\)\.fa\.gz" | head -1)

    if [ -z "$genome_file" ]; then
        echo "[ERROR] Could not find genome file for $species at Ensembl"
        echo "   Tried URL: $fasta_url"
        return 1
    fi

    local download_url="${fasta_url}${genome_file}"
    local output_file="${genome_dir}/${genome_file}"

    echo "[DOWNLOAD]  Downloading: $genome_file"
    echo "   From: $download_url"

    if $DRY_RUN; then
        echo "[dry-run] curl -L -o \"$output_file\" \"$download_url\""
        echo "[dry-run] gunzip \"$output_file\""
    else
        curl -L -o "$output_file" "$download_url" || {
            echo "[ERROR] Download failed"
            return 1
        }
        echo "[PKG] Extracting genome..."
        gunzip "$output_file"
    fi

    echo "[OK] Genome downloaded successfully"
    return 0
}

download_from_ucsc() {
    local species=$1
    local genome_dir=$2

    echo "[NET] Downloading genome for $species from UCSC..."

    # Create genome directory if it doesn't exist
    mkdir -p "$genome_dir"

    # Common UCSC genome database names (bash 3.x compatible - no associative arrays)
    local ucsc_db=""
    case "$species" in
        "Homo_sapiens") ucsc_db="hg38" ;;
        "Mus_musculus") ucsc_db="mm39" ;;
        "Rattus_norvegicus") ucsc_db="rn7" ;;
        "Danio_rerio") ucsc_db="danRer11" ;;
        "Drosophila_melanogaster") ucsc_db="dm6" ;;
        "Caenorhabditis_elegans") ucsc_db="ce11" ;;
        "Saccharomyces_cerevisiae") ucsc_db="sacCer3" ;;
        "Gallus_gallus") ucsc_db="galGal6" ;;
        "Sus_scrofa") ucsc_db="susScr11" ;;
        "Bos_taurus") ucsc_db="bosTau9" ;;
        "Canis_familiaris") ucsc_db="canFam6" ;;
        *)
            echo "[ERROR] Species $species not found in UCSC mapping"
            echo "   Supported species: Homo_sapiens, Mus_musculus, Rattus_norvegicus,"
            echo "   Danio_rerio, Drosophila_melanogaster, Caenorhabditis_elegans,"
            echo "   Saccharomyces_cerevisiae, Gallus_gallus, Sus_scrofa, Bos_taurus, Canis_familiaris"
            return 1
            ;;
    esac

    local base_url="https://hgdownload.soe.ucsc.edu/goldenPath/${ucsc_db}/bigZips"
    local genome_file="${ucsc_db}.fa.gz"
    local download_url="${base_url}/${genome_file}"
    local output_file="${genome_dir}/${genome_file}"

    echo "[DOWNLOAD]  Downloading: $genome_file (UCSC: $ucsc_db)"
    echo "   From: $download_url"

    if $DRY_RUN; then
        echo "[dry-run] curl -L -o \"$output_file\" \"$download_url\""
        echo "[dry-run] gunzip \"$output_file\""
    else
        curl -L -o "$output_file" "$download_url" || {
            echo "[ERROR] Download failed"
            return 1
        }
        echo "[PKG] Extracting genome..."
        gunzip "$output_file"
    fi

    echo "[OK] Genome downloaded successfully"
    return 0
}

# === DOWNLOAD GENOME IF REQUESTED ===
if $DOWNLOAD_GENOME; then
    if [ -z "$SPECIES" ]; then
        echo "[ERROR] Error: Species name required when downloading genome"
        echo "   Use -s or --species flag"
        exit 1
    fi

    echo "=================================="
    echo "Genome Download"
    echo "Species: $SPECIES"
    echo "Database: $DATABASE"
    echo "Output: $GENOME_DIR"
    echo "=================================="

    if [ "$DATABASE" == "ensembl" ]; then
        download_from_ensembl "$SPECIES" "$GENOME_DIR"
        if [ $? -ne 0 ]; then
            echo ""
            echo "[WARN]  Ensembl download failed. Trying UCSC as fallback..."
            download_from_ucsc "$SPECIES" "$GENOME_DIR" || {
                echo "[ERROR] Both Ensembl and UCSC downloads failed"
                exit 1
            }
        fi
    else
        download_from_ucsc "$SPECIES" "$GENOME_DIR"
        if [ $? -ne 0 ]; then
            echo ""
            echo "[WARN]  UCSC download failed. Trying Ensembl as fallback..."
            download_from_ensembl "$SPECIES" "$GENOME_DIR" || {
                echo "[ERROR] Both UCSC and Ensembl downloads failed"
                exit 1
            }
        fi
    fi

    echo ""
fi

# === TOOL SETUP ===
cd "$PROJECT_ROOT" || { echo "[ERROR] Could not cd to project root: $PROJECT_ROOT"; exit 1; }

# Check for bedtools and samtools availability
check_tools() {
    local has_bedtools=false
    local has_samtools=false
    
    if command -v bedtools &> /dev/null; then
        has_bedtools=true
    fi
    if command -v samtools &> /dev/null; then
        has_samtools=true
    fi
    
    if $has_bedtools && $has_samtools; then
        echo "[OK] Found system bedtools and samtools"
        return 0
    fi
    return 1
}

# Try system tools first if --no-micromamba or if micromamba not available
if ! $USE_MICROMAMBA || ! command -v micromamba &> /dev/null; then
    if check_tools; then
        echo "[INFO] Using system bedtools and samtools"
        USE_MICROMAMBA=false
    else
        if ! command -v micromamba &> /dev/null; then
            echo "[WARN] micromamba not found and system bedtools/samtools not available"
            echo "   To install micromamba: curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar xvj"
            echo "   Or install bedtools and samtools via your package manager"
            exit 1
        fi
    fi
fi

# Setup micromamba environment if needed
if $USE_MICROMAMBA; then
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
    
    if command -v micromamba &> /dev/null; then
        eval "$(micromamba shell hook --shell bash)"
    else
        echo "[ERROR] micromamba not found. Use --no-micromamba to use system tools."
        exit 1
    fi

    if ! $DRY_RUN; then
        if ! micromamba env list | grep -q "^$ENV_NAME "; then
            echo "[SETUP] Creating micromamba environment '$ENV_NAME'..."
            micromamba create -n "$ENV_NAME" -c conda-forge -c bioconda bedtools samtools -y || {
                echo "[ERROR] Failed to create micromamba environment: $ENV_NAME"
                exit 1
            }
        fi
        echo "[SETUP] Activating micromamba environment '$ENV_NAME'..."
        micromamba activate "$ENV_NAME" || {
            echo "[ERROR] Failed to activate micromamba environment: $ENV_NAME"
            echo "Try running: micromamba create -n $ENV_NAME -c conda-forge -c bioconda bedtools samtools -y"
            exit 1
        }
    else
        echo "[dry-run] micromamba activate $ENV_NAME"
    fi
fi

# === FUNCTIONS ===
run_cmd() {
    if $DRY_RUN; then
        echo "[dry-run] $*"
    else
        eval "$@"
    fi
}

progress_bar() {
    local progress=$1
    local total=$2
    local percent=$((progress * 100 / total))
    local bar_length=$((percent / 2))
    local spaces=$((50 - bar_length))
    printf "\rProgress: [%-*s%*s] %d%%" "$bar_length" "#" "$spaces" "" "$percent"
}

# === MAIN ===
mkdir -p "$OUTPUT_DIR"

echo ""
echo "=================================="
echo "Genome Window Generator v${VERSION}"
echo "=================================="

# Collect FASTA files to process
GENOMES=()

if [ -n "$DIRECT_GENOME_FILE" ]; then
    # Direct file mode
    if [ ! -f "$DIRECT_GENOME_FILE" ]; then
        echo "[ERROR] Genome file not found: $DIRECT_GENOME_FILE"
        exit 1
    fi
    GENOMES=("$DIRECT_GENOME_FILE")
    echo "Mode: Direct FASTA file"
    echo "Genome file: $DIRECT_GENOME_FILE"
else
    # Directory scanning mode
    echo "Mode: Directory scan"
    echo "Input directory: $GENOME_DIR"
    echo "[SEARCH] Scanning for genome files in: $GENOME_DIR"
    
    # Find genome files (macOS compatible)
    while IFS= read -r file; do
        GENOMES+=("$file")
    done < <(find "$GENOME_DIR" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) 2>/dev/null)
fi

echo "Output directory: $OUTPUT_DIR"
echo "Window sizes (bp): ${WINDOW_SIZES[*]}"
echo "=================================="
echo ""

if [[ ${#GENOMES[@]} -eq 0 ]]; then
    echo "[ERROR] No genome files found"
    echo "   Use --genome to specify a FASTA file directly"
    echo "   Or use --input to specify a directory with .fna/.fa/.fasta files"
    exit 1
fi

echo "[OK] Found ${#GENOMES[@]} genome file(s) to process"
echo ""

TOTAL=${#GENOMES[@]}
COUNT=0

for GENOME in "${GENOMES[@]}"; do
    ((COUNT++))
    
    BASENAME=$(basename "$GENOME")
    BASENAME=${BASENAME%.fna}
    BASENAME=${BASENAME%.fa}
    BASENAME=${BASENAME%.fasta}
    BASENAME=${BASENAME%.gz}

    echo "[PROCESS] Processing: $BASENAME ($COUNT/$TOTAL)"

    GENOME_INDEX="${GENOME}.fai"
    BEDGENOME="${OUTPUT_DIR}/${BASENAME}.genome"

    # Use direct FAI if provided, otherwise check for existing or create
    if [ -n "$DIRECT_FAI_FILE" ]; then
        if [ ! -f "$DIRECT_FAI_FILE" ]; then
            echo "[ERROR] FAI index file not found: $DIRECT_FAI_FILE"
            exit 1
        fi
        GENOME_INDEX="$DIRECT_FAI_FILE"
        echo "   Using provided FAI index: $GENOME_INDEX"
    elif [ ! -f "$GENOME_INDEX" ]; then
        echo "   Creating FASTA index..."
        run_cmd "samtools faidx \"$GENOME\""
    else
        echo "   Using existing FAI index: $GENOME_INDEX"
    fi

    # Create genome length file for bedtools
    echo "   Creating genome length file..."
    run_cmd "cut -f1,2 \"$GENOME_INDEX\" > \"$BEDGENOME\""

    # Generate BED windows for each size
    for WIN in "${WINDOW_SIZES[@]}"; do
        OUTFILE="${OUTPUT_DIR}/${BASENAME}_${WIN}bp.bed"
        echo "   Generating ${WIN}bp windows..."
        run_cmd "bedtools makewindows -g \"$BEDGENOME\" -w $WIN > \"$OUTFILE\""
    done
    
    echo "   [OK] Done with $BASENAME"
    echo ""
done

echo "=================================="
echo "[OK] All done! Windowed BED files written to: $OUTPUT_DIR"
echo ""
echo "Generated files:"
ls -la "$OUTPUT_DIR"/*.bed 2>/dev/null | head -20
echo "=================================="
