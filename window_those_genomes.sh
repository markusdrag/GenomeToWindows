#!/bin/bash

################################################################################
# Genome Window Generator
# Part of the MethylSense package
#
# This script downloads genomes from Ensembl or UCSC and generates genomic
# windows of specified sizes using bedtools.
#
# Author: Markus Drag
# License: MIT
################################################################################

# === DEFAULTS ===
GENOME_DIR="./genomes"
OUTPUT_DIR="./windowed_genomes"
WINDOW_SIZES=(1000 5000 10000 15000 20000 25000)
DRY_RUN=false
ENV_NAME="bedtools_env"
SPECIES=""
DOWNLOAD_GENOME=false
DATABASE="ensembl"  # ensembl or ucsc

# Detect and set project root to script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

# === USAGE INFO ===
usage() {
    cat << EOF
Usage: $0 [options]

Genome Window Generator - Part of the MethylSense package

This script can either process existing genome files or automatically download
genomes from Ensembl or UCSC databases and generate genomic windows.

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

Examples:
  # Download human genome and create windows
  $0 --species Homo_sapiens

  # Download mouse genome with custom window sizes
  $0 --species Mus_musculus --windows 1000,5000,10000

  # Process existing genome files
  $0 --input ./my_genomes --output ./my_windows

  # Use UCSC database instead of Ensembl
  $0 --species Homo_sapiens --database ucsc

EOF
    exit 1
}

# === FLAG PARSING ===
while [[ $# -gt 0 ]]; do
    case "$1" in
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
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate database option
if [[ "$DATABASE" != "ensembl" && "$DATABASE" != "ucsc" ]]; then
    echo "❌ Error: Database must be either 'ensembl' or 'ucsc'"
    exit 1
fi

# === GENOME DOWNLOAD FUNCTIONS ===
download_from_ensembl() {
    local species=$1
    local genome_dir=$2

    echo "🌐 Downloading genome for $species from Ensembl..."

    # Create genome directory if it doesn't exist
    mkdir -p "$genome_dir"

    # Convert species name format (e.g., Homo_sapiens to homo_sapiens for URL)
    local species_lower=$(echo "$species" | tr '[:upper:]' '[:lower:]')

    # Try to get the latest release number
    echo "🔍 Finding latest Ensembl release..."
    local latest_release=$(curl -s "https://ftp.ensembl.org/pub/" | grep -oP 'release-\K[0-9]+' | sort -n | tail -1)

    if [ -z "$latest_release" ]; then
        echo "⚠️  Could not determine latest release, using current_fasta"
        local base_url="https://ftp.ensembl.org/pub/current_fasta"
    else
        echo "📦 Using Ensembl release $latest_release"
        local base_url="https://ftp.ensembl.org/pub/release-${latest_release}/fasta"
    fi

    local fasta_url="${base_url}/${species_lower}/dna/"

    echo "🔍 Searching for primary assembly or toplevel DNA file..."

    # Try to find the primary assembly file first, then toplevel
    local genome_file=$(curl -s "$fasta_url" | grep -oP "${species}[^\"]*\.dna\.(primary_assembly|toplevel)\.fa\.gz" | head -1)

    if [ -z "$genome_file" ]; then
        echo "❌ Could not find genome file for $species at Ensembl"
        echo "   Tried URL: $fasta_url"
        return 1
    fi

    local download_url="${fasta_url}${genome_file}"
    local output_file="${genome_dir}/${genome_file}"

    echo "⬇️  Downloading: $genome_file"
    echo "   From: $download_url"

    if $DRY_RUN; then
        echo "[dry-run] curl -L -o \"$output_file\" \"$download_url\""
        echo "[dry-run] gunzip \"$output_file\""
    else
        curl -L -o "$output_file" "$download_url" || {
            echo "❌ Download failed"
            return 1
        }
        echo "📦 Extracting genome..."
        gunzip "$output_file"
    fi

    echo "✅ Genome downloaded successfully"
    return 0
}

download_from_ucsc() {
    local species=$1
    local genome_dir=$2

    echo "🌐 Downloading genome for $species from UCSC..."

    # Create genome directory if it doesn't exist
    mkdir -p "$genome_dir"

    # Common UCSC genome database names
    declare -A ucsc_names=(
        ["Homo_sapiens"]="hg38"
        ["Mus_musculus"]="mm39"
        ["Rattus_norvegicus"]="rn7"
        ["Danio_rerio"]="danRer11"
        ["Drosophila_melanogaster"]="dm6"
        ["Caenorhabditis_elegans"]="ce11"
        ["Saccharomyces_cerevisiae"]="sacCer3"
        ["Gallus_gallus"]="galGal6"
        ["Sus_scrofa"]="susScr11"
        ["Bos_taurus"]="bosTau9"
        ["Canis_familiaris"]="canFam6"
    )

    local ucsc_db="${ucsc_names[$species]}"

    if [ -z "$ucsc_db" ]; then
        echo "❌ Species $species not found in UCSC mapping"
        echo "   Please specify one of: ${!ucsc_names[@]}"
        return 1
    fi

    local base_url="https://hgdownload.soe.ucsc.edu/goldenPath/${ucsc_db}/bigZips"
    local genome_file="${ucsc_db}.fa.gz"
    local download_url="${base_url}/${genome_file}"
    local output_file="${genome_dir}/${genome_file}"

    echo "⬇️  Downloading: $genome_file (UCSC: $ucsc_db)"
    echo "   From: $download_url"

    if $DRY_RUN; then
        echo "[dry-run] curl -L -o \"$output_file\" \"$download_url\""
        echo "[dry-run] gunzip \"$output_file\""
    else
        curl -L -o "$output_file" "$download_url" || {
            echo "❌ Download failed"
            return 1
        }
        echo "📦 Extracting genome..."
        gunzip "$output_file"
    fi

    echo "✅ Genome downloaded successfully"
    return 0
}

# === DOWNLOAD GENOME IF REQUESTED ===
if $DOWNLOAD_GENOME; then
    if [ -z "$SPECIES" ]; then
        echo "❌ Error: Species name required when downloading genome"
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
            echo "⚠️  Ensembl download failed. Trying UCSC as fallback..."
            download_from_ucsc "$SPECIES" "$GENOME_DIR" || {
                echo "❌ Both Ensembl and UCSC downloads failed"
                exit 1
            }
        fi
    else
        download_from_ucsc "$SPECIES" "$GENOME_DIR"
        if [ $? -ne 0 ]; then
            echo ""
            echo "⚠️  UCSC download failed. Trying Ensembl as fallback..."
            download_from_ensembl "$SPECIES" "$GENOME_DIR" || {
                echo "❌ Both UCSC and Ensembl downloads failed"
                exit 1
            }
        fi
    fi

    echo ""
fi

# === MICROMAMBA SETUP ===
cd "$PROJECT_ROOT" || { echo "❌ Could not cd to project root: $PROJECT_ROOT"; exit 1; }

# Ensure micromamba is initialized
export MAMBA_ROOT_PREFIX="$HOME/micromamba"

# Initialize micromamba for this shell session
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
else
    echo "❌ micromamba not found. Please install micromamba first."
    exit 1
fi

# Check if environment exists, create if it doesn't
if ! $DRY_RUN; then
    if ! micromamba env list | grep -q "^$ENV_NAME "; then
        echo "🔧 Creating micromamba environment '$ENV_NAME'..."
        micromamba create -n "$ENV_NAME" -c conda-forge -c bioconda bedtools samtools -y || {
            echo "❌ Failed to create micromamba environment: $ENV_NAME"
            exit 1
        }
    fi
    echo "🔧 Activating micromamba environment '$ENV_NAME'..."
    micromamba activate "$ENV_NAME" || {
        echo "❌ Failed to activate micromamba environment: $ENV_NAME"
        echo "Try running: micromamba create -n $ENV_NAME -c conda-forge -c bioconda bedtools samtools -y"
        exit 1
    }
else
    echo "[dry-run] micromamba activate $ENV_NAME"
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
echo "Genome Window Generation"
echo "Input directory: $GENOME_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Window sizes (bp): ${WINDOW_SIZES[*]}"
echo "=================================="
echo ""

echo "🔍 Scanning for genome files in: $GENOME_DIR"
mapfile -t GENOMES < <(find "$GENOME_DIR" -type f \( -name "*.fna" -o -name "*.fa" \))

if [[ ${#GENOMES[@]} -eq 0 ]]; then
    echo "❌ No genome files found in $GENOME_DIR"
    exit 1
fi

TOTAL=${#GENOMES[@]}
COUNT=0

for GENOME in "${GENOMES[@]}"; do
    ((COUNT++))
    progress_bar "$COUNT" "$TOTAL"

    BASENAME=$(basename "$GENOME")
    BASENAME=${BASENAME%.fna}
    BASENAME=${BASENAME%.fa}

    GENOME_INDEX="${GENOME}.fai"
    BEDGENOME="${OUTPUT_DIR}/${BASENAME}.genome"

    # Index if needed
    if [ ! -f "$GENOME_INDEX" ]; then
        run_cmd "samtools faidx \"$GENOME\""
    fi

    # Create genome length file
    run_cmd "cut -f1,2 \"$GENOME_INDEX\" > \"$BEDGENOME\""

    # Generate BED windows
    for WIN in "${WINDOW_SIZES[@]}"; do
        OUTFILE="${OUTPUT_DIR}/${BASENAME}_${WIN}bp.bed"
        run_cmd "bedtools makewindows -g \"$BEDGENOME\" -w $WIN > \"$OUTFILE\""
    done
done

echo -e "\n✅ Done. Windowed BED files written to: $OUTPUT_DIR"
