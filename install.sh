#!/bin/bash

################################################################################
# GenomeToWindows Installation Script
# Automatically sets up the required environment with all dependencies
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ENV_NAME="bedtools_env"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${BLUE}"
cat <<"EOF"
+===============================================================+
|                                                               |
|   +===+                                  __      __ _  _  _   |
|   | = | GENOME                          |  |    ||  \ || \ \ |
|   +===+          TO                     |/\| || |\|  || | |  |
|   [###]  [###]  [###]  [###]  [###]     |/\| || | |  ||_/ |  |
|    ~~~    ~~~    ~~~    ~~~    ~~~       \/  || |_/ |__/\_/   |
|                                               WINDOWS         |
+===============================================================+

            GENOME TO WINDOWS - INSTALLATION SCRIPT

EOF
echo -e "${NC}"

echo -e "${GREEN}Starting installation...${NC}\n"

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to detect package manager
detect_package_manager() {
    if command_exists micromamba; then
        echo "micromamba"
    elif command_exists mamba; then
        echo "mamba"
    elif command_exists conda; then
        echo "conda"
    else
        echo "none"
    fi
}

# Check for package manager
echo -e "${BLUE}[1/5] Checking for package manager...${NC}"
PKG_MANAGER=$(detect_package_manager)

if [ "$PKG_MANAGER" = "none" ]; then
    echo -e "${YELLOW}No conda/mamba/micromamba found. Installing micromamba...${NC}"

    # Install micromamba
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command_exists brew; then
            echo -e "${GREEN}Using Homebrew to install micromamba...${NC}"
            brew install micromamba
        else
            echo -e "${GREEN}Installing micromamba via shell script...${NC}"
            "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
        fi
    else
        # Linux
        echo -e "${GREEN}Installing micromamba via shell script...${NC}"
        "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
    fi

    # Re-detect after installation
    PKG_MANAGER=$(detect_package_manager)

    if [ "$PKG_MANAGER" = "none" ]; then
        echo -e "${RED}Failed to install micromamba. Please install manually:${NC}"
        echo -e "${YELLOW}Visit: https://mamba.readthedocs.io/en/latest/installation.html${NC}"
        exit 1
    fi
fi

echo -e "${GREEN}[OK] Found: $PKG_MANAGER${NC}\n"

# Initialize micromamba if needed
if [ "$PKG_MANAGER" = "micromamba" ]; then
    echo -e "${BLUE}[2/5] Initializing micromamba...${NC}"
    export MAMBA_ROOT_PREFIX="${HOME}/micromamba"

    # Check if already initialized
    if [ -f "${HOME}/.bashrc" ] && grep -q "micromamba shell hook" "${HOME}/.bashrc"; then
        echo -e "${GREEN}[OK] Micromamba already initialized${NC}\n"
    else
        eval "$(micromamba shell hook --shell bash)"
        echo -e "${GREEN}[OK] Micromamba initialized${NC}\n"
    fi
else
    echo -e "${BLUE}[2/5] Skipping initialization (using $PKG_MANAGER)${NC}\n"
fi

# Check if environment already exists
echo -e "${BLUE}[3/5] Checking for existing environment...${NC}"
ENV_EXISTS=false

if [ "$PKG_MANAGER" = "micromamba" ]; then
    if micromamba env list | grep -q "^${ENV_NAME} "; then
        ENV_EXISTS=true
    fi
elif [ "$PKG_MANAGER" = "mamba" ]; then
    if mamba env list | grep -q "^${ENV_NAME} "; then
        ENV_EXISTS=true
    fi
else
    if conda env list | grep -q "^${ENV_NAME} "; then
        ENV_EXISTS=true
    fi
fi

if [ "$ENV_EXISTS" = true ]; then
    echo -e "${YELLOW}Environment '${ENV_NAME}' already exists.${NC}"
    read -p "Do you want to recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}Removing existing environment...${NC}"
        $PKG_MANAGER env remove -n "$ENV_NAME" -y
        ENV_EXISTS=false
    else
        echo -e "${GREEN}[OK] Using existing environment${NC}\n"
    fi
fi

# Create environment if needed
if [ "$ENV_EXISTS" = false ]; then
    echo -e "${BLUE}[4/5] Creating conda environment '${ENV_NAME}'...${NC}"
    echo -e "${YELLOW}This may take a few minutes...${NC}"

    $PKG_MANAGER create -n "$ENV_NAME" -c conda-forge -c bioconda \
        bedtools \
        samtools \
        curl \
        -y

    echo -e "${GREEN}[OK] Environment created successfully${NC}\n"
else
    echo -e "${BLUE}[4/5] Environment already configured${NC}\n"
fi

# Make scripts executable
echo -e "${BLUE}[5/5] Making scripts executable...${NC}"
cd "$SCRIPT_DIR"
chmod +x window_those_genomes.sh
if [ -f "example_usage.sh" ]; then
    chmod +x example_usage.sh
fi
echo -e "${GREEN}[OK] Scripts are now executable${NC}\n"

# Final message
echo -e "${GREEN}==========================================================${NC}"
echo -e "${GREEN}Installation complete!${NC}"
echo -e "${GREEN}==========================================================${NC}\n"

echo -e "${BLUE}To use GenomeToWindows:${NC}\n"

if [ "$PKG_MANAGER" = "micromamba" ]; then
    echo -e "  ${YELLOW}1.${NC} Activate the environment:"
    echo -e "     ${GREEN}micromamba activate ${ENV_NAME}${NC}\n"
elif [ "$PKG_MANAGER" = "mamba" ]; then
    echo -e "  ${YELLOW}1.${NC} Activate the environment:"
    echo -e "     ${GREEN}mamba activate ${ENV_NAME}${NC}\n"
else
    echo -e "  ${YELLOW}1.${NC} Activate the environment:"
    echo -e "     ${GREEN}conda activate ${ENV_NAME}${NC}\n"
fi

echo -e "  ${YELLOW}2.${NC} Run the script:"
echo -e "     ${GREEN}./window_those_genomes.sh --species Homo_sapiens${NC}\n"

echo -e "  ${YELLOW}3.${NC} Get help:"
echo -e "     ${GREEN}./window_those_genomes.sh --help${NC}\n"

echo -e "${BLUE}Quick start example:${NC}"
echo -e "  ${GREEN}./window_those_genomes.sh --species Mus_musculus --windows 1000,5000,10000${NC}\n"

echo -e "${BLUE}Note:${NC} The script will automatically activate the environment when run."
echo -e ""
