#!/bin/bash

################################################################################
# Example Usage Script
# Demonstrates various ways to use window_those_genomes.sh
################################################################################

echo "==================================="
echo "Genome Window Generator - Examples"
echo "==================================="
echo ""

# Example 1: Download human genome with default windows
echo "Example 1: Download human genome with default windows (1kb, 5kb, 10kb, 15kb, 20kb, 25kb)"
echo "Command: ./window_those_genomes.sh --species Homo_sapiens"
echo ""
read -p "Press Enter to continue..."
echo ""

# Example 2: Download mouse genome with custom windows
echo "Example 2: Download mouse genome with custom window sizes"
echo "Command: ./window_those_genomes.sh --species Mus_musculus --windows 1000,5000,10000"
echo ""
read -p "Press Enter to continue..."
echo ""

# Example 3: Use UCSC database
echo "Example 3: Download from UCSC instead of Ensembl"
echo "Command: ./window_those_genomes.sh --species Homo_sapiens --database ucsc"
echo ""
read -p "Press Enter to continue..."
echo ""

# Example 4: Process existing files
echo "Example 4: Process existing genome files without downloading"
echo "Command: ./window_those_genomes.sh --input ./genomes --output ./custom_output"
echo ""
read -p "Press Enter to continue..."
echo ""

# Example 5: Dry run
echo "Example 5: Preview commands without executing (dry-run mode)"
echo "Command: ./window_those_genomes.sh --species Homo_sapiens --dry-run"
echo ""
read -p "Press Enter to continue..."
echo ""

echo "==================================="
echo "Ready to run!"
echo "==================================="
echo ""
echo "To run any of these examples, uncomment the corresponding line below:"
echo ""

# Uncomment one of these to run:
# ./window_those_genomes.sh --species Homo_sapiens
# ./window_those_genomes.sh --species Mus_musculus --windows 1000,5000,10000
# ./window_those_genomes.sh --species Homo_sapiens --database ucsc
# ./window_those_genomes.sh --input ./genomes --output ./custom_output
# ./window_those_genomes.sh --species Homo_sapiens --dry-run

echo "Or run with --help to see all options:"
echo "./window_those_genomes.sh --help"
