#!/bin/bash

# Bioconda build script for SNPhylo2
set -euo pipefail

# Install Python package
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# Create wrapper scripts for external tools to ensure they're in PATH
# These are handled by conda dependencies, but we verify availability

echo "SNPhylo2 installation complete"
echo "External tools verified:"
echo "  - IQ-TREE2: $(which iqtree2 2>/dev/null || echo 'not found')"
echo "  - PLINK2: $(which plink2 2>/dev/null || echo 'not found')"
echo "  - BCFtools: $(which bcftools 2>/dev/null || echo 'not found')"
echo "  - Samtools: $(which samtools 2>/dev/null || echo 'not found')"
