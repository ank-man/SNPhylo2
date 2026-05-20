# Installation Guide

This guide covers various methods for installing SNPhylo2 on your system.

## System Requirements

### Minimum Requirements
- Python 3.10 or higher
- 4 GB RAM
- 10 GB free disk space
- Linux or macOS operating system

### Recommended Requirements
- Python 3.11 or 3.12
- 32 GB RAM or more
- SSD storage
- 16+ CPU cores
- 100 GB free disk space for large datasets

### External Dependencies
SNPhylo2 requires several external bioinformatics tools:

- **IQ-TREE2** (≥2.3.0) - Phylogenetic tree inference
- **PLINK2** (≥2.0) - Genotype processing and LD pruning
- **BCFtools** (≥1.19) - VCF/BCF manipulation
- **Samtools** (≥1.19) - Sequence data processing

## Installation Methods

### Option 1: Conda (Recommended)

The easiest way to install SNPhylo2 is through Conda, which handles both Python dependencies and external tools.

```bash
# Install from Bioconda
conda install -c bioconda -c conda-forge snphylo2

# Verify installation
snphylo2 --version
snphylo2 validate  # Check external tools
```

### Option 2: Docker

For reproducible environments or systems where you don't have administrator access:

```bash
# Pull the Docker image
docker pull ghcr.io/snphylo2/snphylo2:latest

# Run SNPhylo2
docker run -v $(pwd):/data ghcr.io/snphylo2/snphylo2:latest \
    run -v /data/input.vcf.gz -o /data/results

# Or run interactively
docker run -it -v $(pwd):/data ghcr.io/snphylo2/snphylo2:latest bash
```

### Option 3: Singularity/Apptainer (HPC)

For HPC environments that use Singularity/Apptainer:

```bash
# Pull the image
singularity pull docker://ghcr.io/snphylo2/snphylo2:latest

# Run
singularity run snphylo2_latest.sif run -v input.vcf.gz -o results/
```

### Option 4: From Source

For development or to use the latest features:

```bash
# Clone the repository
git clone https://github.com/snphylo2/snphylo2.git
cd snphylo2

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Install external tools separately
# See external tool installation section below
```

### Option 5: pip + System Package Manager

If you prefer to manage external tools through your system package manager:

```bash
# Install external tools first
# Ubuntu/Debian:
sudo apt-get install plink2 bcftools samtools

# macOS with Homebrew:
brew install plink2 bcftools samtools

# Install IQ-TREE2 (download from GitHub releases)
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.4/iqtree-2.3.4-Linux.tar.gz
tar -xzf iqtree-2.3.4-Linux.tar.gz
sudo mv iqtree-2.3.4-Linux/bin/iqtree2 /usr/local/bin/

# Then install SNPhylo2
pip install snphylo2
```

## External Tool Installation

### IQ-TREE2

Download from [GitHub releases](https://github.com/iqtree/iqtree2/releases):

```bash
# Linux
wget https://github.com/iqtree/iqtree2/releases/download/v2.3.4/iqtree-2.3.4-Linux.tar.gz
tar -xzf iqtree-2.3.4-Linux.tar.gz
sudo mv iqtree-2.3.4-Linux/bin/iqtree2 /usr/local/bin/

# macOS
brew install iqtree2
```

### PLINK2

```bash
# Linux
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240123.zip
unzip plink2_linux_avx2_20240123.zip
sudo mv plink2 /usr/local/bin/

# macOS
brew install plink2
```

### BCFtools and Samtools

```bash
# Ubuntu/Debian
sudo apt-get install bcftools samtools

# macOS
brew install bcftools samtools

# Conda
conda install -c bioconda bcftools samtools
```

## Verification

After installation, verify that everything is working:

```bash
# Check SNPhylo2 version
snphylo2 --version

# Validate installation
snphylo2 validate

# Should show:
# Core Tools:
#   ✓ python: 3.11.x
#   ✓ tabix: found
#   ✓ bgzip: found
# Tree Engines:
#   ✓ iqtree2: 2.3.x
#   ✓ raxml-ng: 1.2.x (optional)
#   ✓ FastTree: 2.1.x (optional)
# QC/Filtering:
#   ✓ plink2: 2.0.x
#   ✓ bcftools: 1.19
```

## Troubleshooting

### Issue: "Command not found" after Conda installation

```bash
# Ensure conda is properly initialized
conda init bash  # or zsh, fish, etc.
# Restart your shell
```

### Issue: External tools not found

```bash
# Check if tools are in PATH
which iqtree2
which plink2

# If not found, add to PATH
export PATH=$PATH:/path/to/tools
```

### Issue: Permission denied with Docker

```bash
# Add user to docker group (requires logout/login)
sudo usermod -aG docker $USER
```

### Issue: Missing Python packages

```bash
# Reinstall with all dependencies
pip install --force-reinstall -e ".[dev]"
```

## Uninstallation

```bash
# Conda
conda remove snphylo2

# pip
pip uninstall snphylo2

# Docker
docker rmi ghcr.io/snphylo2/snphylo2:latest
```

## Next Steps

- [Quick Start Guide](quickstart.md)
- [Configuration Guide](configuration.md)
- [Tutorials](tutorials/basic-usage.md)
