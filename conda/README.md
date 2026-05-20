# Bioconda Package for SNPhylo2

This directory contains the Conda recipe for building and distributing SNPhylo2 through the Bioconda channel.

## Quick Links

- **Bioconda**: https://bioconda.github.io
- **Package Page**: https://anaconda.org/bioconda/snphylo2
- **Installation**: `conda install -c bioconda -c conda-forge snphylo2`

## Recipe Structure

```
conda/
├── meta.yaml        # Conda package metadata and dependencies
├── build.sh         # Unix build script
├── bld.bat          # Windows build script
└── README.md        # This file
```

## Local Testing

### Test the Recipe Locally

```bash
# Install conda-build
conda install conda-build

# Build the package
cd conda
conda build . --output-folder ../conda-build

# Test the built package
conda install --use-local -c ../conda-build snphylo2
snphylo2 --version
```

### Test with Docker (recommended for Bioconda compatibility)

```bash
# Using Bioconda's testing environment
docker run -v $(pwd):/recipe bioconda/bioconda-utils-build \
    --packages snphylo2 --docker /recipe
```

## Submitting to Bioconda

### Method 1: Manual PR (Recommended for first release)

1. **Fork Bioconda Recipes**
   ```bash
   git clone https://github.com/bioconda/bioconda-recipes.git
   cd bioconda-recipes
   git remote add upstream https://github.com/bioconda/bioconda-recipes.git
   ```

2. **Create Recipe Directory**
   ```bash
   mkdir -p recipes/snphylo2
   cp /path/to/snphylo2/conda/meta.yaml recipes/snphylo2/
   cp /path/to/snphylo2/conda/build.sh recipes/snphylo2/
   cp /path/to/snphylo2/conda/bld.bat recipes/snphylo2/
   ```

3. **Update SHA256**
   ```bash
   # Get SHA256 of release tarball
   curl -sL https://github.com/ank-man/snphylo2/archive/v0.1.0.tar.gz | sha256sum
   
   # Update meta.yaml with the correct SHA256
   ```

4. **Test the Recipe**
   ```bash
   # Using Bioconda's test framework
   ./bootstrap.py
   bioconda-utils build recipes/snphylo2 config.yml
   ```

5. **Submit PR**
   ```bash
   git checkout -b snphylo2
   git add recipes/snphylo2/
   git commit -m "Add snphylo2 v0.1.0"
   git push origin snphylo2
   ```
   Then open PR at https://github.com/bioconda/bioconda-recipes/pulls

### Method 2: Using Bioconda Utils (for updates)

```bash
# Install bioconda-utils
conda install -c conda-forge -c bioconda bioconda-utils

# Update recipe for new version
# Edit recipes/snphylo2/meta.yaml with new version and SHA256

# Test
bioconda-utils build recipes/snphylo2 config.yml

# Lint
bioconda-utils lint recipes/snphylo2 config.yml
```

## Recipe Maintenance

### Version Updates

When releasing a new version:

1. Update `version` in `meta.yaml`
2. Update `sha256` with new release tarball hash
3. Reset `build: number: 0`
4. Test locally
5. Submit PR to Bioconda

### Adding Dependencies

Edit `requirements:` section in `meta.yaml`:

```yaml
requirements:
  host:
    - python >=3.10
    - pip
  run:
    - python >=3.10
    # Add new dependencies here
    - new-package >=1.0
```

### Platform Support

Current recipe supports:
- Linux (primary target for Bioconda)
- macOS (if dependencies available)
- Windows (limited support, some tools may not be available)

## Troubleshooting

### Common Issues

**Issue: SHA256 mismatch**
```bash
# Recalculate SHA256
curl -sL https://github.com/ank-man/snphylo2/archive/vX.X.X.tar.gz | sha256sum
```

**Issue: Missing dependencies**
- Check that all Python packages are available on conda-forge
- External tools (IQ-TREE2, PLINK2) must be on Bioconda

**Issue: Build fails on specific Python version**
- Check Python version constraints in meta.yaml
- Some dependencies may not support latest Python

### Bioconda Guidelines

Follow Bioconda's contribution guidelines:
- https://bioconda.github.io/contributor/guidelines.html
- https://bioconda.github.io/contributor/linting.html

## Contact

For issues with the Bioconda package:
- Open issue at https://github.com/bioconda/bioconda-recipes/issues
- Mention @ank-man for SNPhylo2-specific issues

For general SNPhylo2 issues:
- https://github.com/ank-man/snphylo2/issues
