# Publishing SNPhylo2

This document describes how to publish SNPhylo2 to PyPI and Bioconda.

## PyPI (Python Package Index)

### Prerequisites

1. **PyPI Account**: Register at https://pypi.org/account/register/
2. **TestPyPI Account**: Register at https://test.pypi.org/account/register/
3. **API Token**: Generate at https://pypi.org/manage/account/token/

### Automated Publishing (GitHub Actions)

The easiest method - just create a GitHub release:

1. **Create a new release** on GitHub:
   - Go to https://github.com/ank-man/snphylo2/releases
   - Click "Draft a new release"
   - Choose tag: `v0.1.0` (must match version in `pyproject.toml`)
   - Title: "SNPhylo2 v0.1.0"
   - Add release notes
   - Click "Publish release"

2. **GitHub Actions will automatically**:
   - Build the package
   - Publish to TestPyPI first
   - Then publish to PyPI

### Manual Publishing (if needed)

```bash
# 1. Install build tools
pip install build twine

# 2. Build package
python -m build

# 3. Check package
python -m twine check dist/*

# 4. Upload to TestPyPI first (test)
python -m twine upload --repository testpypi dist/*

# 5. Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ snphylo2

# 6. Upload to PyPI (production)
python -m twine upload dist/*
```

### Installation from PyPI

```bash
pip install snphylo2
```

With optional dependencies:
```bash
pip install snphylo2[plotting,population,dev]
```

## Bioconda

### Prerequisites

1. **GitHub Account** (already have)
2. **Fork Bioconda Recipes**: https://github.com/bioconda/bioconda-recipes

### Submitting Recipe

```bash
# 1. Fork and clone bioconda-recipes
git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git

# 2. Create recipe
cp -r /path/to/snphylo2/conda recipes/snphylo2

# 3. Update SHA256 in recipes/snphylo2/meta.yaml
curl -sL https://github.com/ank-man/snphylo2/archive/v0.1.0.tar.gz | sha256sum

# 4. Test locally (optional)
conda build recipes/snphylo2

# 5. Submit PR
git checkout -b snphylo2
git add recipes/snphylo2
git commit -m "Add snphylo2 v0.1.0"
git push origin snphylo2
# Then open PR at https://github.com/bioconda/bioconda-recipes/pulls
```

### Installation from Bioconda

```bash
conda install -c bioconda -c conda-forge snphylo2
```

## Version Numbering

Follow semantic versioning: `MAJOR.MINOR.PATCH`

- **MAJOR**: Breaking changes (API changes)
- **MINOR**: New features (backwards compatible)
- **PATCH**: Bug fixes

Examples:
- `0.1.0` - Initial release
- `0.1.1` - Bug fix
- `0.2.0` - New features
- `1.0.0` - Stable release

## Release Checklist

Before publishing:

- [ ] Update version in `pyproject.toml`
- [ ] Update `conda/meta.yaml` version and SHA256
- [ ] Update `CHANGELOG.md`
- [ ] Run tests: `pytest`
- [ ] Test build locally: `python -m build`
- [ ] Create git tag: `git tag v0.1.0`
- [ ] Push tag: `git push origin v0.1.0`
- [ ] Create GitHub release
- [ ] Verify PyPI publication
- [ ] Submit Bioconda PR
- [ ] Update documentation

## Troubleshooting

### PyPI Issues

**"File already exists"**
- PyPI doesn't allow overwriting files
- Must increment version number

**"Invalid API Token"**
- Check token hasn't expired
- Verify token has "Upload" scope

### Bioconda Issues

**"Recipe fails linting"**
- Run `bioconda-utils lint recipes/snphylo2`
- Fix all issues before PR

**"Dependency not available"**
- All dependencies must be on conda-forge or bioconda
- Check https://anaconda.org/search

## Contact

For publishing issues:
- PyPI: https://github.com/pypi/support
- Bioconda: https://github.com/bioconda/bioconda-recipes/issues
- SNPhylo2: https://github.com/ank-man/snphylo2/issues
