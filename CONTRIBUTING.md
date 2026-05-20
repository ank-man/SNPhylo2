# Contributing to SNPhylo2

Thank you for your interest in contributing to SNPhylo2! This document provides guidelines for contributing to the project.

## 🚀 Quick Start

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/snphylo2.git`
3. Create a branch: `git checkout -b feature/my-feature`
4. Make changes and commit
5. Push to your fork: `git push origin feature/my-feature`
6. Open a Pull Request

## 📋 Development Setup

```bash
# Clone repository
git clone https://github.com/snphylo2/snphylo2.git
cd snphylo2

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

## 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=snphylo2 --cov-report=html

# Run specific test
pytest tests/unit/test_vcf_reader.py

# Run slow tests
pytest --runslow
```

## 📝 Code Style

- Follow PEP 8 with 100 character line length
- Use type hints where possible
- Write docstrings for all public functions/classes
- Run `black` and `ruff` before committing:

```bash
black snphylo2 tests
ruff check snphylo2 tests
```

## 🐛 Reporting Bugs

Please use GitHub Issues and include:

- Python version
- SNPhylo2 version
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Error messages

## 💡 Suggesting Features

Feature requests are welcome! Please:

- Check if the feature already exists
- Describe the use case
- Explain why it would be valuable
- Provide examples if possible

## 📚 Documentation

- Update docstrings for API changes
- Update README.md for user-facing changes
- Add examples for new features
- Update CHANGELOG.md

## 🏗️ Pull Request Process

1. Update tests for new functionality
2. Ensure all tests pass
3. Update documentation
4. Add entry to CHANGELOG.md
5. Request review from maintainers

## 📋 Commit Message Format

```
type(scope): subject

body (optional)

footer (optional)
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `style`: Formatting
- `refactor`: Code restructuring
- `test`: Tests
- `chore`: Maintenance

Example:
```
feat(filtering): add HWE filter

Implements Hardy-Weinberg equilibrium filtering for population
genetics analyses. Includes p-value threshold configuration.

Closes #123
```

## 🔧 Areas Needing Help

- [ ] Additional tree building engines (PhyML, MrBayes)
- [ ] Visualization improvements (interactive trees)
- [ ] Cloud deployment guides (AWS, GCP, Azure)
- [ ] Additional population genetics analyses
- [ ] Performance optimizations
- [ ] Tutorial notebooks
- [ ] Translation of documentation

## 💬 Community

- GitHub Discussions: Questions and ideas
- GitHub Issues: Bug reports and features
- Email: Ankush Sharma (mr.ank2999@gmail.com)

## 🏆 Recognition

Contributors will be:
- Listed in CONTRIBUTORS.md
- Mentioned in release notes
- Thanked in presentations and publications

Thank you for contributing to SNPhylo2!
