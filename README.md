# SNPhylo2

[![Version](https://img.shields.io/badge/version-0.1.0-blue)](https://github.com/ank-man/snphylo2)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Tests](https://github.com/ank-man/snphylo2/actions/workflows/tests.yml/badge.svg)](https://github.com/ank-man/snphylo2/actions)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://snphylo2.readthedocs.io)

<p align="center">
  <img src="fig1.jpeg" alt="SNPhylo2 Workflow" width="800"/>
  <br>
  <b>Figure 1:</b> SNPhylo2 analysis workflow from raw VCF to phylogenetic tree with population genomics analyses
</p>

> **Next-Generation Phylogenomic Pipeline from SNP Data**

SNPhylo2 is a modern, scalable, and reproducible pipeline for constructing phylogenetic trees from large SNP datasets. It represents a complete reimagining of the original [SNPhylo](https://github.com/thlee/SNPhylo) tool for the era of population-scale resequencing, pangenomics, and cloud computing.

---

## 🚀 Features

- **🔄 Scalable**: Process 10,000+ samples and 100M+ SNPs with chunked streaming
- **📊 Multiple Formats**: VCF, BCF, PLINK, HapMap, GDS, FASTA support
- **🌳 Modern Tree Methods**: IQ-TREE2, RAxML-NG, FastTree, PhyML integration
- **🧬 Population Genomics**: Built-in PCA, FST, IBS, kinship, LD decay (PopLDdecay-style)
- **☁️ Cloud/HPC Ready**: Nextflow, Snakemake, Docker, Singularity workflows
- **📈 Publication-Ready**: Automated HTML reports with high-quality figures
- **🐍 Python-First**: Clean API with optional Rust backend for performance

---

## 📦 Installation

### Option 1: Conda (Recommended)

```bash
conda install -c bioconda -c conda-forge snphylo2
```

### Option 2: From Source

```bash
git clone https://github.com/ank-man/snphylo2.git
cd snphylo2
pip install -e .
```

---

## 🏃 Quick Start

### One-Command Analysis

```bash
snphylo2 run -v input.vcf.gz --threads 16 -o results/
```

This executes the complete pipeline:
1. Quality Control
2. Variant Filtering (MAF, missingness, depth)
3. LD Pruning (PLINK2 backend)
4. Phylogenetic Tree Building (IQ-TREE2)
5. Report Generation

### Step-by-Step Analysis

```bash
# 1. QC and filter variants
snphylo2 filter -v input.vcf.gz -o filtered.vcf.gz --maf 0.05 --max-missing 0.2

# 2. LD pruning
snphylo2 prune -i filtered.vcf.gz -o pruned.vcf.gz --window 50 --r2 0.2

# 3. Build tree
snphylo2 tree -i pruned.vcf.gz -o tree.nwk --engine iqtree2 --bootstrap 1000

# 4. Population genomics
snphylo2 ld-decay -v input.vcf.gz -m metadata.tsv --plot
```

---

## 📊 What SNPhylo2 Does

SNPhylo2 automates the entire phylogenomic analysis pipeline:

| Step | Input | Output | Description |
|------|-------|--------|-------------|
| **QC** | Raw VCF | QC Report | Sample/variant quality metrics |
| **Filter** | VCF | Filtered VCF | MAF, missingness, depth filtering |
| **LD Prune** | Filtered VCF | Pruned VCF | Remove linked SNPs (r² threshold) |
| **Tree** | Pruned SNPs | Newick tree | ML phylogeny with bootstrap |
| **PopGen** | VCF + Metadata | PCA, FST, LD decay | Population structure analysis |

---

## 📈 Comparison with Original SNPhylo

| Feature | SNPhylo (2018) | SNPhylo2 |
|---------|----------------|----------|
| Max Samples | ~100 | ~100,000 |
| Max SNPs | ~1M | ~100M |
| Tree Engines | DNAML only | IQ-TREE2, RAxML-NG, FastTree |
| LD Pruning | SNPRelate | PLINK2 (100x faster) |
| Population Genomics | None | PCA, FST, IBS, LD decay |
| Parallelization | Single-threaded | Chunked, HPC, Cloud |
| Reporting | PNG tree | Full HTML dashboard |
| Containerization | None | Docker, Singularity |

---

## 🎯 Key Capabilities

### Phylogenetic Analysis
- **IQ-TREE2 integration**: Ultrafast bootstrap, model selection
- **Ascertainment bias correction**: Essential for SNP-only data
- **Multiple engines**: RAxML-NG, FastTree for different scales

### Population Genomics
- **PCA**: Principal component analysis on genotype matrix
- **FST**: Weir-Cockerham pairwise population differentiation
- **IBS**: Identity-by-state distance matrices
- **LD Decay**: PopLDdecay-style r² decay curves
- **Kinship**: VanRaden relationship matrix

### Data Handling
- **Streaming VCF/BCF**: Constant memory footprint
- **Chunked processing**: Parallel by chromosome
- **Checkpoint/resume**: Automatic pipeline recovery
- **Multiple formats**: VCF, BCF, PLINK, HapMap, GDS

---

## 🧬 Example: Rice 3K Analysis

```bash
# Download demo data
python examples/download_demo_data.py --dataset rice --samples 50

# Run complete pipeline with population structure
snphylo2 run \
  -v demo_data/rice_3k_50accessions_chr1.vcf.gz \
  --metadata demo_data/rice_3k_50accessions_chr1.metadata.tsv \
  --threads 8 \
  -o rice_results/

# Generate manuscript figures
python examples/generate_manuscript_figures.py \
  demo_data/rice_3k_50accessions_chr1.vcf.gz \
  --metadata demo_data/rice_3k_50accessions_chr1.metadata.tsv \
  -o rice_figures/ \
  --format pdf
```

**Outputs:**
- `rice_results/snphylo2_output.tree.nwk` - Phylogenetic tree
- `rice_results/snphylo2_output_report.html` - Interactive report
- `rice_figures/figure_ld_decay.pdf` - LD decay curves by subpopulation
- `rice_figures/figure_pca_pc1_pc2.pdf` - PCA with indica/japonica coloring

---

## ☁️ HPC/Cloud Deployment

### Nextflow on SLURM

```bash
nextflow run workflows/nextflow/main.nf \
  -profile slurm \
  --input huge_dataset.vcf.gz \
  --output results/
```

### Snakemake with Conda

```bash
snakemake --cores 16 --use-conda \
  -s workflows/snakemake/Snakefile \
  --config input=huge_dataset.vcf.gz
```

### Docker

```bash
docker run -v $(pwd):/data ghcr.io/ank-man/snphylo2:latest \
  run -v /data/input.vcf.gz -o /data/results
```

---

## 📚 Documentation

- [Installation Guide](docs/installation.md)
- [Tutorials](docs/tutorials/)
- [Configuration](docs/configuration.md)
- [API Reference](docs/api/)
- [Benchmarks](docs/benchmarks.md)

---

## 🤝 Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas needing help:
- Additional tree building engines (MrBayes, BEAST)
- Polyploid support
- Pangenome graph compatibility
- Cloud deployment guides

---

## 📄 Citation

If you use SNPhylo2 in your research, please cite:

```bibtex
@article{snphylo2_2026,
  title={SNPhylo2: A Scalable, Reproducible, and Population-Aware Pipeline 
         for Phylogenomic Inference from SNP Data},
  author={Sharma, Ankush and SNPhylo2 Development Team},
  journal={Bioinformatics},
  year={2026},
  publisher={Oxford University Press}
}
```

Also cite the original SNPhylo:

```bibtex
@article{lee2014snphylo,
  title={SNPhylo: a pipeline to construct a phylogenetic tree from huge SNP data},
  author={Lee, Tae-Ho and Guo, Hui and Wang, Xiyin and Kim, Changsoo and Paterson, Andrew H},
  journal={BMC genomics},
  volume={15},
  pages={1--10},
  year={2014},
  publisher={Springer}
}
```

---

## 📬 Contact

- **Issues**: [GitHub Issues](https://github.com/ank-man/snphylo2/issues)
- **Discussions**: [GitHub Discussions](https://github.com/ank-man/snphylo2/discussions)
- **Email**: Ankush Sharma (mr.ank2999@gmail.com)

---

**Made with ❤️ for the phylogenomics community**
