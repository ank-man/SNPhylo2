# SNPhylo2

**Next-Generation Phylogenomic Pipeline from SNP Data**

[![Version](https://img.shields.io/badge/version-0.1.0-blue)](https://github.com/snphylo2/snphylo2)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://snphylo2.readthedocs.io)

SNPhylo2 is a modern, scalable, and reproducible pipeline for constructing phylogenetic trees from large SNP datasets. It represents a complete reimagining of the original [SNPhylo](https://chibba.agtec.uga.edu/snphylo) (Lee et al., 2014) tool for the era of population-scale resequencing, pangenomics, and cloud computing.

**Original SNPhylo**: https://chibba.agtec.uga.edu/snphylo

## 🚀 Features

- **🔄 Scalable**: Process 10,000+ samples and 100M+ SNPs with chunked streaming
- **📊 Multiple Formats**: VCF, BCF, PLINK, HapMap, GDS, FASTA support
- **🌳 Modern Tree Methods**: IQ-TREE2, RAxML-NG, FastTree integration
- **🧬 Population Genomics**: Built-in PCA, FST, IBS, ADMIXTURE support
- **☁️ Cloud/HPC Ready**: Nextflow, Snakemake, Docker, Singularity
- **📈 Publication-Ready**: Automated HTML reports with high-quality figures
- **🐍 Python-First**: Clean API with optional Rust backend for performance

## 🏃 Quick Start

```bash
# Install via Conda (recommended)
conda install -c bioconda snphylo2

# Run complete pipeline
snphylo2 run -v input.vcf.gz --threads 16 -o results/

# Or step by step
snphylo2 filter -v input.vcf.gz -o filtered.vcf.gz --maf 0.05
snphylo2 prune -i filtered.vcf.gz -o pruned.vcf.gz --r2 0.2
snphylo2 tree -i pruned.vcf.gz -o tree.nwk --bootstrap 1000
```

## 📊 Comparison with Original SNPhylo

| Feature | SNPhylo | SNPhylo2 |
|---------|---------|----------|
| Max Samples | ~100 | ~100,000 |
| Max SNPs | ~1M | ~100M |
| Tree Engines | DNAML | IQ-TREE2, RAxML-NG, FastTree |
| Parallelization | None | Chunked, HPC, Cloud |
| Missing Data | Basic | Sensitivity Analysis |
| Polyploidy | No | Yes |
| Reporting | PNG | HTML Dashboard |

## 🎯 Use Cases

- **Crop Diversity Panels**: Analyze thousands of crop accessions
- **Population Genomics**: Integrate phylogeny with population structure
- **Microbial Outbreaks**: Rapid phylogenetic analysis of pathogen genomes
- **Pangenome-Aware**: Handle multiple reference haplotypes

## 📖 Documentation

- [Installation Guide](installation.md)
- [Quick Start](quickstart.md)
- [Tutorials](tutorials/basic-usage.md)
- [API Reference](api/index.md)
- [Benchmarks](benchmarks.md)

## 🤝 Contributing

We welcome contributions! See our [Contributing Guide](contributing.md) for details.

## 📄 Citation

If you use SNPhylo2 in your research, please cite:

```bibtex
@article{snphylo2_2026,
  title={SNPhylo2: A Scalable, Reproducible, and Population-Aware Pipeline 
         for Phylogenomic Inference from SNP Data},
  author={Ankush Sharma and SNPhylo2 Development Team},
  journal={Bioinformatics},
  year={2026}
}
```

## 📬 Support

- GitHub Issues: [github.com/snphylo2/snphylo2/issues](https://github.com/snphylo2/snphylo2/issues)
- GitHub Discussions: [github.com/snphylo2/snphylo2/discussions](https://github.com/snphylo2/snphylo2/discussions)
- Email: Ankush Sharma (mr.ank2999@gmail.com)

---

**Made with ❤️ for the phylogenomics community**
