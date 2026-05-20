# SNPhylo2 Examples and Demo Data

This directory contains scripts and datasets for testing and demonstrating SNPhylo2.

## Quick Start

### 1. Download Demo Datasets

```bash
# Download all demo datasets
python download_demo_data.py --dataset all -o demo_data/

# Or download specific datasets
python download_demo_data.py --dataset human --samples 100
python download_demo_data.py --dataset rice --samples 50
python download_demo_data.py --dataset arabidopsis --samples 80

# Create tiny test dataset for quick testing
python download_demo_data.py --dataset tiny
```

Available datasets:
- **Human (1000 Genomes)**: 100 samples, chr21:10-11Mb (simulated if download fails)
- **Rice (3K)**: 50 accessions with indica/japonica/aus population structure
- **Arabidopsis (1001G)**: 80 accessions with rapid LD decay
- **Tiny test**: 10 samples, 100 SNPs for quick testing

### 2. Generate Manuscript Figures

```bash
# Generate all figures from a dataset
python generate_manuscript_figures.py demo_data/tiny_test.vcf.gz -o figures/

# With metadata for population coloring
python generate_manuscript_figures.py \
    demo_data/rice_3k_50accessions_chr1.vcf.gz \
    --metadata demo_data/rice_3k_50accessions_chr1.metadata.tsv \
    -o rice_figures/ \
    --format pdf
```

Generated figures:
- `figure_workflow.pdf` - SNPhylo2 workflow diagram
- `figure_ld_decay.pdf` - LD decay curves (PopLDdecay-style)
- `figure_pca_pc1_pc2.pdf` - PCA scatter plot
- `figure_pca_pc1_pc3.pdf` - Additional PCA plot
- `figure_ld_heatmap.pdf` - LD heatmap for genomic region
- `figure_fst_heatmap.pdf` - Population differentiation (if metadata provided)
- `figure_combined_overview.pdf` - Multi-panel summary figure

Data tables:
- `ld_decay_{population}.tsv` - LD decay data
- `pca_eigenvalues.tsv` - PCA eigenvalues and variance
- `fst_values.tsv` - Pairwise FST matrix

### 3. Run SNPhylo2 on Demo Data

```bash
# Quick test
snphylo2 run -v demo_data/tiny_test.vcf.gz -o tiny_results/

# With configuration
snphylo2 init --preset plant -o rice_config.yaml
snphylo2 run -v demo_data/rice_3k_50accessions_chr1.vcf.gz \
    --config rice_config.yaml \
    -o rice_results/

# Step-by-step analysis
snphylo2 qc -v demo_data/human_simulated_100samples_chr21_10000000_11000000.vcf.gz
snphylo2 filter -v demo_data/human_simulated_100samples_chr21_10000000_11000000.vcf.gz \
    -o human_filtered.vcf.gz --maf 0.05
snphylo2 prune -i human_filtered.vcf.gz -o human_pruned.vcf.gz --r2 0.2
snphylo2 tree -i human_pruned.vcf.gz -o human_tree.nwk --bootstrap 1000
```

### 4. LD Decay Analysis (PopLDdecay-style)

```python
from snphylo2.population.ld_decay import LDDecayAnalyzer
from snphylo2.visualization.ld_plots import LDPlotter

# Analyze LD decay
analyzer = LDDecayAnalyzer(
    max_distance=200_000,  # 200 kb
    bin_size=5_000,      # 5 kb bins
)

# With population groups
sample_groups = {
    'indica': ['Sample_001', 'Sample_004', ...],
    'japonica': ['Sample_002', 'Sample_005', ...],
    'aus': ['Sample_003', 'Sample_006', ...],
}

results = analyzer.analyze(
    'demo_data/rice_3k_50accessions_chr1.vcf.gz',
    sample_groups=sample_groups,
)

# Create publication-quality plot
plotter = LDPlotter(width=1200, height=900, dpi=300)
plotter.plot_ld_decay(
    results,
    output_path='manuscript_ld_decay.pdf',
    title='Linkage Disequilibrium Decay in Rice Populations',
    show_confidence=True,
    manuscript_ready=True,
)

# Compare populations
comparison = analyzer.compare_populations(results)
print(comparison)
```

## Dataset Characteristics

### Human (1000 Genomes-style)
- **LD decay**: ~100-200 kb (r² = 0.5)
- **Population structure**: CEU, YRI, CHB populations
- **Heterozygosity**: ~0.1%
- **Use case**: Human population genetics, GWAS

### Rice (3K Genomes-style)
- **LD decay**: ~100-200 kb (selfing species)
- **Population structure**: indica, japonica, aus subspecies
- **Heterozygosity**: Moderate
- **Use case**: Crop breeding, diversity analysis

### Arabidopsis (1001 Genomes-style)
- **LD decay**: ~10 kb (rapid decay in outcrossing)
- **Population structure**: European accessions
- **Heterozygosity**: High homozygosity (selfing)
- **Use case**: Adaptation studies, rapid LD mapping

## File Formats

### Metadata TSV Format
```
sample_id	population	subspecies	country
Sample_001	indica	Oryza_sativa	India
Sample_002	japonica	Oryza_sativa	Japan
Sample_003	aus	Oryza_sativa	Bangladesh
```

### Output Data Tables

**LD Decay TSV:**
```
distance	mean_r2	p95_r2	p5_r2	n_pairs
5000	0.823	0.951	0.612	5234
10000	0.654	0.843	0.421	8932
...
```

**PCA TSV:**
```
sample_id	PC1	PC2	PC3	...
Sample_001	-12.3	5.4	-2.1	...
Sample_002	15.2	-8.1	3.5	...
```

## Nextflow Workflow Example

```bash
# Run on SLURM cluster
nextflow run ../workflows/nextflow/main.nf \
    -profile slurm \
    --input demo_data/rice_3k_50accessions_chr1.vcf.gz \
    --output rice_nextflow_results/
```

## Snakemake Workflow Example

```bash
# Edit config.yaml to point to demo data, then run:
snakemake --cores 4 --use-conda \
    -s ../workflows/snakemake/Snakefile \
    --configfile ../workflows/snakemake/config.yaml
```

## Citation

If you use these demo datasets or examples in your research, please cite:

```
SNPhylo2 Development Team (2026). SNPhylo2: A Scalable, Reproducible, 
and Population-Aware Pipeline for Phylogenomic Inference from SNP Data.
Bioinformatics, xx(xx), xxx-xxx.
```

## License

The demo datasets are generated for testing purposes and are released under the same MIT license as SNPhylo2.
