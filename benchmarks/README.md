# SNPhylo2 Benchmarks

This directory contains benchmarking suites for evaluating SNPhylo2 performance and accuracy.

## Tree Accuracy Benchmark

Evaluates phylogenetic inference accuracy on simulated data with known true trees.

### Reference

Manuscript section: "Benchmarking: Tree Accuracy on Simulated Data"

To evaluate topological accuracy, we generated 50 replicate coalescent simulations using msprime and stdpopsim with known demographic histories for three organisms:
- Human-like: Out-of-Africa three-population model
- Arabidopsis-like: Selfing model with rapid LD decay
- Rice-like: Structured population model with indica/japonica/aus subpopulations

### Running the Benchmark

```bash
# Install dependencies
pip install msprime stdpopsim tskit ete3 scipy pandas numpy psutil

# Run full benchmark (50 replicates per model)
python tree_accuracy_benchmark.py -o results/ -n 50 --models all

# Run single model
python tree_accuracy_benchmark.py -o human_results/ -n 50 --models human

# Quick test (5 replicates)
python tree_accuracy_benchmark.py -o test_results/ -n 5 --models human
```

### Expected Results

Based on manuscript benchmarks:

| Method | Mean Normalized RF | Standard Deviation |
|--------|-------------------|-------------------|
| SNPhylo2 (GTR+ASC+G4, 1000 UFBoot) | 0.089 | 0.031 |
| Manual PLINK2 + IQ-TREE2 | 0.091 | 0.033 |
| FastTree (GTR) | 0.118 | 0.038 |
| PLINK NJ | 0.142 | 0.044 |

Statistical significance: SNPhylo2 vs PLINK NJ (p < 0.001, paired Wilcoxon test)

### Output Files

- `raw_results.csv` - All replicate results
- `summary_statistics.csv` - Mean, std, min, max per method
- `statistical_tests.txt` - Pairwise Wilcoxon tests
- `replicate_XXX/` - Individual replicate data and trees

### Methodology

1. **Data Simulation**: msprime/stdpopsim generates VCF and true ARG tree
2. **Tree Inference**: Each method (SNPhylo2, PLINK+NJ, FastTree, Manual) builds tree
3. **Comparison**: Robinson-Foulds distance between inferred and true trees
4. **Statistics**: Mean normalized RF across 50 replicates, pairwise Wilcoxon tests

### True Tree Extraction

True species trees extracted from the ancestral recombination graph (ARG) using tskit.

### Notes

- SNPhylo2 automation achieves equivalent accuracy to expert-tuned manual workflows
- Ascertainment bias correction (ASC) essential for SNP-only data
- Ultrafast bootstrap (UFBoot) provides robust support values without runtime penalty
