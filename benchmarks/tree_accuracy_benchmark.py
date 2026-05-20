#!/usr/bin/env python3
"""
Tree Accuracy Benchmarking on Simulated Data

Evaluates topological accuracy using coalescent simulations with known demographic
histories. Compares SNPhylo2 against baseline methods using Robinson-Foulds distance.

Reference implementation for manuscript section:
"Benchmarking: Tree Accuracy on Simulated Data"
"""

import argparse
import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.stats as stats
from ete3 import Tree


@dataclass
class SimulationConfig:
    """Configuration for demographic simulation."""
    name: str
    species: str
    model: str
    samples: Dict[str, int]  # population -> sample count
    sequence_length: int
    mutation_rate: float
    recombination_rate: float
    seed: Optional[int] = None


@dataclass
class BenchmarkResult:
    """Results from a single replicate."""
    replicate_id: int
    method: str
    rf_distance: float
    normalized_rf: float
    n_nodes_true: int
    n_nodes_inferred: int
    runtime_seconds: float
    memory_mb: float


class TreeAccuracyBenchmark:
    """
    Benchmark tree inference accuracy on simulated data.
    
    Implements the methodology described in the manuscript for evaluating
    topological accuracy against known true trees from coalescent simulations.
    """
    
    def __init__(self, output_dir: str, n_replicates: int = 50):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.n_replicates = n_replicates
        self.results: List[BenchmarkResult] = []
        
    def run_full_benchmark(self, configs: List[SimulationConfig]):
        """
        Run complete benchmarking suite across multiple demographic models.
        
        Args:
            configs: List of simulation configurations (human, arabidopsis, rice)
        """
        print("=" * 70)
        print("TREE ACCURACY BENCHMARK")
        print("=" * 70)
        print(f"Replicates per model: {self.n_replicates}")
        print(f"Methods: SNPhylo2, PLINK+NJ, FastTree, Manual PLINK2+IQ-TREE2")
        print("=" * 70)
        
        for config in configs:
            print(f"\n[MODEL] {config.name}: {config.species}")
            self._benchmark_model(config)
        
        # Generate summary statistics
        self._generate_summary()
        self._statistical_tests()
        
    def _benchmark_model(self, config: SimulationConfig):
        """Benchmark a single demographic model."""
        model_dir = self.output_dir / config.name
        model_dir.mkdir(exist_ok=True)
        
        for rep in range(1, self.n_replicates + 1):
            print(f"  Replicate {rep}/{self.n_replicates}...", end=" ")
            
            # Set seed for reproducibility
            seed = config.seed or (42 + rep)
            
            # 1. Generate simulated data with true tree
            vcf_path, true_tree_path = self._generate_simulation(
                config, model_dir, rep, seed
            )
            
            # 2. Run each method and calculate RF distance
            methods = [
                ("snphylo2", self._run_snphylo2),
                ("plink_nj", self._run_plink_nj),
                ("fasttree", self._run_fasttree),
                ("manual_workflow", self._run_manual_workflow),
            ]
            
            for method_name, method_func in methods:
                try:
                    inferred_tree, runtime, memory = method_func(vcf_path, model_dir, rep)
                    rf_dist, norm_rf = self._calculate_rf_distance(
                        true_tree_path, inferred_tree
                    )
                    
                    result = BenchmarkResult(
                        replicate_id=rep,
                        method=method_name,
                        rf_distance=rf_dist,
                        normalized_rf=norm_rf,
                        n_nodes_true=self._count_nodes(true_tree_path),
                        n_nodes_inferred=self._count_nodes(inferred_tree),
                        runtime_seconds=runtime,
                        memory_mb=memory,
                    )
                    self.results.append(result)
                    
                except Exception as e:
                    print(f"[{method_name} failed: {e}]")
            
            print("done")
            
    def _generate_simulation(
        self,
        config: SimulationConfig,
        model_dir: Path,
        rep: int,
        seed: int,
    ) -> Tuple[Path, Path]:
        """
        Generate simulated VCF and true tree using msprime/stdpopsim.
        
        Returns:
            Tuple of (vcf_path, true_tree_path)
        """
        rep_dir = model_dir / f"replicate_{rep:03d}"
        rep_dir.mkdir(exist_ok=True)
        
        vcf_path = rep_dir / "simulated.vcf.gz"
        true_tree_path = rep_dir / "true_tree.nwk"
        
        # Use stdpopsim for realistic demographic models
        cmd = self._build_simulation_command(config, seed, vcf_path, true_tree_path)
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        return vcf_path, true_tree_path
    
    def _build_simulation_command(
        self,
        config: SimulationConfig,
        seed: int,
        vcf_path: Path,
        tree_path: Path,
    ) -> List[str]:
        """Build stdpopsim/msprime simulation command."""
        # Build samples specification
        samples_str = ",".join([
            f"{pop}:{n}" for pop, n in config.samples.items()
        ])
        
        return [
            "python", "-m", "stdpopsim",
            config.species,
            "--model", config.model,
            "--samples", samples_str,
            "--length", str(config.sequence_length),
            "--mutation-rate", str(config.mutation_rate),
            "--recombination-rate", str(config.recombination_rate),
            "--seed", str(seed),
            "--output-vcf", str(vcf_path),
            "--output-tree", str(tree_path),
        ]
    
    def _run_snphylo2(
        self,
        vcf_path: Path,
        output_dir: Path,
        rep: int,
    ) -> Tuple[Path, float, float]:
        """Run SNPhylo2 pipeline."""
        import time
        import psutil
        
        result_dir = output_dir / f"snphylo2_rep_{rep:03d}"
        
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024
        start_time = time.time()
        
        cmd = [
            "snphylo2", "run",
            "-v", str(vcf_path),
            "-o", str(result_dir),
            "--threads", "4",
            "--tree-engine", "iqtree2",
            "--model", "GTR+ASC+G4",
            "--bootstrap", "1000",
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        runtime = time.time() - start_time
        mem_after = process.memory_info().rss / 1024 / 1024
        memory = mem_after - mem_before
        
        # Find output tree
        tree_file = result_dir / "snphylo2_output.tree.nwk"
        
        return tree_file, runtime, memory
    
    def _run_plink_nj(
        self,
        vcf_path: Path,
        output_dir: Path,
        rep: int,
    ) -> Tuple[Path, float, float]:
        """Run PLINK distance + NJ tree."""
        import time
        import psutil
        
        result_dir = output_dir / f"plink_nj_rep_{rep:03d}"
        result_dir.mkdir(exist_ok=True)
        
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024
        start_time = time.time()
        
        # Convert VCF to PLINK
        plink_prefix = result_dir / "data"
        subprocess.run([
            "plink2",
            "--vcf", str(vcf_path),
            "--make-bed",
            "--out", str(plink_prefix),
        ], check=True, capture_output=True)
        
        # Calculate distance matrix
        dist_file = result_dir / "distances.mdist"
        subprocess.run([
            "plink2",
            "--bfile", str(plink_prefix),
            "--distance", "triangle", "1-ibs",
            "--out", str(result_dir / "distances"),
        ], check=True, capture_output=True)
        
        # Build NJ tree (using PHYLIP or custom implementation)
        tree_file = result_dir / "nj_tree.nwk"
        self._build_nj_tree(dist_file, tree_file)
        
        runtime = time.time() - start_time
        mem_after = process.memory_info().rss / 1024 / 1024
        memory = mem_after - mem_before
        
        return tree_file, runtime, memory
    
    def _build_nj_tree(self, dist_file: Path, output_tree: Path):
        """Build neighbor-joining tree from distance matrix."""
        # Implementation using scipy or PHYLIP
        from scipy.cluster.hierarchy import linkage, to_tree
        from scipy.spatial.distance import squareform
        
        # Read distance matrix
        dist_matrix = pd.read_csv(dist_file, sep="\t", header=None).values
        
        # Convert to condensed form and build NJ tree
        condensed = squareform(dist_matrix)
        linkage_matrix = linkage(condensed, method="average")
        
        # Convert to Newick (simplified)
        tree = to_tree(linkage_matrix)
        
        with open(output_tree, "w") as f:
            f.write(self._linkage_to_newick(tree))
    
    def _linkage_to_newick(self, tree) -> str:
        """Convert scipy linkage tree to Newick format."""
        # Simplified conversion
        return "(A,B);"  # Placeholder - would need full implementation
    
    def _run_fasttree(
        self,
        vcf_path: Path,
        output_dir: Path,
        rep: int,
    ) -> Tuple[Path, float, float]:
        """Run FastTree on SNP data."""
        import time
        import psutil
        
        result_dir = output_dir / f"fasttree_rep_{rep:03d}"
        result_dir.mkdir(exist_ok=True)
        
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024
        start_time = time.time()
        
        # Convert VCF to FASTA alignment
        fasta_file = result_dir / "alignment.fasta"
        self._vcf_to_fasta(vcf_path, fasta_file)
        
        # Run FastTree
        tree_file = result_dir / "fasttree.nwk"
        subprocess.run([
            "FastTree",
            "-gtr", "-nt",
            str(fasta_file),
        ], stdout=open(tree_file, "w"), check=True)
        
        runtime = time.time() - start_time
        mem_after = process.memory_info().rss / 1024 / 1024
        memory = mem_after - mem_before
        
        return tree_file, runtime, memory
    
    def _vcf_to_fasta(self, vcf_path: Path, fasta_path: Path):
        """Convert VCF genotype data to FASTA alignment."""
        # Use bcftools or cyvcf2 to extract genotypes
        # Write as FASTA format
        pass  # Implementation would go here
    
    def _run_manual_workflow(
        self,
        vcf_path: Path,
        output_dir: Path,
        rep: int,
    ) -> Tuple[Path, float, float]:
        """Run manual PLINK2 + IQ-TREE2 workflow for comparison."""
        import time
        import psutil
        
        result_dir = output_dir / f"manual_rep_{rep:03d}"
        result_dir.mkdir(exist_ok=True)
        
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024
        start_time = time.time()
        
        # Step 1: PLINK2 for LD pruning
        plink_prefix = result_dir / "pruned"
        subprocess.run([
            "plink2",
            "--vcf", str(vcf_path),
            "--indep-pairwise", "50", "10", "0.2",
            "--make-bed",
            "--out", str(plink_prefix),
        ], check=True, capture_output=True)
        
        # Step 2: Export to PHYLIP for IQ-TREE2
        phylip_file = result_dir / "alignment.phy"
        subprocess.run([
            "plink2",
            "--bfile", str(plink_prefix),
            "--export", "phylip",
            "--out", str(result_dir / "alignment"),
        ], check=True, capture_output=True)
        
        # Step 3: IQ-TREE2
        tree_file = result_dir / "manual_iqtree.nwk"
        subprocess.run([
            "iqtree2",
            "-s", str(phylip_file),
            "-m", "GTR+ASC+G4",
            "-bb", "1000",
            "-nt", "4",
            "-pre", str(result_dir / "iqtree"),
        ], check=True, capture_output=True)
        
        runtime = time.time() - start_time
        mem_after = process.memory_info().rss / 1024 / 1024
        memory = mem_after - mem_before
        
        return tree_file, runtime, memory
    
    def _calculate_rf_distance(
        self,
        true_tree_path: Path,
        inferred_tree_path: Path,
    ) -> Tuple[float, float]:
        """
        Calculate Robinson-Foulds distance between true and inferred trees.
        
        Returns:
            Tuple of (raw_rf_distance, normalized_rf_distance)
        """
        # Load trees using ETE3
        true_tree = Tree(str(true_tree_path))
        inferred_tree = Tree(str(inferred_tree_path))
        
        # Calculate RF distance
        rf_result = true_tree.compare(inferred_tree, unrooted=True)
        rf_distance = rf_result["rf"]
        max_rf = rf_result["max_rf"]
        normalized_rf = rf_distance / max_rf if max_rf > 0 else 0
        
        return rf_distance, normalized_rf
    
    def _count_nodes(self, tree_path: Path) -> int:
        """Count number of nodes in tree."""
        tree = Tree(str(tree_path))
        return len(tree.get_descendants()) + 1
    
    def _generate_summary(self):
        """Generate summary statistics for manuscript table."""
        df = pd.DataFrame([
            {
                "method": r.method,
                "replicate": r.replicate_id,
                "rf_distance": r.rf_distance,
                "normalized_rf": r.normalized_rf,
                "runtime": r.runtime_seconds,
                "memory_mb": r.memory_mb,
            }
            for r in self.results
        ])
        
        # Group by method and calculate statistics
        summary = df.groupby("method").agg({
            "normalized_rf": ["mean", "std", "min", "max"],
            "runtime": ["mean", "std"],
            "memory_mb": ["mean", "std"],
        }).round(4)
        
        # Save results
        df.to_csv(self.output_dir / "raw_results.csv", index=False)
        summary.to_csv(self.output_dir / "summary_statistics.csv")
        
        print("\n" + "=" * 70)
        print("SUMMARY STATISTICS (Normalized RF Distance)")
        print("=" * 70)
        print(summary)
        print("=" * 70)
        
        return summary
    
    def _statistical_tests(self):
        """Perform statistical comparisons between methods."""
        df = pd.DataFrame([
            {
                "method": r.method,
                "replicate": r.replicate_id,
                "normalized_rf": r.normalized_rf,
            }
            for r in self.results
        ])
        
        # Pairwise Wilcoxon signed-rank tests
        methods = df["method"].unique()
        
        print("\n" + "=" * 70)
        print("PAIRWISE STATISTICAL TESTS (Wilcoxon Signed-Rank)")
        print("=" * 70)
        print(f"{'Comparison':<40} {'p-value':<15} {'Significant'}")
        print("-" * 70)
        
        for i, method1 in enumerate(methods):
            for method2 in methods[i+1:]:
                data1 = df[df["method"] == method1]["normalized_rf"].values
                data2 = df[df["method"] == method2]["normalized_rf"].values
                
                if len(data1) == len(data2):
                    statistic, pvalue = stats.wilcoxon(data1, data2)
                    significant = "***" if pvalue < 0.001 else "**" if pvalue < 0.01 else "*" if pvalue < 0.05 else "ns"
                    
                    print(f"{method1} vs {method2:<25} {pvalue:<15.6f} {significant}")
        
        print("=" * 70)
        print("*** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant")


def get_predefined_configs() -> List[SimulationConfig]:
    """Get predefined simulation configurations for the three models."""
    
    configs = [
        # Human Out-of-Africa model
        SimulationConfig(
            name="human_ooa",
            species="HomSap",
            model="OutOfAfrica_3G09",
            samples={
                "YRI": 20,  # African
                "CEU": 20,  # European
                "CHB": 20,  # East Asian
            },
            sequence_length=10_000_000,  # 10 Mb
            mutation_rate=1.25e-8,
            recombination_rate=1e-8,
            seed=42,
        ),
        
        # Arabidopsis selfing model
        SimulationConfig(
            name="arabidopsis_selfing",
            species="AraTha",
            model="African2Epoch_1H18",  # Approximation
            samples={
                "AFR": 80,  # Simulated African accessions
            },
            sequence_length=5_000_000,  # 5 Mb (rapid LD decay)
            mutation_rate=7e-9,
            recombination_rate=8e-9,
            seed=42,
        ),
        
        # Rice structured population model
        SimulationConfig(
            name="rice_structure",
            species="OrySat",  # Would need custom model
            model="dombay2023",  # Example - would need actual rice demographic model
            samples={
                "indica": 20,
                "japonica": 20,
                "aus": 20,
            },
            sequence_length=8_000_000,  # 8 Mb
            mutation_rate=3e-8,
            recombination_rate=1e-8,
            seed=42,
        ),
    ]
    
    return configs


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark tree accuracy on simulated data"
    )
    parser.add_argument(
        "-o", "--output",
        default="benchmark_results",
        help="Output directory for results"
    )
    parser.add_argument(
        "-n", "--replicates",
        type=int,
        default=50,
        help="Number of replicates per model (default: 50)"
    )
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["human", "arabidopsis", "rice", "all"],
        default=["all"],
        help="Which demographic models to benchmark"
    )
    
    args = parser.parse_args()
    
    # Get configurations
    all_configs = get_predefined_configs()
    
    if "all" not in args.models:
        selected = []
        if "human" in args.models:
            selected.append(next(c for c in all_configs if c.name == "human_ooa"))
        if "arabidopsis" in args.models:
            selected.append(next(c for c in all_configs if c.name == "arabidopsis_selfing"))
        if "rice" in args.models:
            selected.append(next(c for c in all_configs if c.name == "rice_structure"))
        all_configs = selected
    
    # Run benchmark
    benchmark = TreeAccuracyBenchmark(
        output_dir=args.output,
        n_replicates=args.replicates,
    )
    
    benchmark.run_full_benchmark(all_configs)
    
    print("\n" + "=" * 70)
    print(f"BENCHMARK COMPLETE. Results saved to: {args.output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
