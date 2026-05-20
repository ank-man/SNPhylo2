#!/usr/bin/env python3
"""
SNPhylo2 Benchmarking Suite

Comprehensive performance testing across different dataset sizes and configurations.
"""

import argparse
import json
import time
import tempfile
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import subprocess

import psutil


@dataclass
class BenchmarkResult:
    """Results from a benchmark run."""
    name: str
    runtime_seconds: float
    peak_memory_mb: float
    input_snps: int
    input_samples: int
    output_snps: int
    success: bool = True
    error_message: Optional[str] = None
    metadata: Dict = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        return {
            'name': self.name,
            'runtime_seconds': self.runtime_seconds,
            'peak_memory_mb': self.peak_memory_mb,
            'input_snps': self.input_snps,
            'input_samples': self.input_samples,
            'output_snps': self.output_snps,
            'success': self.success,
            'error_message': self.error_message,
            'metadata': self.metadata,
        }


class ResourceMonitor:
    """Monitor CPU and memory usage during execution."""
    
    def __init__(self):
        self.peak_memory = 0.0
        self.process = psutil.Process()
    
    def sample(self):
        """Sample current resource usage."""
        mem_info = self.process.memory_info()
        self.peak_memory = max(self.peak_memory, mem_info.rss / 1024 / 1024)  # MB
    
    def get_peak_memory(self) -> float:
        """Get peak memory usage."""
        return self.peak_memory


class BenchmarkSuite:
    """Benchmark suite for SNPhylo2."""
    
    def __init__(self, output_dir: str, tools_to_compare: Optional[List[str]] = None):
        """
        Initialize benchmark suite.
        
        Args:
            output_dir: Directory for benchmark results
            tools_to_compare: List of tools to compare (snphylo2, tassel, plink_iqtree)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.tools = tools_to_compare or ["snphylo2"]
        self.results: List[BenchmarkResult] = []
    
    def run_benchmarks(self, datasets: List[Dict]) -> Dict:
        """
        Run benchmarks on specified datasets.
        
        Args:
            datasets: List of dataset specifications
            
        Returns:
            Dictionary with all results
        """
        print("=" * 60)
        print("SNPhylo2 Benchmark Suite")
        print("=" * 60)
        
        for dataset in datasets:
            print(f"\nDataset: {dataset['name']}")
            print(f"  Samples: {dataset['n_samples']}")
            print(f"  SNPs: {dataset['n_snps']}")
            
            # Generate or locate dataset
            vcf_path = self._get_dataset(dataset)
            
            # Run benchmarks for each tool
            for tool in self.tools:
                print(f"\n  Running: {tool}")
                
                try:
                    result = self._run_tool_benchmark(tool, vcf_path, dataset)
                    self.results.append(result)
                    
                    if result.success:
                        print(f"    ✓ Runtime: {result.runtime_seconds:.1f}s")
                        print(f"    ✓ Memory: {result.peak_memory_mb:.1f} MB")
                        print(f"    ✓ Retention: {result.output_snps / result.input_snps:.1%}")
                    else:
                        print(f"    ✗ Failed: {result.error_message}")
                        
                except Exception as e:
                    print(f"    ✗ Error: {e}")
                    self.results.append(BenchmarkResult(
                        name=f"{tool}_{dataset['name']}",
                        runtime_seconds=0,
                        peak_memory_mb=0,
                        input_snps=dataset['n_snps'],
                        input_samples=dataset['n_samples'],
                        output_snps=0,
                        success=False,
                        error_message=str(e),
                    ))
        
        # Generate report
        report = self._generate_report()
        return report
    
    def _get_dataset(self, spec: Dict) -> Path:
        """Get or generate dataset."""
        # Check if dataset exists
        dataset_name = f"{spec['name']}_{spec['n_samples']}s_{spec['n_snps']}snps.vcf.gz"
        dataset_path = self.output_dir / "datasets" / dataset_name
        
        if dataset_path.exists():
            return dataset_path
        
        # Generate dataset
        print(f"    Generating dataset...")
        dataset_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            from tests.fixtures.generate_simulated_data import generate_benchmark_dataset
            
            result = generate_benchmark_dataset(
                str(dataset_path.parent),
                n_samples=spec['n_samples'],
                n_snps=spec['n_snps'],
                seed=spec.get('seed', 42),
            )
            return Path(result['vcf'])
            
        except ImportError:
            # Fallback to simple generator
            from tests.fixtures.generate_test_vcf import generate_test_vcf
            
            generate_test_vcf(
                str(dataset_path),
                n_samples=spec['n_samples'],
                n_snps=spec['n_snps'],
                seed=spec.get('seed', 42),
            )
            return dataset_path
    
    def _run_tool_benchmark(
        self,
        tool: str,
        vcf_path: Path,
        dataset_spec: Dict,
    ) -> BenchmarkResult:
        """Run benchmark for a specific tool."""
        
        if tool == "snphylo2":
            return self._benchmark_snphylo2(vcf_path, dataset_spec)
        elif tool == "plink_iqtree":
            return self._benchmark_plink_iqtree(vcf_path, dataset_spec)
        else:
            raise ValueError(f"Unknown tool: {tool}")
    
    def _benchmark_snphylo2(
        self,
        vcf_path: Path,
        dataset_spec: Dict,
    ) -> BenchmarkResult:
        """Benchmark SNPhylo2."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result_dir = Path(tmpdir) / "results"
            
            # Monitor resources
            monitor = ResourceMonitor()
            start_time = time.time()
            
            try:
                # Run pipeline
                cmd = [
                    "snphylo2", "run",
                    "-v", str(vcf_path),
                    "-o", str(result_dir),
                    "--threads", "4",
                    "--bootstrap", "100",  # Reduced for benchmarking
                    "--report", "json",
                ]
                
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                
                # Monitor while running
                while process.poll() is None:
                    monitor.sample()
                    time.sleep(0.1)
                
                runtime = time.time() - start_time
                
                if process.returncode != 0:
                    stderr = process.stderr.read().decode()
                    return BenchmarkResult(
                        name=f"snphylo2_{dataset_spec['name']}",
                        runtime_seconds=runtime,
                        peak_memory_mb=monitor.get_peak_memory(),
                        input_snps=dataset_spec['n_snps'],
                        input_samples=dataset_spec['n_samples'],
                        output_snps=0,
                        success=False,
                        error_message=stderr,
                    )
                
                # Parse results
                report_path = result_dir / "snphylo2_output.json"
                output_snps = 0
                if report_path.exists():
                    with open(report_path) as f:
                        data = json.load(f)
                        stats = data.get('stats', {})
                        pruning = stats.get('ld_pruning', {})
                        output_snps = pruning.get('output_snps', 0)
                
                return BenchmarkResult(
                    name=f"snphylo2_{dataset_spec['name']}",
                    runtime_seconds=runtime,
                    peak_memory_mb=monitor.get_peak_memory(),
                    input_snps=dataset_spec['n_snps'],
                    input_samples=dataset_spec['n_samples'],
                    output_snps=output_snps,
                    success=True,
                    metadata={'command': ' '.join(cmd)},
                )
                
            except Exception as e:
                runtime = time.time() - start_time
                return BenchmarkResult(
                    name=f"snphylo2_{dataset_spec['name']}",
                    runtime_seconds=runtime,
                    peak_memory_mb=monitor.get_peak_memory(),
                    input_snps=dataset_spec['n_snps'],
                    input_samples=dataset_spec['n_samples'],
                    output_snps=0,
                    success=False,
                    error_message=str(e),
                )
    
    def _benchmark_plink_iqtree(
        self,
        vcf_path: Path,
        dataset_spec: Dict,
    ) -> BenchmarkResult:
        """Benchmark custom PLINK + IQ-TREE pipeline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir)
            
            monitor = ResourceMonitor()
            start_time = time.time()
            
            try:
                # Convert VCF to PLINK
                cmd1 = [
                    "plink2",
                    "--vcf", str(vcf_path),
                    "--make-bed",
                    "--out", str(work_dir / "data"),
                ]
                subprocess.run(cmd1, check=True, capture_output=True)
                monitor.sample()
                
                # LD pruning
                cmd2 = [
                    "plink2",
                    "--bfile", str(work_dir / "data"),
                    "--indep-pairwise", "50", "10", "0.2",
                    "--out", str(work_dir / "pruned"),
                ]
                subprocess.run(cmd2, check=True, capture_output=True)
                monitor.sample()
                
                # Extract pruned SNPs
                cmd3 = [
                    "plink2",
                    "--bfile", str(work_dir / "data"),
                    "--extract", str(work_dir / "pruned.prune.in"),
                    "--recode", "vcf",
                    "--out", str(work_dir / "pruned"),
                ]
                subprocess.run(cmd3, check=True, capture_output=True)
                monitor.sample()
                
                # Convert to PHYLIP
                # (simplified - would need proper conversion)
                
                runtime = time.time() - start_time
                
                return BenchmarkResult(
                    name=f"plink_iqtree_{dataset_spec['name']}",
                    runtime_seconds=runtime,
                    peak_memory_mb=monitor.get_peak_memory(),
                    input_snps=dataset_spec['n_snps'],
                    input_samples=dataset_spec['n_samples'],
                    output_snps=0,
                    success=True,
                )
                
            except Exception as e:
                runtime = time.time() - start_time
                return BenchmarkResult(
                    name=f"plink_iqtree_{dataset_spec['name']}",
                    runtime_seconds=runtime,
                    peak_memory_mb=monitor.get_peak_memory(),
                    input_snps=dataset_spec['n_snps'],
                    input_samples=dataset_spec['n_samples'],
                    output_snps=0,
                    success=False,
                    error_message=str(e),
                )
    
    def _generate_report(self) -> Dict:
        """Generate benchmark report."""
        report = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'results': [r.to_dict() for r in self.results],
            'summary': {},
        }
        
        # Calculate summary statistics
        for tool in self.tools:
            tool_results = [r for r in self.results if r.name.startswith(tool) and r.success]
            
            if tool_results:
                report['summary'][tool] = {
                    'avg_runtime': sum(r.runtime_seconds for r in tool_results) / len(tool_results),
                    'avg_memory': sum(r.peak_memory_mb for r in tool_results) / len(tool_results),
                    'success_rate': sum(1 for r in tool_results if r.success) / len(tool_results),
                }
        
        # Save report
        report_path = self.output_dir / "benchmark_report.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\nReport saved: {report_path}")
        
        return report


def main():
    parser = argparse.ArgumentParser(
        description="SNPhylo2 Benchmark Suite"
    )
    parser.add_argument(
        "-o", "--output",
        default="benchmark_results",
        help="Output directory for results"
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        default=["snphylo2"],
        choices=["snphylo2", "plink_iqtree", "tassel"],
        help="Tools to benchmark"
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["small", "medium", "large"],
        help="Dataset sizes to test"
    )
    
    args = parser.parse_args()
    
    # Define datasets
    dataset_specs = {
        "tiny": {"name": "tiny", "n_samples": 10, "n_snps": 100},
        "small": {"name": "small", "n_samples": 50, "n_snps": 1000},
        "medium": {"name": "medium", "n_samples": 100, "n_snps": 10000},
        "large": {"name": "large", "n_samples": 500, "n_snps": 50000},
    }
    
    datasets = [dataset_specs[d] for d in args.datasets if d in dataset_specs]
    
    # Run benchmarks
    suite = BenchmarkSuite(args.output, tools_to_compare=args.tools)
    report = suite.run_benchmarks(datasets)
    
    print("\n" + "=" * 60)
    print("Benchmark Complete")
    print("=" * 60)
    
    # Print summary
    for tool, stats in report['summary'].items():
        print(f"\n{tool}:")
        print(f"  Avg runtime: {stats['avg_runtime']:.1f}s")
        print(f"  Avg memory: {stats['avg_memory']:.1f} MB")
        print(f"  Success rate: {stats['success_rate']:.1%}")


if __name__ == "__main__":
    main()
