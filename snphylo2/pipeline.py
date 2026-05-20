"""
Main SNPhylo2 pipeline orchestration.
"""

import time
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field

from snphylo2.config import SNPhylo2Config
from snphylo2.io.vcf_reader import VCFReader
from snphylo2.filtering.variant_filters import FilterPipeline
from snphylo2.pruning.ld_pruner import LDPruner
from snphylo2.tree.tree_builder import TreeBuilder
from snphylo2.report.html_report import HTMLReport
from snphylo2.utils.logging_utils import get_logger
from snphylo2.exceptions import SNPhylo2Error

logger = get_logger()


@dataclass
class PipelineResults:
    """Results from pipeline execution."""
    output_dir: Path
    filtered_vcf: Optional[Path] = None
    pruned_vcf: Optional[Path] = None
    alignment_file: Optional[Path] = None
    tree_file: Optional[Path] = None
    report_path: Optional[Path] = None
    stats: Dict[str, Any] = field(default_factory=dict)
    runtime_seconds: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'output_dir': str(self.output_dir),
            'filtered_vcf': str(self.filtered_vcf) if self.filtered_vcf else None,
            'pruned_vcf': str(self.pruned_vcf) if self.pruned_vcf else None,
            'alignment_file': str(self.alignment_file) if self.alignment_file else None,
            'tree_file': str(self.tree_file) if self.tree_file else None,
            'report_path': str(self.report_path) if self.report_path else None,
            'stats': self.stats,
            'runtime_seconds': self.runtime_seconds,
        }


class Pipeline:
    """
    Main SNPhylo2 pipeline.
    
    Orchestrates the complete workflow:
    1. Quality control
    2. Variant filtering
    3. LD pruning
    4. Alignment generation
    5. Tree building
    6. Report generation
    """
    
    def __init__(self, config: SNPhylo2Config):
        """
        Initialize pipeline.
        
        Args:
            config: Pipeline configuration
        """
        self.config = config
        self.results = PipelineResults(output_dir=config.output.directory)
        
        # Create output directory
        self.config.output.directory.mkdir(parents=True, exist_ok=True)
    
    def run(self) -> Dict[str, Any]:
        """
        Run the complete pipeline.
        
        Returns:
            Dictionary with pipeline results
        """
        start_time = time.time()
        
        logger.info("=" * 60)
        logger.info("Starting SNPhylo2 Pipeline")
        logger.info("=" * 60)
        
        try:
            # Step 1: QC (optional, can be skipped)
            if self.config.reporting.format in ['html', 'both']:
                self._run_qc()
            
            # Step 2: Filtering
            self._run_filtering()
            
            # Step 3: LD Pruning
            self._run_pruning()
            
            # Step 4: Tree Building
            self._run_tree_building()
            
            # Step 5: Report Generation
            self._generate_report()
            
        except SNPhylo2Error as e:
            logger.error(f"Pipeline failed: {e.message}")
            raise
        except Exception as e:
            logger.error(f"Pipeline failed with unexpected error: {e}")
            raise SNPhylo2Error(f"Pipeline execution failed: {e}")
        
        # Calculate runtime
        self.results.runtime_seconds = time.time() - start_time
        
        logger.info("=" * 60)
        logger.info("Pipeline Complete")
        logger.info(f"Runtime: {self.results.runtime_seconds:.1f} seconds")
        logger.info(f"Output: {self.results.output_dir}")
        logger.info("=" * 60)
        
        return self.results.to_dict()
    
    def _run_qc(self) -> None:
        """Run quality control analysis."""
        logger.info("Step 1: Quality Control")
        
        # QC is primarily informational in the main pipeline
        # Detailed QC can be run separately with 'snphylo2 qc'
        with VCFReader(self.config.input.path) as reader:
            self.results.stats['input_samples'] = reader.n_samples
            
            # Count variants
            variant_count = 0
            for _ in reader:
                variant_count += 1
            
            self.results.stats['input_snps'] = variant_count
        
        logger.info(f"  Input samples: {self.results.stats['input_samples']}")
        logger.info(f"  Input SNPs: {self.results.stats['input_snps']:,}")
    
    def _run_filtering(self) -> None:
        """Run variant and sample filtering."""
        logger.info("Step 2: Variant Filtering")
        
        output_prefix = self.config.output.directory / self.config.output.prefix
        filtered_vcf = output_prefix.with_suffix('.filtered.vcf.gz')
        
        filter_pipeline = FilterPipeline(self.config.filtering)
        stats = filter_pipeline.run(
            self.config.input.path,
            filtered_vcf,
        )
        
        self.results.filtered_vcf = filtered_vcf
        self.results.stats['filtering'] = {
            'input_snps': stats.input_snps,
            'output_snps': stats.output_snps,
            'retention_rate': stats.retention_rate,
            'filters_applied': stats.filters_applied,
        }
        
        logger.info(f"  Filtered SNPs: {stats.output_snps:,} / {stats.input_snps:,}")
    
    def _run_pruning(self) -> None:
        """Run LD pruning."""
        logger.info("Step 3: LD Pruning")
        
        if not self.results.filtered_vcf:
            raise SNPhylo2Error("Cannot run pruning - no filtered VCF")
        
        output_prefix = self.config.output.directory / self.config.output.prefix
        pruned_vcf = output_prefix.with_suffix('.pruned.vcf.gz')
        
        pruner = LDPruner(self.config.ld_pruning, threads=self.config.compute.threads)
        stats = pruner.run(self.results.filtered_vcf, pruned_vcf)
        
        self.results.pruned_vcf = pruned_vcf
        self.results.stats['ld_pruning'] = {
            'input_snps': stats.input_snps,
            'output_snps': stats.output_snps,
            'pruned_count': stats.pruned_count,
            'retention_rate': stats.retention_rate,
        }
        
        logger.info(f"  Pruned SNPs: {stats.output_snps:,} / {stats.input_snps:,}")
    
    def _run_tree_building(self) -> None:
        """Run phylogenetic tree building."""
        logger.info("Step 4: Tree Building")
        
        if not self.results.pruned_vcf:
            raise SNPhylo2Error("Cannot build tree - no pruned VCF")
        
        output_prefix = self.config.output.directory / self.config.output.prefix
        tree_file = output_prefix.with_suffix('.tree.nwk')
        
        builder = TreeBuilder(self.config.tree)
        result = builder.build(self.results.pruned_vcf, tree_file)
        
        self.results.tree_file = tree_file
        self.results.stats['tree'] = {
            'engine': self.config.tree.engine.value,
            'model': result.get('model'),
            'log_likelihood': result.get('log_likelihood'),
            'bootstrap_support': result.get('bootstrap_support'),
        }
        
        logger.info(f"  Tree engine: {self.config.tree.engine.value}")
        logger.info(f"  Model: {result.get('model', 'N/A')}")
        if result.get('bootstrap_support'):
            logger.info(f"  Mean bootstrap: {result['bootstrap_support']:.1f}%")
    
    def _generate_report(self) -> None:
        """Generate HTML report."""
        logger.info("Step 5: Report Generation")
        
        if self.config.reporting.format == 'json':
            # JSON-only output
            import json
            report_path = self.config.output.directory / f"{self.config.output.prefix}.json"
            with open(report_path, 'w') as f:
                json.dump(self.results.to_dict(), f, indent=2)
            self.results.report_path = report_path
        else:
            # HTML report
            report_path = self.config.output.directory / f"{self.config.output.prefix}_report.html"
            
            reporter = HTMLReport(self.config.output.directory)
            reporter.generate(report_path, results=self.results)
            
            self.results.report_path = report_path
        
        logger.info(f"  Report: {self.results.report_path}")
    
    def _cleanup(self) -> None:
        """Clean up intermediate files if requested."""
        if not self.config.output.keep_intermediates:
            logger.info("Cleaning up intermediate files...")
            
            intermediates = [
                self.results.filtered_vcf,
                self.results.pruned_vcf,
                self.results.alignment_file,
            ]
            
            for file in intermediates:
                if file and file.exists():
                    try:
                        file.unlink()
                        logger.debug(f"  Removed: {file}")
                    except Exception as e:
                        logger.warning(f"  Could not remove {file}: {e}")


class PipelineStep:
    """Base class for pipeline steps."""
    
    def __init__(self, name: str, config: Any):
        self.name = name
        self.config = config
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None
    
    def __enter__(self):
        self.start_time = time.time()
        logger.info(f"Starting: {self.name}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        runtime = self.end_time - self.start_time
        
        if exc_type is None:
            logger.info(f"Completed: {self.name} ({runtime:.1f}s)")
        else:
            logger.error(f"Failed: {self.name} ({runtime:.1f}s)")
        
        return False  # Don't suppress exceptions
    
    def run(self) -> Any:
        """Execute the step. Must be implemented by subclasses."""
        raise NotImplementedError
