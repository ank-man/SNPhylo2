"""
Quality control analysis for SNP datasets.
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np

from snphylo2.io.vcf_reader import VCFReader, Variant
from snphylo2.exceptions import InputError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class QCMetrics:
    """Quality control metrics for a dataset."""
    # Overall metrics
    total_snps: int = 0
    total_samples: int = 0
    
    # SNP metrics
    snp_missing_rates: List[float] = field(default_factory=list)
    maf_distribution: List[float] = field(default_factory=list)
    depth_distribution: List[float] = field(default_factory=list)
    qual_distribution: List[float] = field(default_factory=list)
    
    # Sample metrics
    sample_missing_rates: Dict[str, float] = field(default_factory=dict)
    sample_heterozygosity: Dict[str, float] = field(default_factory=dict)
    sample_mean_depth: Dict[str, float] = field(default_factory=dict)
    
    # Type metrics
    transitions: int = 0
    transversions: int = 0
    
    @property
    def overall_missing_rate(self) -> float:
        if not self.snp_missing_rates:
            return 0.0
        return np.mean(self.snp_missing_rates)
    
    @property
    def mean_maf(self) -> float:
        if not self.maf_distribution:
            return 0.0
        return np.mean(self.maf_distribution)
    
    @property
    def ts_tv_ratio(self) -> float:
        if self.transversions == 0:
            return float('inf')
        return self.transitions / self.transversions
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'total_snps': self.total_snps,
            'total_samples': self.total_samples,
            'overall_missing_rate': self.overall_missing_rate,
            'mean_maf': self.mean_maf,
            'ts_tv_ratio': self.ts_tv_ratio,
            'transitions': self.transitions,
            'transversions': self.transversions,
        }


@dataclass
class QCReport:
    """Quality control report."""
    metrics: QCMetrics
    sample_outliers: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    
    def save(self, output_path: Path) -> None:
        """Save report to file."""
        # Placeholder - would generate HTML/JSON
        import json
        with open(output_path, 'w') as f:
            json.dump(self.metrics.to_dict(), f, indent=2)


class QCModule:
    """
    Quality control analysis for VCF datasets.
    """
    
    def __init__(self, vcf_path: Path, chunk_size: int = 10000):
        """
        Initialize QC module.
        
        Args:
            vcf_path: Path to VCF file
            chunk_size: Number of variants to process at once
        """
        self.vcf_path = Path(vcf_path)
        self.chunk_size = chunk_size
        self.metrics = QCMetrics()
    
    def run_analysis(self) -> QCReport:
        """
        Run complete QC analysis.
        
        Returns:
            QCReport with metrics and recommendations
        """
        logger.info(f"Starting QC analysis: {self.vcf_path}")
        
        self._collect_metrics()
        outliers = self._identify_outliers()
        recommendations = self._generate_recommendations()
        
        report = QCReport(
            metrics=self.metrics,
            sample_outliers=outliers,
            recommendations=recommendations,
        )
        
        logger.info(f"QC analysis complete")
        logger.info(f"  SNPs: {self.metrics.total_snps:,}")
        logger.info(f"  Samples: {self.metrics.total_samples}")
        logger.info(f"  Missing rate: {self.metrics.overall_missing_rate:.2%}")
        logger.info(f"  Ts/Tv: {self.metrics.ts_tv_ratio:.3f}")
        
        return report
    
    def _collect_metrics(self) -> None:
        """Collect QC metrics from VCF."""
        with VCFReader(self.vcf_path) as reader:
            self.metrics.total_samples = reader.n_samples
            sample_names = reader.sample_names
            
            # Initialize sample metrics
            sample_missing_counts = defaultdict(int)
            sample_het_counts = defaultdict(int)
            sample_gt_counts = defaultdict(int)
            sample_depth_sums = defaultdict(float)
            sample_depth_counts = defaultdict(int)
            
            # Process variants in chunks
            variant_count = 0
            
            for chunk in reader.iter_chunks(self.chunk_size):
                for variant in chunk:
                    variant_count += 1
                    
                    # SNP-level metrics
                    missing_count = 0
                    alt_allele_count = 0
                    total_alleles = 0
                    
                    for i, sample in enumerate(sample_names):
                        gt = variant.genotypes[i]
                        sample_gt_counts[sample] += 1
                        
                        # Missing
                        if np.any(gt < 0) or np.any(gt == 3):
                            missing_count += 1
                            sample_missing_counts[sample] += 1
                        else:
                            # Count alleles
                            total_alleles += len(gt)
                            alt_allele_count += np.sum(gt > 0)
                            
                            # Heterozygosity
                            if np.any(gt == 1):
                                sample_het_counts[sample] += 1
                        
                        # Depth
                        if variant.depths is not None and i < len(variant.depths):
                            dp = variant.depths[i]
                            if dp >= 0:
                                sample_depth_sums[sample] += dp
                                sample_depth_counts[sample] += 1
                    
                    # Store SNP metrics
                    missing_rate = missing_count / len(sample_names)
                    self.metrics.snp_missing_rates.append(missing_rate)
                    
                    # MAF
                    if total_alleles > 0:
                        maf = min(alt_allele_count / total_alleles, 
                                 1 - alt_allele_count / total_alleles)
                        self.metrics.maf_distribution.append(maf)
                    
                    # Quality
                    if variant.qual is not None:
                        self.metrics.qual_distribution.append(variant.qual)
                    
                    # Ts/Tv
                    if variant.is_snp and variant.is_biallelic:
                        if variant.is_transition:
                            self.metrics.transitions += 1
                        elif variant.is_transversion:
                            self.metrics.transversions += 1
            
            self.metrics.total_snps = variant_count
            
            # Calculate sample metrics
            for sample in sample_names:
                if sample_gt_counts[sample] > 0:
                    self.metrics.sample_missing_rates[sample] = (
                        sample_missing_counts[sample] / sample_gt_counts[sample]
                    )
                    self.metrics.sample_heterozygosity[sample] = (
                        sample_het_counts[sample] / sample_gt_counts[sample]
                    )
                
                if sample_depth_counts[sample] > 0:
                    self.metrics.sample_mean_depth[sample] = (
                        sample_depth_sums[sample] / sample_depth_counts[sample]
                    )
    
    def _identify_outliers(self) -> List[str]:
        """Identify outlier samples."""
        outliers = []
        
        # Missingness outliers (>3 SD from mean)
        if self.metrics.sample_missing_rates:
            missing_rates = list(self.metrics.sample_missing_rates.values())
            mean_missing = np.mean(missing_rates)
            std_missing = np.std(missing_rates)
            
            for sample, rate in self.metrics.sample_missing_rates.items():
                if rate > mean_missing + 3 * std_missing:
                    outliers.append(f"{sample} (high missingness: {rate:.1%})")
        
        # Heterozygosity outliers
        if self.metrics.sample_heterozygosity:
            het_rates = list(self.metrics.sample_heterozygosity.values())
            mean_het = np.mean(het_rates)
            std_het = np.std(het_rates)
            
            for sample, rate in self.metrics.sample_heterozygosity.items():
                if abs(rate - mean_het) > 3 * std_het:
                    outliers.append(f"{sample} (unusual heterozygosity: {rate:.1%})")
        
        return outliers
    
    def _generate_recommendations(self) -> List[str]:
        """Generate filtering recommendations based on QC."""
        recommendations = []
        
        # Missingness recommendations
        if self.metrics.overall_missing_rate > 0.5:
            recommendations.append(
                "Dataset has high overall missingness (>50%). "
                "Consider --max-missing 0.5 or imputation."
            )
        
        # MAF recommendations
        if self.metrics.mean_maf < 0.05:
            recommendations.append(
                "Dataset has low mean MAF. Consider --maf 0.01 for rare variants."
            )
        
        # Ts/Tv recommendations
        ts_tv = self.metrics.ts_tv_ratio
        if ts_tv < 1.0:
            recommendations.append(
                f"Low Ts/Tv ratio ({ts_tv:.2f}). May indicate quality issues."
            )
        elif ts_tv > 3.0:
            recommendations.append(
                f"High Ts/Tv ratio ({ts_tv:.2f}). Expected ~2.0 for genome-wide SNPs."
            )
        
        return recommendations
    
    def get_sample_qc_table(self) -> Dict[str, Dict[str, float]]:
        """
        Get QC metrics table for all samples.
        
        Returns:
            Dictionary mapping sample IDs to metrics
        """
        table = {}
        for sample in self.metrics.sample_missing_rates:
            table[sample] = {
                'missing_rate': self.metrics.sample_missing_rates.get(sample, 0),
                'heterozygosity': self.metrics.sample_heterozygosity.get(sample, 0),
                'mean_depth': self.metrics.sample_mean_depth.get(sample, 0),
            }
        return table
