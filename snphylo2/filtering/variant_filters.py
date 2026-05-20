"""
Variant and sample filtering implementations.
"""

from pathlib import Path
from typing import Dict, List, Optional, Callable, Any, Set
from dataclasses import dataclass, field
from collections import defaultdict
import numpy as np

from snphylo2.io.vcf_reader import VCFReader, Variant
from snphylo2.config import FilteringConfig
from snphylo2.exceptions import FilterError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class FilterStats:
    """Statistics for filtering operations."""
    input_snps: int = 0
    output_snps: int = 0
    input_samples: int = 0
    output_samples: int = 0
    filters_applied: Dict[str, int] = field(default_factory=dict)
    
    @property
    def retention_rate(self) -> float:
        """Calculate SNP retention rate."""
        if self.input_snps == 0:
            return 0.0
        return self.output_snps / self.input_snps


class BaseFilter:
    """Base class for variant filters."""
    
    def __init__(self, name: str):
        self.name = name
        self.removed_count = 0
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        """
        Apply filter to variant.
        
        Args:
            variant: Variant to check
            sample_mask: Boolean mask for samples to include
            
        Returns:
            True if variant passes filter, False otherwise
        """
        raise NotImplementedError
    
    def get_stats(self) -> Dict[str, Any]:
        """Get filter statistics."""
        return {
            'name': self.name,
            'removed': self.removed_count,
        }


class MAFFilter(BaseFilter):
    """Filter by minor allele frequency."""
    
    def __init__(self, min_maf: float = 0.05, max_maf: float = 1.0):
        super().__init__("MAF")
        self.min_maf = min_maf
        self.max_maf = max_maf
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        genotypes = variant.genotypes
        
        if sample_mask is not None:
            genotypes = genotypes[sample_mask]
        
        # Calculate allele frequencies
        # genotypes: (n_samples, ploidy), values: 0, 1, 2, -1 (missing)
        valid_gt = genotypes[genotypes >= 0]
        
        if len(valid_gt) == 0:
            self.removed_count += 1
            return False
        
        # Count alleles
        n_alleles = len(valid_gt)
        alt_count = np.sum(valid_gt > 0)
        
        maf = min(alt_count / n_alleles, 1 - alt_count / n_alleles)
        
        if maf < self.min_maf or maf > self.max_maf:
            self.removed_count += 1
            return False
        
        return True


class MissingnessFilter(BaseFilter):
    """Filter by missing data rate."""
    
    def __init__(self, max_missing_rate: float = 0.2):
        super().__init__("Missingness")
        self.max_missing_rate = max_missing_rate
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        genotypes = variant.genotypes
        
        if sample_mask is not None:
            genotypes = genotypes[sample_mask]
        
        # Count missing (typically -1 or 3)
        missing = np.sum((genotypes < 0) | (genotypes == 3))
        total = genotypes.size
        
        missing_rate = missing / total if total > 0 else 1.0
        
        if missing_rate > self.max_missing_rate:
            self.removed_count += 1
            return False
        
        return True


class DepthFilter(BaseFilter):
    """Filter by read depth."""
    
    def __init__(self, min_depth: int = 5, max_depth: Optional[int] = None):
        super().__init__("Depth")
        self.min_depth = min_depth
        self.max_depth = max_depth
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        if variant.depths is None:
            # No depth info, assume pass
            return True
        
        depths = variant.depths
        
        if sample_mask is not None:
            depths = depths[sample_mask]
        
        # Check that at least some samples meet depth criteria
        valid_depths = depths[depths >= 0]
        
        if len(valid_depths) == 0:
            self.removed_count += 1
            return False
        
        # Require median depth to meet minimum
        median_depth = np.median(valid_depths)
        
        if median_depth < self.min_depth:
            self.removed_count += 1
            return False
        
        if self.max_depth and median_depth > self.max_depth:
            self.removed_count += 1
            return False
        
        return True


class QualityFilter(BaseFilter):
    """Filter by variant quality."""
    
    def __init__(self, min_qual: float = 30.0):
        super().__init__("Quality")
        self.min_qual = min_qual
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        if variant.qual is None or variant.qual < self.min_qual:
            self.removed_count += 1
            return False
        return True


class BiallelicFilter(BaseFilter):
    """Filter to keep only biallelic SNPs."""
    
    def __init__(self):
        super().__init__("Biallelic")
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        if not variant.is_biallelic:
            self.removed_count += 1
            return False
        return True


class TsTvFilter(BaseFilter):
    """Filter by transition/transversion ratio (informational only)."""
    
    def __init__(self):
        super().__init__("TsTv")
        self.transitions = 0
        self.transversions = 0
    
    def apply(self, variant: Variant, sample_mask: Optional[np.ndarray] = None) -> bool:
        # This filter doesn't remove variants, just counts
        if variant.is_snp and variant.is_biallelic:
            if variant.is_transition:
                self.transitions += 1
            elif variant.is_transversion:
                self.transversions += 1
        return True
    
    def get_ts_tv_ratio(self) -> float:
        """Calculate Ts/Tv ratio."""
        if self.transversions == 0:
            return float('inf') if self.transitions > 0 else 0.0
        return self.transitions / self.transversions


class FilterPipeline:
    """Pipeline for applying multiple filters."""
    
    def __init__(self, config: FilteringConfig):
        """
        Initialize filter pipeline.
        
        Args:
            config: Filtering configuration
        """
        self.config = config
        self.filters: List[BaseFilter] = []
        self._build_pipeline()
    
    def _build_pipeline(self) -> None:
        """Build filter chain from configuration."""
        # Add filters in order
        if self.config.biallelic_only:
            self.filters.append(BiallelicFilter())
        
        if self.config.maf.min > 0 or self.config.maf.max < 1.0:
            self.filters.append(MAFFilter(
                min_maf=self.config.maf.min,
                max_maf=self.config.maf.max,
            ))
        
        if self.config.missingness.max_per_snp < 1.0:
            self.filters.append(MissingnessFilter(
                max_missing_rate=self.config.missingness.max_per_snp,
            ))
        
        if self.config.depth.min > 0 or self.config.depth.max:
            self.filters.append(DepthFilter(
                min_depth=self.config.depth.min,
                max_depth=self.config.depth.max,
            ))
        
        if self.config.genotype_quality and self.config.genotype_quality > 0:
            self.filters.append(QualityFilter(min_qual=self.config.genotype_quality))
        
        # TsTv filter for statistics (doesn't remove)
        self.tstv_filter = TsTvFilter()
        self.filters.append(self.tstv_filter)
    
    def run(
        self,
        input_path: Path,
        output_path: Path,
        chunk_size: int = 10000,
    ) -> FilterStats:
        """
        Run filtering pipeline.
        
        Args:
            input_path: Input VCF file
            output_path: Output filtered VCF file
            chunk_size: Number of variants to process at once
            
        Returns:
            Filter statistics
        """
        logger.info(f"Starting filter pipeline")
        logger.info(f"  Input: {input_path}")
        logger.info(f"  Output: {output_path}")
        
        stats = FilterStats()
        stats.input_samples = self._count_samples(input_path)
        
        # For now, we'll write a simplified VCF
        # Full implementation would use pysam or cyvcf2 writer
        
        kept_variants = []
        
        with VCFReader(input_path) as reader:
            stats.input_snps = 0
            
            for chunk in reader.iter_chunks(chunk_size):
                stats.input_snps += len(chunk)
                
                for variant in chunk:
                    if self._apply_all_filters(variant):
                        kept_variants.append(variant)
                
                if stats.input_snps % 100000 == 0:
                    logger.debug(f"Processed {stats.input_snps} variants, kept {len(kept_variants)}")
        
        # Write output
        self._write_filtered_vcf(input_path, output_path, kept_variants)
        
        stats.output_snps = len(kept_variants)
        stats.output_samples = stats.input_samples
        
        for f in self.filters:
            stats.filters_applied[f.name] = f.removed_count
        
        logger.info(f"Filter pipeline complete")
        logger.info(f"  Input SNPs: {stats.input_snps:,}")
        logger.info(f"  Output SNPs: {stats.output_snps:,}")
        logger.info(f"  Retention: {stats.retention_rate:.1%}")
        logger.info(f"  Ts/Tv ratio: {self.tstv_filter.get_ts_tv_ratio():.3f}")
        
        return stats
    
    def _apply_all_filters(self, variant: Variant) -> bool:
        """Apply all filters to a variant."""
        for f in self.filters:
            if not f.apply(variant):
                return False
        return True
    
    def _count_samples(self, path: Path) -> int:
        """Count samples in VCF."""
        with VCFReader(path) as reader:
            return reader.n_samples
    
    def _write_filtered_vcf(
        self,
        input_path: Path,
        output_path: Path,
        variants: List[Variant],
    ) -> None:
        """Write filtered variants to output VCF."""
        # Simplified implementation - would use proper VCF writer
        # For now, just log the count
        logger.info(f"Would write {len(variants)} variants to {output_path}")
        
        # Placeholder for actual implementation
        # This would use pysam.VariantFile to write proper VCF
        pass


class SampleQC:
    """Quality control for samples."""
    
    def __init__(self, max_missing_rate: float = 0.5, min_depth: int = 5):
        self.max_missing_rate = max_missing_rate
        self.min_depth = min_depth
    
    def calculate_metrics(
        self,
        input_path: Path,
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate QC metrics per sample.
        
        Returns:
            Dictionary mapping sample IDs to metrics
        """
        metrics = defaultdict(lambda: {
            'total_snps': 0,
            'missing_snps': 0,
            'mean_depth': 0.0,
            'het_rate': 0.0,
        })
        
        with VCFReader(input_path) as reader:
            sample_names = reader.sample_names
            
            for chunk in reader.iter_chunks(10000):
                for variant in chunk:
                    for i, sample in enumerate(sample_names):
                        metrics[sample]['total_snps'] += 1
                        
                        gt = variant.genotypes[i]
                        if np.any(gt < 0) or np.any(gt == 3):
                            metrics[sample]['missing_snps'] += 1
                        
                        # Additional metrics would be calculated here
        
        # Calculate rates
        for sample in metrics:
            total = metrics[sample]['total_snps']
            if total > 0:
                metrics[sample]['missing_rate'] = metrics[sample]['missing_snps'] / total
            else:
                metrics[sample]['missing_rate'] = 1.0
        
        return dict(metrics)
    
    def identify_outliers(
        self,
        metrics: Dict[str, Dict[str, float]],
    ) -> List[str]:
        """
        Identify outlier samples based on QC metrics.
        
        Returns:
            List of outlier sample IDs
        """
        outliers = []
        
        for sample, m in metrics.items():
            if m.get('missing_rate', 0) > self.max_missing_rate:
                outliers.append(sample)
        
        return outliers
