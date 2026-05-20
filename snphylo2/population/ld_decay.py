"""
LD (Linkage Disequilibrium) decay analysis module.

Inspired by PopLDdecay (BGI), calculates LD decay patterns for population genetics.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from collections import defaultdict
import itertools

import numpy as np
import pandas as pd

try:
    import scipy.stats as stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

from snphylo2.io.vcf_reader import VCFReader, Variant
from snphylo2.exceptions import PopulationError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class LDDecayResult:
    """Results from LD decay analysis."""
    distances: np.ndarray  # Physical distances (bp)
    mean_r2: np.ndarray    # Mean r² at each distance
    percentiles_95: np.ndarray  # 95th percentile
    percentiles_5: np.ndarray   # 5th percentile
    n_pairs: np.ndarray    # Number of SNP pairs at each distance
    half_decay_distance: Optional[float] = None  # Distance where r² = 0.5
    population: Optional[str] = None
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        return pd.DataFrame({
            'distance': self.distances,
            'mean_r2': self.mean_r2,
            'p95_r2': self.percentiles_95,
            'p5_r2': self.percentiles_5,
            'n_pairs': self.n_pairs,
        })
    
    def calculate_half_decay(self) -> float:
        """Calculate distance at which r² decays to 0.5."""
        if np.all(self.mean_r2 < 0.5):
            return self.distances[0]
        if np.all(self.mean_r2 > 0.5):
            return self.distances[-1]
        
        # Interpolate to find where r² = 0.5
        for i in range(len(self.mean_r2) - 1):
            if self.mean_r2[i] >= 0.5 >= self.mean_r2[i + 1]:
                # Linear interpolation
                x1, y1 = self.distances[i], self.mean_r2[i]
                x2, y2 = self.distances[i + 1], self.mean_r2[i + 1]
                x = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
                self.half_decay_distance = x
                return x
        
        return self.distances[-1]


class LDDecayAnalyzer:
    """
    LD decay analysis inspired by PopLDdecay.
    
    Calculates r² decay with physical distance for population structure analysis.
    """
    
    def __init__(
        self,
        max_distance: int = 500_000,  # 500 kb
        bin_size: int = 5_000,        # 5 kb bins
        min_maf: float = 0.05,
        max_missing: float = 0.2,
    ):
        """
        Initialize LD decay analyzer.
        
        Args:
            max_distance: Maximum distance to calculate LD (bp)
            bin_size: Bin size for distance aggregation (bp)
            min_maf: Minimum MAF for included SNPs
            max_missing: Maximum missing rate for included SNPs
        """
        self.max_distance = max_distance
        self.bin_size = bin_size
        self.min_maf = min_maf
        self.max_missing = max_missing
        
        # Distance bins
        self.bins = np.arange(0, max_distance + bin_size, bin_size)
        self.bin_centers = (self.bins[:-1] + self.bins[1:]) / 2
    
    def analyze(
        self,
        vcf_path: Path,
        sample_groups: Optional[Dict[str, List[str]]] = None,
        chromosomes: Optional[List[str]] = None,
    ) -> Dict[str, LDDecayResult]:
        """
        Analyze LD decay for a VCF file.
        
        Args:
            vcf_path: Path to VCF file
            sample_groups: Optional dict mapping group names to sample IDs
            chromosomes: Optional list of chromosomes to analyze
            
        Returns:
            Dictionary mapping group names to LDDecayResult
        """
        logger.info(f"Analyzing LD decay: {vcf_path}")
        logger.info(f"  Max distance: {self.max_distance:,} bp")
        logger.info(f"  Bin size: {self.bin_size:,} bp")
        
        # Load variants
        variants_by_chrom = self._load_variants(vcf_path, chromosomes)
        
        results = {}
        
        if sample_groups is None:
            # Analyze all samples together
            logger.info("Analyzing all samples as single group")
            result = self._analyze_group(variants_by_chrom, group_name="all")
            results["all"] = result
        else:
            # Analyze each group separately
            for group_name, sample_ids in sample_groups.items():
                logger.info(f"Analyzing group: {group_name} ({len(sample_ids)} samples)")
                result = self._analyze_group(
                    variants_by_chrom,
                    group_name=group_name,
                    sample_subset=sample_ids,
                )
                results[group_name] = result
        
        return results
    
    def _load_variants(
        self,
        vcf_path: Path,
        chromosomes: Optional[List[str]] = None,
    ) -> Dict[str, List[Variant]]:
        """Load variants by chromosome."""
        variants_by_chrom = defaultdict(list)
        
        with VCFReader(vcf_path) as reader:
            for variant in reader:
                # Filter by chromosome
                if chromosomes and variant.chrom not in chromosomes:
                    continue
                
                # Filter by MAF and missingness
                if not self._pass_filters(variant):
                    continue
                
                variants_by_chrom[variant.chrom].append(variant)
        
        # Sort by position
        for chrom in variants_by_chrom:
            variants_by_chrom[chrom].sort(key=lambda v: v.pos)
        
        logger.info(f"Loaded {sum(len(v) for v in variants_by_chrom.values())} variants")
        
        return dict(variants_by_chrom)
    
    def _pass_filters(self, variant: Variant) -> bool:
        """Check if variant passes quality filters."""
        # Must be biallelic
        if not variant.is_biallelic:
            return False
        
        # Calculate MAF
        genotypes = variant.genotypes
        valid_gt = genotypes[genotypes >= 0]
        
        if len(valid_gt) == 0:
            return False
        
        n_alleles = len(valid_gt)
        alt_count = np.sum(valid_gt > 0)
        maf = min(alt_count / n_alleles, 1 - alt_count / n_alleles)
        
        if maf < self.min_maf:
            return False
        
        # Calculate missing rate
        missing = np.sum((genotypes < 0) | (genotypes == 3))
        missing_rate = missing / genotypes.size
        
        if missing_rate > self.max_missing:
            return False
        
        return True
    
    def _analyze_group(
        self,
        variants_by_chrom: Dict[str, List[Variant]],
        group_name: str,
        sample_subset: Optional[List[str]] = None,
    ) -> LDDecayResult:
        """Analyze LD decay for a group of samples."""
        # Accumulate r² values by distance bin
        bin_r2_values = defaultdict(list)
        
        # Process each chromosome
        for chrom, variants in variants_by_chrom.items():
            if len(variants) < 2:
                continue
            
            logger.debug(f"Processing {chrom}: {len(variants)} variants")
            
            # Calculate LD for variant pairs
            for i, var1 in enumerate(variants):
                # Only compare with variants within max_distance
                for var2 in variants[i+1:]:
                    distance = var2.pos - var1.pos
                    
                    if distance > self.max_distance:
                        break
                    
                    # Calculate r²
                    r2 = self._calculate_r2(var1, var2, sample_subset)
                    
                    if r2 is not None:
                        bin_idx = int(distance / self.bin_size)
                        bin_r2_values[bin_idx].append(r2)
        
        # Aggregate results
        mean_r2 = np.zeros(len(self.bin_centers))
        p95_r2 = np.zeros(len(self.bin_centers))
        p5_r2 = np.zeros(len(self.bin_centers))
        n_pairs = np.zeros(len(self.bin_centers), dtype=int)
        
        for bin_idx in range(len(self.bin_centers)):
            if bin_idx in bin_r2_values and len(bin_r2_values[bin_idx]) > 0:
                values = np.array(bin_r2_values[bin_idx])
                mean_r2[bin_idx] = np.mean(values)
                p95_r2[bin_idx] = np.percentile(values, 95)
                p5_r2[bin_idx] = np.percentile(values, 5)
                n_pairs[bin_idx] = len(values)
        
        # Filter bins with sufficient data
        valid_bins = n_pairs >= 10
        
        result = LDDecayResult(
            distances=self.bin_centers[valid_bins],
            mean_r2=mean_r2[valid_bins],
            percentiles_95=p95_r2[valid_bins],
            percentiles_5=p5_r2[valid_bins],
            n_pairs=n_pairs[valid_bins],
            population=group_name,
        )
        
        # Calculate half-decay distance
        result.calculate_half_decay()
        
        logger.info(f"Group {group_name}: Half-decay distance = {result.half_decay_distance:.0f} bp")
        
        return result
    
    def _calculate_r2(
        self,
        var1: Variant,
        var2: Variant,
        sample_subset: Optional[List[str]] = None,
    ) -> Optional[float]:
        """
        Calculate r² (correlation coefficient squared) between two variants.
        
        Args:
            var1: First variant
            var2: Second variant
            sample_subset: Optional list of sample indices to include
            
        Returns:
            r² value or None if insufficient data
        """
        # Get genotypes (0, 1, 2)
        g1 = np.sum(var1.genotypes, axis=1)
        g2 = np.sum(var2.genotypes, axis=1)
        
        # Handle missing data
        valid = (g1 >= 0) & (g2 >= 0) & (g1 != 3) & (g2 != 3)
        
        if sample_subset is not None:
            # Filter to subset samples (would need proper implementation)
            pass
        
        g1_valid = g1[valid]
        g2_valid = g2[valid]
        
        if len(g1_valid) < 10:  # Need sufficient sample size
            return None
        
        # Calculate Pearson correlation
        if SCIPY_AVAILABLE:
            r, _ = stats.pearsonr(g1_valid, g2_valid)
        else:
            # Manual calculation
            mean1, mean2 = np.mean(g1_valid), np.mean(g2_valid)
            std1, std2 = np.std(g1_valid), np.std(g2_valid)
            
            if std1 == 0 or std2 == 0:
                return None
            
            r = np.mean((g1_valid - mean1) * (g2_valid - mean2)) / (std1 * std2)
        
        if np.isnan(r):
            return None
        
        return r ** 2
    
    def compare_populations(
        self,
        results: Dict[str, LDDecayResult],
    ) -> pd.DataFrame:
        """
        Compare LD decay across populations.
        
        Args:
            results: Dictionary of LDDecayResult by population
            
        Returns:
            DataFrame with comparison statistics
        """
        comparison = []
        
        for pop_name, result in results.items():
            comparison.append({
                'population': pop_name,
                'half_decay_distance': result.half_decay_distance,
                'mean_r2_first_bin': result.mean_r2[0] if len(result.mean_r2) > 0 else np.nan,
                'max_distance_analyzed': result.distances[-1] if len(result.distances) > 0 else 0,
                'total_pairs': np.sum(result.n_pairs),
            })
        
        return pd.DataFrame(comparison)


class LDHeatmapGenerator:
    """Generate LD heatmaps for genomic regions."""
    
    def __init__(self, window_size: int = 100_000):
        """
        Initialize LD heatmap generator.
        
        Args:
            window_size: Window size for heatmap (bp)
        """
        self.window_size = window_size
    
    def generate_heatmap(
        self,
        vcf_path: Path,
        chrom: str,
        start: int,
        end: int,
        sample_subset: Optional[List[str]] = None,
    ) -> Tuple[np.ndarray, List[int]]:
        """
        Generate LD heatmap for a genomic region.
        
        Args:
            vcf_path: Path to VCF file
            chrom: Chromosome
            start: Start position
            end: End position
            sample_subset: Optional sample subset
            
        Returns:
            Tuple of (r2_matrix, positions)
        """
        # Load variants in region
        variants = []
        positions = []
        
        with VCFReader(vcf_path) as reader:
            for variant in reader.iter_variants(chrom, start, end):
                if not variant.is_biallelic:
                    continue
                
                variants.append(variant)
                positions.append(variant.pos)
        
        if len(variants) < 2:
            raise PopulationError(f"Insufficient variants in region {chrom}:{start}-{end}")
        
        # Calculate r² matrix
        n = len(variants)
        r2_matrix = np.zeros((n, n))
        
        analyzer = LDDecayAnalyzer()
        
        for i in range(n):
            for j in range(i, n):
                r2 = analyzer._calculate_r2(variants[i], variants[j], sample_subset)
                if r2 is not None:
                    r2_matrix[i, j] = r2
                    r2_matrix[j, i] = r2
        
        return r2_matrix, positions
