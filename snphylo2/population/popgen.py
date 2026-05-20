"""
Population genetics analyses for SNPhylo2.

Implements PCA, FST, IBS, and related analyses.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from collections import defaultdict

import numpy as np
import pandas as pd

from snphylo2.io.vcf_reader import VCFReader
from snphylo2.exceptions import PopulationError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class PCAResult:
    """Results from PCA analysis."""
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    explained_variance_ratio: np.ndarray
    sample_ids: List[str]
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        n_components = len(self.eigenvalues)
        
        data = {
            'sample_id': self.sample_ids,
        }
        
        for i in range(n_components):
            data[f'PC{i+1}'] = self.eigenvectors[:, i]
            data[f'PC{i+1}_variance'] = self.explained_variance_ratio[i]
        
        return pd.DataFrame(data)


@dataclass
class FSTResult:
    """Results from FST calculation."""
    population_pairs: List[Tuple[str, str]]
    fst_values: np.ndarray
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        data = []
        for (pop1, pop2), fst in zip(self.population_pairs, self.fst_values):
            data.append({
                'population1': pop1,
                'population2': pop2,
                'fst': fst,
            })
        return pd.DataFrame(data)


class PopulationGenetics:
    """
    Population genetics analysis module.
    
    Implements:
    - PCA on genotype matrix
    - Pairwise FST calculation
    - IBS (Identity-by-State) distance matrix
    - Kinship matrix (VanRaden)
    """
    
    def __init__(self, vcf_path: Path, chunk_size: int = 10000):
        """
        Initialize population genetics module.
        
        Args:
            vcf_path: Path to VCF file
            chunk_size: Number of variants to process at once
        """
        self.vcf_path = Path(vcf_path)
        self.chunk_size = chunk_size
        
        if not self.vcf_path.exists():
            raise PopulationError(f"VCF file not found: {vcf_path}")
    
    def run_pca(self, n_components: int = 10) -> PCAResult:
        """
        Perform PCA on genotype matrix.
        
        Uses allele counts (0, 1, 2) as features, with imputation of
        missing values to the mean.
        
        Args:
            n_components: Number of principal components to compute
            
        Returns:
            PCAResult with eigenvalues, eigenvectors, and sample IDs
        """
        logger.info(f"Running PCA (n_components={n_components})")
        
        # Load genotypes
        genotypes, sample_ids = self._load_genotype_matrix()
        
        # Impute missing values to mean
        genotypes_imputed = self._impute_missing(genotypes)
        
        # Center and scale
        genotypes_centered = genotypes_imputed - np.mean(genotypes_imputed, axis=0)
        
        # Compute covariance matrix
        n_samples = genotypes_centered.shape[0]
        cov_matrix = np.dot(genotypes_centered, genotypes_centered.T) / n_samples
        
        # Eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        
        # Sort by descending eigenvalue
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Keep top components
        eigenvalues = eigenvalues[:n_components]
        eigenvectors = eigenvectors[:, :n_components]
        
        # Calculate explained variance
        total_var = np.sum(eigenvalues)
        explained_var = eigenvalues / total_var if total_var > 0 else np.zeros(n_components)
        
        logger.info(f"PCA complete. Top 3 PCs explain {explained_var[:3].sum():.1%} of variance")
        
        return PCAResult(
            eigenvalues=eigenvalues,
            eigenvectors=eigenvectors,
            explained_variance_ratio=explained_var,
            sample_ids=sample_ids,
        )
    
    def calculate_fst(
        self,
        sample_groups: Dict[str, List[str]],
    ) -> FSTResult:
        """
        Calculate pairwise FST between populations.
        
        Uses Weir-Cockerham FST estimator.
        
        Args:
            sample_groups: Dictionary mapping population names to sample IDs
            
        Returns:
            FSTResult with pairwise FST values
        """
        logger.info(f"Calculating FST for {len(sample_groups)} populations")
        
        # Load genotypes
        genotypes, sample_ids = self._load_genotype_matrix()
        
        # Create sample index mapping
        sample_idx = {s: i for i, s in enumerate(sample_ids)}
        
        # Calculate FST for each pair
        populations = list(sample_groups.keys())
        pairs = []
        fst_values = []
        
        for i, pop1 in enumerate(populations):
            for pop2 in populations[i+1:]:
                idx1 = [sample_idx[s] for s in sample_groups[pop1] if s in sample_idx]
                idx2 = [sample_idx[s] for s in sample_groups[pop2] if s in sample_idx]
                
                if not idx1 or not idx2:
                    logger.warning(f"Skipping FST for {pop1}-{pop2}: missing samples")
                    continue
                
                fst = self._weir_cockerham_fst(genotypes, idx1, idx2)
                
                pairs.append((pop1, pop2))
                fst_values.append(fst)
        
        logger.info(f"FST calculation complete for {len(pairs)} pairs")
        
        return FSTResult(
            population_pairs=pairs,
            fst_values=np.array(fst_values),
        )
    
    def calculate_ibs(self) -> Tuple[np.ndarray, List[str]]:
        """
        Calculate Identity-by-State (IBS) distance matrix.
        
        Returns:
            Tuple of (distance_matrix, sample_ids)
        """
        logger.info("Calculating IBS distance matrix")
        
        genotypes, sample_ids = self._load_genotype_matrix()
        n_samples = len(sample_ids)
        
        # Calculate pairwise IBS
        ibs_matrix = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(i, n_samples):
                # Calculate proportion of matching alleles
                matches = 0
                total = 0
                
                for k in range(genotypes.shape[1]):
                    gi = genotypes[i, k]
                    gj = genotypes[j, k]
                    
                    # Skip missing
                    if gi < 0 or gj < 0:
                        continue
                    
                    # Count matching alleles
                    if gi == gj:
                        matches += 2  # Both match
                    elif abs(gi - gj) == 1:
                        matches += 1  # One matches (heterozygote)
                    
                    total += 2
                
                if total > 0:
                    ibs = matches / total
                else:
                    ibs = 0
                
                ibs_matrix[i, j] = ibs
                ibs_matrix[j, i] = ibs
        
        # Convert to distance (1 - IBS)
        distance_matrix = 1 - ibs_matrix
        
        logger.info(f"IBS matrix complete: {n_samples}x{n_samples}")
        
        return distance_matrix, sample_ids
    
    def calculate_kinship(self) -> Tuple[np.ndarray, List[str]]:
        """
        Calculate VanRaden kinship matrix.
        
        Returns:
            Tuple of (kinship_matrix, sample_ids)
        """
        logger.info("Calculating kinship matrix (VanRaden method)")
        
        genotypes, sample_ids = self._load_genotype_matrix()
        
        # Impute missing to mean
        genotypes_imputed = self._impute_missing(genotypes)
        
        # Center genotypes
        genotypes_centered = genotypes_imputed - np.mean(genotypes_imputed, axis=0)
        
        # Normalize by sqrt of variance
        genotype_std = np.std(genotypes_centered, axis=0)
        genotype_std[genotype_std == 0] = 1  # Avoid division by zero
        genotypes_normalized = genotypes_centered / genotype_std
        
        # Calculate kinship
        n_snps = genotypes_normalized.shape[1]
        kinship = np.dot(genotypes_normalized, genotypes_normalized.T) / n_snps
        
        logger.info(f"Kinship matrix complete: {kinship.shape}")
        
        return kinship, sample_ids
    
    def detect_outliers(
        self,
        pca_result: PCAResult,
        n_sd: float = 3.0,
    ) -> List[str]:
        """
        Detect outlier samples based on PCA.
        
        Samples beyond n_sd standard deviations from mean on any PC are flagged.
        
        Args:
            pca_result: Result from run_pca()
            n_sd: Number of standard deviations threshold
            
        Returns:
            List of outlier sample IDs
        """
        outliers = []
        
        for pc_idx in range(min(3, len(pca_result.eigenvalues))):
            pc_values = pca_result.eigenvectors[:, pc_idx]
            mean = np.mean(pc_values)
            std = np.std(pc_values)
            
            if std == 0:
                continue
            
            for i, sample in enumerate(pca_result.sample_ids):
                z_score = abs(pc_values[i] - mean) / std
                if z_score > n_sd:
                    outliers.append(f"{sample} (PC{pc_idx+1}, z={z_score:.1f})")
        
        return outliers
    
    def _load_genotype_matrix(self) -> Tuple[np.ndarray, List[str]]:
        """
        Load genotype matrix from VCF.
        
        Returns:
            Tuple of (genotypes, sample_ids)
            genotypes shape: (n_samples, n_snps)
        """
        logger.debug("Loading genotype matrix")
        
        genotypes_list = []
        sample_ids = None
        
        with VCFReader(self.vcf_path) as reader:
            sample_ids = reader.sample_names
            
            for chunk in reader.iter_chunks(self.chunk_size):
                for variant in chunk:
                    # Skip non-biallelic
                    if not variant.is_biallelic:
                        continue
                    
                    # Extract genotypes (0, 1, 2, -1 for missing)
                    gts = []
                    for i in range(len(sample_ids)):
                        gt = variant.genotypes[i]
                        if np.any(gt < 0) or np.any(gt == 3):
                            gts.append(-1)  # Missing
                        else:
                            gts.append(int(np.sum(gt)))  # 0, 1, or 2
                    
                    genotypes_list.append(gts)
        
        genotypes = np.array(genotypes_list).T  # Shape: (n_samples, n_snps)
        
        logger.debug(f"Loaded genotype matrix: {genotypes.shape}")
        
        return genotypes, sample_ids
    
    def _impute_missing(self, genotypes: np.ndarray) -> np.ndarray:
        """Impute missing values to column means."""
        genotypes_imputed = genotypes.copy().astype(float)
        
        for j in range(genotypes.shape[1]):
            col = genotypes[:, j]
            mask = col >= 0
            
            if np.any(mask):
                mean_val = np.mean(col[mask])
                genotypes_imputed[~mask, j] = mean_val
            else:
                genotypes_imputed[:, j] = 0
        
        return genotypes_imputed
    
    def _weir_cockerham_fst(
        self,
        genotypes: np.ndarray,
        idx1: List[int],
        idx2: List[int],
    ) -> float:
        """
        Calculate Weir-Cockerham FST between two populations.
        
        Args:
            genotypes: Genotype matrix (n_samples x n_snps)
            idx1: Indices for population 1
            idx2: Indices for population 2
            
        Returns:
            FST value
        """
        # Extract population genotypes
        pop1 = genotypes[idx1, :]
        pop2 = genotypes[idx2, :]
        
        # Calculate allele frequencies
        # Using only non-missing genotypes
        fst_numerator = 0
        fst_denominator = 0
        
        for snp in range(genotypes.shape[1]):
            g1 = pop1[:, snp]
            g2 = pop2[:, snp]
            
            # Skip if all missing in either population
            if np.all(g1 < 0) or np.all(g2 < 0):
                continue
            
            # Filter missing
            g1_valid = g1[g1 >= 0]
            g2_valid = g2[g2 >= 0]
            
            if len(g1_valid) == 0 or len(g2_valid) == 0:
                continue
            
            # Calculate allele frequencies (p = freq of alt allele)
            n1 = len(g1_valid)
            n2 = len(g2_valid)
            
            p1 = np.sum(g1_valid) / (2 * n1)
            p2 = np.sum(g2_valid) / (2 * n2)
            
            # Overall frequency
            p = (np.sum(g1_valid) + np.sum(g2_valid)) / (2 * (n1 + n2))
            
            if p == 0 or p == 1:
                continue
            
            # Variance components
            # Between-population variance
            a = (p1 - p)**2 + (p2 - p)**2 - p * (1 - p) / n1 - p * (1 - p) / n2
            
            # Within-population variance
            b = p * (1 - p)
            
            fst_numerator += a
            fst_denominator += a + b
        
        if fst_denominator > 0:
            return fst_numerator / fst_denominator
        else:
            return 0.0
