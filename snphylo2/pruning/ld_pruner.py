"""
Linkage disequilibrium pruning implementations.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass

import numpy as np

from snphylo2.config import LDPruningConfig
from snphylo2.exceptions import PruningError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class PruningStats:
    """Statistics for LD pruning."""
    input_snps: int = 0
    output_snps: int = 0
    pruned_count: int = 0
    
    @property
    def retention_rate(self) -> float:
        if self.input_snps == 0:
            return 0.0
        return self.output_snps / self.input_snps


class LDPruner:
    """
    LD pruning using multiple backend methods.
    """
    
    def __init__(self, config: LDPruningConfig, threads: int = 1):
        """
        Initialize LD pruner.
        
        Args:
            config: LD pruning configuration
            threads: Number of threads
        """
        self.config = config
        self.threads = threads
    
    def run(self, input_path: Path, output_path: Path) -> PruningStats:
        """
        Run LD pruning.
        
        Args:
            input_path: Input VCF file
            output_path: Output pruned VCF
            
        Returns:
            Pruning statistics
        """
        logger.info(f"Starting LD pruning")
        logger.info(f"  Method: {self.config.method}")
        logger.info(f"  Window: {self.config.window_size}")
        logger.info(f"  Step: {self.config.step_size}")
        logger.info(f"  R² threshold: {self.config.r2_threshold}")
        
        # Use PLINK2 for efficient LD pruning
        stats = self._run_plink_pruning(input_path, output_path)
        
        logger.info(f"LD pruning complete")
        logger.info(f"  Input: {stats.input_snps:,} SNPs")
        logger.info(f"  Output: {stats.output_snps:,} SNPs")
        logger.info(f"  Pruned: {stats.pruned_count:,} SNPs")
        
        return stats
    
    def _run_plink_pruning(
        self,
        input_path: Path,
        output_path: Path,
    ) -> PruningStats:
        """
        Run LD pruning using PLINK2.
        
        PLINK2 is the most efficient method for large datasets.
        """
        stats = PruningStats()
        
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            
            # Convert VCF to PLINK format if needed
            prefix = tmp_path / "data"
            
            try:
                # Import VCF to PLINK
                cmd = [
                    "plink2",
                    "--vcf", str(input_path),
                    "--make-bed",
                    "--out", str(prefix),
                    "--allow-extra-chr",
                ]
                
                logger.debug(f"Running: {' '.join(cmd)}")
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Count input SNPs
                with open(f"{prefix}.bim") as f:
                    stats.input_snps = sum(1 for _ in f)
                
                # Run LD pruning
                prune_prefix = tmp_path / "pruned"
                cmd = [
                    "plink2",
                    "--bfile", str(prefix),
                    "--indep-pairwise",
                    str(self.config.window_size),
                    str(self.config.step_size),
                    str(self.config.r2_threshold),
                    "--out", str(prune_prefix),
                    "--allow-extra-chr",
                    "--threads", str(self.threads),
                ]
                
                logger.debug(f"Running: {' '.join(cmd)}")
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Extract pruned SNPs
                cmd = [
                    "plink2",
                    "--bfile", str(prefix),
                    "--extract", f"{prune_prefix}.prune.in",
                    "--recode", "vcf",
                    "--out", str(tmp_path / "pruned"),
                    "--allow-extra-chr",
                ]
                
                logger.debug(f"Running: {' '.join(cmd)}")
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Move output to final location
                import shutil
                shutil.move(f"{tmp_path}/pruned.vcf", output_path)
                
                # Count output SNPs
                with open(f"{prune_prefix}.prune.in") as f:
                    stats.output_snps = sum(1 for _ in f)
                
                stats.pruned_count = stats.input_snps - stats.output_snps
                
            except subprocess.CalledProcessError as e:
                raise PruningError(
                    f"PLINK2 pruning failed: {e.stderr.decode() if e.stderr else 'Unknown error'}"
                )
            except FileNotFoundError:
                raise PruningError("PLINK2 not found. Please install PLINK2.")
        
        return stats
    
    def _run_snprelate_pruning(
        self,
        input_path: Path,
        output_path: Path,
    ) -> PruningStats:
        """
        Run LD pruning using SNPRelate (R package).
        
        This is slower but doesn't require PLINK.
        """
        # This would call R script with SNPRelate
        # Implementation similar to original SNPhylo
        raise NotImplementedError("SNPRelate pruning not yet implemented")
    
    def get_representative_snps(
        self,
        input_path: Path,
        n_snps: int,
    ) -> List[str]:
        """
        Select a representative set of SNPs maximizing information.
        
        Uses a greedy algorithm to select SNPs that are not in LD
        while maximizing genome coverage.
        
        Args:
            input_path: Input VCF
            n_snps: Target number of SNPs
            
        Returns:
            List of selected SNP IDs
        """
        # This is a placeholder for a more sophisticated algorithm
        # that could be used for very large datasets
        logger.warning("Representative SNP selection not yet fully implemented")
        return []
