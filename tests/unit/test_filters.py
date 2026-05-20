"""
Unit tests for filtering module.
"""

import numpy as np
import pytest

from snphylo2.filtering.variant_filters import (
    MAFFilter,
    MissingnessFilter,
    DepthFilter,
    BiallelicFilter,
)
from snphylo2.io.vcf_reader import Variant


class TestMAFFilter:
    """Test MAF filter."""
    
    def test_maf_calculation(self):
        """Test MAF calculation on known genotype."""
        filter_obj = MAFFilter(min_maf=0.1, max_maf=0.5)
        
        # Create mock variant with MAF = 0.25 (25% alt alleles)
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T"],
            qual=30.0,
            filter=[],
            genotypes=np.array([[0], [0], [1], [2]]),  # 1/4 alt = 0.25 MAF
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == True  # 0.25 is within [0.1, 0.5]
    
    def test_maf_too_low(self):
        """Test filtering out low MAF variants."""
        filter_obj = MAFFilter(min_maf=0.1)
        
        # Create variant with MAF = 0.05 (too low)
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T"],
            qual=30.0,
            filter=[],
            genotypes=np.array([[0], [0], [0], [0], [0], [0], [0], [0], [0], [1]]),
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == False


class TestMissingnessFilter:
    """Test missingness filter."""
    
    def test_low_missingness_passes(self):
        """Test that low missingness passes filter."""
        filter_obj = MissingnessFilter(max_missing_rate=0.2)
        
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T"],
            qual=30.0,
            filter=[],
            genotypes=np.array([[0], [0], [0], [0], [-1]]),  # 20% missing
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == True
    
    def test_high_missingness_fails(self):
        """Test that high missingness fails filter."""
        filter_obj = MissingnessFilter(max_missing_rate=0.2)
        
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T"],
            qual=30.0,
            filter=[],
            genotypes=np.array([[-1], [-1], [0], [0], [0]]),  # 40% missing
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == False


class TestBiallelicFilter:
    """Test biallelic filter."""
    
    def test_biallelic_passes(self):
        """Test that biallelic SNP passes."""
        filter_obj = BiallelicFilter()
        
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T"],  # Single alt = biallelic
            qual=30.0,
            filter=[],
            genotypes=np.array([[0], [1], [2]]),
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == True
    
    def test_multiallelic_fails(self):
        """Test that multiallelic SNP is filtered."""
        filter_obj = BiallelicFilter()
        
        variant = Variant(
            chrom="chr1",
            pos=1000,
            id="rs1",
            ref="A",
            alt=["T", "G"],  # Two alts = multiallelic
            qual=30.0,
            filter=[],
            genotypes=np.array([[0], [1], [2]]),
            genotype_qualities=None,
            depths=None,
            allele_depths=None,
            info={},
        )
        
        assert filter_obj.apply(variant) == False
