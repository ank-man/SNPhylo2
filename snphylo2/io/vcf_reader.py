"""
VCF/BCF file reading using cyvcf2 with streaming support.
"""

from pathlib import Path
from typing import Iterator, Optional, List, Dict, Any, Tuple, Union
from dataclasses import dataclass
import numpy as np

import cyvcf2
from snphylo2.exceptions import InputError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class Variant:
    """Represents a single variant with genotypes."""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: List[str]
    qual: Optional[float]
    filter: List[str]
    genotypes: np.ndarray  # Shape: (n_samples, ploidy)
    genotype_qualities: Optional[np.ndarray]
    depths: Optional[np.ndarray]
    allele_depths: Optional[np.ndarray]
    info: Dict[str, Any]
    
    @property
    def is_snp(self) -> bool:
        """Check if variant is a SNP."""
        return len(self.ref) == 1 and all(len(a) == 1 for a in self.alt)
    
    @property
    def is_biallelic(self) -> bool:
        """Check if variant is biallelic."""
        return len(self.alt) == 1
    
    @property
    def is_transition(self) -> bool:
        """Check if SNP is a transition (A<->G or C<->T)."""
        if not self.is_snp or not self.is_biallelic:
            return False
        alleles = {self.ref, self.alt[0]}
        return alleles in [{'A', 'G'}, {'C', 'T'}]
    
    @property
    def is_transversion(self) -> bool:
        """Check if SNP is a transversion."""
        if not self.is_snp or not self.is_biallelic:
            return False
        return not self.is_transition


class VCFReader:
    """
    Streaming VCF/BCF file reader using cyvcf2.
    
    Supports chunked reading for memory-efficient processing of large files.
    """
    
    def __init__(
        self,
        path: Path,
        region: Optional[str] = None,
        samples: Optional[List[str]] = None,
    ):
        """
        Initialize VCF reader.
        
        Args:
            path: Path to VCF/BCF file (must be indexed for random access)
            region: Optional genomic region (e.g., "chr1:1000-2000")
            samples: Optional list of samples to include
        """
        self.path = Path(path)
        self.region = region
        self.samples = samples
        
        if not self.path.exists():
            raise InputError(f"VCF file not found: {self.path}")
        
        try:
            self._vcf = cyvcf2.VCF(str(self.path), samples=samples)
        except Exception as e:
            raise InputError(f"Failed to open VCF file: {e}")
        
        # Store sample names
        self.sample_names = list(self._vcf.samples)
        self.n_samples = len(self.sample_names)
        
        # Store contig names if available
        self.contigs = list(self._vcf.seqlens.keys()) if hasattr(self._vcf, 'seqlens') else []
        
        logger.debug(f"Opened VCF with {self.n_samples} samples, {len(self.contigs)} contigs")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
    
    def close(self):
        """Close the VCF file."""
        if hasattr(self, '_vcf'):
            self._vcf.close()
    
    def __iter__(self) -> Iterator[Variant]:
        """Iterate over all variants."""
        try:
            for record in self._vcf:
                yield self._variant_from_record(record)
        except Exception as e:
            raise InputError(f"Error reading VCF: {e}")
    
    def iter_variants(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> Iterator[Variant]:
        """
        Iterate over variants in a specific region.
        
        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position
        """
        if chrom:
            region = f"{chrom}:{start or 1}-{end or ''}"
            try:
                for record in self._vcf(region):
                    yield self._variant_from_record(record)
            except Exception as e:
                raise InputError(f"Error reading region {region}: {e}")
        else:
            yield from self
    
    def iter_chunks(self, chunk_size: int = 10000) -> Iterator[List[Variant]]:
        """
        Iterate over variants in chunks for batch processing.
        
        Args:
            chunk_size: Number of variants per chunk
        """
        chunk = []
        for variant in self:
            chunk.append(variant)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk
    
    def get_variant_count(self) -> int:
        """
        Count total variants in file.
        Note: This iterates through the entire file.
        """
        count = 0
        for _ in self:
            count += 1
        # Reset iterator
        self._vcf = cyvcf2.VCF(str(self.path), samples=self.samples)
        return count
    
    def _variant_from_record(self, record: cyvcf2.Variant) -> Variant:
        """Convert cyvcf2 record to Variant dataclass."""
        # Extract genotypes
        genotypes = record.genotype.array()  # Shape: (n_samples, ploidy + 1)
        # Remove phasing info (last column)
        if genotypes.shape[1] > 1:
            genotypes = genotypes[:, :-1]
        
        # Extract genotype qualities if available
        gq = record.gt_quals if hasattr(record, 'gt_quals') else None
        
        # Extract depths if available
        dp = record.gt_depths if hasattr(record, 'gt_depths') else None
        
        # Extract allele depths if available
        ad = None
        try:
            ad = record.format('AD')
        except:
            pass
        
        # Extract INFO fields
        info = {}
        for key in record.INFO:
            try:
                info[key] = record.INFO[key]
            except:
                pass
        
        return Variant(
            chrom=record.CHROM,
            pos=record.POS,
            id=record.ID or f"{record.CHROM}:{record.POS}",
            ref=record.REF,
            alt=list(record.ALT) if record.ALT else [],
            qual=record.QUAL,
            filter=list(record.FILTERS) if record.FILTERS else [],
            genotypes=genotypes,
            genotype_qualities=gq,
            depths=dp,
            allele_depths=ad,
            info=info,
        )
    
    def get_header_info(self) -> Dict[str, Any]:
        """Get VCF header information."""
        return {
            'samples': self.sample_names,
            'contigs': self.contigs,
            'filters': list(self._vcf.filters),
        }


class BCFReader(VCFReader):
    """
    BCF file reader (binary VCF).
    Identical interface to VCFReader.
    """
    pass  # cyvcf2 handles BCF transparently


class VCFFilteredWriter:
    """
    Write filtered variants to a new VCF file.
    """
    
    def __init__(self, output_path: Path, template_reader: VCFReader):
        """
        Initialize writer with template from existing reader.
        
        Args:
            output_path: Output VCF file path
            template_reader: Reader to copy header from
        """
        self.output_path = output_path
        self.template_reader = template_reader
        self._writer = None
    
    def __enter__(self):
        # cyvcf2 doesn't have a direct writer, so we'll use pysam or write manually
        # For now, this is a placeholder for the interface
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._writer:
            self._writer.close()
        return False
    
    def write_variant(self, variant: Variant) -> None:
        """Write a variant to output."""
        # Implementation would use pysam.VariantFile or manual VCF writing
        pass


def create_index(vcf_path: Path) -> Path:
    """
    Create tabix index for a VCF.gz file.
    
    Args:
        vcf_path: Path to bgzipped VCF file
        
    Returns:
        Path to index file
    """
    import subprocess
    
    index_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    
    try:
        subprocess.run(
            ["tabix", "-p", "vcf", str(vcf_path)],
            check=True,
            capture_output=True,
        )
        logger.info(f"Created index: {index_path}")
        return index_path
    except subprocess.CalledProcessError as e:
        raise InputError(f"Failed to create index: {e.stderr.decode()}")
    except FileNotFoundError:
        raise InputError("tabix not found. Please install htslib.")
