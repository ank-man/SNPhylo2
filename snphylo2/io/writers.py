"""
Output format writers for alignments and trees.
"""

from pathlib import Path
from typing import List, Dict, Optional, Union
from dataclasses import dataclass
import numpy as np

from snphylo2.exceptions import InputError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class Alignment:
    """Represents a sequence alignment."""
    sample_names: List[str]
    sequences: Dict[str, str]  # sample_id -> sequence
    variant_positions: Optional[List[int]] = None
    
    def __len__(self) -> int:
        """Return alignment length."""
        if not self.sequences:
            return 0
        return len(next(iter(self.sequences.values())))
    
    def get_consensus(self, threshold: float = 0.5) -> str:
        """Get consensus sequence from alignment."""
        if not self.sequences:
            return ""
        
        length = len(self)
        consensus = []
        
        for i in range(length):
            bases = [seq[i] for seq in self.sequences.values() if seq[i] not in '-N?']
            if not bases:
                consensus.append('N')
                continue
            
            from collections import Counter
            counts = Counter(bases)
            most_common = counts.most_common(1)[0]
            
            if most_common[1] / len(bases) >= threshold:
                consensus.append(most_common[0])
            else:
                # Use IUPAC ambiguity code
                unique_bases = set(bases)
                consensus.append(self._get_ambiguity_code(unique_bases))
        
        return ''.join(consensus)
    
    def _get_ambiguity_code(self, bases: set) -> str:
        """Convert set of bases to IUPAC ambiguity code."""
        ambiguities = {
            frozenset(['A', 'G']): 'R',
            frozenset(['C', 'T']): 'Y',
            frozenset(['A', 'C']): 'M',
            frozenset(['G', 'T']): 'K',
            frozenset(['A', 'T']): 'W',
            frozenset(['C', 'G']): 'S',
            frozenset(['A', 'C', 'G']): 'V',
            frozenset(['A', 'C', 'T']): 'H',
            frozenset(['A', 'G', 'T']): 'D',
            frozenset(['C', 'G', 'T']): 'B',
            frozenset(['A', 'C', 'G', 'T']): 'N',
        }
        return ambiguities.get(frozenset(bases), 'N')


class FASTAWriter:
    """Write alignments in FASTA format."""
    
    def __init__(self, output_path: Path, line_length: int = 60):
        """
        Initialize FASTA writer.
        
        Args:
            output_path: Output file path
            line_length: Characters per line
        """
        self.output_path = Path(output_path)
        self.line_length = line_length
    
    def write(self, alignment: Alignment) -> None:
        """
        Write alignment to FASTA file.
        
        Args:
            alignment: Alignment object to write
        """
        try:
            with open(self.output_path, 'w') as f:
                for sample_id, sequence in alignment.sequences.items():
                    f.write(f">{sample_id}\n")
                    # Write sequence in lines
                    for i in range(0, len(sequence), self.line_length):
                        f.write(sequence[i:i + self.line_length] + '\n')
            
            logger.info(f"Wrote FASTA alignment to {self.output_path}")
            logger.debug(f"  Samples: {len(alignment.sample_names)}, Length: {len(alignment)}")
            
        except Exception as e:
            raise InputError(f"Failed to write FASTA file: {e}")
    
    def write_from_genotypes(
        self,
        sample_names: List[str],
        genotypes: np.ndarray,  # Shape: (n_samples, n_snps)
        ref_alleles: List[str],
        alt_alleles: List[str],
    ) -> None:
        """
        Convert genotype matrix to FASTA and write.
        
        Args:
            sample_names: List of sample IDs
            genotypes: Genotype matrix (0=ref, 1=het, 2=alt, -1=missing)
            ref_alleles: Reference alleles for each SNP
            alt_alleles: Alternate alleles for each SNP
        """
        alignment = self._genotypes_to_alignment(
            sample_names, genotypes, ref_alleles, alt_alleles
        )
        self.write(alignment)
    
    def _genotypes_to_alignment(
        self,
        sample_names: List[str],
        genotypes: np.ndarray,
        ref_alleles: List[str],
        alt_alleles: List[str],
    ) -> Alignment:
        """Convert genotypes to alignment with IUPAC codes."""
        sequences = {}
        
        for i, sample in enumerate(sample_names):
            seq = []
            for j, gt in enumerate(genotypes[i]):
                if gt == -1 or gt == 3:  # Missing
                    seq.append('N')
                elif gt == 0:  # Homozygous reference
                    seq.append(ref_alleles[j])
                elif gt == 2:  # Homozygous alternate
                    seq.append(alt_alleles[j])
                elif gt == 1:  # Heterozygous - use IUPAC ambiguity
                    seq.append(self._get_ambiguity(ref_alleles[j], alt_alleles[j]))
                else:
                    seq.append('N')
            
            sequences[sample] = ''.join(seq)
        
        return Alignment(sample_names=sample_names, sequences=sequences)
    
    def _get_ambiguity(self, ref: str, alt: str) -> str:
        """Get IUPAC ambiguity code for heterozygous genotype."""
        ambiguities = {
            ('A', 'G'): 'R', ('G', 'A'): 'R',
            ('C', 'T'): 'Y', ('T', 'C'): 'Y',
            ('A', 'C'): 'M', ('C', 'A'): 'M',
            ('G', 'T'): 'K', ('T', 'G'): 'K',
            ('A', 'T'): 'W', ('T', 'A'): 'W',
            ('C', 'G'): 'S', ('G', 'C'): 'S',
        }
        return ambiguities.get((ref, alt), 'N')


class PHYLIPWriter:
    """Write alignments in PHYLIP format."""
    
    def __init__(self, output_path: Path, interleaved: bool = False):
        """
        Initialize PHYLIP writer.
        
        Args:
            output_path: Output file path
            interleaved: Use interleaved format (default: sequential)
        """
        self.output_path = Path(output_path)
        self.interleaved = interleaved
    
    def write(self, alignment: Alignment) -> None:
        """
        Write alignment to PHYLIP file.
        
        Args:
            alignment: Alignment object
        """
        try:
            n_samples = len(alignment.sample_names)
            seq_length = len(alignment)
            
            with open(self.output_path, 'w') as f:
                # Header
                f.write(f"{n_samples} {seq_length}\n")
                
                if self.interleaved:
                    # Interleaved format
                    names = list(alignment.sequences.keys())
                    seqs = list(alignment.sequences.values())
                    max_name_len = max(len(n) for n in names)
                    
                    # Write first block with names
                    for name, seq in zip(names, seqs):
                        f.write(f"{name:{max_name_len + 2}}{seq[:60]}\n")
                    
                    # Write remaining blocks
                    for start in range(60, seq_length, 60):
                        f.write('\n')
                        for seq in seqs:
                            f.write(f"{' ' * (max_name_len + 2)}{seq[start:start+60]}\n")
                else:
                    # Sequential format
                    for name, seq in alignment.sequences.items():
                        # PHYLIP restricts names to 10 characters
                        short_name = name[:10]
                        f.write(f"{short_name:10} {seq}\n")
            
            logger.info(f"Wrote PHYLIP alignment to {self.output_path}")
            
        except Exception as e:
            raise InputError(f"Failed to write PHYLIP file: {e}")


class NewickWriter:
    """Write trees in Newick format."""
    
    def __init__(self, output_path: Path):
        """
        Initialize Newick writer.
        
        Args:
            output_path: Output file path
        """
        self.output_path = Path(output_path)
    
    def write(self, newick_string: str) -> None:
        """
        Write Newick tree string to file.
        
        Args:
            newick_string: Tree in Newick format
        """
        try:
            with open(self.output_path, 'w') as f:
                f.write(newick_string)
            
            logger.info(f"Wrote Newick tree to {self.output_path}")
            
        except Exception as e:
            raise InputError(f"Failed to write Newick file: {e}")
    
    def write_with_bootstrap(
        self,
        tree,  # ete3.Tree or Bio.Phylo tree
        output_path: Optional[Path] = None,
    ) -> None:
        """
        Write tree with bootstrap values.
        
        Args:
            tree: Tree object with support values
            output_path: Optional override for output path
        """
        path = output_path or self.output_path
        
        try:
            # Try ete3 first
            if hasattr(tree, 'write'):
                tree.write(outfile=str(path), format=1)
            # Fall back to Bio.Phylo
            elif hasattr(tree, 'root'):
                from Bio import Phylo
                Phylo.write(tree, path, 'newick')
            else:
                raise InputError("Unknown tree format")
            
            logger.info(f"Wrote tree to {path}")
            
        except Exception as e:
            raise InputError(f"Failed to write tree: {e}")


class NexusWriter:
    """Write alignments in NEXUS format."""
    
    def __init__(self, output_path: Path, data_type: str = "DNA"):
        """
        Initialize NEXUS writer.
        
        Args:
            output_path: Output file path
            data_type: Data type (DNA, RNA, Protein)
        """
        self.output_path = Path(output_path)
        self.data_type = data_type
    
    def write(self, alignment: Alignment) -> None:
        """
        Write alignment to NEXUS file.
        
        Args:
            alignment: Alignment object
        """
        try:
            n_samples = len(alignment.sample_names)
            seq_length = len(alignment)
            
            with open(self.output_path, 'w') as f:
                f.write("#NEXUS\n")
                f.write("BEGIN DATA;\n")
                f.write(f"  DIMENSIONS NTAX={n_samples} NCHAR={seq_length};\n")
                f.write(f"  FORMAT DATATYPE={self.data_type} MISSING=? GAP=-;\n")
                f.write("  MATRIX\n")
                
                for name, seq in alignment.sequences.items():
                    # Quote names if they contain special characters
                    if any(c in name for c in [' ', '-', ':', '(', ')', '[', ']']):
                        name = f"'{name}'"
                    f.write(f"    {name:20} {seq}\n")
                
                f.write("  ;\n")
                f.write("END;\n")
            
            logger.info(f"Wrote NEXUS alignment to {self.output_path}")
            
        except Exception as e:
            raise InputError(f"Failed to write NEXUS file: {e}")
