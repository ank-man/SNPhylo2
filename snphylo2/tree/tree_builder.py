"""
Phylogenetic tree building with multiple engines.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass

from snphylo2.config import TreeConfig, TreeEngine, BootstrapMethod
from snphylo2.exceptions import TreeError, EngineError
from snphylo2.io.writers import FASTAWriter, PHYLIPWriter
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class TreeResult:
    """Result from tree building."""
    tree_file: Path
    model: Optional[str] = None
    log_likelihood: Optional[float] = None
    bootstrap_support: Optional[float] = None
    runtime_seconds: Optional[float] = None


class TreeBuilder:
    """
    Unified interface for multiple tree-building engines.
    """
    
    def __init__(self, config: TreeConfig):
        """
        Initialize tree builder.
        
        Args:
            config: Tree building configuration
        """
        self.config = config
        
        # Map engines to implementation methods
        self._builders = {
            TreeEngine.IQTREE2: self._build_iqtree2,
            TreeEngine.RAXML_NG: self._build_raxml_ng,
            TreeEngine.FASTTREE: self._build_fasttree,
            TreeEngine.PHYML: self._build_phyml,
            TreeEngine.DNAML: self._build_dnaml,
        }
    
    def build(self, input_path: Path, output_path: Path) -> Dict[str, Any]:
        """
        Build phylogenetic tree.
        
        Args:
            input_path: Input VCF or alignment
            output_path: Output tree file (Newick)
            
        Returns:
            Dictionary with tree building results
        """
        logger.info(f"Building tree with {self.config.engine.value}")
        
        builder = self._builders.get(self.config.engine)
        if not builder:
            raise TreeError(f"Unknown tree engine: {self.config.engine}")
        
        return builder(input_path, output_path)
    
    def _build_iqtree2(
        self,
        input_path: Path,
        output_path: Path,
    ) -> Dict[str, Any]:
        """
        Build tree using IQ-TREE2.
        
        IQ-TREE2 is the recommended engine for most analyses.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            
            # Convert input to alignment if needed
            alignment_path = self._prepare_alignment(input_path, tmp_path)
            
            # Build IQ-TREE2 command
            cmd = [
                "iqtree2",
                "-s", str(alignment_path),
                "-pre", str(tmp_path / "tree"),
                "-nt", str(self.config.threads or "AUTO"),
            ]
            
            # Model selection
            if self.config.model_selection:
                cmd.extend(["-m", "MFP"])  # ModelFinder Plus
            else:
                # Use specified model
                model = self.config.candidate_models[0] if self.config.candidate_models else "GTR+ASC"
                cmd.extend(["-m", model])
            
            # Bootstrap
            if self.config.bootstrap.method == BootstrapMethod.ULTRAFAST:
                cmd.extend(["-bb", str(self.config.bootstrap.replicates)])
            elif self.config.bootstrap.method == BootstrapMethod.STANDARD:
                cmd.extend(["-b", str(self.config.bootstrap.replicates)])
            elif self.config.bootstrap.method == BootstrapMethod.SHALRT:
                cmd.extend(["-alrt", str(self.config.bootstrap.replicates)])
            
            # Seed
            if self.config.bootstrap.seed:
                cmd.extend(["-seed", str(self.config.bootstrap.seed)])
            
            # Outgroup
            if self.config.outgroup:
                cmd.extend(["-o", self.config.outgroup])
            
            # SNP-specific: ascertainment bias correction
            if self.config.asc_bias_correction:
                # Check if model already has ASC
                if not any("ASC" in m for m in self.config.candidate_models):
                    cmd.extend(["-m", "GTR+ASC"])
            
            # Silent mode
            cmd.append("-quiet")
            
            try:
                logger.debug(f"Running: {' '.join(cmd)}")
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                
                # Parse log for best model and likelihood
                log_file = tmp_path / "tree.log"
                model = None
                log_likelihood = None
                
                if log_file.exists():
                    with open(log_file) as f:
                        for line in f:
                            if "Best-fit model" in line:
                                model = line.split(":")[-1].strip()
                            if "Log-likelihood" in line:
                                try:
                                    log_likelihood = float(line.split(":")[-1].strip())
                                except:
                                    pass
                
                # Copy tree file to output
                tree_file = tmp_path / "tree.treefile"
                if tree_file.exists():
                    import shutil
                    shutil.copy(tree_file, output_path)
                else:
                    raise TreeError("IQ-TREE2 did not produce tree file")
                
                # Calculate mean bootstrap support
                bootstrap_support = self._calculate_mean_support(output_path)
                
                return {
                    'tree_file': output_path,
                    'model': model,
                    'log_likelihood': log_likelihood,
                    'bootstrap_support': bootstrap_support,
                }
                
            except subprocess.CalledProcessError as e:
                raise EngineError(
                    "IQ-TREE2 failed",
                    command=" ".join(cmd),
                    returncode=e.returncode,
                    stderr=e.stderr,
                )
            except FileNotFoundError:
                raise TreeError("iqtree2 not found. Please install IQ-TREE2.")
    
    def _build_raxml_ng(
        self,
        input_path: Path,
        output_path: Path,
    ) -> Dict[str, Any]:
        """Build tree using RAxML-NG."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            
            alignment_path = self._prepare_alignment(input_path, tmp_path)
            
            cmd = [
                "raxml-ng",
                "--all",
                "--msa", str(alignment_path),
                "--model", self.config.candidate_models[0] if self.config.candidate_models else "GTR+G",
                "--tree", "rand{1}",
                "--threads", str(self.config.threads or 1),
                "--seed", str(self.config.bootstrap.seed or 42),
                "--bs-tree", str(self.config.bootstrap.replicates),
                "--prefix", str(tmp_path / "tree"),
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Copy result
                import shutil
                result_tree = tmp_path / "tree.raxml.bestTree"
                if result_tree.exists():
                    shutil.copy(result_tree, output_path)
                
                return {
                    'tree_file': output_path,
                    'model': self.config.candidate_models[0] if self.config.candidate_models else "GTR+G",
                }
                
            except subprocess.CalledProcessError as e:
                raise TreeError(f"RAxML-NG failed: {e.stderr.decode() if e.stderr else 'Unknown'}")
            except FileNotFoundError:
                raise TreeError("raxml-ng not found. Please install RAxML-NG.")
    
    def _build_fasttree(
        self,
        input_path: Path,
        output_path: Path,
    ) -> Dict[str, Any]:
        """
        Build tree using FastTree.
        
        FastTree is the fastest option for very large datasets
        but less accurate than IQ-TREE2.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            
            # FastTree requires FASTA format
            alignment_path = self._prepare_alignment(input_path, tmp_path, format="fasta")
            
            cmd = [
                "FastTree",
                "-gtr",  # Generalized time-reversible model
                "-gamma",  # Gamma rate heterogeneity
                "-nt",  # Nucleotide sequences
                str(alignment_path),
            ]
            
            try:
                with open(output_path, 'w') as f:
                    subprocess.run(cmd, check=True, stdout=f, capture_stderr=True)
                
                return {
                    'tree_file': output_path,
                    'model': 'GTR+G',
                }
                
            except subprocess.CalledProcessError as e:
                raise TreeError(f"FastTree failed: {e.stderr.decode() if e.stderr else 'Unknown'}")
            except FileNotFoundError:
                raise TreeError("FastTree not found. Please install FastTree.")
    
    def _build_phyml(
        self,
        input_path: Path,
        output_path: Path,
    ) -> Dict[str, Any]:
        """Build tree using PhyML."""
        raise NotImplementedError("PhyML builder not yet implemented")
    
    def _build_dnaml(
        self,
        input_path: Path,
        output_path: Path,
    ) -> Dict[str, Any]:
        """
        Build tree using PHYLIP dnaml.
        
        This is for backwards compatibility with original SNPhylo.
        """
        raise NotImplementedError("DNAML builder not yet implemented")
    
    def _prepare_alignment(
        self,
        input_path: Path,
        tmp_path: Path,
        format: str = "phylip",
    ) -> Path:
        """
        Prepare alignment file from input.
        
        Args:
            input_path: Input VCF or alignment
            tmp_path: Temporary directory
            format: Output format (phylip, fasta)
            
        Returns:
            Path to alignment file
        """
        # Check if input is already an alignment
        suffix = input_path.suffix.lower()
        
        if suffix in ['.phy', '.phylip', '.fasta', '.fa']:
            # Already alignment format
            return input_path
        
        # Convert VCF to alignment
        from snphylo2.io.vcf_reader import VCFReader
        from snphylo2.io.writers import FASTAWriter, PHYLIPWriter, Alignment
        
        # Read VCF and convert to alignment
        alignment_data = {
            'sample_names': [],
            'sequences': {},
            'variant_positions': [],
        }
        
        with VCFReader(input_path) as reader:
            alignment_data['sample_names'] = reader.sample_names
            
            # Initialize sequences
            for sample in reader.sample_names:
                alignment_data['sequences'][sample] = []
            
            # Collect variants
            for variant in reader:
                if not variant.is_snp or not variant.is_biallelic:
                    continue
                
                alignment_data['variant_positions'].append(variant.pos)
                
                for i, sample in enumerate(reader.sample_names):
                    gt = variant.genotypes[i]
                    
                    # Convert genotype to base
                    if np.any(gt < 0) or np.any(gt == 3):  # Missing
                        base = 'N'
                    elif np.all(gt == 0):  # Homozygous ref
                        base = variant.ref
                    elif np.all(gt == 2):  # Homozygous alt
                        base = variant.alt[0]
                    else:  # Heterozygous
                        base = self._get_ambiguity(variant.ref, variant.alt[0])
                    
                    alignment_data['sequences'][sample].append(base)
        
        # Join sequences
        for sample in alignment_data['sequences']:
            alignment_data['sequences'][sample] = ''.join(alignment_data['sequences'][sample])
        
        # Create alignment object
        alignment = Alignment(
            sample_names=alignment_data['sample_names'],
            sequences=alignment_data['sequences'],
            variant_positions=alignment_data['variant_positions'],
        )
        
        # Write in requested format
        if format == "fasta":
            output_file = tmp_path / "alignment.fa"
            writer = FASTAWriter(output_file)
            writer.write(alignment)
        else:
            output_file = tmp_path / "alignment.phy"
            writer = PHYLIPWriter(output_file)
            writer.write(alignment)
        
        return output_file
    
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
    
    def _calculate_mean_support(self, tree_path: Path) -> Optional[float]:
        """
        Calculate mean bootstrap support from tree.
        
        Args:
            tree_path: Path to Newick tree with support values
            
        Returns:
            Mean bootstrap support or None if no values
        """
        try:
            from ete3 import Tree
            
            tree = Tree(str(tree_path))
            supports = []
            
            for node in tree.traverse():
                if not node.is_leaf() and node.support is not None:
                    supports.append(node.support)
            
            if supports:
                return sum(supports) / len(supports)
            return None
            
        except Exception as e:
            logger.warning(f"Could not calculate mean support: {e}")
            return None
