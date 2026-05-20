#!/usr/bin/env python3
"""
Generate publication-quality figures for SNPhylo2 manuscripts.

Creates a comprehensive figure set demonstrating:
1. Phylogenetic trees with metadata
2. PCA plots with population structure
3. LD decay curves (PopLDdecay-style)
4. LD heatmaps
5. FST heatmaps
6. Combined multi-panel figures
"""

import argparse
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Add snphylo2 to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from snphylo2.population.ld_decay import LDDecayAnalyzer, LDDecayResult
from snphylo2.population.popgen import PopulationGenetics, PCAResult
from snphylo2.visualization.plots import Visualization, PlotConfig
from snphylo2.visualization.ld_plots import LDPlotter
from snphylo2.io.vcf_reader import VCFReader


class ManuscriptFigureGenerator:
    """Generate publication-ready figures for manuscripts."""
    
    def __init__(self, output_dir: str, figure_format: str = "pdf"):
        """
        Initialize figure generator.
        
        Args:
            output_dir: Directory for output figures
            figure_format: Output format (pdf, png, svg, tiff)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.format = figure_format
        
        # High-resolution settings for publications
        self.plot_config = PlotConfig(
            width=1200,
            height=900,
            dpi=300,
        )
        
        self.viz = Visualization(self.plot_config)
        self.ld_plotter = LDPlotter(width=1200, height=900, dpi=300)
    
    def generate_all_figures(self, vcf_path: Path, metadata_path: Optional[Path] = None):
        """
        Generate complete figure set from a VCF file.
        
        Args:
            vcf_path: Path to input VCF
            metadata_path: Optional metadata file
        """
        print("=" * 60)
        print("Generating Manuscript Figures")
        print("=" * 60)
        
        # Load metadata if provided
        metadata = None
        if metadata_path and metadata_path.exists():
            metadata = self._load_metadata(metadata_path)
            print(f"Loaded metadata for {len(metadata)} samples")
        
        # Figure 1: Phylogenetic analysis workflow
        print("\n[1/6] Generating workflow diagram...")
        self._create_workflow_diagram()
        
        # Figure 2: LD Decay curves
        print("\n[2/6] Analyzing LD decay...")
        self._create_ld_decay_figure(vcf_path, metadata)
        
        # Figure 3: PCA with population structure
        print("\n[3/6] Performing PCA...")
        self._create_pca_figure(vcf_path, metadata)
        
        # Figure 4: LD Heatmap for representative region
        print("\n[4/6] Generating LD heatmap...")
        self._create_ld_heatmap_figure(vcf_path)
        
        # Figure 5: Population differentiation (FST)
        print("\n[5/6] Calculating FST...")
        self._create_fst_figure(vcf_path, metadata)
        
        # Figure 6: Combined overview
        print("\n[6/6] Creating combined figure...")
        self._create_combined_figure(vcf_path, metadata)
        
        print("\n" + "=" * 60)
        print("All figures generated successfully!")
        print(f"Output directory: {self.output_dir}")
        print("=" * 60)
    
    def _load_metadata(self, metadata_path: Path) -> Dict[str, Dict]:
        """Load sample metadata from TSV file."""
        df = pd.read_csv(metadata_path, sep='\t')
        
        # Assume first column is sample_id
        sample_col = df.columns[0]
        
        metadata = {}
        for _, row in df.iterrows():
            sample_id = str(row[sample_col])
            metadata[sample_id] = row.to_dict()
        
        return metadata
    
    def _create_workflow_diagram(self):
        """Create SNPhylo2 workflow diagram."""
        # This would create a schematic diagram of the pipeline
        # For now, we'll create a text-based placeholder
        
        output_path = self.output_dir / f"figure_workflow.{self.format}"
        
        # Create a simple infographic-style figure using matplotlib
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
            from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
            
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.set_xlim(0, 10)
            ax.set_ylim(0, 10)
            ax.axis('off')
            
            # Title
            ax.text(5, 9.5, 'SNPhylo2 Analysis Workflow', 
                   ha='center', va='center', fontsize=16, fontweight='bold')
            
            # Steps
            steps = [
                (2, 8, 'VCF Input', '#3498db'),
                (2, 6.5, 'Quality Control', '#2ecc71'),
                (2, 5, 'Filter Variants', '#f39c12'),
                (2, 3.5, 'LD Pruning', '#e74c3c'),
                (2, 2, 'Tree Building', '#9b59b6'),
                (2, 0.5, 'Report', '#34495e'),
            ]
            
            for x, y, text, color in steps:
                box = FancyBboxPatch((x-1, y-0.4), 2, 0.8,
                                    boxstyle="round,pad=0.1",
                                    facecolor=color, edgecolor='black', linewidth=1.5)
                ax.add_patch(box)
                ax.text(x, y, text, ha='center', va='center', 
                       fontsize=11, color='white', fontweight='bold')
            
            # Arrows
            for i in range(len(steps) - 1):
                arrow = FancyArrowPatch((2, steps[i][1] - 0.4), (2, steps[i+1][1] + 0.4),
                                       arrowstyle='->', mutation_scale=20, 
                                       linewidth=2, color='black')
                ax.add_patch(arrow)
            
            # Side outputs
            outputs = [
                (5, 6.5, 'QC Report', '#2ecc71'),
                (5, 5, 'Filtered VCF', '#f39c12'),
                (5, 3.5, 'Pruned SNPs', '#e74c3c'),
                (5, 2, 'Phylogenetic Tree\n(HTML/PDF)', '#9b59b6'),
                (5, 0.5, 'Final Report', '#34495e'),
            ]
            
            for x, y, text, color in outputs:
                box = FancyBboxPatch((x-1, y-0.4), 2, 0.8,
                                    boxstyle="round,pad=0.1",
                                    facecolor=color, edgecolor='black', 
                                    linewidth=1.5, alpha=0.6)
                ax.add_patch(box)
                ax.text(x, y, text, ha='center', va='center', 
                       fontsize=9, color='black')
                
                # Arrow from step to output
                arrow = FancyArrowPatch((3, y), (4, y),
                                       arrowstyle='->', mutation_scale=15, 
                                       linewidth=1, color='gray', linestyle='--')
                ax.add_patch(arrow)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Saved: {output_path}")
            
        except Exception as e:
            print(f"  Could not create workflow diagram: {e}")
    
    def _create_ld_decay_figure(self, vcf_path: Path, metadata: Optional[Dict]):
        """Create LD decay figure with population comparisons."""
        print("  Analyzing LD decay patterns...")
        
        # Group samples by population if metadata available
        sample_groups = None
        if metadata:
            sample_groups = self._group_by_population(metadata)
        
        # Run LD decay analysis
        analyzer = LDDecayAnalyzer(
            max_distance=200_000,  # 200 kb
            bin_size=5_000,        # 5 kb bins
        )
        
        results = analyzer.analyze(vcf_path, sample_groups=sample_groups)
        
        # Create plot
        output_path = self.output_dir / f"figure_ld_decay.{self.format}"
        
        self.ld_plotter.plot_ld_decay(
            results,
            output_path=output_path,
            title="Linkage Disequilibrium Decay",
            show_confidence=True,
            manuscript_ready=True,
        )
        
        # Save data table
        for pop_name, result in results.items():
            data_path = self.output_dir / f"ld_decay_{pop_name}.tsv"
            result.to_dataframe().to_csv(data_path, sep='\t', index=False)
        
        print(f"  Saved: {output_path}")
    
    def _create_pca_figure(self, vcf_path: Path, metadata: Optional[Dict]):
        """Create PCA figure with population coloring."""
        print("  Performing PCA analysis...")
        
        popgen = PopulationGenetics(vcf_path)
        pca_result = popgen.run_pca(n_components=10)
        
        # Create multiple PCA plots
        plots = [
            (1, 2, "PC1 vs PC2"),
            (1, 3, "PC1 vs PC3"),
            (2, 3, "PC2 vs PC3"),
        ]
        
        for pc_x, pc_y, title in plots:
            output_path = self.output_dir / f"figure_pca_pc{pc_x}_pc{pc_y}.{self.format}"
            
            color_by = None
            if metadata:
                # Try common population columns
                for col in ['population', 'subspecies', 'ecotype', 'country']:
                    if any(col in m for m in metadata.values()):
                        color_by = col
                        break
            
            self.viz.plot_pca(
                pca_result,
                metadata=metadata,
                color_by=color_by,
                pc_x=pc_x,
                pc_y=pc_y,
                output_path=output_path,
            )
            
            print(f"  Saved: {output_path}")
        
        # Save eigenvalues
        eigen_path = self.output_dir / "pca_eigenvalues.tsv"
        eigen_df = pd.DataFrame({
            'PC': range(1, len(pca_result.eigenvalues) + 1),
            'eigenvalue': pca_result.eigenvalues,
            'explained_variance_ratio': pca_result.explained_variance_ratio,
            'cumulative_variance': np.cumsum(pca_result.explained_variance_ratio),
        })
        eigen_df.to_csv(eigen_path, sep='\t', index=False)
    
    def _create_ld_heatmap_figure(self, vcf_path: Path):
        """Create LD heatmap for a representative region."""
        print("  Generating LD heatmap...")
        
        from snphylo2.population.ld_decay import LDHeatmapGenerator
        
        # Use first 100 kb of first chromosome
        chrom = None
        with VCFReader(vcf_path) as reader:
            if reader.contigs:
                chrom = reader.contigs[0]
        
        if not chrom:
            print("  Could not determine chromosome for heatmap")
            return
        
        generator = LDHeatmapGenerator(window_size=100_000)
        
        try:
            r2_matrix, positions = generator.generate_heatmap(
                vcf_path,
                chrom=chrom,
                start=1,
                end=100_000,
            )
            
            output_path = self.output_dir / f"figure_ld_heatmap.{self.format}"
            
            self.ld_plotter.plot_ld_heatmap(
                r2_matrix,
                positions,
                output_path=output_path,
                title=f"LD Heatmap: {chrom}:1-100kb",
                chromosome=chrom,
            )
            
            print(f"  Saved: {output_path}")
            
        except Exception as e:
            print(f"  Could not generate heatmap: {e}")
    
    def _create_fst_figure(self, vcf_path: Path, metadata: Optional[Dict]):
        """Create FST heatmap figure."""
        if not metadata:
            print("  Skipping FST (no metadata)")
            return
        
        print("  Calculating pairwise FST...")
        
        # Get population groups
        sample_groups = self._group_by_population(metadata)
        
        if len(sample_groups) < 2:
            print("  Skipping FST (need 2+ populations)")
            return
        
        popgen = PopulationGenetics(vcf_path)
        fst_result = popgen.calculate_fst(sample_groups)
        
        output_path = self.output_dir / f"figure_fst_heatmap.{self.format}"
        
        self.viz.plot_fst_heatmap(
            fst_result,
            output_path=output_path,
        )
        
        # Save FST table
        fst_df = fst_result.to_dataframe()
        fst_path = self.output_dir / "fst_values.tsv"
        fst_df.to_csv(fst_path, sep='\t', index=False)
        
        print(f"  Saved: {output_path}")
    
    def _create_combined_figure(self, vcf_path: Path, metadata: Optional[Dict]):
        """Create multi-panel combined figure."""
        print("  Creating combined figure...")
        
        output_path = self.output_dir / f"figure_combined_overview.{self.format}"
        
        # This would combine multiple analyses into one figure
        # For now, create a placeholder
        
        try:
            import matplotlib.pyplot as plt
            from matplotlib.gridspec import GridSpec
            
            fig = plt.figure(figsize=(16, 12))
            gs = GridSpec(2, 3, figure=fig)
            
            # Panel A: LD Decay
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.text(0.5, 0.5, 'A) LD Decay', ha='center', va='center', fontsize=14)
            ax1.set_xlabel('Distance (kb)')
            ax1.set_ylabel('r²')
            ax1.axis('off')
            
            # Panel B: PCA
            ax2 = fig.add_subplot(gs[0, 1])
            ax2.text(0.5, 0.5, 'B) PCA', ha='center', va='center', fontsize=14)
            ax2.set_xlabel('PC1')
            ax2.set_ylabel('PC2')
            ax2.axis('off')
            
            # Panel C: FST
            ax3 = fig.add_subplot(gs[0, 2])
            ax3.text(0.5, 0.5, 'C) FST Heatmap', ha='center', va='center', fontsize=14)
            ax3.axis('off')
            
            # Panel D: LD Heatmap
            ax4 = fig.add_subplot(gs[1, :2])
            ax4.text(0.5, 0.5, 'D) LD Heatmap', ha='center', va='center', fontsize=14)
            ax4.axis('off')
            
            # Panel E: Statistics
            ax5 = fig.add_subplot(gs[1, 2])
            ax5.text(0.5, 0.5, 'E) Summary Statistics', ha='center', va='center', fontsize=14)
            ax5.axis('off')
            
            plt.suptitle('SNPhylo2 Population Genomics Analysis', fontsize=16, fontweight='bold')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Saved: {output_path}")
            
        except Exception as e:
            print(f"  Could not create combined figure: {e}")
    
    def _group_by_population(self, metadata: Dict) -> Dict[str, List[str]]:
        """Group samples by population from metadata."""
        groups = {}
        
        for sample_id, info in metadata.items():
            # Try to find population field
            pop = None
            for key in ['population', 'subspecies', 'ecotype', 'Pop', 'pop']:
                if key in info:
                    pop = info[key]
                    break
            
            if not pop:
                pop = 'unknown'
            
            if pop not in groups:
                groups[pop] = []
            groups[pop].append(sample_id)
        
        return groups


def main():
    parser = argparse.ArgumentParser(
        description="Generate manuscript figures for SNPhylo2"
    )
    parser.add_argument(
        "vcf",
        help="Input VCF file"
    )
    parser.add_argument(
        "--metadata",
        help="Sample metadata TSV file (optional)"
    )
    parser.add_argument(
        "-o", "--output",
        default="manuscript_figures",
        help="Output directory (default: manuscript_figures)"
    )
    parser.add_argument(
        "--format",
        choices=["pdf", "png", "svg", "tiff"],
        default="pdf",
        help="Output format (default: pdf)"
    )
    
    args = parser.parse_args()
    
    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        print(f"Error: VCF file not found: {vcf_path}")
        return 1
    
    metadata_path = Path(args.metadata) if args.metadata else None
    
    generator = ManuscriptFigureGenerator(
        args.output,
        figure_format=args.format,
    )
    
    generator.generate_all_figures(vcf_path, metadata_path)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
