"""
Visualization module for SNPhylo2.

Creates publication-quality plots using Plotly and optional R integration.
"""

from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass

import numpy as np

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from snphylo2.exceptions import VisualizationError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class PlotConfig:
    """Configuration for plot generation."""
    width: int = 800
    height: int = 600
    dpi: int = 300
    format: str = "html"  # html, png, pdf, svg
    title_font_size: int = 16
    axis_font_size: int = 12
    color_scheme: str = "default"


class Visualization:
    """
    Visualization module for phylogenetic and population genetics analyses.
    """
    
    def __init__(self, config: Optional[PlotConfig] = None):
        """
        Initialize visualization module.
        
        Args:
            config: Plot configuration
        """
        self.config = config or PlotConfig()
    
    def plot_pca(
        self,
        pca_result: Any,  # PCAResult
        metadata: Optional[Dict[str, Dict]] = None,
        color_by: Optional[str] = None,
        pc_x: int = 1,
        pc_y: int = 2,
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Create PCA scatter plot.
        
        Args:
            pca_result: PCAResult object
            metadata: Optional sample metadata dictionary
            color_by: Metadata field to color by
            pc_x: PC for x-axis (1-based)
            pc_y: PC for y-axis (1-based)
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available, skipping PCA plot")
            return None
        
        # Extract data
        sample_ids = pca_result.sample_ids
        x = pca_result.eigenvectors[:, pc_x - 1]
        y = pca_result.eigenvectors[:, pc_y - 1]
        
        # Get variance explained
        var_x = pca_result.explained_variance_ratio[pc_x - 1] * 100
        var_y = pca_result.explained_variance_ratio[pc_y - 1] * 100
        
        # Create hover text
        hover_text = [f"Sample: {s}<br>PC{pc_x}: {x[i]:.3f}<br>PC{pc_y}: {y[i]:.3f}" 
                      for i, s in enumerate(sample_ids)]
        
        # Determine colors
        if metadata and color_by and color_by in metadata.get(sample_ids[0], {}):
            colors = [metadata[s].get(color_by, "Unknown") for s in sample_ids]
            color_title = color_by
        else:
            colors = ["Sample"] * len(sample_ids)
            color_title = "Sample"
        
        # Create figure
        fig = px.scatter(
            x=x,
            y=y,
            color=colors,
            hover_name=sample_ids,
            hover_data={"PC" + str(pc_x): x, "PC" + str(pc_y): y},
            labels={
                'x': f'PC{pc_x} ({var_x:.1f}%)',
                'y': f'PC{pc_y} ({var_y:.1f}%)',
                'color': color_title,
            },
            title=f'PCA: PC{pc_x} vs PC{pc_y}',
            width=self.config.width,
            height=self.config.height,
        )
        
        fig.update_layout(
            title_font_size=self.config.title_font_size,
            xaxis_title_font_size=self.config.axis_font_size,
            yaxis_title_font_size=self.config.axis_font_size,
        )
        
        # Save if path provided
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def plot_fst_heatmap(
        self,
        fst_result: Any,  # FSTResult
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Create FST heatmap.
        
        Args:
            fst_result: FSTResult object
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        # Get unique populations
        populations = set()
        for pop1, pop2 in fst_result.population_pairs:
            populations.add(pop1)
            populations.add(pop2)
        populations = sorted(list(populations))
        
        # Create FST matrix
        n = len(populations)
        fst_matrix = np.zeros((n, n))
        
        for (pop1, pop2), fst in zip(fst_result.population_pairs, fst_result.fst_values):
            i = populations.index(pop1)
            j = populations.index(pop2)
            fst_matrix[i, j] = fst
            fst_matrix[j, i] = fst
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=fst_matrix,
            x=populations,
            y=populations,
            colorscale='YlOrRd',
            zmin=0,
            zmax=1,
            text=np.round(fst_matrix, 3),
            texttemplate='%{text}',
            textfont={"size": 10},
            hoverongaps=False,
        ))
        
        fig.update_layout(
            title='Pairwise FST Heatmap',
            xaxis_title='Population',
            yaxis_title='Population',
            width=self.config.width,
            height=self.config.height,
        )
        
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def plot_ibs_matrix(
        self,
        ibs_matrix: np.ndarray,
        sample_ids: List[str],
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Create IBS distance matrix heatmap.
        
        Args:
            ibs_matrix: IBS distance matrix
            sample_ids: Sample IDs
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        fig = go.Figure(data=go.Heatmap(
            z=ibs_matrix,
            x=sample_ids,
            y=sample_ids,
            colorscale='Viridis',
            zmin=0,
            zmax=1,
            hoverongaps=False,
        ))
        
        fig.update_layout(
            title='IBS Distance Matrix',
            xaxis_title='Sample',
            yaxis_title='Sample',
            width=self.config.width,
            height=self.config.height,
        )
        
        # Hide axis labels if many samples
        if len(sample_ids) > 20:
            fig.update_xaxes(showticklabels=False)
            fig.update_yaxes(showticklabels=False)
        
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def plot_missingness_heatmap(
        self,
        missingness_data: Dict[str, Dict[str, float]],
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Create missingness heatmap (samples x chromosomes).
        
        Args:
            missingness_data: Dict of sample -> chromosome -> missing_rate
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        # Convert to matrix
        samples = list(missingness_data.keys())
        chromosomes = set()
        for sample_data in missingness_data.values():
            chromosomes.update(sample_data.keys())
        chromosomes = sorted(list(chromosomes))
        
        matrix = np.zeros((len(samples), len(chromosomes)))
        for i, sample in enumerate(samples):
            for j, chrom in enumerate(chromosomes):
                matrix[i, j] = missingness_data[sample].get(chrom, 0)
        
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=chromosomes,
            y=samples,
            colorscale='Reds',
            zmin=0,
            zmax=1,
            hovertemplate='Sample: %{y}<br>Chrom: %{x}<br>Missing: %{z:.1%}',
        ))
        
        fig.update_layout(
            title='Missingness by Sample and Chromosome',
            xaxis_title='Chromosome',
            yaxis_title='Sample',
            width=self.config.width,
            height=max(400, len(samples) * 15),
        )
        
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def plot_maf_distribution(
        self,
        maf_values: List[float],
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Plot MAF distribution histogram.
        
        Args:
            maf_values: List of MAF values
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        fig = go.Figure(data=go.Histogram(
            x=maf_values,
            nbinsx=50,
            marker_color='steelblue',
            opacity=0.7,
        ))
        
        fig.update_layout(
            title='Minor Allele Frequency Distribution',
            xaxis_title='MAF',
            yaxis_title='Count',
            bargap=0.1,
            width=self.config.width,
            height=self.config.height,
        )
        
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def plot_depth_distribution(
        self,
        depth_values: List[int],
        output_path: Optional[Path] = None,
    ) -> Optional[go.Figure]:
        """
        Plot read depth distribution.
        
        Args:
            depth_values: List of depth values
            output_path: Optional output file path
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        # Filter extreme outliers for visualization
        max_depth = min(np.percentile(depth_values, 99), 100)
        filtered = [d for d in depth_values if d <= max_depth]
        
        fig = go.Figure(data=go.Histogram(
            x=filtered,
            nbinsx=50,
            marker_color='green',
            opacity=0.7,
        ))
        
        fig.update_layout(
            title='Read Depth Distribution',
            xaxis_title='Depth',
            yaxis_title='Count',
            bargap=0.1,
            width=self.config.width,
            height=self.config.height,
        )
        
        # Add mean line
        mean_depth = np.mean(depth_values)
        fig.add_vline(x=mean_depth, line_dash="dash", line_color="red",
                      annotation_text=f"Mean: {mean_depth:.1f}")
        
        if output_path:
            self._save_plot(fig, output_path)
        
        return fig
    
    def _save_plot(self, fig: go.Figure, output_path: Path) -> None:
        """
        Save plot to file.
        
        Args:
            fig: Plotly figure
            output_path: Output file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if output_path.suffix == '.html':
            fig.write_html(str(output_path))
        elif output_path.suffix in ['.png', '.jpg', '.jpeg', '.webp']:
            fig.write_image(str(output_path), scale=2)
        elif output_path.suffix == '.pdf':
            fig.write_image(str(output_path))
        elif output_path.suffix == '.svg':
            fig.write_image(str(output_path))
        else:
            # Default to HTML
            fig.write_html(str(output_path.with_suffix('.html')))
        
        logger.info(f"Saved plot: {output_path}")


class RVisualization:
    """
    R-based visualization using rpy2 (optional).
    
    Provides advanced tree visualization via ggtree when available.
    """
    
    def __init__(self):
        self.available = False
        self._check_rpy2()
    
    def _check_rpy2(self):
        """Check if rpy2 is available."""
        try:
            import rpy2.robjects as ro
            self.available = True
        except ImportError:
            pass
    
    def plot_tree_with_ggtree(
        self,
        tree_path: Path,
        metadata: Optional[Dict] = None,
        output_path: Optional[Path] = None,
    ) -> bool:
        """
        Create tree visualization using ggtree.
        
        Args:
            tree_path: Path to Newick tree file
            metadata: Optional sample metadata
            output_path: Output file path
            
        Returns:
            True if successful, False otherwise
        """
        if not self.available:
            logger.warning("rpy2 not available, skipping ggtree plot")
            return False
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.conversion import localconverter
            
            # Activate pandas conversion
            pandas2ri.activate()
            
            # R code for tree plotting
            r_code = f"""
            library(ggtree)
            library(ape)
            
            tree <- read.tree("{tree_path}")
            p <- ggtree(tree) + 
                 geom_tiplab(size=3) +
                 theme_tree2()
            
            ggsave("{output_path}", p, width=10, height=8)
            """
            
            ro.r(r_code)
            
            logger.info(f"Saved ggtree plot: {output_path}")
            return True
            
        except Exception as e:
            logger.warning(f"ggtree plotting failed: {e}")
            return False
