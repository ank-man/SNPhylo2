"""
Manuscript-quality LD decay and heatmap visualization.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from snphylo2.population.ld_decay import LDDecayResult
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


class LDPlotter:
    """
    Create publication-quality LD decay and heatmap plots.
    """
    
    def __init__(
        self,
        width: int = 800,
        height: int = 600,
        dpi: int = 300,
        style: str = "publication",  # "publication", "interactive"
    ):
        """
        Initialize LD plotter.
        
        Args:
            width: Figure width in pixels
            height: Figure height in pixels
            dpi: DPI for raster outputs
            style: Plot style
        """
        self.width = width
        self.height = height
        self.dpi = dpi
        self.style = style
        
        # Publication-style colors
        self.colors = [
            '#1f77b4',  # Blue
            '#ff7f0e',  # Orange
            '#2ca02c',  # Green
            '#d62728',  # Red
            '#9467bd',  # Purple
            '#8c564b',  # Brown
            '#e377c2',  # Pink
            '#7f7f7f',  # Gray
        ]
    
    def plot_ld_decay(
        self,
        results: Dict[str, LDDecayResult],
        output_path: Optional[Path] = None,
        title: str = "LD Decay",
        xlim: Optional[Tuple[float, float]] = None,
        ylim: Optional[Tuple[float, float]] = None,
        show_confidence: bool = True,
        manuscript_ready: bool = True,
    ) -> Optional[go.Figure]:
        """
        Create LD decay plot with multiple populations.
        
        Args:
            results: Dictionary of LDDecayResult by population
            output_path: Output file path
            title: Plot title
            xlim: X-axis limits (bp)
            ylim: Y-axis limits (r²)
            show_confidence: Show 5-95% confidence intervals
            manuscript_ready: Use manuscript-quality styling
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly not available")
            return None
        
        fig = go.Figure()
        
        # Add traces for each population
        for i, (pop_name, result) in enumerate(results.items()):
            color = self.colors[i % len(self.colors)]
            
            # Main line
            fig.add_trace(go.Scatter(
                x=result.distances / 1000,  # Convert to kb
                y=result.mean_r2,
                mode='lines',
                name=pop_name,
                line=dict(color=color, width=2.5 if manuscript_ready else 2),
                hovertemplate=f'<b>{pop_name}</b><br>' +
                            'Distance: %{x:.1f} kb<br>' +
                            'Mean r²: %{y:.3f}<br>' +
                            'N pairs: %{{customdata:,}}<extra></extra>',
                customdata=result.n_pairs.reshape(-1, 1),
            ))
            
            # Confidence interval
            if show_confidence:
                fig.add_trace(go.Scatter(
                    x=np.concatenate([result.distances / 1000, 
                                     result.distances[::-1] / 1000]),
                    y=np.concatenate([result.percentiles_95, 
                                     result.percentiles_5[::-1]]),
                    fill='toself',
                    fillcolor=f'rgba{tuple(list(int(color.lstrip("#")[i:i+2], 16) for i in (0, 2, 4)) + [0.2])}',
                    line=dict(width=0),
                    showlegend=False,
                    name=f'{pop_name} (5-95%)',
                    hoverinfo='skip',
                ))
        
        # Add half-decay reference line
        fig.add_hline(
            y=0.5,
            line_dash="dash",
            line_color="gray",
            annotation_text="r² = 0.5",
            annotation_position="right",
        )
        
        # Layout
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=16 if manuscript_ready else 14),
                x=0.5,
            ) if title else None,
            xaxis=dict(
                title='Distance (kb)',
                titlefont=dict(size=14),
                tickfont=dict(size=12),
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                range=[0, xlim[1]/1000] if xlim else None,
            ),
            yaxis=dict(
                title='Linkage Disequilibrium (r²)',
                titlefont=dict(size=14),
                tickfont=dict(size=12),
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                range=ylim if ylim else [0, 1],
            ),
            legend=dict(
                x=0.98,
                y=0.98,
                xanchor='right',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='gray',
                borderwidth=1,
            ),
            plot_bgcolor='white',
            width=self.width,
            height=self.height,
            margin=dict(l=70, r=40, t=60, b=60),
        )
        
        # Add annotation with half-decay distances
        if manuscript_ready:
            annotation_text = "Half-decay distances:<br>"
            for pop_name, result in results.items():
                if result.half_decay_distance:
                    annotation_text += f"{pop_name}: {result.half_decay_distance/1000:.1f} kb<br>"
            
            fig.add_annotation(
                x=0.02,
                y=0.98,
                xref='paper',
                yref='paper',
                text=annotation_text,
                showarrow=False,
                font=dict(size=10),
                align='left',
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='gray',
                borderwidth=1,
            )
        
        # Save
        if output_path:
            self._save_figure(fig, output_path)
        
        return fig
    
    def plot_ld_heatmap(
        self,
        r2_matrix: np.ndarray,
        positions: List[int],
        output_path: Optional[Path] = None,
        title: str = "LD Heatmap",
        chromosome: Optional[str] = None,
        highlight_snps: Optional[List[int]] = None,
    ) -> Optional[go.Figure]:
        """
        Create LD heatmap for a genomic region.
        
        Args:
            r2_matrix: r² matrix (n_snps x n_snps)
            positions: SNP positions
            output_path: Output file path
            title: Plot title
            chromosome: Chromosome name
            highlight_snps: List of SNP indices to highlight
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        # Convert positions to Mb
        positions_mb = np.array(positions) / 1_000_000
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=r2_matrix,
            x=positions_mb,
            y=positions_mb,
            colorscale=[
                [0, 'white'],
                [0.2, 'lightyellow'],
                [0.4, 'orange'],
                [0.6, 'red'],
                [0.8, 'darkred'],
                [1, 'black'],
            ],
            zmin=0,
            zmax=1,
            colorbar=dict(
                title='r²',
                titleside='right',
            ),
            hovertemplate='Pos X: %{x:.3f} Mb<br>' +
                         'Pos Y: %{y:.3f} Mb<br>' +
                         'r²: %{z:.3f}<extra></extra>',
        ))
        
        # Add diagonal line
        fig.add_trace(go.Scatter(
            x=[positions_mb[0], positions_mb[-1]],
            y=[positions_mb[0], positions_mb[-1]],
            mode='lines',
            line=dict(color='gray', width=1, dash='dash'),
            showlegend=False,
            hoverinfo='skip',
        ))
        
        # Highlight specific SNPs
        if highlight_snps:
            for idx in highlight_snps:
                if 0 <= idx < len(positions_mb):
                    fig.add_vline(x=positions_mb[idx], line_color='blue', line_width=1)
                    fig.add_hline(y=positions_mb[idx], line_color='blue', line_width=1)
        
        # Layout
        chrom_text = f" ({chromosome})" if chromosome else ""
        fig.update_layout(
            title=dict(
                text=f"{title}{chrom_text}",
                x=0.5,
            ),
            xaxis=dict(
                title='Position (Mb)',
                scaleanchor='y',
                scaleratio=1,
            ),
            yaxis=dict(
                title='Position (Mb)',
            ),
            width=self.height,  # Square aspect
            height=self.height,
        )
        
        if output_path:
            self._save_figure(fig, output_path)
        
        return fig
    
    def plot_ld_decay_comparison(
        self,
        results_list: List[Dict[str, LDDecayResult]],
        labels: List[str],
        output_path: Optional[Path] = None,
        title: str = "LD Decay Comparison",
    ) -> Optional[go.Figure]:
        """
        Compare LD decay across multiple datasets/experiments.
        
        Args:
            results_list: List of result dictionaries (one per dataset)
            labels: Labels for each dataset
            output_path: Output file path
            title: Plot title
            
        Returns:
            Plotly figure if available
        """
        if not PLOTLY_AVAILABLE:
            return None
        
        n_plots = len(results_list)
        fig = make_subplots(
            rows=1,
            cols=n_plots,
            subplot_titles=labels,
            shared_yaxes=True,
            horizontal_spacing=0.05,
        )
        
        for col_idx, (results, label) in enumerate(zip(results_list, labels), 1):
            for pop_idx, (pop_name, result) in enumerate(results.items()):
                color = self.colors[pop_idx % len(self.colors)]
                
                fig.add_trace(
                    go.Scatter(
                        x=result.distances / 1000,
                        y=result.mean_r2,
                        mode='lines',
                        name=pop_name,
                        line=dict(color=color, width=2),
                        showlegend=(col_idx == 1),
                        legendgroup=pop_name,
                    ),
                    row=1,
                    col=col_idx,
                )
            
            # Add reference line
            fig.add_hline(
                y=0.5,
                line_dash="dash",
                line_color="gray",
                row=1,
                col=col_idx,
            )
        
        fig.update_layout(
            title=title,
            height=self.height,
            width=self.width * n_plots,
        )
        
        fig.update_xaxes(title_text='Distance (kb)')
        fig.update_yaxes(title_text='r²', col=1)
        
        if output_path:
            self._save_figure(fig, output_path)
        
        return fig
    
    def create_manuscript_figure(
        self,
        ld_decay_results: Dict[str, LDDecayResult],
        output_path: Path,
        include_heatmap: bool = True,
        heatmap_data: Optional[Tuple[np.ndarray, List[int]]] = None,
    ) -> None:
        """
        Create a comprehensive manuscript figure combining multiple LD analyses.
        
        Args:
            ld_decay_results: LD decay results
            output_path: Output file path
            include_heatmap: Include LD heatmap subplot
            heatmap_data: Optional (r2_matrix, positions) for heatmap
        """
        if not PLOTLY_AVAILABLE:
            return
        
        if include_heatmap and heatmap_data:
            # Combined figure with decay and heatmap
            fig = make_subplots(
                rows=1,
                cols=2,
                subplot_titles=('A) LD Decay', 'B) LD Heatmap'),
                column_widths=[0.6, 0.4],
            )
            
            # Add decay plot to first subplot
            for pop_idx, (pop_name, result) in enumerate(ld_decay_results.items()):
                color = self.colors[pop_idx % len(self.colors)]
                
                fig.add_trace(
                    go.Scatter(
                        x=result.distances / 1000,
                        y=result.mean_r2,
                        mode='lines',
                        name=pop_name,
                        line=dict(color=color, width=2.5),
                    ),
                    row=1,
                    col=1,
                )
            
            fig.add_hline(y=0.5, line_dash="dash", line_color="gray", row=1, col=1)
            
            # Add heatmap to second subplot
            r2_matrix, positions = heatmap_data
            positions_mb = np.array(positions) / 1_000_000
            
            fig.add_trace(
                go.Heatmap(
                    z=r2_matrix,
                    x=positions_mb,
                    y=positions_mb,
                    colorscale='YlOrRd',
                    zmin=0,
                    zmax=1,
                    colorbar=dict(title='r²', x=0.97),
                ),
                row=1,
                col=2,
            )
            
            fig.update_xaxes(title_text='Distance (kb)', row=1, col=1)
            fig.update_yaxes(title_text='r²', row=1, col=1)
            fig.update_xaxes(title_text='Position (Mb)', row=1, col=2)
            fig.update_yaxes(title_text='Position (Mb)', row=1, col=2)
            
            fig.update_layout(
                height=self.height,
                width=int(self.width * 1.5),
                title_text="Linkage Disequilibrium Analysis",
            )
            
        else:
            # Just decay plot
            fig = self.plot_ld_decay(ld_decay_results)
        
        if fig and output_path:
            self._save_figure(fig, output_path)
    
    def _save_figure(self, fig: go.Figure, output_path: Path) -> None:
        """
        Save figure to file.
        
        Args:
            fig: Plotly figure
            output_path: Output file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if output_path.suffix == '.html':
            fig.write_html(str(output_path))
        elif output_path.suffix in ['.png', '.jpg']:
            fig.write_image(str(output_path), scale=3, width=self.width, height=self.height)
        elif output_path.suffix == '.pdf':
            fig.write_image(str(output_path), width=self.width, height=self.height)
        elif output_path.suffix == '.svg':
            fig.write_image(str(output_path), width=self.width, height=self.height)
        elif output_path.suffix in ['.eps', '.tiff', '.tif']:
            # For publication formats, use matplotlib
            self._save_with_matplotlib(fig, output_path)
        else:
            # Default to high-res PNG
            output_path = output_path.with_suffix('.png')
            fig.write_image(str(output_path), scale=3, width=self.width, height=self.height)
        
        logger.info(f"Saved LD plot: {output_path}")
    
    def _save_with_matplotlib(self, fig: go.Figure, output_path: Path) -> None:
        """Convert Plotly figure to matplotlib for publication formats."""
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("Matplotlib not available for publication format")
            return
        
        # This is a simplified conversion
        # Full implementation would extract data and recreate
        logger.warning(f"Publication format {output_path.suffix} requires manual conversion")
