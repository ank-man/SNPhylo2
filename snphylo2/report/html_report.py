"""
HTML report generation for SNPhylo2 results.
"""

from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime

from snphylo2.exceptions import ReportError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


class HTMLReport:
    """
    Generate HTML reports from pipeline results.
    """
    
    def __init__(self, results_dir: Path):
        """
        Initialize report generator.
        
        Args:
            results_dir: Directory containing results
        """
        self.results_dir = Path(results_dir)
    
    def generate(
        self,
        output_path: Path,
        results: Optional[Any] = None,
        format: str = "html",
    ) -> Path:
        """
        Generate report.
        
        Args:
            output_path: Output file path
            results: Pipeline results (optional)
            format: Report format (html, json, both)
            
        Returns:
            Path to generated report
        """
        logger.info(f"Generating report: {output_path}")
        
        if format == "json":
            return self._generate_json(output_path, results)
        else:
            return self._generate_html(output_path, results)
    
    def _generate_html(self, output_path: Path, results: Optional[Any]) -> Path:
        """Generate HTML report."""
        try:
            html = self._build_html(results)
            
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            
            logger.info(f"HTML report saved: {output_path}")
            return output_path
            
        except Exception as e:
            raise ReportError(f"Failed to generate HTML report: {e}")
    
    def _generate_json(self, output_path: Path, results: Optional[Any]) -> Path:
        """Generate JSON report."""
        import json
        
        data = results.to_dict() if results else {}
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        return output_path
    
    def _build_html(self, results: Optional[Any]) -> str:
        """Build HTML content."""
        stats = results.stats if results else {}
        
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SNPhylo2 Report</title>
    <style>
        :root {{
            --primary: #2c3e50;
            --secondary: #3498db;
            --success: #27ae60;
            --warning: #f39c12;
            --danger: #e74c3c;
            --light: #ecf0f1;
            --dark: #2c3e50;
        }}
        
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            color: var(--dark);
            background: #f5f6fa;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        header {{
            background: var(--primary);
            color: white;
            padding: 2rem;
            border-radius: 8px;
            margin-bottom: 2rem;
        }}
        
        header h1 {{
            font-size: 2rem;
            margin-bottom: 0.5rem;
        }}
        
        header .meta {{
            opacity: 0.8;
            font-size: 0.9rem;
        }}
        
        .card {{
            background: white;
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        .card h2 {{
            color: var(--primary);
            margin-bottom: 1rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid var(--light);
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
        }}
        
        .stat-box {{
            background: var(--light);
            padding: 1rem;
            border-radius: 6px;
            text-align: center;
        }}
        
        .stat-box .value {{
            font-size: 2rem;
            font-weight: bold;
            color: var(--secondary);
        }}
        
        .stat-box .label {{
            font-size: 0.9rem;
            color: #666;
            margin-top: 0.25rem;
        }}
        
        .progress-bar {{
            background: var(--light);
            height: 20px;
            border-radius: 10px;
            overflow: hidden;
            margin-top: 0.5rem;
        }}
        
        .progress-bar .fill {{
            background: var(--secondary);
            height: 100%;
            transition: width 0.3s ease;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 1rem;
        }}
        
        th, td {{
            padding: 0.75rem;
            text-align: left;
            border-bottom: 1px solid var(--light);
        }}
        
        th {{
            background: var(--light);
            font-weight: 600;
        }}
        
        tr:hover {{
            background: #f8f9fa;
        }}
        
        .badge {{
            display: inline-block;
            padding: 0.25rem 0.5rem;
            border-radius: 4px;
            font-size: 0.85rem;
            font-weight: 500;
        }}
        
        .badge-success {{
            background: #d4edda;
            color: #155724;
        }}
        
        .badge-warning {{
            background: #fff3cd;
            color: #856404;
        }}
        
        .footer {{
            text-align: center;
            padding: 2rem;
            color: #666;
            font-size: 0.9rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>SNPhylo2 Analysis Report</h1>
            <div class="meta">
                Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
                Version: 0.1.0
            </div>
        </header>
        
        <div class="card">
            <h2>Summary</h2>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="value">{stats.get('input_samples', 'N/A')}</div>
                    <div class="label">Samples</div>
                </div>
                <div class="stat-box">
                    <div class="value">{stats.get('input_snps', 'N/A'):,}</div>
                    <div class="label">Input SNPs</div>
                </div>
                <div class="stat-box">
                    <div class="value">{self._get_filtered_snps(stats):,}</div>
                    <div class="label">After Filtering</div>
                </div>
                <div class="stat-box">
                    <div class="value">{self._get_pruned_snps(stats):,}</div>
                    <div class="label">After LD Pruning</div>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>Filtering Results</h2>
            {self._build_filtering_section(stats)}
        </div>
        
        <div class="card">
            <h2>Tree Building</h2>
            {self._build_tree_section(stats)}
        </div>
        
        <div class="card">
            <h2>Output Files</h2>
            {self._build_files_section(results)}
        </div>
        
        <div class="footer">
            SNPhylo2: Next-Generation Phylogenomics Pipeline<br>
            <a href="https://github.com/snphylo2/snphylo2">https://github.com/snphylo2/snphylo2</a>
        </div>
    </div>
</body>
</html>"""
        
        return html
    
    def _get_filtered_snps(self, stats: Dict) -> int:
        """Get number of SNPs after filtering."""
        filtering = stats.get('filtering', {})
        return filtering.get('output_snps', 0)
    
    def _get_pruned_snps(self, stats: Dict) -> int:
        """Get number of SNPs after pruning."""
        pruning = stats.get('ld_pruning', {})
        return pruning.get('output_snps', 0)
    
    def _build_filtering_section(self, stats: Dict) -> str:
        """Build filtering results HTML."""
        filtering = stats.get('filtering', {})
        
        if not filtering:
            return "<p>No filtering performed.</p>"
        
        input_snps = filtering.get('input_snps', 0)
        output_snps = filtering.get('output_snps', 0)
        retention = filtering.get('retention_rate', 0) * 100
        
        filters = filtering.get('filters_applied', {})
        filter_rows = ""
        for name, count in filters.items():
            filter_rows += f"<tr><td>{name}</td><td>{count:,}</td></tr>"
        
        return f"""
        <div class="stats-grid">
            <div class="stat-box">
                <div class="value">{input_snps:,}</div>
                <div class="label">Input SNPs</div>
            </div>
            <div class="stat-box">
                <div class="value">{output_snps:,}</div>
                <div class="label">Retained SNPs</div>
            </div>
            <div class="stat-box">
                <div class="value">{retention:.1f}%</div>
                <div class="label">Retention Rate</div>
            </div>
        </div>
        <h3>Filters Applied</h3>
        <table>
            <thead>
                <tr><th>Filter</th><th>Variants Removed</th></tr>
            </thead>
            <tbody>{filter_rows}</tbody>
        </table>
        """
    
    def _build_tree_section(self, stats: Dict) -> str:
        """Build tree building results HTML."""
        tree = stats.get('tree', {})
        
        if not tree:
            return "<p>No tree building performed.</p>"
        
        engine = tree.get('engine', 'Unknown')
        model = tree.get('model', 'N/A')
        bootstrap = tree.get('bootstrap_support')
        
        bootstrap_str = f"{bootstrap:.1f}%" if bootstrap else "N/A"
        
        return f"""
        <table>
            <tr><th>Setting</th><th>Value</th></tr>
            <tr><td>Engine</td><td><span class="badge badge-success">{engine}</span></td></tr>
            <tr><td>Model</td><td>{model}</td></tr>
            <tr><td>Mean Bootstrap Support</td><td>{bootstrap_str}</td></tr>
        </table>
        """
    
    def _build_files_section(self, results: Optional[Any]) -> str:
        """Build output files section HTML."""
        if not results:
            return "<p>No output files.</p>"
        
        files = []
        if results.filtered_vcf:
            files.append(("Filtered VCF", results.filtered_vcf))
        if results.pruned_vcf:
            files.append(("Pruned VCF", results.pruned_vcf))
        if results.tree_file:
            files.append(("Tree (Newick)", results.tree_file))
        if results.report_path:
            files.append(("Report", results.report_path))
        
        rows = ""
        for name, path in files:
            rows += f"<tr><td>{name}</td><td><code>{path}</code></td></tr>"
        
        return f"""
        <table>
            <thead>
                <tr><th>File Type</th><th>Path</th></tr>
            </thead>
            <tbody>{rows}</tbody>
        </table>
        """
