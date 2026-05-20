"""
Command-line interface for SNPhylo2 using Click and Rich.
"""

import sys
from pathlib import Path
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich import box

from snphylo2 import __version__
from snphylo2.config import (
    SNPhylo2Config,
    load_config,
    save_config,
    create_default_config,
    InputConfig,
    MetadataConfig,
    FilteringConfig,
    LDPruningConfig,
    TreeConfig,
    OutputConfig,
    ComputeConfig,
    TreeEngine,
)
from snphylo2.exceptions import SNPhylo2Error
from snphylo2.pipeline import Pipeline
from snphylo2.utils.logging_utils import setup_logging

console = Console()


def print_header():
    """Print SNPhylo2 header."""
    header = Text()
    header.append("SNPhylo2 ", style="bold cyan")
    header.append(f"v{__version__}", style="cyan")
    header.append(" - Next-Generation Phylogenomics Pipeline", style="dim")
    console.print(Panel(header, box=box.ROUNDED))


@click.group()
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose logging")
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output")
@click.version_option(version=__version__, prog_name="snphylo2")
@click.pass_context
def main(ctx, verbose, quiet):
    """
    SNPhylo2: Next-Generation Phylogenomics Pipeline
    
    A modern, scalable, and reproducible pipeline for constructing phylogenetic
    trees from large SNP datasets.
    
    Use 'snphylo2 COMMAND --help' for detailed help on each command.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    
    setup_logging(verbose=verbose, quiet=quiet)
    
    if not quiet:
        print_header()


@main.command()
@click.option("-p", "--preset", 
              type=click.Choice(["default", "human", "plant", "microbial"]),
              default="default",
              help="Configuration preset")
@click.option("-o", "--output", type=click.Path(), default="snphylo2_config.yaml",
              help="Output configuration file path")
def init(preset: str, output: str):
    """
    Create a default configuration template.
    
    Example:
        snphylo2 init --preset plant -o my_config.yaml
    """
    try:
        output_path = Path(output)
        config = create_default_config(output_path, preset=preset)
        
        console.print(f"[green]Created configuration file:[/green] {output_path}")
        console.print(f"[dim]Preset: {preset}[/dim]")
        console.print("\n[dim]Edit this file to customize your analysis, then run:[/dim]")
        console.print(f"[cyan]snphylo2 run --config {output_path}[/cyan]")
        
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
@click.option("-v", "--vcf", type=click.Path(exists=True), help="Input VCF/BCF file")
@click.option("-c", "--config", type=click.Path(exists=True), help="Configuration file")
@click.option("-o", "--output", type=click.Path(), help="Output directory")
@click.option("--maf", type=float, help="Minimum minor allele frequency")
@click.option("--max-missing", type=float, help="Maximum missingness per SNP")
@click.option("--min-depth", type=int, help="Minimum read depth")
@click.option("--ld-window", type=int, help="LD pruning window size")
@click.option("--ld-r2", type=float, help="LD pruning R² threshold")
@click.option("--tree-engine", type=click.Choice([e.value for e in TreeEngine]),
              help="Tree building engine")
@click.option("--bootstrap", type=int, help="Number of bootstrap replicates")
@click.option("--threads", type=int, help="Number of threads")
@click.option("--report", type=click.Choice(["html", "json", "both"]),
              default="html", help="Report format")
@click.option("--keep-intermediates", is_flag=True, help="Keep intermediate files")
@click.pass_context
def run(ctx, vcf, config, output, maf, max_missing, min_depth, ld_window, ld_r2,
        tree_engine, bootstrap, threads, report, keep_intermediates):
    """
    Run the complete SNPhylo2 pipeline.
    
    This is the main command that runs quality control, filtering, LD pruning,
    tree building, and generates the final report.
    
    Examples:
        snphylo2 run -v input.vcf.gz --threads 16
        
        snphylo2 run -v input.vcf.gz -c config.yaml -o results/
        
        snphylo2 run -v input.vcf.gz --maf 0.05 --tree-engine iqtree2 --bootstrap 1000
    """
    try:
        # Load or create configuration
        if config:
            cfg = load_config(config)
            console.print(f"[dim]Loaded configuration from {config}[/dim]")
        else:
            cfg = SNPhylo2Config(
                input=InputConfig(path=Path(vcf) if vcf else Path("input.vcf.gz")),
                output=OutputConfig(directory=Path(output) if output else Path(".")),
                compute=ComputeConfig(threads=threads or 1),
            )
        
        # Override with command-line options
        if vcf:
            cfg.input.path = Path(vcf)
        if output:
            cfg.output.directory = Path(output)
            cfg.output.prefix = f"{Path(output).name}"
        if maf is not None:
            cfg.filtering.maf.min = maf
        if max_missing is not None:
            cfg.filtering.missingness.max_per_snp = max_missing
        if min_depth is not None:
            cfg.filtering.depth.min = min_depth
        if ld_window is not None:
            cfg.ld_pruning.window_size = ld_window
        if ld_r2 is not None:
            cfg.ld_pruning.r2_threshold = ld_r2
        if tree_engine:
            cfg.tree.engine = TreeEngine(tree_engine)
        if bootstrap is not None:
            cfg.tree.bootstrap.replicates = bootstrap
        if threads is not None:
            cfg.compute.threads = threads
        if report:
            cfg.reporting.format = report  # type: ignore
        if keep_intermediates:
            cfg.output.keep_intermediates = True
        
        # Validate input
        if not cfg.input.path.exists():
            raise click.BadParameter(f"Input file not found: {cfg.input.path}")
        
        # Run pipeline
        if not ctx.obj.get("quiet"):
            console.print("\n[bold]Starting SNPhylo2 pipeline...[/bold]")
            console.print(f"[dim]Input: {cfg.input.path}[/dim]")
            console.print(f"[dim]Output: {cfg.output.directory}[/dim]")
            console.print(f"[dim]Threads: {cfg.compute.threads}[/dim]")
        
        pipeline = Pipeline(cfg)
        results = pipeline.run()
        
        if not ctx.obj.get("quiet"):
            console.print(f"\n[green]Pipeline completed successfully![/green]")
            console.print(f"[dim]Results saved to: {results['output_dir']}[/dim]")
            if "report_path" in results:
                console.print(f"[dim]Report: {results['report_path']}[/dim]")
        
    except SNPhylo2Error as e:
        console.print(f"\n[red]Error: {e.message}[/red]")
        if e.details:
            for key, value in e.details.items():
                console.print(f"[dim]  {key}: {value}[/dim]")
        sys.exit(1)
    except Exception as e:
        console.print(f"\n[red]Unexpected error: {e}[/red]")
        if ctx.obj.get("verbose"):
            console.print_exception()
        sys.exit(1)


@main.command()
@click.option("-v", "--vcf", type=click.Path(exists=True), required=True,
              help="Input VCF/BCF file")
@click.option("-o", "--output", type=click.Path(), help="Output QC report file")
@click.pass_context
def qc(ctx, vcf: str, output: Optional[str]):
    """
    Run quality control analysis only.
    
    Generates QC metrics and visualizations without running the full pipeline.
    
    Example:
        snphylo2 qc -v input.vcf.gz -o qc_report.html
    """
    from snphylo2.qc.qc_module import QCModule
    
    try:
        if not ctx.obj.get("quiet"):
            console.print("[bold]Running QC analysis...[/bold]")
        
        qc_module = QCModule(Path(vcf))
        report = qc_module.run_analysis()
        
        output_path = Path(output) if output else Path("qc_report.html")
        report.save(output_path)
        
        if not ctx.obj.get("quiet"):
            console.print(f"[green]QC report saved to:[/green] {output_path}")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
@click.option("-v", "--vcf", type=click.Path(exists=True), required=True,
              help="Input VCF/BCF file")
@click.option("-o", "--output", type=click.Path(), required=True,
              help="Output filtered VCF file")
@click.option("--maf", type=float, default=0.05, help="Minimum MAF")
@click.option("--max-missing", type=float, default=0.2, help="Maximum missingness")
@click.option("--min-depth", type=int, default=5, help="Minimum depth")
@click.option("--biallelic-only", is_flag=True, default=True, help="Biallelic SNPs only")
@click.pass_context
def filter(ctx, vcf: str, output: str, maf: float, max_missing: float, 
           min_depth: int, biallelic_only: bool):
    """
    Apply variant and sample filters.
    
    Example:
        snphylo2 filter -v input.vcf.gz -o filtered.vcf.gz --maf 0.05 --max-missing 0.2
    """
    from snphylo2.filtering.variant_filters import FilterPipeline
    from snphylo2.config import FilteringConfig, MAFFilter, MissingnessFilter, DepthFilter
    
    try:
        if not ctx.obj.get("quiet"):
            console.print("[bold]Running filters...[/bold]")
        
        config = FilteringConfig(
            maf=MAFFilter(min=maf),
            missingness=MissingnessFilter(max_per_snp=max_missing),
            depth=DepthFilter(min=min_depth),
            biallelic_only=biallelic_only,
        )
        
        pipeline = FilterPipeline(config)
        stats = pipeline.run(Path(vcf), Path(output))
        
        if not ctx.obj.get("quiet"):
            console.print(f"[green]Filtering complete:[/green]")
            console.print(f"[dim]  Input SNPs: {stats['input_snps']:,}[/dim]")
            console.print(f"[dim]  Output SNPs: {stats['output_snps']:,}[/dim]")
            console.print(f"[dim]  Retained: {stats['retention_rate']:.1%}[/dim]")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
@click.option("-i", "--input", "input_file", type=click.Path(exists=True), required=True,
              help="Input filtered VCF file")
@click.option("-o", "--output", type=click.Path(), required=True,
              help="Output pruned VCF file")
@click.option("--window", type=int, default=50, help="Window size")
@click.option("--step", type=int, default=10, help="Step size")
@click.option("--r2", type=float, default=0.2, help="R² threshold")
@click.option("-t", "--threads", type=int, default=1, help="Number of threads")
@click.pass_context
def prune(ctx, input_file: str, output: str, window: int, step: int, r2: float, threads: int):
    """
    Perform LD pruning on filtered variants.
    
    Example:
        snphylo2 prune -i filtered.vcf.gz -o pruned.vcf.gz --window 50 --r2 0.2
    """
    from snphylo2.pruning.ld_pruner import LDPruner
    
    try:
        if not ctx.obj.get("quiet"):
            console.print("[bold]Running LD pruning...[/bold]")
        
        config = LDPruningConfig(
            window_size=window,
            step_size=step,
            r2_threshold=r2,
        )
        
        pruner = LDPruner(config, threads=threads)
        stats = pruner.run(Path(input_file), Path(output))
        
        if not ctx.obj.get("quiet"):
            console.print(f"[green]LD pruning complete:[/green]")
            console.print(f"[dim]  Input SNPs: {stats['input_snps']:,}[/dim]")
            console.print(f"[dim]  Output SNPs: {stats['output_snps']:,}[/dim]")
            console.print(f"[dim]  Pruned: {stats['pruned_count']:,}[/dim]")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
@click.option("-i", "--input", "input_file", type=click.Path(exists=True), required=True,
              help="Input VCF or alignment file")
@click.option("-o", "--output", type=click.Path(), required=True,
              help="Output tree file (Newick format)")
@click.option("--engine", type=click.Choice([e.value for e in TreeEngine]),
              default="iqtree2", help="Tree building engine")
@click.option("--model", type=str, default="GTR+ASC", help="Substitution model")
@click.option("--bootstrap", type=int, default=1000, help="Bootstrap replicates")
@click.option("-t", "--threads", type=int, default=1, help="Number of threads")
@click.option("--outgroup", type=str, help="Outgroup sample name")
@click.pass_context
def tree(ctx, input_file: str, output: str, engine: str, model: str,
         bootstrap: int, threads: int, outgroup: Optional[str]):
    """
    Build phylogenetic tree from filtered variants.
    
    Example:
        snphylo2 tree -i pruned.vcf.gz -o tree.nwk --engine iqtree2 --bootstrap 1000
    """
    try:
        if not ctx.obj.get("quiet"):
            console.print(f"[bold]Building tree with {engine}...[/bold]")
        
        config = TreeConfig(
            engine=TreeEngine(engine),
            model_selection=False,
            candidate_models=[model],
            bootstrap=LDPruningConfig.__pydantic_fields__["bootstrap"].annotation,  # type: ignore
            threads=threads,
            outgroup=outgroup,
        )
        
        from snphylo2.tree.tree_builder import TreeBuilder
        builder = TreeBuilder(config)
        result = builder.build(Path(input_file), Path(output))
        
        if not ctx.obj.get("quiet"):
            console.print(f"[green]Tree building complete:[/green]")
            console.print(f"[dim]  Tree file: {result['tree_file']}[/dim]")
            if result.get('model'):
                console.print(f"[dim]  Model: {result['model']}[/dim]")
            if result.get('bootstrap_support'):
                console.print(f"[dim]  Mean bootstrap: {result['bootstrap_support']:.1f}%[/dim]")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
@click.option("-v", "--vcf", type=click.Path(exists=True), required=True,
              help="Input VCF file")
@click.option("-o", "--output", type=click.Path(), default="ld_decay",
              help="Output directory/prefix")
@click.option("-d", "--max-distance", type=int, default=200000,
              help="Maximum distance (bp)")
@click.option("--bin-size", type=int, default=5000,
              help="Bin size (bp)")
@click.option("-m", "--metadata", type=click.Path(exists=True),
              help="Metadata file with population info")
@click.option("--plot", is_flag=True, default=True,
              help="Generate plots")
@click.pass_context
def ld_decay(ctx, vcf: str, output: str, max_distance: int, bin_size: int,
             metadata: Optional[str], plot: bool):
    """
    Analyze LD decay patterns (PopLDdecay-style).
    
    Calculates r² decay with physical distance for population genetics analysis.
    
    Example:
        snphylo2 ld-decay -v input.vcf.gz -o ld_results/
        
        snphylo2 ld-decay -v input.vcf.gz -m metadata.tsv --plot
    """
    try:
        from snphylo2.population.ld_decay import LDDecayAnalyzer
        from snphylo2.visualization.ld_plots import LDPlotter
        import pandas as pd
        
        if not ctx.obj.get("quiet"):
            console.print("[bold]Analyzing LD decay...[/bold]")
        
        # Load metadata if provided
        sample_groups = None
        if metadata:
            df = pd.read_csv(metadata, sep='\t')
            sample_col = df.columns[0]
            pop_col = None
            for col in ['population', 'subspecies', 'ecotype']:
                if col in df.columns:
                    pop_col = col
                    break
            
            if pop_col:
                sample_groups = {}
                for _, row in df.iterrows():
                    pop = row[pop_col]
                    if pop not in sample_groups:
                        sample_groups[pop] = []
                    sample_groups[pop].append(row[sample_col])
                
                if not ctx.obj.get("quiet"):
                    console.print(f"[dim]Found {len(sample_groups)} populations[/dim]")
        
        # Analyze
        analyzer = LDDecayAnalyzer(
            max_distance=max_distance,
            bin_size=bin_size,
        )
        
        results = analyzer.analyze(
            Path(vcf),
            sample_groups=sample_groups,
        )
        
        # Save results
        output_path = Path(output)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for pop_name, result in results.items():
            data_file = output_path / f"ld_decay_{pop_name}.tsv"
            result.to_dataframe().to_csv(data_file, sep='\t', index=False)
            
            if not ctx.obj.get("quiet"):
                console.print(f"[dim]  Saved: {data_file}[/dim]")
        
        # Generate plots
        if plot and results:
            plotter = LDPlotter(width=1200, height=900, dpi=300)
            
            plot_file = output_path / f"ld_decay_plot.pdf"
            plotter.plot_ld_decay(
                results,
                output_path=plot_file,
                title="Linkage Disequilibrium Decay",
                show_confidence=True,
                manuscript_ready=True,
            )
            
            if not ctx.obj.get("quiet"):
                console.print(f"[green]Plot saved:[/green] {plot_file}")
        
        # Summary
        if not ctx.obj.get("quiet"):
            console.print("\n[bold]LD Decay Summary:[/bold]")
            for pop_name, result in results.items():
                if result.half_decay_distance:
                    console.print(f"  {pop_name}: half-decay = {result.half_decay_distance/1000:.1f} kb")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        sys.exit(1)


@main.command()
@click.option("-d", "--directory", type=click.Path(exists=True), required=True,
              help="Results directory")
@click.option("-o", "--output", type=click.Path(), default="report.html",
              help="Output report file")
@click.option("-f", "--format", "fmt", type=click.Choice(["html", "json", "both"]),
              default="html", help="Report format")
@click.pass_context
def report(ctx, directory: str, output: str, fmt: str):
    """
    Generate HTML report from existing results.
    
    Example:
        snphylo2 report -d results/ -o final_report.html
    """
    try:
        if not ctx.obj.get("quiet"):
            console.print("[bold]Generating report...[/bold]")
        
        from snphylo2.report.html_report import HTMLReport
        
        reporter = HTMLReport(Path(directory))
        report_path = reporter.generate(Path(output), format=fmt)
        
        if not ctx.obj.get("quiet"):
            console.print(f"[green]Report generated:[/green] {report_path}")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


@main.command()
def validate():
    """
    Check installation and dependencies.
    
    Verifies that all required external tools are installed and accessible.
    """
    from snphylo2.utils.validators import validate_installation
    
    console.print("[bold]Validating SNPhylo2 installation...[/bold]\n")
    
    results = validate_installation()
    
    all_ok = True
    for category, tools in results.items():
        console.print(f"[bold]{category}:[/bold]")
        for tool, status in tools.items():
            if status["available"]:
                console.print(f"  [green]✓[/green] {tool}: {status['version'] or 'found'}")
            else:
                console.print(f"  [red]✗[/red] {tool}: {status['error'] or 'not found'}")
                all_ok = False
    
    if all_ok:
        console.print("\n[green]All dependencies are available![/green]")
    else:
        console.print("\n[yellow]Some dependencies are missing. See documentation for installation.[/yellow]")


@main.command()
@click.option("-t1", "--tree1", type=click.Path(exists=True), required=True,
              help="First tree (Newick format)")
@click.option("-t2", "--tree2", type=click.Path(exists=True), required=True,
              help="Second tree (Newick format)")
@click.option("-o", "--output", type=click.Path(), help="Output comparison file")
def compare(tree1: str, tree2: str, output: Optional[str]):
    """
    Compare two phylogenetic trees using Robinson-Foulds distance.
    
    Example:
        snphylo2 compare -t1 tree1.nwk -t2 tree2.nwk
    """
    try:
        from snphylo2.tree.tree_comparison import compare_trees
        
        console.print("[bold]Comparing trees...[/bold]")
        
        result = compare_trees(Path(tree1), Path(tree2))
        
        console.print(f"\n[green]Comparison results:[/green]")
        console.print(f"[dim]  Robinson-Foulds distance: {result['rf_distance']}[/dim]")
        console.print(f"[dim]  Normalized RF: {result['normalized_rf']:.4f}[/dim]")
        console.print(f"[dim]  Shared splits: {result['shared_splits']}[/dim]")
        console.print(f"[dim]  Unique to tree1: {result['unique_tree1']}[/dim]")
        console.print(f"[dim]  Unique to tree2: {result['unique_tree2']}[/dim]")
        
        if output:
            import json
            with open(output, "w") as f:
                json.dump(result, f, indent=2)
            console.print(f"\n[dim]Results saved to: {output}[/dim]")
            
    except SNPhylo2Error as e:
        console.print(f"[red]Error: {e.message}[/red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
