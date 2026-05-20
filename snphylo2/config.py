"""
Configuration management for SNPhylo2 using Pydantic.
"""

from pathlib import Path
from typing import Literal, List, Optional, Dict, Any, Union
from enum import Enum

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

from snphylo2.exceptions import ConfigurationError


class InputFormat(str, Enum):
    """Supported input formats."""
    AUTO = "auto"
    VCF = "vcf"
    BCF = "bcf"
    PLINK = "plink"
    HAPMAP = "hapmap"
    GDS = "gds"
    FASTA = "fasta"
    PHYLIP = "phylip"
    SIMPLE_SNP = "simple_snp"


class TreeEngine(str, Enum):
    """Supported tree-building engines."""
    IQTREE2 = "iqtree2"
    RAXML_NG = "raxml_ng"
    FASTTREE = "fasttree"
    PHYML = "phyml"
    DNAML = "dnaml"  # For backwards compatibility


class BootstrapMethod(str, Enum):
    """Bootstrap methods."""
    STANDARD = "standard"
    ULTRAFAST = "ultrafast"
    SHALRT = "shalrt"


class MissingDataStrategy(str, Enum):
    """Strategies for handling missing data."""
    FILTER = "filter"
    AMBIGUITY = "ambiguity"
    IMPUTE_MAJOR = "impute_major"
    IMPUTE_REFERENCE = "impute_reference"


class ComputeMode(str, Enum):
    """Compute execution modes."""
    LOCAL = "local"
    SLURM = "slurm"
    PBS = "pbs"
    SGE = "sge"
    CLOUD = "cloud"


class MAFFilter(BaseModel):
    """MAF filtering configuration."""
    min: float = Field(default=0.05, ge=0.0, le=1.0)
    max: float = Field(default=1.0, ge=0.0, le=1.0)


class MissingnessFilter(BaseModel):
    """Missingness filtering configuration."""
    max_per_snp: float = Field(default=0.2, ge=0.0, le=1.0)
    max_per_sample: float = Field(default=0.5, ge=0.0, le=1.0)


class DepthFilter(BaseModel):
    """Depth filtering configuration."""
    min: int = Field(default=5, ge=0)
    max: Optional[int] = Field(default=None, ge=0)


class FilteringConfig(BaseModel):
    """Variant and sample filtering configuration."""
    maf: MAFFilter = Field(default_factory=MAFFilter)
    missingness: MissingnessFilter = Field(default_factory=MissingnessFilter)
    depth: DepthFilter = Field(default_factory=DepthFilter)
    genotype_quality: Optional[int] = Field(default=20, ge=0)
    biallelic_only: bool = True
    ts_tv_ratio: Optional[float] = Field(default=None, ge=0)
    hwe_pvalue: Optional[float] = Field(default=1e-6, ge=0, le=1)
    exclude_singletons: bool = False
    mask_bed: Optional[Path] = None
    include_chromosomes: Optional[List[str]] = None
    exclude_chromosomes: Optional[List[str]] = Field(default_factory=list)
    thin_distance: Optional[int] = None  # bp
    chromosome_balance: bool = False
    asc_bias_aware: bool = False


class LDPruningConfig(BaseModel):
    """LD pruning configuration."""
    method: Literal["windowed", "pairwise", "clumping"] = "windowed"
    window_size: int = Field(default=50, ge=1)
    step_size: int = Field(default=10, ge=1)
    r2_threshold: float = Field(default=0.2, ge=0.0, le=1.0)
    by_chromosome: bool = True
    n_representative: Optional[int] = None


class MissingDataConfig(BaseModel):
    """Missing data handling configuration."""
    strategy: MissingDataStrategy = MissingDataStrategy.FILTER
    sensitivity_analysis: bool = False
    thresholds: List[float] = Field(default_factory=lambda: [0.1, 0.2, 0.3])


class BootstrapConfig(BaseModel):
    """Bootstrap configuration."""
    method: BootstrapMethod = BootstrapMethod.ULTRAFAST
    replicates: int = Field(default=1000, ge=1)
    seed: Optional[int] = None


class TreeConfig(BaseModel):
    """Tree building configuration."""
    engine: TreeEngine = TreeEngine.IQTREE2
    model_selection: bool = True
    candidate_models: List[str] = Field(default_factory=lambda: ["GTR+ASC", "GTR+ASC+G"])
    model_test_criterion: Literal["AIC", "AICc", "BIC"] = "BIC"
    asc_bias_correction: bool = True
    invariant_sites: Optional[Literal["estimate", "empirical"]] = "estimate"
    bootstrap: BootstrapConfig = Field(default_factory=BootstrapConfig)
    outgroup: Optional[str] = None
    threads: Optional[int] = None  # Auto-detect if None


class PopulationConfig(BaseModel):
    """Population genetics analysis configuration."""
    run_pca: bool = True
    pca_components: int = Field(default=10, ge=1)
    run_admixture: bool = False
    admixture_k: List[int] = Field(default_factory=lambda: [2, 3, 4, 5])
    calculate_fst: bool = True
    calculate_ibs: bool = True
    calculate_kinship: bool = True


class VisualizationConfig(BaseModel):
    """Visualization configuration."""
    tree_formats: List[str] = Field(default_factory=lambda: ["pdf", "svg", "html"])
    tree_layout: Literal["rectangular", "circular", "radial"] = "rectangular"
    show_bootstrap: bool = True
    color_by: Optional[str] = None
    pca_color_by: Optional[str] = None


class ReportingConfig(BaseModel):
    """Report generation configuration."""
    format: Literal["html", "json", "both"] = "html"
    include_command: bool = True
    include_config: bool = True
    include_versions: bool = True


class ComputeConfig(BaseModel):
    """Compute resource configuration."""
    mode: ComputeMode = ComputeMode.LOCAL
    threads: int = Field(default=1, ge=1)
    memory: str = "16G"
    tmp_dir: Path = Field(default_factory=lambda: Path("/tmp"))
    
    # HPC-specific settings
    slurm_partition: Optional[str] = None
    slurm_nodes: int = 1
    slurm_time: str = "24:00:00"
    
    # Cloud-specific settings
    cloud_provider: Optional[Literal["aws", "gcp", "azure"]] = None
    cloud_instance_type: Optional[str] = None
    cloud_spot_instances: bool = True


class InputConfig(BaseModel):
    """Input file configuration."""
    format: InputFormat = InputFormat.AUTO
    path: Path
    index: Optional[Path] = None
    region: Optional[str] = None  # e.g., "chr1:1-1000000"
    
    @field_validator("path")
    @classmethod
    def path_must_exist(cls, v: Path) -> Path:
        if not v.exists():
            raise ConfigurationError(f"Input file does not exist: {v}")
        return v


class MetadataConfig(BaseModel):
    """Metadata file configuration."""
    path: Optional[Path] = None
    sample_id_column: str = "sample_id"
    group_column: Optional[str] = None


class OutputConfig(BaseModel):
    """Output configuration."""
    prefix: str = "snphylo2_output"
    directory: Path = Field(default_factory=lambda: Path("."))
    keep_intermediates: bool = False
    compress_outputs: bool = True


class SNPhylo2Config(BaseModel):
    """Main SNPhylo2 configuration model."""
    
    version: str = "2.0"
    
    input: InputConfig
    metadata: MetadataConfig = Field(default_factory=MetadataConfig)
    filtering: FilteringConfig = Field(default_factory=FilteringConfig)
    ld_pruning: LDPruningConfig = Field(default_factory=LDPruningConfig)
    missing_data: MissingDataConfig = Field(default_factory=MissingDataConfig)
    tree: TreeConfig = Field(default_factory=TreeConfig)
    population: PopulationConfig = Field(default_factory=PopulationConfig)
    visualization: VisualizationConfig = Field(default_factory=VisualizationConfig)
    reporting: ReportingConfig = Field(default_factory=ReportingConfig)
    compute: ComputeConfig = Field(default_factory=ComputeConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    
    @model_validator(mode="after")
    def validate_consistency(self) -> "SNPhylo2Config":
        """Validate cross-field consistency."""
        # Check that MAF min < max
        if self.filtering.maf.min >= self.filtering.maf.max:
            raise ConfigurationError("MAF min must be less than MAF max")
        
        # Check that LD pruning window > step
        if self.ld_pruning.window_size <= self.ld_pruning.step_size:
            raise ConfigurationError("LD pruning window_size must be greater than step_size")
        
        return self


def load_config(path: Union[str, Path]) -> SNPhylo2Config:
    """
    Load configuration from YAML file.
    
    Args:
        path: Path to YAML configuration file
        
    Returns:
        SNPhylo2Config instance
        
    Raises:
        ConfigurationError: If file is invalid or cannot be loaded
    """
    path = Path(path)
    
    if not path.exists():
        raise ConfigurationError(f"Configuration file not found: {path}")
    
    try:
        with open(path, "r") as f:
            data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Invalid YAML syntax: {e}")
    except Exception as e:
        raise ConfigurationError(f"Failed to read configuration file: {e}")
    
    try:
        return SNPhylo2Config(**data)
    except Exception as e:
        raise ConfigurationError(f"Invalid configuration: {e}")


def save_config(config: SNPhylo2Config, path: Union[str, Path]) -> None:
    """
    Save configuration to YAML file.
    
    Args:
        config: Configuration to save
        path: Output file path
    """
    path = Path(path)
    
    # Convert to dict, handling Path objects
    def convert_paths(obj: Any) -> Any:
        if isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, dict):
            return {k: convert_paths(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_paths(item) for item in obj]
        elif isinstance(obj, BaseModel):
            return convert_paths(obj.model_dump())
        return obj
    
    data = convert_paths(config.model_dump())
    
    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def create_default_config(path: Union[str, Path], preset: str = "default") -> SNPhylo2Config:
    """
    Create a default configuration file.
    
    Args:
        path: Output file path
        preset: Preset name (default, human, plant, microbial)
        
    Returns:
        SNPhylo2Config instance
    """
    presets = {
        "default": SNPhylo2Config(
            input=InputConfig(path=Path("input.vcf.gz")),
        ),
        "human": SNPhylo2Config(
            input=InputConfig(path=Path("input.vcf.gz")),
            filtering=FilteringConfig(
                exclude_chromosomes=["chrM", "chrY"],
                maf=MAFFilter(min=0.01),  # Lower MAF for humans
            ),
        ),
        "plant": SNPhylo2Config(
            input=InputConfig(path=Path("input.vcf.gz")),
            filtering=FilteringConfig(
                maf=MAFFilter(min=0.05),
                depth=DepthFilter(min=3),  # Lower depth for plants
            ),
            ld_pruning=LDPruningConfig(
                window_size=100,  # Larger windows for plants
                r2_threshold=0.1,
            ),
        ),
        "microbial": SNPhylo2Config(
            input=InputConfig(path=Path("input.vcf.gz")),
            filtering=FilteringConfig(
                maf=MAFFilter(min=0.0),  # Include rare variants
                missingness=MissingnessFilter(max_per_snp=0.1),
            ),
            tree=TreeConfig(
                engine=TreeEngine.IQTREE2,
                model_selection=False,
                candidate_models=["GTR+G"],
            ),
        ),
    }
    
    if preset not in presets:
        raise ConfigurationError(f"Unknown preset: {preset}. Available: {list(presets.keys())}")
    
    config = presets[preset]
    save_config(config, path)
    return config
