"""
SNPhylo2: Next-Generation Phylogenomic Pipeline from SNP Data

A modern, scalable, and reproducible pipeline for constructing phylogenetic
trees from large SNP datasets with integrated population-genomics analyses.
"""

__version__ = "0.1.0"
__author__ = "Ankush Sharma"
__license__ = "MIT"

from snphylo2.config import SNPhylo2Config, load_config
from snphylo2.exceptions import (
    SNPhylo2Error,
    ConfigurationError,
    InputError,
    FilterError,
    TreeError,
)

__all__ = [
    "__version__",
    "SNPhylo2Config",
    "load_config",
    "SNPhylo2Error",
    "ConfigurationError",
    "InputError",
    "FilterError",
    "TreeError",
]
