"""
Installation validation utilities.
"""

import shutil
import subprocess
from typing import Dict, Any

from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


def validate_installation() -> Dict[str, Dict[str, Any]]:
    """
    Validate that all required external tools are installed.
    
    Returns:
        Dictionary with validation results for each tool
    """
    results = {
        "Core Tools": {},
        "Tree Engines": {},
        "QC/Filtering": {},
    }
    
    # Core tools
    results["Core Tools"]["python"] = _check_tool("python", ["--version"])
    results["Core Tools"]["tabix"] = _check_tool("tabix")
    results["Core Tools"]["bgzip"] = _check_tool("bgzip")
    
    # Tree engines
    results["Tree Engines"]["iqtree2"] = _check_tool("iqtree2", ["-version"])
    results["Tree Engines"]["raxml-ng"] = _check_tool("raxml-ng", ["--version"])
    results["Tree Engines"]["FastTree"] = _check_tool("FastTree")
    
    # QC/Filtering tools
    results["QC/Filtering"]["plink2"] = _check_tool("plink2", ["--version"])
    results["QC/Filtering"]["bcftools"] = _check_tool("bcftools", ["--version"])
    
    return results


def _check_tool(name: str, version_args: list = None) -> Dict[str, Any]:
    """
    Check if a tool is available.
    
    Args:
        name: Tool executable name
        version_args: Arguments to get version
        
    Returns:
        Dictionary with availability status and version info
    """
    result = {
        "available": False,
        "version": None,
        "error": None,
    }
    
    # Check if executable exists
    executable = shutil.which(name)
    if not executable:
        result["error"] = "Not found in PATH"
        return result
    
    result["available"] = True
    
    # Try to get version
    if version_args:
        try:
            cmd = [name] + version_args
            output = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=5,
            )
            
            # Parse version from stdout or stderr
            version_text = output.stdout or output.stderr
            if version_text:
                # Extract first line as version
                result["version"] = version_text.strip().split('\n')[0][:50]
                
        except subprocess.TimeoutExpired:
            result["error"] = "Version check timed out"
        except Exception as e:
            result["error"] = str(e)
    
    return result


def check_vcf_index(vcf_path: str) -> bool:
    """
    Check if VCF file has a tabix index.
    
    Args:
        vcf_path: Path to VCF file
        
    Returns:
        True if indexed, False otherwise
    """
    from pathlib import Path
    
    vcf = Path(vcf_path)
    
    # Check for .tbi or .csi index
    tbi = vcf.with_suffix(vcf.suffix + ".tbi")
    csi = vcf.with_suffix(vcf.suffix + ".csi")
    
    return tbi.exists() or csi.exists()


def check_file_format(file_path: str) -> str:
    """
    Detect file format from extension and content.
    
    Args:
        file_path: Path to file
        
    Returns:
        Detected format string
    """
    from pathlib import Path
    
    path = Path(file_path)
    suffix = path.suffix.lower()
    
    # Check compressed extensions first
    if suffix == ".gz":
        inner_suffix = path.stem.lower().split('.')[-1] if '.' in path.stem else ''
        if inner_suffix in ["vcf"]:
            return "vcf.gz"
    
    format_map = {
        ".vcf": "vcf",
        ".bcf": "bcf",
        ".bed": "plink",
        ".bim": "plink",
        ".fam": "plink",
        ".gds": "gds",
        ".fasta": "fasta",
        ".fa": "fasta",
        ".phy": "phylip",
        ".phylip": "phylip",
        ".nwk": "newick",
        ".newick": "newick",
        ".nex": "nexus",
        ".nexus": "nexus",
    }
    
    return format_map.get(suffix, "unknown")
