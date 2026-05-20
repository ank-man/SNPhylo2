"""
Tree comparison utilities including Robinson-Foulds distance.
"""

from pathlib import Path
from typing import Dict, Any, Tuple

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False

try:
    from Bio import Phylo
    from Bio.Phylo.TreeDistance import robinson_foulds
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from snphylo2.exceptions import TreeError
from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


def compare_trees(tree1_path: Path, tree2_path: Path) -> Dict[str, Any]:
    """
    Compare two phylogenetic trees using Robinson-Foulds distance.
    
    Args:
        tree1_path: First tree (Newick format)
        tree2_path: Second tree (Newick format)
        
    Returns:
        Dictionary with comparison metrics
    """
    logger.info(f"Comparing trees: {tree1_path} vs {tree2_path}")
    
    if not ETE3_AVAILABLE and not BIOPYTHON_AVAILABLE:
        raise TreeError("Either ete3 or Bio.Phylo required for tree comparison")
    
    try:
        if BIOPYTHON_AVAILABLE:
            return _compare_biopython(tree1_path, tree2_path)
        else:
            return _compare_ete3(tree1_path, tree2_path)
    except Exception as e:
        raise TreeError(f"Tree comparison failed: {e}")


def _compare_biopython(tree1_path: Path, tree2_path: Path) -> Dict[str, Any]:
    """Compare trees using Bio.Phylo."""
    tree1 = Phylo.read(tree1_path, 'newick')
    tree2 = Phylo.read(tree2_path, 'newick')
    
    # Robinson-Foulds distance
    rf_result = robinson_foulds(tree1, tree2)
    rf_distance = rf_result[0]
    
    # Get total number of splits
    n_splits_tree1 = len(list(tree1.get_nonterminals()))
    n_splits_tree2 = len(list(tree2.get_nonterminals()))
    max_possible_rf = n_splits_tree1 + n_splits_tree2
    
    normalized_rf = rf_distance / max_possible_rf if max_possible_rf > 0 else 0.0
    
    # Shared and unique splits
    shared = (max_possible_rf - rf_distance) // 2
    unique1 = n_splits_tree1 - shared
    unique2 = n_splits_tree2 - shared
    
    return {
        'rf_distance': int(rf_distance),
        'normalized_rf': float(normalized_rf),
        'max_possible_rf': int(max_possible_rf),
        'shared_splits': int(shared),
        'unique_tree1': int(unique1),
        'unique_tree2': int(unique2),
        'tree1_splits': int(n_splits_tree1),
        'tree2_splits': int(n_splits_tree2),
    }


def _compare_ete3(tree1_path: Path, tree2_path: Path) -> Dict[str, Any]:
    """Compare trees using ete3."""
    tree1 = Tree(str(tree1_path))
    tree2 = Tree(str(tree2_path))
    
    # Robinson-Foulds distance
    rf_result = tree1.robinson_foulds(tree2, unrooted_trees=True)
    rf_distance = rf_result[0]
    max_possible_rf = rf_result[1]
    
    normalized_rf = rf_distance / max_possible_rf if max_possible_rf > 0 else 0.0
    
    # Get splits info
    n_splits_tree1 = len([n for n in tree1.traverse() if not n.is_leaf()])
    n_splits_tree2 = len([n for n in tree2.traverse() if not n.is_leaf()])
    
    shared = (max_possible_rf - rf_distance) // 2
    unique1 = n_splits_tree1 - shared
    unique2 = n_splits_tree2 - shared
    
    return {
        'rf_distance': int(rf_distance),
        'normalized_rf': float(normalized_rf),
        'max_possible_rf': int(max_possible_rf),
        'shared_splits': int(shared),
        'unique_tree1': int(unique1),
        'unique_tree2': int(unique2),
        'tree1_splits': int(n_splits_tree1),
        'tree2_splits': int(n_splits_tree2),
    }


def get_tree_statistics(tree_path: Path) -> Dict[str, Any]:
    """
    Calculate basic statistics for a tree.
    
    Args:
        tree_path: Path to Newick tree
        
    Returns:
        Dictionary with tree statistics
    """
    if not ETE3_AVAILABLE:
        raise TreeError("ete3 required for tree statistics")
    
    try:
        tree = Tree(str(tree_path))
        
        n_leaves = len(tree)
        n_internal = len([n for n in tree.traverse() if not n.is_leaf()])
        total_nodes = n_leaves + n_internal
        
        # Tree height (max distance from root to leaf)
        heights = [tree.get_distance(leaf) for leaf in tree.get_leaves()]
        max_height = max(heights) if heights else 0
        mean_height = sum(heights) / len(heights) if heights else 0
        
        # Bootstrap support statistics
        supports = [n.support for n in tree.traverse() 
                   if not n.is_leaf() and n.support is not None]
        
        support_stats = {}
        if supports:
            support_stats = {
                'mean_bootstrap': sum(supports) / len(supports),
                'median_bootstrap': sorted(supports)[len(supports) // 2],
                'min_bootstrap': min(supports),
                'max_bootstrap': max(supports),
                'bootstrap_>70': sum(1 for s in supports if s > 70) / len(supports),
                'bootstrap_>90': sum(1 for s in supports if s > 90) / len(supports),
            }
        
        return {
            'n_leaves': n_leaves,
            'n_internal_nodes': n_internal,
            'total_nodes': total_nodes,
            'max_height': max_height,
            'mean_height': mean_height,
            **support_stats,
        }
        
    except Exception as e:
        raise TreeError(f"Failed to calculate tree statistics: {e}")


def find_long_branches(tree_path: Path, threshold: float = 0.1) -> list:
    """
    Identify samples on unusually long branches (potential outliers).
    
    Args:
        tree_path: Path to Newick tree
        threshold: Branch length threshold (as multiple of mean)
        
    Returns:
        List of (sample_name, branch_length) tuples
    """
    if not ETE3_AVAILABLE:
        raise TreeError("ete3 required for branch length analysis")
    
    try:
        tree = Tree(str(tree_path))
        
        # Get all terminal branch lengths
        branch_lengths = []
        for leaf in tree.get_leaves():
            dist = leaf.dist
            branch_lengths.append((leaf.name, dist))
        
        if not branch_lengths:
            return []
        
        # Calculate threshold
        lengths = [bl[1] for bl in branch_lengths]
        mean_length = sum(lengths) / len(lengths)
        threshold_value = mean_length * threshold
        
        # Find long branches
        long_branches = [(name, length) for name, length in branch_lengths 
                        if length > threshold_value]
        long_branches.sort(key=lambda x: x[1], reverse=True)
        
        return long_branches
        
    except Exception as e:
        raise TreeError(f"Failed to analyze branch lengths: {e}")
