#!/usr/bin/env python3
"""
Generate simulated SNP data with known phylogeny for validation.
Uses msprime for coalescent simulations.
"""

import argparse
import tempfile
from pathlib import Path
from typing import List, Tuple

import numpy as np


def simulate_phylogeny(
    n_samples: int = 10,
    sequence_length: int = 10000,
    mutation_rate: float = 1e-8,
    recombination_rate: float = 1e-8,
    effective_pop_size: float = 1e4,
    seed: int = 42,
) -> Tuple['tskit.TreeSequence', str]:
    """
    Simulate SNP data with known tree using msprime.
    
    Returns:
        Tree sequence and Newick string of true tree
    """
    try:
        import msprime
    except ImportError:
        raise ImportError("msprime is required for simulations. Install with: pip install msprime")
    
    # Simulate ancestry
    ts = msprime.sim_ancestry(
        samples=n_samples,
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        population_size=effective_pop_size,
        random_seed=seed,
    )
    
    # Add mutations
    ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        random_seed=seed,
    )
    
    # Extract true tree (first tree in sequence)
    true_tree = ts.first().newick()
    
    return ts, true_tree


def ts_to_vcf(
    ts: 'tskit.TreeSequence',
    output_path: str,
    sample_names: List[str] = None,
):
    """
    Write tree sequence to VCF file.
    
    Args:
        ts: Tree sequence from msprime
        output_path: Output VCF file path
        sample_names: Optional list of sample names
    """
    import gzip
    
    if sample_names is None:
        sample_names = [f"Sample_{i}" for i in range(ts.num_samples)]
    
    opener = gzip.open if str(output_path).endswith('.gz') else open
    
    with opener(output_path, 'wt') as f:
        # Write header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=msprime_simulation\n")
        f.write("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write("\t".join(sample_names) + "\n")
        
        # Write variants
        for variant in ts.variants():
            pos = int(variant.site.position)
            alleles = variant.alleles
            
            # Skip if more than 2 alleles (not biallelic)
            if len(alleles) > 2:
                continue
            
            ref = alleles[0]
            alt = alleles[1] if len(alleles) > 1 else "."
            
            # Convert genotypes to VCF format
            genotypes = []
            for g in variant.genotypes:
                a1, a2 = g
                if a1 == -1 or a2 == -1:
                    genotypes.append("./.")
                else:
                    genotypes.append(f"{a1}/{a2}")
            
            # Write variant
            info = f"AN={len(sample_names) * 2}"
            f.write(f"1\t{pos}\t.\t{ref}\t{alt}\t30\tPASS\t{info}\tGT\t")
            f.write("\t".join(genotypes) + "\n")


def generate_benchmark_dataset(
    output_dir: str,
    n_samples: int,
    n_snps: int,
    seed: int = 42,
) -> dict:
    """
    Generate a benchmark dataset with known phylogeny.
    
    Args:
        output_dir: Output directory
        n_samples: Number of samples
        n_snps: Approximate number of SNPs
        seed: Random seed
        
    Returns:
        Dictionary with paths to generated files and metrics
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate sequence length to get approximately n_snps
    mutation_rate = 1e-8
    sequence_length = int(n_snps / (mutation_rate * 10000))  # Rough estimate
    
    # Simulate
    ts, true_tree = simulate_phylogeny(
        n_samples=n_samples,
        sequence_length=sequence_length,
        mutation_rate=mutation_rate,
        seed=seed,
    )
    
    # Write VCF
    vcf_path = output_dir / f"simulated_{n_samples}s_{n_snps}snps.vcf.gz"
    sample_names = [f"Sample_{i:03d}" for i in range(n_samples)]
    ts_to_vcf(ts, str(vcf_path), sample_names)
    
    # Write true tree
    tree_path = output_dir / f"simulated_{n_samples}s_{n_snps}snps_true_tree.nwk"
    with open(tree_path, 'w') as f:
        f.write(true_tree)
    
    # Generate metadata
    metadata_path = output_dir / f"simulated_{n_samples}s_{n_snps}snps_metadata.tsv"
    with open(metadata_path, 'w') as f:
        f.write("sample_id\tpopulation\tlatitude\tlongitude\n")
        for i, sample in enumerate(sample_names):
            pop = f"Pop{(i % 3) + 1}"
            lat = 30.0 + (i % 10)
            lon = -90.0 + (i % 20)
            f.write(f"{sample}\t{pop}\t{lat}\t{lon}\n")
    
    # Count actual SNPs
    actual_snps = sum(1 for _ in ts.variants() if len(_.alleles) <= 2)
    
    return {
        'vcf': str(vcf_path),
        'true_tree': str(tree_path),
        'metadata': str(metadata_path),
        'n_samples': n_samples,
        'n_snps_requested': n_snps,
        'n_snps_actual': actual_snps,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate simulated SNP datasets with known phylogeny"
    )
    parser.add_argument(
        "output_dir",
        help="Output directory"
    )
    parser.add_argument(
        "--samples", "-s",
        type=int,
        default=10,
        help="Number of samples"
    )
    parser.add_argument(
        "--snps", "-n",
        type=int,
        default=1000,
        help="Approximate number of SNPs"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed"
    )
    
    args = parser.parse_args()
    
    result = generate_benchmark_dataset(
        args.output_dir,
        n_samples=args.samples,
        n_snps=args.snps,
        seed=args.seed,
    )
    
    print(f"Generated benchmark dataset:")
    print(f"  VCF: {result['vcf']}")
    print(f"  True tree: {result['true_tree']}")
    print(f"  Metadata: {result['metadata']}")
    print(f"  Samples: {result['n_samples']}")
    print(f"  SNPs: {result['n_snps_actual']} (requested {result['n_snps_requested']})")


if __name__ == "__main__":
    main()
