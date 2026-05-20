#!/usr/bin/env python3
"""
Generate synthetic test VCF files for SNPhylo2 testing.
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np


def generate_test_vcf(
    output_path: str,
    n_samples: int = 10,
    n_snps: int = 100,
    seed: int = 42,
    missing_rate: float = 0.05,
    maf_range: tuple = (0.1, 0.4),
):
    """
    Generate a synthetic VCF file for testing.
    
    Args:
        output_path: Output VCF file path
        n_samples: Number of samples
        n_snps: Number of SNPs
        seed: Random seed
        missing_rate: Rate of missing genotypes
        maf_range: Range of minor allele frequencies
    """
    random.seed(seed)
    np.random.seed(seed)
    
    output = Path(output_path)
    
    # Generate sample names
    samples = [f"Sample_{i:03d}" for i in range(1, n_samples + 1)]
    
    # Open output
    import gzip
    opener = gzip.open if str(output).endswith('.gz') else open
    
    with opener(output, 'wt') as f:
        # Write header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=snphylo2_test_generator\n")
        f.write(f"##random_seed={seed}\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
        
        # Write column header
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write("\t".join(samples) + "\n")
        
        # Generate SNPs
        bases = ['A', 'C', 'G', 'T']
        
        for i in range(n_snps):
            chrom = f"chr{(i % 10) + 1}"
            pos = (i // 10) * 1000 + 100
            snp_id = f"rs{i+1}"
            
            # Choose alleles
            ref, alt = random.sample(bases, 2)
            
            # Generate MAF
            maf = random.uniform(*maf_range)
            
            # Generate genotypes
            genotypes = []
            for _ in range(n_samples):
                if random.random() < missing_rate:
                    genotypes.append("./.")
                else:
                    # Generate genotype based on MAF
                    p_alt = maf
                    n_alleles = np.random.binomial(2, p_alt)
                    
                    if n_alleles == 0:
                        gt = "0/0"
                    elif n_alleles == 1:
                        gt = "0/1"
                    else:
                        gt = "1/1"
                    
                    # Add depth and quality
                    dp = random.randint(5, 50)
                    gq = random.randint(20, 99)
                    genotypes.append(f"{gt}:{dp}:{gq}")
            
            # Write variant line
            qual = random.randint(20, 100)
            info = f"DP={sum(int(gt.split(':')[1]) for gt in genotypes if ':' in gt)}"
            format_str = "GT:DP:GQ"
            
            f.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t{qual}\tPASS\t{info}\t{format_str}\t")
            f.write("\t".join(genotypes) + "\n")
    
    print(f"Generated test VCF: {output}")
    print(f"  Samples: {n_samples}")
    print(f"  SNPs: {n_snps}")
    print(f"  Missing rate: {missing_rate}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic test VCF files"
    )
    parser.add_argument(
        "output",
        help="Output VCF file path (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "--samples", "-s",
        type=int,
        default=10,
        help="Number of samples (default: 10)"
    )
    parser.add_argument(
        "--snps", "-n",
        type=int,
        default=100,
        help="Number of SNPs (default: 100)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--missing-rate",
        type=float,
        default=0.05,
        help="Missing genotype rate (default: 0.05)"
    )
    
    args = parser.parse_args()
    
    generate_test_vcf(
        args.output,
        n_samples=args.samples,
        n_snps=args.snps,
        seed=args.seed,
        missing_rate=args.missing_rate,
    )


if __name__ == "__main__":
    main()
