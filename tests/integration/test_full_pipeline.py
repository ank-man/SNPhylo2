"""
Integration tests for the complete SNPhylo2 pipeline.
"""

import tempfile
from pathlib import Path

import pytest

from snphylo2.config import SNPhylo2Config, InputConfig, OutputConfig, ComputeConfig
from snphylo2.pipeline import Pipeline


@pytest.mark.integration
@pytest.mark.slow
class TestFullPipeline:
    """Test complete pipeline execution."""
    
    def test_simple_pipeline(self, tmp_path):
        """Test simple end-to-end pipeline execution."""
        # This test requires test data
        # For now, it's a placeholder that shows the structure
        
        # Generate test data
        import subprocess
        test_vcf = tmp_path / "test.vcf.gz"
        
        # Generate test VCF
        script = Path(__file__).parent.parent / "fixtures" / "generate_test_vcf.py"
        subprocess.run([
            "python", str(script), str(test_vcf),
            "--samples", "5",
            "--snps", "50",
            "--seed", "42",
        ], check=True)
        
        # Create config
        config = SNPhylo2Config(
            input=InputConfig(path=test_vcf),
            output=OutputConfig(
                directory=tmp_path / "results",
                prefix="test",
                keep_intermediates=True,
            ),
            compute=ComputeConfig(threads=1),
            tree={"bootstrap": {"replicates": 100}},  # Reduced for speed
        )
        
        # Run pipeline
        pipeline = Pipeline(config)
        results = pipeline.run()
        
        # Verify outputs
        assert results['tree_file'] is not None
        assert Path(results['tree_file']).exists()
        assert results['report_path'] is not None
        assert Path(results['report_path']).exists()
        
        # Check stats
        assert results['stats']['input_samples'] == 5
        assert results['stats']['input_snps'] == 50
    
    def test_pipeline_with_metadata(self, tmp_path):
        """Test pipeline with metadata file."""
        # Generate test data with metadata
        test_vcf = tmp_path / "test.vcf.gz"
        test_meta = tmp_path / "metadata.tsv"
        
        # Generate VCF
        script = Path(__file__).parent.parent / "fixtures" / "generate_test_vcf.py"
        subprocess.run([
            "python", str(script), str(test_vcf),
            "--samples", "5",
            "--snps", "50",
        ], check=True)
        
        # Generate metadata
        with open(test_meta, 'w') as f:
            f.write("sample_id\tpopulation\n")
            for i in range(5):
                f.write(f"Sample_{i+1:03d}\tPop{(i % 2) + 1}\n")
        
        config = SNPhylo2Config(
            input=InputConfig(path=test_vcf),
            metadata={"path": test_meta, "sample_id_column": "sample_id"},
            output=OutputConfig(directory=tmp_path / "results", prefix="test"),
            compute=ComputeConfig(threads=1),
        )
        
        pipeline = Pipeline(config)
        results = pipeline.run()
        
        assert results['tree_file'] is not None


@pytest.mark.integration
class TestWithSimulatedData:
    """Tests using simulated data with known phylogeny."""
    
    @pytest.fixture(scope="class")
    def simulated_dataset(self, tmp_path_factory):
        """Generate simulated dataset once for all tests."""
        tmp_path = tmp_path_factory.mktemp("simulated")
        
        try:
            from tests.fixtures.generate_simulated_data import generate_benchmark_dataset
            
            dataset = generate_benchmark_dataset(
                str(tmp_path),
                n_samples=10,
                n_snps=100,
                seed=42,
            )
            return dataset
        except ImportError:
            pytest.skip("msprime not available for simulation")
    
    def test_accuracy_on_simulated_data(self, simulated_dataset, tmp_path):
        """Test that inferred tree approximates true tree."""
        vcf_path = Path(simulated_dataset['vcf'])
        true_tree_path = Path(simulated_dataset['true_tree'])
        
        # Run pipeline
        config = SNPhylo2Config(
            input=InputConfig(path=vcf_path),
            output=OutputConfig(directory=tmp_path / "results", prefix="sim"),
            compute=ComputeConfig(threads=2),
            tree={"bootstrap": {"replicates": 100}},
        )
        
        pipeline = Pipeline(config)
        results = pipeline.run()
        
        # Compare trees
        from snphylo2.tree.tree_comparison import compare_trees
        
        comparison = compare_trees(
            Path(results['tree_file']),
            true_tree_path,
        )
        
        # Check that trees are similar (normalized RF < 0.5)
        assert comparison['normalized_rf'] < 0.5, \
            f"Inferred tree differs too much from true tree (RF={comparison['normalized_rf']:.3f})"
