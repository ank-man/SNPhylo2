"""
Unit tests for configuration module.
"""

import tempfile
from pathlib import Path

import pytest

from snphylo2.config import (
    SNPhylo2Config,
    InputConfig,
    FilteringConfig,
    TreeConfig,
    create_default_config,
    load_config,
    save_config,
)
from snphylo2.exceptions import ConfigurationError


class TestConfigValidation:
    """Test configuration validation."""
    
    def test_maf_validation(self):
        """Test that MAF min < max is enforced."""
        with pytest.raises(ConfigurationError):
            SNPhylo2Config(
                input=InputConfig(path=Path("test.vcf")),
                filtering=FilteringConfig(
                    maf={"min": 0.5, "max": 0.3}  # Invalid: min > max
                ),
            )
    
    def test_ld_pruning_window_validation(self):
        """Test that window > step is enforced."""
        with pytest.raises(ConfigurationError):
            SNPhylo2Config(
                input=InputConfig(path=Path("test.vcf")),
                ld_pruning={"window_size": 10, "step_size": 20},  # Invalid: window < step
            )


class TestConfigIO:
    """Test configuration file I/O."""
    
    def test_save_and_load_config(self, tmp_path):
        """Test saving and loading configuration."""
        config = SNPhylo2Config(
            input=InputConfig(path=Path("test.vcf")),
            output={"prefix": "test_output"},
        )
        
        config_path = tmp_path / "test_config.yaml"
        save_config(config, config_path)
        
        assert config_path.exists()
        
        loaded = load_config(config_path)
        assert loaded.output.prefix == "test_output"
    
    def test_load_nonexistent_config(self):
        """Test loading non-existent config file."""
        with pytest.raises(ConfigurationError):
            load_config("/nonexistent/config.yaml")


class TestDefaultConfigs:
    """Test default configuration presets."""
    
    def test_default_preset(self, tmp_path):
        """Test default preset."""
        config_path = tmp_path / "default.yaml"
        config = create_default_config(config_path, preset="default")
        
        assert config.filtering.maf.min == 0.05
    
    def test_human_preset(self, tmp_path):
        """Test human preset."""
        config_path = tmp_path / "human.yaml"
        config = create_default_config(config_path, preset="human")
        
        assert config.filtering.maf.min == 0.01
        assert "chrM" in config.filtering.exclude_chromosomes
    
    def test_plant_preset(self, tmp_path):
        """Test plant preset."""
        config_path = tmp_path / "plant.yaml"
        config = create_default_config(config_path, preset="plant")
        
        assert config.filtering.depth.min == 3
        assert config.ld_pruning.window_size == 100
    
    def test_invalid_preset(self, tmp_path):
        """Test invalid preset raises error."""
        with pytest.raises(ConfigurationError):
            create_default_config(tmp_path / "invalid.yaml", preset="invalid")
