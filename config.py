"""
Configuration management for chimera simulator.
"""

import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, Any, Optional
import logging


@dataclass
class AssemblyConfig:
    """Configuration for synthetic assembly generation."""
    total_contigs: int = 1000
    chimera_percentage: float = 0.15
    normal_contig_length_range: tuple = (500, 5000)
    chimeric_contig_length_range: tuple = (1000, 8000)


@dataclass
class ChimeraTypesConfig:
    """Configuration for different chimera types."""
    gap_mediated: float = 0.4      # Clear insert size jumps
    orientation_flip: float = 0.3   # Improper pair orientations
    distant_join: float = 0.3       # Extreme insert sizes
    
    def __post_init__(self):
        total = self.gap_mediated + self.orientation_flip + self.distant_join
        if abs(total - 1.0) > 0.001:
            raise ValueError(f"Chimera type proportions must sum to 1.0, got {total}")


@dataclass
class CoverageConfig:
    """Configuration for coassembly-like coverage simulation."""
    samples: int = 5
    base_coverage_range: tuple = (0.1, 500.0)
    distribution: str = "log_normal"
    coverage_correlation: float = 0.0  # Independence from chimera status


@dataclass
class ReadSimulationConfig:
    """Configuration for read simulation."""
    coverage_target: float = 100.0
    read_length: int = 150
    insert_size_mean: int = 300
    insert_size_std: int = 50
    error_rate: float = 0.01
    quality_profile: str = "illumina"


@dataclass
class MappingConfig:
    """Configuration for read mapping."""
    aligner: str = "bwa_mem"  # or "minimap2"
    threads: int = 4
    min_mapping_quality: int = 10
    

@dataclass
class ValidationConfig:
    """Configuration for PE artifact validation."""
    min_mapping_rate: float = 0.85
    insert_size_outlier_threshold: float = 3.0  # Standard deviations
    orientation_artifact_threshold: float = 0.05  # Fraction of improper pairs


class ChimeraConfig:
    """Main configuration class for chimera simulator."""
    
    def __init__(self, config_path: Optional[Path] = None, seed: int = 42):
        self.seed = seed
        self.source = "defaults"
        
        # Load default configuration
        self._load_defaults()
        
        # Override with file configuration if provided
        if config_path:
            self._load_from_file(config_path)
            self.source = str(config_path)
    
    def _load_defaults(self):
        """Load default configuration values."""
        self.assembly = AssemblyConfig()
        self.chimera_types = ChimeraTypesConfig()
        self.coverage = CoverageConfig()
        self.read_simulation = ReadSimulationConfig()
        self.mapping = MappingConfig()
        self.validation = ValidationConfig()
    
    def _load_from_file(self, config_path: Path):
        """Load configuration from YAML file."""
        logger = logging.getLogger(__name__)
        
        try:
            with open(config_path, 'r') as f:
                config_data = yaml.safe_load(f)
            
            # Update configuration sections
            if 'assembly' in config_data:
                self._update_config(self.assembly, config_data['assembly'])
            
            if 'chimera_types' in config_data:
                self._update_config(self.chimera_types, config_data['chimera_types'])
            
            if 'coverage' in config_data:
                self._update_config(self.coverage, config_data['coverage'])
            
            if 'read_simulation' in config_data:
                self._update_config(self.read_simulation, config_data['read_simulation'])
            
            if 'mapping' in config_data:
                self._update_config(self.mapping, config_data['mapping'])
            
            if 'validation' in config_data:
                self._update_config(self.validation, config_data['validation'])
            
            logger.info(f"Configuration loaded from {config_path}")
            
        except Exception as e:
            logger.error(f"Failed to load configuration from {config_path}: {e}")
            raise
    
    def _update_config(self, config_obj, updates: Dict[str, Any]):
        """Update configuration object with new values."""
        for key, value in updates.items():
            if hasattr(config_obj, key):
                setattr(config_obj, key, value)
            else:
                logging.getLogger(__name__).warning(
                    f"Unknown configuration parameter: {key}"
                )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary for serialization."""
        return {
            'assembly': self.assembly.__dict__,
            'chimera_types': self.chimera_types.__dict__,
            'coverage': self.coverage.__dict__,
            'read_simulation': self.read_simulation.__dict__,
            'mapping': self.mapping.__dict__,
            'validation': self.validation.__dict__,
            'seed': self.seed,
            'source': self.source
        }
    
    def save_to_file(self, output_path: Path):
        """Save current configuration to YAML file."""
        with open(output_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, indent=2)
    
    def validate(self):
        """Validate configuration parameters."""
        # Validate chimera type proportions
        self.chimera_types.__post_init__()
        
        # Validate coverage parameters
        if self.coverage.base_coverage_range[0] >= self.coverage.base_coverage_range[1]:
            raise ValueError("Coverage range minimum must be less than maximum")
        
        # Validate read simulation parameters
        if self.read_simulation.insert_size_mean <= 0:
            raise ValueError("Insert size mean must be positive")
        
        if self.read_simulation.insert_size_std <= 0:
            raise ValueError("Insert size standard deviation must be positive")
        
        # Validate mapping parameters
        if self.mapping.aligner not in ['bwa_mem', 'minimap2']:
            raise ValueError(f"Unsupported aligner: {self.mapping.aligner}")
        
        logging.getLogger(__name__).info("Configuration validation passed")