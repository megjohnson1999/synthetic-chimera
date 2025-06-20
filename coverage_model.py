"""
Coverage modeling for coassembly-like simulations.

Generates realistic uneven coverage patterns across multiple samples.
"""

import random
import numpy as np
import logging
from typing import List, Dict, Tuple
from dataclasses import dataclass

from contig_generators import ContigInfo


@dataclass
class SampleCoverage:
    """Coverage information for a single sample."""
    sample_id: str
    base_coverage: float
    coverage_multiplier: float = 1.0
    

@dataclass
class ContigCoverage:
    """Coverage information for a contig across all samples."""
    contig_id: str
    sample_coverages: Dict[str, float]
    total_coverage: float
    average_coverage: float
    coverage_std: float


class CoverageModel:
    """
    Generates coassembly-like coverage patterns.
    
    Simulates multiple samples with different coverage profiles,
    creating realistic uneven coverage independent of chimera status.
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Generate sample profiles
        self.samples = self._generate_sample_profiles()
        
    def _generate_sample_profiles(self) -> List[SampleCoverage]:
        """Generate coverage profiles for multiple samples."""
        samples = []
        
        for i in range(self.config.coverage.samples):
            sample_id = f"sample_{i+1:02d}"
            
            # Generate base coverage using log-normal distribution
            if self.config.coverage.distribution == "log_normal":
                # Log-normal parameters to achieve desired range
                min_cov, max_cov = self.config.coverage.base_coverage_range
                
                # Use log-normal distribution
                mu = np.log(np.sqrt(min_cov * max_cov))  # Geometric mean
                sigma = (np.log(max_cov) - np.log(min_cov)) / 4  # Span ~4 standard deviations
                
                base_coverage = np.random.lognormal(mu, sigma)
                
                # Clamp to range
                base_coverage = max(min_cov, min(max_cov, base_coverage))
                
            else:
                # Uniform distribution fallback
                min_cov, max_cov = self.config.coverage.base_coverage_range
                base_coverage = random.uniform(min_cov, max_cov)
            
            samples.append(SampleCoverage(
                sample_id=sample_id,
                base_coverage=base_coverage
            ))
        
        self.logger.info(f"Generated {len(samples)} sample profiles")
        self.logger.info(f"Coverage range: {min(s.base_coverage for s in samples):.2f}x - "
                        f"{max(s.base_coverage for s in samples):.2f}x")
        
        return samples
    
    def generate_contig_coverage(self, contigs: List[ContigInfo]) -> List[ContigCoverage]:
        """
        Generate coverage for all contigs across all samples.
        
        Args:
            contigs: List of generated contigs
            
        Returns:
            List of ContigCoverage objects
        """
        contig_coverages = []
        
        for contig in contigs:
            sample_coverages = {}
            
            for sample in self.samples:
                # Generate per-contig coverage variation
                # Independent of chimera status (realistic for coassemblies)
                coverage_variation = self._generate_coverage_variation(contig, sample)
                
                final_coverage = sample.base_coverage * coverage_variation
                sample_coverages[sample.sample_id] = final_coverage
            
            # Calculate statistics
            coverage_values = list(sample_coverages.values())
            total_coverage = sum(coverage_values)
            average_coverage = np.mean(coverage_values)
            coverage_std = np.std(coverage_values)
            
            contig_coverage = ContigCoverage(
                contig_id=contig.id,
                sample_coverages=sample_coverages,
                total_coverage=total_coverage,
                average_coverage=average_coverage,
                coverage_std=coverage_std
            )
            
            contig_coverages.append(contig_coverage)
        
        self._log_coverage_statistics(contig_coverages)
        return contig_coverages
    
    def _generate_coverage_variation(self, contig: ContigInfo, sample: SampleCoverage) -> float:
        """
        Generate per-contig coverage variation for a sample.
        
        This simulates biological and technical variation in coverage
        that's independent of chimera status.
        """
        # Base variation factor (log-normal around 1.0)
        base_variation = np.random.lognormal(0, 0.5)  # ~2-fold variation typical
        
        # Length-dependent bias (shorter contigs often have more variable coverage)
        length_bias = 1.0
        if contig.length < 1000:
            length_bias = np.random.lognormal(0, 0.3)  # Additional variation for short contigs
        
        # Sample-specific bias (some samples may be systematically different)
        sample_bias = sample.coverage_multiplier
        
        # Combine factors
        total_variation = base_variation * length_bias * sample_bias
        
        # Prevent extremely high or low coverage
        total_variation = max(0.01, min(100, total_variation))
        
        return total_variation
    
    def _log_coverage_statistics(self, contig_coverages: List[ContigCoverage]):
        """Log coverage statistics for monitoring."""
        total_coverages = [cc.total_coverage for cc in contig_coverages]
        avg_coverages = [cc.average_coverage for cc in contig_coverages]
        
        self.logger.info(f"Coverage statistics:")
        self.logger.info(f"  Total coverage range: {min(total_coverages):.2f}x - {max(total_coverages):.2f}x")
        self.logger.info(f"  Average coverage range: {min(avg_coverages):.2f}x - {max(avg_coverages):.2f}x")
        self.logger.info(f"  Mean total coverage: {np.mean(total_coverages):.2f}x")
        
        # Check for realistic distribution
        low_coverage_contigs = sum(1 for cc in contig_coverages if cc.total_coverage < 10)
        high_coverage_contigs = sum(1 for cc in contig_coverages if cc.total_coverage > 500)
        
        self.logger.info(f"  Low coverage contigs (<10x): {low_coverage_contigs} ({100*low_coverage_contigs/len(contig_coverages):.1f}%)")
        self.logger.info(f"  High coverage contigs (>500x): {high_coverage_contigs} ({100*high_coverage_contigs/len(contig_coverages):.1f}%)")
    
    def calculate_read_counts(self, contig_coverages: List[ContigCoverage], 
                            contigs: List[ContigInfo]) -> Dict[str, int]:
        """
        Calculate number of read pairs to generate for each contig.
        
        Args:
            contig_coverages: Coverage information
            contigs: Contig information
            
        Returns:
            Dictionary mapping contig_id to number of read pairs
        """
        read_counts = {}
        read_length = self.config.read_simulation.read_length
        
        for contig_coverage, contig in zip(contig_coverages, contigs):
            # Calculate reads needed for total coverage
            total_bases_needed = contig_coverage.total_coverage * contig.length
            
            # Convert to read pairs (each pair contributes 2 * read_length bases)
            read_pairs_needed = int(total_bases_needed / (2 * read_length))
            
            # Ensure minimum read pairs for very low coverage
            read_pairs_needed = max(1, read_pairs_needed)
            
            read_counts[contig.id] = read_pairs_needed
        
        total_read_pairs = sum(read_counts.values())
        self.logger.info(f"Total read pairs to generate: {total_read_pairs:,}")
        
        return read_counts
    
    def save_coverage_info(self, contig_coverages: List[ContigCoverage], 
                          output_file: str):
        """Save coverage information to file."""
        with open(output_file, 'w') as f:
            # Header
            sample_headers = [f"{s.sample_id}_coverage" for s in self.samples]
            header = ["contig_id", "total_coverage", "average_coverage", "coverage_std"] + sample_headers
            f.write('\t'.join(header) + '\n')
            
            # Data
            for cc in contig_coverages:
                row = [
                    cc.contig_id,
                    f"{cc.total_coverage:.2f}",
                    f"{cc.average_coverage:.2f}",
                    f"{cc.coverage_std:.2f}"
                ]
                
                # Add per-sample coverages
                for sample in self.samples:
                    coverage = cc.sample_coverages.get(sample.sample_id, 0.0)
                    row.append(f"{coverage:.2f}")
                
                f.write('\t'.join(row) + '\n')
        
        self.logger.info(f"Coverage information saved to {output_file}")


def apply_coverage_model(contigs: List[ContigInfo], config) -> Tuple[List[ContigCoverage], Dict[str, int]]:
    """
    Apply coverage model to contigs and return coverage info and read counts.
    
    Args:
        contigs: List of generated contigs
        config: Configuration object
        
    Returns:
        Tuple of (contig_coverages, read_counts)
    """
    logger = logging.getLogger(__name__)
    logger.info("Applying coassembly coverage model...")
    
    # Initialize coverage model
    coverage_model = CoverageModel(config)
    
    # Generate coverage for all contigs
    contig_coverages = coverage_model.generate_contig_coverage(contigs)
    
    # Calculate read counts needed
    read_counts = coverage_model.calculate_read_counts(contig_coverages, contigs)
    
    return contig_coverages, read_counts