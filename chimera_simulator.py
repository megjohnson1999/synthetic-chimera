#!/usr/bin/env python3
"""
Metagenomic Chimera Simulator

Generates synthetic metagenomic assemblies containing chimeric contigs with
realistic paired-end mapping artifacts for testing chimera detection algorithms.
"""

import argparse
import sys
import os
import random
from pathlib import Path
from datetime import datetime
import logging

from config import ChimeraConfig
from utils import create_output_structure, setup_logging, validate_fasta
from contig_generators import GenomePool, generate_contigs
from coverage_model import apply_coverage_model
from read_simulator import simulate_reads_from_contigs
from mapping import map_reads_to_contigs
from ground_truth import generate_ground_truth_files


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate synthetic metagenomic assemblies with chimeric contigs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  python chimera_simulator.py --genomes viral_genomes.fasta
  
  # Custom configuration
  python chimera_simulator.py --genomes viral_genomes.fasta --config custom_params.yaml --output experiment_1/
  
  # Pipeline integration
  for dataset in datasets/*.fasta; do
      python chimera_simulator.py --genomes "$dataset" --output "results/$(basename "$dataset" .fasta)"
  done
        """
    )
    
    parser.add_argument(
        "--genomes",
        required=True,
        type=Path,
        help="Path to input viral genomes FASTA file"
    )
    
    parser.add_argument(
        "--config",
        type=Path,
        help="YAML configuration file (falls back to defaults)"
    )
    
    parser.add_argument(
        "--output",
        type=Path,
        help="Output directory (falls back to timestamped directory)"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    return parser.parse_args()


def validate_inputs(args):
    """Validate input arguments and files."""
    if not args.genomes.exists():
        raise FileNotFoundError(f"Input genomes file not found: {args.genomes}")
    
    if not args.genomes.suffix.lower() in ['.fasta', '.fa', '.fas']:
        raise ValueError(f"Input file must be FASTA format: {args.genomes}")
    
    if args.config and not args.config.exists():
        raise FileNotFoundError(f"Configuration file not found: {args.config}")
    
    # Set default output directory if not provided
    if not args.output:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        genome_name = args.genomes.stem
        args.output = Path(f"chimera_simulation_{genome_name}_{timestamp}")
    
    return args


def main():
    try:
        # Parse and validate arguments
        args = parse_arguments()
        args = validate_inputs(args)
        
        # Create output directory structure
        output_dirs = create_output_structure(args.output)
        
        # Setup logging
        setup_logging(output_dirs['reports'] / 'simulation.log', args.verbose)
        logger = logging.getLogger(__name__)
        
        logger.info(f"Starting chimera simulation")
        logger.info(f"Input genomes: {args.genomes}")
        logger.info(f"Output directory: {args.output}")
        logger.info(f"Random seed: {args.seed}")
        
        # Load configuration
        config = ChimeraConfig(args.config, args.seed)
        config.validate()
        logger.info(f"Configuration loaded: {config.source}")
        
        # Save configuration to output directory
        config.save_to_file(output_dirs['config'] / 'used_config.yaml')
        
        # Validate input genomes
        genome_stats = validate_fasta(args.genomes)
        logger.info(f"Input validation: {genome_stats['num_sequences']} sequences, "
                   f"avg length: {genome_stats['avg_length']:.0f} bp")
        
        # Set random seed
        random.seed(args.seed)
        
        # Phase 2: Generate contigs
        logger.info("Phase 2: Generating contigs...")
        genome_pool = GenomePool(args.genomes)
        contigs = generate_contigs(genome_pool, config)
        
        # Save contigs to FASTA
        contigs_file = output_dirs['contigs'] / 'synthetic_assembly.fasta'
        with open(contigs_file, 'w') as f:
            for contig in contigs:
                f.write(f">{contig.id}\n{contig.sequence}\n")
        
        logger.info(f"Generated {len(contigs)} contigs saved to {contigs_file}")
        
        # Phase 3: Coverage modeling and read simulation
        logger.info("Phase 3: Applying coverage model and simulating reads...")
        contig_coverages, read_counts = apply_coverage_model(contigs, config)
        
        # Save coverage information
        coverage_file = output_dirs['reports'] / 'coverage_info.tsv'
        from coverage_model import CoverageModel
        coverage_model = CoverageModel(config)
        coverage_model.save_coverage_info(contig_coverages, coverage_file)
        
        # Simulate reads
        r1_file, r2_file = simulate_reads_from_contigs(
            contigs, contig_coverages, read_counts, config, output_dirs['reads']
        )
        
        logger.info(f"Read simulation completed:")
        logger.info(f"  R1: {r1_file}")
        logger.info(f"  R2: {r2_file}")
        
        # Phase 4: Mapping & validation
        logger.info("Phase 4: Mapping reads and analyzing PE artifacts...")
        bam_file, pe_analysis = map_reads_to_contigs(
            contigs_file, r1_file, r2_file, contigs, config, output_dirs['mappings']
        )
        
        logger.info(f"Mapping completed: {bam_file}")
        
        # Phase 5: Ground truth generation
        logger.info("Phase 5: Generating ground truth data...")
        ground_truth_files = generate_ground_truth_files(
            contigs, contig_coverages, pe_analysis, config, output_dirs['ground_truth']
        )
        
        logger.info("Ground truth generation completed")
        
        # Generate final summary
        logger.info("Simulation pipeline completed successfully!")
        logger.info("Generated files:")
        logger.info(f"  Synthetic assembly: {contigs_file}")
        logger.info(f"  Reads: {r1_file}, {r2_file}")
        logger.info(f"  Mappings: {bam_file}")
        logger.info(f"  Ground truth: {len(ground_truth_files)} files")
        logger.info(f"  Reports: {coverage_file} + mapping stats")
        
        logger.info("Chimera simulation completed successfully")
        print(f"Results saved to: {args.output}")
        
    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()