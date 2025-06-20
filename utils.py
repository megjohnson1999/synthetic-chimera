"""
Utility functions for chimera simulator.
"""

import logging
import sys
from pathlib import Path
from typing import Dict
from datetime import datetime


def create_output_structure(output_dir: Path) -> Dict[str, Path]:
    """
    Create organized output directory structure.
    
    Returns:
        Dictionary mapping directory names to Path objects
    """
    output_dir = Path(output_dir)
    
    directories = {
        'root': output_dir,
        'contigs': output_dir / 'contigs',
        'reads': output_dir / 'reads', 
        'mappings': output_dir / 'mappings',
        'ground_truth': output_dir / 'ground_truth',
        'reports': output_dir / 'reports',
        'config': output_dir / 'config'
    }
    
    # Create all directories
    for dir_path in directories.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    return directories


def setup_logging(log_file: Path, verbose: bool = False):
    """
    Configure logging for the application.
    
    Args:
        log_file: Path to log file
        verbose: Enable debug level logging
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_formatter = logging.Formatter(
        '%(levelname)s: %(message)s'
    )
    
    # Setup file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(file_formatter)
    
    # Setup console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(console_formatter)
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        handlers=[file_handler, console_handler]
    )


def parse_fasta(fasta_path: Path):
    """
    Simple FASTA parser that yields (header, sequence) tuples.
    
    Args:
        fasta_path: Path to FASTA file
        
    Yields:
        Tuple of (header, sequence)
    """
    with open(fasta_path, 'r') as f:
        header = None
        sequence = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Yield previous record if exists
                if header and sequence:
                    yield header, ''.join(sequence)
                
                # Start new record
                header = line[1:]  # Remove '>'
                sequence = []
            else:
                sequence.append(line)
        
        # Yield final record
        if header and sequence:
            yield header, ''.join(sequence)


def validate_fasta(fasta_path: Path) -> Dict[str, int]:
    """
    Validate FASTA file and return statistics.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary with validation statistics
    """
    stats = {
        'num_sequences': 0,
        'total_length': 0,
        'min_length': float('inf'),
        'max_length': 0,
        'avg_length': 0
    }
    
    lengths = []
    
    for header, sequence in parse_fasta(fasta_path):
        # Validate sequence characters
        valid_chars = set('ATGCNatgcn')
        if not set(sequence).issubset(valid_chars):
            invalid_chars = set(sequence) - valid_chars
            raise ValueError(f"Invalid characters in sequence '{header}': {invalid_chars}")
        
        length = len(sequence)
        lengths.append(length)
        
        stats['num_sequences'] += 1
        stats['total_length'] += length
        stats['min_length'] = min(stats['min_length'], length)
        stats['max_length'] = max(stats['max_length'], length)
    
    if stats['num_sequences'] == 0:
        raise ValueError("No sequences found in FASTA file")
    
    stats['avg_length'] = stats['total_length'] / stats['num_sequences']
    
    # Reset min_length if no sequences found
    if stats['min_length'] == float('inf'):
        stats['min_length'] = 0
    
    return stats


def format_size(size_bytes: int) -> str:
    """Format file size in human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def get_timestamp() -> str:
    """Get current timestamp string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def reverse_complement(sequence: str) -> str:
    """
    Return reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n'
    }
    
    return ''.join(complement.get(base, base) for base in reversed(sequence))