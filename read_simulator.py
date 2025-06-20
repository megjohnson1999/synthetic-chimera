"""
Read simulation module for chimera simulator.

Generates realistic paired-end reads from synthetic contigs with proper
insert size distributions and error profiles.
"""

import subprocess
import random
import logging
import tempfile
import gzip
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from contig_generators import ContigInfo
from coverage_model import ContigCoverage


class ReadSimulator:
    """
    Simulates paired-end reads from synthetic contigs.
    
    Uses external tools (ART or BBTools) for realistic error profiles
    and insert size distributions.
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Check available tools
        self.simulator_tool = self._detect_simulator_tool()
    
    def _detect_simulator_tool(self) -> str:
        """Detect which read simulation tool is available."""
        # Check for ART (recommended)
        try:
            result = subprocess.run(['art_illumina', '--help'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                self.logger.info("Using ART (Illumina) for read simulation")
                return 'art'
        except FileNotFoundError:
            pass
        
        # Check for BBTools randomreads.sh
        try:
            result = subprocess.run(['randomreads.sh', '--help'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                self.logger.info("Using BBTools randomreads.sh for read simulation")
                return 'bbtools'
        except FileNotFoundError:
            pass
        
        # Fallback to simple simulation
        self.logger.warning("No external read simulator found. Using simple built-in simulator.")
        self.logger.warning("For more realistic reads, install ART or BBTools.")
        return 'simple'
    
    def simulate_reads(self, contigs: List[ContigInfo], 
                      contig_coverages: List[ContigCoverage],
                      read_counts: Dict[str, int],
                      output_dir: Path) -> Tuple[Path, Path]:
        """
        Simulate paired-end reads from all contigs.
        
        Args:
            contigs: List of contig information
            contig_coverages: Coverage information
            read_counts: Number of read pairs per contig
            output_dir: Output directory for reads
            
        Returns:
            Tuple of (R1_file_path, R2_file_path)
        """
        self.logger.info("Starting read simulation...")
        
        # Create temporary files for individual contig reads
        temp_files = []
        
        try:
            # Simulate reads for each contig
            for contig in contigs:
                if contig.id not in read_counts:
                    continue
                
                num_pairs = read_counts[contig.id]
                if num_pairs == 0:
                    continue
                
                self.logger.debug(f"Simulating {num_pairs} read pairs for {contig.id}")
                
                # Create temporary files for this contig
                temp_r1 = tempfile.NamedTemporaryFile(mode='w', suffix='_R1.fastq', delete=False)
                temp_r2 = tempfile.NamedTemporaryFile(mode='w', suffix='_R2.fastq', delete=False)
                temp_files.extend([temp_r1.name, temp_r2.name])
                
                # Simulate reads for this contig
                self._simulate_contig_reads(contig, num_pairs, temp_r1.name, temp_r2.name)
                
                temp_r1.close()
                temp_r2.close()
            
            # Merge all temporary files
            r1_output = output_dir / 'reads_R1.fastq.gz'
            r2_output = output_dir / 'reads_R2.fastq.gz'
            
            self._merge_read_files(temp_files, r1_output, r2_output)
            
        finally:
            # Clean up temporary files
            for temp_file in temp_files:
                try:
                    Path(temp_file).unlink()
                except FileNotFoundError:
                    pass
        
        return r1_output, r2_output
    
    def _simulate_contig_reads(self, contig: ContigInfo, num_pairs: int, 
                              r1_output: str, r2_output: str):
        """Simulate reads for a single contig."""
        if self.simulator_tool == 'art':
            self._simulate_with_art(contig, num_pairs, r1_output, r2_output)
        elif self.simulator_tool == 'bbtools':
            self._simulate_with_bbtools(contig, num_pairs, r1_output, r2_output)
        else:
            self._simulate_simple(contig, num_pairs, r1_output, r2_output)
    
    def _simulate_with_art(self, contig: ContigInfo, num_pairs: int,
                          r1_output: str, r2_output: str):
        """Simulate reads using ART."""
        # Create temporary reference file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_file:
            ref_file.write(f">{contig.id}\n{contig.sequence}\n")
            ref_path = ref_file.name
        
        try:
            # Calculate coverage for ART
            read_length = self.config.read_simulation.read_length
            coverage = (num_pairs * 2 * read_length) / contig.length
            
            # ART command
            output_prefix = r1_output.replace('_R1.fastq', '')
            
            cmd = [
                'art_illumina',
                '-ss', 'HS25',  # HiSeq 2500 profile
                '-i', ref_path,
                '-p',  # Paired-end
                '-l', str(read_length),
                '-f', str(coverage),
                '-m', str(self.config.read_simulation.insert_size_mean),
                '-s', str(self.config.read_simulation.insert_size_std),
                '-o', output_prefix
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"ART failed: {result.stderr}")
            
            # ART creates files with different naming convention
            art_r1 = f"{output_prefix}1.fq"
            art_r2 = f"{output_prefix}2.fq"
            
            # Rename to expected names
            if Path(art_r1).exists():
                Path(art_r1).rename(r1_output)
            if Path(art_r2).exists():
                Path(art_r2).rename(r2_output)
            
        finally:
            Path(ref_path).unlink()
    
    def _simulate_with_bbtools(self, contig: ContigInfo, num_pairs: int,
                              r1_output: str, r2_output: str):
        """Simulate reads using BBTools randomreads.sh."""
        # Create temporary reference file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_file:
            ref_file.write(f">{contig.id}\n{contig.sequence}\n")
            ref_path = ref_file.name
        
        try:
            cmd = [
                'randomreads.sh',
                f'ref={ref_path}',
                f'out1={r1_output}',
                f'out2={r2_output}',
                f'length={self.config.read_simulation.read_length}',
                f'reads={num_pairs}',
                f'insertsize={self.config.read_simulation.insert_size_mean}',
                f'insertdev={self.config.read_simulation.insert_size_std}',
                f'errorrate={self.config.read_simulation.error_rate}',
                'paired=t',
                'gaussian=t'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"BBTools randomreads.sh failed: {result.stderr}")
            
        finally:
            Path(ref_path).unlink()
    
    def _simulate_simple(self, contig: ContigInfo, num_pairs: int,
                        r1_output: str, r2_output: str):
        """Simple built-in read simulation."""
        read_length = self.config.read_simulation.read_length
        insert_mean = self.config.read_simulation.insert_size_mean
        insert_std = self.config.read_simulation.insert_size_std
        error_rate = self.config.read_simulation.error_rate
        
        with open(r1_output, 'w') as f1, open(r2_output, 'w') as f2:
            for i in range(num_pairs):
                # Generate insert size
                insert_size = max(read_length * 2, 
                                int(random.normalvariate(insert_mean, insert_std)))
                
                # Choose random start position
                max_start = len(contig.sequence) - insert_size
                if max_start <= 0:
                    continue
                
                start_pos = random.randint(0, max_start)
                end_pos = start_pos + insert_size
                
                # Extract fragment
                fragment = contig.sequence[start_pos:end_pos]
                
                # Generate R1 (forward)
                r1_seq = fragment[:read_length]
                r1_qual = 'I' * read_length  # High quality
                
                # Generate R2 (reverse complement of end)
                r2_fragment = fragment[-read_length:]
                r2_seq = self._reverse_complement(r2_fragment)
                r2_qual = 'I' * read_length
                
                # Add errors
                r1_seq = self._add_errors(r1_seq, error_rate)
                r2_seq = self._add_errors(r2_seq, error_rate)
                
                # Write FASTQ records
                read_id = f"{contig.id}_{i+1}"
                
                f1.write(f"@{read_id}/1\n{r1_seq}\n+\n{r1_qual}\n")
                f2.write(f"@{read_id}/2\n{r2_seq}\n+\n{r2_qual}\n")
    
    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def _add_errors(self, seq: str, error_rate: float) -> str:
        """Add sequencing errors to sequence."""
        if error_rate <= 0:
            return seq
        
        bases = ['A', 'T', 'G', 'C']
        seq_list = list(seq)
        
        for i, base in enumerate(seq_list):
            if random.random() < error_rate:
                # Introduce error
                error_bases = [b for b in bases if b != base]
                seq_list[i] = random.choice(error_bases)
        
        return ''.join(seq_list)
    
    def _merge_read_files(self, temp_files: List[str], r1_output: Path, r2_output: Path):
        """Merge temporary read files and compress."""
        # Separate R1 and R2 files
        r1_files = [f for f in temp_files if '_R1.fastq' in f]
        r2_files = [f for f in temp_files if '_R2.fastq' in f]
        
        # Merge R1 files
        with gzip.open(r1_output, 'wt') as out_f1:
            for r1_file in r1_files:
                try:
                    with open(r1_file, 'r') as in_f:
                        out_f1.write(in_f.read())
                except FileNotFoundError:
                    continue
        
        # Merge R2 files
        with gzip.open(r2_output, 'wt') as out_f2:
            for r2_file in r2_files:
                try:
                    with open(r2_file, 'r') as in_f:
                        out_f2.write(in_f.read())
                except FileNotFoundError:
                    continue
        
        # Count reads
        r1_count = self._count_reads(r1_output)
        r2_count = self._count_reads(r2_output)
        
        self.logger.info(f"Generated {r1_count:,} R1 reads and {r2_count:,} R2 reads")
        
        if r1_count != r2_count:
            self.logger.warning(f"Mismatch in read counts: R1={r1_count}, R2={r2_count}")
    
    def _count_reads(self, fastq_file: Path) -> int:
        """Count reads in compressed FASTQ file."""
        count = 0
        with gzip.open(fastq_file, 'rt') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
        return count


def simulate_reads_from_contigs(contigs: List[ContigInfo],
                               contig_coverages: List[ContigCoverage],
                               read_counts: Dict[str, int],
                               config, output_dir: Path) -> Tuple[Path, Path]:
    """
    Main function to simulate reads from contigs.
    
    Args:
        contigs: List of contig information
        contig_coverages: Coverage information
        read_counts: Number of read pairs per contig
        config: Configuration object
        output_dir: Output directory for reads
        
    Returns:
        Tuple of (R1_file_path, R2_file_path)
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting read simulation from synthetic contigs...")
    
    simulator = ReadSimulator(config)
    r1_file, r2_file = simulator.simulate_reads(
        contigs, contig_coverages, read_counts, output_dir
    )
    
    logger.info(f"Read simulation completed:")
    logger.info(f"  R1 reads: {r1_file}")
    logger.info(f"  R2 reads: {r2_file}")
    
    return r1_file, r2_file