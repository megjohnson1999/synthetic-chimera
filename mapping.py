"""
Read mapping module for chimera simulator.

Maps simulated reads back to synthetic contigs using BWA-MEM or minimap2,
preserving realistic PE mapping artifacts for chimera detection testing.
"""

import subprocess
import logging
import tempfile
import pysam
from pathlib import Path
from typing import Tuple, Dict, List, Optional

from contig_generators import ContigInfo


class ReadMapper:
    """
    Maps paired-end reads to synthetic contigs.
    
    Uses BWA-MEM or minimap2 to generate realistic BAM files with
    proper PE flags and mapping artifacts at chimeric junctions.
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Detect available mapping tool
        self.mapper_tool = self._detect_mapper_tool()
    
    def _detect_mapper_tool(self) -> str:
        """Detect which mapping tool is available."""
        preferred_tool = self.config.mapping.aligner
        
        # Check for BWA-MEM
        if preferred_tool == "bwa_mem":
            try:
                result = subprocess.run(['bwa'], capture_output=True, text=True)
                if result.returncode == 0 or "Usage:" in result.stderr:
                    self.logger.info("Using BWA-MEM for read mapping")
                    return 'bwa_mem'
            except FileNotFoundError:
                pass
        
        # Check for minimap2
        if preferred_tool == "minimap2":
            try:
                result = subprocess.run(['minimap2', '--version'], capture_output=True, text=True)
                if result.returncode == 0:
                    self.logger.info("Using minimap2 for read mapping")
                    return 'minimap2'
            except FileNotFoundError:
                pass
        
        # Try alternative if preferred not available
        if preferred_tool == "bwa_mem":
            try:
                result = subprocess.run(['minimap2', '--version'], capture_output=True, text=True)
                if result.returncode == 0:
                    self.logger.warning("BWA-MEM not found, using minimap2 instead")
                    return 'minimap2'
            except FileNotFoundError:
                pass
        else:
            try:
                result = subprocess.run(['bwa'], capture_output=True, text=True)
                if result.returncode == 0 or "Usage:" in result.stderr:
                    self.logger.warning("minimap2 not found, using BWA-MEM instead")
                    return 'bwa_mem'
            except FileNotFoundError:
                pass
        
        raise RuntimeError("Neither BWA-MEM nor minimap2 found. Please install one of these mapping tools.")
    
    def map_reads(self, contigs_file: Path, r1_file: Path, r2_file: Path, 
                  output_dir: Path) -> Path:
        """
        Map paired-end reads to synthetic contigs.
        
        Args:
            contigs_file: Path to synthetic assembly FASTA
            r1_file: Path to R1 reads (compressed FASTQ)
            r2_file: Path to R2 reads (compressed FASTQ)
            output_dir: Output directory for BAM file
            
        Returns:
            Path to sorted, indexed BAM file
        """
        self.logger.info("Starting read mapping...")
        
        output_bam = output_dir / 'reads_to_contigs.bam'
        
        if self.mapper_tool == 'bwa_mem':
            bam_file = self._map_with_bwa(contigs_file, r1_file, r2_file, output_bam)
        else:
            bam_file = self._map_with_minimap2(contigs_file, r1_file, r2_file, output_bam)
        
        # Sort and index BAM file
        sorted_bam = self._sort_and_index_bam(bam_file)
        
        # Generate mapping statistics
        self._generate_mapping_stats(sorted_bam, output_dir)
        
        return sorted_bam
    
    def _map_with_bwa(self, contigs_file: Path, r1_file: Path, r2_file: Path,
                      output_bam: Path) -> Path:
        """Map reads using BWA-MEM."""
        self.logger.info("Mapping reads with BWA-MEM...")
        
        # Index reference
        self.logger.info("Indexing reference contigs...")
        index_cmd = ['bwa', 'index', str(contigs_file)]
        result = subprocess.run(index_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BWA index failed: {result.stderr}")
        
        # Map reads
        self.logger.info("Mapping paired-end reads...")
        map_cmd = [
            'bwa', 'mem',
            '-t', str(self.config.mapping.threads),
            '-M',  # Mark shorter split hits as secondary
            str(contigs_file),
            str(r1_file),
            str(r2_file)
        ]
        
        # Convert SAM to BAM
        with open(output_bam, 'wb') as bam_out:
            map_process = subprocess.Popen(map_cmd, stdout=subprocess.PIPE)
            sam_to_bam_cmd = ['samtools', 'view', '-bS', '-']
            subprocess.run(sam_to_bam_cmd, stdin=map_process.stdout, stdout=bam_out)
            map_process.wait()
        
        if map_process.returncode != 0:
            raise RuntimeError("BWA-MEM mapping failed")
        
        return output_bam
    
    def _map_with_minimap2(self, contigs_file: Path, r1_file: Path, r2_file: Path,
                          output_bam: Path) -> Path:
        """Map reads using minimap2."""
        self.logger.info("Mapping reads with minimap2...")
        
        map_cmd = [
            'minimap2',
            '-ax', 'sr',  # Short read mode
            '-t', str(self.config.mapping.threads),
            str(contigs_file),
            str(r1_file),
            str(r2_file)
        ]
        
        # Convert SAM to BAM
        with open(output_bam, 'wb') as bam_out:
            map_process = subprocess.Popen(map_cmd, stdout=subprocess.PIPE)
            sam_to_bam_cmd = ['samtools', 'view', '-bS', '-']
            subprocess.run(sam_to_bam_cmd, stdin=map_process.stdout, stdout=bam_out)
            map_process.wait()
        
        if map_process.returncode != 0:
            raise RuntimeError("minimap2 mapping failed")
        
        return output_bam
    
    def _sort_and_index_bam(self, bam_file: Path) -> Path:
        """Sort and index BAM file."""
        self.logger.info("Sorting and indexing BAM file...")
        
        sorted_bam = bam_file.with_suffix('.sorted.bam')
        
        # Sort BAM
        sort_cmd = [
            'samtools', 'sort',
            '-@', str(self.config.mapping.threads),
            '-o', str(sorted_bam),
            str(bam_file)
        ]
        
        result = subprocess.run(sort_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BAM sorting failed: {result.stderr}")
        
        # Index BAM
        index_cmd = ['samtools', 'index', str(sorted_bam)]
        result = subprocess.run(index_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BAM indexing failed: {result.stderr}")
        
        # Remove unsorted BAM
        bam_file.unlink()
        
        return sorted_bam
    
    def _generate_mapping_stats(self, bam_file: Path, output_dir: Path):
        """Generate mapping statistics."""
        self.logger.info("Generating mapping statistics...")
        
        stats_file = output_dir / 'mapping_stats.txt'
        
        # Use samtools flagstat
        flagstat_cmd = ['samtools', 'flagstat', str(bam_file)]
        result = subprocess.run(flagstat_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            with open(stats_file, 'w') as f:
                f.write("=== Mapping Statistics ===\n\n")
                f.write(result.stdout)
                f.write("\n=== Additional Statistics ===\n")
                
                # Parse basic stats
                lines = result.stdout.strip().split('\n')
                total_reads = 0
                mapped_reads = 0
                
                for line in lines:
                    if 'in total' in line:
                        total_reads = int(line.split()[0])
                    elif 'mapped (' in line and 'paired in sequencing' not in line:
                        mapped_reads = int(line.split()[0])
                
                if total_reads > 0:
                    mapping_rate = mapped_reads / total_reads
                    f.write(f"Mapping rate: {mapping_rate:.3f} ({mapping_rate*100:.1f}%)\n")
                    
                    if mapping_rate < self.config.validation.min_mapping_rate:
                        self.logger.warning(f"Low mapping rate: {mapping_rate:.3f} "
                                          f"(minimum: {self.config.validation.min_mapping_rate:.3f})")
                    else:
                        self.logger.info(f"Mapping rate: {mapping_rate:.3f}")
        
        # Generate insert size statistics
        self._generate_insert_size_stats(bam_file, output_dir)
    
    def _generate_insert_size_stats(self, bam_file: Path, output_dir: Path):
        """Generate insert size statistics using samtools."""
        self.logger.info("Generating insert size statistics...")
        
        insert_stats_file = output_dir / 'insert_size_stats.txt'
        
        # Use samtools stats for insert size distribution
        stats_cmd = ['samtools', 'stats', str(bam_file)]
        result = subprocess.run(stats_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            with open(insert_stats_file, 'w') as f:
                f.write("=== Insert Size Statistics ===\n\n")
                
                # Extract insert size information
                for line in result.stdout.split('\n'):
                    if line.startswith('IS\t'):  # Insert size histogram
                        f.write(line + '\n')
                    elif line.startswith('SN\t') and 'insert size' in line.lower():
                        f.write(line + '\n')
        
        self.logger.info(f"Insert size statistics saved to {insert_stats_file}")


class PEArtifactAnalyzer:
    """
    Analyzes BAM file for paired-end artifacts at chimeric junctions.
    
    Identifies insert size anomalies and orientation problems that
    should be detectable by chimera detection algorithms.
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def analyze_pe_artifacts(self, bam_file: Path, contigs: List[ContigInfo],
                           output_dir: Path) -> Dict:
        """
        Analyze PE mapping artifacts in BAM file.
        
        Args:
            bam_file: Path to sorted, indexed BAM file
            contigs: List of contig information including chimera junctions
            output_dir: Output directory for analysis results
            
        Returns:
            Dictionary with artifact analysis results
        """
        self.logger.info("Analyzing PE mapping artifacts...")
        
        results = {
            'insert_size_anomalies': {},
            'orientation_anomalies': {},
            'chimera_detection_signals': {}
        }
        
        # Create contig lookup
        contig_lookup = {contig.id: contig for contig in contigs}
        
        # Analyze BAM file
        with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
            # Calculate global insert size statistics
            insert_sizes = self._calculate_insert_size_stats(bam)
            results['global_insert_stats'] = insert_sizes
            
            # Analyze chimeric contigs specifically
            for contig in contigs:
                if contig.is_chimeric:
                    contig_results = self._analyze_chimeric_contig(
                        bam, contig, insert_sizes
                    )
                    results['chimera_detection_signals'][contig.id] = contig_results
        
        # Save analysis results
        self._save_artifact_analysis(results, output_dir)
        
        return results
    
    def _calculate_insert_size_stats(self, bam: pysam.AlignmentFile) -> Dict:
        """Calculate global insert size statistics."""
        insert_sizes = []
        
        for read in bam.fetch():
            if (read.is_proper_pair and not read.is_secondary and 
                not read.is_supplementary and read.template_length > 0):
                insert_sizes.append(abs(read.template_length))
        
        if not insert_sizes:
            return {'mean': 0, 'std': 0, 'median': 0}
        
        import numpy as np
        return {
            'mean': np.mean(insert_sizes),
            'std': np.std(insert_sizes),
            'median': np.median(insert_sizes),
            'count': len(insert_sizes)
        }
    
    def _analyze_chimeric_contig(self, bam: pysam.AlignmentFile, 
                               contig: ContigInfo, global_stats: Dict) -> Dict:
        """Analyze PE artifacts for a specific chimeric contig."""
        results = {
            'junction_coords': contig.junction_coords,
            'expected_artifacts': contig.expected_pe_artifacts,
            'observed_artifacts': {
                'large_inserts': [],
                'small_inserts': [],
                'improper_orientations': [],
                'discordant_pairs': []
            }
        }
        
        # Analyze reads mapping to this contig
        for read in bam.fetch(contig.id):
            if read.is_secondary or read.is_supplementary:
                continue
            
            # Check for insert size anomalies
            if read.is_paired and read.template_length != 0:
                insert_size = abs(read.template_length)
                
                # Check if insert size is anomalous
                if global_stats['count'] > 0:
                    z_score = (insert_size - global_stats['mean']) / max(global_stats['std'], 1)
                    
                    if z_score > self.config.validation.insert_size_outlier_threshold:
                        results['observed_artifacts']['large_inserts'].append({
                            'read_name': read.query_name,
                            'position': read.reference_start,
                            'insert_size': insert_size,
                            'z_score': z_score
                        })
                    elif z_score < -self.config.validation.insert_size_outlier_threshold:
                        results['observed_artifacts']['small_inserts'].append({
                            'read_name': read.query_name,
                            'position': read.reference_start,
                            'insert_size': insert_size,
                            'z_score': z_score
                        })
            
            # Check for orientation anomalies
            if read.is_paired and not read.is_proper_pair:
                results['observed_artifacts']['improper_orientations'].append({
                    'read_name': read.query_name,
                    'position': read.reference_start,
                    'flags': read.flag
                })
        
        return results
    
    def _save_artifact_analysis(self, results: Dict, output_dir: Path):
        """Save artifact analysis results."""
        analysis_file = output_dir / 'pe_artifact_analysis.tsv'
        
        with open(analysis_file, 'w') as f:
            f.write("=== PE Artifact Analysis Results ===\n\n")
            
            # Global statistics
            f.write("Global Insert Size Statistics:\n")
            stats = results['global_insert_stats']
            f.write(f"Mean: {stats['mean']:.2f}\n")
            f.write(f"Std: {stats['std']:.2f}\n")
            f.write(f"Median: {stats['median']:.2f}\n")
            f.write(f"Count: {stats['count']}\n\n")
            
            # Per-chimera analysis
            f.write("Chimeric Contig Analysis:\n")
            f.write("contig_id\tjunction_coords\tlarge_inserts\tsmall_inserts\timproper_orientations\n")
            
            for contig_id, analysis in results['chimera_detection_signals'].items():
                f.write(f"{contig_id}\t")
                f.write(f"{';'.join(map(str, analysis['junction_coords']))}\t")
                f.write(f"{len(analysis['observed_artifacts']['large_inserts'])}\t")
                f.write(f"{len(analysis['observed_artifacts']['small_inserts'])}\t")
                f.write(f"{len(analysis['observed_artifacts']['improper_orientations'])}\n")
        
        self.logger.info(f"PE artifact analysis saved to {analysis_file}")


def map_reads_to_contigs(contigs_file: Path, r1_file: Path, r2_file: Path,
                        contigs: List[ContigInfo], config, output_dir: Path) -> Tuple[Path, Dict]:
    """
    Main function to map reads to contigs and analyze PE artifacts.
    
    Args:
        contigs_file: Path to synthetic assembly FASTA
        r1_file: Path to R1 reads
        r2_file: Path to R2 reads
        contigs: List of contig information
        config: Configuration object
        output_dir: Output directory
        
    Returns:
        Tuple of (bam_file_path, pe_artifact_analysis)
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting read mapping and PE artifact analysis...")
    
    # Map reads
    mapper = ReadMapper(config)
    bam_file = mapper.map_reads(contigs_file, r1_file, r2_file, output_dir)
    
    # Analyze PE artifacts
    analyzer = PEArtifactAnalyzer(config)
    pe_analysis = analyzer.analyze_pe_artifacts(bam_file, contigs, output_dir)
    
    logger.info(f"Mapping and analysis completed:")
    logger.info(f"  BAM file: {bam_file}")
    logger.info(f"  PE artifacts detected in {len(pe_analysis['chimera_detection_signals'])} chimeric contigs")
    
    return bam_file, pe_analysis