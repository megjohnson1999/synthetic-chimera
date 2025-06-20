"""
Contig generation modules for chimera simulator.

Generates normal and chimeric contigs with PE mapping artifacts.
"""

import random
import logging
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from pathlib import Path

from utils import parse_fasta, reverse_complement


@dataclass
class ContigInfo:
    """Information about a generated contig."""
    id: str
    sequence: str
    length: int
    is_chimeric: bool
    chimera_type: Optional[str] = None
    junction_coords: List[int] = None
    source_segments: List[Dict] = None
    expected_pe_artifacts: Dict = None
    
    def __post_init__(self):
        if self.junction_coords is None:
            self.junction_coords = []
        if self.source_segments is None:
            self.source_segments = []
        if self.expected_pe_artifacts is None:
            self.expected_pe_artifacts = {}


class GenomePool:
    """Pool of input genomes for contig generation."""
    
    def __init__(self, fasta_path: Path):
        self.genomes = {}
        self.genome_list = []
        self.logger = logging.getLogger(__name__)
        
        self._load_genomes(fasta_path)
    
    def _load_genomes(self, fasta_path: Path):
        """Load genomes from FASTA file."""
        for header, sequence in parse_fasta(fasta_path):
            genome_id = header.split()[0]  # Use first part of header as ID
            self.genomes[genome_id] = sequence.upper()
            self.genome_list.append(genome_id)
        
        self.logger.info(f"Loaded {len(self.genomes)} genomes from {fasta_path}")
    
    def get_random_fragment(self, min_length: int, max_length: int) -> Tuple[str, Dict]:
        """
        Get random fragment from random genome.
        
        Returns:
            Tuple of (sequence, source_info)
        """
        genome_id = random.choice(self.genome_list)
        genome_seq = self.genomes[genome_id]
        
        # Choose fragment length
        fragment_length = random.randint(min_length, max_length)
        fragment_length = min(fragment_length, len(genome_seq))
        
        # Choose random start position
        max_start = len(genome_seq) - fragment_length
        if max_start <= 0:
            start = 0
            fragment_length = len(genome_seq)
        else:
            start = random.randint(0, max_start)
        
        end = start + fragment_length
        fragment_seq = genome_seq[start:end]
        
        source_info = {
            'genome_id': genome_id,
            'start': start,
            'end': end,
            'length': fragment_length,
            'strand': '+'
        }
        
        return fragment_seq, source_info
    
    def get_distant_fragments(self, min_length: int, max_length: int, 
                            min_distance: int = 10000) -> Tuple[str, str, Dict, Dict]:
        """
        Get two fragments that are distant from each other or from different genomes.
        
        Returns:
            Tuple of (seq1, seq2, source1_info, source2_info)
        """
        # Decide whether to use same genome or different genomes
        use_different_genomes = random.choice([True, False])
        
        if use_different_genomes or len(self.genome_list) == 1:
            # Use different genomes
            seq1, source1 = self.get_random_fragment(min_length, max_length)
            seq2, source2 = self.get_random_fragment(min_length, max_length)
        else:
            # Use same genome but distant locations
            genome_id = random.choice(self.genome_list)
            genome_seq = self.genomes[genome_id]
            
            # Get first fragment
            frag1_length = random.randint(min_length, max_length)
            frag1_length = min(frag1_length, len(genome_seq) // 2)
            
            max_start1 = len(genome_seq) - frag1_length - min_distance
            if max_start1 <= 0:
                # Genome too small for distant fragments, use different genomes
                seq1, source1 = self.get_random_fragment(min_length, max_length)
                seq2, source2 = self.get_random_fragment(min_length, max_length)
                return seq1, seq2, source1, source2
            
            start1 = random.randint(0, max_start1)
            end1 = start1 + frag1_length
            seq1 = genome_seq[start1:end1]
            
            # Get second fragment at least min_distance away
            frag2_length = random.randint(min_length, max_length)
            frag2_length = min(frag2_length, len(genome_seq) - end1 - min_distance)
            
            if frag2_length < min_length:
                # Not enough space, use different genomes
                seq1, source1 = self.get_random_fragment(min_length, max_length)
                seq2, source2 = self.get_random_fragment(min_length, max_length)
                return seq1, seq2, source1, source2
            
            start2 = random.randint(end1 + min_distance, len(genome_seq) - frag2_length)
            end2 = start2 + frag2_length
            seq2 = genome_seq[start2:end2]
            
            source1 = {
                'genome_id': genome_id, 'start': start1, 'end': end1,
                'length': frag1_length, 'strand': '+'
            }
            source2 = {
                'genome_id': genome_id, 'start': start2, 'end': end2,
                'length': frag2_length, 'strand': '+'
            }
        
        return seq1, seq2, source1, source2


class ContigGenerator:
    """Base class for contig generators."""
    
    def __init__(self, genome_pool: GenomePool, config):
        self.genome_pool = genome_pool
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def generate(self, contig_id: str) -> ContigInfo:
        """Generate a contig. To be implemented by subclasses."""
        raise NotImplementedError


class NormalContigGenerator(ContigGenerator):
    """Generator for normal (non-chimeric) contigs."""
    
    def generate(self, contig_id: str) -> ContigInfo:
        """Generate a normal contig from a single genome fragment."""
        min_len, max_len = self.config.assembly.normal_contig_length_range
        
        sequence, source_info = self.genome_pool.get_random_fragment(min_len, max_len)
        
        return ContigInfo(
            id=contig_id,
            sequence=sequence,
            length=len(sequence),
            is_chimeric=False,
            source_segments=[source_info]
        )


class GapMediatedChimeraGenerator(ContigGenerator):
    """
    Generator for gap-mediated chimeras.
    
    Creates chimeras with Ns between fragments, causing clear insert size jumps.
    """
    
    def generate(self, contig_id: str) -> ContigInfo:
        """Generate gap-mediated chimera."""
        min_len, max_len = self.config.assembly.chimeric_contig_length_range
        
        # Reserve space for gap
        gap_size = random.randint(50, 500)  # Gap of Ns
        available_length = random.randint(min_len, max_len) - gap_size
        
        # Split available length between two fragments
        frag1_length = random.randint(available_length // 4, 3 * available_length // 4)
        frag2_length = available_length - frag1_length
        
        # Get fragments
        seq1, source1 = self.genome_pool.get_random_fragment(frag1_length, frag1_length)
        seq2, source2 = self.genome_pool.get_random_fragment(frag2_length, frag2_length)
        
        # Create gap
        gap = 'N' * gap_size
        
        # Assemble chimeric sequence
        chimeric_seq = seq1 + gap + seq2
        junction_coord = len(seq1) + gap_size // 2  # Middle of gap
        
        # Expected PE artifacts: abnormally large insert sizes spanning gap
        expected_artifacts = {
            'insert_size_jump': {
                'location': junction_coord,
                'expected_increase': gap_size,
                'artifact_type': 'large_insert'
            }
        }
        
        return ContigInfo(
            id=contig_id,
            sequence=chimeric_seq,
            length=len(chimeric_seq),
            is_chimeric=True,
            chimera_type='gap_mediated',
            junction_coords=[junction_coord],
            source_segments=[source1, source2],
            expected_pe_artifacts=expected_artifacts
        )


class OrientationFlipChimeraGenerator(ContigGenerator):
    """
    Generator for orientation-flip chimeras.
    
    Creates chimeras where second fragment is reverse complemented,
    causing improper pair orientations.
    """
    
    def generate(self, contig_id: str) -> ContigInfo:
        """Generate orientation-flip chimera."""
        min_len, max_len = self.config.assembly.chimeric_contig_length_range
        total_length = random.randint(min_len, max_len)
        
        # Split length between two fragments
        frag1_length = random.randint(total_length // 4, 3 * total_length // 4)
        frag2_length = total_length - frag1_length
        
        # Get fragments
        seq1, source1 = self.genome_pool.get_random_fragment(frag1_length, frag1_length)
        seq2, source2 = self.genome_pool.get_random_fragment(frag2_length, frag2_length)
        
        # Reverse complement second fragment
        seq2_rc = reverse_complement(seq2)
        source2['strand'] = '-'  # Mark as reverse strand
        
        # Assemble chimeric sequence
        chimeric_seq = seq1 + seq2_rc
        junction_coord = len(seq1)
        
        # Expected PE artifacts: improper pair orientations at junction
        expected_artifacts = {
            'orientation_flip': {
                'location': junction_coord,
                'expected_orientation': 'improper',
                'artifact_type': 'orientation_mismatch'
            }
        }
        
        return ContigInfo(
            id=contig_id,
            sequence=chimeric_seq,
            length=len(chimeric_seq),
            is_chimeric=True,
            chimera_type='orientation_flip',
            junction_coords=[junction_coord],
            source_segments=[source1, source2],
            expected_pe_artifacts=expected_artifacts
        )


class DistantJoinChimeraGenerator(ContigGenerator):
    """
    Generator for distant-join chimeras.
    
    Creates chimeras by joining fragments from distant genomic locations,
    causing extreme insert size anomalies.
    """
    
    def generate(self, contig_id: str) -> ContigInfo:
        """Generate distant-join chimera."""
        min_len, max_len = self.config.assembly.chimeric_contig_length_range
        total_length = random.randint(min_len, max_len)
        
        # Split length between two fragments
        frag1_length = random.randint(total_length // 4, 3 * total_length // 4)
        frag2_length = total_length - frag1_length
        
        # Get distant fragments
        seq1, seq2, source1, source2 = self.genome_pool.get_distant_fragments(
            frag1_length, frag2_length
        )
        
        # Assemble chimeric sequence
        chimeric_seq = seq1 + seq2
        junction_coord = len(seq1)
        
        # Calculate expected distance between fragments
        if source1['genome_id'] == source2['genome_id']:
            genomic_distance = abs(source2['start'] - source1['end'])
        else:
            genomic_distance = float('inf')  # Different genomes
        
        # Expected PE artifacts: extreme insert sizes for reads spanning junction
        expected_artifacts = {
            'distant_join': {
                'location': junction_coord,
                'genomic_distance': genomic_distance,
                'artifact_type': 'extreme_insert'
            }
        }
        
        return ContigInfo(
            id=contig_id,
            sequence=chimeric_seq,
            length=len(chimeric_seq),
            is_chimeric=True,
            chimera_type='distant_join',
            junction_coords=[junction_coord],
            source_segments=[source1, source2],
            expected_pe_artifacts=expected_artifacts
        )


def generate_contigs(genome_pool: GenomePool, config) -> List[ContigInfo]:
    """
    Generate complete set of normal and chimeric contigs.
    
    Args:
        genome_pool: Pool of input genomes
        config: Configuration object
        
    Returns:
        List of ContigInfo objects
    """
    logger = logging.getLogger(__name__)
    
    total_contigs = config.assembly.total_contigs
    chimera_percentage = config.assembly.chimera_percentage
    num_chimeric = int(total_contigs * chimera_percentage)
    num_normal = total_contigs - num_chimeric
    
    logger.info(f"Generating {total_contigs} contigs ({num_normal} normal, {num_chimeric} chimeric)")
    
    # Initialize generators
    normal_gen = NormalContigGenerator(genome_pool, config)
    gap_gen = GapMediatedChimeraGenerator(genome_pool, config)
    orient_gen = OrientationFlipChimeraGenerator(genome_pool, config)
    distant_gen = DistantJoinChimeraGenerator(genome_pool, config)
    
    # Determine chimera type distribution
    chimera_types = config.chimera_types
    num_gap = int(num_chimeric * chimera_types.gap_mediated)
    num_orient = int(num_chimeric * chimera_types.orientation_flip)
    num_distant = num_chimeric - num_gap - num_orient  # Remainder
    
    logger.info(f"Chimera distribution - Gap: {num_gap}, Orientation: {num_orient}, Distant: {num_distant}")
    
    contigs = []
    contig_counter = 1
    
    # Generate normal contigs
    for i in range(num_normal):
        contig_id = f"contig_{contig_counter:06d}"
        contig = normal_gen.generate(contig_id)
        contigs.append(contig)
        contig_counter += 1
    
    # Generate chimeric contigs
    generators_and_counts = [
        (gap_gen, num_gap),
        (orient_gen, num_orient),
        (distant_gen, num_distant)
    ]
    
    for generator, count in generators_and_counts:
        for i in range(count):
            contig_id = f"contig_{contig_counter:06d}"
            contig = generator.generate(contig_id)
            contigs.append(contig)
            contig_counter += 1
    
    # Shuffle to mix normal and chimeric contigs
    random.shuffle(contigs)
    
    logger.info(f"Generated {len(contigs)} contigs successfully")
    return contigs