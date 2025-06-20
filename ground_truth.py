"""
Ground truth generation for chimera simulator.

Creates comprehensive ground truth data focused on junction coordinates
and expected paired-end mapping artifacts for algorithm testing.
"""

import logging
from pathlib import Path
from typing import List, Dict, Any
import json

from contig_generators import ContigInfo
from coverage_model import ContigCoverage


class GroundTruthGenerator:
    """
    Generates ground truth data for chimera detection algorithm testing.
    
    Focuses on junction coordinates, PE artifact predictions, and
    detection difficulty classification.
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def generate_ground_truth(self, contigs: List[ContigInfo],
                            contig_coverages: List[ContigCoverage],
                            pe_analysis: Dict, output_dir: Path) -> Dict[str, Path]:
        """
        Generate comprehensive ground truth files.
        
        Args:
            contigs: List of contig information
            contig_coverages: Coverage information
            pe_analysis: PE artifact analysis results
            output_dir: Output directory for ground truth files
            
        Returns:
            Dictionary mapping file types to file paths
        """
        self.logger.info("Generating ground truth data...")
        
        output_files = {}
        
        # Generate chimera annotations
        chimera_file = self._generate_chimera_annotations(contigs, output_dir)
        output_files['chimera_annotations'] = chimera_file
        
        # Generate junction details
        junction_file = self._generate_junction_details(contigs, contig_coverages, output_dir)
        output_files['junction_details'] = junction_file
        
        # Generate PE artifact predictions
        pe_artifacts_file = self._generate_pe_artifact_predictions(contigs, pe_analysis, output_dir)
        output_files['pe_artifacts'] = pe_artifacts_file
        
        # Generate detection difficulty assessment
        difficulty_file = self._generate_difficulty_assessment(contigs, contig_coverages, output_dir)
        output_files['difficulty_assessment'] = difficulty_file
        
        # Generate comprehensive JSON summary
        json_file = self._generate_json_summary(contigs, contig_coverages, pe_analysis, output_dir)
        output_files['json_summary'] = json_file
        
        self.logger.info(f"Ground truth generation completed. Files created:")
        for file_type, file_path in output_files.items():
            self.logger.info(f"  {file_type}: {file_path}")
        
        return output_files
    
    def _generate_chimera_annotations(self, contigs: List[ContigInfo], 
                                    output_dir: Path) -> Path:
        """Generate basic chimera status annotations."""
        output_file = output_dir / 'chimera_annotations.tsv'
        
        with open(output_file, 'w') as f:
            f.write("contig_id\tis_chimeric\tchimera_type\tlength\tnum_junctions\n")
            
            for contig in contigs:
                f.write(f"{contig.id}\t")
                f.write(f"{contig.is_chimeric}\t")
                f.write(f"{contig.chimera_type or 'normal'}\t")
                f.write(f"{contig.length}\t")
                f.write(f"{len(contig.junction_coords)}\n")
        
        return output_file
    
    def _generate_junction_details(self, contigs: List[ContigInfo],
                                 contig_coverages: List[ContigCoverage],
                                 output_dir: Path) -> Path:
        """Generate detailed junction information."""
        output_file = output_dir / 'junction_details.tsv'
        
        # Create coverage lookup
        coverage_lookup = {cc.contig_id: cc for cc in contig_coverages}
        
        with open(output_file, 'w') as f:
            f.write("contig_id\tjunction_coord\tchimera_type\ttotal_coverage\t")
            f.write("avg_coverage\tleft_segment_info\tright_segment_info\n")
            
            for contig in contigs:
                if not contig.is_chimeric:
                    continue
                
                coverage_info = coverage_lookup.get(contig.id)
                total_cov = coverage_info.total_coverage if coverage_info else 0
                avg_cov = coverage_info.average_coverage if coverage_info else 0
                
                for junction_coord in contig.junction_coords:
                    f.write(f"{contig.id}\t{junction_coord}\t{contig.chimera_type}\t")
                    f.write(f"{total_cov:.2f}\t{avg_cov:.2f}\t")
                    
                    # Add segment information
                    if len(contig.source_segments) >= 2:
                        left_seg = contig.source_segments[0]
                        right_seg = contig.source_segments[1]
                        
                        left_info = f"{left_seg['genome_id']}:{left_seg['start']}-{left_seg['end']}({left_seg['strand']})"
                        right_info = f"{right_seg['genome_id']}:{right_seg['start']}-{right_seg['end']}({right_seg['strand']})"
                        
                        f.write(f"{left_info}\t{right_info}\n")
                    else:
                        f.write("unknown\tunknown\n")
        
        return output_file
    
    def _generate_pe_artifact_predictions(self, contigs: List[ContigInfo],
                                        pe_analysis: Dict, output_dir: Path) -> Path:
        """Generate PE artifact predictions and observed results."""
        output_file = output_dir / 'pe_artifact_predictions.tsv'
        
        with open(output_file, 'w') as f:
            f.write("contig_id\tjunction_coord\texpected_artifact_type\t")
            f.write("predicted_signal\tobserved_large_inserts\tobserved_small_inserts\t")
            f.write("observed_improper_orientations\tdetection_confidence\n")
            
            for contig in contigs:
                if not contig.is_chimeric:
                    continue
                
                # Get observed artifacts from PE analysis
                observed = pe_analysis.get('chimera_detection_signals', {}).get(contig.id, {})
                observed_artifacts = observed.get('observed_artifacts', {})
                
                for junction_coord in contig.junction_coords:
                    # Determine expected artifact type
                    expected_type = self._get_expected_artifact_type(contig, junction_coord)
                    predicted_signal = self._get_predicted_signal_strength(contig, junction_coord)
                    
                    # Count observed artifacts
                    large_inserts = len(observed_artifacts.get('large_inserts', []))
                    small_inserts = len(observed_artifacts.get('small_inserts', []))
                    improper_orient = len(observed_artifacts.get('improper_orientations', []))
                    
                    # Calculate detection confidence
                    confidence = self._calculate_detection_confidence(
                        expected_type, large_inserts, small_inserts, improper_orient
                    )
                    
                    f.write(f"{contig.id}\t{junction_coord}\t{expected_type}\t")
                    f.write(f"{predicted_signal}\t{large_inserts}\t{small_inserts}\t")
                    f.write(f"{improper_orient}\t{confidence}\n")
        
        return output_file
    
    def _generate_difficulty_assessment(self, contigs: List[ContigInfo],
                                      contig_coverages: List[ContigCoverage],
                                      output_dir: Path) -> Path:
        """Generate detection difficulty assessment."""
        output_file = output_dir / 'difficulty_assessment.tsv'
        
        # Create coverage lookup
        coverage_lookup = {cc.contig_id: cc for cc in contig_coverages}
        
        with open(output_file, 'w') as f:
            f.write("contig_id\tdifficulty_level\tdifficulty_factors\ttotal_coverage\t")
            f.write("coverage_variability\trecommended_for_testing\n")
            
            for contig in contigs:
                if not contig.is_chimeric:
                    continue
                
                coverage_info = coverage_lookup.get(contig.id)
                difficulty_level, factors = self._assess_detection_difficulty(contig, coverage_info)
                
                total_cov = coverage_info.total_coverage if coverage_info else 0
                cov_var = coverage_info.coverage_std if coverage_info else 0
                
                # Recommend for testing based on difficulty and coverage
                recommended = self._recommend_for_testing(difficulty_level, total_cov)
                
                f.write(f"{contig.id}\t{difficulty_level}\t{';'.join(factors)}\t")
                f.write(f"{total_cov:.2f}\t{cov_var:.2f}\t{recommended}\n")
        
        return output_file
    
    def _generate_json_summary(self, contigs: List[ContigInfo],
                             contig_coverages: List[ContigCoverage],
                             pe_analysis: Dict, output_dir: Path) -> Path:
        """Generate comprehensive JSON summary."""
        output_file = output_dir / 'ground_truth_summary.json'
        
        # Create coverage lookup
        coverage_lookup = {cc.contig_id: cc for cc in contig_coverages}
        
        summary = {
            'metadata': {
                'total_contigs': len(contigs),
                'chimeric_contigs': sum(1 for c in contigs if c.is_chimeric),
                'normal_contigs': sum(1 for c in contigs if not c.is_chimeric),
                'chimera_types': {}
            },
            'contigs': {},
            'global_stats': pe_analysis.get('global_insert_stats', {}),
            'detection_summary': {}
        }
        
        # Count chimera types
        for contig in contigs:
            if contig.is_chimeric and contig.chimera_type:
                chimera_type = contig.chimera_type
                summary['metadata']['chimera_types'][chimera_type] = \
                    summary['metadata']['chimera_types'].get(chimera_type, 0) + 1
        
        # Add detailed contig information
        for contig in contigs:
            coverage_info = coverage_lookup.get(contig.id)
            
            contig_data = {
                'is_chimeric': contig.is_chimeric,
                'chimera_type': contig.chimera_type,
                'length': contig.length,
                'junction_coords': contig.junction_coords,
                'source_segments': contig.source_segments,
                'expected_pe_artifacts': contig.expected_pe_artifacts,
                'coverage': {
                    'total': coverage_info.total_coverage if coverage_info else 0,
                    'average': coverage_info.average_coverage if coverage_info else 0,
                    'std': coverage_info.coverage_std if coverage_info else 0,
                    'per_sample': coverage_info.sample_coverages if coverage_info else {}
                }
            }
            
            # Add observed artifacts if available
            if contig.id in pe_analysis.get('chimera_detection_signals', {}):
                contig_data['observed_artifacts'] = pe_analysis['chimera_detection_signals'][contig.id]
            
            summary['contigs'][contig.id] = contig_data
        
        # Add detection summary
        chimeric_contigs = [c for c in contigs if c.is_chimeric]
        if chimeric_contigs:
            summary['detection_summary'] = {
                'easy_detection': sum(1 for c in chimeric_contigs 
                                    if self._assess_detection_difficulty(c, coverage_lookup.get(c.id))[0] == 'easy'),
                'medium_detection': sum(1 for c in chimeric_contigs 
                                      if self._assess_detection_difficulty(c, coverage_lookup.get(c.id))[0] == 'medium'),
                'hard_detection': sum(1 for c in chimeric_contigs 
                                    if self._assess_detection_difficulty(c, coverage_lookup.get(c.id))[0] == 'hard')
            }
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        return output_file
    
    def _get_expected_artifact_type(self, contig: ContigInfo, junction_coord: int) -> str:
        """Determine expected PE artifact type for a junction."""
        if contig.chimera_type == 'gap_mediated':
            return 'large_insert_size'
        elif contig.chimera_type == 'orientation_flip':
            return 'improper_orientation'
        elif contig.chimera_type == 'distant_join':
            return 'extreme_insert_size'
        else:
            return 'unknown'
    
    def _get_predicted_signal_strength(self, contig: ContigInfo, junction_coord: int) -> str:
        """Predict signal strength for detection."""
        if contig.chimera_type == 'gap_mediated':
            # Gap size affects signal strength
            return 'strong'  # Gaps create clear signals
        elif contig.chimera_type == 'orientation_flip':
            return 'strong'  # Orientation flips are usually clear
        elif contig.chimera_type == 'distant_join':
            return 'medium'  # Depends on distance and coverage
        else:
            return 'unknown'
    
    def _calculate_detection_confidence(self, expected_type: str, large_inserts: int,
                                      small_inserts: int, improper_orient: int) -> str:
        """Calculate detection confidence based on observed artifacts."""
        if expected_type == 'large_insert_size' and large_inserts > 0:
            return 'high'
        elif expected_type == 'improper_orientation' and improper_orient > 0:
            return 'high'
        elif expected_type == 'extreme_insert_size' and (large_inserts > 0 or small_inserts > 0):
            return 'medium'
        else:
            return 'low'
    
    def _assess_detection_difficulty(self, contig: ContigInfo, 
                                   coverage_info: ContigCoverage) -> tuple:
        """Assess detection difficulty for a chimeric contig."""
        factors = []
        
        # Coverage-based factors
        if coverage_info:
            if coverage_info.total_coverage < 20:
                factors.append('low_coverage')
            elif coverage_info.total_coverage > 200:
                factors.append('high_coverage')
            
            if coverage_info.coverage_std > coverage_info.average_coverage:
                factors.append('high_coverage_variability')
        
        # Length-based factors
        if contig.length < 1000:
            factors.append('short_contig')
        elif contig.length > 5000:
            factors.append('long_contig')
        
        # Chimera type factors
        if contig.chimera_type == 'gap_mediated':
            factors.append('clear_gap_signal')
        elif contig.chimera_type == 'orientation_flip':
            factors.append('clear_orientation_signal')
        elif contig.chimera_type == 'distant_join':
            factors.append('distance_dependent_signal')
        
        # Determine overall difficulty
        if 'low_coverage' in factors:
            difficulty = 'hard'
        elif 'clear_gap_signal' in factors or 'clear_orientation_signal' in factors:
            difficulty = 'easy'
        else:
            difficulty = 'medium'
        
        return difficulty, factors
    
    def _recommend_for_testing(self, difficulty_level: str, total_coverage: float) -> str:
        """Recommend whether contig is good for testing."""
        if total_coverage < 10:
            return 'no_insufficient_coverage'
        elif difficulty_level == 'easy' and total_coverage > 30:
            return 'yes_positive_control'
        elif difficulty_level == 'hard' and total_coverage > 50:
            return 'yes_challenging_case'
        elif difficulty_level == 'medium':
            return 'yes_standard_case'
        else:
            return 'maybe_edge_case'


def generate_ground_truth_files(contigs: List[ContigInfo],
                               contig_coverages: List[ContigCoverage],
                               pe_analysis: Dict, config,
                               output_dir: Path) -> Dict[str, Path]:
    """
    Main function to generate all ground truth files.
    
    Args:
        contigs: List of contig information
        contig_coverages: Coverage information
        pe_analysis: PE artifact analysis results
        config: Configuration object
        output_dir: Output directory
        
    Returns:
        Dictionary mapping file types to file paths
    """
    logger = logging.getLogger(__name__)
    logger.info("Generating comprehensive ground truth data...")
    
    generator = GroundTruthGenerator(config)
    ground_truth_files = generator.generate_ground_truth(
        contigs, contig_coverages, pe_analysis, output_dir
    )
    
    logger.info("Ground truth generation completed successfully")
    return ground_truth_files