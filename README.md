# Metagenomic Chimera Simulator

A Python tool that generates realistic synthetic metagenomic assemblies containing chimeric contigs with corresponding reads and mappings for testing chimera detection algorithms.

## Overview

This simulator creates synthetic datasets specifically designed for testing paired-end mapping-based chimera detection algorithms. It generates:

- **Synthetic contig assemblies** with realistic normal and chimeric contigs
- **Paired-end reads** simulated from synthetic contigs with authentic insert size distributions
- **BAM mappings** with realistic PE artifacts at chimeric junctions
- **Comprehensive ground truth** data for algorithm benchmarking

## Key Features

### PE Artifact-Focused Design
- **Gap-mediated chimeras**: Create insert size jumps when reads span N-gaps
- **Orientation-flip chimeras**: Generate improper pair orientations (FR→RF/FF/RR)
- **Distant-join chimeras**: Produce extreme insert sizes from distant genomic fragments

### Coassembly Realism
- **Multi-sample coverage simulation**: 5 samples with independent coverage profiles
- **Log-normal coverage distribution**: Realistic 0.1x-500x coverage range
- **Coverage independence**: Coverage patterns don't reveal chimera locations

### Algorithm Testing Ready
- **Multiple difficulty levels**: Easy/medium/hard detection classification
- **Precise ground truth**: Exact junction coordinates and expected PE artifacts
- **Benchmarking support**: Comprehensive metadata for systematic evaluation

## Installation

### Requirements

**Python dependencies:**
```bash
pip install -r requirements.txt
```

**External tools (choose one mapping tool):**
- **BWA-MEM** (recommended): `conda install bwa` or compile from source
- **minimap2** (alternative): `conda install minimap2`

**Required for all installations:**
- **samtools**: `conda install samtools`

**Optional read simulators (for higher quality):**
- **ART** (recommended): Install from [ART website](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
- **BBTools** (alternative): `conda install bbtools`

*Note: If external read simulators are not available, the tool includes a built-in simulator.*

## Usage

### Basic Usage
```bash
python chimera_simulator.py --genomes viral_genomes.fasta
```

### Custom Configuration
```bash
python chimera_simulator.py --genomes viral_genomes.fasta --config custom_params.yaml --output experiment_1/
```

### Pipeline Integration
```bash
for dataset in datasets/*.fasta; do
    python chimera_simulator.py --genomes "$dataset" --output "results/$(basename "$dataset" .fasta)"
done
```

## Command Line Options

| Option | Required | Description |
|--------|----------|-------------|
| `--genomes` | Yes | Path to input viral genomes FASTA file |
| `--config` | No | YAML configuration file (uses defaults if not provided) |
| `--output` | No | Output directory (auto-generated timestamp if not provided) |
| `--seed` | No | Random seed for reproducibility (default: 42) |
| `--verbose` | No | Enable verbose logging |

## Configuration

The tool accepts a YAML configuration file. See `example_config.yaml` for all available options:

```yaml
assembly:
  total_contigs: 1000
  chimera_percentage: 0.15
  normal_contig_length_range: [500, 5000]
  chimeric_contig_length_range: [1000, 8000]

chimera_types:
  gap_mediated: 0.4      # Insert size jumps
  orientation_flip: 0.3   # Improper pair orientations
  distant_join: 0.3       # Extreme insert sizes

coverage:
  samples: 5
  base_coverage_range: [0.1, 500.0]
  distribution: "log_normal"

read_simulation:
  coverage_target: 100.0
  read_length: 150
  insert_size_mean: 300
  insert_size_std: 50

mapping:
  aligner: "bwa_mem"  # or "minimap2"
  threads: 4
```

## Output Structure

```
output/
├── contigs/
│   └── synthetic_assembly.fasta
├── reads/
│   ├── reads_R1.fastq.gz
│   └── reads_R2.fastq.gz
├── mappings/
│   ├── reads_to_contigs.sorted.bam
│   ├── reads_to_contigs.sorted.bam.bai
│   ├── mapping_stats.txt
│   └── insert_size_stats.txt
├── ground_truth/
│   ├── chimera_annotations.tsv
│   ├── junction_details.tsv
│   ├── pe_artifact_predictions.tsv
│   ├── difficulty_assessment.tsv
│   └── ground_truth_summary.json
├── reports/
│   ├── coverage_info.tsv
│   ├── pe_artifact_analysis.tsv
│   └── simulation.log
└── config/
    └── used_config.yaml
```

## Ground Truth Files

### `chimera_annotations.tsv`
Basic binary classification of contigs as chimeric or normal.

### `junction_details.tsv`
Precise breakpoint coordinates with source genome information.

### `pe_artifact_predictions.tsv`
Expected vs observed PE mapping artifacts for each chimeric junction.

### `difficulty_assessment.tsv`
Detection difficulty classification and recommended test cases.

### `ground_truth_summary.json`
Comprehensive metadata including all contig information, coverage patterns, and detection statistics.

## Expected PE Artifacts

| Chimera Type | Expected Artifact | Detection Signal |
|--------------|-------------------|------------------|
| Gap-mediated | Large insert sizes | Reads spanning N-gaps show insert size jumps |
| Orientation-flip | Improper orientations | FR pairs become RF/FF/RR at junctions |
| Distant-join | Extreme insert sizes | Massive insert sizes for cross-genome joins |

## Testing Your Detection Algorithm

1. **Start with easy cases**: Use contigs marked as "easy" detection difficulty
2. **Validate expected artifacts**: Check that your algorithm detects the predicted PE anomalies
3. **Systematic evaluation**: Test across all difficulty levels and chimera types
4. **Coverage independence**: Ensure detection doesn't rely on coverage patterns

## Performance Metrics

The simulator generates datasets with:
- **>85% mapping rate**: Realistic for challenging metagenomic data
- **Coassembly-like coverage**: 0.1x-500x range with log-normal distribution
- **Clear PE signals**: Detectable artifacts at all chimeric junctions
- **Balanced difficulty**: Mix of easy, medium, and hard detection cases

## Troubleshooting

### Low mapping rates
- Check that external mapping tools (BWA/minimap2) are installed
- Verify read simulation completed successfully
- Consider adjusting error rates in configuration

### No PE artifacts detected
- Ensure BAM file contains properly paired reads
- Check insert size distribution in mapping statistics
- Verify chimeric contigs were generated (check logs)

### Memory issues
- Reduce `total_contigs` in configuration
- Lower `coverage_target` for large datasets
- Use fewer `samples` in coverage model

## Citation

If you use this simulator in your research, please cite:
[Citation information to be added]

## License

[License information to be added]

## Contributing

Issues and pull requests welcome! Please see the GitHub repository for details.