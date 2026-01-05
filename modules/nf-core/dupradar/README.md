# nf-core/modules: dupradar

## Description

Assessment of duplication rates in RNA-Seq datasets using optimized featureCounts-based implementation with cloud storage compatibility.

## Usage

```nextflow
include { DUPRADAR } from './modules/nf-core/dupradar/main'

workflow {
    ch_bam = Channel.fromPath("*.bam").map { file -> 
        [[id: file.baseName, single_end: false, strandedness: 'unstranded'], file] 
    }
    ch_gtf = Channel.fromPath("genome.gtf").map { file -> [[:], file] }
    
    DUPRADAR(ch_bam, ch_gtf)
}
```

## Inputs

- **meta** (map): Groovy map containing sample information
  - Required keys: `id`, `single_end`, `strandedness`
  - Example: `[id: 'sample1', single_end: false, strandedness: 'forward']`
- **bam** (file): BAM file containing read alignments (duplicate-marked recommended)
- **gtf** (file): GTF file with genomic feature annotations

## Outputs

- **scatter2d** (file): `*_duprateExpDens.pdf` - 2D density scatter plot
- **boxplot** (file): `*_duprateExpBoxplot.pdf` - Boxplot of duplication vs expression
- **hist** (file): `*_expressionHist.pdf` - Expression histogram
- **dupmatrix** (file): `*_dupMatrix.txt` - Raw duplication matrix data
- **intercept_slope** (file): `*_intercept_slope.txt` - GLM model parameters
- **multiqc** (file): `*_mqc.txt` - MultiQC-compatible files
- **session_info** (file): `*.R_sessionInfo.log` - Session information
- **versions** (file): `versions.yml` - Software versions

## Parameters

- **ext.args** (string): Additional arguments
  - `--feature_type <type>`: Feature type to count (default: 'exon')
  - Example: `ext.args = '--feature_type CDS'`

### Meta Parameters

- **strandedness**: `'unstranded'`, `'forward'`, or `'reverse'`
- **single_end**: `true` for single-end, `false` for paired-end

## Technical Details

This implementation uses:
1. **featureCounts** (Subread): Direct BAM/GTF processing with/without duplicate reads
2. **Python analysis**: Duplication rate calculations and GLM fitting
3. **Plot generation**: Visualization using matplotlib

Performance: ~10-60 minutes, 4-8 GB memory, 2-4 cores.

## Troubleshooting

### Memory Issues
Increase memory allocation:
```nextflow
process {
    withName: 'DUPRADAR' {
        memory = '16.GB'
    }
}
```

### Container Issues
```bash
nextflow run pipeline.nf --profile docker,wave
```

## Integration with MultiQC

dupRadar outputs are automatically detected by MultiQC:
- `*_dup_intercept_mqc.txt`: General statistics table
- `*_duprateExpDensCurve_mqc.txt`: Line graph for the dupRadar plot

## Citations

- [dupRadar](https://doi.org/10.1186/s12859-016-1276-2): Sayols S, Scherzinger D, Klein H. dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data. BMC Bioinformatics. 2016;17:428.

- [Subread](https://doi.org/10.1093/nar/gkt214): Liao Y, Smyth GK, Shi W. The subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research. 2013;41(10):e108.

## Authors

- Original dupRadar module: [@pinin4fjords](https://github.com/pinin4fjords)
- Optimized implementation: [@edmundmiller](https://github.com/edmundmiller)
