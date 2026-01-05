# nf-core/modules: dupradar

## Description

Assessment of duplication rates in RNA-Seq datasets using [dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html).

This module provides two modes of operation:

1. **Traditional Mode** (default): Uses the original R-based dupRadar implementation
2. **Fast Mode**: High-performance implementation using native tools, optimized for cloud storage systems like Fusion S3

## Usage

### Basic Usage

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

### Enabling Fast Mode

#### Method 1: Process Configuration
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = true
        
        // Optional: Reduce resource requirements for fast mode
        memory = '8.GB'
        cpus = 2
        time = '2.h'
    }
}
```

#### Method 2: Conditional Enabling
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = { meta ->
            // Enable for specific samples or conditions
            meta.id.contains('large') || meta.reads_count > 100000000
        }
    }
}
```

#### Method 3: Configuration File
Create a `dupradar_fast.config` file:
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = true
        ext.args = '--feature_type exon'
    }
}
```

Then run: `nextflow run pipeline.nf -c dupradar_fast.config`

## Inputs

- **meta** (map): Groovy map containing sample information
  - Required keys: `id`, `single_end`, `strandedness`
  - Example: `[id: 'sample1', single_end: false, strandedness: 'forward']`
- **bam** (file): BAM file containing read alignments (duplicate-marked recommended)
- **gtf** (file): GTF file with genomic feature annotations

## Outputs

All outputs are identical between traditional and fast modes:

- **scatter2d** (file): `*_duprateExpDens.pdf` - 2D density scatter plot
- **boxplot** (file): `*_duprateExpBoxplot.pdf` - Boxplot of duplication vs expression
- **hist** (file): `*_expressionHist.pdf` - Expression histogram
- **dupmatrix** (file): `*_dupMatrix.txt` - Raw duplication matrix data
- **intercept_slope** (file): `*_intercept_slope.txt` - GLM model parameters
- **multiqc** (file): `*_mqc.txt` - MultiQC-compatible files
- **session_info** (file): `*.R_sessionInfo.log` - Session information
- **versions** (file): `versions.yml` - Software versions

## Parameters

### Process-specific Parameters

- **ext.use_fast_dupradar** (boolean): Enable fast mode (default: false)
- **ext.args** (string): Additional arguments
  - `--feature_type <type>`: Feature type to count (default: 'exon')
  - Example: `ext.args = '--feature_type CDS'`

### Meta Parameters

- **strandedness**: `'unstranded'`, `'forward'`, or `'reverse'`
- **single_end**: `true` for single-end, `false` for paired-end

## Mode Comparison

| Feature | Traditional Mode | Fast Mode |
|---------|------------------|-----------|
| **Runtime** | 2-12 hours | 10-60 minutes |
| **Memory** | 16-32 GB | 4-8 GB |
| **CPU** | 4-8 cores | 2-4 cores |
| **Storage I/O** | High latency sensitive | Optimized patterns |
| **Implementation** | R/dupRadar package | Native tools + minimal R |
| **Plot Quality** | Full dupRadar plots | Identical using dupRadar functions |
| **MultiQC Output** | ✅ Identical | ✅ Identical |
| **Fusion S3 Performance** | ⚠️ Slow (10+ hours) | ✅ Fast (30-60 min) |

## Fast Mode Technical Details

### Architecture
Fast mode uses a modular approach with specialized scripts:

1. **featureCounts** (native binary): Direct BAM/GTF processing
2. **duplication_rates.awk**: High-speed duplication calculations
3. **intercept_and_slope.py**: GLM model fitting
4. **generate_multiqc.py**: MultiQC file generation
5. **generate_plots.r**: Plot generation using original dupRadar functions

### Performance Optimization
- Bypasses R/Rsubread wrapper overhead
- Uses native binaries with optimal I/O patterns
- Eliminates memory-intensive R data structures
- Streams data processing where possible

### Compatibility Guarantee
- **100% MultiQC compatibility**: Identical output formats
- **Identical calculations**: Same algorithms as dupRadar
- **Same visualizations**: Uses extracted dupRadar plotting functions
- **API compatibility**: Drop-in replacement

## Troubleshooting

### Common Issues

#### 1. Fast Mode Not Working
```bash
Error: command not found: duplication_rates.awk
```
**Solution**: Ensure the bin scripts are executable and in PATH. This usually indicates a container/conda environment issue.

#### 2. Memory Issues in Traditional Mode
```bash
Error: Cannot allocate vector of size X GB
```
**Solution**: Switch to fast mode or increase memory allocation:
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = true  // or increase memory
        memory = '32.GB'
    }
}
```

#### 3. Fusion S3 Timeouts
```bash
Process DUPRADAR terminated with exit status (124) - timeout
```
**Solution**: Enable fast mode specifically for Fusion environments:
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = true
        time = '2.h'  // Reduced from default 12.h
    }
}
```

#### 4. Container Issues
```bash
Unable to find image 'community.wave.seqera.io/...'
```
**Solution**: Use Wave profile or pre-built containers:
```bash
nextflow run pipeline.nf --profile docker,wave
```

### Performance Debugging

Enable verbose logging to diagnose performance issues:

```nextflow
process {
    withName: 'DUPRADAR' {
        debug = true
        beforeScript = 'echo "Starting dupRadar analysis at $(date)"'
        afterScript = 'echo "Completed dupRadar analysis at $(date)"'
    }
}
```

## Examples

### Example 1: Basic Usage with Fast Mode
```nextflow
#!/usr/bin/env nextflow

include { DUPRADAR } from './modules/nf-core/dupradar/main'

workflow {
    // Input channels
    ch_bam = Channel.of([
        [[id: 'sample1', single_end: false, strandedness: 'forward'], 
         file('sample1.markdup.bam')]
    ])
    ch_gtf = Channel.of([[[:], file('genes.gtf')]])
    
    // Run dupRadar in fast mode
    DUPRADAR(ch_bam, ch_gtf)
    
    // Access outputs
    DUPRADAR.out.scatter2d.view { "Scatter plot: ${it[1]}" }
    DUPRADAR.out.multiqc.view { "MultiQC files: ${it[1]}" }
}
```

### Example 2: Conditional Fast Mode Based on File Size
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = { meta ->
            // Enable fast mode for files larger than 1GB
            def bamFile = task.inputFiles.find { it.name.endsWith('.bam') }
            return bamFile?.size() > 1_000_000_000
        }
    }
}
```

### Example 3: Different Configurations for Different Samples
```nextflow
process {
    withName: 'DUPRADAR' {
        ext.use_fast_dupradar = { meta ->
            switch(meta.sample_type) {
                case 'large_dataset': return true
                case 'fusion_storage': return true  
                default: return false
            }
        }
        
        memory = { meta ->
            task.ext.use_fast_dupradar ? '8.GB' : '16.GB'
        }
        
        time = { meta ->
            task.ext.use_fast_dupradar ? '2.h' : '8.h'
        }
    }
}
```

## Integration with MultiQC

dupRadar outputs are automatically detected by MultiQC. The generated files include:

- `*_dup_intercept_mqc.txt`: General statistics table
- `*_duprateExpDensCurve_mqc.txt`: Line graph for the dupRadar plot

No additional MultiQC configuration is required.

## Citations

If you use this module, please cite:

- [dupRadar](https://doi.org/10.1186/s12859-016-1276-2): Sayols S, Scherzinger D, Klein H. dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data. BMC Bioinformatics. 2016;17:428.

- [Subread](https://doi.org/10.1093/nar/gkt214): Liao Y, Smyth GK, Shi W. The subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research. 2013;41(10):e108.

## Authors

- Original dupRadar module: [@pinin4fjords](https://github.com/pinin4fjords)
- Fast mode implementation: [@edmundmiller](https://github.com/edmundmiller)

## Version History

- v1.0.0: Original dupRadar implementation
- v2.0.0: Added fast mode for improved performance with cloud storage systems