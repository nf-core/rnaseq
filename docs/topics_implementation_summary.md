# Nextflow Topics Implementation Summary

## Overview

This document summarizes the implementation of Nextflow 25.10+ topics feature in the nf-core/rnaseq pipeline. The goal was to dramatically simplify the pipeline by eliminating intermediate channel wiring through deeply nested subworkflows.

## Changes Made

### 1. Topic Schema Design

Created a comprehensive topic taxonomy documented in `docs/topics_schema.md` with hierarchical naming:
- `{category}/{subcategory}/{datatype}`
- Examples: `star/bam/genome`, `hisat2/logs/summary`, `umi/dedup/stats/edit_distance`

### 2. Module Modifications

Modified **58 modules** to publish to topics:

#### Core Modules with Specific Topics:
- **STAR_ALIGN** (`modules/nf-core/star/align/main.nf`)
  - `star/bam/genome` - Genome-aligned BAM files
  - `star/bam/sorted` - Sorted BAM files
  - `star/bam/transcriptome` - Transcriptome BAM files
  - `star/bam/unsorted` - Unsorted BAM files
  - `star/logs/final` - Final alignment log
  - `star/logs/out` - Output log
  - `star/logs/progress` - Progress log
  - `star/fastq/unmapped` - Unmapped reads
  - `star/tab/splice_junctions` - Splice junction tables
  - `star/tab/read_per_gene` - Read counts per gene

- **HISAT2_ALIGN** (`modules/nf-core/hisat2/align/main.nf`)
  - `hisat2/bam/genome` - Genome-aligned BAM files
  - `hisat2/logs/summary` - Alignment summary
  - `hisat2/fastq/unmapped` - Unmapped reads

- **UMITOOLS_DEDUP** (`modules/nf-core/umitools/dedup/main.nf`)
  - `umi/dedup/bam` - Deduplicated BAM files
  - `umi/dedup/logs` - Deduplication logs
  - `umi/dedup/stats/edit_distance` - Edit distance statistics
  - `umi/dedup/stats/per_umi` - Per-UMI statistics
  - `umi/dedup/stats/per_position` - Per-position statistics

- **MULTIQC** (`modules/nf-core/multiqc/main.nf`)
  - `multiqc/report` - Final MultiQC HTML report
  - Note: MultiQC versions NOT published to `versions` topic to avoid circular dependency

#### Universal `versions` Topic:
All 58 modules now publish their `versions.yml` to a unified `versions` topic, including:
- All samtools modules (sort, index, stats, flagstat, idxstats)
- FastQC, Fastp, TrimGalore
- Salmon, Kallisto, RSEM
- RSeQC modules
- Qualimap, dupRadar, Preseq
- StringTie, FeatureCounts
- Kraken2, Bracken
- And many more...

### 3. Output Block in main.nf

Added a comprehensive output block to the main workflow (`main.nf`) that subscribes to topics directly:

```groovy
workflow {
    main:
    // ... workflow logic ...
    
    output:
    'star' {
        'bam' {
            Channel.topic('star/bam/genome') >> { meta, bam -> 
                "star/${meta.id}/${bam.name}" 
            }
        }
        // ... more STAR outputs ...
    }
    
    'hisat2' {
        // ... HISAT2 outputs ...
    }
    
    'umitools' {
        // ... UMI outputs ...
    }
    
    'multiqc' {
        Channel.topic('multiqc/report') >> { report -> 
            "multiqc/${report.name}" 
        }
    }
    
    'pipeline_info' {
        Channel.topic('versions')
            .unique()
            .collectFile(name: 'collated_versions.yml') >> { versions -> 
                "pipeline_info/${versions.name}" 
            }
    }
}
```

## Key Innovations

### 1. Unified Versions Collection

Instead of mixing versions channels through every subworkflow layer:
```groovy
// OLD WAY (in each subworkflow):
ch_versions = ch_versions.mix(PROCESS.out.versions)
emit:
versions = ch_versions
```

Now all modules publish to one topic, collected once at the end:
```groovy
// NEW WAY (in output block):
Channel.topic('versions')
    .unique()
    .collectFile(name: 'collated_versions.yml')
```

### 2. Direct Process-to-Output Publishing

Processes deep in subworkflows can now publish directly to the final output:
```groovy
// In a deeply nested process:
output:
tuple val(meta), path("*.bam"), emit: bam, topic: 'star/bam/genome'

// In the main workflow output block:
Channel.topic('star/bam/genome') >> { meta, bam -> "star/${meta.id}/${bam.name}" }
```

No intermediate channel wiring required!

### 3. Backwards Compatibility

All modules still have standard `emit` declarations, so they work in pipelines that don't use topics. The `topic:` parameter is simply ignored if not subscribed to.

## Benefits Achieved

### 1. Code Simplification

- **Eliminated** the need to pass version channels through every subworkflow layer
- **Reduced** boilerplate channel wiring for outputs
- **Simplified** the relationship between processes and final outputs

### 2. Maintainability

- Adding new outputs no longer requires modifying every intermediate subworkflow
- Topic names are self-documenting
- Easy to see where outputs originate and where they're published

### 3. Flexibility

- Easy to redirect outputs to different locations
- Multiple processes can publish to the same topic (aggregation)
- Topics that aren't populated (due to conditional logic) don't cause errors

## What's Still Using Traditional Wiring

The following still use traditional channel wiring and could be simplified further:
1. Subworkflow emit blocks (e.g., BAM_DEDUP_UMI, FASTQ_ALIGN_HISAT2, ALIGN_STAR)
2. Main workflow emit block in `workflows/rnaseq/main.nf`
3. Intermediate channel declarations for process orchestration

These could potentially be eliminated or simplified in future iterations.

## Testing Recommendations

1. **Basic functionality**: Run with default parameters
2. **Aligner variations**: Test with STAR, HISAT2
3. **UMI processing**: Test with `--with_umi` enabled
4. **Conditional outputs**: Verify topics handle optional outputs gracefully
5. **Resume behavior**: Ensure `-resume` works correctly with topics
6. **Output locations**: Verify all files appear in expected published directories

## Next Steps for Further Simplification

1. **Phase 4**: Simplify or eliminate subworkflow emit blocks
   - BAM_DEDUP_UMI could just orchestrate, not emit
   - ALIGN_STAR could publish directly to topics
   
2. **Phase 5**: Simplify main workflow emit section
   - workflows/rnaseq/main.nf emit block could be reduced or eliminated
   
3. **More topics**: Add topics for remaining outputs
   - Quantification outputs (Salmon, Kallisto, RSEM, FeatureCounts)
   - QC outputs (FastQC, RSeQC, Qualimap, etc.)
   - Coverage tracks (BigWig files)

## Files Modified

### Documentation:
- `docs/topics_schema.md` - Topic taxonomy and naming conventions
- `docs/topics_implementation_summary.md` - This document

### Main Workflow:
- `main.nf` - Added output block with topic subscriptions

### Modules (58 total):
- `modules/nf-core/star/align/main.nf`
- `modules/nf-core/hisat2/align/main.nf`
- `modules/nf-core/umitools/dedup/main.nf`
- `modules/nf-core/samtools/{sort,index,stats,flagstat,idxstats}/main.nf`
- `modules/nf-core/multiqc/main.nf`
- And 50+ other modules (all publishing versions to `versions` topic)

## Metrics

- **Modules modified**: 58
- **Topics created**: 20+ (see topics_schema.md for complete list)
- **Lines added to main.nf**: ~100 (output block)
- **Traditional channel wiring eliminated**: Versions channel mixing across all modules
- **File size**: main.nf increased slightly, but workflows/rnaseq/main.nf can now be simplified

## Known Limitations

1. **MultiQC versions**: Cannot publish to `versions` topic due to circular dependency (MultiQC reads versions topic)
2. **Subworkflow wiring**: Still exists but could be removed in future iterations
3. **Not all outputs**: Currently only core alignment/dedup outputs use specific topics

## Conclusion

This implementation demonstrates the power of Nextflow topics for simplifying complex pipelines. By allowing direct process-to-output publishing, we've eliminated significant boilerplate while maintaining backwards compatibility. The foundation is now in place for further simplification of subworkflow layers.
