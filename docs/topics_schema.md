# Nextflow Topics Schema for nf-core/rnaseq

## Overview
This document defines the topic naming schema for the nf-core/rnaseq pipeline using Nextflow 25.10+ topics feature.

## Topic Naming Convention

Topics follow a hierarchical naming scheme: `{category}/{subcategory}/{datatype}`

## Topic Taxonomy

### 1. Alignment Outputs

#### STAR Aligner
- `star/bam/genome` - Genome-aligned BAM files
- `star/bam/transcriptome` - Transcriptome-aligned BAM files
- `star/bam/sorted` - Coordinate-sorted genome BAM files
- `star/bam/unsorted` - Unsorted BAM files
- `star/index/bai` - BAM index files
- `star/index/csi` - CSI index files
- `star/logs/final` - STAR final alignment log
- `star/logs/out` - STAR stdout log
- `star/logs/progress` - STAR progress log
- `star/stats/samtools` - Samtools stats for STAR alignments
- `star/stats/flagstat` - Samtools flagstat for STAR alignments
- `star/stats/idxstats` - Samtools idxstats for STAR alignments
- `star/stats/percent_mapped` - Percent mapped reads from STAR
- `star/fastq/unmapped` - Unmapped reads from STAR
- `star/tab/splice_junctions` - Splice junction tables
- `star/tab/read_per_gene` - Read counts per gene

#### HISAT2 Aligner
- `hisat2/bam/genome` - Genome-aligned BAM files
- `hisat2/bam/sorted` - Coordinate-sorted BAM files
- `hisat2/index/bai` - BAM index files
- `hisat2/index/csi` - CSI index files
- `hisat2/logs/summary` - HISAT2 alignment summary
- `hisat2/stats/samtools` - Samtools stats for HISAT2 alignments
- `hisat2/stats/flagstat` - Samtools flagstat for HISAT2 alignments
- `hisat2/stats/idxstats` - Samtools idxstats for HISAT2 alignments
- `hisat2/fastq/unmapped` - Unmapped reads from HISAT2

### 2. UMI Deduplication Outputs

- `umi/dedup/bam/genome` - Deduplicated genome BAM files
- `umi/dedup/bam/transcriptome` - Deduplicated transcriptome BAM files
- `umi/dedup/index/bai` - Index for deduplicated BAMs
- `umi/dedup/index/csi` - CSI index for deduplicated BAMs
- `umi/dedup/logs` - UMI deduplication logs
- `umi/dedup/stats/samtools` - Samtools stats for deduplicated BAMs
- `umi/dedup/stats/flagstat` - Samtools flagstat for deduplicated BAMs
- `umi/dedup/stats/idxstats` - Samtools idxstats for deduplicated BAMs
- `umi/dedup/stats/edit_distance` - UMI edit distance statistics
- `umi/dedup/stats/per_umi` - Per-UMI statistics
- `umi/dedup/stats/per_position` - Per-position UMI statistics

### 3. Quantification Outputs

- `salmon/counts` - Salmon transcript quantification
- `salmon/logs` - Salmon quantification logs
- `kallisto/counts` - Kallisto transcript quantification
- `kallisto/logs` - Kallisto quantification logs
- `rsem/counts/genes` - RSEM gene counts
- `rsem/counts/transcripts` - RSEM transcript counts
- `rsem/logs` - RSEM quantification logs
- `featurecounts/counts` - featureCounts gene counts
- `featurecounts/summary` - featureCounts summary

### 4. Quality Control Outputs

- `qc/fastqc/raw` - FastQC reports for raw reads
- `qc/fastqc/trimmed` - FastQC reports for trimmed reads
- `qc/rseqc` - RSeQC output files
- `qc/qualimap` - Qualimap RNA-seq QC
- `qc/dupradar` - dupRadar duplication analysis
- `qc/preseq` - Preseq complexity curves
- `qc/deseq2` - DESeq2 QC plots

### 5. Read Processing Outputs

- `trimming/fastq` - Trimmed FASTQ files
- `trimming/logs` - Trimming logs (fastp, TrimGalore, etc.)
- `rrna_removal/fastq` - Reads after rRNA removal
- `rrna_removal/logs` - rRNA removal logs (SortMeRNA, BBSplit)

### 6. Coverage and Visualization

- `coverage/bedgraph/forward` - Forward strand coverage bedGraph
- `coverage/bedgraph/reverse` - Reverse strand coverage bedGraph
- `coverage/bigwig/forward` - Forward strand coverage bigWig
- `coverage/bigwig/reverse` - Reverse strand coverage bigWig

### 7. MultiQC Aggregation

- `multiqc/files` - All files to be included in MultiQC report
- `multiqc/report` - Final MultiQC HTML report

### 8. Versions and Metadata

- `versions` - All software version files (versions.yml)
- `metadata/trim_status` - Sample trim status
- `metadata/map_status` - Sample mapping status
- `metadata/strand_status` - Sample strandedness

### 9. Additional Analyses

- `stringtie/gtf` - StringTie assembled transcripts
- `stringtie/coverage` - StringTie coverage tables
- `kraken2/report` - Kraken2 classification reports
- `bracken/report` - Bracken abundance estimates

## Implementation Notes

1. **Universal `versions` topic**: All modules emitting `versions.yml` should publish to the single `versions` topic for centralized collection.

2. **MultiQC aggregation**: Instead of passing files through subworkflow chains, publish directly to `multiqc/files` topic.

3. **Conditional outputs**: Topics that are only populated when certain parameters are enabled (e.g., UMI processing) will simply remain empty - the output block can handle this gracefully.

4. **Backwards compatibility**: Modules modified to use topics will still work in pipelines that don't subscribe to those topics - the `emit` declarations are standard Nextflow outputs.

5. **Topic namespacing**: The hierarchical naming prevents collisions and makes the output block self-documenting.

## Migration Strategy

1. Start with the most-used outputs: BAM files, indices, and logs
2. Add `versions` topic to all modules
3. Add MultiQC-related topics
4. Gradually eliminate intermediate channel wiring in subworkflows
5. Simplify the main workflow emit block
6. Add output block to main.nf that subscribes to topics

## Expected Benefits

- **File size reduction**: Eliminate thousands of lines of channel wiring
- **Clarity**: Output routing is explicit in topic names
- **Maintainability**: Adding new outputs doesn't require modifying every subworkflow layer
- **Flexibility**: Easy to redirect outputs or add new publishing destinations
