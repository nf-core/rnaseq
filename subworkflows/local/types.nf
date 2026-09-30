// Record types shared across the pipeline's top-level scripts. Documentation only: nothing casts a
// record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { SamtoolsStatsFiles } from '../nf-core/bam_stats_samtools/types'
include { StarLogs           } from './align_star/types'
include { Bowtie2Logs        } from './align_bowtie2/types'
include { Hisat2Logs         } from '../nf-core/fastq_align_hisat2/types'
include { RsemMerge          } from '../nf-core/quantify_rsem/types'
include { BigwigFiles        } from '../nf-core/bedgraph_bedclip_bedgraphtobigwig/types'
include { StringtieAssembly  } from '../nf-core/bam_stringtie_merge/types'

record Sample {
    id:   String
    meta: Map
}

// Alignment stage result of STAR, Bowtie2 or HISAT2: the aligner-specific log group is set for
// the aligner that ran, and null for the others.
record AlignedSample {
    id:                String
    meta:              Map
    aligner:           String
    orig_bam:          List<Path>
    bam:               Path
    bai:               Path
    transcriptome_bam: Path?
    unmapped:          List<Path>?
    percent_mapped:    Float?
    samtools:          SamtoolsStatsFiles
    star:              StarLogs?
    hisat2:            Hisat2Logs?
    bowtie2:           Bowtie2Logs?
}

// Per-sample quantification of the alignment-based quantifier (RSEM or Salmon) or a pseudo-aligner.
record QuantSample {
    id:                String
    meta:              Map
    counts_gene:       Path?
    counts_transcript: Path?
    stat:              Path?
    quant_dir:         Path?
    log:               Path?
}

record RsemMergeResult {
    id:         String
    rsem_merge: RsemMerge
}

record Deseq2Results {
    rdata:        Path
    pca_vals:     Path
    plots_pdf:    Path
    sample_dists: Path
    size_factors: Path
    log:          Path
}

record ContaminantsKraken2 {
    report:                       Path?
    classified_reads_fastq:       Path?
    unclassified_reads_fastq:     Path?
    classified_reads_assignment:  Path?
}

record ContaminantsBracken {
    abundance: Path
    report:    Path
}

record ContaminantsSylph {
    profile: Path
}

record ContaminantsSylphtax {
    taxprof: Path
}

// Exactly one tool branch is populated per run.
record ContaminantsSample {
    id:        String
    meta:      Map
    kraken2:   ContaminantsKraken2?
    bracken:   ContaminantsBracken?
    sylph:     ContaminantsSylph?
    sylphtax:  ContaminantsSylphtax?
}

record StringtieSample {
    id:             String
    meta:           Map
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
    denovo:         StringtieAssembly?
}

record BigwigSample {
    id:      String
    meta:    Map
    combined: BigwigFiles?
    forward:  BigwigFiles?
    reverse:  BigwigFiles?
}

record RustqcSamtools {
    stats:    Path?
    flagstat: Path?
    idxstats: Path?
}

record RustqcPreseq {
    lc_extrap: Path?
}

record RustqcDupradar {
    scatter2d:       List<Path>
    boxplot:         List<Path>
    hist:            List<Path>
    dupmatrix:       Path?
    intercept_slope: Path?
    multiqc:         List<Path>
}

record RustqcFeaturecounts {
    counts:  Path?
    summary: Path?
}

record RustqcBiotype {
    tsv:  Path?
    mqc:  Path?
    rrna: Path?
}

record RustqcTin {
    txt: Path?
    xls: Path?
}

record RustqcInnerdistance {
    distance: Path?
    freq:     Path?
    mean:     Path?
    summary:  Path?
    plot:     List<Path>
    rscript:  Path?
}

record RustqcJunctionannotation {
    bed:          Path?
    interact_bed: Path?
    xls:          Path?
    log:          Path?
    plot:         List<Path>
    rscript:      Path?
}

record RustqcJunctionsaturation {
    summary: Path?
    plot:    List<Path>
    rscript: Path?
}

record RustqcReadduplication {
    seq_xls: Path?
    pos_xls: Path?
    plot:    List<Path>
    rscript: Path?
}

record RustqcRseqc {
    bamstat:            Path?
    inferexperiment:    Path?
    readdistribution:   Path?
    tin:                RustqcTin
    innerdistance:      RustqcInnerdistance
    junctionannotation: RustqcJunctionannotation
    junctionsaturation: RustqcJunctionsaturation
    readduplication:    RustqcReadduplication
}

record RustqcSample {
    id:            String
    meta:          Map
    samtools:      RustqcSamtools
    preseq:        RustqcPreseq
    dupradar:      RustqcDupradar
    featurecounts: RustqcFeaturecounts
    biotype:       RustqcBiotype
    rseqc:         RustqcRseqc
    qualimap:      Path?
    all_files:     List<Path>
}

// One FQ_LINT result. Raw, trimmed, BBSplit-filtered and rRNA-removed reads each publish through
// their own output because the four results share a basename.
record LintFile {
    id:   String
    file: Path
}

record PipelineInfo {
    versions: Path
}

// One row of samplesheet_with_bams.csv; field order is the CSV column order.
record SamplesheetRow {
    sample:            String
    fastq_1:           Path
    fastq_2:           Path?
    strandedness:      String
    seq_platform:      String?
    seq_center:        String?
    genome_bam:        String?
    percent_mapped:    Float?
    transcriptome_bam: String?
}
