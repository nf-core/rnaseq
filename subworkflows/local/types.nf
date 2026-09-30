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

// One FQ_LINT result. Raw, trimmed, BBSplit-filtered and rRNA-removed reads each publish through
// their own output because the four results share a basename.
record LintFile {
    id:   String
    file: Path
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

// A genome-aligned BAM with its index, from an aligner or a pre-aligned samplesheet entry.
// percent_mapped is null when the samplesheet does not provide it.
record GenomeBam {
    id:             String
    meta:           Map
    bam:            Path
    bai:            Path
    percent_mapped: Float?
}

// CUSTOM_RSEMMERGECOUNTS outputs, a single 'all_samples' row
record RsemMergeSample {
    id:         String
    meta:       Map
    rsem_merge: RsemMerge
}

record Kraken2Result {
    id:                          String
    meta:                        Map
    report:                      Path
    classified_reads_fastq:      List<Path>
    unclassified_reads_fastq:    List<Path>
    classified_reads_assignment: Path?
}

record BrackenResult {
    id:        String
    meta:      Map
    abundance: Path
    report:    Path
}

record SylphProfile {
    profile: Path
}

record SylphtaxTaxprof {
    taxprof: Path
}

// Exactly one screening tool branch is populated per run: kraken2 (plus bracken), or sylph (plus sylphtax
// when the profile was not empty).
record Contaminants {
    id:       String
    meta:     Map
    kraken2:  Kraken2Result?
    bracken:  BrackenResult?
    sylph:    SylphProfile?
    sylphtax: SylphtaxTaxprof?
}

// StringTie assembly with the per-sample de novo assembly that fed the merged GTF (--stringtie_ignore_gtf only)
record StringtieSample {
    id:             String
    meta:           Map
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
    denovo:         StringtieAssembly?
}

// Coverage tracks; forward and reverse are null for unstranded samples
record BigwigSample {
    id:       String
    meta:     Map
    combined: BigwigFiles
    forward:  BigwigFiles?
    reverse:  BigwigFiles?
}

record Deseq2Qc {
    id:            String
    meta:          Map
    rdata:         Path?
    pca_vals:      Path?
    plots_pdf:     Path?
    sample_dists:  Path?
    size_factors:  Path?
    log:           Path?
    pca_multiqc:   Path?
    dists_multiqc: Path?
}

record PipelineInfo {
    versions: Path
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
    scatter2d:       Set<Path>
    boxplot:         Set<Path>
    hist:            Set<Path>
    dupmatrix:       Path?
    intercept_slope: Path?
    multiqc:         Set<Path>
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
    plot:     Set<Path>
    rscript:  Path?
}

record RustqcJunctionannotation {
    bed:          Path?
    interact_bed: Path?
    xls:          Path?
    log:          Path?
    plot:         Set<Path>
    rscript:      Path?
}

record RustqcJunctionsaturation {
    summary: Path?
    plot:    Set<Path>
    rscript: Path?
}

record RustqcReadduplication {
    seq_xls: Path?
    pos_xls: Path?
    plot:    Set<Path>
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

record RustqcResult {
    id:            String
    meta:          Map
    samtools:      RustqcSamtools
    preseq:        RustqcPreseq
    dupradar:      RustqcDupradar
    featurecounts: RustqcFeaturecounts
    biotype:       RustqcBiotype
    rseqc:         RustqcRseqc
    qualimap:      Path?
    all_files:     Set<Path>
}
