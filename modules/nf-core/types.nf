// Record types shared by modules, subworkflows and the pipeline.
// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
//
// A module takes its record input as `sample: <Type>` with a type from this file (bound as one
// variable, so scripts and task.ext closures read sample.meta, sample.bam ...). The *Result types
// describe module outputs, for annotating the channels a process call returns.

// ============================================================================
// Shared vocabulary
// ============================================================================

record Sample {
    id:   String
    meta: Map
}

// FASTQ files of one sample
record ReadsInput {
    id:    String
    meta:  Map
    reads: List<Path>
}

// A BAM with its index
record BamBaiInput {
    id:   String
    meta: Map
    bam:  Path
    bai:  Path
}

// A BAM without an index
record BamInput {
    id:   String
    meta: Map
    bam:  Path
}

// One FASTA
record FastaInput {
    id:    String
    meta:  Map
    fasta: Path
}

// A FASTA with its annotation
record FastaGtfInput {
    id:    String
    meta:  Map
    fasta: Path
    gtf:   Path
}

// One GTF
record GtfInput {
    id:   String
    meta: Map
    gtf:  Path
}

// Quantification outputs of one sample or run
record QuantsInput {
    id:     String
    meta:   Map
    quants: List<Path>
}

// One compressed archive
record ArchiveInput {
    id:      String
    meta:    Map
    archive: Path
}

// One bedGraph
record BedgraphInput {
    id:       String
    meta:     Map
    bedgraph: Path
}

// ============================================================================
// bowtie2
// ============================================================================

record Bowtie2Logs {
    log: Path
}

// ============================================================================
// bracken
// ============================================================================

record BrackenResult {
    id:        String
    meta:      Map
    abundance: Path
    report:    Path
}

// ============================================================================
// dupradar
// ============================================================================

record BamQcDupradar {
    id:              String
    meta:            Map
    scatter2d:       Path
    boxplot:         Path
    hist:            Path
    dupmatrix:       Path
    intercept_slope: Path
    multiqc:         List<Path>
    session_info:    Path
}

// ============================================================================
// hisat2
// ============================================================================

record Hisat2Logs {
    summary: Path
}

// ============================================================================
// kraken2
// ============================================================================

record Kraken2Result {
    id:                          String
    meta:                        Map
    report:                      Path
    classified_reads_fastq:      Set<Path>
    unclassified_reads_fastq:    Set<Path>
    classified_reads_assignment: Path?
}

// ============================================================================
// preseq
// ============================================================================

record BamQcPreseq {
    id:        String
    meta:      Map
    lc_extrap: Path
    log:       Path
}

// ============================================================================
// rseqc
// ============================================================================

record RseqcInnerDistance {
    id:       String
    meta:     Map
    distance: Path
    freq:     Path?
    mean:     Path?
    pdf:      Path?
    rscript:  Path?
}

record RseqcJunctionAnnotation {
    id:           String
    meta:         Map
    bed:          Path?
    interact_bed: Path?
    xls:          Path
    pdf:          Path?
    events_pdf:   Path?
    rscript:      Path
    log:          Path
}

record RseqcJunctionSaturation {
    id:      String
    meta:    Map
    pdf:     Path
    rscript: Path
}

record RseqcReadDuplication {
    id:      String
    meta:    Map
    seq_xls: Path
    pos_xls: Path
    pdf:     Path
    rscript: Path
}

record RseqcTin {
    id:   String
    meta: Map
    txt:  Path
    xls:  Path
}

// ============================================================================
// samtools
// ============================================================================

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record SamtoolsStats {
    id:       String
    meta:     Map
    samtools: SamtoolsStatsFiles
}

record RawBams {
    id:       String
    meta:     Map
    raw_bams: List<Path>
}

// ============================================================================
// star
// ============================================================================

record StarLogs {
    log_final:    Path
    log_out:      Path
    log_progress: Path
    tab:          List<Path>
}

// ============================================================================
// stringtie
// ============================================================================

record StringtieInput {
    id:    String
    meta:  Map
    bam:   Path?
    lrbam: Path?
}

record StringtieAssembly {
    id:             String
    meta:           Map
    transcript_gtf: Path
    abundance:      Path
    coverage_gtf:   Path?
    ballgown:       List<Path>?
}

// ============================================================================
// subread
// ============================================================================

record BamQcFeaturecounts {
    id:      String
    meta:    Map
    counts:  Path
    summary: Path
}

// ============================================================================
// Subworkflow and pipeline records
// ============================================================================

record Bowtie2Aligned {
    id:                String
    meta:              Map
    aligner:           String
    raw_bams:          List<Path>
    transcriptome_bam: Path
    unmapped:          List<Path>?
    percent_mapped:    Float?
    bowtie2:           Bowtie2Logs
    bam:               Path
    bai:               Path
    samtools:          SamtoolsStatsFiles
}

// StarAlignResult plus aligner, percent_mapped and the sorted BAM with its samtools stats.
record StarAligned {
    id:                 String
    meta:               Map
    aligner:            String
    raw_bams:           List<Path>
    bam_sorted:         Path?
    bam_sorted_aligned: Path?
    bam_unsorted:       Path?
    transcriptome_bam:  Path?
    unmapped:           List<Path>
    sam:                Path?
    junction:           Path?
    spl_junc_tab:       Path?
    read_per_gene_tab:  Path?
    wig:                List<Path>
    bedgraph:           List<Path>
    orig_bai:           Path?
    qc_metrics:         Path?
    duplicate_metrics:  Path?
    percent_mapped:     Float?
    star:               StarLogs
    bam:                Path
    bai:                Path
    samtools:           SamtoolsStatsFiles
}

// MultiQC input files one stage contributes for one sample
record MultiqcFiles {
    id:    String
    files: List<Path>
}

// Alignment stage result of STAR, Bowtie2 or HISAT2: the aligner-specific log group is set for
// the aligner that ran, and null for the others.
record AlignedSample {
    id:                String
    meta:              Map
    aligner:           String
    raw_bams:          List<Path>
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

// A coordinate-sorted BAM, from an aligner, a stage that rewrites it, or a pre-aligned samplesheet
// entry. Stage records that add fields of their own (metrics, dedup logs, aligner logs) keep these
// names and types. percent_mapped is null when the samplesheet does not provide it, and
// transcriptome_bam is set only where a transcriptome-aligned BAM exists.
record Bam {
    id:                String
    meta:              Map
    bam:               Path
    bai:               Path?
    samtools:          SamtoolsStatsFiles?
    percent_mapped:    Float?
    transcriptome_bam: Path?
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

record PipelineInfo {
    versions: Path
}

// One independently-produced, optionally-user-supplied publish artifact. Used to stream
// genome reference/index/intermediate files as records instead of fusing them into one
// wide row via a fake join - each producer mixes in its own record, or none at all.
record GenomeArtifact {
    kind: String
    file: Path
}

record UmicollapseDedupBam {
    id:       String
    meta:     Map
    bam:      Path
    bai:      Path
    log:      Path
    samtools: SamtoolsStatsFiles
}

// umi_tools writes all three stats tables together, so the tsv fields are either all set or all null.
record UmitoolsDedupBam {
    id:                   String
    meta:                 Map
    bam:                  Path
    bai:                  Path
    log:                  Path
    tsv_edit_distance:    Path?
    tsv_per_umi:          Path?
    tsv_umi_per_position: Path?
    samtools:             SamtoolsStatsFiles
}

record UmitoolsDedupStats {
    edit_distance:    Path
    per_umi:          Path
    umi_per_position: Path
}

record UmiDedupTranscriptome {
    transcriptome_bam:      Path
    dedup_bam:              Path
    sorted_bam:             Path
    sorted_bam_index:       Path
    filtered_bam:           Path?
    samtools:               SamtoolsStatsFiles
    tsv:                    UmitoolsDedupStats?
    coord_sorted_bam:       Path?
    coord_sorted_bam_index: Path?
    coord_sorted_samtools:  SamtoolsStatsFiles?
}

record UmiDedupBam {
    id:                       String
    meta:                     Map
    bam:                      Path
    bai:                      Path
    samtools:                 SamtoolsStatsFiles
    genomic_dedup_log:        Path
    transcriptomic_dedup_log: Path?
    prepare_for_rsem_log:     Path?
    transcriptome_bam:        Path?
    transcriptome:            UmiDedupTranscriptome?
    tsv:                      UmitoolsDedupStats?
}

record MarkdupBam {
    id:       String
    meta:     Map
    bam:      Path?
    cram:     Path?
    bai:      Path
    metrics:  Path
    samtools: SamtoolsStatsFiles?
}

record BamQcBiotype {
    tsv:  Path
    rrna: Path
}

record BamQcRnaseq {
    id:            String
    meta:          Map
    preseq:        BamQcPreseq?
    featurecounts: BamQcFeaturecounts?
    biotype:       BamQcBiotype?
    qualimap:      Path?
    dupradar:      BamQcDupradar?
    rseqc:         Rseqc?
    mqc_files:     List<Path>
}

record Rseqc {
    id:                 String
    meta:               Map
    bamstat:            Path?
    inferexperiment:    Path?
    innerdistance:      RseqcInnerDistance?
    junctionannotation: RseqcJunctionAnnotation?
    junctionsaturation: RseqcJunctionSaturation?
    readdistribution:   Path?
    readduplication:    RseqcReadDuplication?
    tin:                RseqcTin?
}

// One record for the whole run: the merged GTF and the per-sample assemblies that fed it.
record StringtieMerged {
    id:         String
    meta:       Map
    gtf:        List<Path>
    assemblies: Bag<StringtieAssembly>
    merged_gtf: Path
}

record BigwigFiles {
    id:       String
    meta:     Map
    bedgraph: Path
    bigwig:   Path
}

record Hisat2Aligned {
    id:       String
    meta:     Map
    aligner:  String
    raw_bams: List<Path>
    unmapped: List<Path>
    hisat2:   Hisat2Logs
    bam:      Path
    bai:      Path
    samtools: SamtoolsStatsFiles
}

record FastpReads {
    id:            String
    meta:          Map
    reads:         List<Path>
    adapter_fasta: Path?
}

record FastpTrim {
    html:         Path
    json:         Path
    log:          Path
    reads_fail:   List<Path>?
    reads_merged: Path?
}

record UmitoolsExtractFiles {
    log:   Path
    reads: List<Path>
}

// `reads` and `meta` are the reads handed to the next stage: UMI-extracted and
// with R1 or R2 discarded when requested, then trimmed. `reads` is null for
// samples below min_trimmed_reads or without any trimmed reads.
record FastqFastqcUmitoolsFastp {
    id:                String
    meta:              Map
    reads:             List<Path>?
    fastqc_raw_html:   List<Path>?
    fastqc_raw_zip:    List<Path>?
    fastqc_trim_html:  List<Path>?
    fastqc_trim_zip:   List<Path>?
    umi:               UmitoolsExtractFiles?
    trim:              FastpTrim?
    adapter_seq:       String?
    num_trimmed_reads: Long?
}

record TrimgaloreTrim {
    html:     List<Path>?
    zip:      List<Path>?
    log:      List<Path>?
    json:     List<Path>?
    unpaired: List<Path>?
}

// `reads` and `meta` are the reads handed to the next stage: UMI-extracted and
// with R1 or R2 discarded when requested, then trimmed. `reads` is null for
// samples below min_trimmed_reads. TrimgaloreTrim.html and zip are FastQC
// reports on the trimmed reads.
record FastqFastqcUmitoolsTrimgalore {
    id:                String
    meta:              Map
    reads:             List<Path>?
    fastqc_raw_html:   List<Path>?
    fastqc_raw_zip:    List<Path>?
    umi:               UmitoolsExtractFiles?
    trim:              TrimgaloreTrim?
    num_trimmed_reads: Float?
}

// TrimGalore and fastp report the same kinds of file with different
// cardinality, so every multi-file field is a list and tool-specific fields
// are nullable. TrimGalore's own html/zip are FastQC reports on the trimmed
// reads and are held in fastqc.trim_html/trim_zip.
record PreprocessedFastqc {
    raw_html:      List<Path>?
    raw_zip:       List<Path>?
    trim_html:     List<Path>?
    trim_zip:      List<Path>?
    filtered_html: List<Path>?
    filtered_zip:  List<Path>?
}

record PreprocessedTrim {
    html:         List<Path>?
    log:          List<Path>?
    json:         List<Path>?
    unpaired:     List<Path>?
    reads_fail:   List<Path>?
    reads_merged: Path?
}

record PreprocessedUmi {
    log:   Path
    reads: List<Path>
}

record PreprocessedBbsplit {
    stats:              Path
    primary_reads:      List<Path>?
    other_genome_reads: List<Path>?
}

record PreprocessedLint {
    raw:     Path?
    trimmed: Path?
    bbsplit: Path?
    ribo:    Path?
}

// The tool logs are set only for the selected tool.
record PreprocessedRrna {
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
}

// Run-level references, each set only when it is built here.
record RrnaReferences {
    sortmerna_index:  Path?
    bowtie2_index:    Path?
    seqkit_prefixed:  List<Path>?
    seqkit_converted: List<Path>?
}

// Samples that fail min_trimmed_reads keep their record with null reads and
// reads_trimmed, and a meta that lacks the inferred strandedness.
record FastqQcTrimFilterSetstrandedness {
    id:                String
    meta:              Map
    reads:             List<Path>?
    reads_cat:         List<Path>
    reads_trimmed:     List<Path>?
    num_trimmed_reads: Long?  // Float for TrimGalore, Long for fastp
    fastqc:            PreprocessedFastqc?
    trim:              PreprocessedTrim?
    umi:               PreprocessedUmi?
    bbsplit:           PreprocessedBbsplit?
    lint:              PreprocessedLint?
    rrna:              PreprocessedRrna?
}

// `reads` is null for samples with no reads left after rRNA removal. The tool
// logs are set only for the selected tool.
record FastqRemoveRrna {
    id:               String
    meta:             Map
    reads:            List<Path>?
    sortmerna_log:    Path?
    ribodetector_log: Path?
    seqkit_stats:     Path?
    bowtie2_log:      Path?
}

// reads are the sub-sampled FASTQ files.
record SalmonSubsampled {
    id:                String
    meta:              Map
    reads:             List<Path>
    quant_dir:         Path
    json_info:         Path?
    lib_format_counts: Path?
}

record QuantMerged {
    id:                        String
    meta:                      Map
    tpm_gene:                  Path
    counts_gene:               Path
    lengths_gene:              Path
    counts_gene_length_scaled: Path
    counts_gene_scaled:        Path
    tpm_transcript:            Path
    counts_transcript:         Path
    lengths_transcript:        Path
    tx2gene:                   Path
    tx2gene_augmented:         Path
    merged_gene_rds:           Path?
    merged_transcript_rds:     Path?
}
