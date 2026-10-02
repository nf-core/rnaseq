nextflow.enable.types = true

//
// Build or load aligner / pseudo-aligner / filtering indices (STAR, RSEM, HISAT2, Bowtie2, Salmon, Kallisto, BBSplit, SortMeRNA)
//

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_BOWTIE2_RRNA_INDEX } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_BOWTIE2_INDEX      } from '../../../modules/nf-core/untar'

include { BOWTIE2_BUILD                     } from '../../../modules/nf-core/bowtie2/build'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_GENOMEGENERATE as PARABRICKS_STARGENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { KALLISTO_INDEX                    } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { SENTIEON_RSEMPREPAREREFERENCE as SENTIEON_RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/sentieon/rsempreparereference'

include { STAR_GENOMEPARAMS_UPGRADE         } from '../../../modules/local/star_genomeparams_upgrade'

include { taskOutputOrNull                  } from '../utils_nfcore_rnaseq_pipeline'
include { GenomeArtifact } from '../../../modules/nf-core/types'
include { BbmapBbsplitResult } from '../../../modules/nf-core/bbmap/bbsplit/main'
include { KallistoIndexResult } from '../../../modules/nf-core/kallisto/index/main'
include { SalmonIndexResult } from '../../../modules/nf-core/salmon/index/main'
include { SortmernaResult } from '../../../modules/nf-core/sortmerna/main'
include { StarGenomegenerateResult } from '../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME_INDICES {

    take:
    ch_fasta: Value<Path?>                      // genome.fasta - emitted from PREPARE_GENOME_REFERENCES
    ch_gtf: Value<Path?>                        // genome.gtf - emitted from PREPARE_GENOME_REFERENCES
    ch_transcript_fasta: Value<Path?>           // transcript.fasta - emitted from PREPARE_GENOME_REFERENCES
    ch_rrna_fastas: Channel<Path>               // rRNA fastas - emitted from PREPARE_GENOME_REFERENCES
    fasta_provided: Boolean                     // whether a genome FASTA was provided
    splicesites: String?                        // file: /path/to/splicesites.txt
    bbsplit_fasta_list: String?                 // file: /path/to/bbsplit_fasta_list.txt
    star_index: String?                         // directory: /path/to/star/index/
    rsem_index: String?                         // directory: /path/to/rsem/index/
    salmon_index: String?                       // directory: /path/to/salmon/index/
    kallisto_index: String?                     // directory: /path/to/kallisto/index/
    hisat2_index: String?                       // directory: /path/to/hisat2/index/
    bowtie2_index: String?                      // directory: /path/to/bowtie2/index/
    bbsplit_index: String?                      // directory: /path/to/bbsplit/index/
    sortmerna_index: String?                    // directory: /path/to/sortmerna/index/
    bowtie2_rrna_index: String?                 // directory: /path/to/bowtie2/index/
    aligner: String?                            // Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem', 'hisat2', and 'bowtie2_salmon'
    pseudo_aligner: String?                     // Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_bbsplit: Boolean                       // Skip BBSplit for removal of non-reference genome reads
    ribo_removal_tool: String?                  // Tool for rRNA removal - 'sortmerna', 'ribodetector', or 'bowtie2' (null if skip)
    skip_alignment: Boolean                     // Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment: Boolean              // Skip all of the pseudoalignment-based processes within the pipeline
    use_sentieon_star: Boolean                  // whether to use sentieon STAR version
    use_parabricks_star: Boolean                // whether to use parabricks STAR version
    star_index_legacy: Boolean                  // whether the supplied star_index was built with STAR 2.6.x and needs genomeParameters.txt upgraded to the 2.7.4a metadata schema
    hisat2_build_memory: String?                // memory threshold for HISAT2 index building with splice sites
    any_auto_strandedness: Boolean              // whether any sample in the input samplesheet declares strandedness 'auto', requiring a Salmon index for strandedness inference

    main:
    // Absent indices are Values holding null so every stream keeps one type.
    ch_no_path = channel.value(null as Path)

    //------------------------------------------------
    // 1) Determine which indices we actually want built
    //------------------------------------------------
    def prepare_tool_indices = (!skip_bbsplit ? ['bbsplit'] : []) +
        (ribo_removal_tool == 'sortmerna' ? ['sortmerna'] : []) +
        // If no index is provided, this subworkflow does not need to build an index as that is handled by the fastq_remove_rrna subworkflow.
        (ribo_removal_tool == 'bowtie2' && bowtie2_rrna_index ? ['bowtie2_rrna'] : []) +
        ((!skip_alignment && aligner) || aligner == 'star_rsem' ? [aligner as String] : []) +
        (!skip_pseudo_alignment && pseudo_aligner ? [pseudo_aligner as String] : []) +
        // needed to infer strandedness even without --pseudo_aligner salmon
        (any_auto_strandedness ? ['salmon'] : [])

    //---------------------------------------------------------
    // 2) BBSplit index: uses FASTA only if we generate from scratch
    //---------------------------------------------------------
    if ('bbsplit' in prepare_tool_indices && bbsplit_index && bbsplit_index.endsWith('.tar.gz')) {
        ch_bbsplit_index = UNTAR_BBSPLIT_INDEX(record(id: 'bbsplit_index', meta: [:], archive: file(bbsplit_index, checkIfExists: true))).map { r -> r.dir }
        ch_bbsplit_log   = ch_no_path
    } else if ('bbsplit' in prepare_tool_indices && bbsplit_index) {
        ch_bbsplit_index = channel.value(file(bbsplit_index, checkIfExists: true))
        ch_bbsplit_log   = ch_no_path
    } else if ('bbsplit' in prepare_tool_indices && fasta_provided) {
        // bbsplit_fasta_list is a 2 column csv: short_name,path_to_fasta
        def bbsplit_rows = file(bbsplit_fasta_list, checkIfExists: true).splitCsv()
        ch_bbsplit_fasta_list = channel.value(
            tuple(
                bbsplit_rows.collect { row -> row[0] },
                bbsplit_rows.collect { row -> file(row[1], checkIfExists: true) }
            )
        )

        // Index-only run: no reads
        def ch_bbsplit: Value<BbmapBbsplitResult> = BBMAP_BBSPLIT(
            channel.value(record(id: 'bbsplit_index', meta: [:], reads: [])),
            ch_no_path,
            ch_fasta,
            ch_bbsplit_fasta_list,
            true
        )
        ch_bbsplit_index = ch_bbsplit.map { r -> r.index }
        ch_bbsplit_log   = ch_bbsplit.map { r -> r.log }
    } else {
        ch_bbsplit_index = ch_no_path
        ch_bbsplit_log   = ch_no_path
    }

    //-------------------------------------------------------------
    // 3) SortMeRNA index
    //-------------------------------------------------------------
    if ('sortmerna' in prepare_tool_indices && sortmerna_index && sortmerna_index.endsWith('.tar.gz')) {
        ch_sortmerna_index = UNTAR_SORTMERNA_INDEX(record(id: 'sortmerna_index', meta: [:], archive: file(sortmerna_index, checkIfExists: true))).map { r -> r.dir }
    } else if ('sortmerna' in prepare_tool_indices && sortmerna_index) {
        ch_sortmerna_index = channel.value(file(sortmerna_index, checkIfExists: true))
    } else if ('sortmerna' in prepare_tool_indices) {
        // Index-only run: no reads
        def ch_sortmerna_built: Value<SortmernaResult> = SORTMERNA_INDEX(
            channel.value(record(id: 'rrna_refs', meta: [:], reads: [])),
            ch_rrna_fastas.collect().map { refs -> refs.toList() },
            ch_no_path
        )
        ch_sortmerna_index = ch_sortmerna_built.map { r -> r.index }
    } else {
        ch_sortmerna_index = ch_no_path
    }

    //-------------------------------------------------------------
    // 3b) Bowtie2 rRNA index - only handles untar
    //-------------------------------------------------------------
    // No need to build as that is handled in the fastq_remove_rrna subworkflow from nf-core
    if ('bowtie2_rrna' in prepare_tool_indices && bowtie2_rrna_index.endsWith('.tar.gz')) {
        ch_bowtie2_rrna_index = UNTAR_BOWTIE2_RRNA_INDEX(record(id: 'bowtie2_rrna_index', meta: [:], archive: file(bowtie2_rrna_index, checkIfExists: true))).map { r -> r.dir }
    } else if ('bowtie2_rrna' in prepare_tool_indices) {
        ch_bowtie2_rrna_index = channel.value(file(bowtie2_rrna_index, checkIfExists: true))
    } else {
        ch_bowtie2_rrna_index = ch_no_path
    }

    //----------------------------------------------------
    // 4) STAR index (e.g. for 'star_salmon') -> needs FASTA if built
    //----------------------------------------------------
    // ch_star_index_publish is the publishable form; for a legacy index this is the raw untarred index, not the upgraded copy STAR_ALIGN uses.
    // Pre-built STAR indices supplied by the user with star_index_legacy set (genomes-map opt-in for indices
    // built with STAR 2.6.x, e.g. AWS iGenomes) go through STAR_GENOMEPARAMS_UPGRADE to rewrite
    // `versionGenome 20201` and add the genomeType / genomeTransformType / genomeTransformVCF fields
    // that STAR 2.7.4a+ requires. Modern indices skip the upgrade entirely.
    def build_star = prepare_tool_indices.intersect(['star_salmon', 'star_rsem']) ? true : false
    if (build_star && use_parabricks_star && fasta_provided) {
        // Parabricks needs its own STAR index built with its bundled STAR version
        def ch_star_generated: Value<StarGenomegenerateResult> = PARABRICKS_STARGENOMEGENERATE(
            ch_fasta.map { item -> record(id: 'genome', meta: [:], fasta: item) },
            ch_gtf
        )
        ch_star_index         = ch_star_generated.map { r -> r.index }
        ch_star_index_publish = ch_star_index
    } else if (build_star && star_index) {
        if (star_index.endsWith('.tar.gz')) {
            ch_star_raw = UNTAR_STAR_INDEX(record(id: 'star_index', meta: [:], archive: file(star_index, checkIfExists: true)))
                .map { r -> record(id: r.id, meta: r.meta, index: r.dir) }
        } else {
            ch_star_raw = channel.value(record(id: 'star_index', meta: [:], index: file(star_index, checkIfExists: true)))
        }
        ch_star_index_publish = ch_star_raw.map { r -> r.index }
        if (star_index_legacy) {
            ch_star_index = STAR_GENOMEPARAMS_UPGRADE(ch_star_raw).map { r -> r.index }
        } else {
            ch_star_index = ch_star_index_publish
        }
    } else if (build_star && fasta_provided) {
        ch_star_generated = STAR_GENOMEGENERATE(
            ch_fasta.map { item -> record(id: 'genome', meta: [:], fasta: item) },
            ch_gtf
        )
        ch_star_index         = ch_star_generated.map { r -> r.index }
        ch_star_index_publish = ch_star_index
    } else {
        ch_star_index         = ch_no_path
        ch_star_index_publish = ch_no_path
    }

    //------------------------------------------------
    // 5) RSEM index -> needs FASTA & GTF if built
    //------------------------------------------------
    // ch_rsem_transcript_fasta is the incidental *transcripts.fa the index step also emits; unused by the pipeline, published under --save_reference.
    def build_rsem = 'star_rsem' in prepare_tool_indices
    if (build_rsem && rsem_index && rsem_index.endsWith('.tar.gz')) {
        ch_rsem_index            = UNTAR_RSEM_INDEX(record(id: 'rsem_index', meta: [:], archive: file(rsem_index, checkIfExists: true))).map { r -> r.dir }
        ch_rsem_transcript_fasta = ch_no_path
    } else if (build_rsem && rsem_index) {
        ch_rsem_index            = channel.value(file(rsem_index, checkIfExists: true))
        ch_rsem_transcript_fasta = ch_no_path
    } else if (build_rsem && fasta_provided && use_sentieon_star) {
        ch_rsem_reference        = SENTIEON_RSEM_PREPAREREFERENCE_GENOME(
            ch_fasta.map { fasta_file -> record(id: 'genome', meta: [id: 'genome'], fasta: fasta_file) }.combine(gtf: ch_gtf)
        )
        ch_rsem_index            = ch_rsem_reference.map { r -> r.index }
        ch_rsem_transcript_fasta = ch_rsem_reference.map { r -> r.transcript_fasta }
    } else if (build_rsem && fasta_provided) {
        ch_rsem_reference        = RSEM_PREPAREREFERENCE_GENOME(
            ch_fasta.map { fasta_file -> record(id: 'genome', meta: [id: 'genome'], fasta: fasta_file) }.combine(gtf: ch_gtf)
        )
        ch_rsem_index            = ch_rsem_reference.map { r -> r.index }
        ch_rsem_transcript_fasta = ch_rsem_reference.map { r -> r.transcript_fasta }
    } else {
        ch_rsem_index            = ch_no_path
        ch_rsem_transcript_fasta = ch_no_path
    }

    //---------------------------------------------------------
    // 6) HISAT2 index -> needs FASTA & GTF if built
    //---------------------------------------------------------
    def build_hisat2 = 'hisat2' in prepare_tool_indices
    if (build_hisat2 && splicesites) {
        ch_splicesites = channel.value(file(splicesites, checkIfExists: true))
    } else if (build_hisat2 && fasta_provided) {
        ch_splicesites = HISAT2_EXTRACTSPLICESITES(ch_gtf.map { item -> record(id: 'genome', meta: [:], gtf: item) }).map { r -> r.splicesites }
    } else {
        ch_splicesites = ch_no_path
    }
    if (build_hisat2 && hisat2_index && hisat2_index.endsWith('.tar.gz')) {
        ch_hisat2_index = UNTAR_HISAT2_INDEX(record(id: 'hisat2_index', meta: [:], archive: file(hisat2_index, checkIfExists: true))).map { r -> r.dir }
    } else if (build_hisat2 && hisat2_index) {
        ch_hisat2_index = channel.value(file(hisat2_index, checkIfExists: true))
    } else if (build_hisat2 && fasta_provided) {
        ch_hisat2_index = HISAT2_BUILD(
            ch_fasta
                .combine(ch_gtf)
                .combine(ch_splicesites)
                .map { fasta_file, gtf_file, ss_file -> record(id: 'genome', meta: [:], fasta: fasta_file, gtf: gtf_file, splicesites: ss_file) },
            hisat2_build_memory
        ).map { r -> r.index }
    } else {
        ch_hisat2_index = ch_no_path
    }

    //---------------------------------------------------------
    // 7) Bowtie2 index -> built from transcript FASTA for Salmon alignment mode
    //---------------------------------------------------------
    def build_bowtie2 = 'bowtie2_salmon' in prepare_tool_indices
    if (build_bowtie2 && bowtie2_index && bowtie2_index.endsWith('.tar.gz')) {
        ch_bowtie2_index = UNTAR_BOWTIE2_INDEX(record(id: 'bowtie2_index', meta: [:], archive: file(bowtie2_index, checkIfExists: true))).map { r -> r.dir }
    } else if (build_bowtie2 && bowtie2_index) {
        ch_bowtie2_index = channel.value(file(bowtie2_index, checkIfExists: true))
    } else if (build_bowtie2) {
        ch_bowtie2_index = BOWTIE2_BUILD(
            ch_transcript_fasta.map { fasta_file -> record(id: 'transcripts', meta: [id: 'transcripts'], fasta: fasta_file) }
        ).map { r -> r.index }
    } else {
        ch_bowtie2_index = ch_no_path
    }

    //------------------------------------------------------
    // 8) Salmon index -> can skip genome if transcript_fasta is enough
    //------------------------------------------------------
    if (salmon_index && salmon_index.endsWith('.tar.gz')) {
        ch_salmon_index = UNTAR_SALMON_INDEX(record(id: 'salmon_index', meta: [:], archive: file(salmon_index))).map { r -> r.dir }
    } else if (salmon_index) {
        ch_salmon_index = channel.value(file(salmon_index))
    } else if ('salmon' in prepare_tool_indices && fasta_provided) {
        // genome_fasta may be null (no decoys)
        def ch_salmon_built: Value<SalmonIndexResult> = SALMON_INDEX(
            ch_transcript_fasta
                .map { transcript_fasta_file -> record(id: 'salmon_index', meta: [:], transcript_fasta: transcript_fasta_file) }
                .combine(genome_fasta: ch_fasta)
        )
        ch_salmon_index = ch_salmon_built.map { r -> r.index }
    } else if ('salmon' in prepare_tool_indices) {
        ch_salmon_built = SALMON_INDEX(
            ch_transcript_fasta.map { item -> record(id: 'salmon_index', meta: [:], transcript_fasta: item, genome_fasta: null) }
        )
        ch_salmon_index = ch_salmon_built.map { r -> r.index }
    } else {
        ch_salmon_index = ch_no_path
    }

    //--------------------------------------------------
    // 9) Kallisto index -> only needs transcript FASTA
    //--------------------------------------------------
    if (kallisto_index && kallisto_index.endsWith('.tar.gz')) {
        ch_kallisto_index = UNTAR_KALLISTO_INDEX(record(id: 'kallisto_index', meta: [:], archive: file(kallisto_index))).map { r -> r.dir }
    } else if (kallisto_index) {
        ch_kallisto_index = channel.value(file(kallisto_index))
    } else if ('kallisto' in prepare_tool_indices) {
        def ch_kallisto_built: Value<KallistoIndexResult> = KALLISTO_INDEX(ch_transcript_fasta.map { item -> record(id: 'kallisto_index', meta: [:], fasta: item) })
        ch_kallisto_index = ch_kallisto_built.map { r -> r.index }
    } else {
        ch_kallisto_index = ch_no_path
    }

    //--------------------------------------------------
    // 10) Streamed index artifacts, one record per producer
    //--------------------------------------------------
    // Each index is an independent, optionally user-supplied producer with no shared
    // key, so they stream as tagged records via mix() rather than fusing into one row.
    // Absent indices hold null and are dropped by the trailing filter.
    ch_indices = ch_star_index_publish.flatMap { p -> [ record(kind: 'star', file: p) ] }
        .mix(ch_rsem_index.flatMap              { p -> [ record(kind: 'rsem', file: p) ] })
        .mix(ch_rsem_transcript_fasta.flatMap   { p -> [ record(kind: 'rsem_transcript_fasta', file: p) ] })
        .mix(ch_hisat2_index.flatMap            { p -> [ record(kind: 'hisat2', file: p) ] })
        .mix(ch_splicesites.flatMap             { p -> [ record(kind: 'hisat2_splicesites', file: p) ] })
        .mix(ch_bowtie2_index.flatMap           { p -> [ record(kind: 'bowtie2', file: p) ] })
        .mix(ch_salmon_index.flatMap            { p -> [ record(kind: 'salmon', file: p) ] })
        .mix(ch_kallisto_index.flatMap          { p -> [ record(kind: 'kallisto', file: p) ] })
        .mix(ch_bbsplit_index.flatMap           { p -> [ record(kind: 'bbsplit', file: p) ] })
        .mix(ch_bbsplit_log.flatMap             { p -> [ record(kind: 'bbsplit_log', file: p) ] })
        .mix(ch_sortmerna_index.flatMap         { p -> [ record(kind: 'sortmerna', file: p) ] })
        .mix(ch_bowtie2_rrna_index.flatMap      { p -> [ record(kind: 'bowtie2_rrna', file: p) ] })
        .filter { r -> taskOutputOrNull(r.file) != null }

    // The index streams stay separate named emits: RNASEQ consumes each index independently,
    // and no per-sample key exists to fuse them into one record.
    emit:
    splicesites:        Value<Path?>            = ch_splicesites        // genome.splicesites.txt, null when absent
    bbsplit_index:      Value<Path?>            = ch_bbsplit_index      // bbsplit/index/, null when absent
    sortmerna_index:    Value<Path?>            = ch_sortmerna_index    // sortmerna/index/, null when absent
    bowtie2_rrna_index: Value<Path?>            = ch_bowtie2_rrna_index // bowtie2/index/, null when absent
    star_index:         Value<Path?>            = ch_star_index         // star/index/, null when absent
    rsem_index:         Value<Path?>            = ch_rsem_index         // rsem/index/, null when absent
    hisat2_index:       Value<Path?>            = ch_hisat2_index       // hisat2/index/, null when absent
    bowtie2_index:      Value<Path?>            = ch_bowtie2_index      // bowtie2/index/, null when absent
    salmon_index:       Value<Path?>            = ch_salmon_index       // salmon/index/, null when absent
    kallisto_index:     Value<Path?>            = ch_kallisto_index     // kallisto/index/, null when absent
    indices:            Channel<GenomeArtifact> = ch_indices            // one record per index/log actually built or supplied
}
