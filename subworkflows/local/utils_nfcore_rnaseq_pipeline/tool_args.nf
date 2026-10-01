//
// Tool arguments that follow from pipeline params. The workflows pass them to the processes that need
// them, so that a pipeline that includes this one controls them through its params record and a change
// only affects the tasks that use the value. Plain functions of values, so they can be shared by the
// processes of a tool.
//

// UMI-tools extract options
def umiExtractArgs(params) {
    return [
        params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : '',
        params.umitools_bc_pattern     ? "--bc-pattern='${params.umitools_bc_pattern}'" : '',
        params.umitools_bc_pattern2    ? "--bc-pattern2='${params.umitools_bc_pattern2}'" : '',
        params.umitools_umi_separator  ? "--umi-separator='${params.umitools_umi_separator}'" : ''
    ].join(' ').trim()
}

// Merges extra options into default options, so that the extra value of an option replaces the default one
// (a repeated option is an error for some tools). Options start where a token begins with `prefix`.
def mergeOptions(List defaults, String extras, String option_start) {
    if (!extras) {
        return defaults
    }
    return (defaults.join(' ') + ' ' + extras)
        .split("\\s(?=${option_start})")
        .collectEntries { String token ->
            def parts = token.trim().split(/\s+/, 2)
            [(parts[0]): parts.size() > 1 ? parts[1] : '']
        }
        .collect { k, v -> v ? "${k} ${v}" : k }
}

// Read group values: the per-sample ones of the samplesheet take priority over those of the params
def readGroupCenter(params, meta) {
    return meta.seq_center ?: params.seq_center
}
def readGroupPlatform(params, meta) {
    return meta.seq_platform ?: params.seq_platform
}

// STAR options of the process that runs the alignment: 'star', 'sentieon' or 'parabricks'
def starAlignArgs(params, tool, meta) {
    def isPbrun = tool == 'parabricks'
    def isSentieon = tool == 'sentieon'
    def quantifier = params.aligner == 'star_rsem' ? 'rsem' : 'salmon'

    def args = []

    // Common args - pbrun uses kebab-case equivalents
    if (isPbrun) {
        args += [
            '--quantMode TranscriptomeSAM',
            '--out-sam-attributes NH HI AS NM MD',
            '--read-files-command zcat',
            "--read-group-sm ${meta.id}",
            "--read-group-id-prefix ${meta.id}"
        ]
    } else {
        args += [
            '--quantMode TranscriptomeSAM',
            '--outSAMtype BAM Unsorted',
            '--outSAMattributes NH HI AS NM MD',
            '--readFilesCommand zcat'
        ]
    }

    // Quantifier-specific args
    //
    // NOTE: pbrun rna_fq2bam is based on STAR 2.7.2a and does NOT support:
    //   --outFilterType BySJout  : no pbrun equivalent; pbrun uses its own filtering defaults
    //   --sjdbScore 1            : no pbrun equivalent; affects junction scoring priority
    //   --quantTranscriptomeBan Singleend / --quantTranscriptomeSAMoutput BanSingleEnd :
    //       no pbrun equivalent; Salmon handles mixed single/paired records gracefully
    //   --runRNGseed 0           : no pbrun equivalent; pbrun uses deterministic alignment selection
    //
    if (quantifier == 'rsem') {
        if (isPbrun) {
            args += [
                '--out-sam-unmapped Within',
                '--max-out-filter-multimap 20',
                '--max-out-filter-mismatch 999',
                '--max-out-filter-mismatch-ratio 0.04',
                '--min-intron-size 20',
                '--max-intron-size 1000000',
                '--max-align-mates-gap 1000000',
                '--min-align-sj-overhang 8',
                '--min-align-sjdb-overhang 1'
            ]
        } else {
            args += [
                '--outSAMunmapped Within',
                '--outFilterType BySJout',
                '--outFilterMultimapNmax 20',
                '--outFilterMismatchNmax 999',
                '--outFilterMismatchNoverLmax 0.04',
                '--alignIntronMin 20',
                '--alignIntronMax 1000000',
                '--alignMatesGapMax 1000000',
                '--alignSJoverhangMin 8',
                '--alignSJDBoverhangMin 1',
                '--sjdbScore 1'
            ]
        }
    } else {
        if (isPbrun) {
            args += [
                '--two-pass-mode Basic',
                '--max-out-filter-multimap 20',
                '--min-align-sjdb-overhang 1',
                '--out-sam-strand-field intronMotif'
            ]
        } else {
            args += [
                '--twopassMode Basic',
                '--runRNGseed 0',
                '--outFilterMultimapNmax 20',
                '--alignSJDBoverhangMin 1',
                '--outSAMstrandField intronMotif',
                // Sentieon's STAR build predates 2.7.7a so it expects the
                // pre-2.7.7a --quantTranscriptomeBan Singleend spelling. Drop this
                // branch (and isSentieon) once Sentieon's bundled STAR is >= 2.7.7a.
                isSentieon ? '--quantTranscriptomeBan Singleend' : '--quantTranscriptomeSAMoutput BanSingleEnd'
            ]
        }
    }

    // Unmapped reads output
    if (params.save_unaligned || (params.contaminant_screening && params.contaminant_screening_input == 'unmapped')) {
        args += isPbrun
            ? ['--out-reads-unmapped Fastx']
            : ['--outReadsUnmapped Fastx']
    }

    // Prokaryotic: disable spliced alignment
    if (params.prokaryotic) {
        if (isPbrun) {
            args += ['--max-intron-size 1']
        } else {
            args += ['--sjdbGTFfeatureExon CDS', '--alignIntronMax 1']
        }
    }

    // Merge user-supplied extra args, deduplicating so user values
    // override pipeline defaults (prevents STAR "duplicate parameter" error)
    args = mergeOptions(args, params.extra_star_align_args, '--')

    // Read group tags
    def seqCenter   = readGroupCenter(params, meta)
    def seqPlatform = readGroupPlatform(params, meta)
    if (isPbrun) {
        // pbrun supports --read-group-pl but has no --read-group-cn flag
        if (seqPlatform) args << "--read-group-pl ${seqPlatform}"
    } else {
        def rgParts = ["'ID:${meta.id}'", "'SM:${meta.id}'"]
        if (seqPlatform) rgParts << "'PL:${seqPlatform}'"
        if (seqCenter) rgParts << "'CN:${seqCenter}'"
        args += ["--outSAMattrRGline ${rgParts.join(' ')}"]
    }

    return args.join(' ')
}

// HISAT2 alignment options
def hisat2AlignArgs(params, meta) {
    def argsList = ['--met-stderr', '--new-summary', '--dta']

    // Strandedness
    if (meta.strandedness == 'forward') {
        argsList << (meta.single_end ? '--rna-strandness F' : '--rna-strandness FR')
    } else if (meta.strandedness == 'reverse') {
        argsList << (meta.single_end ? '--rna-strandness R' : '--rna-strandness RF')
    }

    // Merge extras so user values override defaults (last occurrence of a flag wins).
    // Must run before RG args since `--rg` repeats and would collapse.
    argsList = mergeOptions(argsList, params.extra_hisat2_align_args, '--')

    // Read group tags
    def seqCenter   = readGroupCenter(params, meta)
    def seqPlatform = readGroupPlatform(params, meta)
    argsList << "--rg-id ${meta.id}"
    argsList << "--rg SM:${meta.id}"
    if (seqPlatform) argsList << "--rg PL:${seqPlatform}"
    if (seqCenter) argsList << "--rg CN:${seqCenter}"

    return argsList.join(' ')
}

// Bowtie2 alignment options
def bowtie2AlignArgs(params, meta) {
    def args = [
        '--very-sensitive',
        '--no-discordant',
        '-k 200'
    ]

    // Merge extras so user values override defaults (e.g. `-k N`).
    // Must run before RG args since `--rg` repeats and would collapse.
    args = mergeOptions(args, params.extra_bowtie2_align_args, '-')

    // Read group tags
    def seqCenter   = readGroupCenter(params, meta)
    def seqPlatform = readGroupPlatform(params, meta)
    def rgParts = ["--rg-id ${meta.id}", "--rg SM:${meta.id}"]
    if (seqPlatform) rgParts << "--rg PL:${seqPlatform}"
    if (seqCenter) rgParts << "--rg CN:${seqCenter}"
    args += rgParts

    return args.join(' ')
}

// samtools index options
def samtoolsIndexArgs(params) {
    return params.bam_csi_index ? '-c' : ''
}
