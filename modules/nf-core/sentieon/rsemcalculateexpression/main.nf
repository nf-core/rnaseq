nextflow.enable.types = true

include { ReadsInput } from '../../types'
include { RsemQuantSample } from '../../rsem/calculateexpression/main'

process SENTIEON_RSEMCALCULATEEXPRESSION {
    tag "$sample.meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39a3e1a85912520836ad054c8ac0497b463bb5170e0e907183dbd08509dad997/data' :
        'community.wave.seqera.io/library/rsem_sentieon:3e4315fa0b636313' }"

    input:
    sample: ReadsInput  // FASTQ files or BAM file for --alignments mode
    index: Path

    output:
    record(
        id:                sample.id,
        meta:              sample.meta,
        counts_gene:       file("*.genes.results"),
        counts_transcript: file("*.isoforms.results"),
        stat:              file("*.stat"),
        log:               file("*.log", optional: true),
        bam_star:          file("*.STAR.genome.bam", optional: true),
        bam_genome:        file("${prefix}.genome.bam", optional: true),
        bam_transcript:    file("${prefix}.transcript.bam", optional: true)
    ) as RsemQuantSample

    topic:
    tuple(task.process, 'rsem', eval('rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g"')) >> 'versions'
    tuple(task.process, 'star', eval('STAR --version | sed -e "s/STAR_//g"')) >> 'versions'
    tuple(task.process, 'sentieon', eval('sentieon driver --version 2>&1 | sed -e "s/sentieon-genomics-//g"')) >> 'versions'

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${sample.meta.id}"

    def strandedness = ''
    if (sample.meta.strandedness == 'forward') {
        strandedness = '--strandedness forward'
    } else if (sample.meta.strandedness == 'reverse') {
        strandedness = '--strandedness reverse'
    }

    // Detect if input is BAM file(s); a lone file arrives as a Path, which iterates over its name components
    def reads_names = sample.reads instanceof Path ? "${sample.reads}" : sample.reads.join(' ')
    def reads_count = sample.reads instanceof Path ? 1 : sample.reads.size()
    def is_bam = reads_names.toLowerCase().endsWith('.bam')
    def alignment_mode = is_bam ? '--alignments' : ''

    // Use metadata for paired-end detection if available, otherwise empty (auto-detect)
    def paired_end = sample.meta.containsKey('single_end') ? (sample.meta.single_end ? "" : "--paired-end") : "unknown"

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""

    """
    $sentieonLicense

    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`

    # Create symlink to sentieon in PATH
    ln -sf \$(which sentieon) ./STAR
    export PATH=".:\$PATH"

    # Use metadata-based paired-end detection, or auto-detect if no metadata provided
    PAIRED_END_FLAG="$paired_end"
    if [ "${paired_end}" == "unknown" ]; then
        # Auto-detect only if no metadata provided
        if [ "${is_bam}" == "true" ]; then
            samtools flagstat $reads_names | grep -q 'paired in sequencing' && PAIRED_END_FLAG="--paired-end"
        else
            [ ${reads_count} -gt 1 ] && PAIRED_END_FLAG="--paired-end"
        fi
    fi

    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --temporary-folder ./tmp/ \\
        $alignment_mode \\
        $strandedness \\
        \$PAIRED_END_FLAG \\
        $args \\
        $reads_names \\
        \$INDEX \\
        $prefix
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def is_bam = (sample.reads instanceof Path ? "${sample.reads}" : sample.reads.join(' ')).toLowerCase().endsWith('.bam')
    """
    # Create symlink to sentieon in PATH for version detection
    ln -sf \$(which sentieon) ./STAR
    export PATH=".:\$PATH"

    touch ${prefix}.genes.results
    touch ${prefix}.isoforms.results
    touch ${prefix}.stat
    touch ${prefix}.log

    # Only create STAR BAM output when not in alignment mode
    if [ "${is_bam}" == "false" ]; then
        touch ${prefix}.STAR.genome.bam
    fi

    touch ${prefix}.genome.bam
    touch ${prefix}.transcript.bam
    """
}
